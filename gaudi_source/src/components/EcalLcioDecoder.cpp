/*
 * EcalLcioDecoder: a second front-end to the SAME decode backend as
 * EcalRawDecoder, reading the EUDAQ LCIO files in Data/eudaq/ instead of the raw
 * SL byte stream, and writing the identical `siwecaldecoded` tree.
 *
 * WHY. 13 of the 14 th210 muon runs (eudaq 146-178) exist only as the EUDAQ
 * container the raw decoder cannot parse -- but they also exist, already decoded,
 * as LCIO in Data/eudaq/ROC_run_<N>_tp.slcio. This component reads those so th210
 * (and any other LCIO-only run) enters the existing event-building/calibration
 * pipeline unchanged.
 *
 * NOT A REIMPLEMENTATION. It fills the same k4siwecal::Acquisition struct the raw
 * decoder fills, then calls the same computeCorrectedBcid + tagBadBcid and binds the
 * tree through the same bindAcquisitionBranches. The only new logic is reading the
 * LCGenericObject records; corrected_bcid, badbcid and the tree schema are the
 * already-validated backend, so an LCIO-decoded run cannot drift from a raw-decoded
 * one.
 *
 * FIELD MAP, pinned by exact integer match against the raw decode of run_000254
 * (which exists in both formats: 71052 frames, 100% agreement):
 *
 *   collection EUDAQDataSiECAL, one LCGenericObject per (CycleNr, Layer, SkirocID,
 *   sca); header ints [CycleNr, BunchXID, sca, Layer, SkirocID, NChannels] then, per
 *   channel, 6 interleaved ints:
 *       [0] adc_low   [1] hitbit   [2] autogainbit   [3] adc_high   [4] hitbit(dup)   [5] autogainbit(dup)
 *   CycleNr -> acqNumber,  BunchXID -> bcid,  Layer -> slab,  SkirocID -> chip.
 *   nhits per (slab,chip,sca) = sum of the hit bits (== raw nhits, verified exactly).
 */
#include "k4SiWEcalReco/AcquisitionTreeIO.h"
#include "k4SiWEcalReco/RunSettings.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TFile.h"
#include "TTree.h"

#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCGenericObject.h"
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"

#include <memory>
#include <string>
#include <vector>

using k4siwecal::kChannelsInSkiroc;
using k4siwecal::kScasInSkiroc;
using k4siwecal::kSkirocsPerAsu;
using k4siwecal::kSlbDepth;

struct EcalLcioDecoder final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    if (m_inputFile.value().empty()) {
      return error() << "InputFile (a .slcio path) not set" << endmsg, StatusCode::FAILURE;
    }
    if (m_outputFile.value().empty()) {
      return error() << "OutputFile not set" << endmsg, StatusCode::FAILURE;
    }

    // Same ZSTD-5 as the raw decoder -- see k4siwecal::kDecodedFileCompression.
    auto fout = std::unique_ptr<TFile>(
        TFile::Open(m_outputFile.value().c_str(), "RECREATE", "", k4siwecal::kDecodedFileCompression));
    if (!fout || fout->IsZombie()) {
      return error() << "Cannot create output file: " << m_outputFile.value() << endmsg, StatusCode::FAILURE;
    }
    fout->cd();
    auto* tree = new TTree(m_treeName.value().c_str(), m_treeName.value().c_str());

    auto acq = std::make_unique<k4siwecal::Acquisition>();
    k4siwecal::RunSettings runSettings;
    if (!m_runSettingsFile.value().empty()) {
      runSettings = k4siwecal::parseRunSettings(m_runSettingsFile.value());
    }
    k4siwecal::bindAcquisitionBranches(*tree, *acq, runSettings);

    std::unique_ptr<IO::LCReader> reader(IOIMPL::LCFactory::getInstance()->createLCReader());
    try {
      reader->open(m_inputFile.value());
    } catch (const std::exception& e) {
      return error() << "Cannot open LCIO file " << m_inputFile.value() << ": " << e.what() << endmsg,
             StatusCode::FAILURE;
    }

    // Finalise one acquisition: fill the slab bookkeeping the LCIO records do not
    // carry per-object, then run the same overflow-correction + retrigger tagging
    // the raw decoder runs, and stream the entry out.
    long long nWritten = 0;
    const auto flush = [&]() {
      for (int slab = 0; slab < kSlbDepth; ++slab) {
        bool present = false;
        for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
          if (acq->chipId[slab][chip] >= 0) {
            present = true;
            break;
          }
        }
        if (present) {
          acq->slboardId[slab] = slab;  // Layer IS the slab id; identity, as the raw
          acq->slot[slab] = -1;         // decoder leaves slot unknown (SlbFrameDecoder.h)
        }
      }
      acq->nSlboards = kSlbDepth;  // recordFrame sets this to kSlbDepth on every frame
      k4siwecal::computeCorrectedBcid(*acq);
      k4siwecal::tagBadBcid(*acq);
      tree->Fill();
      ++nWritten;
    };

    int currentCycle = -1;
    bool haveCycle = false;
    long long nEvents = 0;
    EVENT::LCEvent* evt = nullptr;
    while ((evt = reader->readNextEvent()) != nullptr) {
      ++nEvents;
      EVENT::LCCollection* col = nullptr;
      try {
        col = evt->getCollection(m_collectionName.value());
      } catch (const std::exception&) {
        continue;  // event with no SiECAL data (e.g. a slow-control-only record)
      }
      const int nObj = col->getNumberOfElements();
      for (int k = 0; k < nObj; ++k) {
        auto* obj = dynamic_cast<EVENT::LCGenericObject*>(col->getElementAt(k));
        if (obj == nullptr) continue;
        const int ni = obj->getNInt();
        if (ni < 6) continue;
        const int cyc = obj->getIntVal(0);
        const int bx = obj->getIntVal(1);
        const int sca = obj->getIntVal(2);
        const int layer = obj->getIntVal(3);
        const int chip = obj->getIntVal(4);
        const int nc = obj->getIntVal(5);
        // Only well-formed full-channel records; the header must account for exactly
        // 6 + 6*NChannels ints, matching the validated interleaved layout.
        if (nc != kChannelsInSkiroc || ni != 6 + 6 * nc) continue;
        if (layer < 0 || layer >= kSlbDepth || chip < 0 || chip >= kSkirocsPerAsu || sca < 0 ||
            sca >= kScasInSkiroc)
          continue;

        if (!haveCycle) {
          k4siwecal::initAcquisition(*acq);
          currentCycle = cyc;
          acq->acqNumber = cyc;
          haveCycle = true;
        } else if (cyc != currentCycle) {
          flush();
          if (m_maxAcq.value() > 0 && nWritten >= m_maxAcq.value()) {
            haveCycle = false;
            break;
          }
          k4siwecal::initAcquisition(*acq);
          currentCycle = cyc;
          acq->acqNumber = cyc;
        }

        acq->chipId[layer][chip] = chip;
        if (sca + 1 > acq->numCol[layer][chip]) acq->numCol[layer][chip] = sca + 1;
        acq->bcid[layer][chip][sca] = bx;
        int nhit = 0;
        for (int ch = 0; ch < nc; ++ch) {
          const int base = 6 + ch * 6;
          const int adcLow = obj->getIntVal(base + 0);
          const int hit = obj->getIntVal(base + 1);
          const int gain = obj->getIntVal(base + 2);
          const int adcHigh = obj->getIntVal(base + 3);
          acq->adcLow[layer][chip][sca][ch] = adcLow;
          acq->adcHigh[layer][chip][sca][ch] = adcHigh;
          acq->hitbitLow[layer][chip][sca][ch] = hit;
          acq->hitbitHigh[layer][chip][sca][ch] = hit;
          acq->autogainbitLow[layer][chip][sca][ch] = gain;
          acq->autogainbitHigh[layer][chip][sca][ch] = gain;
          nhit += hit;
        }
        acq->nhits[layer][chip][sca] = nhit;
      }
      if (m_maxAcq.value() > 0 && nWritten >= m_maxAcq.value()) break;
    }
    if (haveCycle && !(m_maxAcq.value() > 0 && nWritten >= m_maxAcq.value())) flush();
    reader->close();

    fout->cd();
    fout->Write(nullptr, TObject::kOverwrite);
    fout->Close();
    info() << "EcalLcioDecoder: wrote " << nWritten << " acquisition cycles from " << nEvents
           << " LCIO event(s) to " << m_outputFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override { return StatusCode::SUCCESS; }

 private:
  Gaudi::Property<std::string> m_inputFile{this, "InputFile", "", "Input EUDAQ .slcio file for one run"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "", "Output ROOT file (siwecaldecoded tree)"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Output TTree name"};
  Gaudi::Property<std::string> m_collectionName{this, "CollectionName", "EUDAQDataSiECAL",
                                                 "LCGenericObject collection holding the SiECAL frames"};
  Gaudi::Property<std::string> m_runSettingsFile{
      this, "RunSettingsFile", "",
      "Path to Run_Settings.txt (optional); supplies thresholdDac/holdDelay/fsPeakTime/... exactly as for "
      "a raw-decoded run. Empty leaves those branches at their -1 sentinel."};
  Gaudi::Property<int> m_maxAcq{this, "MaxAcq", -1, "Stop after this many acquisitions (<=0 = all); for testing"};
};

DECLARE_COMPONENT(EcalLcioDecoder)
