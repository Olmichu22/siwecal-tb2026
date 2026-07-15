/*
 * AcquisitionTreeIO: the single definition of the `siwecaldecoded` tree schema.
 *
 * Both decoders bind their output tree through bindAcquisitionBranches -- the raw
 * SL decoder (EcalRawDecoder) and the EUDAQ-LCIO decoder (EcalLcioDecoder) -- so the
 * two front-ends cannot drift into producing subtly different trees. A branch added
 * or renamed here reaches both at once; that shared schema is exactly why an event
 * built from an LCIO-converted th210 run and one built from a raw-decoded run are
 * read by the same downstream code with no special-casing.
 *
 * Extracted verbatim from EcalRawDecoder::bindBranches; the only change is that it is
 * now a free function in a header instead of a private static, so it can be reused.
 */
#pragma once

#include "k4SiWEcalReco/RunSettings.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"  // Acquisition, kSlbDepth, kSkirocsPerAsu, ...

#include "TString.h"  // Form
#include "TTree.h"

namespace k4siwecal {

inline void bindAcquisitionBranches(TTree& tree, Acquisition& acq, RunSettings& runSettings) {
  tree.Branch("acqNumber", &acq.acqNumber, "acqNumber/I");
  tree.Branch("n_slboards", &acq.nSlboards, "n_slboards/I");
  tree.Branch("slot", acq.slot, Form("slot[%d]/I", kSlbDepth));
  tree.Branch("slboard_id", acq.slboardId, Form("slboard_id[%d]/I", kSlbDepth));
  tree.Branch("chipid", acq.chipId, Form("chipid[%d][%d]/I", kSlbDepth, kSkirocsPerAsu));
  tree.Branch("nColumns", acq.numCol, Form("nColumns[%d][%d]/I", kSlbDepth, kSkirocsPerAsu));
  tree.Branch("startACQ", acq.startAcq, Form("startACQ[%d]/F", kSlbDepth));
  tree.Branch("rawTSD", acq.rawTsd, Form("rawTSD[%d]/I", kSlbDepth));
  tree.Branch("TSD", acq.tsd, Form("TSD[%d]/F", kSlbDepth));
  tree.Branch("rawAVDD0", acq.rawAvdd0, Form("rawAVDD0[%d]/I", kSlbDepth));
  tree.Branch("rawAVDD1", acq.rawAvdd1, Form("rawAVDD1[%d]/I", kSlbDepth));
  tree.Branch("AVDD0", acq.avdd0, Form("AVDD0[%d]/F", kSlbDepth));
  tree.Branch("AVDD1", acq.avdd1, Form("AVDD1[%d]/F", kSlbDepth));
  tree.Branch("bcid", acq.bcid, Form("bcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
  tree.Branch("corrected_bcid", acq.correctedBcid,
              Form("corrected_bcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
  tree.Branch("badbcid", acq.badbcid, Form("badbcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
  tree.Branch("nhits", acq.nhits, Form("nhits[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
  tree.Branch("adc_low", acq.adcLow,
              Form("adc_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
  tree.Branch("adc_high", acq.adcHigh,
              Form("adc_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
  tree.Branch("autogainbit_low", acq.autogainbitLow,
              Form("autogainbit_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
  tree.Branch("autogainbit_high", acq.autogainbitHigh,
              Form("autogainbit_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
  tree.Branch("hitbit_low", acq.hitbitLow,
              Form("hitbit_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
  tree.Branch("hitbit_high", acq.hitbitHigh,
              Form("hitbit_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));

  tree.Branch("acqWindowMs", &runSettings.acqWindowMs, "acqWindowMs/F");
  tree.Branch("delayBetweenCycleMs", &runSettings.delayBetweenCycleMs, "delayBetweenCycleMs/F");
  tree.Branch("thresholdDac", &runSettings.thresholdDac, "thresholdDac/I");
  tree.Branch("holdDelay", &runSettings.holdDelay, "holdDelay/I");
  tree.Branch("fsPeakTime", &runSettings.fsPeakTime, "fsPeakTime/I");
  tree.Branch("gainSelectionThreshold", &runSettings.gainSelectionThreshold, "gainSelectionThreshold/I");
}

}  // namespace k4siwecal
