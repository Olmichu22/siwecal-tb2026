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

/* ROOT compression setting for every `siwecaldecoded` file: ZSTD level 5
 * (algorithm 5, level 5 -> 505). Both decoders pass it to TFile::Open, so the
 * two front-ends stay byte-comparable.
 *
 * NOT a micro-optimisation. An Acquisition is 17.4 GB of uncompressed tree data
 * per ~3000-entry chunk, because the arrays are fixed [15][16][15][64] and go
 * out mostly empty -- and ROOT's default (ZLIB level 1, setting 101) is
 * remarkably bad at that shape. Measured on one real chunk of
 * TB2026CERN_eudaq_run_000285 (2026-07-25):
 *
 *     ZLIB-1 (the old default) 185.4 MB   <- what production used to write
 *     ZLIB-6                    25.5 MB
 *     LZ4-4                     78.4 MB
 *     ZSTD-5                    10.0 MB   <- 18.6x smaller, ~27 s to write
 *     ZSTD-9                     8.3 MB   (22.3x, but ~69 s)
 *     LZMA-5                     8.7 MB   (21.2x, but ~214 s)
 *
 * Control: rewriting the same tree with ZLIB-1 again gives 184.8 MB (1.00x), so
 * the gain is the algorithm, not the rewrite. ZSTD-5 is the knee of the curve --
 * ZSTD-9 buys 17% more space for 2.5x the CPU. This took a full campaign from
 * ~2.4 TB to ~130 GB, which is the difference between fitting in the EOS quota
 * and aborting the DAG halfway through it (it did, on 2026-07-25).
 *
 * Transparent to every reader: ROOT decompresses on open, and ZSTD has been in
 * ROOT since 6.20, far below the key4hep stack this builds against.
 */
inline constexpr int kDecodedFileCompression = 505;

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
