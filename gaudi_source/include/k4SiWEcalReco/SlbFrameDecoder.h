/*
 * SlbFrameDecoder: pure C++ (no ROOT/Gaudi dependency) port of the SLB raw
 * binary decoder from the external reference tool
 * /eos/experiment/drdcalo/siw-ecal/suehara/converter_SLB/SLBraw2ROOT.cc --
 * the version that actually produced the official `Data/rundata_converted`
 * files. (Older public checkouts of SiWECAL-TB-analysis carry a
 * `cycleIDDecoding` that indexes at `2*n+1` instead of `2*n+3`; do not use
 * them as the reference.) Kept ROOT/Gaudi-independent so it can be
 * unit-tested with synthetic byte sequences, same separation principle as
 * EcalShowerVars.h.
 *
 * This is a byte-for-byte-faithful port (same index arithmetic, same bit
 * masks, same channel-order reversal) of a delicate, already-validated
 * bit-level protocol decoder. Do not "clean up" the arithmetic without
 * re-verifying against the reference tool's output on a real raw run.
 *
 * One deliberate deviation: the reference's startAcqTimeStamp decode loop
 * (SLBraw2ROOT.cc:291-297) shifts by `(30-2*n)` with `n` running 16..31,
 * i.e. by negative amounts -2..-32 -- undefined behaviour in C++, and not
 * reproducible deterministically across compilers. Since it's a slow-control
 * diagnostic (not read by AcquisitionReader or any calibration/physics code),
 * this port resets the shift counter for that field the same way the
 * TSD/AVDD0/AVDD1/transmitID loops already do in the reference, instead of
 * replicating undefined behaviour.
 */
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <istream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace k4siwecal {

constexpr int kSlbDepth = 15;
constexpr int kSkirocsPerAsu = 16;
constexpr int kScasInSkiroc = 15;
constexpr int kChannelsInSkiroc = 64;
constexpr int kSingleSkirocEventSize = 129 * 2;  // bytes
constexpr int kBcidThresDefault = 3;

// ---------------------------------------------------------------------------
// Convert_FromGrayToBinary (SLBraw2ROOT.cc:41-66), ported verbatim.
inline int convertGrayToBinary(int grayValue, int nbOfBits) {
  int grayBit = 1 << (nbOfBits - 1);
  int binary = grayValue & grayBit;
  int binBit = binary;
  while ((grayBit >>= 1) != 0) {
    binBit >>= 1;
    binary |= binBit ^ (grayValue & grayBit);
    binBit = binary & grayBit;
  }
  return binary;
}

// ---------------------------------------------------------------------------
// One raw frame = one SKIROC chip's readout for one acquisition cycle,
// containing up to kScasInSkiroc SCA slots. Port of the per-frame scalar and
// per-SCA/channel outputs of SLBraw2ROOT::DecodeRawFrame.
struct DecodedFrame {
  int coreDaughterIndex = -1;
  int slabAdd = -1;
  int chipId = -1;
  int asuIndex = -1;
  int skirocIndex = -1;
  int nSca = 0;  // nbOfSingleSkirocEventsInFrame

  int cycleId = 0;
  double startAcqTimeStamp = 0.;
  int rawTsd = 0;
  int rawAvdd0 = 0;
  int rawAvdd1 = 0;
  int transmitId = 0;
  float temperature = 0.f;
  float avdd0 = 0.f;
  float avdd1 = 0.f;

  int bcid[kScasInSkiroc] = {0};
  int nhits[kScasInSkiroc] = {0};
  int adcLow[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
  int adcHigh[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
  int gainLow[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
  int gainHigh[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
  int hitLow[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
  int hitHigh[kScasInSkiroc][kChannelsInSkiroc] = {{0}};
};

// Cycle-id key used to GROUP frames while buffering (SLBraw2ROOT.cc's
// cycleIDDecoding, called at ReadFile:570 on the frame bytes BEFORE the
// 2-byte core/slab prefix is stripped). The `+3` offset skips that prefix, so
// this reads exactly the same cycleId field `decodeFrame` reads on the
// prefix-stripped bytes -- reading at `2*n+1` instead would return
// `cycleId >> 2`, collapsing four consecutive acquisition cycles into one
// buffered Acquisition (whose per-(slab,chip) slots then overwrite each
// other). `frameBytes` must be the full frame including its 2-byte prefix,
// exactly as buffered.
inline int decodeCycleId(const std::vector<unsigned char>& frameBytes) {
  int result = 0;
  for (int n = 0; n < 16; ++n) {
    result += static_cast<unsigned int>(((frameBytes.at(2 * n + 3) & 0xC0) >> 6) << (30 - 2 * n));
  }
  return result;
}

// Port of SLBraw2ROOT::DecodeRawFrame (SLBraw2ROOT.cc:246-425). `frameBytes`
// is the full frame including its 2-byte core/slab prefix (same buffer given
// to decodeCycleId); this function makes its own copy to erase the prefix,
// exactly as the reference does internally.
inline bool decodeFrame(std::vector<unsigned char> frameBytes, DecodedFrame& out) {
  out = DecodedFrame{};
  out.coreDaughterIndex = frameBytes.at(0);
  out.slabAdd = frameBytes.at(1);

  frameBytes.erase(frameBytes.begin());
  frameBytes.erase(frameBytes.begin());

  const int actualDataSize = static_cast<int>(frameBytes.size());
  out.chipId = frameBytes.at(actualDataSize - 2 - 2);
  out.asuIndex = out.chipId / kSkirocsPerAsu;
  out.skirocIndex = out.chipId - out.asuIndex * kSkirocsPerAsu;
  out.nSca = (actualDataSize - 2 - 2) / kSingleSkirocEventSize;

  // Metadata block (DecodeRawFrame:283-345). cycleId re-derived here (on the
  // prefix-stripped bytes) is the value actually stored in the output tree;
  // it differs by design from decodeCycleId()'s buffering key (see above).
  int cycleId = 0;
  for (int n = 0; n < 16; ++n) {
    cycleId += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (30 - 2 * n));
  }
  out.cycleId = cycleId;

  // Deliberate deviation from the reference's UB-laden loop -- see file
  // header note. Uses the same well-formed pattern as the TSD/AVDD0/AVDD1
  // loops immediately below (reset shift counter, 16 terms, non-negative
  // shifts), applied to the same byte range (n = 16..31) the reference reads.
  double startAcq = 0.;
  {
    int i = 0;
    for (int n = 16; n < 32; ++n) {
      startAcq += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (30 - 2 * i));
      ++i;
    }
  }
  out.startAcqTimeStamp = startAcq;

  int rawTsd = 0;
  {
    int i = 0;
    for (int n = 32; n < 38; ++n) {
      rawTsd += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (10 - 2 * i));
      ++i;
    }
  }
  out.rawTsd = rawTsd;
  out.temperature = -0.00035207f * static_cast<float>(rawTsd) * static_cast<float>(rawTsd) +
                     2.102526f * static_cast<float>(rawTsd) - 2946.38f;

  int rawAvdd0 = 0;
  {
    int i = 0;
    for (int n = 38; n < 44; ++n) {
      rawAvdd0 += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (10 - 2 * i));
      ++i;
    }
  }
  out.rawAvdd0 = rawAvdd0;
  out.avdd0 = 1.212766f * (static_cast<float>(rawAvdd0) * 3.3f) / 4095.0f;

  int rawAvdd1 = 0;
  {
    int i = 0;
    for (int n = 44; n < 50; ++n) {
      rawAvdd1 += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (10 - 2 * i));
      ++i;
    }
  }
  out.rawAvdd1 = rawAvdd1;
  out.avdd1 = 1.212766f * (static_cast<float>(rawAvdd1) * 3.3f) / 4095.0f;

  int transmitId = 0;
  {
    int i = 0;
    for (int n = 50; n < 54; ++n) {
      transmitId += static_cast<unsigned int>(((frameBytes.at(2 * n + 1) & 0xC0) >> 6) << (6 - 2 * i));
      ++i;
    }
  }
  out.transmitId = transmitId;

  // Per-SCA decode (DecodeRawFrame:347-421): BCID is read BACKWARDS from the
  // end of the payload (one 12-bit Gray value per SCA), while per-channel
  // ADC/hit/gain data is read FORWARDS from `index=0`, accumulating 2*64
  // bytes (high-gain then low-gain) per SCA.
  int index = 0;
  for (int n = 0; n < out.nSca; ++n) {
    const int rawValueBcid = static_cast<int>(frameBytes.at(actualDataSize - 2 * (n + 1) - 2 - 2)) +
                              ((static_cast<int>(frameBytes.at(actualDataSize - 1 - 2 * (n + 1) - 2)) & 0x0F) << 8);
    const int scaAscii = out.nSca - n - 1;
    const int sca = out.nSca - (scaAscii + 1);
    if (sca > kScasInSkiroc - 1) {
      return false;  // matches DecodeRawFrame:356-358
    }

    out.bcid[sca] = convertGrayToBinary(rawValueBcid, 12);

    for (int channel = 0; channel < kChannelsInSkiroc; ++channel) {
      const unsigned short rawData = static_cast<unsigned short>(frameBytes.at(index + 2 * channel)) +
                                      (static_cast<unsigned short>(frameBytes.at(index + 1 + 2 * channel)) << 8);
      const int rawValue = static_cast<int>(rawData & 0xFFF);
      out.hitHigh[sca][channel] = (rawData & 0x1000) >> 12;
      out.gainHigh[sca][channel] = (rawData & 0x2000) >> 13;
      out.adcHigh[sca][channel] = convertGrayToBinary(rawValue, 12);
    }
    index += kChannelsInSkiroc * 2;

    out.nhits[sca] = 0;
    for (int channel = 0; channel < kChannelsInSkiroc; ++channel) {
      const unsigned short rawData = static_cast<unsigned short>(frameBytes.at(index + 2 * channel)) +
                                      (static_cast<unsigned short>(frameBytes.at(index + 1 + 2 * channel)) << 8);
      const int rawValue = static_cast<int>(rawData & 0xFFF);
      out.hitLow[sca][channel] = (rawData & 0x1000) >> 12;
      out.gainLow[sca][channel] = (rawData & 0x2000) >> 13;
      out.adcLow[sca][channel] = convertGrayToBinary(rawValue, 12);
      if (out.hitLow[sca][channel] > 0) {
        ++out.nhits[sca];
      }
    }
    index += kChannelsInSkiroc * 2;
  }
  return true;
}

// ---------------------------------------------------------------------------
// One flushed acquisition cycle, matching the `siwecaldecoded` tree layout
// (SiWECAL-TB-analysis/include/siwecaldecoded.h). All arrays indexed
// [slab][chip][sca] or [slab][chip][sca][channel].
struct Acquisition {
  int acqNumber = -1;
  int nSlboards = kSlbDepth;
  int slot[kSlbDepth];
  int slboardId[kSlbDepth];
  int chipId[kSlbDepth][kSkirocsPerAsu];
  int numCol[kSlbDepth][kSkirocsPerAsu];
  float startAcq[kSlbDepth];
  int rawTsd[kSlbDepth];
  float tsd[kSlbDepth];
  int rawAvdd0[kSlbDepth];
  int rawAvdd1[kSlbDepth];
  float avdd0[kSlbDepth];
  float avdd1[kSlbDepth];

  int bcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  int correctedBcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  int badbcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  int nhits[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  int adcLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  int adcHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  int autogainbitLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  int autogainbitHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  int hitbitLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  int hitbitHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
};

// Port of SLBraw2ROOT::treeInit (SLBraw2ROOT.cc:755-792): reset one
// Acquisition to the "nothing recorded yet" sentinel state before frames for
// a new cycle are recorded into it.
inline void initAcquisition(Acquisition& acq) {
  for (int isl = 0; isl < kSlbDepth; ++isl) {
    for (int k = 0; k < kSkirocsPerAsu; ++k) {
      for (int i = 0; i < kScasInSkiroc; ++i) {
        acq.bcid[isl][k][i] = -999;
        acq.badbcid[isl][k][i] = -999;
        acq.correctedBcid[isl][k][i] = -999;
        acq.nhits[isl][k][i] = -999;
        for (int j = 0; j < kChannelsInSkiroc; ++j) {
          acq.adcLow[isl][k][i][j] = -999;
          acq.adcHigh[isl][k][i][j] = -999;
          acq.autogainbitLow[isl][k][i][j] = -999;
          acq.autogainbitHigh[isl][k][i][j] = -999;
          acq.hitbitLow[isl][k][i][j] = -999;
          acq.hitbitHigh[isl][k][i][j] = -999;
        }
      }
      acq.chipId[isl][k] = -999;
      acq.numCol[isl][k] = 0;
      acq.startAcq[isl] = -1;
      acq.rawTsd[isl] = -1;
      acq.tsd[isl] = -1;
      acq.rawAvdd0[isl] = -1;
      acq.rawAvdd1[isl] = -1;
      acq.avdd0[isl] = -1;
      acq.avdd1[isl] = -1;
    }
    acq.slot[isl] = -1;
    acq.slboardId[isl] = -1;
  }
  acq.nSlboards = -1;
  acq.acqNumber = -1;
}

// Port of SLBraw2ROOT::RecordCycle (SLBraw2ROOT.cc:794-838). Records one
// decoded frame (one slab/chip's data) into the shared Acquisition. Channel
// order is reversed on write (kChannelsInSkiroc-ichn-1), matching the
// reference exactly.
inline void recordFrame(Acquisition& acq, const DecodedFrame& frame, bool zeroSuppression = false) {
  acq.acqNumber = frame.cycleId;
  acq.nSlboards = kSlbDepth;
  int previousBcid = -1000;
  int loopBcid = 0;

  const int slab = frame.slabAdd;
  const int chip = frame.chipId;

  acq.startAcq[slab] = static_cast<float>(frame.startAcqTimeStamp);
  acq.rawTsd[slab] = frame.rawTsd;
  acq.tsd[slab] = frame.temperature;
  acq.rawAvdd0[slab] = frame.rawAvdd0;
  acq.rawAvdd1[slab] = frame.rawAvdd1;
  acq.avdd0[slab] = frame.avdd0;
  acq.avdd1[slab] = frame.avdd1;
  acq.slot[slab] = -1;  // unknown from raw data, as in the reference
  acq.slboardId[slab] = slab;

  if (chip < 0 || chip >= kSkirocsPerAsu) {
    return;  // matches RecordCycle's "Wrong chipId" guard
  }
  acq.chipId[slab][chip] = chip;

  for (int isca = 0; isca < frame.nSca; ++isca) {
    acq.bcid[slab][chip][isca] = frame.bcid[isca];
    acq.nhits[slab][chip][isca] = frame.nhits[isca];
    acq.numCol[slab][chip] = isca + 1;
    if (acq.bcid[slab][chip][isca] > 0 && acq.bcid[slab][chip][isca] - previousBcid < 0) {
      ++loopBcid;
    }
    if (acq.bcid[slab][chip][isca] > 0) {
      acq.correctedBcid[slab][chip][isca] = acq.bcid[slab][chip][isca] + loopBcid * 4096;
    }
    previousBcid = frame.bcid[isca];

    const int nchn = zeroSuppression ? frame.nhits[isca] : kChannelsInSkiroc;
    for (int ichn = 0; ichn < nchn; ++ichn) {
      const int src = kChannelsInSkiroc - ichn - 1;
      acq.adcLow[slab][chip][isca][ichn] = frame.adcLow[isca][src];
      acq.autogainbitLow[slab][chip][isca][ichn] = frame.gainLow[isca][src];
      acq.hitbitLow[slab][chip][isca][ichn] = frame.hitLow[isca][src];
      acq.adcHigh[slab][chip][isca][ichn] = frame.adcHigh[isca][src];
      acq.hitbitHigh[slab][chip][isca][ichn] = frame.hitHigh[isca][src];
      acq.autogainbitHigh[slab][chip][isca][ichn] = frame.gainHigh[isca][src];
    }
  }
}

// Port of SLBraw2ROOT::GetBadBCID (SLBraw2ROOT.cc:840-963), verbatim
// translation of the per-(slab,chip) consecutive-SCA state machine that
// classifies each SCA slot as good(0)/empty-before(1)/empty-after(2)/
// retrigger(3). This is exactly the kind of delicate, already-validated
// logic the implementation plan calls for porting as-is rather than
// re-deriving; do not "simplify" without re-verifying against the reference.
inline void tagBadBcid(Acquisition& acq, int bcidThres = kBcidThresDefault) {
  for (int i = 0; i < kSlbDepth; ++i) {
    for (int k = 0; k < kSkirocsPerAsu; ++k) {
      if (acq.chipId[i][k] < 0) continue;  // only valid chips in this spill
      for (int ibc = 0; ibc < acq.numCol[i][k]; ++ibc) {
        if (ibc == 0) {
          acq.badbcid[i][k][ibc] = 0;
          const int corrBcid = acq.correctedBcid[i][k][ibc];
          int corrBcid1 = 0;
          int corrBcid2 = 0;

          if (acq.correctedBcid[i][k][ibc + 1] > 0 && acq.correctedBcid[i][k][ibc] > 0 &&
              (acq.correctedBcid[i][k][ibc + 1] - acq.correctedBcid[i][k][ibc]) > 0) {
            corrBcid1 = acq.correctedBcid[i][k][ibc + 1];
          }
          if (acq.correctedBcid[i][k][ibc + 2] > 0 &&
              (acq.correctedBcid[i][k][ibc + 2] - acq.correctedBcid[i][k][ibc + 1]) > 0) {
            corrBcid2 = acq.correctedBcid[i][k][ibc + 2];
          }

          if (corrBcid2 > 0) {
            if ((corrBcid2 - corrBcid1) > (bcidThres - 1) && (corrBcid1 - corrBcid) == 1) {
              if (acq.nhits[i][k][ibc] > 0 && acq.nhits[i][k][ibc + 1] == 0) {
                acq.badbcid[i][k][ibc] = 0;
                acq.badbcid[i][k][ibc + 1] = 2;
              }
              if (acq.nhits[i][k][ibc] == 0 && acq.nhits[i][k][ibc + 1] > 0) {
                acq.badbcid[i][k][ibc] = 1;
                acq.badbcid[i][k][ibc + 1] = 0;
              }
            }
            if ((corrBcid2 - corrBcid1) < bcidThres && (corrBcid1 - corrBcid) < bcidThres) {
              acq.badbcid[i][k][ibc] = 3;
              acq.badbcid[i][k][ibc + 1] = 3;
              acq.badbcid[i][k][ibc + 2] = 3;
            }
          }

          if (corrBcid1 > 0 && (corrBcid1 - corrBcid) > 1 && (corrBcid1 - corrBcid) < bcidThres) {
            acq.badbcid[i][k][ibc] = 3;
            acq.badbcid[i][k][ibc + 1] = 3;
          }
        }  // ibc==0

        if (ibc > 0 && acq.badbcid[i][k][ibc] < 0 && acq.correctedBcid[i][k][ibc] > 0 &&
            (acq.correctedBcid[i][k][ibc] - acq.correctedBcid[i][k][ibc - 1]) > 0) {
          acq.badbcid[i][k][ibc] = 0;
          const int corrBcid = acq.correctedBcid[i][k][ibc];
          const int corrBcidMinus = acq.correctedBcid[i][k][ibc - 1];

          if (acq.correctedBcid[i][k][ibc + 1] > 0 &&
              (acq.correctedBcid[i][k][ibc + 1] - acq.correctedBcid[i][k][ibc]) > 0) {
            const int corrBcid1 = acq.correctedBcid[i][k][ibc + 1];

            if (acq.correctedBcid[i][k][ibc + 2] > 0 &&
                (acq.correctedBcid[i][k][ibc + 2] - acq.correctedBcid[i][k][ibc + 1]) > 0) {
              const int corrBcid2 = acq.correctedBcid[i][k][ibc + 2];
              if ((corrBcid2 - corrBcid1) < bcidThres && (corrBcid1 - corrBcid) < bcidThres) {
                acq.badbcid[i][k][ibc] = 3;
                acq.badbcid[i][k][ibc + 1] = 3;
                acq.badbcid[i][k][ibc + 2] = 3;
              }
              if ((corrBcid1 - corrBcid) < bcidThres && (corrBcid - corrBcidMinus) < bcidThres) {
                acq.badbcid[i][k][ibc] = 3;
              }
              if (acq.badbcid[i][k][ibc] != 3 && (corrBcid2 - corrBcid1) > (bcidThres - 1) &&
                  (corrBcid1 - corrBcid) == 1) {
                if (acq.nhits[i][k][ibc] > 0 && acq.nhits[i][k][ibc + 1] == 0) {
                  acq.badbcid[i][k][ibc] = 0;
                  acq.badbcid[i][k][ibc + 1] = 2;
                }
                if (acq.nhits[i][k][ibc] == 0 && acq.nhits[i][k][ibc + 1] > 0) {
                  acq.badbcid[i][k][ibc] = 1;
                  acq.badbcid[i][k][ibc + 1] = 0;
                }
              }
              if (acq.badbcid[i][k][ibc] != 3 && (corrBcid2 - corrBcid1) > (bcidThres - 1) &&
                  (corrBcid1 - corrBcid) > 1 && (corrBcid1 - corrBcid) < bcidThres) {
                acq.badbcid[i][k][ibc] = 3;
                acq.badbcid[i][k][ibc + 1] = 3;
              }
              if ((corrBcid - corrBcidMinus) < bcidThres) {
                acq.badbcid[i][k][ibc] = 3;
              }
            } else {
              if ((corrBcid1 - corrBcid) < bcidThres) acq.badbcid[i][k][ibc] = 3;
              if ((corrBcid - corrBcidMinus) < bcidThres) acq.badbcid[i][k][ibc] = 3;
            }
          } else {
            if ((corrBcid - corrBcidMinus) < bcidThres) acq.badbcid[i][k][ibc] = 3;
          }
        }  // ibc>0
      }  // ibc
    }  // k
  }  // i
}

// ---------------------------------------------------------------------------
// Byte-stream framing: reproduces SLBraw2ROOT::ReadFile's frame-sync state
// machine (SLBraw2ROOT.cc:471-626) up to (but not including) cycle buffering,
// which is CycleAssembler's job below. Yields one raw frame's bytes
// (including its 2-byte core/slab prefix, same layout `decodeCycleId` and
// `decodeFrame` expect) per successful call; frames with a bad trailer are
// silently discarded and scanning continues, matching the reference's
// "WARNING NO TRAILER FOUND" branch.
class RawFrameReader {
 public:
  explicit RawFrameReader(std::istream& in, bool eudaqFormat = false) : m_in(in), m_eudaq(eudaqFormat) {}

  bool nextFrame(std::vector<unsigned char>& frameBytes) {
    if (m_eof) return false;
    for (;;) {
      bool firstframe = false;
      if (m_eudaq) {
        unsigned char c;
        if (!readByte(c)) return false;
        if (c == 0xAB) {
          if (!readByte(c)) return false;
          if (c == 0xCD) {
            firstframe = true;
            m_initFileFound = true;
          }
        }
      } else if (!m_initFileFound) {
        unsigned char c;
        if (!readByte(c)) return false;
        if (c == 0xEE) {
          for (int i = 0; i < 3; ++i) {
            if (!readByte(c)) return false;
            if (c != 0xEE) continue;
            if (i == 2) {
              firstframe = true;
              m_initFileFound = true;
            }
          }
        }
      }

      const bool headerHit = firstframe || (!m_eudaq && m_dataResult == 0xEEEEEEEEu) ||
                              (m_eudaq && (0xFFFFu & m_dataResult) == 0xABCDu);

      if (!headerHit) {
        if (!m_eudaq && m_initFileFound) {
          if (!readRaw(&m_dataResult, sizeof(m_dataResult))) return false;
        }
        continue;
      }

      frameBytes.clear();
      int datasize = 0;
      if (!m_eudaq) {
        int skirocEventNumber = 0;
        if (!readRaw(&skirocEventNumber, sizeof(skirocEventNumber))) return false;
        for (int j = 0; j < 2; ++j) {
          unsigned char c;
          if (!readByte(c)) return false;
          frameBytes.push_back(c);
        }
        unsigned int rawDatasize = 0;
        if (!readRaw(&rawDatasize, sizeof(rawDatasize))) return false;
        datasize = static_cast<int>(rawDatasize);
      } else {
        unsigned char c0, c1, c2;
        if (!readByte(c0) || !readByte(c1) || !readByte(c2)) return false;
        datasize = (static_cast<unsigned short>(c1) << 8) + c0;
        for (int ibuf = 0; ibuf < 5; ++ibuf) {
          unsigned char c;
          if (!readByte(c)) return false;
          if (ibuf == 0) frameBytes.push_back(c);
          if (ibuf == 3) frameBytes.push_back(c);
        }
      }
      for (int j = 0; j < datasize; ++j) {
        unsigned char c;
        if (!readByte(c)) return false;
        frameBytes.push_back(c);
      }
      const unsigned short trailer =
          (static_cast<unsigned short>(frameBytes.at(datasize + 2 - 1)) << 8) + frameBytes.at(datasize + 2 - 2);

      // Unconditional read-ahead for the next iteration's header check
      // (ReadFile:619-621), regardless of whether this frame's trailer is
      // valid.
      if (!m_eudaq && m_initFileFound) {
        if (!readRaw(&m_dataResult, sizeof(m_dataResult))) {
          m_eof = true;
        }
      }

      if (trailer == 0x9697) {
        return true;
      }
      // bad trailer: discard silently and keep scanning (or stop if the
      // read-ahead above already hit EOF).
      if (m_eof) return false;
    }
  }

 private:
  bool readByte(unsigned char& c) { return static_cast<bool>(m_in.read(reinterpret_cast<char*>(&c), 1)); }
  bool readRaw(void* dst, std::size_t n) { return static_cast<bool>(m_in.read(reinterpret_cast<char*>(dst), n)); }

  std::istream& m_in;
  bool m_eudaq;
  bool m_initFileFound = false;
  bool m_eof = false;
  unsigned int m_dataResult = 0;
};

// ---------------------------------------------------------------------------
// Cycle buffering (ReadFile:462-652): frames arrive in arbitrary
// slab/chip interleaving; group them by cycle-id key and only flush (yield)
// the earliest-keyed cycle once more than `maxReadoutCycleJump` distinct
// cycle keys are buffered, so slower/out-of-order slabs' frames for the same
// cycle have a chance to arrive first. `repeatedCycleCount()` mirrors the
// reference's dumped-cycles sort+unique check: if a cycle key reappears as a
// "new" map entry after having already been flushed, the buffer window was
// too small and the caller should retry the whole file with a larger one
// (ConvertDirectorySL_Raw.cc:66-71).
class CycleAssembler {
 public:
  explicit CycleAssembler(int maxReadoutCycleJump = 10, int bcidThres = kBcidThresDefault)
      : m_maxReadoutCycleJump(maxReadoutCycleJump), m_bcidThres(bcidThres) {}

  // Feed one raw frame (as yielded by RawFrameReader::nextFrame). Appends any
  // acquisitions that became ready to flush into `flushed`.
  void addFrame(const std::vector<unsigned char>& frameBytes, std::vector<Acquisition>& flushed) {
    const int key = decodeCycleId(frameBytes);
    auto it = m_cycles.find(key);
    if (it == m_cycles.end()) {
      m_cycles.emplace(key, std::vector<std::vector<unsigned char>>{frameBytes});
      m_seenKeys.push_back(key);
    } else {
      it->second.push_back(frameBytes);
    }

    if (static_cast<int>(m_cycles.size()) > m_maxReadoutCycleJump) {
      flushOldest(flushed);
    }
  }

  // Flush all remaining buffered cycles (call once at EOF).
  void drain(std::vector<Acquisition>& flushed) {
    while (!m_cycles.empty()) {
      flushOldest(flushed);
    }
  }

  // Count of cycle keys that were flushed and later reappeared as a "new"
  // entry -- signals the buffer window was too small (see class comment).
  std::size_t repeatedCycleCount() const {
    std::vector<int> sorted = m_seenKeys;
    std::sort(sorted.begin(), sorted.end());
    const auto last = std::unique(sorted.begin(), sorted.end());
    return m_seenKeys.size() - static_cast<std::size_t>(std::distance(sorted.begin(), last));
  }

 private:
  void flushOldest(std::vector<Acquisition>& flushed) {
    // std::map keeps keys sorted, so begin() is the smallest (earliest) key,
    // matching the reference's "sort into a vector, take element 0" pattern.
    auto it = m_cycles.begin();
    // Acquisition is ~5.5MB (fixed [15][16][15][64] int arrays) -- too large
    // for a stack frame on the smaller worker-thread stacks Gaudi/k4run can
    // use (observed to segfault there despite working fine on a plain
    // process's default stack). Heap-allocate it explicitly.
    auto acq = std::make_unique<Acquisition>();
    bool any = false;
    for (std::size_t jframes = 0; jframes < it->second.size(); ++jframes) {
      if (jframes == 0) {
        initAcquisition(*acq);
      }
      DecodedFrame frame;
      if (decodeFrame(it->second[jframes], frame)) {
        recordFrame(*acq, frame, m_zeroSuppression);
        any = true;
      }
    }
    if (any) {
      tagBadBcid(*acq, m_bcidThres);
      flushed.push_back(*acq);
    }
    m_cycles.erase(it);
  }

  std::map<int, std::vector<std::vector<unsigned char>>> m_cycles;
  std::vector<int> m_seenKeys;
  int m_maxReadoutCycleJump;
  int m_bcidThres;
  bool m_zeroSuppression = false;
};

}  // namespace k4siwecal
