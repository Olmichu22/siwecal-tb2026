"""Can badbcid / hitbit purge the multi-peak pedestals?

plot_common_mode.py established the mechanism: the second pedestal population is
a beam-induced, CHIP-LOCAL common-mode sag, and it scales with how many hits the
SAME chip took in the same acquisition (|dev|>4 ADC goes 4.65% -> 53.22% from 1
hit on the chip to 50+). That is an activity variable, and the decoded tree
already carries three flags that measure activity: badbcid, hitbit and nhits. So
the question this script answers is whether the calibration can simply REFUSE the
disturbed samples instead of correcting them -- and what that costs.

WHAT THE FLAGS ACTUALLY ARE (measured on run_000060, not assumed):

  * hitbit_low and hitbit_high are THE SAME BIT. 0 disagreements in 623,923,200
    channel-slots. The SKIROC writes one trigger bit per channel and the decoder
    copies it into both gain arrays, so "also require the other gain to be
    quiet" is not a cut at all -- it is a no-op, and this script includes it as
    a variant only to show that it changes literally nothing.
  * nhits == count(hitbit==1) in the slice, exactly, in every filled slice. The
    production MaxNhit=1 is therefore ALREADY a hitbit cut: at most one channel
    of the 64 may be over threshold. The only hitbit tightening left is nhits==0.
  * badbcid is the SLBraw2ROOT state machine (SlbFrameDecoder.h::tagBadBcid):
    0 good, 1 empty-before, 2 empty-after, 3 retrigger. On filled slices it is
    74.2% / 0.5% / 0.7% / 24.7%, i.e. the production badbcid==0 cut already
    throws away a quarter of the sample, essentially all of it retriggers.

The flags say nothing about the OTHER SCAs of the same chip, and that is where
the mechanism lives: the production fill happily accepts a quiet SCA slice from a
chip that took 30 hits in its other slices in the same acquisition. Only 8.6% of
the accepted sample has a chip with <=1 hit total. Hence the variants:

  baseline       production: badbcid==0, bcid>=0, nhits<=MaxNhit, coincidence
  hitbit_bothg   + hitbit_high==0 too (the no-op control described above)
  keep_badbcid   baseline WITHOUT the badbcid==0 cut -- the superset
  only_badbcid   baseline with badbcid!=0 INSTEAD of ==0: the quarter of the
                 sample the production cut already throws away, on its own. If
                 badbcid purges anything, this is where it has to show.
  badbcid_chip   + no OTHER filled SCA of this chip carries badbcid!=0
  nhits0         + nhits==0: no channel of the chip over threshold in the slice
  nhits1         + nhits==1, the complement of nhits0 inside the baseline --
                 the pair isolates what the one allowed hitbit==1 channel does
  chipbusy_le1   + <=1 hit on the whole chip, all its good SCAs, this acquisition
  combo          + badbcid_chip AND chipbusy<=1

TWO METRICS, AND THEY ARE NOT THE SAME QUESTION:

  * PRODUCTION REJECT rate -- TSpectrum(h, 2, "", 0.8) on cells with >50 entries,
    i.e. verbatim what fitPedestalSca does. A cell it calls multi-peak gets no
    fit; the table falls back to the channel average with the -5 sentinel. This
    is the number the calibration actually pays, and it is quoted RAW, at each
    variant's own statistics, because losing entries to a cut is a real cost.
  * DIAGNOSTIC multi-peak fraction -- TSpectrum(h, 2, "nobackground", 0.30) on
    cells with >=MIN_ENTRIES, the setting of plot_pedestal_multipeak.py, so the
    numbers sit on the same scale as the tables. Quoted BOTH raw and with the
    statistics equalised (intersection of cells + multinomial resampling to
    exactly N, the method of compare_multipeak_stat_matched.py) -- a tighter cut
    has fewer entries per cell, and fewer entries hide a shoulder, so a raw
    comparison between variants would credit every cut with a purge it did not
    make.

MEASURED (2026-09-07): 7 th220 runs with beam (60, 61, 62, 86, 87, 88, 254),
1,130 chunks, low gain, 1.534e9 baseline entries. Matched pairwise at N=200:

    variant          cells   multi-peak (base)   reject (base)    RMS (base)   kept
    hitbit_bothg   115,805   21.03% (21.03%)   24.44% (24.38%)   2.79 (2.79)   100%
    keep_badbcid   115,805   21.82% (21.08%)   25.04% (24.44%)   2.80 (2.79)   109%
    only_badbcid    92,407   23.26% (21.21%)   28.00% (24.71%)   2.81 (2.80)   8.8%
    badbcid_chip   111,845   21.83% (21.20%)   23.85% (24.89%)   2.68 (2.81)  64.7%
    nhits0          84,670   20.25% (21.32%)   24.30% (26.53%)   2.74 (2.86)   5.0%
    nhits1         115,497   20.41% (21.06%)   23.41% (24.28%)   2.78 (2.80)  95.0%
    chipbusy_le1    24,582    0.13% (22.54%)    0.90% (27.22%)   1.56 (2.88)   7.6%
    combo           14,976    0.03% (23.00%)    0.55% (28.00%)   0.85 (2.94)   6.1%

  hitbit_bothg is identical to the baseline in every field, as the 0 mismatches
  guarantee it must be. nhits==0 buys 1 pp for 95% of the statistics. badbcid
  does purge something and the rejected quarter is where it shows -- 23.26 vs
  21.21 multi-peak, 28.00 vs 24.71 reject -- but that is 2-3 pp of a 21 pp
  problem, so the cut is worth keeping and is not the mechanism.

  chipbusy<=1 is a different order of thing: the multi-peak fraction goes to
  0.13% and the RMS to 1.56 ADC, BELOW the 1.70 out-of-spill floor. Confirmed
  without any peak finder by the pooled deviation: |dev|>4 ADC goes 14.12% ->
  4.69% -> 0.27% for baseline -> chipbusy_le1 -> combo.

THE CUT SCAN, AND WHY THERE IS NO COMPROMISE WORKING POINT:

    chipbusy <=    1       2       3       5      10    (baseline)
    multi-peak   0.13%  14.71%  13.73%  17.82%  22.36%    ~21%
    reject       0.90%  15.74%  12.96%  16.39%  23.67%    ~24%
    RMS [ADC]     1.56    1.82    2.04    2.33    2.69     2.79
    kept          7.6%   17.1%   27.1%   44.8%   75.9%     100%

  ONE extra hit anywhere on the chip takes the multi-peak fraction from 0.13%
  to 14.71%. The satellite structure appears abruptly between one hit and two,
  while the WIDTH degrades smoothly (1.56 -> 1.82 -> 2.04 -> 2.33) -- two
  different things, and only the first is what the multi-peak count measures.
  The le2/le3/le5 ordering is not monotonic (14.71/13.73/17.82) because each row
  matches a different cell set, which the moving baseline column shows; do not
  read a trend into it. What is unambiguous is the cliff at one hit.

  So the only clean working point is chipbusy<=1, and it costs 92.4% of the
  entries. That is not the real cost. The real cost is CELLS: it leaves 29,118
  cells with >50 entries against the baseline's 118,270, so three quarters of
  the detector would swap the -5 (multi-peak) fallback for the -10 (low
  statistics) one and land on the same channel average either way. The full 89
  th220 runs are only ~1.2x the statistics used here, so pooling everything
  does not rescue it.

  AND THE PEDESTAL MOVES. chipbusy<=1 raises the pedestal mean by +0.20 ADC
  (IQR [-1.33, +1.76]), combo by +0.43 -- the sign the common-mode sag predicts.
  A higher pedestal makes every hit smaller, on a low-gain MIP of ~20 ADC.
  => Purging is not a drop-in fix. Correcting is the route (plot_common_mode.py
  section 10c already gets 46% of the RMS back with the chip almost quiet).

SCOPE: SCA 0-7 only (91% of the accepted sample; nothing above SCA7 ever reaches
the entry counts a matched comparison needs -- see section 3 of the multipeak
tables). One gain, one threshold, whatever runs are given.

Usage:
  RUN=<dir with chunks/>[,<dir>,...] [GAIN=low] [MAXCHUNKS=0] [NS=200,500] [MIN_ENTRIES=200]
  [NSLABS_HIT=8] [MAXNHIT=1] [REF=<out-of-spill grid>] [SEED=0]
  [CACHE=<grids.npz>] \

The grids these numbers come from, and the logs holding the tables above:
  <DEFAULT_CACHE_DIR>/grids_th220_low.npz        the nine flag variants
  <DEFAULT_CACHE_DIR>/grids_th220_low_scan.npz   the chipbusy cut scan
  <DEFAULT_CACHE_DIR>/{purge_low2,pairwise_low,scan_low}.log
  <DEFAULT_CACHE_DIR>/flag_purge_low.png
Re-analysing costs minutes with CACHE set; refilling costs ~20 min.

      python3 diagnostics/plot_flag_purge.py [outdir]
"""
import glob
import os
import re
import sys
import time

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUTDIR, exist_ok=True)

RUN = os.environ["RUN"]
GAIN = os.environ.get("GAIN", "low")
if GAIN not in ("low", "high"):
    raise SystemExit(f"GAIN must be low|high, got {GAIN}")
MAXCHUNKS = int(os.environ.get("MAXCHUNKS", "0"))
NS = [int(x) for x in os.environ.get("NS", "100,200").split(",")]
MIN_ENTRIES = int(os.environ.get("MIN_ENTRIES", "200"))
NSLABS_HIT = int(os.environ.get("NSLABS_HIT", "8"))
MAXNHIT = int(os.environ.get("MAXNHIT", "1"))
SEED = int(os.environ.get("SEED", "0"))
REF = os.environ.get("REF", "")
# Filling 1,130 chunks takes ~20 min; the scan that follows is seconds. CACHE
# stores the filled grids (non-empty cells only) so the analysis can be redone
# without re-reading a terabyte of chunks. Delete the file to force a refill.
#
# Point it at EOS, NOT at a session scratchpad. The 2026-08-06 grids behind
# section 7 of the multipeak tables were written to a scratchpad and went with
# it, which cost a refill on 08-09 to make that section reproducible again. The
# grids this study's numbers come from live in DEFAULT_CACHE_DIR below.
DEFAULT_CACHE_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_flag_purge"
CACHE = os.environ.get("CACHE", "")

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
NSCA_KEEP = 8                      # see SCOPE
NBIN, ADC_LO, ADC_HI = 400, 100.5, 500.5     # the production pedestal booking
NCELL = NSL * NCHIP * NCHN * NSCA_KEEP
SIZE = NCELL * NBIN

# Production fit (fitPedestalSca): >50 entries and TSpectrum(2, "", 0.8)==1 peak.
PROD_MIN_ENTRIES, PROD_SIGMA, PROD_THRESH = 50, 2.0, 0.8
# Diagnostic scan (plot_pedestal_multipeak.py).
DIAG_SIGMA, DIAG_THRESH = 2.0, 0.30

# Any "chipbusy_leN" is understood, so the cut can be SCANNED rather than
# asserted: the interesting question is not whether <=1 hit on the chip is
# clean (it is) but where between 1 and "no cut at all" the cleanliness goes
# away, because that is what decides whether any usable working point exists.
# VARIANTS overrides the list; every name must be one this file defines.
VARIANTS = os.environ.get(
    "VARIANTS",
    "baseline,hitbit_bothg,keep_badbcid,only_badbcid,badbcid_chip,"
    "nhits0,nhits1,chipbusy_le1,combo").split(",")
_CHIPBUSY = re.compile(r"^chipbusy_le(\d+)$")
_KNOWN = {"baseline", "hitbit_bothg", "keep_badbcid", "only_badbcid",
          "badbcid_chip", "nhits0", "nhits1", "combo"}
for _v in VARIANTS:
    if _v not in _KNOWN and not _CHIPBUSY.match(_v):
        raise SystemExit(f"unknown variant {_v!r}")
if "baseline" not in VARIANTS:
    raise SystemExit("VARIANTS must include baseline -- everything is matched to it")

rng = np.random.default_rng(SEED)


# --------------------------------------------------------------------------
# Fill: one pass over the chunks, all variants at once.

class Grid:
    """Per-cell pedestal histograms, filled by flat (cell*NBIN + bin) index.

    Indices are buffered and flushed through np.bincount rather than np.add.at:
    the buffer costs one 8-byte int per entry and the flush is a single pass,
    where np.add.at pays a scatter per element. FLUSH_AT is a memory/speed
    trade, not a correctness knob.
    """
    FLUSH_AT = 8_000_000

    def __init__(self):
        self.h = np.zeros(SIZE, dtype=np.uint32)
        self.buf = []
        self.nbuf = 0
        self.slices = 0
        self.entries = 0

    def add(self, flat):
        self.buf.append(flat)
        self.nbuf += flat.size
        self.entries += flat.size
        if self.nbuf >= self.FLUSH_AT:
            self.flush()

    def flush(self):
        if not self.buf:
            return
        idx = np.concatenate(self.buf)
        self.h += np.bincount(idx, minlength=SIZE).astype(np.uint32)
        self.buf, self.nbuf = [], 0

    def matrix(self):
        """(cells, NBIN) of the non-empty cells only; frees the full grid.

        The full grid is NCELL*NBIN uint32 = ~197 MB per variant, and there are
        nine of them. Handing back only the non-empty rows and dropping the
        dense array keeps the peak at roughly two grids instead of eighteen.
        """
        self.flush()
        m = self.h.reshape(NCELL, NBIN)
        nz = np.flatnonzero(m.sum(axis=1) > 0)
        out = m[nz].copy()
        self.h = None
        return nz, out


cached = CACHE and os.path.exists(CACHE)
grids = {v: Grid() for v in VARIANTS}

bcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
bad = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
nh = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
adc = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.int32)
hitb = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.int32)
hitb_o = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.int32)
nslb = np.zeros(1, dtype=np.int32)

# RUN takes a comma-separated list: the tight variants keep <10% of the sample,
# so one run cannot carry them to the entries-per-cell a matched comparison
# needs. Runs are pooled into the SAME grids, exactly as the production fill
# pools its InputFiles.
runs = [r for r in RUN.split(",") if r]
paths = []
if not cached:
    for r in runs:
        p = sorted(glob.glob(os.path.join(r, "chunks", "chunk_*.root")))
        if not p:
            raise SystemExit(f"no chunks under {r}")
        paths += p[:MAXCHUNKS] if MAXCHUNKS else p
run_name = ",".join(os.path.basename(r.rstrip("/")) for r in runs)
print(f"{run_name}: {len(paths)} chunks, gain={GAIN}, "
      f"MaxNhit={MAXNHIT}, NSlabsHit={NSLABS_HIT}, SCA 0-{NSCA_KEEP - 1}")

# chn index broadcast against the accepted-slice axis, once.
CHN = np.arange(NCHN, dtype=np.int64)
hb_mismatch = 0
t0 = time.time()

for ip, path in enumerate(paths):
    f = ROOT.TFile.Open(path)
    t = f.Get("siwecaldecoded")
    other = "high" if GAIN == "low" else "low"
    t.SetBranchStatus("*", 0)
    for b in ("n_slboards", "bcid", "badbcid", "nhits",
              f"adc_{GAIN}", f"hitbit_{GAIN}", f"hitbit_{other}"):
        t.SetBranchStatus(b, 1)
    t.SetBranchAddress("n_slboards", nslb)
    t.SetBranchAddress("bcid", bcid)
    t.SetBranchAddress("badbcid", bad)
    t.SetBranchAddress("nhits", nh)
    t.SetBranchAddress(f"adc_{GAIN}", adc)
    t.SetBranchAddress(f"hitbit_{GAIN}", hitb)
    t.SetBranchAddress(f"hitbit_{other}", hitb_o)

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        n = int(nslb[0])
        b_, bc, nh_ = bad[:n], bcid[:n], nh[:n]
        hb_mismatch += int((hitb[:n] != hitb_o[:n]).sum())

        filled = bc >= 0
        good = filled & (b_ == 0)
        # Every candidate the LOOSEST variant could want: keep_badbcid drops the
        # badbcid==0 requirement, so the coincidence has to be evaluated for
        # those slices too.
        cand = filled & (nh_ <= MAXNHIT)
        cand[:, :, NSCA_KEEP:] = False
        if not cand.any():
            continue

        # SimpleCoincidenceTagger, vectorised. Per slab, the sorted bcids of its
        # GOOD slices (the tagger's own definition, unchanged by the variant);
        # a candidate is matched by a slab that holds any of them within +-1.
        si, ci, ki = np.nonzero(cand)
        ref = bc[si, ci, ki].astype(np.int64)
        cnt = np.zeros((n, ref.size), dtype=np.int64)
        for s in range(n):
            B = np.sort(bc[s][good[s]])
            if B.size:
                cnt[s] = (np.searchsorted(B, ref + 1, "right")
                          - np.searchsorted(B, ref - 1, "left"))
        # The C++ skips the candidate's own slab entirely, so subtract its own
        # indicator rather than assuming it is 1 -- for a badbcid!=0 candidate
        # its slab may hold no good slice near ref at all.
        seen = (cnt > 0).sum(axis=0) - (cnt[si, np.arange(ref.size)] > 0)
        keep = seen >= (NSLABS_HIT - 1)
        if not keep.any():
            continue
        si, ci, ki = si[keep], ci[keep], ki[keep]

        chip_hits = np.where(good, nh_, 0).sum(axis=2)        # (n, NCHIP)
        chip_bad = (filled & (b_ != 0)).sum(axis=2)           # (n, NCHIP)
        is_good = b_[si, ci, ki] == 0
        nhv = nh_[si, ci, ki]
        cbv = chip_hits[si, ci]
        cbadv = chip_bad[si, ci]

        masks = {
            "baseline": is_good,
            "hitbit_bothg": is_good,
            "keep_badbcid": np.ones(si.size, dtype=bool),
            "only_badbcid": ~is_good,
            "badbcid_chip": is_good & (cbadv == 0),
            "nhits0": is_good & (nhv == 0),
            "nhits1": is_good & (nhv == 1),
            "combo": is_good & (cbadv == 0) & (cbv <= 1),
        }
        for v in VARIANTS:
            mm = _CHIPBUSY.match(v)
            if mm:
                masks[v] = is_good & (cbv <= int(mm.group(1)))

        for v, m in ((v, masks[v]) for v in VARIANTS):
            if not m.any():
                continue
            s_, c_, k_ = si[m], ci[m], ki[m]
            A = adc[s_, c_, k_, :]
            H = hitb[s_, c_, k_, :]
            ok = (H == 0) & (A > ADC_LO) & (A < ADC_HI)
            if v == "hitbit_bothg":
                ok &= hitb_o[s_, c_, k_, :] == 0
            if not ok.any():
                continue
            cell = (((s_.astype(np.int64) * NCHIP + c_) * NCHN)[:, None] + CHN
                    ) * NSCA_KEEP + k_[:, None]
            flat = cell * NBIN + (A - 101).astype(np.int64)
            grids[v].add(flat[ok])
            grids[v].slices += int(m.sum())

    f.Close()
    if (ip + 1) % 10 == 0 or ip + 1 == len(paths):
        print(f"  {ip + 1}/{len(paths)} chunks, {time.time() - t0:.0f}s, "
              f"baseline entries {grids['baseline'].entries:,}")

if not cached:
    print(f"\nhitbit_{GAIN} != hitbit_{other}: {hb_mismatch:,} channel-slots\n")


# --------------------------------------------------------------------------
# Scan.

h1 = ROOT.TH1F("scan", "", NBIN, ADC_LO, ADC_HI)
h1.SetDirectory(0)
sp = ROOT.TSpectrum(10)
buf = np.zeros(NBIN + 2, dtype=np.float64)
CENTRES = ADC_LO + 0.5 + np.arange(NBIN)


def npeaks(counts, sigma, thresh, opt):
    buf[1:NBIN + 1] = counts
    h1.SetContent(buf)
    return sp.Search(h1, sigma, opt, thresh)


def rms(counts):
    tot = counts.sum()
    mu = (counts * CENTRES).sum() / tot
    return float(np.sqrt(max((counts * (CENTRES - mu) ** 2).sum() / tot, 0.0)))


mats, raw, rms_dist, pooled = {}, {}, {}, {}
if cached:
    z = np.load(CACHE)
    for v in VARIANTS:
        mats[v] = (z[f"{v}_nz"], z[f"{v}_m"])
    slices = {v: int(z[f"{v}_slices"]) for v in VARIANTS}
    print(f"loaded grids from {CACHE}")
else:
    for v in VARIANTS:
        mats[v] = grids[v].matrix()
    slices = {v: grids[v].slices for v in VARIANTS}
    if CACHE:
        np.savez_compressed(
            CACHE, **{f"{v}_{k}": val for v in VARIANTS
                      for k, val in (("nz", mats[v][0]), ("m", mats[v][1]),
                                     ("slices", np.int64(slices[v])))})
        print(f"grids cached to {CACHE}")

for v in VARIANTS:
    nz, m = mats[v]
    ent = m.sum(axis=1, dtype=np.int64)
    prod = np.flatnonzero(ent > PROD_MIN_ENTRIES)
    diag = np.flatnonzero(ent >= MIN_ENTRIES)
    nrej = sum(npeaks(m[c], PROD_SIGMA, PROD_THRESH, "") != 1 for c in prod)
    nmul = sum(npeaks(m[c], DIAG_SIGMA, DIAG_THRESH, "nobackground") >= 2
               for c in diag)
    cell_rms = np.array([rms(m[c].astype(np.float64)) for c in diag]) if diag.size \
        else np.array([])
    med_rms = float(np.median(cell_rms)) if cell_rms.size else float("nan")
    rms_dist[v] = cell_rms
    raw[v] = dict(slices=slices[v], entries=int(ent.sum()),
                  cells=len(nz), prod_n=len(prod), prod_rej=nrej,
                  diag_n=len(diag), diag_multi=nmul, rms=med_rms)

# ------------------------------------------------------------------------
# The same comparison WITHOUT TSpectrum. Each cell's histogram is shifted by
# its own mean and the shifted histograms are summed, so the satellite comb
# shows up directly as the tails of one pooled distribution. A peak finder has
# thresholds and a detection power that depend on statistics; this does not,
# and it is the check that the near-zero multi-peak fractions of the tight
# variants are a real narrowing rather than TSpectrum losing sensitivity.
DEV_HALF = 30
DEV_N = 2 * DEV_HALF + 1
COLS = np.arange(NBIN)
for v in VARIANTS:
    nz, m = mats[v]
    ent = m.sum(axis=1, dtype=np.int64)
    sel = np.flatnonzero(ent >= MIN_ENTRIES)
    acc = np.zeros(DEV_N, dtype=np.float64)
    for lo in range(0, sel.size, 8192):
        blk = sel[lo:lo + 8192]
        mm = m[blk].astype(np.float64)
        mu = (mm * CENTRES).sum(axis=1) / mm.sum(axis=1)
        shift = np.rint(mu - (ADC_LO + 0.5)).astype(np.int64)
        tgt = COLS[None, :] - shift[:, None] + DEV_HALF
        ok = (tgt >= 0) & (tgt < DEV_N)
        acc += np.bincount(tgt[ok], weights=mm[ok], minlength=DEV_N)
    pooled[v] = acc

base_e = max(raw["baseline"]["entries"], 1)
print("RAW, each variant at its own statistics")
print(f"{'variant':<14}{'slices':>12}{'entries':>14}{'kept':>7}"
      f"{'cells>50':>10}{'REJECT':>9}{'cells>=' + str(MIN_ENTRIES):>11}"
      f"{'multi':>8}{'RMS':>7}")
for v in VARIANTS:
    r = raw[v]
    pr = 100 * r["prod_rej"] / max(r["prod_n"], 1)
    dm = 100 * r["diag_multi"] / max(r["diag_n"], 1)
    print(f"{v:<14}{r['slices']:>12,}{r['entries']:>14,}"
          f"{100 * r['entries'] / base_e:>6.1f}%{r['prod_n']:>10,}{pr:>8.2f}%"
          f"{r['diag_n']:>11,}{dm:>7.2f}%{r['rms']:>7.2f}")

# --------------------------------------------------------------------------
# Statistics-equalised: identical cells, identical detection power.

def resample(counts, n):
    tot = counts.sum()
    return rng.multinomial(n, counts / tot)


matched = {}
print("\nPOOLED deviation from each cell's own mean "
      f"(cells >= {MIN_ENTRIES} entries, no peak finder)")
print(f"{'variant':<14}{'entries':>15}{'RMS':>8}{'|dev|>4':>10}{'|dev|>2':>10}"
      f"{'cell RMS quartiles':>26}")
DEV_X = np.arange(-DEV_HALF, DEV_HALF + 1, dtype=np.float64)
for v in VARIANTS:
    a = pooled[v]
    tot = a.sum()
    if tot <= 0:
        continue
    mu = (a * DEV_X).sum() / tot
    sd = float(np.sqrt((a * (DEV_X - mu) ** 2).sum() / tot))
    f4 = a[np.abs(DEV_X) > 4].sum() / tot
    f2 = a[np.abs(DEV_X) > 2].sum() / tot
    q = np.percentile(rms_dist[v], [25, 50, 75]) if rms_dist[v].size else [np.nan] * 3
    print(f"{v:<14}{int(tot):>15,}{sd:>8.2f}{100 * f4:>9.2f}%{100 * f2:>9.2f}%"
          f"   {q[0]:>6.2f} {q[1]:>6.2f} {q[2]:>6.2f}")

# PAIRWISE against the baseline, not an intersection of all the variants at
# once. Intersecting them all is dominated by the tightest cut -- it left 1,012
# cells, a fine sample for a 22%-vs-0.3% gap and a poor one for the 2-pp
# differences among the badbcid variants. Each variant is therefore matched to
# the baseline on the cells THOSE TWO share, and the baseline is re-measured on
# each such set, so every row compares two numbers taken on identical cells at
# identical detection power -- the only comparison that means anything here.
def scan_cells(m, pos, cells, N):
    multi, rej, rmss = 0, 0, []
    for c in cells:
        rs = resample(m[pos[c]].astype(np.float64), N).astype(np.float64)
        if npeaks(rs, DIAG_SIGMA, DIAG_THRESH, "nobackground") >= 2:
            multi += 1
        if npeaks(rs, PROD_SIGMA, PROD_THRESH, "") != 1:
            rej += 1
        rmss.append(rms(rs))
    n = max(len(cells), 1)
    return multi / n, rej / n, float(np.median(rmss)) if rmss else float("nan")


print("\nMATCHED pairwise against the baseline "
      "(cells the pair shares, each resampled to exactly N)")
ok_cells = {v: {N: set(mats[v][0][mats[v][1].sum(axis=1, dtype=np.int64) >= N].tolist())
                for N in NS} for v in VARIANTS}
positions = {v: {c: i for i, c in enumerate(mats[v][0].tolist())} for v in VARIANTS}
for N in NS:
    row = {}
    print(f"  N={N}")
    print(f"      {'variant':<14}{'cells':>9}{'multi-peak':>24}"
          f"{'production reject':>26}{'RMS [ADC]':>21}")
    for v in VARIANTS:
        common = sorted(ok_cells[v][N] & ok_cells["baseline"][N])
        if not common:
            print(f"      {v:<14}{0:>9}   no cell clears {N} entries in both")
            continue
        fv = scan_cells(mats[v][1], positions[v], common, N)
        fb = scan_cells(mats["baseline"][1], positions["baseline"], common, N)
        row[v] = (fv, fb, len(common))
        print(f"      {v:<14}{len(common):>9,}"
              f"   {100 * fv[0]:6.2f}% (base {100 * fb[0]:5.2f}%)"
              f"   {100 * fv[1]:6.2f}% (base {100 * fb[1]:5.2f}%)"
              f"   {fv[2]:5.2f} (base {fb[2]:5.2f})")
    matched[N] = row

# --------------------------------------------------------------------------
# What the purge costs in the pedestal VALUE, not just its width.
#
# Every cut here selects quiet chips, and the disturbance it removes is a
# common-mode SAG: plot_common_mode.py measures the median deviation going
# +0.14 ADC at one hit on the chip to -2.44 ADC at fifty. So a purged pedestal
# necessarily sits HIGHER than the production one, and subtracting a higher
# pedestal makes every hit smaller. That shift is a systematic on the energy
# scale and has to be quoted next to any gain in cleanliness. Means are
# unbiased at any statistics, so this is measured on the full samples, on the
# cells the two variants share.
base_nz, base_m = mats["baseline"]
base_ent = base_m.sum(axis=1, dtype=np.int64)
base_mean = (base_m * CENTRES).sum(axis=1) / np.maximum(base_ent, 1)
base_pos = {c: i for i, c in enumerate(base_nz.tolist())}
print(f"\npedestal MEAN shift vs the production baseline "
      f"(cells with >={MIN_ENTRIES} entries in both)")
shifts = {}
for v in VARIANTS:
    nz, m = mats[v]
    ent = m.sum(axis=1, dtype=np.int64)
    mean = (m * CENTRES).sum(axis=1) / np.maximum(ent, 1)
    d = [mean[i] - base_mean[base_pos[c]]
         for i, c in enumerate(nz.tolist())
         if ent[i] >= MIN_ENTRIES and c in base_pos
         and base_ent[base_pos[c]] >= MIN_ENTRIES]
    if not d:
        print(f"  {v:<14}      -- no shared cell")
        continue
    d = np.array(d)
    shifts[v] = float(np.median(d))
    print(f"  {v:<14} {np.median(d):+6.2f} ADC   "
          f"IQR [{np.percentile(d, 25):+.2f}, {np.percentile(d, 75):+.2f}]   "
          f"on {len(d):,} cells")

# The out-of-spill pedestal width of the same cells: the floor any purge is
# trying to reach. Without it a "narrower" variant has nothing to be narrow
# against.
ref_rms = float("nan")
if REF:
    fr_ = ROOT.TFile.Open(REF)
    if fr_ and not fr_.IsZombie():
        cell_re = re.compile(rf"^ped_{GAIN}_s(\d+)_c(\d+)_ch(\d+)_sca(\d+)$")
        vals = []
        for k in fr_.GetListOfKeys():
            mm = cell_re.match(k.GetName())
            if not mm:
                continue
            sl, ch, cn, sc = (int(x) for x in mm.groups())
            if sc >= NSCA_KEEP:
                continue
            h = k.ReadObj()
            if h.GetEntries() >= MIN_ENTRIES:
                vals.append(h.GetRMS())
        fr_.Close()
        if vals:
            ref_rms = float(np.median(vals))
        print(f"\nout-of-spill reference RMS on {len(vals):,} cells: {ref_rms:.2f} ADC")

# --------------------------------------------------------------------------
# Figure: what each cut buys, and what it costs.

c = ROOT.TCanvas("c", "flag purge", 2100, 1250)
c.Divide(3, 2)
keep_alive = []


def bars(pad, title, ylab, values, logy=False, ref=None):
    c.cd(pad)
    ROOT.gPad.SetMargin(0.12, 0.04, 0.22, 0.10)
    if logy:
        ROOT.gPad.SetLogy()
    hb = ROOT.TH1F(f"h{pad}", f"{title};;{ylab}", len(VARIANTS), 0, len(VARIANTS))
    for i, v in enumerate(VARIANTS, start=1):
        hb.SetBinContent(i, values[i - 1])
        hb.GetXaxis().SetBinLabel(i, v)
    hb.GetXaxis().LabelsOption("v")
    hb.SetFillColor(ROOT.kAzure - 9)
    hb.SetBarWidth(0.7)
    hb.SetBarOffset(0.15)
    hb.SetMinimum(0 if not logy else max(min(v for v in values if v > 0) * 0.5, 1e-6))
    hb.SetMaximum(max(values) * 1.25)
    hb.Draw("bar")
    keep_alive.append(hb)
    if ref is not None and np.isfinite(ref):
        ln = ROOT.TLine(0, ref, len(VARIANTS), ref)
        ln.SetLineColor(ROOT.kRed + 1)
        ln.SetLineStyle(2)
        ln.SetLineWidth(2)
        ln.Draw()
        keep_alive.append(ln)
    return hb


N0 = NS[0] if NS[0] in matched else (sorted(matched)[0] if matched else None)
if N0 is not None and matched[N0]:
    row = matched[N0]
    def col(i):
        return [row[v][0][i] if v in row else 0.0 for v in VARIANTS]
    bars(1, f"multi-peak fraction, matched to the baseline at N={N0} "
            f"({GAIN} gain)", "multi-peak fraction", col(0))
    bars(3, f"cells the production fit REFUSES (-5 sentinel), matched at "
            f"N={N0}", "reject fraction", col(1))
    # The floor goes in the axis title, not the pad title: a second line of
    # text under the pad collides with the rotated variant labels.
    bars(4, f"pedestal RMS at N={N0}",
         "median RMS [ADC]  (red dashes = out-of-spill floor)", col(2),
         ref=ref_rms)
bars(2, "statistics kept, relative to the production baseline",
     "entries / baseline entries",
     [raw[v]["entries"] / base_e for v in VARIANTS], logy=True)

# Panels 5-6: the peak-finder-free view. Overlaid curves, one per variant.
c.cd(5)
ROOT.gPad.SetMargin(0.12, 0.04, 0.12, 0.10)
ROOT.gPad.SetLogy()
COLORS = [ROOT.kBlack, ROOT.kGray + 1, ROOT.kGray + 2, ROOT.kRed + 1,
          ROOT.kOrange + 7, ROOT.kSpring - 6, ROOT.kGreen + 3, ROOT.kAzure + 1,
          ROOT.kViolet + 1]
leg = ROOT.TLegend(0.60, 0.55, 0.96, 0.88)
leg.SetBorderSize(0)
leg.SetTextSize(0.030)
first = True
for i, v in enumerate(VARIANTS):
    a = pooled[v]
    if a.sum() <= 0:
        continue
    hp = ROOT.TH1F(f"dev{i}", "pedestal ADC minus that cell's own mean;"
                              "deviation [ADC];fraction of entries",
                   DEV_N, -DEV_HALF - 0.5, DEV_HALF + 0.5)
    for b in range(DEV_N):
        hp.SetBinContent(b + 1, a[b] / a.sum())
    hp.SetLineColor(COLORS[i % len(COLORS)])
    hp.SetLineWidth(2)
    hp.SetMinimum(1e-6)
    hp.SetMaximum(1.0)
    hp.Draw("hist" if first else "hist same")
    first = False
    leg.AddEntry(hp, v, "l")
    keep_alive.append(hp)
leg.Draw()
keep_alive.append(leg)

c.cd(6)
ROOT.gPad.SetMargin(0.12, 0.04, 0.12, 0.10)
leg2 = ROOT.TLegend(0.60, 0.55, 0.96, 0.88)
leg2.SetBorderSize(0)
leg2.SetTextSize(0.030)
first = True
for i, v in enumerate(VARIANTS):
    d = rms_dist[v]
    if d.size == 0:
        continue
    hr = ROOT.TH1F(f"rms{i}", f"per-cell pedestal RMS (cells >= {MIN_ENTRIES} "
                              "entries);RMS [ADC];fraction of cells", 60, 0, 6)
    for x in d:
        hr.Fill(x)
    if hr.Integral() > 0:
        hr.Scale(1.0 / hr.Integral())
    hr.SetLineColor(COLORS[i % len(COLORS)])
    hr.SetLineWidth(2)
    hr.SetMaximum(0.35)
    hr.Draw("hist" if first else "hist same")
    first = False
    leg2.AddEntry(hr, v, "l")
    keep_alive.append(hr)
if np.isfinite(ref_rms):
    ln2 = ROOT.TLine(ref_rms, 0, ref_rms, 0.35)
    ln2.SetLineColor(ROOT.kRed + 1)
    ln2.SetLineStyle(2)
    ln2.SetLineWidth(2)
    ln2.Draw()
    keep_alive.append(ln2)
leg2.Draw()
keep_alive.append(leg2)

out = os.path.join(OUTDIR, f"flag_purge_{GAIN}.png")
c.SaveAs(out)
print(f"\nsaved {out}")
