"""Stage 1 of the trains_out_of_spill pipeline: pull selected acquisitions out
of a run's decoded chunks into one small siwecaldecoded file.

The three acquisitions of interest on run_000060 sit ~13,000 apart in a 28,860
acquisition run spread over 59 chunks. Event-building the whole run to look at
three of them would be absurd, so this copies just those entries -- the output
is a normal siwecaldecoded tree that EcalEventBuilder reads unchanged.

TWO THINGS THIS MUST GET RIGHT, both about chunk boundaries:

  * An acquisition split across a chunk boundary is written TWICE, as two
    COMPLEMENTARY halves sharing an acqNumber. Both halves must be copied or
    the event loses half its hits.
  * EcalEventBuilder detects that split by comparing acqNumber with the
    PREVIOUS entry (EcalEventBuilder.cpp:234) and merging. So the two halves
    must land in the output as CONSECUTIVE entries. Copying chunk by chunk in
    entry order preserves that; reordering or deduplicating would break it.

Note the `spill` branch the viewer shows is the builder's ENTRY INDEX in its
input, not acqNumber -- so in the output the acquisitions come out as spill
0, 1, 2 in the order listed here. The mapping is printed at the end.

Usage:
  RUN=<dir with chunks/> ACQS=2631,16038,26053 OUT=<file.root> \
      python3 gaudi_jobs/trains_out_of_spill/extract_acquisitions.py
"""
import glob
import os

import numpy as np
import ROOT

RUN = os.environ["RUN"]
OUT = os.environ["OUT"]
ACQS = [int(x) for x in os.environ["ACQS"].split(",")]
WANTED = set(ACQS)

paths = sorted(glob.glob(os.path.join(RUN, "chunks", "chunk_*.root")))
if not paths:
    raise SystemExit(f"no chunks under {RUN}")
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)

# Everything goes through a TChain, and that is not a stylistic choice.
# Cloning per file does NOT work: TTree::CloneTree(0) binds the output tree's
# branch addresses to the SOURCE tree's buffers, so as soon as the next chunk is
# opened those addresses no longer point at the tree being read and Fill() silently
# writes stale/zero entries. The first version of this script did exactly that and
# produced a file whose 2nd and 3rd entries were empty (acqNumber=0, no hits) --
# which then looked like a chunk-boundary splice to the event builder, because two
# consecutive entries shared acqNumber 0. TChain::LoadTree re-copies addresses to
# registered clones on every file switch, so cloning from the chain is correct.

# Pass 1: a separate acqNumber-only chain to locate the entries. Reading full
# ~4 MB entries for all 28,860 acquisitions just to find three would dominate
# the runtime, and branch status cannot be narrowed on the copy chain (the clone
# only carries branches that are active when it is made).
scan = ROOT.TChain("siwecaldecoded")
for path in paths:
    scan.Add(path)
scan.SetBranchStatus("*", 0)
scan.SetBranchStatus("acqNumber", 1)
acqn = np.zeros(1, dtype=np.int32)
scan.SetBranchAddress("acqNumber", acqn)
found = []                                   # (global entry index, acqNumber)
for i in range(scan.GetEntries()):
    scan.GetEntry(i)
    if int(acqn[0]) in WANTED:
        found.append((i, int(acqn[0])))
print(f"scanned {scan.GetEntries():,} entries over {len(paths)} chunks, "
      f"{len(found)} match")
if not found:
    raise SystemExit(f"none of {ACQS} found in {RUN}")

# Pass 2: copy those entries, in chain order so that the two halves of an
# acquisition split across a chunk boundary stay CONSECUTIVE.
chain = ROOT.TChain("siwecaldecoded")
for path in paths:
    chain.Add(path)
out_file = ROOT.TFile.Open(OUT, "RECREATE")
out_tree = chain.CloneTree(0)
copied = {}          # acqNumber -> how many entries (2 = split across a boundary)
order = []           # acqNumbers in output order, for the spill-index mapping
for i, a in found:
    chain.GetEntry(i)
    out_tree.Fill()
    copied[a] = copied.get(a, 0) + 1
    if a not in order:
        order.append(a)
    print(f"  entry {i}: acq {a}")

out_file.cd()
out_tree.Write()
n = out_tree.GetEntries()
out_file.Close()

print(f"\nwrote {OUT}: {n} entries for {len(copied)} acquisitions")
missing = WANTED - set(copied)
if missing:
    print(f"  WARNING: not found: {sorted(missing)}")
for a in order:
    half = "  (split across a chunk boundary, both halves copied)" if copied[a] == 2 else ""
    print(f"  acq {a}: {copied[a]} entry/entries{half}")
# The builder sets `spill` to the ENTRY INDEX of an acquisition in its input
# (EcalEventBuilder.cpp:271), and when it splices two halves it keeps the index
# of the FIRST half. So the mapping is a running sum over copied entry counts,
# not 0,1,2 -- it only looks sequential when nothing was split.
print("\nspill index in the viewer -> acquisition:")
idx = 0
for a in order:
    print(f"  spill {idx} -> acq {a}")
    idx += copied[a]
