"""Global calibration diagnostics: MPV map, fit-status map, MPV and chi2/ndf
distributions, straight from the DiagnosticsFile written by the calibrator.

  python3 plot_summary.py <outdir> <diagnostics.root> <tag> <mip|ped>
"""
import ROOT, sys, os
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
OUT, ROOTF, TAG, KIND = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
f=ROOT.TFile.Open(ROOTF)
names=[k.GetName() for k in f.GetListOfKeys()]
c=ROOT.TCanvas("c","c",1800,1300); c.Divide(2,2)
for i,n in enumerate(names[:4]):
    c.cd(i+1); ROOT.gPad.SetMargin(0.12,0.14,0.11,0.09)
    h=f.Get(n)
    if h.InheritsFrom("TH2"): h.Draw("colz")
    else: ROOT.gPad.SetLogy(); h.SetLineColor(ROOT.kBlue+1); h.SetLineWidth(2); h.Draw("hist")
    ROOT.SetOwnership(h,False)
os.makedirs(OUT,exist_ok=True)
png="%s/%s_%s_diagnostics.png"%(OUT,TAG,KIND)
c.SaveAs(png); print("saved",png)
