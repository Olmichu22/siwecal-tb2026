import ROOT, sys, numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptFit(0)
OUT=sys.argv[1]
MERGED="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th230/merged_th230.root"
TXT="/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th230/MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt"

# ---- read new txt ----
rows=[l.split() for l in open(TXT) if not l.startswith('#') and l.strip()]
M=np.array([[float(x) for x in r] for r in rows])
def rec(lay,chip,ch):
    m=M[(M[:,0]==lay)&(M[:,1]==chip)&(M[:,2]==ch)]
    return m[0] if len(m) else None

# langaus TF1 (standard ROOT langaus tutorial form)
ROOT.gInterpreter.Declare(r'''
double langaufun(double *x, double *par){
  double invsq2pi=0.3989422804014, mpshift=-0.22278298;
  double np=100.0, sc=5.0;
  double summ=0.0, xlow, xupp, step, i, mpc, fland, xx;
  mpc=par[1]-mpshift*par[0];
  xlow=x[0]-sc*par[3]; xupp=x[0]+sc*par[3]; step=(xupp-xlow)/np;
  for(i=1.0;i<=np/2;i++){
    xx=xlow+(i-0.5)*step;
    fland=TMath::Landau(xx,mpc,par[0])/par[0];
    summ+=fland*TMath::Gaus(x[0],xx,par[3]);
    xx=xupp-(i-0.5)*step;
    fland=TMath::Landau(xx,mpc,par[0])/par[0];
    summ+=fland*TMath::Gaus(x[0],xx,par[3]);
  }
  return par[2]*step*summ*invsq2pi/par[3];
}''')

def build_hmips(f,lay,chip,ch):
    h=ROOT.TH1F("hm_%d_%d_%d"%(lay,chip,ch),"",900,-100.5,799.5)
    for sca in range(15):
        mh=f.Get("mip_high_s%d_c%d_ch%d_sca%d"%(lay,chip,ch,sca))
        ph=f.Get("ped_high_s%d_c%d_ch%d_sca%d"%(lay,chip,ch,sca))
        if not mh or not ph: continue
        pedMean=0.
        if ph.GetEntries()>0:
            ax=ph.GetXaxis(); ax.SetRangeUser(ph.GetMean()-20,ph.GetMean()+20); pedMean=ph.GetMean()
        if pedMean<=0: continue
        for k in range(900):
            y=mh.GetBinContent(k)
            if y>0: h.Fill(int(mh.GetXaxis().GetBinCenter(k)-pedMean), y)
    return h

FIT_LO, FIT_HI = 0., 100.   # forced fit range (experiment)
def calibfit(h):
    # k4siwecal::fitLangau (highGain) but with a FORCED fit range [FIT_LO,FIT_HI]
    pllo=[0.1,14.,1.0,0.5]; plhi=[20.,220.,1e8,10.]
    fr0, fr1 = FIT_LO, FIT_HI
    h.GetXaxis().SetRangeUser(fr0,fr1)
    fl=ROOT.TF1("fl_%s"%h.GetName(),"landau",fr0,fr1); h.Fit(fl,"QRE")
    sv=[min(max(fl.GetParameter(2),pllo[0]),plhi[0]), min(max(fl.GetParameter(1),pllo[1]),plhi[1]),
        h.Integral("width")*1.2, min(max(h.GetRMS()*0.1,pllo[3]),plhi[3])]
    lg=ROOT.TF1("lg_%s"%h.GetName(),ROOT.langaufun,fr0,fr1,4); lg.SetParameters(*sv)
    for i in range(4): lg.SetParLimits(i,pllo[i],plhi[i])
    h.Fit(lg,"RBQM"); h.GetXaxis().SetRangeUser(-20,120)
    return lg

# channel list: the 9 borderline + 3 high-stat reference
borderline=[(4,2,57),(4,2,59),(4,5,60),(4,5,63),(4,7,44),(4,7,61),(4,10,13),(4,10,31),(13,5,35)]
ref=[(0,4,0),(0,4,10),(7,3,20)]
chans=borderline+ref

f=ROOT.TFile.Open(MERGED)
c=ROOT.TCanvas("c","th230 MIP fits",1800,1200); c.Divide(4,3)
# empv sentinel -> MPV source label (English)
def src_label(empv):
    if empv>=0: return ("fit curve MPV", ROOT.kGreen+2)
    if empv==-2: return ("hist-peak fallback", ROOT.kOrange+7)
    if empv==-3: return ("chip-fit fallback", ROOT.kAzure+1)
    if empv==-6: return ("slab-fit fallback", ROOT.kCyan+2)
    if empv==-4: return ("global-avg fallback", ROOT.kOrange+7)
    return ("masked", ROOT.kRed)
for i,(lay,chip,ch) in enumerate(chans,1):
    c.cd(i)
    h=build_hmips(f,lay,chip,ch)
    h.SetTitle("layer %d  chip %d  ch %d;ADC (pedestal-subtracted);entries"%(lay,chip,ch))
    h.GetXaxis().SetRangeUser(-20,120)
    h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1); h.Draw("hist")
    row=rec(lay,chip,ch)
    mpv=row[3] if row is not None else 0
    empv=row[4] if row is not None else -10
    chi2=row[6] if row is not None else 0
    nent=row[7] if row is not None else 0
    lbl,col=src_label(empv)
    histmax=h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    # Overlay the calibrator's fit (range [0,100]) for visual reference.
    if h.Integral()>50:
        lg=calibfit(h)
        lg.SetLineColor(ROOT.kRed); lg.SetLineWidth(2); lg.SetNpx(400)
        lg.SetRange(0,120); lg.Draw("same")
        ROOT.SetOwnership(lg,False)
    ymax=h.GetMaximum()
    # STORED MPV from the txt (blue dashed): fit-curve MPV if the fit was kept,
    # or the histogram-max if the histmax-consistency filter sent it to fallback
    if mpv>0:
        ln=ROOT.TLine(mpv,0,mpv,ymax*1.05); ln.SetLineColor(ROOT.kBlue+1); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw()
        ROOT.SetOwnership(ln,False)
    # histogram max bin (magenta dotted) -- reference
    lm=ROOT.TLine(histmax,0,histmax,ymax*1.05); lm.SetLineColor(ROOT.kMagenta+1); lm.SetLineStyle(3); lm.SetLineWidth(1); lm.Draw()
    ROOT.SetOwnership(lm,False)
    t=ROOT.TLatex(); t.SetNDC()
    t.SetTextSize(0.055); t.SetTextColor(col)
    t.DrawLatex(0.40,0.83,"MPV=%.1f"%mpv)
    t.SetTextSize(0.042); t.DrawLatex(0.40,0.76,lbl)
    t.SetTextColor(ROOT.kBlack); t.SetTextSize(0.044)
    t.DrawLatex(0.40,0.69,"#chi^{2}/ndf=%.2f"%chi2)
    t.DrawLatex(0.40,0.63,"hist-max=%.0f"%histmax)
    t.DrawLatex(0.40,0.57,"N=%.0f"%nent)
    ROOT.SetOwnership(t,False)
c.cd(0)
tt=ROOT.TLatex(); tt.SetNDC(); tt.SetTextSize(0.018)
tt.DrawLatex(0.20,0.985,"th230 MIP high gain  --  black: histogram   red: calibrator Langau fit   blue dashed: STORED MPV   magenta dotted: histogram-max bin   (fit range 0-100, histmax25%, minN200, slab-fb)")
ROOT.SetOwnership(tt,False)
c.SaveAs(OUT+"/th230_mip_fits.png")
print("saved",OUT+"/th230_mip_fits.png")
