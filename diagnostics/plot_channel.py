"""Zoom on a single channel: pedestal-subtracted MIP spectrum at several
rebinnings + the per-SCA decomposition, to judge whether an apparent
double peak survives rebinning/statistics.

  TH=th220 python3 plot_channel.py <outdir> <layer> <chip> <channel>
"""
import ROOT, sys, os, numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptFit(0)
OUT=sys.argv[1]; LAYER=int(sys.argv[2]); CHIP=int(sys.argv[3]); CH=int(sys.argv[4])
TH=os.environ.get("TH","th220")
_D={"th230":("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th230/merged_th230.root",
             "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th230/MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt"),
    "th220":("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th220/merged_th220.root",
             "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th220/MIP_pedestalsubmode1_TB2026CERN_run_000th220_highgain.txt")}
MERGED, TXT = _D[TH]
M=np.array([[float(x) for x in l.split()] for l in open(TXT) if not l.startswith('#') and l.strip()])
m=M[(M[:,0]==LAYER)&(M[:,1]==CHIP)&(M[:,2]==CH)]
row=m[0] if len(m) else None

ROOT.gInterpreter.Declare(r'''
double langaufunC(double *x, double *par){
  double invsq2pi=0.3989422804014, mpshift=-0.22278298; double np=100.0, sc=5.0;
  double summ=0.0,xlow,xupp,step,i,mpc,fland,xx;
  mpc=par[1]-mpshift*par[0]; xlow=x[0]-sc*par[3]; xupp=x[0]+sc*par[3]; step=(xupp-xlow)/np;
  for(i=1.0;i<=np/2;i++){ xx=xlow+(i-0.5)*step; fland=TMath::Landau(xx,mpc,par[0])/par[0]; summ+=fland*TMath::Gaus(x[0],xx,par[3]);
    xx=xupp-(i-0.5)*step; fland=TMath::Landau(xx,mpc,par[0])/par[0]; summ+=fland*TMath::Gaus(x[0],xx,par[3]); }
  return par[2]*step*summ*invsq2pi/par[3]; }''')

f=ROOT.TFile.Open(MERGED)

def build(name, scas):
    h=ROOT.TH1F(name,"",900,-100.5,799.5)
    for sca in scas:
        mh=f.Get("mip_high_s%d_c%d_ch%d_sca%d"%(LAYER,CHIP,CH,sca))
        ph=f.Get("ped_high_s%d_c%d_ch%d_sca%d"%(LAYER,CHIP,CH,sca))
        if not mh or not ph: continue
        if ph.GetEntries()<=0: continue
        ph.GetXaxis().SetRangeUser(ph.GetMean()-20,ph.GetMean()+20); pm=ph.GetMean()
        if pm<=0: continue
        for k in range(900):
            y=mh.GetBinContent(k)
            if y>0: h.Fill(int(mh.GetXaxis().GetBinCenter(k)-pm),y)
    return h

def calibfit(h,tag):
    pllo=[0.1,14.,1.0,0.5]; plhi=[20.,220.,1e8,10.]; fr0,fr1=0.,100.
    h.GetXaxis().SetRangeUser(fr0,fr1)
    fl=ROOT.TF1("fl_"+tag,"landau",fr0,fr1); h.Fit(fl,"QRE")
    sv=[min(max(fl.GetParameter(2),pllo[0]),plhi[0]), min(max(fl.GetParameter(1),pllo[1]),plhi[1]),
        h.Integral("width")*1.2, min(max(h.GetRMS()*0.1,pllo[3]),plhi[3])]
    lg=ROOT.TF1("lg_"+tag,ROOT.langaufunC,fr0,fr1,4); lg.SetParameters(*sv)
    for i in range(4): lg.SetParLimits(i,pllo[i],plhi[i])
    r=h.Fit(lg,"RBQMS"); h.GetXaxis().SetRangeUser(-10,90)
    chi2ndf=lg.GetChisquare()/lg.GetNDF() if lg.GetNDF()>0 else -1
    return lg, chi2ndf

c=ROOT.TCanvas("c","ch",2000,1300); c.Divide(2,2)
allsca=list(range(15))

# --- pads 1-3: same spectrum at 1 / 2 / 4 ADC binning ---
for i,rb in enumerate([1,2,4]):
    c.cd(i+1); ROOT.gPad.SetMargin(0.10,0.03,0.11,0.09)
    h=build("h_rb%d"%rb, allsca)
    if rb>1: h.Rebin(rb)
    h.SetTitle("%s  L%d chip%d ch%d -- %d ADC/bin;ADC - pedestal;entries"%(TH,LAYER,CHIP,CH,rb))
    h.GetXaxis().SetRangeUser(-10,90); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1); h.Draw("hist e")
    lg,chi2=calibfit(h,"rb%d"%rb)
    lg.SetLineColor(ROOT.kRed); lg.SetLineWidth(2); lg.SetNpx(400); lg.SetRange(0,100); lg.Draw("same")
    ROOT.SetOwnership(h,False); ROOT.SetOwnership(lg,False)
    t=ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.038)
    t.DrawLatex(0.55,0.84,"N = %.0f"%h.Integral())
    t.DrawLatex(0.55,0.78,"curve MPV = %.2f"%lg.GetParameter(1))
    t.DrawLatex(0.55,0.72,"#chi^{2}/ndf = %.2f"%chi2)
    if row is not None:
        t.SetTextColor(ROOT.kBlue+1)
        t.DrawLatex(0.55,0.66,"stored MPV = %.2f (status %s)"%(row[3], "fit" if row[4]>=0 else "fb %.0f"%row[4]))
    ROOT.SetOwnership(t,False)

# --- pad 4: per-SCA decomposition ---
c.cd(4); ROOT.gPad.SetMargin(0.10,0.03,0.11,0.09)
cols=[ROOT.kRed+1,ROOT.kBlue+1,ROOT.kGreen+2,ROOT.kMagenta+1,ROOT.kOrange+7,ROOT.kCyan+2,ROOT.kGray+2]
leg=ROOT.TLegend(0.60,0.45,0.96,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
hsca0=build("h_sca0",[0]); hsca0.Rebin(2)
hrest=build("h_scarest",list(range(1,15))); hrest.Rebin(2)
hsca0.SetTitle("%s  L%d chip%d ch%d -- SCA0 vs SCA1-14 (2 ADC/bin, area-normalised);ADC - pedestal;a.u."%(TH,LAYER,CHIP,CH))
for h,cl,lab in [(hsca0,ROOT.kRed+1,"SCA0  (N=%.0f)"%hsca0.Integral()),(hrest,ROOT.kBlue+1,"SCA1-14  (N=%.0f)"%hrest.Integral())]:
    if h.Integral()>0: h.Scale(1.0/h.Integral())
    h.SetLineColor(cl); h.SetLineWidth(2); leg.AddEntry(h,lab,"l"); ROOT.SetOwnership(h,False)
hsca0.GetXaxis().SetRangeUser(-10,90)
hsca0.SetMaximum(1.25*max(hsca0.GetMaximum(),hrest.GetMaximum()))
hsca0.Draw("hist"); hrest.Draw("hist same"); leg.Draw(); ROOT.SetOwnership(leg,False)

os.makedirs(OUT,exist_ok=True)
png=OUT+"/%s_L%dc%dch%d_channel.png"%(TH,LAYER,CHIP,CH)
c.SaveAs(png); print("saved",png)
