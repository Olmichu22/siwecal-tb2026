import ROOT, sys, os, numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptFit(0)
OUT=sys.argv[1]; LAYER=int(sys.argv[2]); CHIP=int(sys.argv[3])
TH=os.environ.get("TH","th230")
_D={"th230":("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th230/merged_th230.root",
             "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th230/MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt"),
    "th220":("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th220/merged_th220.root",
             "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026/calibration/MuonCalib_gaudi/mips/th220/MIP_pedestalsubmode1_TB2026CERN_run_000th220_highgain.txt")}
MERGED, TXT = _D[TH]
M=np.array([[float(x) for x in l.split()] for l in open(TXT) if not l.startswith('#') and l.strip()])
def rec(l,c,ch):
    m=M[(M[:,0]==l)&(M[:,1]==c)&(M[:,2]==ch)]
    return m[0] if len(m) else None

ROOT.gInterpreter.Declare(r'''
double langaufunC(double *x, double *par){
  double invsq2pi=0.3989422804014, mpshift=-0.22278298; double np=100.0, sc=5.0;
  double summ=0.0,xlow,xupp,step,i,mpc,fland,xx;
  mpc=par[1]-mpshift*par[0]; xlow=x[0]-sc*par[3]; xupp=x[0]+sc*par[3]; step=(xupp-xlow)/np;
  for(i=1.0;i<=np/2;i++){ xx=xlow+(i-0.5)*step; fland=TMath::Landau(xx,mpc,par[0])/par[0]; summ+=fland*TMath::Gaus(x[0],xx,par[3]);
    xx=xupp-(i-0.5)*step; fland=TMath::Landau(xx,mpc,par[0])/par[0]; summ+=fland*TMath::Gaus(x[0],xx,par[3]); }
  return par[2]*step*summ*invsq2pi/par[3]; }''')

def build(f,l,c,ch):
    h=ROOT.TH1F("h_%d_%d_%d"%(l,c,ch),"",900,-100.5,799.5)
    for sca in range(15):
        mh=f.Get("mip_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca)); ph=f.Get("ped_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca))
        if not mh or not ph: continue
        pm=0.
        if ph.GetEntries()>0: ph.GetXaxis().SetRangeUser(ph.GetMean()-20,ph.GetMean()+20); pm=ph.GetMean()
        if pm<=0: continue
        for k in range(900):
            y=mh.GetBinContent(k)
            if y>0: h.Fill(int(mh.GetXaxis().GetBinCenter(k)-pm),y)
    return h

def calibfit(h):
    pllo=[0.1,14.,1.0,0.5]; plhi=[20.,220.,1e8,10.]; fr0,fr1=0.,100.
    h.GetXaxis().SetRangeUser(fr0,fr1)
    fl=ROOT.TF1("fl_%s"%h.GetName(),"landau",fr0,fr1); h.Fit(fl,"QRE")
    sv=[min(max(fl.GetParameter(2),pllo[0]),plhi[0]), min(max(fl.GetParameter(1),pllo[1]),plhi[1]),
        h.Integral("width")*1.2, min(max(h.GetRMS()*0.1,pllo[3]),plhi[3])]
    lg=ROOT.TF1("lg_%s"%h.GetName(),ROOT.langaufunC,fr0,fr1,4); lg.SetParameters(*sv)
    for i in range(4): lg.SetParLimits(i,pllo[i],plhi[i])
    h.Fit(lg,"RBQM"); h.GetXaxis().SetRangeUser(-20,120)
    return lg

def status_col(empv):
    if empv>=0: return ROOT.kGreen+2     # genuine fit
    if empv==-2: return ROOT.kOrange+7   # hist-peak fallback
    if empv==-3: return ROOT.kAzure+1    # chip fallback
    if empv==-6: return ROOT.kCyan+2     # slab fallback
    if empv==-4: return ROOT.kViolet     # global fallback
    return ROOT.kRed                     # masked

f=ROOT.TFile.Open(MERGED)
c=ROOT.TCanvas("c","chip",2600,2600); c.Divide(8,8,0.001,0.001)
for ch in range(64):
    c.cd(ch+1); ROOT.gPad.SetMargin(0.12,0.02,0.10,0.10)
    h=build(f,LAYER,CHIP,ch)
    h.SetTitle("ch %d"%ch); h.GetXaxis().SetRangeUser(-10,90)
    h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1); h.Draw("hist")
    row=rec(LAYER,CHIP,ch)
    mpv=row[3] if row is not None else 0; empv=row[4] if row is not None else -10; nent=row[7] if row is not None else 0
    if h.Integral()>50:
        lg=calibfit(h); lg.SetLineColor(ROOT.kRed); lg.SetLineWidth(2); lg.SetNpx(300); lg.SetRange(0,120); lg.Draw("same"); ROOT.SetOwnership(lg,False)
    ymax=h.GetMaximum()
    if mpv>0:
        ln=ROOT.TLine(mpv,0,mpv,ymax*1.05); ln.SetLineColor(status_col(empv)); ln.SetLineStyle(2); ln.SetLineWidth(2); ln.Draw(); ROOT.SetOwnership(ln,False)
    t=ROOT.TLatex(); t.SetNDC(); t.SetTextSize(0.11); t.SetTextColor(status_col(empv))
    t.DrawLatex(0.45,0.80,"MPV=%.0f"%mpv); ROOT.SetOwnership(t,False)
    t2=ROOT.TLatex(); t2.SetNDC(); t2.SetTextSize(0.09); t2.SetTextColor(ROOT.kGray+2)
    t2.DrawLatex(0.45,0.68,"N=%.0f"%nent); ROOT.SetOwnership(t2,False)
c.cd(0)
tt=ROOT.TLatex(); tt.SetNDC(); tt.SetTextSize(0.013); tt.SetTextColor(ROOT.kBlack)
tt.DrawLatex(0.13,0.992,TH+" MIP high gain -- layer %d chip %d, all 64 channels.  green=genuine fit  orange=hist-max fallback  blue=chip-fb  cyan=slab-fb  violet=global-fb  red=masked"%(LAYER,CHIP))
ROOT.SetOwnership(tt,False)
c.SaveAs(OUT+"/%s_chip_L%d_c%d.png"%(TH,LAYER,CHIP))
print("saved",OUT+"/%s_chip_L%d_c%d.png"%(TH,LAYER,CHIP))
