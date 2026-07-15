import ROOT, sys
ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
OUT=sys.argv[1]
MERGED="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch/hist/th230/merged_th230.root"
f=ROOT.TFile.Open(MERGED)

def per_sca(l,c,ch):
    """return list of (sca, TH1F pedestal-subtracted) for this channel."""
    out=[]
    for sca in range(15):
        mh=f.Get("mip_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca)); ph=f.Get("ped_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca))
        if not mh or not ph: continue
        pm=0.
        if ph.GetEntries()>0: ph.GetXaxis().SetRangeUser(ph.GetMean()-20,ph.GetMean()+20); pm=ph.GetMean()
        if pm<=0: continue
        h=ROOT.TH1F("sca_%d_%d_%d_%d"%(l,c,ch,sca),"",140,-20.5,119.5)
        n=0
        for k in range(900):
            y=mh.GetBinContent(k)
            if y>0:
                h.Fill(mh.GetXaxis().GetBinCenter(k)-pm, y); n+=y
        if n>0: out.append((sca,h))
    return out

def combined(l,c,ch,exclude=()):
    h=ROOT.TH1F("comb_%d_%d_%d_%s"%(l,c,ch,"_".join(map(str,exclude))),"",140,-20.5,119.5)
    for sca in range(15):
        if sca in exclude: continue
        mh=f.Get("mip_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca)); ph=f.Get("ped_high_s%d_c%d_ch%d_sca%d"%(l,c,ch,sca))
        if not mh or not ph: continue
        pm=0.
        if ph.GetEntries()>0: ph.GetXaxis().SetRangeUser(ph.GetMean()-20,ph.GetMean()+20); pm=ph.GetMean()
        if pm<=0: continue
        for k in range(900):
            y=mh.GetBinContent(k)
            if y>0: h.Fill(mh.GetXaxis().GetBinCenter(k)-pm, y)
    return h

cols=[ROOT.kBlack,ROOT.kRed,ROOT.kBlue,ROOT.kGreen+2,ROOT.kMagenta,ROOT.kOrange+7,ROOT.kCyan+2,ROOT.kGray+1]
chans=[(4,2,59),(4,2,57),(4,10,31),(0,4,0)]
c=ROOT.TCanvas("c","sca",1800,2000); c.Divide(2,4)
keep=[]
for i,(l,cp,ch) in enumerate(chans):
    # left: per-SCA overlay
    c.cd(2*i+1); ROOT.gPad.SetGrid()
    scas=per_sca(l,cp,ch)
    leg=ROOT.TLegend(0.6,0.5,0.98,0.92); leg.SetTextSize(0.035)
    first=True
    for j,(sca,h) in enumerate(scas[:8]):
        h.SetLineColor(cols[j%len(cols)]); h.SetLineWidth(2 if sca in (0,1) else 1)
        h.SetTitle("layer %d chip %d ch %d -- per-SCA (ped-subtracted);ADC;entries"%(l,cp,ch))
        h.Draw("hist" if first else "hist same"); first=False
        leg.AddEntry(h,"SCA %d%s"%(sca,"  <== first cell" if sca==0 else ""),"l"); keep.append(h)
    leg.Draw(); keep.append(leg)
    # right: combined all vs excl SCA0 vs excl SCA0+1
    c.cd(2*i+2); ROOT.gPad.SetGrid()
    hall=combined(l,cp,ch); h0=combined(l,cp,ch,exclude=(0,)); h01=combined(l,cp,ch,exclude=(0,1))
    hall.SetLineColor(ROOT.kBlack); hall.SetLineWidth(2)
    h0.SetLineColor(ROOT.kRed); h0.SetLineWidth(2)
    h01.SetLineColor(ROOT.kBlue); h01.SetLineWidth(2)
    hall.SetTitle("layer %d chip %d ch %d -- combined;ADC;entries"%(l,cp,ch))
    hall.Draw("hist"); h0.Draw("hist same"); h01.Draw("hist same")
    leg2=ROOT.TLegend(0.55,0.68,0.98,0.92); leg2.SetTextSize(0.035)
    leg2.AddEntry(hall,"all SCAs","l"); leg2.AddEntry(h0,"excluding SCA0","l"); leg2.AddEntry(h01,"excluding SCA0+1","l")
    leg2.Draw(); keep+= [hall,h0,h01,leg2]
c.SaveAs(OUT+"/th230_sca_compare.png")
print("saved",OUT+"/th230_sca_compare.png")
