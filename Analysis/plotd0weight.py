import ROOT as rt

fileelecms = rt.TFile("./Analysis/HEPData-ins1317640-v1-Table_5.root", "READ")
histelecms = fileelecms.Get("Table 5/Hist1D_y1")
filemuocms = rt.TFile("./Analysis/HEPData-ins1317640-v1-Table_6.root", "READ")
histmuocms = filemuocms.Get("Table 6/Hist1D_y1")
filefre = rt.TFile("./Analysis/d0weightFreya.root", "READ")
histelefre = filefre.Get("d0EffEleFreya")
histmuofre = filefre.Get("d0EffMuoFreya")
filenis = rt.TFile("./Analysis/d0weightNishita.root", "READ")
histelenis = filenis.Get("d0EffEleNishita")
histmuonis = filenis.Get("d0EffMuoNishita")
filekam = rt.TFile("./Analysis/d0weightKamal.root", "READ")
histelekam = filekam.Get("d0EffEleKamal")
histmuokam = filekam.Get("d0EffMuoKamal")

def plotCosmetics(histlist,
                  labellist,
                  xaxistitle,
                  yaxistitle,
                  outplotname,
                  logy):

    c1 = rt.TCanvas("","",10,32,782,552);

    ctr = 0
    for hist in histlist:
        hist.GetXaxis().SetTitle(xaxistitle)
        hist.GetYaxis().SetTitle(yaxistitle)
        if ctr<4:
            hist.SetLineColorAlpha(ctr+1,0.5)
        if ctr>=4:
            hist.SetLineColorAlpha(ctr+2,0.5)
        ctr = ctr+1
        hist.SetLineWidth(2)

    c1.SetLogy(logy)
    rt.gStyle.SetOptStat(0)

    for hist in histlist:
        hist.Draw("hist same")

    legc1 = rt.TLegend(0.5, 0.9, 0.89, 1.0, "", "brNDC")
    ctr = 0
    for hist in histlist:
        legc1.AddEntry(hist, labellist[ctr], "l")
        ctr = ctr+1
    legc1.SetTextSize(0.03)
    legc1.SetBorderSize(0)
    legc1.Draw()

    c1.SaveAs(outplotname+".pdf")
    return

histlist = [histelekam, histelefre, histelenis, histelecms]
labellist = ["Kamal weight scheme", "Freya weight scheme", "Nishita weight scheme", "CMS HEP Database"]
plotCosmetics(histlist, labellist, "d0ele", "efficiency", "weightele", False)

histlist = [histmuokam, histmuofre, histmuonis, histmuocms]
labellist = ["Kamal weight scheme", "Freya weight scheme", "Nishita weight scheme", "CMS HEP Database"]
plotCosmetics(histlist, labellist, "d0muo", "efficiency", "weightmuo", False)
