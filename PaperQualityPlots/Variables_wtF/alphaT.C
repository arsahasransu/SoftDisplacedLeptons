void alphaT()
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Jun  4 12:25:37 2020) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,782,600);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: 
   TPad *pad = new TPad("", "",0.05,0.05,1,0.95);
   TPad *pad2 = new TPad("", "",0.001,0.08,0.05,0.85);
   TPad *pad3 = new TPad("", "",0.08,0.001,1,0.075);
   TPad *pad4 = new TPad("", "",0.001,0.88,1,0.98);

   pad->Draw();
   pad->cd();
   pad->Range(0,0,1,1);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetLogy();
   pad->SetFrameBorderMode(0);
   
   TH1F *alphaT0__21 = new TH1F("alphaT0__21","",25,0,1.25);
   alphaT0__21->SetBinContent(3,0.0009545021);
   alphaT0__21->SetBinContent(4,0.0009545021);
   alphaT0__21->SetBinContent(5,0.006999682);
   alphaT0__21->SetBinContent(6,0.01972638);
   alphaT0__21->SetBinContent(7,0.05217944);
   alphaT0__21->SetBinContent(8,0.1244034);
   alphaT0__21->SetBinContent(9,0.2577156);
   alphaT0__21->SetBinContent(10,0.4129812);
   alphaT0__21->SetBinContent(11,0.08685969);
   alphaT0__21->SetBinContent(12,0.009545021);
   alphaT0__21->SetBinContent(13,0.007317849);
   alphaT0__21->SetBinContent(14,0.004454343);
   alphaT0__21->SetBinContent(15,0.001909004);
   alphaT0__21->SetBinContent(16,0.002227172);
   alphaT0__21->SetBinContent(17,0.001272669);
   alphaT0__21->SetBinContent(18,0.002227172);
   alphaT0__21->SetBinContent(19,0.0006363347);
   alphaT0__21->SetBinContent(20,0.001590837);
   alphaT0__21->SetBinContent(21,0.001272669);
   alphaT0__21->SetBinContent(22,0.0009545021);
   alphaT0__21->SetBinContent(23,0.0009545021);
   alphaT0__21->SetBinContent(24,0.001272669);
   alphaT0__21->SetBinContent(25,0.001590837);
   alphaT0__21->SetBinContent(26,0.3070315);
   alphaT0__21->SetBinError(3,0.000551082);
   alphaT0__21->SetBinError(4,0.000551082);
   alphaT0__21->SetBinError(5,0.001492337);
   alphaT0__21->SetBinError(6,0.002505252);
   alphaT0__21->SetBinError(7,0.00407453);
   alphaT0__21->SetBinError(8,0.006291352);
   alphaT0__21->SetBinError(9,0.009055202);
   alphaT0__21->SetBinError(10,0.01146286);
   alphaT0__21->SetBinError(11,0.005256987);
   alphaT0__21->SetBinError(12,0.001742674);
   alphaT0__21->SetBinError(13,0.001525877);
   alphaT0__21->SetBinError(14,0.001190473);
   alphaT0__21->SetBinError(15,0.0007793477);
   alphaT0__21->SetBinError(16,0.0008417917);
   alphaT0__21->SetBinError(17,0.0006363347);
   alphaT0__21->SetBinError(18,0.0008417917);
   alphaT0__21->SetBinError(19,0.0004499566);
   alphaT0__21->SetBinError(20,0.0007114438);
   alphaT0__21->SetBinError(21,0.0006363347);
   alphaT0__21->SetBinError(22,0.000551082);
   alphaT0__21->SetBinError(23,0.000551082);
   alphaT0__21->SetBinError(24,0.0006363347);
   alphaT0__21->SetBinError(25,0.0007114438);
   alphaT0__21->SetBinError(26,0.009883694);
   alphaT0__21->SetMinimum(0.0001);
   alphaT0__21->SetMaximum(0.5);
   alphaT0__21->SetEntries(4108);
   alphaT0__21->SetLineColor(1);
   alphaT0__21->SetLineWidth(3);
   alphaT0__21->GetXaxis()->SetRange(1,26);
   alphaT0__21->GetXaxis()->SetTicks("-");
   alphaT0__21->GetXaxis()->SetLabelFont(42);
   alphaT0__21->GetXaxis()->SetLabelOffset(-0.08);
   alphaT0__21->GetXaxis()->SetLabelSize(0.055);
   alphaT0__21->GetXaxis()->SetTitleSize(0.035);
   alphaT0__21->GetXaxis()->SetTitleOffset(1);
   alphaT0__21->GetXaxis()->SetTitleFont(42);
   alphaT0__21->GetYaxis()->SetTicks("+");
   alphaT0__21->GetYaxis()->SetLabelFont(42);
   alphaT0__21->GetYaxis()->SetLabelOffset(-0.04);
   alphaT0__21->GetYaxis()->SetLabelSize(0.055);
   alphaT0__21->GetYaxis()->SetTitleSize(0.035);
   alphaT0__21->GetYaxis()->SetTitleFont(42);
   alphaT0__21->GetZaxis()->SetLabelFont(42);
   alphaT0__21->GetZaxis()->SetLabelSize(0.035);
   alphaT0__21->GetZaxis()->SetTitleSize(0.035);
   alphaT0__21->GetZaxis()->SetTitleOffset(1);
   alphaT0__21->GetZaxis()->SetTitleFont(42);
   alphaT0__21->Draw("hist  E");
   
   TH1F *alphaT1__22 = new TH1F("alphaT1__22","",25,0,1.25);
   alphaT1__22->SetBinContent(3,6.703743e-05);
   alphaT1__22->SetBinContent(4,0.0002414081);
   alphaT1__22->SetBinContent(5,0.001889003);
   alphaT1__22->SetBinContent(6,0.01328201);
   alphaT1__22->SetBinContent(7,0.05272788);
   alphaT1__22->SetBinContent(8,0.1330014);
   alphaT1__22->SetBinContent(9,0.270748);
   alphaT1__22->SetBinContent(10,0.3901241);
   alphaT1__22->SetBinContent(11,0.1046047);
   alphaT1__22->SetBinContent(12,0.0161931);
   alphaT1__22->SetBinContent(13,0.005490457);
   alphaT1__22->SetBinContent(14,0.002252116);
   alphaT1__22->SetBinContent(15,0.001773848);
   alphaT1__22->SetBinContent(16,0.001170696);
   alphaT1__22->SetBinContent(17,0.001613413);
   alphaT1__22->SetBinContent(18,0.0008005109);
   alphaT1__22->SetBinContent(19,0.0008355223);
   alphaT1__22->SetBinContent(20,0.0006353742);
   alphaT1__22->SetBinContent(21,0.0005448821);
   alphaT1__22->SetBinContent(22,0.0006053544);
   alphaT1__22->SetBinContent(23,0.0002537007);
   alphaT1__22->SetBinContent(24,0.0005559261);
   alphaT1__22->SetBinContent(25,0.0005894126);
   alphaT1__22->SetBinContent(26,0.008697438);
   alphaT1__22->SetBinError(3,6.703743e-05);
   alphaT1__22->SetBinError(4,0.0001212132);
   alphaT1__22->SetBinError(5,0.0003122881);
   alphaT1__22->SetBinError(6,0.0008448906);
   alphaT1__22->SetBinError(7,0.001693289);
   alphaT1__22->SetBinError(8,0.002696148);
   alphaT1__22->SetBinError(9,0.003860315);
   alphaT1__22->SetBinError(10,0.00464071);
   alphaT1__22->SetBinError(11,0.00236891);
   alphaT1__22->SetBinError(12,0.0009232786);
   alphaT1__22->SetBinError(13,0.0005324385);
   alphaT1__22->SetBinError(14,0.0003573298);
   alphaT1__22->SetBinError(15,0.0003139775);
   alphaT1__22->SetBinError(16,0.0002524756);
   alphaT1__22->SetBinError(17,0.0003019916);
   alphaT1__22->SetBinError(18,0.0002047313);
   alphaT1__22->SetBinError(19,0.0002099951);
   alphaT1__22->SetBinError(20,0.0001925149);
   alphaT1__22->SetBinError(21,0.0001736828);
   alphaT1__22->SetBinError(22,0.0001849424);
   alphaT1__22->SetBinError(23,0.0001209053);
   alphaT1__22->SetBinError(24,0.0001789274);
   alphaT1__22->SetBinError(25,0.0001797627);
   alphaT1__22->SetBinError(26,0.0006960215);
   alphaT1__22->SetMinimum(0.0005);
   alphaT1__22->SetMaximum(0.2);
   alphaT1__22->SetEntries(19133);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   alphaT1__22->SetLineColor(ci);
   alphaT1__22->SetLineWidth(3);
   alphaT1__22->GetXaxis()->SetRange(1,26);
   alphaT1__22->GetXaxis()->SetLabelFont(42);
   alphaT1__22->GetXaxis()->SetLabelSize(0.035);
   alphaT1__22->GetXaxis()->SetTitleSize(0.035);
   alphaT1__22->GetXaxis()->SetTitleOffset(1);
   alphaT1__22->GetXaxis()->SetTitleFont(42);
   alphaT1__22->GetYaxis()->SetLabelFont(42);
   alphaT1__22->GetYaxis()->SetLabelSize(0.035);
   alphaT1__22->GetYaxis()->SetTitleSize(0.035);
   alphaT1__22->GetYaxis()->SetTitleFont(42);
   alphaT1__22->GetZaxis()->SetLabelFont(42);
   alphaT1__22->GetZaxis()->SetLabelSize(0.035);
   alphaT1__22->GetZaxis()->SetTitleSize(0.035);
   alphaT1__22->GetZaxis()->SetTitleOffset(1);
   alphaT1__22->GetZaxis()->SetTitleFont(42);
   alphaT1__22->Draw("hist SAME E");
   
   TH1F *alphaT2__23 = new TH1F("alphaT2__23","",25,0,1.25);
   alphaT2__23->SetBinContent(3,0.0001370897);
   alphaT2__23->SetBinContent(4,0.0005898003);
   alphaT2__23->SetBinContent(5,0.001348767);
   alphaT2__23->SetBinContent(6,0.00882848);
   alphaT2__23->SetBinContent(7,0.04275706);
   alphaT2__23->SetBinContent(8,0.1294783);
   alphaT2__23->SetBinContent(9,0.2672014);
   alphaT2__23->SetBinContent(10,0.4137777);
   alphaT2__23->SetBinContent(11,0.1073028);
   alphaT2__23->SetBinContent(12,0.01400538);
   alphaT2__23->SetBinContent(13,0.004720047);
   alphaT2__23->SetBinContent(14,0.002144516);
   alphaT2__23->SetBinContent(15,0.002677866);
   alphaT2__23->SetBinContent(16,0.0008251562);
   alphaT2__23->SetBinContent(17,0.0006110577);
   alphaT2__23->SetBinContent(18,0.0006741749);
   alphaT2__23->SetBinContent(19,0.0006575893);
   alphaT2__23->SetBinContent(20,0.0002719565);
   alphaT2__23->SetBinContent(21,0.0003560615);
   alphaT2__23->SetBinContent(22,0.0003528697);
   alphaT2__23->SetBinContent(23,0.00038495);
   alphaT2__23->SetBinContent(24,0.0005396196);
   alphaT2__23->SetBinContent(25,0.0003573751);
   alphaT2__23->SetBinContent(26,0.007160985);
   alphaT2__23->SetBinError(3,9.693796e-05);
   alphaT2__23->SetBinError(4,0.0001971922);
   alphaT2__23->SetBinError(5,0.0002771757);
   alphaT2__23->SetBinError(6,0.0007086156);
   alphaT2__23->SetBinError(7,0.001562816);
   alphaT2__23->SetBinError(8,0.002732024);
   alphaT2__23->SetBinError(9,0.003931258);
   alphaT2__23->SetBinError(10,0.004905767);
   alphaT2__23->SetBinError(11,0.002460045);
   alphaT2__23->SetBinError(12,0.0008825964);
   alphaT2__23->SetBinError(13,0.0005116658);
   alphaT2__23->SetBinError(14,0.0003519121);
   alphaT2__23->SetBinError(15,0.0003921107);
   alphaT2__23->SetBinError(16,0.0002218229);
   alphaT2__23->SetBinError(17,0.0001865277);
   alphaT2__23->SetBinError(18,0.0001969981);
   alphaT2__23->SetBinError(19,0.0002006877);
   alphaT2__23->SetBinError(20,0.0001219617);
   alphaT2__23->SetBinError(21,0.0001466063);
   alphaT2__23->SetBinError(22,0.0001448287);
   alphaT2__23->SetBinError(23,0.0001495404);
   alphaT2__23->SetBinError(24,0.0001804839);
   alphaT2__23->SetBinError(25,0.0001469868);
   alphaT2__23->SetBinError(26,0.0006536583);
   alphaT2__23->SetMinimum(0.0005);
   alphaT2__23->SetMaximum(0.2);
   alphaT2__23->SetEntries(17847);

   ci = TColor::GetColor("#6666ff");
   alphaT2__23->SetLineColor(ci);
   alphaT2__23->SetLineWidth(3);
   alphaT2__23->GetXaxis()->SetRange(1,26);
   alphaT2__23->GetXaxis()->SetLabelFont(42);
   alphaT2__23->GetXaxis()->SetLabelSize(0.035);
   alphaT2__23->GetXaxis()->SetTitleSize(0.035);
   alphaT2__23->GetXaxis()->SetTitleOffset(1);
   alphaT2__23->GetXaxis()->SetTitleFont(42);
   alphaT2__23->GetYaxis()->SetLabelFont(42);
   alphaT2__23->GetYaxis()->SetLabelSize(0.035);
   alphaT2__23->GetYaxis()->SetTitleSize(0.035);
   alphaT2__23->GetYaxis()->SetTitleFont(42);
   alphaT2__23->GetZaxis()->SetLabelFont(42);
   alphaT2__23->GetZaxis()->SetLabelSize(0.035);
   alphaT2__23->GetZaxis()->SetTitleSize(0.035);
   alphaT2__23->GetZaxis()->SetTitleOffset(1);
   alphaT2__23->GetZaxis()->SetTitleFont(42);
   alphaT2__23->Draw("hist SAME E");
   
   TH1F *alphaT3__24 = new TH1F("alphaT3__24","",25,0,1.25);
   alphaT3__24->SetBinContent(2,9.400307e-06);
   alphaT3__24->SetBinContent(3,3.792914e-05);
   alphaT3__24->SetBinContent(4,0.0005285988);
   alphaT3__24->SetBinContent(5,0.005855191);
   alphaT3__24->SetBinContent(6,0.02424989);
   alphaT3__24->SetBinContent(7,0.07084686);
   alphaT3__24->SetBinContent(8,0.1439648);
   alphaT3__24->SetBinContent(9,0.2191919);
   alphaT3__24->SetBinContent(10,0.27946);
   alphaT3__24->SetBinContent(11,0.1173923);
   alphaT3__24->SetBinContent(12,0.04587062);
   alphaT3__24->SetBinContent(13,0.0247272);
   alphaT3__24->SetBinContent(14,0.01650823);
   alphaT3__24->SetBinContent(15,0.01188433);
   alphaT3__24->SetBinContent(16,0.008286222);
   alphaT3__24->SetBinContent(17,0.006107347);
   alphaT3__24->SetBinContent(18,0.005051979);
   alphaT3__24->SetBinContent(19,0.004710097);
   alphaT3__24->SetBinContent(20,0.003752603);
   alphaT3__24->SetBinContent(21,0.002861466);
   alphaT3__24->SetBinContent(22,0.002680203);
   alphaT3__24->SetBinContent(23,0.002210675);
   alphaT3__24->SetBinContent(24,0.002042358);
   alphaT3__24->SetBinContent(25,0.001769815);
   alphaT3__24->SetBinContent(26,0.03425211);
   alphaT3__24->SetBinError(2,9.400307e-06);
   alphaT3__24->SetBinError(3,1.89651e-05);
   alphaT3__24->SetBinError(4,6.569208e-05);
   alphaT3__24->SetBinError(5,0.0002204871);
   alphaT3__24->SetBinError(6,0.0004501103);
   alphaT3__24->SetBinError(7,0.0007730437);
   alphaT3__24->SetBinError(8,0.001104245);
   alphaT3__24->SetBinError(9,0.00136527);
   alphaT3__24->SetBinError(10,0.001542752);
   alphaT3__24->SetBinError(11,0.0009962292);
   alphaT3__24->SetBinError(12,0.0006198256);
   alphaT3__24->SetBinError(13,0.0004555785);
   alphaT3__24->SetBinError(14,0.0003727982);
   alphaT3__24->SetBinError(15,0.000316164);
   alphaT3__24->SetBinError(16,0.0002641731);
   alphaT3__24->SetBinError(17,0.0002271948);
   alphaT3__24->SetBinError(18,0.0002071642);
   alphaT3__24->SetBinError(19,0.0002008242);
   alphaT3__24->SetBinError(20,0.0001792606);
   alphaT3__24->SetBinError(21,0.0001561917);
   alphaT3__24->SetBinError(22,0.0001515054);
   alphaT3__24->SetBinError(23,0.0001381873);
   alphaT3__24->SetBinError(24,0.00013225);
   alphaT3__24->SetBinError(25,0.0001237462);
   alphaT3__24->SetBinError(26,0.0005439501);
   alphaT3__24->SetMinimum(0.0005);
   alphaT3__24->SetMaximum(0.2);
   alphaT3__24->SetEntries(123638);

   ci = TColor::GetColor("#cc3399");
   alphaT3__24->SetLineColor(kPink+7);
   alphaT3__24->SetLineWidth(3);
   alphaT3__24->GetXaxis()->SetRange(1,26);
   alphaT3__24->GetXaxis()->SetLabelFont(42);
   alphaT3__24->GetXaxis()->SetLabelSize(0.035);
   alphaT3__24->GetXaxis()->SetTitleSize(0.035);
   alphaT3__24->GetXaxis()->SetTitleOffset(1);
   alphaT3__24->GetXaxis()->SetTitleFont(42);
   alphaT3__24->GetYaxis()->SetLabelFont(42);
   alphaT3__24->GetYaxis()->SetLabelSize(0.035);
   alphaT3__24->GetYaxis()->SetTitleSize(0.035);
   alphaT3__24->GetYaxis()->SetTitleFont(42);
   alphaT3__24->GetZaxis()->SetLabelFont(42);
   alphaT3__24->GetZaxis()->SetLabelSize(0.035);
   alphaT3__24->GetZaxis()->SetTitleSize(0.035);
   alphaT3__24->GetZaxis()->SetTitleOffset(1);
   alphaT3__24->GetZaxis()->SetTitleFont(42);
   alphaT3__24->Draw("hist SAME E");
   TText *text = new TText(1.275,0.00011,"OVERFLOW");
   text->SetTextAlign(12);
   text->SetTextSize(0.035);
   text->SetTextAngle(90);
   text->Draw();
   pad->Modified();
   c1->cd();
  
// ------------>Primitives in pad2: 
   pad2->Draw();
   pad2->cd();
   pad2->Range(0,0,1,1);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetFrameBorderMode(0);
   TLatex *   tex = new TLatex(0.8,0.4,"normalized no. of events");
   tex->SetTextFont(42);
   tex->SetTextSize(0.65);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
   pad2->Modified();
   c1->cd();
  
// ------------>Primitives in pad3: 
   pad3->Draw();
   pad3->cd();
   pad3->Range(0,0,1,1);
   pad3->SetFillColor(0);
   pad3->SetBorderMode(0);
   pad3->SetBorderSize(2);
   pad3->SetFrameBorderMode(0);
      tex = new TLatex(0.85,0.4,"#scale[1.2]{#alpha}_{T}");
   tex->SetTextFont(42);
   tex->SetTextSize(0.6);
   tex->SetLineWidth(2);
   tex->Draw();
   pad3->Modified();
   c1->cd();
  
   // ------------>Primitives in pad4: 
   pad4->Draw();
   pad4->cd();
   pad4->Range(0,0,1,1);
   pad4->SetFillColor(0);
   pad4->SetBorderMode(0);
   pad4->SetBorderSize(2);
   pad4->SetFrameBorderMode(0);
   
   TLegend *leg = new TLegend(0.02,0.362,0.33,1.0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.4);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(3);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("MET0","HF background","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   tex = new TLatex(0.025,0.15,"signal (m_{c} [GeV], #Deltam [GeV])");
   tex->SetTextFont(42);
   tex->SetTextSize(0.4);
   tex->SetLineWidth(2);
   tex->Draw();
      
   float legStart = 0.4;
   float legDiff = 0.36;

   leg = new TLegend(legStart,0,legStart+legDiff,1.0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.4);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(3);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("MET1","DM: (220, 20)","l");
   entry->SetLineColor(kBlue);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("MET2","DM: (324, 20)","l");
   entry->SetLineColor(kBlue-7);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();

   float correct = -0.05;
   leg = new TLegend(legStart+legDiff+correct,0,legStart+2*legDiff+correct,1.0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.4);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(3);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("MET4","(220, 40)","l");
   entry->SetLineColor(kPink+7);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();

   pad4->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->SaveAs("alphaT_manuallyEdited.pdf");
}
