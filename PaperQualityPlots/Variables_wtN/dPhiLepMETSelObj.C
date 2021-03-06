void dPhiLepMETSelObj()
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Jun  4 12:26:19 2020) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,782,600);
   gStyle->SetOptStat(0);
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: 
   TPad *pad = new TPad("", "",0.09,0.08,1,0.95);
   pad->Draw();
   pad->cd();
   pad->Range(0,0,1,1);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetLogy();
   pad->SetFrameBorderMode(0);
   
   TH1F *dPhiLepMETSelObj0__9 = new TH1F("dPhiLepMETSelObj0__9","",50,0,3.141593);
   dPhiLepMETSelObj0__9->SetBinContent(1,0.2453749);
   dPhiLepMETSelObj0__9->SetBinContent(2,0.04941577);
   dPhiLepMETSelObj0__9->SetBinContent(3,0.04576436);
   dPhiLepMETSelObj0__9->SetBinContent(4,0.04479065);
   dPhiLepMETSelObj0__9->SetBinContent(5,0.03675755);
   dPhiLepMETSelObj0__9->SetBinContent(6,0.03675755);
   dPhiLepMETSelObj0__9->SetBinContent(7,0.03919182);
   dPhiLepMETSelObj0__9->SetBinContent(8,0.03383642);
   dPhiLepMETSelObj0__9->SetBinContent(9,0.03602726);
   dPhiLepMETSelObj0__9->SetBinContent(10,0.030185);
   dPhiLepMETSelObj0__9->SetBinContent(11,0.02775073);
   dPhiLepMETSelObj0__9->SetBinContent(12,0.03164557);
   dPhiLepMETSelObj0__9->SetBinContent(13,0.031889);
   dPhiLepMETSelObj0__9->SetBinContent(14,0.03261928);
   dPhiLepMETSelObj0__9->SetBinContent(15,0.02775073);
   dPhiLepMETSelObj0__9->SetBinContent(16,0.03140214);
   dPhiLepMETSelObj0__9->SetBinContent(17,0.03115872);
   dPhiLepMETSelObj0__9->SetBinContent(18,0.02775073);
   dPhiLepMETSelObj0__9->SetBinContent(19,0.03091529);
   dPhiLepMETSelObj0__9->SetBinContent(20,0.02604674);
   dPhiLepMETSelObj0__9->SetBinContent(21,0.02629017);
   dPhiLepMETSelObj0__9->SetBinContent(22,0.01898734);
   dPhiLepMETSelObj0__9->SetBinContent(23,0.01265823);
   dPhiLepMETSelObj0__9->SetBinContent(24,0.007059396);
   dPhiLepMETSelObj0__9->SetBinContent(25,0.00194742);
   dPhiLepMETSelObj0__9->SetBinContent(26,0.001460565);
   dPhiLepMETSelObj0__9->SetBinContent(27,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(28,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(29,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(30,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(31,0.0002434275);
   dPhiLepMETSelObj0__9->SetBinContent(32,0.001460565);
   dPhiLepMETSelObj0__9->SetBinContent(33,0.00194742);
   dPhiLepMETSelObj0__9->SetBinContent(34,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(35,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinContent(36,0.002434275);
   dPhiLepMETSelObj0__9->SetBinContent(37,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(38,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(39,0.002434275);
   dPhiLepMETSelObj0__9->SetBinContent(40,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(41,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(42,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(43,0.00194742);
   dPhiLepMETSelObj0__9->SetBinContent(44,0.0009737099);
   dPhiLepMETSelObj0__9->SetBinContent(45,0.002434275);
   dPhiLepMETSelObj0__9->SetBinContent(46,0.002190847);
   dPhiLepMETSelObj0__9->SetBinContent(47,0.001217137);
   dPhiLepMETSelObj0__9->SetBinContent(48,0.002677702);
   dPhiLepMETSelObj0__9->SetBinContent(49,0.001460565);
   dPhiLepMETSelObj0__9->SetBinContent(50,0.001703992);
   dPhiLepMETSelObj0__9->SetBinError(1,0.007728582);
   dPhiLepMETSelObj0__9->SetBinError(2,0.003468307);
   dPhiLepMETSelObj0__9->SetBinError(3,0.003337709);
   dPhiLepMETSelObj0__9->SetBinError(4,0.003302011);
   dPhiLepMETSelObj0__9->SetBinError(5,0.002991287);
   dPhiLepMETSelObj0__9->SetBinError(6,0.002991287);
   dPhiLepMETSelObj0__9->SetBinError(7,0.003088748);
   dPhiLepMETSelObj0__9->SetBinError(8,0.002869967);
   dPhiLepMETSelObj0__9->SetBinError(9,0.002961423);
   dPhiLepMETSelObj0__9->SetBinError(10,0.002710693);
   dPhiLepMETSelObj0__9->SetBinError(11,0.002599094);
   dPhiLepMETSelObj0__9->SetBinError(12,0.0027755);
   dPhiLepMETSelObj0__9->SetBinError(13,0.002786155);
   dPhiLepMETSelObj0__9->SetBinError(14,0.002817877);
   dPhiLepMETSelObj0__9->SetBinError(15,0.002599094);
   dPhiLepMETSelObj0__9->SetBinError(16,0.002764804);
   dPhiLepMETSelObj0__9->SetBinError(17,0.002754067);
   dPhiLepMETSelObj0__9->SetBinError(18,0.002599094);
   dPhiLepMETSelObj0__9->SetBinError(19,0.002743288);
   dPhiLepMETSelObj0__9->SetBinError(20,0.002518033);
   dPhiLepMETSelObj0__9->SetBinError(21,0.002529772);
   dPhiLepMETSelObj0__9->SetBinError(22,0.002149893);
   dPhiLepMETSelObj0__9->SetBinError(23,0.00175538);
   dPhiLepMETSelObj0__9->SetBinError(24,0.001310897);
   dPhiLepMETSelObj0__9->SetBinError(25,0.0006885168);
   dPhiLepMETSelObj0__9->SetBinError(26,0.0005962731);
   dPhiLepMETSelObj0__9->SetBinError(27,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(28,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(29,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(30,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(31,0.0002434275);
   dPhiLepMETSelObj0__9->SetBinError(32,0.0005962731);
   dPhiLepMETSelObj0__9->SetBinError(33,0.0006885168);
   dPhiLepMETSelObj0__9->SetBinError(34,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(35,0.0003442584);
   dPhiLepMETSelObj0__9->SetBinError(36,0.0007697852);
   dPhiLepMETSelObj0__9->SetBinError(37,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(38,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(39,0.0007697852);
   dPhiLepMETSelObj0__9->SetBinError(40,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(41,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(42,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(43,0.0006885168);
   dPhiLepMETSelObj0__9->SetBinError(44,0.0004868549);
   dPhiLepMETSelObj0__9->SetBinError(45,0.0007697852);
   dPhiLepMETSelObj0__9->SetBinError(46,0.0007302824);
   dPhiLepMETSelObj0__9->SetBinError(47,0.0005443203);
   dPhiLepMETSelObj0__9->SetBinError(48,0.0008073575);
   dPhiLepMETSelObj0__9->SetBinError(49,0.0005962731);
   dPhiLepMETSelObj0__9->SetBinError(50,0.0006440485);
   dPhiLepMETSelObj0__9->SetMinimum(0.0005);
   dPhiLepMETSelObj0__9->SetMaximum(1);
   dPhiLepMETSelObj0__9->SetEntries(4108);
   dPhiLepMETSelObj0__9->SetLineWidth(3);
   dPhiLepMETSelObj0__9->GetXaxis()->SetLabelFont(42);
   dPhiLepMETSelObj0__9->GetXaxis()->SetLabelOffset(-0.07);
   dPhiLepMETSelObj0__9->GetXaxis()->SetLabelSize(0.055);
   dPhiLepMETSelObj0__9->GetXaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj0__9->GetXaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj0__9->GetXaxis()->SetTitleFont(42);
   dPhiLepMETSelObj0__9->GetYaxis()->SetLabelFont(42);
   dPhiLepMETSelObj0__9->GetYaxis()->SetLabelOffset(-0.05);
   dPhiLepMETSelObj0__9->GetYaxis()->SetLabelSize(0.055);
   dPhiLepMETSelObj0__9->GetYaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj0__9->GetYaxis()->SetTitleFont(42);
   dPhiLepMETSelObj0__9->GetZaxis()->SetLabelFont(42);
   dPhiLepMETSelObj0__9->GetZaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj0__9->GetZaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj0__9->GetZaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj0__9->GetZaxis()->SetTitleFont(42);
   dPhiLepMETSelObj0__9->Draw("hist  E");
   
   TH1F *dPhiLepMETSelObj1__10 = new TH1F("dPhiLepMETSelObj1__10","",50,0,3.141593);
   dPhiLepMETSelObj1__10->SetBinContent(1,0.05657154);
   dPhiLepMETSelObj1__10->SetBinContent(2,0.05351218);
   dPhiLepMETSelObj1__10->SetBinContent(3,0.05212447);
   dPhiLepMETSelObj1__10->SetBinContent(4,0.0487801);
   dPhiLepMETSelObj1__10->SetBinContent(5,0.0481527);
   dPhiLepMETSelObj1__10->SetBinContent(6,0.04466635);
   dPhiLepMETSelObj1__10->SetBinContent(7,0.04379641);
   dPhiLepMETSelObj1__10->SetBinContent(8,0.04573928);
   dPhiLepMETSelObj1__10->SetBinContent(9,0.04251098);
   dPhiLepMETSelObj1__10->SetBinContent(10,0.04054043);
   dPhiLepMETSelObj1__10->SetBinContent(11,0.04256396);
   dPhiLepMETSelObj1__10->SetBinContent(12,0.04444501);
   dPhiLepMETSelObj1__10->SetBinContent(13,0.04415809);
   dPhiLepMETSelObj1__10->SetBinContent(14,0.04269494);
   dPhiLepMETSelObj1__10->SetBinContent(15,0.04789648);
   dPhiLepMETSelObj1__10->SetBinContent(16,0.04682392);
   dPhiLepMETSelObj1__10->SetBinContent(17,0.04213271);
   dPhiLepMETSelObj1__10->SetBinContent(18,0.04282993);
   dPhiLepMETSelObj1__10->SetBinContent(19,0.0402153);
   dPhiLepMETSelObj1__10->SetBinContent(20,0.03440918);
   dPhiLepMETSelObj1__10->SetBinContent(21,0.02903221);
   dPhiLepMETSelObj1__10->SetBinContent(22,0.02154991);
   dPhiLepMETSelObj1__10->SetBinContent(23,0.01159414);
   dPhiLepMETSelObj1__10->SetBinContent(24,0.00651244);
   dPhiLepMETSelObj1__10->SetBinContent(25,0.002301537);
   dPhiLepMETSelObj1__10->SetBinContent(26,0.0009240628);
   dPhiLepMETSelObj1__10->SetBinContent(27,0.0007050614);
   dPhiLepMETSelObj1__10->SetBinContent(28,0.001389985);
   dPhiLepMETSelObj1__10->SetBinContent(29,0.0009381303);
   dPhiLepMETSelObj1__10->SetBinContent(30,0.0007621202);
   dPhiLepMETSelObj1__10->SetBinContent(31,0.0008262109);
   dPhiLepMETSelObj1__10->SetBinContent(32,0.001451418);
   dPhiLepMETSelObj1__10->SetBinContent(33,0.001053702);
   dPhiLepMETSelObj1__10->SetBinContent(34,0.0008267105);
   dPhiLepMETSelObj1__10->SetBinContent(35,0.0006076437);
   dPhiLepMETSelObj1__10->SetBinContent(36,0.001074907);
   dPhiLepMETSelObj1__10->SetBinContent(37,0.001282965);
   dPhiLepMETSelObj1__10->SetBinContent(38,0.001331162);
   dPhiLepMETSelObj1__10->SetBinContent(39,0.00105306);
   dPhiLepMETSelObj1__10->SetBinContent(40,0.0008519923);
   dPhiLepMETSelObj1__10->SetBinContent(41,0.0009340009);
   dPhiLepMETSelObj1__10->SetBinContent(42,0.001236651);
   dPhiLepMETSelObj1__10->SetBinContent(43,0.001113855);
   dPhiLepMETSelObj1__10->SetBinContent(44,0.001000644);
   dPhiLepMETSelObj1__10->SetBinContent(45,0.001368591);
   dPhiLepMETSelObj1__10->SetBinContent(46,0.0008896691);
   dPhiLepMETSelObj1__10->SetBinContent(47,0.0005075901);
   dPhiLepMETSelObj1__10->SetBinContent(48,0.0009192156);
   dPhiLepMETSelObj1__10->SetBinContent(49,0.000815213);
   dPhiLepMETSelObj1__10->SetBinContent(50,0.0005812243);
   dPhiLepMETSelObj1__10->SetBinError(1,0.002002033);
   dPhiLepMETSelObj1__10->SetBinError(2,0.001933163);
   dPhiLepMETSelObj1__10->SetBinError(3,0.0018965);
   dPhiLepMETSelObj1__10->SetBinError(4,0.001838031);
   dPhiLepMETSelObj1__10->SetBinError(5,0.001822561);
   dPhiLepMETSelObj1__10->SetBinError(6,0.001741324);
   dPhiLepMETSelObj1__10->SetBinError(7,0.001718265);
   dPhiLepMETSelObj1__10->SetBinError(8,0.001756641);
   dPhiLepMETSelObj1__10->SetBinError(9,0.001692648);
   dPhiLepMETSelObj1__10->SetBinError(10,0.001659275);
   dPhiLepMETSelObj1__10->SetBinError(11,0.001701972);
   dPhiLepMETSelObj1__10->SetBinError(12,0.001731716);
   dPhiLepMETSelObj1__10->SetBinError(13,0.00172041);
   dPhiLepMETSelObj1__10->SetBinError(14,0.001684417);
   dPhiLepMETSelObj1__10->SetBinError(15,0.001787942);
   dPhiLepMETSelObj1__10->SetBinError(16,0.001768369);
   dPhiLepMETSelObj1__10->SetBinError(17,0.001670708);
   dPhiLepMETSelObj1__10->SetBinError(18,0.001696281);
   dPhiLepMETSelObj1__10->SetBinError(19,0.001642343);
   dPhiLepMETSelObj1__10->SetBinError(20,0.001523421);
   dPhiLepMETSelObj1__10->SetBinError(21,0.001405475);
   dPhiLepMETSelObj1__10->SetBinError(22,0.001213605);
   dPhiLepMETSelObj1__10->SetBinError(23,0.0008975695);
   dPhiLepMETSelObj1__10->SetBinError(24,0.0006790142);
   dPhiLepMETSelObj1__10->SetBinError(25,0.0004010051);
   dPhiLepMETSelObj1__10->SetBinError(26,0.0002590509);
   dPhiLepMETSelObj1__10->SetBinError(27,0.0002194916);
   dPhiLepMETSelObj1__10->SetBinError(28,0.0003049192);
   dPhiLepMETSelObj1__10->SetBinError(29,0.000255989);
   dPhiLepMETSelObj1__10->SetBinError(30,0.0002240844);
   dPhiLepMETSelObj1__10->SetBinError(31,0.0002362859);
   dPhiLepMETSelObj1__10->SetBinError(32,0.0003103238);
   dPhiLepMETSelObj1__10->SetBinError(33,0.000268958);
   dPhiLepMETSelObj1__10->SetBinError(34,0.0002351261);
   dPhiLepMETSelObj1__10->SetBinError(35,0.0002036215);
   dPhiLepMETSelObj1__10->SetBinError(36,0.0002753243);
   dPhiLepMETSelObj1__10->SetBinError(37,0.0003037051);
   dPhiLepMETSelObj1__10->SetBinError(38,0.0003079614);
   dPhiLepMETSelObj1__10->SetBinError(39,0.0002732019);
   dPhiLepMETSelObj1__10->SetBinError(40,0.0002471607);
   dPhiLepMETSelObj1__10->SetBinError(41,0.0002527735);
   dPhiLepMETSelObj1__10->SetBinError(42,0.0002935722);
   dPhiLepMETSelObj1__10->SetBinError(43,0.0002754371);
   dPhiLepMETSelObj1__10->SetBinError(44,0.0002641856);
   dPhiLepMETSelObj1__10->SetBinError(45,0.0003190846);
   dPhiLepMETSelObj1__10->SetBinError(46,0.0002488685);
   dPhiLepMETSelObj1__10->SetBinError(47,0.0001811494);
   dPhiLepMETSelObj1__10->SetBinError(48,0.0002546929);
   dPhiLepMETSelObj1__10->SetBinError(49,0.0002388339);
   dPhiLepMETSelObj1__10->SetBinError(50,0.0001955195);
   dPhiLepMETSelObj1__10->SetMinimum(0.0005);
   dPhiLepMETSelObj1__10->SetMaximum(0.2);
   dPhiLepMETSelObj1__10->SetEntries(19133);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   dPhiLepMETSelObj1__10->SetLineColor(ci);
   dPhiLepMETSelObj1__10->SetLineWidth(3);
   dPhiLepMETSelObj1__10->GetXaxis()->SetLabelFont(42);
   dPhiLepMETSelObj1__10->GetXaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj1__10->GetXaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj1__10->GetXaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj1__10->GetXaxis()->SetTitleFont(42);
   dPhiLepMETSelObj1__10->GetYaxis()->SetLabelFont(42);
   dPhiLepMETSelObj1__10->GetYaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj1__10->GetYaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj1__10->GetYaxis()->SetTitleFont(42);
   dPhiLepMETSelObj1__10->GetZaxis()->SetLabelFont(42);
   dPhiLepMETSelObj1__10->GetZaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj1__10->GetZaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj1__10->GetZaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj1__10->GetZaxis()->SetTitleFont(42);
   dPhiLepMETSelObj1__10->Draw("hist SAME E");
   
   TH1F *dPhiLepMETSelObj2__11 = new TH1F("dPhiLepMETSelObj2__11","",50,0,3.141593);
   dPhiLepMETSelObj2__11->SetBinContent(1,0.04693467);
   dPhiLepMETSelObj2__11->SetBinContent(2,0.04919487);
   dPhiLepMETSelObj2__11->SetBinContent(3,0.04869206);
   dPhiLepMETSelObj2__11->SetBinContent(4,0.04573451);
   dPhiLepMETSelObj2__11->SetBinContent(5,0.04362616);
   dPhiLepMETSelObj2__11->SetBinContent(6,0.04100103);
   dPhiLepMETSelObj2__11->SetBinContent(7,0.04288911);
   dPhiLepMETSelObj2__11->SetBinContent(8,0.0433395);
   dPhiLepMETSelObj2__11->SetBinContent(9,0.04076611);
   dPhiLepMETSelObj2__11->SetBinContent(10,0.04420128);
   dPhiLepMETSelObj2__11->SetBinContent(11,0.04329301);
   dPhiLepMETSelObj2__11->SetBinContent(12,0.04616644);
   dPhiLepMETSelObj2__11->SetBinContent(13,0.04375425);
   dPhiLepMETSelObj2__11->SetBinContent(14,0.04583018);
   dPhiLepMETSelObj2__11->SetBinContent(15,0.047338);
   dPhiLepMETSelObj2__11->SetBinContent(16,0.04832962);
   dPhiLepMETSelObj2__11->SetBinContent(17,0.0470207);
   dPhiLepMETSelObj2__11->SetBinContent(18,0.04576373);
   dPhiLepMETSelObj2__11->SetBinContent(19,0.04354521);
   dPhiLepMETSelObj2__11->SetBinContent(20,0.03777711);
   dPhiLepMETSelObj2__11->SetBinContent(21,0.02943534);
   dPhiLepMETSelObj2__11->SetBinContent(22,0.02317625);
   dPhiLepMETSelObj2__11->SetBinContent(23,0.0131967);
   dPhiLepMETSelObj2__11->SetBinContent(24,0.006789755);
   dPhiLepMETSelObj2__11->SetBinContent(25,0.002454981);
   dPhiLepMETSelObj2__11->SetBinContent(26,0.001099405);
   dPhiLepMETSelObj2__11->SetBinContent(27,0.00187303);
   dPhiLepMETSelObj2__11->SetBinContent(28,0.001096848);
   dPhiLepMETSelObj2__11->SetBinContent(29,0.001735942);
   dPhiLepMETSelObj2__11->SetBinContent(30,0.001801363);
   dPhiLepMETSelObj2__11->SetBinContent(31,0.0003975682);
   dPhiLepMETSelObj2__11->SetBinContent(32,0.001228946);
   dPhiLepMETSelObj2__11->SetBinContent(33,0.0008666629);
   dPhiLepMETSelObj2__11->SetBinContent(34,0.001544445);
   dPhiLepMETSelObj2__11->SetBinContent(35,0.001454313);
   dPhiLepMETSelObj2__11->SetBinContent(36,0.001055666);
   dPhiLepMETSelObj2__11->SetBinContent(37,0.001086517);
   dPhiLepMETSelObj2__11->SetBinContent(38,0.0006897373);
   dPhiLepMETSelObj2__11->SetBinContent(39,0.001475266);
   dPhiLepMETSelObj2__11->SetBinContent(40,0.001034462);
   dPhiLepMETSelObj2__11->SetBinContent(41,0.00130463);
   dPhiLepMETSelObj2__11->SetBinContent(42,0.001169473);
   dPhiLepMETSelObj2__11->SetBinContent(43,0.001047885);
   dPhiLepMETSelObj2__11->SetBinContent(44,0.001153831);
   dPhiLepMETSelObj2__11->SetBinContent(45,0.001104896);
   dPhiLepMETSelObj2__11->SetBinContent(46,0.0008332531);
   dPhiLepMETSelObj2__11->SetBinContent(47,0.0008502761);
   dPhiLepMETSelObj2__11->SetBinContent(48,0.001235738);
   dPhiLepMETSelObj2__11->SetBinContent(49,0.0007964575);
   dPhiLepMETSelObj2__11->SetBinContent(50,0.001812797);
   dPhiLepMETSelObj2__11->SetBinError(1,0.001753002);
   dPhiLepMETSelObj2__11->SetBinError(2,0.001793276);
   dPhiLepMETSelObj2__11->SetBinError(3,0.001765957);
   dPhiLepMETSelObj2__11->SetBinError(4,0.001716403);
   dPhiLepMETSelObj2__11->SetBinError(5,0.001672944);
   dPhiLepMETSelObj2__11->SetBinError(6,0.001620672);
   dPhiLepMETSelObj2__11->SetBinError(7,0.001659432);
   dPhiLepMETSelObj2__11->SetBinError(8,0.001665955);
   dPhiLepMETSelObj2__11->SetBinError(9,0.001612767);
   dPhiLepMETSelObj2__11->SetBinError(10,0.001679109);
   dPhiLepMETSelObj2__11->SetBinError(11,0.00166014);
   dPhiLepMETSelObj2__11->SetBinError(12,0.001710622);
   dPhiLepMETSelObj2__11->SetBinError(13,0.001666583);
   dPhiLepMETSelObj2__11->SetBinError(14,0.00170709);
   dPhiLepMETSelObj2__11->SetBinError(15,0.001723988);
   dPhiLepMETSelObj2__11->SetBinError(16,0.001743767);
   dPhiLepMETSelObj2__11->SetBinError(17,0.001723757);
   dPhiLepMETSelObj2__11->SetBinError(18,0.001706742);
   dPhiLepMETSelObj2__11->SetBinError(19,0.00166955);
   dPhiLepMETSelObj2__11->SetBinError(20,0.001548322);
   dPhiLepMETSelObj2__11->SetBinError(21,0.001371975);
   dPhiLepMETSelObj2__11->SetBinError(22,0.001227249);
   dPhiLepMETSelObj2__11->SetBinError(23,0.0009284185);
   dPhiLepMETSelObj2__11->SetBinError(24,0.0006657375);
   dPhiLepMETSelObj2__11->SetBinError(25,0.0003939562);
   dPhiLepMETSelObj2__11->SetBinError(26,0.0002566896);
   dPhiLepMETSelObj2__11->SetBinError(27,0.0003524857);
   dPhiLepMETSelObj2__11->SetBinError(28,0.0002612308);
   dPhiLepMETSelObj2__11->SetBinError(29,0.0003391996);
   dPhiLepMETSelObj2__11->SetBinError(30,0.0003446605);
   dPhiLepMETSelObj2__11->SetBinError(31,0.0001571248);
   dPhiLepMETSelObj2__11->SetBinError(32,0.0002806431);
   dPhiLepMETSelObj2__11->SetBinError(33,0.0002285413);
   dPhiLepMETSelObj2__11->SetBinError(34,0.0003118091);
   dPhiLepMETSelObj2__11->SetBinError(35,0.0003071851);
   dPhiLepMETSelObj2__11->SetBinError(36,0.0002515015);
   dPhiLepMETSelObj2__11->SetBinError(37,0.0002648692);
   dPhiLepMETSelObj2__11->SetBinError(38,0.0002132235);
   dPhiLepMETSelObj2__11->SetBinError(39,0.0003078541);
   dPhiLepMETSelObj2__11->SetBinError(40,0.0002546648);
   dPhiLepMETSelObj2__11->SetBinError(41,0.0002873062);
   dPhiLepMETSelObj2__11->SetBinError(42,0.0002721531);
   dPhiLepMETSelObj2__11->SetBinError(43,0.0002587756);
   dPhiLepMETSelObj2__11->SetBinError(44,0.0002727906);
   dPhiLepMETSelObj2__11->SetBinError(45,0.0002694792);
   dPhiLepMETSelObj2__11->SetBinError(46,0.0002385835);
   dPhiLepMETSelObj2__11->SetBinError(47,0.0002367465);
   dPhiLepMETSelObj2__11->SetBinError(48,0.0002818653);
   dPhiLepMETSelObj2__11->SetBinError(49,0.0002177434);
   dPhiLepMETSelObj2__11->SetBinError(50,0.0003449253);
   dPhiLepMETSelObj2__11->SetMinimum(0.0005);
   dPhiLepMETSelObj2__11->SetMaximum(0.2);
   dPhiLepMETSelObj2__11->SetEntries(17847);

   ci = TColor::GetColor("#6666ff");
   dPhiLepMETSelObj2__11->SetLineColor(ci);
   dPhiLepMETSelObj2__11->SetLineWidth(3);
   dPhiLepMETSelObj2__11->GetXaxis()->SetLabelFont(42);
   dPhiLepMETSelObj2__11->GetXaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj2__11->GetXaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj2__11->GetXaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj2__11->GetXaxis()->SetTitleFont(42);
   dPhiLepMETSelObj2__11->GetYaxis()->SetLabelFont(42);
   dPhiLepMETSelObj2__11->GetYaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj2__11->GetYaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj2__11->GetYaxis()->SetTitleFont(42);
   dPhiLepMETSelObj2__11->GetZaxis()->SetLabelFont(42);
   dPhiLepMETSelObj2__11->GetZaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj2__11->GetZaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj2__11->GetZaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj2__11->GetZaxis()->SetTitleFont(42);
   dPhiLepMETSelObj2__11->Draw("hist SAME E");
   
   TH1F *dPhiLepMETSelObj3__12 = new TH1F("dPhiLepMETSelObj3__12","",50,0,3.141593);
   dPhiLepMETSelObj3__12->SetBinContent(1,0.04626237);
   dPhiLepMETSelObj3__12->SetBinContent(2,0.04346067);
   dPhiLepMETSelObj3__12->SetBinContent(3,0.04484943);
   dPhiLepMETSelObj3__12->SetBinContent(4,0.04531748);
   dPhiLepMETSelObj3__12->SetBinContent(5,0.04808355);
   dPhiLepMETSelObj3__12->SetBinContent(6,0.04985348);
   dPhiLepMETSelObj3__12->SetBinContent(7,0.0526091);
   dPhiLepMETSelObj3__12->SetBinContent(8,0.05599964);
   dPhiLepMETSelObj3__12->SetBinContent(9,0.056127);
   dPhiLepMETSelObj3__12->SetBinContent(10,0.05885225);
   dPhiLepMETSelObj3__12->SetBinContent(11,0.05909697);
   dPhiLepMETSelObj3__12->SetBinContent(12,0.05935093);
   dPhiLepMETSelObj3__12->SetBinContent(13,0.05611105);
   dPhiLepMETSelObj3__12->SetBinContent(14,0.05294034);
   dPhiLepMETSelObj3__12->SetBinContent(15,0.05006793);
   dPhiLepMETSelObj3__12->SetBinContent(16,0.04437699);
   dPhiLepMETSelObj3__12->SetBinContent(17,0.03921176);
   dPhiLepMETSelObj3__12->SetBinContent(18,0.03310524);
   dPhiLepMETSelObj3__12->SetBinContent(19,0.02766237);
   dPhiLepMETSelObj3__12->SetBinContent(20,0.02133344);
   dPhiLepMETSelObj3__12->SetBinContent(21,0.01506845);
   dPhiLepMETSelObj3__12->SetBinContent(22,0.01050835);
   dPhiLepMETSelObj3__12->SetBinContent(23,0.006314442);
   dPhiLepMETSelObj3__12->SetBinContent(24,0.003138533);
   dPhiLepMETSelObj3__12->SetBinContent(25,0.001502816);
   dPhiLepMETSelObj3__12->SetBinContent(26,0.001066455);
   dPhiLepMETSelObj3__12->SetBinContent(27,0.001059278);
   dPhiLepMETSelObj3__12->SetBinContent(28,0.001000778);
   dPhiLepMETSelObj3__12->SetBinContent(29,0.0009702494);
   dPhiLepMETSelObj3__12->SetBinContent(30,0.0009740559);
   dPhiLepMETSelObj3__12->SetBinContent(31,0.000698045);
   dPhiLepMETSelObj3__12->SetBinContent(32,0.0007175894);
   dPhiLepMETSelObj3__12->SetBinContent(33,0.0008109839);
   dPhiLepMETSelObj3__12->SetBinContent(34,0.0006590003);
   dPhiLepMETSelObj3__12->SetBinContent(35,0.0008250245);
   dPhiLepMETSelObj3__12->SetBinContent(36,0.0006485906);
   dPhiLepMETSelObj3__12->SetBinContent(37,0.0008041514);
   dPhiLepMETSelObj3__12->SetBinContent(38,0.0007594578);
   dPhiLepMETSelObj3__12->SetBinContent(39,0.0006700605);
   dPhiLepMETSelObj3__12->SetBinContent(40,0.0006817272);
   dPhiLepMETSelObj3__12->SetBinContent(41,0.0006483002);
   dPhiLepMETSelObj3__12->SetBinContent(42,0.0006784866);
   dPhiLepMETSelObj3__12->SetBinContent(43,0.0006713576);
   dPhiLepMETSelObj3__12->SetBinContent(44,0.000702826);
   dPhiLepMETSelObj3__12->SetBinContent(45,0.000535069);
   dPhiLepMETSelObj3__12->SetBinContent(46,0.0006304521);
   dPhiLepMETSelObj3__12->SetBinContent(47,0.0005633726);
   dPhiLepMETSelObj3__12->SetBinContent(48,0.0007628028);
   dPhiLepMETSelObj3__12->SetBinContent(49,0.0006241706);
   dPhiLepMETSelObj3__12->SetBinContent(50,0.0006331438);
   dPhiLepMETSelObj3__12->SetBinError(1,0.0006324084);
   dPhiLepMETSelObj3__12->SetBinError(2,0.0006114606);
   dPhiLepMETSelObj3__12->SetBinError(3,0.0006199118);
   dPhiLepMETSelObj3__12->SetBinError(4,0.000622168);
   dPhiLepMETSelObj3__12->SetBinError(5,0.0006405568);
   dPhiLepMETSelObj3__12->SetBinError(6,0.0006507543);
   dPhiLepMETSelObj3__12->SetBinError(7,0.0006680745);
   dPhiLepMETSelObj3__12->SetBinError(8,0.0006882562);
   dPhiLepMETSelObj3__12->SetBinError(9,0.0006887278);
   dPhiLepMETSelObj3__12->SetBinError(10,0.0007054407);
   dPhiLepMETSelObj3__12->SetBinError(11,0.0007066553);
   dPhiLepMETSelObj3__12->SetBinError(12,0.0007085263);
   dPhiLepMETSelObj3__12->SetBinError(13,0.0006898013);
   dPhiLepMETSelObj3__12->SetBinError(14,0.0006705099);
   dPhiLepMETSelObj3__12->SetBinError(15,0.0006530521);
   dPhiLepMETSelObj3__12->SetBinError(16,0.0006139228);
   dPhiLepMETSelObj3__12->SetBinError(17,0.0005781721);
   dPhiLepMETSelObj3__12->SetBinError(18,0.0005310062);
   dPhiLepMETSelObj3__12->SetBinError(19,0.0004861696);
   dPhiLepMETSelObj3__12->SetBinError(20,0.0004278948);
   dPhiLepMETSelObj3__12->SetBinError(21,0.0003593051);
   dPhiLepMETSelObj3__12->SetBinError(22,0.0003000621);
   dPhiLepMETSelObj3__12->SetBinError(23,0.0002325631);
   dPhiLepMETSelObj3__12->SetBinError(24,0.0001641715);
   dPhiLepMETSelObj3__12->SetBinError(25,0.000112868);
   dPhiLepMETSelObj3__12->SetBinError(26,9.54087e-05);
   dPhiLepMETSelObj3__12->SetBinError(27,9.473869e-05);
   dPhiLepMETSelObj3__12->SetBinError(28,9.215147e-05);
   dPhiLepMETSelObj3__12->SetBinError(29,9.151143e-05);
   dPhiLepMETSelObj3__12->SetBinError(30,9.166586e-05);
   dPhiLepMETSelObj3__12->SetBinError(31,7.751997e-05);
   dPhiLepMETSelObj3__12->SetBinError(32,7.777534e-05);
   dPhiLepMETSelObj3__12->SetBinError(33,8.381774e-05);
   dPhiLepMETSelObj3__12->SetBinError(34,7.52418e-05);
   dPhiLepMETSelObj3__12->SetBinError(35,8.391433e-05);
   dPhiLepMETSelObj3__12->SetBinError(36,7.502851e-05);
   dPhiLepMETSelObj3__12->SetBinError(37,8.311964e-05);
   dPhiLepMETSelObj3__12->SetBinError(38,8.087558e-05);
   dPhiLepMETSelObj3__12->SetBinError(39,7.52133e-05);
   dPhiLepMETSelObj3__12->SetBinError(40,7.623714e-05);
   dPhiLepMETSelObj3__12->SetBinError(41,7.42611e-05);
   dPhiLepMETSelObj3__12->SetBinError(42,7.721566e-05);
   dPhiLepMETSelObj3__12->SetBinError(43,7.598894e-05);
   dPhiLepMETSelObj3__12->SetBinError(44,7.764323e-05);
   dPhiLepMETSelObj3__12->SetBinError(45,6.812782e-05);
   dPhiLepMETSelObj3__12->SetBinError(46,7.371446e-05);
   dPhiLepMETSelObj3__12->SetBinError(47,6.9847e-05);
   dPhiLepMETSelObj3__12->SetBinError(48,8.143953e-05);
   dPhiLepMETSelObj3__12->SetBinError(49,7.300414e-05);
   dPhiLepMETSelObj3__12->SetBinError(50,7.354375e-05);
   dPhiLepMETSelObj3__12->SetMinimum(0.0005);
   dPhiLepMETSelObj3__12->SetMaximum(0.2);
   dPhiLepMETSelObj3__12->SetEntries(123638);

   ci = TColor::GetColor("#cc3399");
   dPhiLepMETSelObj3__12->SetLineColor(ci);
   dPhiLepMETSelObj3__12->SetLineWidth(3);
   dPhiLepMETSelObj3__12->GetXaxis()->SetLabelFont(42);
   dPhiLepMETSelObj3__12->GetXaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj3__12->GetXaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj3__12->GetXaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj3__12->GetXaxis()->SetTitleFont(42);
   dPhiLepMETSelObj3__12->GetYaxis()->SetLabelFont(42);
   dPhiLepMETSelObj3__12->GetYaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj3__12->GetYaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj3__12->GetYaxis()->SetTitleFont(42);
   dPhiLepMETSelObj3__12->GetZaxis()->SetLabelFont(42);
   dPhiLepMETSelObj3__12->GetZaxis()->SetLabelSize(0.035);
   dPhiLepMETSelObj3__12->GetZaxis()->SetTitleSize(0.035);
   dPhiLepMETSelObj3__12->GetZaxis()->SetTitleOffset(1);
   dPhiLepMETSelObj3__12->GetZaxis()->SetTitleFont(42);
   dPhiLepMETSelObj3__12->Draw("hist SAME E");
   pad->Modified();
   c1->cd();
  
// ------------>Primitives in pad: 
   TPad *pad = new TPad("", "",0.001,0.08,0.08,0.85);
   pad->Draw();
   pad->cd();
   pad->Range(0,0,1,1);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetFrameBorderMode(0);
   TLatex *   tex = new TLatex(0.7,0.3,"normalized no. of events");
   tex->SetTextFont(42);
   tex->SetTextSize(0.45);
   tex->SetTextAngle(90);
   tex->SetLineWidth(2);
   tex->Draw();
   pad->Modified();
   c1->cd();
  
// ------------>Primitives in pad: 
   TPad *pad = new TPad("", "",0.08,0.001,1,0.08);
   pad->Draw();
   pad->cd();
   pad->Range(0,0,1,1);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetFrameBorderMode(0);
      tex = new TLatex(0.5,0.4,"#Delta#phi(l_{1}, #left|#Sigma_{lep+jet}p_{T}#right|)");
   tex->SetTextFont(42);
   tex->SetTextSize(0.6);
   tex->SetLineWidth(2);
   tex->Draw();
   pad->Modified();
   c1->cd();
  
// ------------>Primitives in pad: 
   TPad *pad = new TPad("", "",0.001,0.88,1,0.98);
   pad->Draw();
   pad->cd();
   pad->Range(0,0,1,1);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetFrameBorderMode(0);
   
   TLegend *leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.4);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("dPhiLepMETSelObj0","HF Background","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
      tex = new TLatex(0.12,0.15,"Signal (m_{c}, #Deltam, c#scale[1.2]{#tau}_{c})");
   tex->SetTextFont(42);
   tex->SetTextSize(0.4);
   tex->SetLineWidth(2);
   tex->Draw();
   
   leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.4);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("dPhiLepMETSelObj1","(220, 20, DM)","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("dPhiLepMETSelObj2","(324, 20, DM)","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   pad->Modified();
   c1->cd();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
