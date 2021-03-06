void dPhiLepMETSelObj()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed May 20 12:40:57 2020) by ROOT version 6.18/04
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
   dPhiLepMETSelObj1__10->SetBinContent(1,0.05174306);
   dPhiLepMETSelObj1__10->SetBinContent(2,0.04991376);
   dPhiLepMETSelObj1__10->SetBinContent(3,0.04986149);
   dPhiLepMETSelObj1__10->SetBinContent(4,0.0460461);
   dPhiLepMETSelObj1__10->SetBinContent(5,0.04750954);
   dPhiLepMETSelObj1__10->SetBinContent(6,0.04494852);
   dPhiLepMETSelObj1__10->SetBinContent(7,0.04390321);
   dPhiLepMETSelObj1__10->SetBinContent(8,0.04568024);
   dPhiLepMETSelObj1__10->SetBinContent(9,0.04280562);
   dPhiLepMETSelObj1__10->SetBinContent(10,0.04113312);
   dPhiLepMETSelObj1__10->SetBinContent(11,0.04379867);
   dPhiLepMETSelObj1__10->SetBinContent(12,0.04426906);
   dPhiLepMETSelObj1__10->SetBinContent(13,0.04578477);
   dPhiLepMETSelObj1__10->SetBinContent(14,0.04526211);
   dPhiLepMETSelObj1__10->SetBinContent(15,0.05001829);
   dPhiLepMETSelObj1__10->SetBinContent(16,0.04865938);
   dPhiLepMETSelObj1__10->SetBinContent(17,0.04573251);
   dPhiLepMETSelObj1__10->SetBinContent(18,0.04458266);
   dPhiLepMETSelObj1__10->SetBinContent(19,0.04097632);
   dPhiLepMETSelObj1__10->SetBinContent(20,0.03548842);
   dPhiLepMETSelObj1__10->SetBinContent(21,0.02817122);
   dPhiLepMETSelObj1__10->SetBinContent(22,0.02116762);
   dPhiLepMETSelObj1__10->SetBinContent(23,0.01108033);
   dPhiLepMETSelObj1__10->SetBinContent(24,0.005592432);
   dPhiLepMETSelObj1__10->SetBinContent(25,0.00219516);
   dPhiLepMETSelObj1__10->SetBinContent(26,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(27,0.0008362515);
   dPhiLepMETSelObj1__10->SetBinContent(28,0.00146344);
   dPhiLepMETSelObj1__10->SetBinContent(29,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(30,0.0008885172);
   dPhiLepMETSelObj1__10->SetBinContent(31,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(32,0.001306643);
   dPhiLepMETSelObj1__10->SetBinContent(33,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(34,0.0007839858);
   dPhiLepMETSelObj1__10->SetBinContent(35,0.0007839858);
   dPhiLepMETSelObj1__10->SetBinContent(36,0.001045314);
   dPhiLepMETSelObj1__10->SetBinContent(37,0.001254377);
   dPhiLepMETSelObj1__10->SetBinContent(38,0.001254377);
   dPhiLepMETSelObj1__10->SetBinContent(39,0.0008885172);
   dPhiLepMETSelObj1__10->SetBinContent(40,0.0007839858);
   dPhiLepMETSelObj1__10->SetBinContent(41,0.0009930487);
   dPhiLepMETSelObj1__10->SetBinContent(42,0.001202111);
   dPhiLepMETSelObj1__10->SetBinContent(43,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(44,0.0009407829);
   dPhiLepMETSelObj1__10->SetBinContent(45,0.001045314);
   dPhiLepMETSelObj1__10->SetBinContent(46,0.0008362515);
   dPhiLepMETSelObj1__10->SetBinContent(47,0.0005226572);
   dPhiLepMETSelObj1__10->SetBinContent(48,0.0008362515);
   dPhiLepMETSelObj1__10->SetBinContent(49,0.0007839858);
   dPhiLepMETSelObj1__10->SetBinContent(50,0.0005226572);
   dPhiLepMETSelObj1__10->SetBinError(1,0.001644502);
   dPhiLepMETSelObj1__10->SetBinError(2,0.001615171);
   dPhiLepMETSelObj1__10->SetBinError(3,0.001614326);
   dPhiLepMETSelObj1__10->SetBinError(4,0.001551332);
   dPhiLepMETSelObj1__10->SetBinError(5,0.001575792);
   dPhiLepMETSelObj1__10->SetBinError(6,0.001532732);
   dPhiLepMETSelObj1__10->SetBinError(7,0.001514804);
   dPhiLepMETSelObj1__10->SetBinError(8,0.001545157);
   dPhiLepMETSelObj1__10->SetBinError(9,0.00149575);
   dPhiLepMETSelObj1__10->SetBinError(10,0.001466237);
   dPhiLepMETSelObj1__10->SetBinError(11,0.001513);
   dPhiLepMETSelObj1__10->SetBinError(12,0.001521103);
   dPhiLepMETSelObj1__10->SetBinError(13,0.001546924);
   dPhiLepMETSelObj1__10->SetBinError(14,0.001538069);
   dPhiLepMETSelObj1__10->SetBinError(15,0.001616862);
   dPhiLepMETSelObj1__10->SetBinError(16,0.001594747);
   dPhiLepMETSelObj1__10->SetBinError(17,0.001546041);
   dPhiLepMETSelObj1__10->SetBinError(18,0.001526481);
   dPhiLepMETSelObj1__10->SetBinError(19,0.00146344);
   dPhiLepMETSelObj1__10->SetBinError(20,0.001361921);
   dPhiLepMETSelObj1__10->SetBinError(21,0.00121342);
   dPhiLepMETSelObj1__10->SetBinError(22,0.001051827);
   dPhiLepMETSelObj1__10->SetBinError(23,0.0007610004);
   dPhiLepMETSelObj1__10->SetBinError(24,0.0005406408);
   dPhiLepMETSelObj1__10->SetBinError(25,0.0003387206);
   dPhiLepMETSelObj1__10->SetBinError(26,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(27,0.0002090629);
   dPhiLepMETSelObj1__10->SetBinError(28,0.0002765642);
   dPhiLepMETSelObj1__10->SetBinError(29,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(30,0.0002154971);
   dPhiLepMETSelObj1__10->SetBinError(31,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(32,0.0002613286);
   dPhiLepMETSelObj1__10->SetBinError(33,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(34,0.0002024243);
   dPhiLepMETSelObj1__10->SetBinError(35,0.0002024243);
   dPhiLepMETSelObj1__10->SetBinError(36,0.0002337394);
   dPhiLepMETSelObj1__10->SetBinError(37,0.0002560487);
   dPhiLepMETSelObj1__10->SetBinError(38,0.0002560487);
   dPhiLepMETSelObj1__10->SetBinError(39,0.0002154971);
   dPhiLepMETSelObj1__10->SetBinError(40,0.0002024243);
   dPhiLepMETSelObj1__10->SetBinError(41,0.000227821);
   dPhiLepMETSelObj1__10->SetBinError(42,0.0002506576);
   dPhiLepMETSelObj1__10->SetBinError(43,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(44,0.0002217447);
   dPhiLepMETSelObj1__10->SetBinError(45,0.0002337394);
   dPhiLepMETSelObj1__10->SetBinError(46,0.0002090629);
   dPhiLepMETSelObj1__10->SetBinError(47,0.0001652787);
   dPhiLepMETSelObj1__10->SetBinError(48,0.0002090629);
   dPhiLepMETSelObj1__10->SetBinError(49,0.0002024243);
   dPhiLepMETSelObj1__10->SetBinError(50,0.0001652787);
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
   dPhiLepMETSelObj2__11->SetBinContent(1,0.04443324);
   dPhiLepMETSelObj2__11->SetBinContent(2,0.0468426);
   dPhiLepMETSelObj2__11->SetBinContent(3,0.04801928);
   dPhiLepMETSelObj2__11->SetBinContent(4,0.04465736);
   dPhiLepMETSelObj2__11->SetBinContent(5,0.0434807);
   dPhiLepMETSelObj2__11->SetBinContent(6,0.0408472);
   dPhiLepMETSelObj2__11->SetBinContent(7,0.04342467);
   dPhiLepMETSelObj2__11->SetBinContent(8,0.04286435);
   dPhiLepMETSelObj2__11->SetBinContent(9,0.04123943);
   dPhiLepMETSelObj2__11->SetBinContent(10,0.04409705);
   dPhiLepMETSelObj2__11->SetBinContent(11,0.04370482);
   dPhiLepMETSelObj2__11->SetBinContent(12,0.04706673);
   dPhiLepMETSelObj2__11->SetBinContent(13,0.04499356);
   dPhiLepMETSelObj2__11->SetBinContent(14,0.04617022);
   dPhiLepMETSelObj2__11->SetBinContent(15,0.04942007);
   dPhiLepMETSelObj2__11->SetBinContent(16,0.04970023);
   dPhiLepMETSelObj2__11->SetBinContent(17,0.04846753);
   dPhiLepMETSelObj2__11->SetBinContent(18,0.04560991);
   dPhiLepMETSelObj2__11->SetBinContent(19,0.04353673);
   dPhiLepMETSelObj2__11->SetBinContent(20,0.03832577);
   dPhiLepMETSelObj2__11->SetBinContent(21,0.03003306);
   dPhiLepMETSelObj2__11->SetBinContent(22,0.0221886);
   dPhiLepMETSelObj2__11->SetBinContent(23,0.0124951);
   dPhiLepMETSelObj2__11->SetBinContent(24,0.006611756);
   dPhiLepMETSelObj2__11->SetBinContent(25,0.002409369);
   dPhiLepMETSelObj2__11->SetBinContent(26,0.0012327);
   dPhiLepMETSelObj2__11->SetBinContent(27,0.001736987);
   dPhiLepMETSelObj2__11->SetBinContent(28,0.0012327);
   dPhiLepMETSelObj2__11->SetBinContent(29,0.001624923);
   dPhiLepMETSelObj2__11->SetBinContent(30,0.001624923);
   dPhiLepMETSelObj2__11->SetBinContent(31,0.0004482546);
   dPhiLepMETSelObj2__11->SetBinContent(32,0.001176668);
   dPhiLepMETSelObj2__11->SetBinContent(33,0.001120636);
   dPhiLepMETSelObj2__11->SetBinContent(34,0.001624923);
   dPhiLepMETSelObj2__11->SetBinContent(35,0.001456828);
   dPhiLepMETSelObj2__11->SetBinContent(36,0.001176668);
   dPhiLepMETSelObj2__11->SetBinContent(37,0.001008573);
   dPhiLepMETSelObj2__11->SetBinContent(38,0.0006723819);
   dPhiLepMETSelObj2__11->SetBinContent(39,0.001400796);
   dPhiLepMETSelObj2__11->SetBinContent(40,0.001176668);
   dPhiLepMETSelObj2__11->SetBinContent(41,0.001288732);
   dPhiLepMETSelObj2__11->SetBinContent(42,0.001064605);
   dPhiLepMETSelObj2__11->SetBinContent(43,0.001064605);
   dPhiLepMETSelObj2__11->SetBinContent(44,0.001064605);
   dPhiLepMETSelObj2__11->SetBinContent(45,0.001064605);
   dPhiLepMETSelObj2__11->SetBinContent(46,0.0007284138);
   dPhiLepMETSelObj2__11->SetBinContent(47,0.0007844456);
   dPhiLepMETSelObj2__11->SetBinContent(48,0.001120636);
   dPhiLepMETSelObj2__11->SetBinContent(49,0.0008965092);
   dPhiLepMETSelObj2__11->SetBinContent(50,0.001568891);
   dPhiLepMETSelObj2__11->SetBinError(1,0.001577871);
   dPhiLepMETSelObj2__11->SetBinError(2,0.001620085);
   dPhiLepMETSelObj2__11->SetBinError(3,0.001640307);
   dPhiLepMETSelObj2__11->SetBinError(4,0.001581845);
   dPhiLepMETSelObj2__11->SetBinError(5,0.001560866);
   dPhiLepMETSelObj2__11->SetBinError(6,0.001512859);
   dPhiLepMETSelObj2__11->SetBinError(7,0.00155986);
   dPhiLepMETSelObj2__11->SetBinError(8,0.001549764);
   dPhiLepMETSelObj2__11->SetBinError(9,0.001520105);
   dPhiLepMETSelObj2__11->SetBinError(10,0.00157189);
   dPhiLepMETSelObj2__11->SetBinError(11,0.001564884);
   dPhiLepMETSelObj2__11->SetBinError(12,0.001623957);
   dPhiLepMETSelObj2__11->SetBinError(13,0.001587788);
   dPhiLepMETSelObj2__11->SetBinError(14,0.001608416);
   dPhiLepMETSelObj2__11->SetBinError(15,0.00166406);
   dPhiLepMETSelObj2__11->SetBinError(16,0.00166877);
   dPhiLepMETSelObj2__11->SetBinError(17,0.001647945);
   dPhiLepMETSelObj2__11->SetBinError(18,0.001598626);
   dPhiLepMETSelObj2__11->SetBinError(19,0.001561871);
   dPhiLepMETSelObj2__11->SetBinError(20,0.001465422);
   dPhiLepMETSelObj2__11->SetBinError(21,0.001297231);
   dPhiLepMETSelObj2__11->SetBinError(22,0.001115019);
   dPhiLepMETSelObj2__11->SetBinError(23,0.0008367336);
   dPhiLepMETSelObj2__11->SetBinError(24,0.0006086614);
   dPhiLepMETSelObj2__11->SetBinError(25,0.0003674253);
   dPhiLepMETSelObj2__11->SetBinError(26,0.0002628126);
   dPhiLepMETSelObj2__11->SetBinError(27,0.000311972);
   dPhiLepMETSelObj2__11->SetBinError(28,0.0002628126);
   dPhiLepMETSelObj2__11->SetBinError(29,0.0003017406);
   dPhiLepMETSelObj2__11->SetBinError(30,0.0003017406);
   dPhiLepMETSelObj2__11->SetBinError(31,0.0001584819);
   dPhiLepMETSelObj2__11->SetBinError(32,0.0002567701);
   dPhiLepMETSelObj2__11->SetBinError(33,0.0002505819);
   dPhiLepMETSelObj2__11->SetBinError(34,0.0003017406);
   dPhiLepMETSelObj2__11->SetBinError(35,0.0002857074);
   dPhiLepMETSelObj2__11->SetBinError(36,0.0002567701);
   dPhiLepMETSelObj2__11->SetBinError(37,0.0002377229);
   dPhiLepMETSelObj2__11->SetBinError(38,0.0001940999);
   dPhiLepMETSelObj2__11->SetBinError(39,0.0002801591);
   dPhiLepMETSelObj2__11->SetBinError(40,0.0002567701);
   dPhiLepMETSelObj2__11->SetBinError(41,0.0002687192);
   dPhiLepMETSelObj2__11->SetBinError(42,0.0002442371);
   dPhiLepMETSelObj2__11->SetBinError(43,0.0002442371);
   dPhiLepMETSelObj2__11->SetBinError(44,0.0002442371);
   dPhiLepMETSelObj2__11->SetBinError(45,0.0002442371);
   dPhiLepMETSelObj2__11->SetBinError(46,0.0002020256);
   dPhiLepMETSelObj2__11->SetBinError(47,0.0002096519);
   dPhiLepMETSelObj2__11->SetBinError(48,0.0002505819);
   dPhiLepMETSelObj2__11->SetBinError(49,0.0002241273);
   dPhiLepMETSelObj2__11->SetBinError(50,0.0002964926);
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
   dPhiLepMETSelObj3__12->SetBinContent(1,0.0450347);
   dPhiLepMETSelObj3__12->SetBinContent(2,0.04286708);
   dPhiLepMETSelObj3__12->SetBinContent(3,0.04452515);
   dPhiLepMETSelObj3__12->SetBinContent(4,0.0452369);
   dPhiLepMETSelObj3__12->SetBinContent(5,0.04807584);
   dPhiLepMETSelObj3__12->SetBinContent(6,0.05031624);
   dPhiLepMETSelObj3__12->SetBinContent(7,0.05302577);
   dPhiLepMETSelObj3__12->SetBinContent(8,0.0566816);
   dPhiLepMETSelObj3__12->SetBinContent(9,0.0568191);
   dPhiLepMETSelObj3__12->SetBinContent(10,0.05956906);
   dPhiLepMETSelObj3__12->SetBinContent(11,0.05979553);
   dPhiLepMETSelObj3__12->SetBinContent(12,0.05976318);
   dPhiLepMETSelObj3__12->SetBinContent(13,0.05653602);
   dPhiLepMETSelObj3__12->SetBinContent(14,0.05305812);
   dPhiLepMETSelObj3__12->SetBinContent(15,0.04990375);
   dPhiLepMETSelObj3__12->SetBinContent(16,0.04431485);
   dPhiLepMETSelObj3__12->SetBinContent(17,0.03879066);
   dPhiLepMETSelObj3__12->SetBinContent(18,0.03302383);
   dPhiLepMETSelObj3__12->SetBinContent(19,0.02720846);
   dPhiLepMETSelObj3__12->SetBinContent(20,0.02096443);
   dPhiLepMETSelObj3__12->SetBinContent(21,0.01480127);
   dPhiLepMETSelObj3__12->SetBinContent(22,0.01039325);
   dPhiLepMETSelObj3__12->SetBinContent(23,0.006292564);
   dPhiLepMETSelObj3__12->SetBinContent(24,0.003105841);
   dPhiLepMETSelObj3__12->SetBinContent(25,0.001496304);
   dPhiLepMETSelObj3__12->SetBinContent(26,0.001027192);
   dPhiLepMETSelObj3__12->SetBinContent(27,0.001043369);
   dPhiLepMETSelObj3__12->SetBinContent(28,0.001011016);
   dPhiLepMETSelObj3__12->SetBinContent(29,0.0009301347);
   dPhiLepMETSelObj3__12->SetBinContent(30,0.0009382229);
   dPhiLepMETSelObj3__12->SetBinContent(31,0.0006874909);
   dPhiLepMETSelObj3__12->SetBinContent(32,0.0007036672);
   dPhiLepMETSelObj3__12->SetBinContent(33,0.0008007247);
   dPhiLepMETSelObj3__12->SetBinContent(34,0.000630874);
   dPhiLepMETSelObj3__12->SetBinContent(35,0.0008007247);
   dPhiLepMETSelObj3__12->SetBinContent(36,0.0006146978);
   dPhiLepMETSelObj3__12->SetBinContent(37,0.0008007247);
   dPhiLepMETSelObj3__12->SetBinContent(38,0.0007602841);
   dPhiLepMETSelObj3__12->SetBinContent(39,0.0006794028);
   dPhiLepMETSelObj3__12->SetBinContent(40,0.0006632265);
   dPhiLepMETSelObj3__12->SetBinContent(41,0.0006470503);
   dPhiLepMETSelObj3__12->SetBinContent(42,0.0006551384);
   dPhiLepMETSelObj3__12->SetBinContent(43,0.0006632265);
   dPhiLepMETSelObj3__12->SetBinContent(44,0.0006874909);
   dPhiLepMETSelObj3__12->SetBinContent(45,0.0005338165);
   dPhiLepMETSelObj3__12->SetBinContent(46,0.0006146978);
   dPhiLepMETSelObj3__12->SetBinContent(47,0.0005419046);
   dPhiLepMETSelObj3__12->SetBinContent(48,0.0007360197);
   dPhiLepMETSelObj3__12->SetBinContent(49,0.0006066096);
   dPhiLepMETSelObj3__12->SetBinContent(50,0.0006227859);
   dPhiLepMETSelObj3__12->SetBinError(1,0.0006035283);
   dPhiLepMETSelObj3__12->SetBinError(2,0.0005888246);
   dPhiLepMETSelObj3__12->SetBinError(3,0.0006001042);
   dPhiLepMETSelObj3__12->SetBinError(4,0.0006048817);
   dPhiLepMETSelObj3__12->SetBinError(5,0.0006235732);
   dPhiLepMETSelObj3__12->SetBinError(6,0.0006379375);
   dPhiLepMETSelObj3__12->SetBinError(7,0.0006548887);
   dPhiLepMETSelObj3__12->SetBinError(8,0.0006770879);
   dPhiLepMETSelObj3__12->SetBinError(9,0.0006779087);
   dPhiLepMETSelObj3__12->SetBinError(10,0.0006941198);
   dPhiLepMETSelObj3__12->SetBinError(11,0.0006954379);
   dPhiLepMETSelObj3__12->SetBinError(12,0.0006952498);
   dPhiLepMETSelObj3__12->SetBinError(13,0.0006762178);
   dPhiLepMETSelObj3__12->SetBinError(14,0.0006550885);
   dPhiLepMETSelObj3__12->SetBinError(15,0.0006353172);
   dPhiLepMETSelObj3__12->SetBinError(16,0.0005986854);
   dPhiLepMETSelObj3__12->SetBinError(17,0.0005601284);
   dPhiLepMETSelObj3__12->SetBinError(18,0.0005168181);
   dPhiLepMETSelObj3__12->SetBinError(19,0.0004691114);
   dPhiLepMETSelObj3__12->SetBinError(20,0.0004117803);
   dPhiLepMETSelObj3__12->SetBinError(21,0.000345998);
   dPhiLepMETSelObj3__12->SetBinError(22,0.0002899343);
   dPhiLepMETSelObj3__12->SetBinError(23,0.0002255993);
   dPhiLepMETSelObj3__12->SetBinError(24,0.0001584943);
   dPhiLepMETSelObj3__12->SetBinError(25,0.0001100104);
   dPhiLepMETSelObj3__12->SetBinError(26,9.114858e-05);
   dPhiLepMETSelObj3__12->SetBinError(27,9.186348e-05);
   dPhiLepMETSelObj3__12->SetBinError(28,9.042802e-05);
   dPhiLepMETSelObj3__12->SetBinError(29,8.673551e-05);
   dPhiLepMETSelObj3__12->SetBinError(30,8.711181e-05);
   dPhiLepMETSelObj3__12->SetBinError(31,7.456886e-05);
   dPhiLepMETSelObj3__12->SetBinError(32,7.544104e-05);
   dPhiLepMETSelObj3__12->SetBinError(33,8.047586e-05);
   dPhiLepMETSelObj3__12->SetBinError(34,7.143241e-05);
   dPhiLepMETSelObj3__12->SetBinError(35,8.047586e-05);
   dPhiLepMETSelObj3__12->SetBinError(36,7.051067e-05);
   dPhiLepMETSelObj3__12->SetBinError(37,8.047586e-05);
   dPhiLepMETSelObj3__12->SetBinError(38,7.841731e-05);
   dPhiLepMETSelObj3__12->SetBinError(39,7.412892e-05);
   dPhiLepMETSelObj3__12->SetBinError(40,7.324112e-05);
   dPhiLepMETSelObj3__12->SetBinError(41,7.234242e-05);
   dPhiLepMETSelObj3__12->SetBinError(42,7.279315e-05);
   dPhiLepMETSelObj3__12->SetBinError(43,7.324112e-05);
   dPhiLepMETSelObj3__12->SetBinError(44,7.456886e-05);
   dPhiLepMETSelObj3__12->SetBinError(45,6.570826e-05);
   dPhiLepMETSelObj3__12->SetBinError(46,7.051067e-05);
   dPhiLepMETSelObj3__12->SetBinError(47,6.620418e-05);
   dPhiLepMETSelObj3__12->SetBinError(48,7.715583e-05);
   dPhiLepMETSelObj3__12->SetBinError(49,7.004525e-05);
   dPhiLepMETSelObj3__12->SetBinError(50,7.097304e-05);
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
