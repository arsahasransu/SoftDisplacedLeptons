void Phi_El()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed May 20 11:43:05 2020) by ROOT version 6.18/04
   TCanvas *c1 = new TCanvas("c1", "c1",10,32,782,600);
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
   
   TH1F *Phi_El0__29 = new TH1F("Phi_El0__29","",51,-3.2,3.2);
   Phi_El0__29->SetBinContent(1,0.008519961);
   Phi_El0__29->SetBinContent(2,0.02142162);
   Phi_El0__29->SetBinContent(3,0.01996105);
   Phi_El0__29->SetBinContent(4,0.02312561);
   Phi_El0__29->SetBinContent(5,0.02044791);
   Phi_El0__29->SetBinContent(6,0.02044791);
   Phi_El0__29->SetBinContent(7,0.01898734);
   Phi_El0__29->SetBinContent(8,0.02020448);
   Phi_El0__29->SetBinContent(9,0.02117819);
   Phi_El0__29->SetBinContent(10,0.02190847);
   Phi_El0__29->SetBinContent(11,0.02361246);
   Phi_El0__29->SetBinContent(12,0.01630964);
   Phi_El0__29->SetBinContent(13,0.01971762);
   Phi_El0__29->SetBinContent(14,0.0248296);
   Phi_El0__29->SetBinContent(15,0.02020448);
   Phi_El0__29->SetBinContent(16,0.02044791);
   Phi_El0__29->SetBinContent(17,0.01898734);
   Phi_El0__29->SetBinContent(18,0.02312561);
   Phi_El0__29->SetBinContent(19,0.02312561);
   Phi_El0__29->SetBinContent(20,0.01825706);
   Phi_El0__29->SetBinContent(21,0.01923077);
   Phi_El0__29->SetBinContent(22,0.0194742);
   Phi_El0__29->SetBinContent(23,0.02044791);
   Phi_El0__29->SetBinContent(24,0.01996105);
   Phi_El0__29->SetBinContent(25,0.01898734);
   Phi_El0__29->SetBinContent(26,0.01825706);
   Phi_El0__29->SetBinContent(27,0.02093476);
   Phi_El0__29->SetBinContent(28,0.01484908);
   Phi_El0__29->SetBinContent(29,0.0194742);
   Phi_El0__29->SetBinContent(30,0.01923077);
   Phi_El0__29->SetBinContent(31,0.0194742);
   Phi_El0__29->SetBinContent(32,0.02020448);
   Phi_El0__29->SetBinContent(33,0.02288218);
   Phi_El0__29->SetBinContent(34,0.0221519);
   Phi_El0__29->SetBinContent(35,0.02093476);
   Phi_El0__29->SetBinContent(36,0.01679649);
   Phi_El0__29->SetBinContent(37,0.01898734);
   Phi_El0__29->SetBinContent(38,0.02093476);
   Phi_El0__29->SetBinContent(39,0.01971762);
   Phi_El0__29->SetBinContent(40,0.02288218);
   Phi_El0__29->SetBinContent(41,0.01898734);
   Phi_El0__29->SetBinContent(42,0.01582278);
   Phi_El0__29->SetBinContent(43,0.02239533);
   Phi_El0__29->SetBinContent(44,0.02190847);
   Phi_El0__29->SetBinContent(45,0.02117819);
   Phi_El0__29->SetBinContent(46,0.01436222);
   Phi_El0__29->SetBinContent(47,0.01655307);
   Phi_El0__29->SetBinContent(48,0.01971762);
   Phi_El0__29->SetBinContent(49,0.01630964);
   Phi_El0__29->SetBinContent(50,0.02117819);
   Phi_El0__29->SetBinContent(51,0.01095424);
   Phi_El0__29->SetBinError(1,0.001440136);
   Phi_El0__29->SetBinError(2,0.002283552);
   Phi_El0__29->SetBinError(3,0.002204329);
   Phi_El0__29->SetBinError(4,0.002372637);
   Phi_El0__29->SetBinError(5,0.00223105);
   Phi_El0__29->SetBinError(6,0.00223105);
   Phi_El0__29->SetBinError(7,0.002149893);
   Phi_El0__29->SetBinError(8,0.00221773);
   Phi_El0__29->SetBinError(9,0.00227054);
   Phi_El0__29->SetBinError(10,0.002309356);
   Phi_El0__29->SetBinError(11,0.002397482);
   Phi_El0__29->SetBinError(12,0.00199254);
   Phi_El0__29->SetBinError(13,0.002190847);
   Phi_El0__29->SetBinError(14,0.002458497);
   Phi_El0__29->SetBinError(15,0.00221773);
   Phi_El0__29->SetBinError(16,0.00223105);
   Phi_El0__29->SetBinError(17,0.002149893);
   Phi_El0__29->SetBinError(18,0.002372637);
   Phi_El0__29->SetBinError(19,0.002372637);
   Phi_El0__29->SetBinError(20,0.002108144);
   Phi_El0__29->SetBinError(21,0.002163631);
   Phi_El0__29->SetBinError(22,0.002177281);
   Phi_El0__29->SetBinError(23,0.00223105);
   Phi_El0__29->SetBinError(24,0.002204329);
   Phi_El0__29->SetBinError(25,0.002149893);
   Phi_El0__29->SetBinError(26,0.002108144);
   Phi_El0__29->SetBinError(27,0.002257453);
   Phi_El0__29->SetBinError(28,0.001901229);
   Phi_El0__29->SetBinError(29,0.002177281);
   Phi_El0__29->SetBinError(30,0.002163631);
   Phi_El0__29->SetBinError(31,0.002177281);
   Phi_El0__29->SetBinError(32,0.00221773);
   Phi_El0__29->SetBinError(33,0.002360117);
   Phi_El0__29->SetBinError(34,0.00232215);
   Phi_El0__29->SetBinError(35,0.002257453);
   Phi_El0__29->SetBinError(36,0.00202206);
   Phi_El0__29->SetBinError(37,0.002149893);
   Phi_El0__29->SetBinError(38,0.002257453);
   Phi_El0__29->SetBinError(39,0.002190847);
   Phi_El0__29->SetBinError(40,0.002360117);
   Phi_El0__29->SetBinError(41,0.002149893);
   Phi_El0__29->SetBinError(42,0.001962575);
   Phi_El0__29->SetBinError(43,0.002334874);
   Phi_El0__29->SetBinError(44,0.002309356);
   Phi_El0__29->SetBinError(45,0.00227054);
   Phi_El0__29->SetBinError(46,0.001869802);
   Phi_El0__29->SetBinError(47,0.002007354);
   Phi_El0__29->SetBinError(48,0.002190847);
   Phi_El0__29->SetBinError(49,0.00199254);
   Phi_El0__29->SetBinError(50,0.00227054);
   Phi_El0__29->SetBinError(51,0.001632961);
   Phi_El0__29->SetMinimum(0.002);
   Phi_El0__29->SetMaximum(0.0496592);
   Phi_El0__29->SetEntries(4108);
   Phi_El0__29->SetLineWidth(3);
   Phi_El0__29->GetXaxis()->SetLabelFont(42);
   Phi_El0__29->GetXaxis()->SetLabelOffset(-0.07);
   Phi_El0__29->GetXaxis()->SetLabelSize(0.055);
   Phi_El0__29->GetXaxis()->SetTitleSize(0.035);
   Phi_El0__29->GetXaxis()->SetTitleOffset(1);
   Phi_El0__29->GetXaxis()->SetTitleFont(42);
   Phi_El0__29->GetYaxis()->SetLabelFont(42);
   Phi_El0__29->GetYaxis()->SetLabelOffset(-0.05);
   Phi_El0__29->GetYaxis()->SetLabelSize(0.055);
   Phi_El0__29->GetYaxis()->SetTitleSize(0.035);
   Phi_El0__29->GetYaxis()->SetTitleFont(42);
   Phi_El0__29->GetZaxis()->SetLabelFont(42);
   Phi_El0__29->GetZaxis()->SetLabelSize(0.035);
   Phi_El0__29->GetZaxis()->SetTitleSize(0.035);
   Phi_El0__29->GetZaxis()->SetTitleOffset(1);
   Phi_El0__29->GetZaxis()->SetTitleFont(42);
   Phi_El0__29->Draw("hist  E");
   
   TH1F *Phi_El1__30 = new TH1F("Phi_El1__30","",51,-3.2,3.2);
   Phi_El1__30->SetBinContent(1,0.01123713);
   Phi_El1__30->SetBinContent(2,0.0186066);
   Phi_El1__30->SetBinContent(3,0.01965191);
   Phi_El1__30->SetBinContent(4,0.02080176);
   Phi_El1__30->SetBinContent(5,0.02054043);
   Phi_El1__30->SetBinContent(6,0.01980871);
   Phi_El1__30->SetBinContent(7,0.01818847);
   Phi_El1__30->SetBinContent(8,0.02163801);
   Phi_El1__30->SetBinContent(9,0.01865886);
   Phi_El1__30->SetBinContent(10,0.02033136);
   Phi_El1__30->SetBinContent(11,0.02007004);
   Phi_El1__30->SetBinContent(12,0.01970418);
   Phi_El1__30->SetBinContent(13,0.01850206);
   Phi_El1__30->SetBinContent(14,0.02080176);
   Phi_El1__30->SetBinContent(15,0.01980871);
   Phi_El1__30->SetBinContent(16,0.01944285);
   Phi_El1__30->SetBinContent(17,0.01881566);
   Phi_El1__30->SetBinContent(18,0.02048816);
   Phi_El1__30->SetBinContent(19,0.02054043);
   Phi_El1__30->SetBinContent(20,0.02048816);
   Phi_El1__30->SetBinContent(21,0.01928605);
   Phi_El1__30->SetBinContent(22,0.018293);
   Phi_El1__30->SetBinContent(23,0.01902472);
   Phi_El1__30->SetBinContent(24,0.02069722);
   Phi_El1__30->SetBinContent(25,0.01886792);
   Phi_El1__30->SetBinContent(26,0.01944285);
   Phi_El1__30->SetBinContent(27,0.02095855);
   Phi_El1__30->SetBinContent(28,0.02048816);
   Phi_El1__30->SetBinContent(29,0.02121988);
   Phi_El1__30->SetBinContent(30,0.01750902);
   Phi_El1__30->SetBinContent(31,0.01771808);
   Phi_El1__30->SetBinContent(32,0.01986097);
   Phi_El1__30->SetBinContent(33,0.02054043);
   Phi_El1__30->SetBinContent(34,0.0204359);
   Phi_El1__30->SetBinContent(35,0.01902472);
   Phi_El1__30->SetBinContent(36,0.02069722);
   Phi_El1__30->SetBinContent(37,0.01928605);
   Phi_El1__30->SetBinContent(38,0.02242199);
   Phi_El1__30->SetBinContent(39,0.02054043);
   Phi_El1__30->SetBinContent(40,0.02033136);
   Phi_El1__30->SetBinContent(41,0.01975644);
   Phi_El1__30->SetBinContent(42,0.01912925);
   Phi_El1__30->SetBinContent(43,0.02121988);
   Phi_El1__30->SetBinContent(44,0.02038363);
   Phi_El1__30->SetBinContent(45,0.01986097);
   Phi_El1__30->SetBinContent(46,0.01991324);
   Phi_El1__30->SetBinContent(47,0.02069722);
   Phi_El1__30->SetBinContent(48,0.02200387);
   Phi_El1__30->SetBinContent(49,0.01996551);
   Phi_El1__30->SetBinContent(50,0.02022683);
   Phi_El1__30->SetBinContent(51,0.01207338);
   Phi_El1__30->SetBinError(1,0.0007663659);
   Phi_El1__30->SetBinError(2,0.0009861476);
   Phi_El1__30->SetBinError(3,0.00101347);
   Phi_El1__30->SetBinError(4,0.001042698);
   Phi_El1__30->SetBinError(5,0.001036128);
   Phi_El1__30->SetBinError(6,0.001017505);
   Phi_El1__30->SetBinError(7,0.0009750043);
   Phi_El1__30->SetBinError(8,0.00106345);
   Phi_El1__30->SetBinError(9,0.0009875317);
   Phi_El1__30->SetBinError(10,0.001030841);
   Phi_El1__30->SetBinError(11,0.001024195);
   Phi_El1__30->SetBinError(12,0.001014817);
   Phi_El1__30->SetBinError(13,0.0009833736);
   Phi_El1__30->SetBinError(14,0.001042698);
   Phi_El1__30->SetBinError(15,0.001017505);
   Phi_El1__30->SetBinError(16,0.001008065);
   Phi_El1__30->SetBinError(17,0.0009916723);
   Phi_El1__30->SetBinError(18,0.001034808);
   Phi_El1__30->SetBinError(19,0.001036128);
   Phi_El1__30->SetBinError(20,0.001034808);
   Phi_El1__30->SetBinError(21,0.001003992);
   Phi_El1__30->SetBinError(22,0.0009778021);
   Phi_El1__30->SetBinError(23,0.0009971664);
   Phi_El1__30->SetBinError(24,0.001040075);
   Phi_El1__30->SetBinError(25,0.0009930487);
   Phi_El1__30->SetBinError(26,0.001008065);
   Phi_El1__30->SetBinError(27,0.00104662);
   Phi_El1__30->SetBinError(28,0.001034808);
   Phi_El1__30->SetBinError(29,0.001053125);
   Phi_El1__30->SetBinError(30,0.0009566197);
   Phi_El1__30->SetBinError(31,0.0009623139);
   Phi_El1__30->SetBinError(32,0.001018846);
   Phi_El1__30->SetBinError(33,0.001036128);
   Phi_El1__30->SetBinError(34,0.001033488);
   Phi_El1__30->SetBinError(35,0.0009971664);
   Phi_El1__30->SetBinError(36,0.001040075);
   Phi_El1__30->SetBinError(37,0.001003992);
   Phi_El1__30->SetBinError(38,0.001082544);
   Phi_El1__30->SetBinError(39,0.001036128);
   Phi_El1__30->SetBinError(40,0.001030841);
   Phi_El1__30->SetBinError(41,0.001016162);
   Phi_El1__30->SetBinError(42,0.0009999021);
   Phi_El1__30->SetBinError(43,0.001053125);
   Phi_El1__30->SetBinError(44,0.001032165);
   Phi_El1__30->SetBinError(45,0.001018846);
   Phi_El1__30->SetBinError(46,0.001020186);
   Phi_El1__30->SetBinError(47,0.001040075);
   Phi_El1__30->SetBinError(48,0.001072403);
   Phi_El1__30->SetBinError(49,0.001021524);
   Phi_El1__30->SetBinError(50,0.001028188);
   Phi_El1__30->SetBinError(51,0.0007943702);
   Phi_El1__30->SetMinimum(0.002);
   Phi_El1__30->SetMaximum(0.04484399);
   Phi_El1__30->SetEntries(19133);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#0000ff");
   Phi_El1__30->SetLineColor(ci);
   Phi_El1__30->SetLineWidth(3);
   Phi_El1__30->GetXaxis()->SetLabelFont(42);
   Phi_El1__30->GetXaxis()->SetLabelSize(0.035);
   Phi_El1__30->GetXaxis()->SetTitleSize(0.035);
   Phi_El1__30->GetXaxis()->SetTitleOffset(1);
   Phi_El1__30->GetXaxis()->SetTitleFont(42);
   Phi_El1__30->GetYaxis()->SetLabelFont(42);
   Phi_El1__30->GetYaxis()->SetLabelSize(0.035);
   Phi_El1__30->GetYaxis()->SetTitleSize(0.035);
   Phi_El1__30->GetYaxis()->SetTitleFont(42);
   Phi_El1__30->GetZaxis()->SetLabelFont(42);
   Phi_El1__30->GetZaxis()->SetLabelSize(0.035);
   Phi_El1__30->GetZaxis()->SetTitleSize(0.035);
   Phi_El1__30->GetZaxis()->SetTitleOffset(1);
   Phi_El1__30->GetZaxis()->SetTitleFont(42);
   Phi_El1__30->Draw("hist SAME E");
   
   TH1F *Phi_El2__31 = new TH1F("Phi_El2__31","",51,-3.2,3.2);
   Phi_El2__31->SetBinContent(1,0.009861602);
   Phi_El2__31->SetBinContent(2,0.01994733);
   Phi_El2__31->SetBinContent(3,0.02067574);
   Phi_El2__31->SetBinContent(4,0.02089987);
   Phi_El2__31->SetBinContent(5,0.01944304);
   Phi_El2__31->SetBinContent(6,0.02056368);
   Phi_El2__31->SetBinContent(7,0.02140416);
   Phi_El2__31->SetBinContent(8,0.02202051);
   Phi_El2__31->SetBinContent(9,0.01977924);
   Phi_El2__31->SetBinContent(10,0.02045162);
   Phi_El2__31->SetBinContent(11,0.01916288);
   Phi_El2__31->SetBinContent(12,0.02045162);
   Phi_El2__31->SetBinContent(13,0.02157225);
   Phi_El2__31->SetBinContent(14,0.01826638);
   Phi_El2__31->SetBinContent(15,0.01916288);
   Phi_El2__31->SetBinContent(16,0.01955511);
   Phi_El2__31->SetBinContent(17,0.01882669);
   Phi_El2__31->SetBinContent(18,0.02022749);
   Phi_El2__31->SetBinContent(19,0.01921892);
   Phi_El2__31->SetBinContent(20,0.01916288);
   Phi_El2__31->SetBinContent(21,0.01961114);
   Phi_El2__31->SetBinContent(22,0.01888273);
   Phi_El2__31->SetBinContent(23,0.01961114);
   Phi_El2__31->SetBinContent(24,0.01899479);
   Phi_El2__31->SetBinContent(25,0.01837844);
   Phi_El2__31->SetBinContent(26,0.0184905);
   Phi_El2__31->SetBinContent(27,0.02061971);
   Phi_El2__31->SetBinContent(28,0.02269289);
   Phi_El2__31->SetBinContent(29,0.01882669);
   Phi_El2__31->SetBinContent(30,0.02207654);
   Phi_El2__31->SetBinContent(31,0.02106797);
   Phi_El2__31->SetBinContent(32,0.02146019);
   Phi_El2__31->SetBinContent(33,0.02005939);
   Phi_El2__31->SetBinContent(34,0.01905082);
   Phi_El2__31->SetBinContent(35,0.0209559);
   Phi_El2__31->SetBinContent(36,0.01821034);
   Phi_El2__31->SetBinContent(37,0.02050765);
   Phi_El2__31->SetBinContent(38,0.01916288);
   Phi_El2__31->SetBinContent(39,0.01837844);
   Phi_El2__31->SetBinContent(40,0.01860257);
   Phi_El2__31->SetBinContent(41,0.02241273);
   Phi_El2__31->SetBinContent(42,0.02022749);
   Phi_El2__31->SetBinContent(43,0.02056368);
   Phi_El2__31->SetBinContent(44,0.01927495);
   Phi_El2__31->SetBinContent(45,0.02022749);
   Phi_El2__31->SetBinContent(46,0.02162828);
   Phi_El2__31->SetBinContent(47,0.01882669);
   Phi_El2__31->SetBinContent(48,0.01921892);
   Phi_El2__31->SetBinContent(49,0.01977924);
   Phi_El2__31->SetBinContent(50,0.01933098);
   Phi_El2__31->SetBinContent(51,0.01221494);
   Phi_El2__31->SetBinError(1,0.0007433462);
   Phi_El2__31->SetBinError(2,0.001057206);
   Phi_El2__31->SetBinError(3,0.001076336);
   Phi_El2__31->SetBinError(4,0.001082154);
   Phi_El2__31->SetBinError(5,0.001043757);
   Phi_El2__31->SetBinError(6,0.001073415);
   Phi_El2__31->SetBinError(7,0.001095132);
   Phi_El2__31->SetBinError(8,0.001110788);
   Phi_El2__31->SetBinError(9,0.001052742);
   Phi_El2__31->SetBinError(10,0.001070487);
   Phi_El2__31->SetBinError(11,0.00103621);
   Phi_El2__31->SetBinError(12,0.001070487);
   Phi_El2__31->SetBinError(13,0.001099424);
   Phi_El2__31->SetBinError(14,0.001011681);
   Phi_El2__31->SetBinError(15,0.00103621);
   Phi_El2__31->SetBinError(16,0.001046761);
   Phi_El2__31->SetBinError(17,0.00102708);
   Phi_El2__31->SetBinError(18,0.001064605);
   Phi_El2__31->SetBinError(19,0.001037724);
   Phi_El2__31->SetBinError(20,0.00103621);
   Phi_El2__31->SetBinError(21,0.001048259);
   Phi_El2__31->SetBinError(22,0.001028608);
   Phi_El2__31->SetBinError(23,0.001048259);
   Phi_El2__31->SetBinError(24,0.001031655);
   Phi_El2__31->SetBinError(25,0.00101478);
   Phi_El2__31->SetBinError(26,0.001017869);
   Phi_El2__31->SetBinError(27,0.001074877);
   Phi_El2__31->SetBinError(28,0.001127619);
   Phi_El2__31->SetBinError(29,0.00102708);
   Phi_El2__31->SetBinError(30,0.0011122);
   Phi_El2__31->SetBinError(31,0.001086497);
   Phi_El2__31->SetBinError(32,0.001096564);
   Phi_El2__31->SetBinError(33,0.001060172);
   Phi_El2__31->SetBinError(34,0.001033176);
   Phi_El2__31->SetBinError(35,0.001083604);
   Phi_El2__31->SetBinError(36,0.001010128);
   Phi_El2__31->SetBinError(37,0.001071952);
   Phi_El2__31->SetBinError(38,0.00103621);
   Phi_El2__31->SetBinError(39,0.00101478);
   Phi_El2__31->SetBinError(40,0.001020948);
   Phi_El2__31->SetBinError(41,0.001120637);
   Phi_El2__31->SetBinError(42,0.001064605);
   Phi_El2__31->SetBinError(43,0.001073415);
   Phi_El2__31->SetBinError(44,0.001039236);
   Phi_El2__31->SetBinError(45,0.001064605);
   Phi_El2__31->SetBinError(46,0.001100851);
   Phi_El2__31->SetBinError(47,0.00102708);
   Phi_El2__31->SetBinError(48,0.001037724);
   Phi_El2__31->SetBinError(49,0.001052742);
   Phi_El2__31->SetBinError(50,0.001040745);
   Phi_El2__31->SetBinError(51,0.0008273);
   Phi_El2__31->SetMinimum(0.002);
   Phi_El2__31->SetMaximum(0.04538578);
   Phi_El2__31->SetEntries(17847);

   ci = TColor::GetColor("#6666ff");
   Phi_El2__31->SetLineColor(ci);
   Phi_El2__31->SetLineWidth(3);
   Phi_El2__31->GetXaxis()->SetLabelFont(42);
   Phi_El2__31->GetXaxis()->SetLabelSize(0.035);
   Phi_El2__31->GetXaxis()->SetTitleSize(0.035);
   Phi_El2__31->GetXaxis()->SetTitleOffset(1);
   Phi_El2__31->GetXaxis()->SetTitleFont(42);
   Phi_El2__31->GetYaxis()->SetLabelFont(42);
   Phi_El2__31->GetYaxis()->SetLabelSize(0.035);
   Phi_El2__31->GetYaxis()->SetTitleSize(0.035);
   Phi_El2__31->GetYaxis()->SetTitleFont(42);
   Phi_El2__31->GetZaxis()->SetLabelFont(42);
   Phi_El2__31->GetZaxis()->SetLabelSize(0.035);
   Phi_El2__31->GetZaxis()->SetTitleSize(0.035);
   Phi_El2__31->GetZaxis()->SetTitleOffset(1);
   Phi_El2__31->GetZaxis()->SetTitleFont(42);
   Phi_El2__31->Draw("hist SAME E");
   
   TH1F *Phi_El3__32 = new TH1F("Phi_El3__32","",51,-3.2,3.2);
   Phi_El3__32->SetBinContent(1,0.01047413);
   Phi_El3__32->SetBinContent(2,0.02012326);
   Phi_El3__32->SetBinContent(3,0.0197593);
   Phi_El3__32->SetBinContent(4,0.01955709);
   Phi_El3__32->SetBinContent(5,0.01962989);
   Phi_El3__32->SetBinContent(6,0.01983209);
   Phi_El3__32->SetBinContent(7,0.01996959);
   Phi_El3__32->SetBinContent(8,0.01959753);
   Phi_El3__32->SetBinContent(9,0.01984827);
   Phi_El3__32->SetBinContent(10,0.02027694);
   Phi_El3__32->SetBinContent(11,0.0196865);
   Phi_El3__32->SetBinContent(12,0.02027694);
   Phi_El3__32->SetBinContent(13,0.02029311);
   Phi_El3__32->SetBinContent(14,0.02007473);
   Phi_El3__32->SetBinContent(15,0.01984018);
   Phi_El3__32->SetBinContent(16,0.02079458);
   Phi_El3__32->SetBinContent(17,0.02034164);
   Phi_El3__32->SetBinContent(18,0.02046296);
   Phi_El3__32->SetBinContent(19,0.02013944);
   Phi_El3__32->SetBinContent(20,0.01988062);
   Phi_El3__32->SetBinContent(21,0.01950856);
   Phi_El3__32->SetBinContent(22,0.02018797);
   Phi_El3__32->SetBinContent(23,0.02013135);
   Phi_El3__32->SetBinContent(24,0.01950856);
   Phi_El3__32->SetBinContent(25,0.02068943);
   Phi_El3__32->SetBinContent(26,0.02001812);
   Phi_El3__32->SetBinContent(27,0.02026076);
   Phi_El3__32->SetBinContent(28,0.01994532);
   Phi_El3__32->SetBinContent(29,0.01970268);
   Phi_El3__32->SetBinContent(30,0.01946812);
   Phi_El3__32->SetBinContent(31,0.02027694);
   Phi_El3__32->SetBinContent(32,0.01966224);
   Phi_El3__32->SetBinContent(33,0.020099);
   Phi_El3__32->SetBinContent(34,0.01933063);
   Phi_El3__32->SetBinContent(35,0.02010709);
   Phi_El3__32->SetBinContent(36,0.02007473);
   Phi_El3__32->SetBinContent(37,0.02006665);
   Phi_El3__32->SetBinContent(38,0.02061664);
   Phi_El3__32->SetBinContent(39,0.01929827);
   Phi_El3__32->SetBinContent(40,0.0204387);
   Phi_El3__32->SetBinContent(41,0.01974312);
   Phi_El3__32->SetBinContent(42,0.0197593);
   Phi_El3__32->SetBinContent(43,0.019824);
   Phi_El3__32->SetBinContent(44,0.01931445);
   Phi_El3__32->SetBinContent(45,0.02008282);
   Phi_El3__32->SetBinContent(46,0.0201637);
   Phi_El3__32->SetBinContent(47,0.0197593);
   Phi_El3__32->SetBinContent(48,0.02017179);
   Phi_El3__32->SetBinContent(49,0.02033355);
   Phi_El3__32->SetBinContent(50,0.019824);
   Phi_El3__32->SetBinContent(51,0.01077339);
   Phi_El3__32->SetBinError(1,0.0002910603);
   Phi_El3__32->SetBinError(2,0.0004034347);
   Phi_El3__32->SetBinError(3,0.0003997696);
   Phi_El3__32->SetBinError(4,0.0003977189);
   Phi_El3__32->SetBinError(5,0.0003984583);
   Phi_El3__32->SetBinError(6,0.0004005053);
   Phi_El3__32->SetBinError(7,0.0004018913);
   Phi_El3__32->SetBinError(8,0.0003981298);
   Phi_El3__32->SetBinError(9,0.0004006686);
   Phi_El3__32->SetBinError(10,0.0004049722);
   Phi_El3__32->SetBinError(11,0.0003990325);
   Phi_El3__32->SetBinError(12,0.0004049722);
   Phi_El3__32->SetBinError(13,0.0004051337);
   Phi_El3__32->SetBinError(14,0.0004029479);
   Phi_El3__32->SetBinError(15,0.000400587);
   Phi_El3__32->SetBinError(16,0.0004101088);
   Phi_El3__32->SetBinError(17,0.0004056178);
   Phi_El3__32->SetBinError(18,0.0004068256);
   Phi_El3__32->SetBinError(19,0.0004035968);
   Phi_El3__32->SetBinError(20,0.000400995);
   Phi_El3__32->SetBinError(21,0.0003972251);
   Phi_El3__32->SetBinError(22,0.0004040828);
   Phi_El3__32->SetBinError(23,0.0004035157);
   Phi_El3__32->SetBinError(24,0.0003972251);
   Phi_El3__32->SetBinError(25,0.0004090706);
   Phi_El3__32->SetBinError(26,0.0004023793);
   Phi_El3__32->SetBinError(27,0.0004048106);
   Phi_El3__32->SetBinError(28,0.000401647);
   Phi_El3__32->SetBinError(29,0.0003991965);
   Phi_El3__32->SetBinError(30,0.0003968132);
   Phi_El3__32->SetBinError(31,0.0004049722);
   Phi_El3__32->SetBinError(32,0.0003987866);
   Phi_El3__32->SetBinError(33,0.0004031914);
   Phi_El3__32->SetBinError(34,0.0003954094);
   Phi_El3__32->SetBinError(35,0.0004032725);
   Phi_El3__32->SetBinError(36,0.0004029479);
   Phi_El3__32->SetBinError(37,0.0004028667);
   Phi_El3__32->SetBinError(38,0.0004083504);
   Phi_El3__32->SetBinError(39,0.0003950784);
   Phi_El3__32->SetBinError(40,0.0004065843);
   Phi_El3__32->SetBinError(41,0.0003996059);
   Phi_El3__32->SetBinError(42,0.0003997696);
   Phi_El3__32->SetBinError(43,0.0004004236);
   Phi_El3__32->SetBinError(44,0.0003952439);
   Phi_El3__32->SetBinError(45,0.0004030291);
   Phi_El3__32->SetBinError(46,0.0004038398);
   Phi_El3__32->SetBinError(47,0.0003997696);
   Phi_El3__32->SetBinError(48,0.0004039208);
   Phi_El3__32->SetBinError(49,0.0004055372);
   Phi_El3__32->SetBinError(50,0.0004004236);
   Phi_El3__32->SetBinError(51,0.000295189);
   Phi_El3__32->SetMinimum(0.002);
   Phi_El3__32->SetMaximum(0.04158916);
   Phi_El3__32->SetEntries(123638);

   ci = TColor::GetColor("#cc3399");
   Phi_El3__32->SetLineColor(ci);
   Phi_El3__32->SetLineWidth(3);
   Phi_El3__32->GetXaxis()->SetLabelFont(42);
   Phi_El3__32->GetXaxis()->SetLabelSize(0.035);
   Phi_El3__32->GetXaxis()->SetTitleSize(0.035);
   Phi_El3__32->GetXaxis()->SetTitleOffset(1);
   Phi_El3__32->GetXaxis()->SetTitleFont(42);
   Phi_El3__32->GetYaxis()->SetLabelFont(42);
   Phi_El3__32->GetYaxis()->SetLabelSize(0.035);
   Phi_El3__32->GetYaxis()->SetTitleSize(0.035);
   Phi_El3__32->GetYaxis()->SetTitleFont(42);
   Phi_El3__32->GetZaxis()->SetLabelFont(42);
   Phi_El3__32->GetZaxis()->SetLabelSize(0.035);
   Phi_El3__32->GetZaxis()->SetTitleSize(0.035);
   Phi_El3__32->GetZaxis()->SetTitleOffset(1);
   Phi_El3__32->GetZaxis()->SetTitleFont(42);
   Phi_El3__32->Draw("hist SAME E");
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
      tex = new TLatex(0.6,0.3,"electron #phi");
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
   TLegendEntry *entry=leg->AddEntry("Phi_El0","HF Background","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
      tex = new TLatex(0.12,0.15,"signal (m_{c} [GeV], #Deltam [GeV])");
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
   entry=leg->AddEntry("Phi_El1","DM: (220, 20)","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Phi_El2","DM: (324, 20)","l");
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
