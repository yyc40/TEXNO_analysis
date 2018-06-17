#include "TCutG.h"
#include "TCut.h"
#include "TFile.h"

void maincode()
{
 
  TFile *f = new TFile("/RAID_2/p93_rootfiles/160422_1252_phys_on_minus100_FADC_RAW_Data_1_1270.root");
  TTree *tr=(TTree*)f->Get("tr");
  TFormula::SetMaxima(1000000,1000000,1000000);

  TCut ran = Form("random_trig_on_off==1");
  TCut aran = Form("random_trig_on_off==0");

  //////////////Basic Cuts/////////////

 
  TCut os_cut="(max_60m[0]<32000) ";

  TCut min_cut = Form("min_60m[0]>-900&&min_60m[0]<-300");
  TCut ped_cut = Form("ped_60m[0]>-700&&ped_60m[0]<300");
  TCut pedt_cut = Form("pedt_60m[0]>-800&&pedt_60m[0]<250");
  TCut diff_pedt0_cut=" ped_60m[0]-pedt_60m[0] >-700&&ped_60m[0]-pedt_60m[0]<700";
 
  TCut min1_cut = Form("min_60m[1]>-850&&min_60m[1]<-200");
  TCut ped1_cut="ped_60m[1]>-700&& ped_60m[1]<100";
  TCut pedt1_cut="pedt_60m[1]>-750&& pedt_60m[1]<150";
  TCut diff_pedt1_cut=" ped_60m[1]-pedt_60m[1] >-550&&ped_60m[1]-pedt_60m[1]<550";
 
 
  TCut allcut0=os_cut&&min_cut&&ped_cut&&pedt_cut;
  TCut allcut1=min1_cut&&ped1_cut&&pedt1_cut;
  TCut ped_diff_cut= diff_pedt0_cut&&diff_pedt1_cut;
  TCut Basic_cut=allcut0&&allcut1&&ped_diff_cut;
 
  //////////////Linear Cuts/////////////
  TCut time_binch0_cut = Form("max_time_bin_60m[0]>1250&&(max_time_bin_60m[0]<(2150.0-(600.0)*exp(-8.0E-04*max_60m[0])))");
  TCut time_binch1_cut = Form("max_time_bin_60m[1]>1950&&(max_time_bin_60m[1]<(3050.0-(200.0)*exp(-4.0E-04*max_60m[1])))");

  TCut timebin_cut = time_binch0_cut&&time_binch1_cut;
  /////AQCut/////
  ///////upper level up part
  Double_t qax1 =2.32244e+07; 
  Double_t qay1 =33040.9;
  Double_t qax2 =-2.23747e+06;
  Double_t qay2 =-471.573;
  Double_t  slope_s1_ua = (qay2-qay1)/(qax2-qax1);
  Double_t  intercept_s1_ua = ((qay2-qay1)/(qax2-qax1))*(-1*qax1)+qay1;
  TCut qa_cutup1a= Form("max_60m[0]<(%f*q_60m[0]+%f)",slope_s1_ua,intercept_s1_ua);
  TCut qa_cutup1b=Form("max_60m[0]>700");
  TCut Uppera=qa_cutup1a&&qa_cutup1b;

  ///////upper level down part
  Double_t qa2x1 =-1.12804e+06;
  Double_t qa2y1= 906.425;
  Double_t qa2x2 =-1.90441e+06;
  Double_t qa2y2=-153.756; 
  Double_t  slope_s1_ub = (qa2y2-qa2y1)/(qa2x2-qa2x1);
  Double_t  intercept_s1_ub = ((qa2y2-qa2y1)/(qa2x2-qa2x1))*(-1*qa2x1)+qa2y1;
  TCut  qa_cutup2a= Form("max_60m[0]<(%f*q_60m[0]+%f)",slope_s1_ub,intercept_s1_ub);
  TCut  qa_cutup2b=Form("max_60m[0]<=700");
  TCut  Upperb=qa_cutup2a&& qa_cutup2b;
  /////////////******* Lower*********///////////////  
  Double_t qax3 =2.66978e+07;
  Double_t qay3 =32990.3;
  Double_t qax4 =-1.26148e+06;
  Double_t qay4 =-269.081;
  Double_t  slope_s1_la = (qay4-qay3)/(qax4-qax3);
  Double_t  intercept_s1_la = ((qay4-qay3)/(qax4-qax3))*(-1*qax3)+qay3;
  TCut qa_cut3a= Form("max_60m[0]>(%f*q_60m[0]+%f)",slope_s1_la,intercept_s1_la);
  TCut  qa_cut3b=Form("max_60m[0]>700");
  TCut Lowera=qa_cut3a&& qa_cut3b;
  ///////---------------------------------	   
  Double_t qbx3 =-366742;
  Double_t qby3 =786.968;
  Double_t qbx4 =-1.25618e+06;
  Double_t qby4 =-138.824;
  Double_t  slope_s1_lb = (qby4-qby3)/(qbx4-qbx3);
  Double_t  intercept_s1_lb = ((qby4-qby3)/(qbx4-qbx3))*(-1*qbx3)+qby3;
  TCut qa_cut4a= Form("max_60m[0]>(%f*q_60m[0]+%f)",slope_s1_lb,intercept_s1_lb);
  TCut qa_cut4b=Form("max_60m[0]<=700");
  TCut  Lowerb=qa_cut4a&&qa_cut4b;
  TCut qa_cut = (Uppera&&Lowera)||(Upperb&&Lowerb);

  ////////////S12 Cut /////////////
  Double_t sx1 =634.508;
  Double_t sy1 =402.055;
  Double_t sx2 =103.028;
  Double_t sy2=-275.805;
  Double_t  slope_s2 = (sy2-sy1)/(sx2-sx1);
  Double_t  intercept_s2 = ((sy2-sy1)/(sx2-sx1))*(-1*sx1)+sy1;
  TCut S12= Form("(max_60m[1]>(%f*max_60m[0]+%f)&&max_60m[0]<=600&&max_60m[0]>150)",slope_s2,intercept_s2);
  TCut s12_cut= (S12&&"max_60m[0]>150.0")||"max_60m[0]>600";  
   
  //////200M Cut///////
  
  Double_t Mx1 =243812;
  Double_t My1 =7097.24;
  Double_t Mx2 =-550630;
  Double_t My2=4529.73;
  Double_t  slopeM = (My2-My1)/(Mx2-Mx1);
  Double_t  interceptM = ((My2-My1)/(Mx2-Mx1))*(-1*Mx1)+My1;
  TCut ML;
  ML = Form("(max_200m1[0]<(%f*optq_60m[2]+%f)&&optq_60m[2]<=100000)||optq_60m[2]>100000",slopeM,interceptM);
 
  Double_t Mx3 =-569417;
  Double_t My3 =2315.11;
  Double_t Mx4 =177524;
  Double_t My4=2476.03;
  Double_t  slopeM1 = (My4-My3)/(Mx4-Mx3);
  Double_t  interceptM1 = ((My4-My3)/(Mx4-Mx3))*(-1*Mx3)+My3;
	
  TCut MU = Form("(max_200m1[0]>(%f*optq_60m[2]+%f)&&optq_60m[2]<=100000)||optq_60m[2]>100000",slopeM1,interceptM1);
  TCut thd_cut =(ML&&MU);
  //////optimal Q cut///////////////
  //////////-----------------
  Double_t uTRx1 =35641.6;
  Double_t  uTRy1 =702.166;
  Double_t uTRx2 =-508762;
  Double_t  uTRy2=145.013;
  Double_t  uTRslope = (uTRy2-uTRy1)/(uTRx2-uTRx1);
  Double_t  uTRintercept = ((uTRy2-uTRy1)/(uTRx2-uTRx1))*(-1*uTRx1)+uTRy1;
  TCut pq_cut_u;
  pq_cut_u= Form("(max_60m[0]<(%f*optq_60m[2]+%f)&&optq_60m[2]<=158500)||optq_60m[2]>158500",uTRslope,uTRintercept);
  Double_t lTRx1 =207875;
  Double_t  lTRy1 =322.493;
  Double_t lTRx2 =-61917.5;
  Double_t  lTRy2=185.451;
  Double_t  lTRslope = (lTRy2-lTRy1)/(lTRx2-lTRx1);
  Double_t  lTRintercept = ((lTRy2-lTRy1)/(lTRx2-lTRx1))*(-1*lTRx1)+lTRy1;
  TCut pq_cut_l,pq_cut;
  pq_cut_l= Form("(max_60m[0]>(%f*optq_60m[2]+%f)&&optq_60m[2]<=158500)||optq_60m[2]>158500",lTRslope,lTRintercept);
  pq_cut=pq_cut_u&&pq_cut_l;
  /////------------------------------------

  TCut Linear_cut = timebin_cut&&qa_cut&&s12_cut&&thd_cut&&pq_cut;
  /////////////ACV &&CRV Cut/////////////
  
  TCut acv_cut = Form("min_20m1[1]>=-2600");
  TCut act_cut = !acv_cut;

  //////////////////////crv cut
  ////////////upper veto
  TCut crv_tdc0_cut = Form("(veto_tdc[0]< 425 || veto_tdc[0] > 950)");
  TCut crv_tdc1_cut = Form("(veto_tdc[1]< 425 || veto_tdc[1] > 950)");
  TCut crv_tdc2_cut = Form("(veto_tdc[2]< 425 || veto_tdc[2] > 950)");
  TCut crv_tdc3_cut = Form("(veto_tdc[3]< 425 || veto_tdc[3] > 950)");
  TCut crv_tdc4_cut = Form("(veto_tdc[4]< 425 || veto_tdc[4] > 950)");
  TCut crv_tdc5_cut = Form("(veto_tdc[5]< 425 || veto_tdc[5] > 950)");
  TCut crv_tdc6_cut = Form("(veto_tdc[6]< 425 || veto_tdc[6] > 950)");
  TCut crv_tdc7_cut = Form("(veto_tdc[7]< 425 || veto_tdc[7] > 950)");
  TCut crv_tdc8_cut = Form("(veto_tdc[8]< 425 || veto_tdc[8] > 950)");
  TCut crv_tdc9_cut = Form("(veto_tdc[9]< 425 || veto_tdc[9] > 950)");
  TCut crv_tdc10_cut = Form("(veto_tdc[10]< 425 || veto_tdc[10] > 950)");   
  TCut crv_tdc11_cut = Form("(veto_tdc[11]< 425 || veto_tdc[11] > 950)");
  TCut crv_tdc12_cut = Form("(veto_tdc[12]< 425 || veto_tdc[12] > 950)");
  TCut crv_tdc13_cut = Form("(veto_tdc[13]< 425 || veto_tdc[13] > 950)");
  TCut crv_tdc14_cut = Form("(veto_tdc[14]< 425 || veto_tdc[14] > 950)");
  TCut crv_tdc15_cut = Form("(veto_tdc[15]< 425 || veto_tdc[15] > 950)");
  TCut crv_tdc16_cut = Form("(veto_tdc[16]< 425 || veto_tdc[16] > 950)");
  TCut crv_tdc17_cut = Form("(veto_tdc[17]< 425 || veto_tdc[17] > 950)");
  TCut crv_tdc18_cut = Form("(veto_tdc[18]< 425 || veto_tdc[18] > 950)");
  TCut crv_tdc19_cut = Form("(veto_tdc[19]< 425 || veto_tdc[19] > 950)");
  TCut crv_tdc21_cut = Form("(veto_tdc[21]< 425 || veto_tdc[21] > 950)");
  TCut crv_tdc22_cut = Form("(veto_tdc[22]< 425 || veto_tdc[22] > 950)");
  TCut crv_tdc23_cut = Form("(veto_tdc[23]< 425 || veto_tdc[23] > 950)");
  TCut crv_tdc24_cut = Form("(veto_tdc[24]< 425 || veto_tdc[24] > 950)");
  TCut crv_tdc25_cut = Form("(veto_tdc[25]< 425 || veto_tdc[25] > 950)");
  TCut crv_tdc26_cut = Form("(veto_tdc[26]< 425 || veto_tdc[26] > 950)");
  TCut crv_tdc27_cut = Form("(veto_tdc[27]< 425 || veto_tdc[27] > 950)");
 
  ////////////////////
  TCut crv_cut = crv_tdc0_cut && crv_tdc1_cut && crv_tdc2_cut && crv_tdc3_cut && crv_tdc4_cut && crv_tdc5_cut && crv_tdc6_cut && crv_tdc7_cut && crv_tdc8_cut && crv_tdc9_cut && crv_tdc10_cut && crv_tdc11_cut && crv_tdc13_cut && crv_tdc14_cut && crv_tdc15_cut && crv_tdc16_cut && crv_tdc17_cut && crv_tdc18_cut && crv_tdc19_cut && crv_tdc21_cut&& crv_tdc22_cut && crv_tdc23_cut && crv_tdc24_cut && crv_tdc25_cut && crv_tdc26_cut && crv_tdc27_cut;
  
  TCut crt_cut = !crv_cut;
  TCut vv_cut = crv_cut&&acv_cut;
  TCut tt_cut=crt_cut&&act_cut;
   

  /////////////////////
  TCut NAI_RAW = Form("max_20m1[0]>1500");
  TCut NAI_Time = Form("max_time_bin_20m1[0]<430&&max_time_bin_20m1[0]>100");
  TCut HE_muon = NAI_RAW&&NAI_Time;

  //trigger time bin
  Double_t trg_bin_up_x1 = 2025.74;
  Double_t trg_bin_up_y1 = 352.19;
  Double_t trg_bin_up_x2 = 1195.03;
  Double_t trg_bin_up_y2 = 111.602;
  Double_t trg_bin_up_slope;
  Double_t trg_bin_up_intercept;
  trg_bin_up_slope = (trg_bin_up_y2-trg_bin_up_y1)/(trg_bin_up_x2-trg_bin_up_x1);
  trg_bin_up_intercept = ((trg_bin_up_y2-trg_bin_up_y1)/(trg_bin_up_x2-trg_bin_up_x1))*(-1*trg_bin_up_x1)+trg_bin_up_y1;
  TCut trg_bin_up = Form("max_time_bin_20m1[0]<(%f*(max_time_bin_60m[0])+%f)",trg_bin_up_slope,trg_bin_up_intercept);
	  
  Double_t trg_bin_low_x1 =2123.72;
  Double_t trg_bin_low_y1 =227.804;
  Double_t trg_bin_low_x2 =1156.69;
  Double_t trg_bin_low_y2 =6.85593;
  Double_t trg_bin_low_slope;
  Double_t trg_bin_low_intercept;
  trg_bin_low_slope = (trg_bin_low_y2-trg_bin_low_y1)/(trg_bin_low_x2-trg_bin_low_x1);
  trg_bin_low_intercept = ((trg_bin_low_y2-trg_bin_low_y1)/(trg_bin_low_x2-trg_bin_low_x1))*(-1*trg_bin_low_x1)+trg_bin_low_y1;
  TCut trg_bin_low = Form("max_time_bin_20m1[0]>(%f*(max_time_bin_60m[0])+%f)",trg_bin_low_slope,trg_bin_low_intercept);
	  
  TCut trg_bin = trg_bin_low && trg_bin_up;	      

  TCut act_sel = Form("max_20m1[0]>1500");

  TCut ac_tag = crt_cut && act_sel && trg_bin;
      
  ///////
  /////////////////Basic cut start/////////////////
  TH2F *min_hist = new TH2F("min_hist","",3400,-1000,33000,100,-2000,500);
  tr->Project("min_hist","min_60m[0]:max_60m[0]",aran&&os_cut);
  min_hist->SetMarkerStyle(21);
  min_hist->SetMarkerSize(0.6);

  TH2F *min_hist_ran = new TH2F("min_hist_ran","",3400,-1000,33000,100,-2000,500);
  tr->Project("min_hist_ran","min_60m[0]:max_60m[0]",ran&&os_cut);
  min_hist_ran->SetMarkerStyle(21);
  min_hist_ran->SetMarkerSize(0.6);
  min_hist_ran->SetMarkerColor(2);
	 
  TH2F *min_hist_cut = new TH2F("min_hist_cut","",3400,-1000,33000,100,-2000,500);
  tr->Project("min_hist_cut","min_60m[0]:max_60m[0]",aran&&os_cut&&min_cut);
  min_hist_cut->SetMarkerStyle(21);
  min_hist_cut->SetMarkerSize(0.6);
  min_hist_cut->SetMarkerColor(3);
  ////////////---------------------------------------
  TH2F *ped_hist = new TH2F("ped_hist","",3400,-1000,33000,300,-1000,2000);
  tr->Project("ped_hist","ped_60m[0]:max_60m[0]",aran&&os_cut&&min_cut);
  ped_hist->SetMarkerStyle(21);
  ped_hist->SetMarkerSize(0.6);
  ped_hist->SetMarkerColor(1);

  TH2F *ped_hist_ran = new TH2F("ped_hist_ran","",3400,-1000,33000,300,-1000,2000);
  tr->Project("ped_hist_ran","ped_60m[0]:max_60m[0]",ran&&os_cut&&min_cut);
  ped_hist_ran->SetMarkerStyle(21);
  ped_hist_ran->SetMarkerSize(0.6);
  ped_hist_ran->SetMarkerColor(2);

  TH2F *ped_hist_cut = new TH2F("ped_hist_cut","",3400,-1000,33000,300,-1000,2000);
  tr->Project("ped_hist_cut","ped_60m[0]:max_60m[0]",aran&&os_cut&&min_cut&&ped_cut);
  ped_hist_cut->SetMarkerStyle(21);
  ped_hist_cut->SetMarkerSize(0.6);
  ped_hist_cut->SetMarkerColor(3);

  /////////////////////////
  TH2F *pedt_hist = new TH2F("pedt_hist","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt_hist","pedt_60m[0]:max_60m[0]",aran&&os_cut&&min_cut&&ped_cut);
  pedt_hist->SetMarkerStyle(21);
  pedt_hist->SetMarkerSize(0.6);
  pedt_hist->SetMarkerColor(1);

  TH2F *pedt_hist_ran = new TH2F("pedt_hist_ran","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt_hist_ran","pedt_60m[0]:max_60m[0]",ran&&os_cut&&min_cut&&ped_cut);
  pedt_hist_ran->SetMarkerStyle(21);
  pedt_hist_ran->SetMarkerSize(0.6);
  pedt_hist_ran->SetMarkerColor(2);

  TH2F *pedt_hist_cut = new TH2F("pedt_hist_cut","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt_hist_cut","pedt_60m[0]:max_60m[0]",aran&&os_cut&&min_cut&&ped_cut&&pedt_cut);
  pedt_hist_cut->SetMarkerStyle(21);
  pedt_hist_cut->SetMarkerSize(0.6);
  pedt_hist_cut->SetMarkerColor(3);

  
  ///////////-----------------
	 
  TH2F *min1_hist = new TH2F("min1_hist","",3400,-1000,33000,400,-1000,3000);
  tr->Project("min1_hist","min_60m[1]:max_60m[1]",aran&&allcut0);
  min1_hist->SetMarkerStyle(21);
  min1_hist->SetMarkerSize(0.6);
  min1_hist->SetMarkerColor(1);

  TH2F *min1_hist_ran = new TH2F("min1_hist_ran","",3400,-1000,33000,400,-1000,3000);
  tr->Project("min1_hist_ran","min_60m[1]:max_60m[1]",ran&&allcut0);
  min1_hist_ran->SetMarkerStyle(21);
  min1_hist_ran->SetMarkerSize(0.6);
  min1_hist_ran->SetMarkerColor(2);

  TH2F *min1_hist_cut = new TH2F("min1_hist_cut","",3400,-1000,33000,400,-1000,3000);
  tr->Project("min1_hist_cut","min_60m[1]:max_60m[1]",aran&&allcut0&&min1_cut);
  min1_hist_cut->SetMarkerStyle(21);
  min1_hist_cut->SetMarkerSize(0.6);
  min1_hist_cut->SetMarkerColor(3);
 
    
    
  TH2F *ped1_hist = new TH2F("ped1_hist","",3400,-1000,33000,460,-1000,3600);
  tr->Project("ped1_hist","ped_60m[1]:max_60m[1]",aran&&allcut0&&min1_cut);
  ped1_hist->SetMarkerStyle(21);
  ped1_hist->SetMarkerSize(0.6);
  ped1_hist->SetMarkerColor(1);
     
  TH2F *ped1_hist_ran = new TH2F("ped1_hist_ran","",3400,-1000,33000,460,-1000,3600);
  tr->Project("ped1_hist_ran","ped_60m[1]:max_60m[1]",ran&&allcut0&&min1_cut);
  ped1_hist_ran->SetMarkerStyle(21);
  ped1_hist_ran->SetMarkerSize(0.6);
  ped1_hist_ran->SetMarkerColor(2);
      	 
  TH2F *ped1_hist_cut = new TH2F("ped1_hist_cut","",3400,-1000,33000,460,-1000,3600);
  tr->Project("ped1_hist_cut","ped_60m[1]:max_60m[1]",aran&&allcut0&&min1_cut&&ped1_cut);
  ped1_hist_cut->SetMarkerStyle(21);
  ped1_hist_cut->SetMarkerSize(0.6);
  ped1_hist_cut->SetMarkerColor(3);
   
  /////////////-
  TH2F *pedt1_hist = new TH2F("pedt1_hist","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt1_hist","pedt_60m[1]:max_60m[1]",aran&&allcut0&&min1_cut&&ped1_cut);
  pedt1_hist->SetMarkerStyle(21);
  pedt1_hist->SetMarkerSize(0.6);
  pedt1_hist->SetMarkerColor(1);
 
  TH2F *pedt1_hist_ran = new TH2F("pedt1_hist_ran","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt1_hist_ran","pedt_60m[1]:max_60m[1]",ran&&allcut0&&min1_cut&&ped1_cut);
  pedt1_hist_ran->SetMarkerStyle(21);
  pedt1_hist_ran->SetMarkerSize(0.6);
  pedt1_hist_ran->SetMarkerColor(2);

  TH2F *pedt1_hist_cut = new TH2F("pedt1_hist_cut","",3400,-1000,33000,300,-1000,2000);
  tr->Project("pedt1_hist_cut","pedt_60m[1]:max_60m[1]",aran&&allcut0&&allcut1);
  pedt1_hist_cut->SetMarkerStyle(21);
  pedt1_hist_cut->SetMarkerSize(0.6);
  pedt1_hist_cut->SetMarkerColor(3);
  ///////----------------------

  TH2F *diff_pedt0_hist = new TH2F("diff_pedt0_hist","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt0_hist"," ped_60m[0]-pedt_60m[0]:max_60m[0]",aran&&allcut0&&allcut1);
  diff_pedt0_hist->SetMarkerStyle(21);
  diff_pedt0_hist->SetMarkerSize(0.6);
  diff_pedt0_hist->SetMarkerColor(1);

  TH2F *diff_pedt0_hist_ran = new TH2F("diff_pedt0_hist_ran","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt0_hist_ran","ped_60m[0]-pedt_60m[0]:max_60m[0]",ran&&allcut0&&allcut1);
  diff_pedt0_hist_ran->SetMarkerStyle(21);
  diff_pedt0_hist_ran->SetMarkerSize(0.6);
  diff_pedt0_hist_ran->SetMarkerColor(2);

  TH2F *diff_pedt0_hist_cut = new TH2F("diff_pedt0_hist_cut","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt0_hist_cut","ped_60m[0]-pedt_60m[0]:max_60m[0]",aran&&allcut0&&allcut1&&diff_pedt0_cut);
  diff_pedt0_hist_cut->SetMarkerStyle(21);
  diff_pedt0_hist_cut->SetMarkerSize(0.6);
  diff_pedt0_hist_cut->SetMarkerColor(3);

  //////////--------------------

  TH2F *diff_pedt1_hist = new TH2F("diff_pedt1_hist","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt1_hist"," ped_60m[1]-pedt_60m[1]:max_60m[1]",aran&&allcut0&&allcut1&&diff_pedt0_cut);
  diff_pedt1_hist->SetMarkerStyle(21);
  diff_pedt1_hist->SetMarkerSize(0.6);
  diff_pedt1_hist->SetMarkerColor(1);

  TH2F *diff_pedt1_hist_ran = new TH2F("diff_pedt1_hist_ran","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt1_hist_ran","ped_60m[1]-pedt_60m[1]:max_60m[1]",ran&&allcut0&&allcut1&&diff_pedt0_cut);
  diff_pedt1_hist_ran->SetMarkerStyle(21);
  diff_pedt1_hist_ran->SetMarkerSize(0.6);
  diff_pedt1_hist_ran->SetMarkerColor(2);

  TH2F *diff_pedt1_hist_cut = new TH2F("diff_pedt1_hist_cut","",3400,-1000,33000,300,-1000,2000);
  tr->Project("diff_pedt1_hist_cut","ped_60m[1]-pedt_60m[1]:max_60m[1]",aran&&Basic_cut);
  diff_pedt1_hist_cut->SetMarkerStyle(21);
  diff_pedt1_hist_cut->SetMarkerSize(0.6);
  diff_pedt1_hist_cut->SetMarkerColor(3);

  /////////////////Basic cut end/////////////////

  /////////////////AC Tag cut start/////////////////
  TH2F *NAI_hist = new TH2F("NAI_hist","",6000,0,6000,1000,0,1000);
  tr->Project("NAI_hist","max_time_bin_20m1[0]:max_time_bin_60m[0]",aran&&Basic_cut&&crt_cut&&act_sel);
  NAI_hist->SetMarkerStyle(21);
  NAI_hist->SetMarkerSize(0.6);
  NAI_hist->SetMarkerColor(1);

  TH2F *NAI_tag = new TH2F("NAI_tag","",6000,0,6000,1000,0,1000);
  tr->Project("NAI_tag","max_time_bin_20m1[0]:max_time_bin_60m[0]",aran&&Basic_cut&&crt_cut&&act_sel&&trg_bin);
  NAI_tag->SetMarkerStyle(21);
  NAI_tag->SetMarkerSize(0.6);
  NAI_tag->SetMarkerColor(2);
	 
  TH1F *NAI_1D= new TH1F("NAI_1D","",12000,-1000,5000);
  tr->Project("NAI_1D","max_20m1[0]",aran&&Basic_cut);
  NAI_1D->SetMarkerStyle(21);
  NAI_1D->SetMarkerSize(0.6);
  NAI_1D->SetMarkerColor(3); 

  /////////////////AC Tag cut end/////////////////

  /////////////////Linear cut start/////////////////
  /////////////////////--------------------  
  TH2F *aq_plot = new TH2F("aq_plot","",1950,-4e+06,40e+06,1950,-4e+03,35e+03);
  tr->Project("aq_plot","max_60m[0]:q_60m[0]",aran&&Basic_cut);
  aq_plot->SetMarkerStyle(21);
  aq_plot->SetMarkerSize(0.3);
  aq_plot->SetMarkerColor(1);

  TH2F *aq_plot_tag = new TH2F("aq_plot_tag","",1950,-4e+06,40e+06,1950,-4e+03,35e+03);
  tr->Project("aq_plot_tag","max_60m[0]:q_60m[0]",aran&&Basic_cut&&ac_tag);
  aq_plot_tag->SetMarkerStyle(21);
  aq_plot_tag->SetMarkerSize(0.3);
  aq_plot_tag->SetMarkerColor(2);
	 

  TH2F *aq_plot_cut= new TH2F("aq_plot_cut","",1950,-4e+06,40e+06,1950,-4e+03,35e+03);
  tr->Project("aq_plot_cut","max_60m[0]:q_60m[0]",aran&&Basic_cut&&qa_cut);
  aq_plot_cut->SetMarkerStyle(21);
  aq_plot_cut->SetMarkerSize(0.3);
  aq_plot_cut->SetMarkerColor(3);

  /////////------------------------
  TH2F *maxtime_bin_hist = new TH2F("maxtime_bin_hist","",3400,-1000,33000,300,-1000,2500);
  tr->Project("maxtime_bin_hist","max_time_bin_60m[0]:max_60m[0]",aran&&Basic_cut&&qa_cut);
  maxtime_bin_hist->SetMarkerStyle(21);
  maxtime_bin_hist->SetMarkerSize(0.6);
  maxtime_bin_hist->SetMarkerColor(1);

  TH2F *maxtime_bin_hist_tag = new TH2F("maxtime_bin_hist_tag","",3400,-1000,33000,300,-1000,2500);
  tr->Project("maxtime_bin_hist_tag","max_time_bin_60m[0]:max_60m[0]",aran&&Basic_cut&&qa_cut&&ac_tag);
  maxtime_bin_hist_tag->SetMarkerStyle(21);
  maxtime_bin_hist_tag->SetMarkerSize(0.6);
  maxtime_bin_hist_tag->SetMarkerColor(2);

  TH2F *maxtime_bin_hist_cut = new TH2F("maxtime_bin_hist_cut","",3400,-1000,33000,300,-1000,2500);
  tr->Project("maxtime_bin_hist_cut","max_time_bin_60m[0]:max_60m[0]",aran&&Basic_cut&&qa_cut&&time_binch0_cut);
  maxtime_bin_hist_cut->SetMarkerStyle(21);
  maxtime_bin_hist_cut->SetMarkerSize(0.6);
  maxtime_bin_hist_cut->SetMarkerColor(3);
  ////////////------------
  TH2F *maxtime_bin_hist1 = new TH2F("maxtime_bin_hist1","",3400,-1000,33000,700,-1000,6000);
  tr->Project("maxtime_bin_hist1","max_time_bin_60m[1]:max_60m[1]",aran&&Basic_cut&&qa_cut&&time_binch0_cut);
  maxtime_bin_hist1->SetMarkerStyle(21);
  maxtime_bin_hist1->SetMarkerSize(0.6);
  maxtime_bin_hist1->SetMarkerColor(1);

  TH2F *maxtime_bin_hist_tag1= new TH2F("maxtime_bin_hist_tag1","",3400,-1000,33000,700,-1000,6000);
  tr->Project("maxtime_bin_hist_tag1","max_time_bin_60m[1]:max_60m[1]",aran&&Basic_cut&&qa_cut&&time_binch0_cut&&ac_tag);
  maxtime_bin_hist_tag1->SetMarkerStyle(21);
  maxtime_bin_hist_tag1->SetMarkerSize(0.6);
  maxtime_bin_hist_tag1->SetMarkerColor(2);

  TH2F *maxtime_bin_hist_cut1 = new TH2F("maxtime_bin_hist_cut1","",3400,-1000,33000,700,-1000,6000);
  tr->Project("maxtime_bin_hist_cut1","max_time_bin_60m[1]:max_60m[1]",aran&&Basic_cut&&qa_cut&&timebin_cut);
  maxtime_bin_hist_cut1->SetMarkerStyle(21);
  maxtime_bin_hist_cut1->SetMarkerSize(0.6);
  maxtime_bin_hist_cut1->SetMarkerColor(3);

  ///////////////--------------
  TH2F *s12_hist = new TH2F("s12_hist","",1950,-4e+03,40e+03,1950,-4e+03,35e+03);
  tr->Project("s12_hist","max_60m[1]:max_60m[0]",aran&&Basic_cut&&qa_cut&&timebin_cut);
  s12_hist->SetMarkerStyle(21);
  s12_hist->SetMarkerSize(0.6);
  s12_hist->SetMarkerColor(1);
 
  TH2F *s12_ac_tag = new TH2F("s12_ac_tag","",1950,-4e+03,40e+03,1950,-4e+03,35e+03);
  tr->Project("s12_ac_tag","max_60m[1]:max_60m[0]",aran&&Basic_cut&&qa_cut&&timebin_cut&&ac_tag);
  s12_ac_tag->SetMarkerStyle(21);
  s12_ac_tag->SetMarkerSize(0.6);
  s12_ac_tag->SetMarkerColor(2);
	 
  TH2F *s12_hist_cut= new TH2F("s12_hist_cut","",1950,-4e+03,40e+03,1950,-4e+03,35e+03);
  tr->Project("s12_hist_cut","max_60m[1]:max_60m[0]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut);
  s12_hist_cut->SetMarkerStyle(21);
  s12_hist_cut->SetMarkerSize(0.6);
  s12_hist_cut->SetMarkerColor(3); 
  
 
  /////**********
  TH2F *M_hist = new TH2F("M_hist","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("M_hist","max_200m1[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut);
  M_hist->SetMarkerStyle(21);
  M_hist->SetMarkerSize(0.6);
  M_hist->SetMarkerColor(1);

  TH2F *M_ac_tag = new TH2F("M_ac_tag","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("M_ac_tag","max_200m1[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut&&ac_tag);
  M_ac_tag->SetMarkerStyle(21);
  M_ac_tag->SetMarkerSize(0.6);
  M_ac_tag->SetMarkerColor(2);
	 

  TH2F *M_hist_cut= new TH2F("M_hist_cut","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("M_hist_cut","max_200m1[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut&&thd_cut);
  M_hist_cut->SetMarkerStyle(21);
  M_hist_cut->SetMarkerSize(0.6);
  M_hist_cut->SetMarkerColor(3);
  ///////-----------------------------------

  TH2F *pq_hist = new TH2F("pq_hist","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("pq_hist","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut&&thd_cut);
  pq_hist->SetMarkerStyle(21);
  pq_hist->SetMarkerSize(0.6);
  pq_hist->SetMarkerColor(1);

  TH2F *pq_ac_tag = new TH2F("pq_ac_tag","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("pq_ac_tag","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut&&thd_cut&&ac_tag);
  pq_ac_tag->SetMarkerStyle(21);
  pq_ac_tag->SetMarkerSize(0.6);
  pq_ac_tag->SetMarkerColor(2);
	 

  TH2F *pq_hist_cut= new TH2F("pq_hist_cut","",1950,-8e+06,40e+06,1950,-10e+03,35e+03);
  tr->Project("pq_hist_cut","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&qa_cut&&timebin_cut&&s12_cut&&thd_cut&&pq_cut);
  pq_hist_cut->SetMarkerStyle(21);
  pq_hist_cut->SetMarkerSize(0.6);
  pq_hist_cut->SetMarkerColor(3);
  /////////////////Linear cut end/////////////////
  /*
/////////////////global edge cut start/////////////////

TH2F *Final_2d_Plot = new TH2F("Final_2d_Plot","",20000,-5e+06,35e+06,3400,-10e+03,33e+03);
tr->Project("Final_2d_Plot","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&Linear_cut);
Final_2d_Plot->SetMarkerStyle(21);
Final_2d_Plot->SetMarkerSize(0.6);
Final_2d_Plot->SetMarkerColor(1);

TH2F *Final_2d_Plot_vv = new TH2F("Final_2d_Plot_vv","",20000,-5e+06,35e+06,3400,-10e+03,33e+03);
tr->Project("Final_2d_Plot_vv","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&Linear_cut&&vv_cut);
Final_2d_Plot_vv->SetMarkerStyle(21);
Final_2d_Plot_vv->SetMarkerSize(0.6);
Final_2d_Plot_vv->SetMarkerColor(2);

TH2F *Final_2d_Plot_tt = new TH2F("Final_2d_Plot_tt","",20000,-5e+06,35e+06,3400,-10e+03,33e+03);
tr->Project("Final_2d_Plot_tt","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&Linear_cut&&ac_tag);
Final_2d_Plot_tt->SetMarkerStyle(21);
Final_2d_Plot_tt->SetMarkerSize(0.6);
Final_2d_Plot_tt->SetMarkerColor(4);


TH2F *Final_2d_Plot_psd_cut = new TH2F("Final_2d_Plot_psd_cut","",20000,-5e+06,35e+06,3400,-10e+03,33e+03);
tr->Project("Final_2d_Plot_psd_cut","max_60m[0]:optq_60m[2]",aran&&Basic_cut&&Linear_cut&&vv_cut&&globle_edge_cut);
Final_2d_Plot_psd_cut->SetMarkerStyle(21);
Final_2d_Plot_psd_cut->SetMarkerSize(0.6);
Final_2d_Plot_psd_cut->SetMarkerColor(3);



TH2F *VV_log_base_10 = new TH2F("VV_log_base_10","",1800,0,18,2000,-10,10);
tr->Project("VV_log_base_10","log10(log(19)/(200.0*fslope_200m)):((2.59651e+03+(5.71669e-03*optq_60m[2]))/10000)",aran&&Basic_cut&&Linear_cut&&vv_cut&&globle_edge_cut);
VV_log_base_10->SetMarkerStyle(21);
VV_log_base_10->SetMarkerSize(0.6);
VV_log_base_10->SetMarkerColor(1);

TH2F *VV_log_baseE = new TH2F("VV_log_baseE","",1800,0,18,2000,-10,10);
tr->Project("VV_log_baseE","(log(19)/(200.0*fslope_200m)):((2.59651e+03+(5.71669e-03*optq_60m[2]))/10000)",aran&&Basic_cut&&Linear_cut&&vv_cut&&globle_edge_cut);
VV_log_baseE->SetMarkerStyle(21);
VV_log_baseE->SetMarkerSize(0.6);
VV_log_baseE->SetMarkerColor(1);

TH2F *VV_linear = new TH2F("VV_linear","",1800,0,18,2000,-1000,1000);
tr->Project("VV_linear","(log(19)/(fslope_200m)):((2.59651e+03+(5.71669e-03*optq_60m[2]))/10000)",aran&&Basic_cut&&Linear_cut&&vv_cut&&globle_edge_cut);
VV_linear->SetMarkerStyle(21);
VV_linear->SetMarkerSize(0.6);
VV_linear->SetMarkerColor(1);

 

/////////////////global edge cut end/////////////////
*/
  TFile f1("Basic_cuts_160415.root","recreate");
  min_hist->Write();
  min_hist_ran->Write();
  min_hist_cut->Write();
  ped_hist->Write();
  ped_hist_ran->Write();
  ped_hist_cut->Write();
  pedt_hist->Write();
  pedt_hist_ran->Write();
  pedt_hist_cut->Write();
  min1_hist->Write();
  min1_hist_ran->Write();
  min1_hist_cut->Write();
  ped1_hist->Write();
  ped1_hist_ran->Write();
  ped1_hist_cut->Write();
  pedt1_hist->Write();
  pedt1_hist_ran->Write();
  pedt1_hist_cut->Write();
  diff_pedt0_hist->Write();
  diff_pedt0_hist_ran->Write();
  diff_pedt0_hist_cut->Write();
  diff_pedt1_hist->Write();
  diff_pedt1_hist_ran->Write();
  diff_pedt1_hist_cut->Write();
      
  min_hist->Delete();
  min_hist_ran->Delete();
  min_hist_cut->Delete();
  ped_hist->Delete();
  ped_hist_ran->Delete();
  ped_hist_cut->Delete();
  pedt_hist->Delete();
  pedt_hist_ran->Delete();
  pedt_hist_cut->Delete();
  min1_hist->Delete();
  min1_hist_ran->Delete();
  min1_hist_cut->Delete();
  ped1_hist->Delete();
  ped1_hist_ran->Delete();
  ped1_hist_cut->Delete();
  pedt1_hist->Delete();
  pedt1_hist_ran->Delete();
  pedt1_hist_cut->Delete();
  diff_pedt0_hist->Delete();
  diff_pedt0_hist_ran->Delete();
  diff_pedt0_hist_cut->Delete();
  diff_pedt1_hist->Delete();
  diff_pedt1_hist_ran->Delete();
  diff_pedt1_hist_cut->Delete();


  TFile f2("AC_tag_selection_160422.root","recreate");
  NAI_hist->Write(); NAI_hist->Delete();
  NAI_tag->Write();NAI_tag->Delete();
  NAI_1D->Write(); NAI_1D->Delete();

      
  TFile f3("Linear_cuts_160422.root","recreate");
  aq_plot->Write(); 
  aq_plot_tag->Write();
  aq_plot_cut->Write() ; 
  maxtime_bin_hist->Write();
  maxtime_bin_hist_tag->Write();
  maxtime_bin_hist_cut->Write();
  maxtime_bin_hist1->Write();
  maxtime_bin_hist_tag1->Write();
  maxtime_bin_hist_cut1->Write();
  s12_hist->Write(); 
  s12_ac_tag->Write(); 
  s12_hist_cut->Write(); 
  M_hist->Write(); 
  M_ac_tag->Write(); 
  M_hist_cut->Write();
  pq_hist->Write(); 
  pq_ac_tag->Write(); 
  pq_hist_cut->Write();
 
  aq_plot->Delete(); 
  aq_plot_tag->Delete();
  aq_plot_cut->Delete() ; 
  maxtime_bin_hist->Delete();
  maxtime_bin_hist_tag->Delete();
  maxtime_bin_hist_cut->Delete();
  maxtime_bin_hist1->Delete();
  maxtime_bin_hist_tag1->Delete();
  maxtime_bin_hist_cut1->Delete();
  s12_hist->Delete(); 
  s12_ac_tag->Delete(); 
  s12_hist_cut->Delete(); 
  M_hist->Delete(); 
  M_ac_tag->Delete(); 
  M_hist_cut->Delete();
  pq_hist->Delete(); 
  pq_ac_tag->Delete(); 
  pq_hist_cut->Delete();

  /* Final_2d_Plot->Write();
     Final_2d_Plot_vv->Write();
     Final_2d_Plot_tt->Write();
     Final_2d_Plot_psd_cut->Write();
     
     TFile f4("BS_cuts_160422.root","recreate");
     VV_log_base_10->Write();
     VV_log_baseE->Write();
     VV_linear->Write();
  */
  
          
  // Efficiency of the Cut
  //[1]Number of RAN
  int n_ran;
  TH1F *ran_hist = new TH1F("ran_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("ran_hist","max_60m[0]",ran);
  n_ran = ran_hist->GetEntries();
  ran_hist->Delete(); 
  printf("Ran event number: %d\n", n_ran);
  
  //[2]Basic_cut Eff
  UInt_t n_Basic_ran;
  TH1F *Basic_hist = new TH1F("Basic_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("Basic_hist","max_60m[0]",ran&&Basic_cut);
  n_Basic_ran = Basic_hist->GetEntries();
  Basic_hist->Delete();
  float eff_Basic;
  eff_Basic = (float)n_Basic_ran/(float)n_ran;
  printf("Basic eff:%f\n",eff_Basic);
     
    
  //[12]NaI-HE Muon
  UInt_t n_nai_he;
  TH1F *nai_he = new TH1F("nai_he", "nai_he", 6, -1000, 35000);
  tr->Project("nai_he","max_60m[0]",aran&&Basic_cut&&HE_muon);
  n_nai_he=nai_he->GetEntries();
  printf("total H.E. mu: %d\n",n_nai_he);
  
  
  //[14](NaI-HE) Muon + CRT Cut
  UInt_t n_nai_he_mu_tag;
  TH1F *nai_he_mu_tag = new TH1F("nai_he_mu_tag", "nai_he_mu_tag", 6, -1000, 35000);
  tr->Project("nai_he_mu_tag","max_60m[0]",aran&&Basic_cut&&HE_muon&&crt_cut);
  n_nai_he_mu_tag=nai_he_mu_tag->GetEntries();
  
  float lambda_value = (float)n_nai_he_mu_tag/(float)n_nai_he;
  printf("total H.E. mu tag: %d\n",n_nai_he_mu_tag);
  printf("\nLambda  %f \n",(float)n_nai_he_mu_tag/(float)n_nai_he);

  //////////////////////////////////////////////////
  char EA[1000];
  sprintf(EA,"%.8e+(%.8e*max_60m[0])",9.92217014e-02, 4.01297019e-04);

  char outfilename[200];
  sprintf(outfilename,"spec_rootfile/E_spec_%d_50eV.root",timeinfo[jj]);  
  
  int min=0;
  int max=14;
  int bin=280;
 
  // Efficiency of the Cut
  //[1]Number of RAN
  int n_ran;
  TH1F *ran_hist = new TH1F("ran_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("ran_hist","max_60m[0]",gran);
  n_ran = ran_hist->GetEntries();
  ran_hist->Delete(); 
  printf("Ran event number: %d\n", n_ran);
      
  //[2]Offset_cut Eff
  UInt_t n_offset_ran;
  TH1F *offset_hist = new TH1F("offset_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("offset_hist","max_60m[0]",gran&&goffset_cut);
  n_offset_ran = offset_hist->GetEntries();
  offset_hist->Delete();
  float eff_offset;
  eff_offset = (float)n_offset_ran/(float)n_ran;
  printf("Offset eff:%f\n",eff_offset);
      
  //[3]ACV Cut Eff
  UInt_t n_acv_ran;
  TH1F *ran_acv_hist = new TH1F("ran_acv_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("ran_acv_hist","max_60m[0]",gran&&acv_cut);
  n_acv_ran = ran_acv_hist->GetEntries();
  ran_acv_hist->Delete();
  float eff_acv;
  eff_acv = (float)n_acv_ran/(float)n_ran;
  printf("ACV eff:%f\n",eff_acv);
      
  //[4]CRV Cut Eff
  UInt_t n_crv_ran;
  TH1F *ran_crv_hist = new TH1F("ran_crv_hist", "Random Trigger", 340, -1000, 33000);
  tr->Project("ran_crv_hist","max_60m[0]",gran&&crv_cut);
  n_crv_ran = ran_crv_hist->GetEntries();
  ran_crv_hist->Delete();
  float eff_crv;
  eff_crv = (float)n_crv_ran/(float)n_ran;
  printf("CRV eff:%f\n",eff_crv);
      
  printf("Live time: %f\n",livetime[jj]);

        	
  //[12]NaI-HE Muon
  UInt_t n_nai_he;
  TH1F *nai_he = new TH1F("nai_he", "nai_he", 6, -1000, 35000);
  tr->Project("nai_he","max_60m[0]",garan&&goffset_cut&&glinear_cut&&HE_muon);
  n_nai_he=nai_he->GetEntries();
  printf("total H.E. mu: %d\n",n_nai_he);
     
      
  //[14](NaI-HE) Muon + CRT Cut
  UInt_t n_nai_he_mu_tag;
  TH1F *nai_he_mu_tag = new TH1F("nai_he_mu_tag", "nai_he_mu_tag", 6, -1000, 35000);
  tr->Project("nai_he_mu_tag","max_60m[0]",garan&&goffset_cut&&glinear_cut&&HE_muon&&crt_cut);
  n_nai_he_mu_tag=nai_he_mu_tag->GetEntries();
  printf("total H.E. mu tag: %d\n",n_nai_he_mu_tag);
  printf("\nLambda  %f \n",(float)n_nai_he_mu_tag/(float)n_nai_he);


  TFile *outfile = new TFile(outfilename,"recreate");
  
  TH1F *N_RAN = new TH1F("N_RAN", "N_RAN", 10, 0,0);
  N_RAN->Fill(n_ran);
  N_RAN->Write();
  N_RAN->Delete();

  TH1F *N_NaI_HE = new TH1F("N_NaI_HE", "N_NaI_HE", 10,0,0);
  N_NaI_HE->Fill(n_nai_he);
  N_NaI_HE->Write();
  N_NaI_HE->Delete();

  TH1F *N_NaI_HE_MU_TAG = new TH1F("N_NaI_HE_MU_TAG", "N_NaI_HE_MU_TAG", 10, 0,0);
  N_NaI_HE_MU_TAG->Fill(n_nai_he_mu_tag);
  N_NaI_HE_MU_TAG->Write();
  N_NaI_HE_MU_TAG->Delete();

  TH1F *livetime_hist = new TH1F("livetime_hist","",10,0,0);
  TH1F *eff_Only_OFFSET = new TH1F("eff_Only_OFFSET","",10,0,0);
  TH1F *eff_Only_ACV = new TH1F("eff_Only_ACV","",10,0,0);
  TH1F *eff_Only_CRV = new TH1F("eff_Only_CRV","",10,0,0);
      
  livetime_hist->Fill(livetime[jj]);
  eff_Only_OFFSET->Fill(eff_offset);
  eff_Only_ACV->Fill(eff_acv);
  eff_Only_CRV->Fill(eff_crv);
     
  livetime_hist->Write();
  eff_Only_OFFSET->Write();
  eff_Only_ACV->Write();
  eff_Only_CRV->Write();
      
  livetime_hist->Delete();
  eff_Only_OFFSET->Delete();
  eff_Only_ACV->Delete();
  eff_Only_CRV->Delete();
  ///////////////////////////////// 
  TH1F *E_RAN = new TH1F("E_RAN","",bin,min,max);
  tr->Project("E_RAN",EA,gran&&goffset_cut);
  E_RAN->Write();
  E_RAN->Delete();
  
  TH1F *E_RAW = new TH1F("E_RAW","",bin,min,max);
  tr->Project("E_RAW",EA,garan);
  E_RAW->Write();
  E_RAW->Delete();
  
  TH1F *E_minRAW = new TH1F("E_minRAW","",bin,min,max);
  tr->Project("E_minRAW",EA,garan&& gminch0_cut&&"max_60m[0]<32750");
  E_minRAW->Write();
  E_minRAW->Delete();
  
  TH1F *E_ORAW = new TH1F("E_ORAW","",bin,min,max);
  tr->Project("E_ORAW",EA,garan&&goffset_cut);
  E_ORAW->Write();
  E_ORAW->Delete();
  
  TH1F *E_OLRAW = new TH1F("E_OLRAW","",bin,min,max);
  tr->Project("E_OLRAW",EA,garan&&goffset_cut&&glinear_cut);
  E_OLRAW->Write();
  E_OLRAW->Delete();

  
  TH1F *E_CRV = new TH1F("E_CRV","",bin,min,max);
  tr->Project("E_CRV",EA,garan&&goffset_cut&&glinear_cut&&crv_cut);
  E_CRV->Write();
  E_CRV->Delete();
  
  TH1F *E_ACV = new TH1F("E_ACV","",bin,min,max);
  tr->Project("E_ACV",EA,garan&&goffset_cut&&glinear_cut&&acv_cut); 
  E_ACV->Write();
  E_ACV->Delete();
  
  

  TH1F *E_CRVACV = new TH1F("E_CRVACV","",bin,min,max);
  tr->Project("E_CRVACV",EA,garan&&goffset_cut&&glinear_cut&&crv_cut&&acv_cut);
  E_CRVACV->Write();
  E_CRVACV->Delete();
  
  TH1F *E_CRVACT = new TH1F("E_CRVACT","",bin,min,max);
  tr->Project("E_CRVACT",EA,garan&&goffset_cut&&glinear_cut&&crv_cut&&act_cut);
  E_CRVACT->Write();
  E_CRVACT->Delete();

  TH1F *E_CRTACV = new TH1F("E_CRTACV","",bin,min,max);
  tr->Project("E_CRTACV",EA,garan&&goffset_cut&&glinear_cut&&crt_cut&&acv_cut);
  E_CRTACV->Write();
  E_CRTACV->Delete();

  TH1F *E_CRTACT = new TH1F("E_CRTACT","",bin,min,max);
  tr->Project("E_CRTACT",EA,garan&&goffset_cut&&glinear_cut&&crt_cut&&act_cut);
  E_CRTACT->Write();
  E_CRTACT->Delete();
 

      
  //Energy Depended PSD
  TH1F *E_actag_af_offset_cut = new TH1F("E_actag_af_offset_cut","",bin,min,max);
  tr->Project("E_actag_af_offset_cut",EA,garan&&goffset_cut&&ac_tag);
  E_actag_af_offset_cut->Write();
  E_actag_af_offset_cut->Delete();

  
  TH1F *E_actag_af_linear_cut = new TH1F("E_actag_af_linear_cut","",bin,min,max);
  tr->Project("E_actag_af_linear_cut",EA,garan&&goffset_cut&&ac_tag&&glinear_cut);
  E_actag_af_linear_cut->Write();
  E_actag_af_linear_cut->Write();
  E_actag_af_linear_cut->Delete();
  outfile->Close();
 
  printf("Histogram written \n");
      
}







} 

