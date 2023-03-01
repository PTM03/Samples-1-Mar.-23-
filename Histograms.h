

#include <cmath>
#include <math.h>


void nano9Ana::BookHistograms()
{
  //The histograms are booked here.
  //Binning etc are done here.
  //These histograms are stored in the hst_<process name>.root file in the same order.
  
  //Example : new TH1F ("hst_name", "hst title", NBins, startVal, EndVal);
  

  h.nmu = new TH1F("nmuons", "Number of Muons", 10, 0, 10);
  h.mupt = new TH1F("muonpt", "Muon pT", 200, 0, 200);
  
  h.muprop[0] = new TH1F("mu0_pt","Leading muon pT",200,0,200);
  h.muprop[1] = new TH1F("mu0_eta","Leading muon Eta",120,-3.,3.);
  h.muprop[2] = new TH1F("mu0_isol","Leading muon PFRelIso",100,0,2);  

  h.nElectron =  new TH1F("nElectrons", "Number of Electrons", 10, 0, 10);


  // h.nLepton[0] =  new TH1F("n_1Lp_evt", "Number of Leptons", 10, 0, 10);
  // h.nLepton[1] =  new TH1F("n_2Lp_evt", "Number of Leptons", 10, 0, 10);

  h.nLep[2] =  new TH1F("n_3Lp_evt", "Number of Leptons", 10, 0, 10);
  h.nLepton = new TH1F("nLeptons","number of Leptons",10,0,10);

  h.LT[0] = new TH1F ("LT_1lp_evt","LT(Min 1 Lepton events)",150,0,1500);
  h.LT[1] = new TH1F ("LT_2lp_evt","LT(Min 2 Lepton events)",150,0,1500);
  h.LT[2] = new TH1F ("LT_3lp_evt","LT(Min 3 Lepton events)",150,0,1500);
  
  h.HT[0] = new TH1F ("HT_1lp_evt","HT(Min 1 Lepton events)",150,0,1500);
  h.HT[1] = new TH1F ("HT_2lp_evt","HT(Min 2 Lepton events)",150,0,1500);
  h.HT[2] = new TH1F ("HT_3lp_evt","HT(Min 3 Lepton events)",150,0,1500);

  //Jets
 

  h.Jet[0] = new TH1F ("Jet Pt_1lp_evt","JetPT_1Lp_evt",1000,0,1000);
  h.Jet[1] = new TH1F ("LeadingJet_Pt","leading_Jet_Pt",1000,0,1000);
  h.Jet[2] = new TH1F ("b_tag_score(min_3lp_evt)","btag_score_(min_3lp_evt)",200,0,20);
  h.Jet[3] = new TH1F ("Number of Jets(min_1lp_evt)","number of Jets", 10,0,15);
  h.Jet[4] = new TH1F ("Number of Jets(min_2lp_evt)","number of Jets", 10,0,15);
  h.Jet[5] = new TH1F ("Number of Jets(min_3lp_evt)","number of Jets", 10,0,15);
  h.Jet[6] = new TH1F ("leading_b_tag_score(min_3lp_evt)","leading_btag_score_(min_3lp_evt)",200,0,20);
  h.Jet[7] = new TH1F ("subleading_b_tag_score(min_3lp_evt)","s-leading_btag_score_(min_3lp_evt)",200,0,20);
  h.Jet[8] = new TH1F ("sub-subleading_b_tag_score(min_3lp_evt)","s-s-leading_btag_score_(min_3lp_evt)",200,0,20);
  

  //MET 
  h.dPhi_MET = new  TH1F ("dPHi_MET", "Transverse_Mass_3lp",100,-4,4);
  h.MET[0] = new TH1F("1_lp_evt_MET_Pt"," MET Pt for 1 Lepton events",200,0,200);
  h.MET[1] = new TH1F("2_lp_evt_MET_Pt"," MET Pt for 2 Lepton events",200,0,200);
  h.MET[2] = new TH1F("3_lp_evt_MET_Pt"," MET Pt for 3 Lepton events",200,0,200);
  h.MET[3] = new TH1F("MET_phi_3lp_evt" ,"MET_phi_3lp_evt",100,-4,4);
  h.MET[4] = new TH1F ("MET_ll_dphi_3lp_evt" ,"MET and ll_dphi_3lp_evt",100,0,4);
  h.MET[5] = new TH1F ("MET_sll_dphi_3lp_evt" ,"MET and sll_dphi_3lp_evt",100,0,4);
  h.MET[6] = new TH1F ("MET_ssll_dphi_3lp_evt" ,"MET and ssll_dphi_3lp_evt",100,0,4);


  //Transverse Mass
  h.mT[0] = new TH1F ("Transv_Mass_1Lp_0", "Transverse_Mass_1lp(leading lepton)", 100,0,1000);
  h.mT[1] = new TH1F ("Transv_Mass_2Lp_0", "Transverse_Mass_2lp(leading lepton)", 100,0,1000);
  h.mT[2] = new TH1F ("Transv_Mass_2Lp_1", "Transverse_Mass_2lp(Sub-leading lepton)", 100,0,1000);
  h.mT[3] = new TH1F ("Transv_Mass_3Lp_0", "Transverse_Mass_3lp(leading lepton)",  100,0,1000);
  h.mT[4] = new TH1F ("Transv_Mass_3Lp_1", "Transverse_Mass_3lp(subleading lepton)",100,0,1000);
  h.mT[5] = new TH1F ("Transv_Mass_3Lp_2", "Transverse_Mass_3lp(ssubleading lepton)",100,0,1000);
  h.mT[6] = new TH1F ("square_Transv_Mass_3Lp_0", "Transverse_Mass_3lp(leading lepton)",  100,0,100000);
  h.mT[7] = new TH1F ("mT(Geometric_Mean)", "Geometric mean of Transverse Mass_3lp",  100,0,1000);
  h.mT[8] = new TH1F ("Avg_Transv_Mass_3Lp", "Sum of Transverse_Mass_3lp",  100,0,1000);
  h.mT[9] = new TH1F ("cube_Transv_Mass_3Lp_0", "Transverse_Mass_3lp(leading lepton)",  100,0,1000000);
  h.mT[10] = new TH1F ("FourthPower_Transv_Mass_3Lp_0", "Transverse_Mass_3lp(leading lepton)",  100,0,1000000);
   h.mT[11] = new TH1F ("FifthPower_Transv_Mass_3Lp_0", "Transverse_Mass_3lp(leading lepton)",  100,0,1000000);
   h.mT[12] = new TH1F ("max_transv_mass(min3levt)", "Max_transverse_mass(min_3lp_evt)", 100,0,500);
  //Delta Phi MET
 
  

  //Invariant Mass
  h.Invar_Mass[0] = new TH1F ("Invar_Mass_2Lp_0_1", "Invar_Mass_2lp_0_1",100,0,1000);
  h.Invar_Mass[1] = new TH1F ("Invar_Mass_3Lp_0_1", "Invar_Mass_3lp_0_1",100,0,1000);
  h.Invar_Mass[2] = new TH1F ("Invar_Mass_3Lp_0_2", "Invar_Mass_3lp_0_2",100,0,1000);
  h.Invar_Mass[3] = new TH1F ("Invar_Mass_3Lp_1_2", "Invar_Mass_3lp_1_2",100,0,1000);
  h.Invar_Mass[4] = new TH1F ("Tri_Invar_Mass_3Lp", "Tri_Invar_Mass_3lp",100,0,2000);
  h.Invar_Mass[5] = new TH1F ("Tetra_Invar_Mass_3Lp", "Tetra_Invar_Mass_3lp",100,0,2000);
  h.Invar_Mass[6] = new TH1F ("GM_Invar_Mass_3Lp", "Tri_Invar_Mass_3lp",100,0,2000);
  h.Invar_Mass[7] = new TH1F ("AVG_Invar_Mass_3Lp", "AVG_Invar_Mass_3lp",100,0,2000);
  h.Invar_Mass[8] = new TH1F("max_invarmass_min_3lp_evt","max invar mass(min 3lp evt)",100,0,500);

  //Delta_plots
  h.Delta[0] = new TH1F ("Delta_eta(min_3lp_evt)", "delta_eta(Min_3lp_evt)",100,0,4);
  h.Delta[1] = new TH1F ("Delta_Phi_01(min_3lp_evt)", "delta_phi_01(Min_3lp_evt)",100,0,4);
  h.Delta[2] = new TH1F ("Delta_R(min_3lp_evt)", "delta_R(Min_3lp_evt)",100,0,4);  
  h.Delta[3] = new TH1F ("Delta_Phi_02(min_3lp_evt)", "delta_phi_02(Min_3lp_evt)",100,0,4);
  h.Delta[4] = new TH1F ("Delta_Phi_12(min_3lp_evt)", "delta_phi_12(Min_3lp_evt)",100,0,4);
  
  
}
