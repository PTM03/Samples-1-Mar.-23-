#define nano9Ana_cxx
// The class definition in nano9Ana.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived

// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

#include "nano9Ana.h"
#include <TH2.h>
#include <TStyle.h>
#include <cmath>
#include <math.h>
#include "Functions.C"
#include "Histograms.h"

void nano9Ana::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
}

void nano9Ana::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  //Initialization of the counters:
  nEvtRan     = 0;
  nEvtTotal   = 0;
 
  //Other custom counters can be initialized here.
  
  _HstFile = new TFile(_HstFileName,"recreate");
  BookHistograms();
  
  Events1 =-1;
  Events2 =-1;
  Events3 =-1;

  btag_score = 0;
 
  
 
}

void nano9Ana::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
  _HstFile->Write();
  _HstFile->Close();
  
  //The following lines are displayed on the root prompt.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total good events = "<<nEvtTotal<<endl;
  
  //The following lines are written on the sum_<process name>.txt file
  ofstream fout(_SumFileName);
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total good events  = "<<nEvtTotal<<endl;
  

  fout<<"Number of Events having only one Lepton" <<Events1<<endl;
  fout<<"Number of Events having Two Leptons" <<Events2<<endl;
  fout<<"Number of Events having Three Leptons" <<Events3<<endl;     


 
  
}

void nano9Ana::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file. 
}

Bool_t nano9Ana::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  





  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC.SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);
  
  //Verbosity determines the number of processed events after which the root prompt is supposed to display a status update.
  if(_verbosity==0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
  else if(_verbosity>0 && nEvtTotal%10000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
  
  //The following flags throws away some events based on unwanted properties (such as detector problems)
  GoodEvt2018 = (_year==2018 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2017 = (_year==2017 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  GoodEvt2016 = (_year==2016 ? *Flag_goodVertices && *Flag_globalSuperTightHalo2016Filter && *Flag_HBHENoiseFilter && *Flag_HBHENoiseIsoFilter && *Flag_EcalDeadCellTriggerPrimitiveFilter && *Flag_BadPFMuonFilter && (_data ? *Flag_eeBadScFilter : 1) : 1);
  
  GoodEvt = GoodEvt2018 && GoodEvt2017 && GoodEvt2016;
  
  nEvtRan++;                             //Total number of events containing everything (including the trash events).
  
  
  
  
  if(GoodEvt){
    nEvtTotal++;                         //Total number of events containing goodEvents
    
    goodLep.clear();
    //Construction of the arrays:
    
    //goodMu array :
    int nmu = 0;                         // This counts the number of muons in each event.
    goodMu.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nMuon); i++){
                                         // This loop runs over all the muon candidates. Some of them will pass our selection criteria.
                                         // These will be stored in the goodMu array.
      Lepton temp;                       // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105); //the muon mass in GeV is 0.105
      temp.id = -13*Muon_charge[i];      //pdgID for mu- = 13, pdgID for mu+ = -13  
      temp.ind = i;

      //These are the flags the 'temp' object i.e. the muon candidate has to pass.
      bool passCuts = temp.v.Pt()>15 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
      //passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
      passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;
      
      if(passCuts){
	goodMu.push_back(temp);          // If 'temp' satisfies all the conditions, it is pushed back into goodMu
	goodLep.push_back(temp);

	
      }
    }                                    // This 'for' loop has created a goodMu array.
    
    //Now we sort the goodMu in decreasing pT order
    Sort(0);                           

    //Other arrays, such as RecoEle, GenMu, GenEle can be constructed here.

         

 //goodElectron array :     
    int nEle = 0;                         // This counts the number of electrons in each event.
    goodElectron.clear();                      // Make sure to empty the array from previous event.
    for(unsigned int i=0; i<(*nElectron); i++){
      // This loop runs over all the Electron candidates. Some of them will pass our selection criteria.



      // These will be stored in the goodElectron array.
      Lepton temp;                                                              // 'temp' is the i-th candidate.
      temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.005);//Electron mass in GeV is 0.000000511
      temp.id = -11*Electron_charge[i];                                         //pdgID for mu- = 11, pdgID for mu+ = -11  
      temp.ind = i;
      
      
      //  h.Electron_cutBased->Fill(Electron_cutBased[i]);
      

      //These are the flags the 'temp' object i.e. the Electron candidate has to pass.
      bool passCuts = temp.v.Pt() > 15 && fabs(temp.v.Eta()) < 2.4  && Electron_cutBased[i]>2 ;
      bool Well_Separated = true;
      //Separations 
      for(int j=0; j<(int)goodMu.size(); j++){
	float Electron_Muon_Separation = temp.v.DeltaR(goodMu.at(j).v);
	if(Electron_Muon_Separation < 0.4) Well_Separated = false;
      }
      
      bool isprompt = false;
      if(fabs(temp.v.Eta())<=1.479){
	if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	  isprompt = true;
      }
      
      
      //   if(fabs(temp.v.Eta())>1.479){
      //if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
      //isprompt = true;
      //  }
      

      passCuts = passCuts && fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1;
      passCuts = passCuts && Well_Separated && isprompt;
      

      if(passCuts ){
	goodElectron.push_back(temp);
	goodLep.push_back(temp);
      }  
    } // This 'for' loop has created a goodElectron array.
    

    //Now we sort the goodElectron in decreasing pT order

    Sort(1);                           


    Sort(6);


 ////////////////////////////////////////////goodJet Array/////////////////////////////////////////


    
    int njet = 0;                   //counts no. of Photons in each events
    goodJet.clear();                //emptying the array
    for(unsigned int i=0; i<(*nJet);i++){
      
      
      Lepton temp,temp2;
      temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],0);
      // temp.id= ;                //pdg ID for Jet is .....since Photon has no Charge keeping the id constant over the entire array

      temp.ind =i;                                                          

      bool Separated_Jets = true;
      /*
      for(unsigned int j=0;j!=i && j<(*nJet); j++){
	temp2.v.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],0);
	temp2.ind =j;
	if(temp.v.DeltaR(temp2.v)<2) Separated_Jets=false;
	}*/
      
      bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<5.0;
      bool Well_Separated_Mu = true;
      bool Well_Separated_Ele = true;

      // Separating fake Jets      
      for(int j=0;j<(int)goodMu.size();j++){
	float Jet_Muon_Separation = temp.v.DeltaR(goodMu.at(j).v);
	if(Jet_Muon_Separation < 0.4) Well_Separated_Mu =false;
      }

      for(int j=0;j<(int)goodElectron.size();j++){
	float Jet_Electron_Separation = temp.v.DeltaR(goodElectron.at(j).v);
	if(Jet_Electron_Separation < 0.4) Well_Separated_Ele =false;
	
	}
      passCuts = passCuts &&  Well_Separated_Mu &&  Well_Separated_Ele;

      if(passCuts){
	goodJet.push_back(temp);
	
	if(Jet_btagDeepB[i]>0.4184){
	  good_bJet.push_back(temp);
	}
	njet++;
      }  // If template satisfies the conditions in cuts then jet is pushed back into goodJet
    }


    
    // sorting the goodJets
    Sort(3);



   
  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    h.nLepton-> Fill((int)goodLep.size());

    h.nElectron->Fill((int)goodElectron.size());
    
    //##############
    // Analysis
    //##############


    //Fill the size of goodMu array
    h.nmu->Fill((int)goodMu.size());


    //Plot pT of all muons in same plot
    for(int i=0; i<(int)goodMu.size(); i++){
      h.mupt->Fill(goodMu.at(i).v.Pt());
    }


    //Plotting the leading muon pT,Rta and relative isolation in each event.
    if((int)goodMu.size()>0){             
      h.muprop[0]->Fill(goodMu.at(0).v.Pt());
      h.muprop[1]->Fill(goodMu.at(0).v.Eta());
      h.muprop[2]->Fill(Muon_pfRelIso04_all[goodMu.at(0).ind]);
    }
    
    //Number of Leptons(Electrons +Muons)
    h.nLepton-> Fill((int)goodLep.size());
    h.nElectron->Fill((int)goodElectron.size());
    
    
    
    //Initializing some required variables		      

    LeptonSumPT1 = 0;
    LeptonSumPT2 = 0;
    LeptonSumPT3 = 0;
    
    JetSumPT1 = 0;
    JetSumPT2 = 0;
    JetSumPT3 = 0;	      
    



    // Events containing at least one Lepton
    if(goodLep.size()>0){
      if(goodLep.at(0).v.Pt()>25){
	Events1++;
	//h.nLepton[0]->Fill((int)goodLep.size());
	h.MET[0]->Fill(*PuppiMET_pt);
	
	//This for loop Sums up the Pts of goodLeps for all events containing at least one Lepton
	for(int i=0; i<(int)goodLep.size();i++){
	  float  Pt = goodLep.at(i).v.Pt();
	  LeptonSumPT1 = LeptonSumPT1 + Pt;
	}
	h.LT[0]->Fill(LeptonSumPT1);
	
	for(int i=0; i<(int)goodJet.size();i++){
	  float Pt =( goodJet.at(i).v.Pt());
	  JetSumPT1 = JetSumPT1 + Pt;
	}
	float E_Lep_leading = sqrt(pow(goodLep.at(0).v.M(),2) +pow(goodLep.at(0).v.Pt(),2));
	float dPhi_MET_leading = delta_phi(goodLep.at(0).v.Phi() ,*PuppiMET_phi);
	float Transv_Mass_lead = transv_mass(goodLep.at(0).v.Pt(),*PuppiMET_pt,dPhi_MET_leading);
	h.mT[0]->Fill(Transv_Mass_lead);
	
      }
      h.HT[0]->Fill(JetSumPT1);
    }
    


    //Events containing at least Two Leptons
    if((int)goodLep.size()>1){
      if(goodLep.at(0).v.Pt()>25 && goodLep.at(1).v.Pt()> 15){
	Events2++;
	


	//  h.nLepton[1]->Fill((int)goodLep.size());
	h.MET[1]->Fill(*PuppiMET_pt);
	


	//This for loop Sums up the Pts of goodLeps for all events containing at least Two Leptons
	for(int i=0; i<(int)goodLep.size();i++){
	  
	  float  Pt = goodLep.at(i).v.Pt();
	  LeptonSumPT2 = LeptonSumPT2 + Pt;
	}
	h.LT[1]->Fill(LeptonSumPT2);
	


	for(int i=0; i<(int)goodJet.size();i++){
	  float Pt = goodJet.at(i).v.Pt();
	  JetSumPT2 = JetSumPT2 + Pt;
	}
	

	float E_Lep_leading = sqrt(pow(goodLep.at(0).v.M(),2) +pow(goodLep.at(0).v.Pt(),2));
	float dPhi_MET_leading = delta_phi(goodLep.at(0).v.Phi() ,*PuppiMET_phi);
	float Transv_Mass_lead = transv_mass(goodLep.at(0).v.Pt(),*PuppiMET_pt,dPhi_MET_leading);
	h.mT[1]->Fill(Transv_Mass_lead);
	

	float E_Lep_subleading = sqrt(pow(goodLep.at(1).v.M(),2) +pow(goodLep.at(1).v.Pt(),2));
	float dPhi_MET_subleading = delta_phi(goodLep.at(1).v.Phi() ,*PuppiMET_phi);
	float Transv_Mass_sublead = transv_mass(goodLep.at(1).v.Pt(),*PuppiMET_pt,dPhi_MET_leading);
	h.mT[2]->Fill(Transv_Mass_sublead);


	float di_mass_0_1 =(goodLep.at(0).v + goodLep.at(1).v).M();
	h.Invar_Mass[0]->Fill( di_mass_0_1);
      }
       h.HT[1]->Fill(JetSumPT2);
    }
    
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////Events having at least Three Leptons/////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool METcont = *PuppiMET_pt<90 && *PuppiMET_pt>50;
    if((int)goodLep.size()>2){



      bool passCuts = goodLep.at(0).v.Pt()>25 && goodLep.at(1).v.Pt()> 15 && goodLep.at(2).v.Pt()>10;
      if(passCuts){

	Events3++;
	float HT = 0;
	float LT = 0;
	h.nLep[2]->Fill((int)goodLep.size());

	

	//Filling JetSumPT3
	for(int i=0; i<(int)goodJet.size();i++){
	  float Pt = goodJet.at(i).v.Pt();
	 
	  HT = HT+Pt;
	}
	
	//Filling Number of Jets	
	for(int i=0;i<(int)goodJet.size();i++){
	  h.Jet[0]->Fill(goodJet.at(i).v.Pt());
	}
       	
	//Filling MET
	h.MET[3]->Fill(*PuppiMET_phi);     
	//	if(LT>150 && LT<450 && *PuppiMET_pt>50 && *PuppiMET_pt<90 )


	//This for loop Sums up the Pts of goodLeps for all events containing at least Two Leptons
	for(int i=0; i<(int)goodLep.size();i++){
	  float  Pt = goodLep.at(i).v.Pt();
	  //LeptonSumPT3 = LeptonSumPT3 + Pt;
	  LT = LT + Pt;
	}
	

	// bool METcut =*PuppiMET_pt>33 && *PuppiMET_pt<91;
	// bool LTcut = LT<160 && LT>110;
	// bool HTcut = HT<100 && HT>30;
	// //	if(METcut && LTcut){// HTcut){
	  h.MET[2]->Fill(*PuppiMET_pt);
	
	  h.LT[2]->Fill(LT);
	  h.HT[2]->Fill(HT);
	  //	}
	//float E_Lep_leading = pow(pow(goodLep.at(0).v.M(),2) +pow(goodLep.at(0).v.Pt(),2),1/2);
	float dPhi_MET_lead = delta_phi( goodLep.at(0).v.Phi() ,*PuppiMET_phi);
	float Transv_Mass_lead = transv_mass(goodLep.at(0).v.Pt(),*PuppiMET_pt,dPhi_MET_lead);
	
	
	//float E_Lep_subleading = pow(pow(goodLep.at(1).v.M(),2) +pow(goodLep.at(1).v.Pt(),2),1/2);
	float dPhi_MET_sublead = delta_phi( goodLep.at(1).v.Phi(),*PuppiMET_phi);
	float Transv_Mass_sublead = transv_mass(goodLep.at(1).v.Pt(),*PuppiMET_pt,dPhi_MET_sublead);
       	
	//float E_Lep_ssubleading = pow(pow(goodLep.at(2).v.M(),2) + pow(goodLep.at(2).v.Pt(),2),1/2);
	float dPhi_MET_ssublead = delta_phi(goodLep.at(2).v.Phi(), *PuppiMET_phi);
	float Transv_Mass_ssublead = transv_mass(goodLep.at(2).v.Pt(), *PuppiMET_pt, dPhi_MET_ssublead);
	float max_transvmass =0;
	if((int)goodLep.size()>0 && (float)*PuppiMET_pt>0){
	  h.MET[4]->Fill(dPhi_MET_lead);
	  h.MET[5]->Fill(dPhi_MET_sublead);
	  h.MET[6]->Fill(dPhi_MET_ssublead);



	  //For loop filling transverse mass of pair giving maximum transverse mass in each event
	
	  float dPhi_MET =0;
	  float transverse_m =0;
	  
	  for(int i=0;i<3;i++){
	    dPhi_MET =delta_phi( goodLep.at(i).v.Phi(),*PuppiMET_phi);
	    transverse_m= transv_mass(goodLep.at(i).v.Pt(),*PuppiMET_pt,dPhi_MET);
	    if(transverse_m > max_transvmass){
	      max_transvmass= transverse_m;
	    }
	  }
	  h.mT[12]->Fill(max_transvmass);  
	}
	
	
	//Transverse Mass Plots
	bool mT1 = Transv_Mass_lead>60 && Transv_Mass_lead<100;
	bool mT2 = Transv_Mass_sublead>60 && Transv_Mass_sublead<120;
	bool mT3 = Transv_Mass_ssublead>70 && Transv_Mass_ssublead<145;
	float mT_c_root = pow(Transv_Mass_lead*Transv_Mass_sublead*Transv_Mass_ssublead,0.33333333333333);
	//	if( mT2 && mT3 && mT1 ){
 	h.mT[3]->Fill(Transv_Mass_lead);  
	h.mT[4]->Fill(Transv_Mass_sublead);
	h.mT[5]->Fill(Transv_Mass_ssublead);
	h.mT[6]->Fill(pow(Transv_Mass_lead,2));
	
	h.mT[7]->Fill(mT_c_root);
	float mT_avg =(Transv_Mass_lead+Transv_Mass_sublead+Transv_Mass_ssublead)/3;
	
	

	// if(zzCUT)
	h.dPhi_MET->Fill(dPhi_MET_lead);
	


	h.mT[8]->Fill(mT_avg);
	//	h.mT[9]->Fill(pow(Transv_Mass_lead,3));
	//	h.mT[10]->Fill(pow(Transv_Mass_lead,4));
	//	h.mT[11]->Fill(pow(Transv_Mass_lead,5));
	
	//	}	

	// The Delta Plots
	
 
	float d_eta =fabs(goodLep.at(0).v.Eta() - goodLep.at(1).v.Eta());
	h.Delta[0]->Fill(d_eta);
	
	float d_phi_01 =fabs(goodLep.at(0).v.Phi() - goodLep.at(1).v.Phi());
	if(d_phi_01>TMath::Pi()){d_phi_01 = 2*(TMath::Pi())-d_phi_01;}
	h.Delta[1]->Fill(d_phi_01);

	float delta_R= goodLep.at(0).v.DeltaR(goodLep.at(1).v);
	h.Delta[2]->Fill(delta_R);
	
	float d_phi_02 =fabs(goodLep.at(0).v.Phi() - goodLep.at(2).v.Phi());
	if(d_phi_02>TMath::Pi()){d_phi_02 = 2*(TMath::Pi())-d_phi_02;}
	h.Delta[3]->Fill(d_phi_02);
	
	float d_phi_12 =fabs(goodLep.at(1).v.Phi() - goodLep.at(2).v.Phi());
	if(d_phi_12>TMath::Pi()){d_phi_12 = 2*(TMath::Pi())-d_phi_12;}
	h.Delta[4]->Fill(d_phi_12);


	
	//Filling of the Invariant_Mass
	float di_mass_0_1 =(goodLep.at(0).v + goodLep.at(1).v).M();
	h.Invar_Mass[1]->Fill( di_mass_0_1);
	
	float di_mass_0_2 =(goodLep.at(0).v + goodLep.at(2).v).M();
	h.Invar_Mass[2]->Fill( di_mass_0_2);
	float di_mass_1_2 =(goodLep.at(1).v + goodLep.at(2).v).M();
	float Avg_invar_mass =(di_mass_0_1 +di_mass_0_2+ di_mass_1_2)/3;
	float GM_Invar_mass = pow(di_mass_0_1*di_mass_0_2*di_mass_1_2,0.333333333333);

	h.Invar_Mass[3]->Fill( di_mass_1_2);
	h.Invar_Mass[6]->Fill(pow(di_mass_0_1*di_mass_0_2*di_mass_1_2,0.333333333333));
	h.Invar_Mass[7]->Fill((di_mass_0_1 +di_mass_0_2+ di_mass_1_2)/3);


	






	//For loop filling invar mass of pairgiving maximum invariant mass in each event
	float max_invarmass =0;
	TLorentzVector ith_vector(0,0,0,0);
	TLorentzVector jth_vector(0,0,0,0);
	float Invarmass = 0;
	for(int i=0;i<3;i++)
	  {
	    ith_vector =goodLep.at(i).v;
	    for(int j=0;j<3;j++)
	      {
		jth_vector =goodLep.at(j).v;
	      }
	    Invarmass = (ith_vector+jth_vector).M();
	    if(Invarmass>max_invarmass){
	      max_invarmass = Invarmass;
	    }
	  }
	h.Invar_Mass[8]->Fill(max_invarmass) ;
	

	


	


	//btag_score
	btag_score =0;
	for(int i=0;i<(int)goodJet.size(); i++){
	  btag_score = btag_score + Jet_btagDeepB[goodJet.at(i).ind];
	}
	h.Jet[2]->Fill(btag_score);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Testing the combined cuts here
	bool Cut_ttz = btag_score>0.4 && HT>160;
	bool Cut_zz = *PuppiMET_pt>40 && mT_avg>60 && max_invarmass>100 ;
	//	bool Cut_ttz =	btag_score<0.4 && HT<160 && Avg_invar_mass>80 &&Avg_invar_mass<250;
	if(Cut_zz &&  Cut_ttz){
	  h.Jet[5]->Fill(goodJet.size());	     
	}
	// }
	//}

	float tri_mass = (goodLep.at(0).v+goodLep.at(1).v+goodLep.at(2).v).M();
	h.Invar_Mass[4]->Fill(tri_mass);

	//btag score individuals
	if((int)good_bJet.size()>0){
	  h.Jet[6]-> Fill(Jet_btagDeepB[good_bJet.at(0).ind]);
	  if((int)good_bJet.size()>1){ h.Jet[7]-> Fill(Jet_btagDeepB[good_bJet.at(1).ind]);}
	  
	  if((int)good_bJet.size()>2){ h.Jet[8]-> Fill(Jet_btagDeepB[good_bJet.at(2).ind]);}
	}
      }
     
    }
    if((int)goodLep.size()>3/* && METcont*/){
       bool passCuts = goodLep.at(0).v.Pt()>25 && goodLep.at(1).v.Pt()> 15 && goodLep.at(2).v.Pt()>10;
       if(passCuts){
       float tetra_mass= (goodLep.at(0).v+goodLep.at(1).v+goodLep.at(2).v+goodLep.at(3).v).M();
       h.Invar_Mass[5]->Fill(tetra_mass);
       }
    }
    // Events containing Min 3 Leptons
    //Events Containing at least three leptons

 
    // MET plot
    //dphi between MET and Lepton
   
    
    


//########### ANALYSIS ENDS HERE ##############
  }//GoodEvt


	return kTRUE;
}



