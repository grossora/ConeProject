#ifndef LARLITE_PI0CONTAINED_CXX
#define LARLITE_PI0CONTAINED_CXX

#include "Pi0Contained.h"


namespace larlite {

  bool Pi0Contained::initialize() {

	// Initialize tree and hists 
        InitializeAnaTree();

   return true;
  }
  
  bool Pi0Contained::analyze(storage_manager* storage) {
//=================================================
//Define some things 
//=================================================
        unsigned int nplanes = larutil::Geometry::GetME()->Nplanes();
	//== Truth Photon Start Position  
		TLorentzVector ConePos_A;
		TLorentzVector ConePos_B;
	//== Truth Photon Start Dir  
		TLorentzVector ConeDir_A;
		TLorentzVector ConeDir_B;
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define some variables-----------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================




//=========================================
//========== Bring in info  ===============
//=========================================
//        auto hits = storage->get_data<event_hit>("gaushit");
  //      if(!hits){print(msg::kERROR,__FUNCTION__,"No DBCluster or associated hits found!");
    //            throw std::exception();
      //          return false;}
//  This uses the truth info to get the cone
		// Get this info from pi0 decay
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
		if(mcshower->size()==2){ // This means we have a gamma-gamma decay
			unsigned int showercount=0;
			for(auto const& mcs : *mcshower)
				{
					//Set Shower B
					if(showercount==1){
					auto ShowerDetProf =  mcs.DetProfile();
					ConePos_B = ShowerDetProf.Position();
					auto pos = ShowerDetProf.Position();
					ConeDir_B = ShowerDetProf.Momentum();
					auto dir = ShowerDetProf.Momentum();
					// Check if this cone is inside of the tpc. look at edges. 
					try{
				     	Energy_B = ShowerDetProf.E();
					AxisLength_B = fconeprofile.Length(Energy_B);
					OpeningAngle_B = fconeprofile.OpeningAngle(Energy_B)*180.0/PI;
					coneintpc_B = fgeoconic.ConeInTPC(pos,dir,AxisLength_B,OpeningAngle_B, smoothness);
					}catch(const ::larutil::LArUtilException& e){
					coneintpc_B = false;
					continue;
					}
					break;
					}// if showercount==0
					//Set Shower A
					if(showercount==0){
					auto ShowerDetProf =  mcs.DetProfile();
					ConePos_A = ShowerDetProf.Position();
					auto pos = ShowerDetProf.Position();
					ConeDir_A = ShowerDetProf.Momentum();
					auto dir = ShowerDetProf.Momentum();
					// Check if this cone is inside of the tpc. look at edges. 
					try{
				     	Energy_A = ShowerDetProf.E();
					AxisLength_A = fconeprofile.Length(Energy_A);
					OpeningAngle_A = fconeprofile.OpeningAngle(Energy_A)*180.0/PI;
					coneintpc_A = fgeoconic.ConeInTPC(pos,dir,AxisLength_A,OpeningAngle_A, smoothness);
					}catch(const ::larutil::LArUtilException& e){
					coneintpc_A = false;
					continue;
					}
					showercount++;
					}// if showercount==0
				}//mcshower 
			}//if mcshower==2
		if(mcshower->size()!=2){// This means a bad mcshower or dalitz decay	
		 coneintpc_A = false;
		 coneintpc_B = false;}
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================

//========================================================
//==== Look to see if there are overlap in the cones =====
//========================================================
		overlap = true;
		// Definition is odd.... if return true that means that there is not overlap and cones are good
		// if false then that means that there is overlap and don't proceed
	if(coneintpc_A && coneintpc_B){
		for( unsigned int a = 0 ; a<nplanes; a++)
		{
		auto polyprojA = fgeoconic.ConicalFeatures(ConePos_A, ConeDir_A ,AxisLength_A , OpeningAngle_A, a ,smoothness);
		auto polyprojB = fgeoconic.ConicalFeatures(ConePos_B, ConeDir_B ,AxisLength_B , OpeningAngle_B, a ,smoothness);
		auto ABOverlap = fgeoconic.ConicalOverlap(polyprojA , polyprojB);
			if(ABOverlap) overlap = false;
		}
	}
	
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------:Looking for overlap  ----------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================


//========================================================
//==== Fill Tree for Non-Dalitz Decay ============== =====
//========================================================
	if(mcshower->size()==2) ConeTree->Fill();// This means a  mcshower is not a dalitz decay    
		

	

    return true;
  }

  bool Pi0Contained::finalize() {

	if(_fout)
	_fout->cd();
	ConeTree->Write();

    return true;
  }


        void Pi0Contained::InitializeAnaTree()
        {
        ConeTree = new TTree("ConeTree","ConeTree");
                ConeTree->Branch("Energy_A",&Energy_A,"Energy_A/D");
                ConeTree->Branch("Energy_B",&Energy_B,"EShower_B/D");
                ConeTree->Branch("AxisLength_A",&AxisLength_A,"AxisLength_A/D");
                ConeTree->Branch("AxisLength_B",&AxisLength_B,"AxisLength_B/D");
                ConeTree->Branch("OpeningAngle_A",&OpeningAngle_A,"OpeningAngle_A/D");
                ConeTree->Branch("OpeningAngle_B",&OpeningAngle_B,"OpeningAngle_B/D");
                ConeTree->Branch("coneintpc_A",&coneintpc_A,"coneintpc_A/B");
                ConeTree->Branch("coneintpc_B",&coneintpc_B,"coneintpc_B/B");
                ConeTree->Branch("overlap",&overlap,"overlap/B");
	
		
	}// initialize tree






}
#endif
