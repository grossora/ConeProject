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
//=========================================
//========== Bring in info  ===============
//=========================================

        auto mctruth = storage->get_data<event_mctruth>("generator");
        auto mcpart = mctruth->at(0).GetParticles();

		for( auto const& mcp : mcpart){
                        auto traj = mcp.Trajectory();
                        pionenergy = traj[0].E();
                        P_x = traj[0].X();
                        P_y = traj[0].Y();
                        P_z = traj[0].Z();
                        P_px = traj[0].Px();
                        P_py = traj[0].Py();
                        P_pz = traj[0].Pz();
                        P_pmag = sqrt(P_px*P_px+P_py*P_py+ P_pz*P_pz);
			}





	// Get this info from pi0 decay
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
		if(mcshower->size()==2){ // This means we have a gamma-gamma decay
			unsigned int showercount=0;
			for(auto const& mcs : *mcshower)
				{
					//Set Shower B
					if(showercount==1){
					auto ShowerDetProf =  mcs.DetProfile();
					auto gammamc = mcs.Start();
					ConePos_B = ShowerDetProf.Position();
					auto pos = ShowerDetProf.Position();
					ConeDir_B = ShowerDetProf.Momentum();
					auto dir = ShowerDetProf.Momentum();
					// Check if this cone is inside of the tpc. look at edges. 
					try{
				     	//Energy_B = ShowerDetProf.E();
				     	Energy_B = gammamc.E();
					//AxisLength_B = fconeprofile.Length(Energy_B);
					//OpeningAngle_B = fconeprofile.OpeningAngle(Energy_B)*180.0/PI;
					//AxisLength_B = fconeprofile.GammaLength(Energy_B);
					AxisLength_B = 15;
					OpeningAngle_B = 33.0;
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
					auto gammamc = mcs.Start();
					ConePos_A = ShowerDetProf.Position();
					auto pos = ShowerDetProf.Position();
					ConeDir_A = ShowerDetProf.Momentum();
					auto dir = ShowerDetProf.Momentum();
					// Check if this cone is inside of the tpc. look at edges. 
					try{
				     	//Energy_A = ShowerDetProf.E();
				     	Energy_A = gammamc.E();
					//AxisLength_A = fconeprofile.Length(Energy_A);
					//OpeningAngle_A = fconeprofile.OpeningAngle(Energy_A)*180.0/PI;
					//AxisLength_A = fconeprofile.GammaLength(Energy_A);
					AxisLength_A = 15;
					OpeningAngle_A = 33.0;
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
	if(coneintpc_A && coneintpc_B){
		for( unsigned int a = 0 ; a<nplanes; a++)
		{
		auto polyprojA = fgeoconic.ConicalFeatures(ConePos_A, ConeDir_A ,AxisLength_A , OpeningAngle_A, a ,smoothness);
		auto polyprojB = fgeoconic.ConicalFeatures(ConePos_B, ConeDir_B ,AxisLength_B , OpeningAngle_B, a ,smoothness);
		auto ABOverlap = fgeoconic.ConicalOverlap(polyprojA , polyprojB);
			if(ABOverlap) overlap = false;
		}
	}





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
                ConeTree->Branch("Energy_B",&Energy_B,"Energy_B/D");
                ConeTree->Branch("AxisLength_A",&AxisLength_A,"AxisLength_A/D");
                ConeTree->Branch("AxisLength_B",&AxisLength_B,"AxisLength_B/D");
                ConeTree->Branch("OpeningAngle_A",&OpeningAngle_A,"OpeningAngle_A/D");
                ConeTree->Branch("OpeningAngle_B",&OpeningAngle_B,"OpeningAngle_B/D");
                ConeTree->Branch("coneintpc_A",&coneintpc_A,"coneintpc_A/B");
                ConeTree->Branch("coneintpc_B",&coneintpc_B,"coneintpc_B/B");
                ConeTree->Branch("overlap",&overlap,"overlap/B");


                ConeTree->Branch("energy",&pionenergy,"pionenergy/D");
                ConeTree->Branch("PposX",&P_x,"P_x/D");
                ConeTree->Branch("PposY",&P_y,"P_y/D");
                ConeTree->Branch("PposZ",&P_z,"P_z/D");
                ConeTree->Branch("PdirX",&P_px,"P_px/D");
                ConeTree->Branch("PdirY",&P_py,"P_py/D");
                ConeTree->Branch("PdirZ",&P_pz,"P_pz/D");
                ConeTree->Branch("Ppmag",&P_pmag,"P_pmag/D");

	
		
	}// initialize tree






}
#endif
