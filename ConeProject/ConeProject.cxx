#ifndef LARLITE_CONEPROJECT_CXX
#define LARLITE_CONEPROJECT_CXX

#include "ConeProject.h"
#include "DataFormat/cluster.h"


namespace larlite {

  bool ConeProject::initialize() {

	// Initialize tree and hists 
        InitializeAnaTree();
        cratio0 = new TH2D("cratio0","Ratio of Contained Charge / Total ; Length of Shower Axis; Charge: Contained/Deposit ",200,0,150,200,0,1);
        cratio1 = new TH2D("cratio1","Ratio of Contained Charge / Total ; Length of Shower Axis; Charge: Contained/Deposit",200,0,150,200,0,1);
        cratio2 = new TH2D("cratio2","Ratio of Contained Charge / Total ; Length of Shower Axis; Charge: Contained/Deposit",200,0,150,200,0,1);

        dangle0 = new TH1D("dangle0","bsaha",100,-1,1);
        dangle1 = new TH1D("dangle1","bsaha",100,-1,1);
        dangle2 = new TH1D("dangle2","bsaha",100,-1,1);

        dangledeg0 = new TH1D("dangledeg0","bsaha",201,-180,180);
        dangledeg1 = new TH1D("dangledeg1","bsaha",201,-180,180);
        dangledeg2 = new TH1D("dangledeg2","bsaha",201,-180,180);


   return true;
  }
  
  bool ConeProject::analyze(storage_manager* storage) {
//=================================================
//Define some things 
//=================================================
        auto geom = larutil::GeometryUtilities::GetME();
        auto tservice = larutil::TimeService::GetME();
	auto const& tpc_clock = tservice->TPCClock();
	double tick_offset = tservice->TriggerOffsetTPC() * tpc_clock.Frequency();// Kazu's fix
        unsigned int nplanes = larutil::Geometry::GetME()->Nplanes();

        // these things will be filled and used 
	//== LarLight vector of hits
		std::vector<larlite::hit> hitsvect;
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> AxisSEpt(nplanes);
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::vector<larutil::PxPoint>> ConePolygonProjection(nplanes);
	//== Slope and Cept of the cone axis
		std::vector<std::pair<double,double>> ConeAxisSlopeCept;
	//== If the cone is contained in TPC
		bool coneintpc =false;
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<larutil::PxHit>> PxHitsVect(nplanes);
	//== RecoFit Based on weights 
		std::vector<std::pair<double,double>> recofitvec(nplanes);
	//== Truth Photon Start Position  
		TLorentzVector StartConePos;
	//== Truth Photon Start Dir  
		TLorentzVector StartConeDir;
        //== Make the pxhit  vector by plane for now... using hits
		std::vector<std::vector<larutil::PxHit>> ContainedPxHitsVect(nplanes);
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<unsigned int>> ContainedHitsVect(nplanes);
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define some variables-----------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================






//=========================================
//========== Bring in info  ===============
//=========================================
        auto hits = storage->get_data<event_hit>("gaushit");
        if(!hits){print(msg::kERROR,__FUNCTION__,"No DBCluster or associated hits found!");
                throw std::exception();
                return false;}
        for(auto const& h : *hits) hitsvect.push_back(h);
                        for(auto const& hit : *hits){ 
                              ::larutil::PxHit h;
                              h.t = (hit.PeakTime() + tick_offset )* geom->TimeToCm() ;
                              h.w = hit.WireID().Wire     * geom->WireToCm();
                              h.charge = hit.Integral();
                              h.plane  = hit.View();
                                if( (int)hit.View() ==0) PxHitsVect[0].push_back(h);
                                if( (int)hit.View() ==1) PxHitsVect[1].push_back(h);
                                if( (int)hit.View() ==2) PxHitsVect[2].push_back(h);
                                  }
//  This uses the truth info to get the cone
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
			for(auto const& mcs : *mcshower){
                                auto ShowerDetProf =  mcs.DetProfile();
				StartConePos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				StartConeDir = ShowerDetProf.Momentum();
				auto dir = ShowerDetProf.Momentum();
		// Check if this cone is inside of the tpc. look at edges. 
	try{
                coneintpc = fgeoconic.ConeInTPC(pos,dir,AxisLength,openingangle, smoothness);
	}catch(const ::larutil::LArUtilException& e){
        continue;
		}
			//Info for the tree 
				_Energy=ShowerDetProf.E();
		auto ShowerLength = fconeprofile.Length(_Energy);
//		auto ShowerOpenAngle = fconeprofile.OpeningAngle(_Energy);
		//AxisLength = ShowerLength;
		//openingangle = ShowerOpenAngle*180/PI;
		AxisLength = 50;
		openingangle = 60;
	
				_Vertexx = pos.X();
				_Vertexy = pos.Y();
				_Vertexz = pos.Z();
	}//mcshower 


//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================





//=========================================
//========== Define the output  ===========
//=========================================
    //auto Output_cluster = storage->get_data<event_cluster>("ncfilter");
    //auto Output_cluster = storage->get_data<event_cluster>("ShowerCone");
    //Output_cluster->clear_data();
    //Output_cluster->set_event_id(hits->event_id());
    //Output_cluster->set_id();
    //Output_cluster->set_run(hits->run());
    //Output_cluster->set_subrun(hits->subrun());
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define The Output  -------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================







//=================================================
//========== Find hits that are in the cone  ======
//=================================================
	if(coneintpc){
	//Loop over the planes 
		for( unsigned int a = 0 ; a<nplanes; a++){
	// $$1  Need first find the cone polygon. 
		// $$ CONICALFEATURES returns the 2d polygon for the 3d cone volume
	try{
		auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,AxisLength , openingangle, a ,smoothness);
		std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
		ConePolygonProjection[a] = polyproj;
		ContainedPxHitsVect[a] = contp;
		std::vector<unsigned int> conth = fgeoconic.PolyContainHit(hitsvect, polyproj, a);
		ContainedHitsVect[a] = conth;
		}catch(const ::larutil::LArUtilException& e)
			{
		ConePolygonProjection[a].clear();
		ContainedPxHitsVect[a].clear();
		ContainedHitsVect[a].clear();
			}
		}// loop over the planes
	}// cone in tpc

//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------- Find hits that are in the cone ---------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================





// Let's see what we can do... Everything abover here gets me all the info that I need. 
	// Now I am just looking at somethings
		//=================================================
		//======  LETS TRY TO LOOK AT ENERGYCONTAIN  ======
		//=================================================

		// This is just playing around
	if(coneintpc){
		for(unsigned int r=0 ; r<RDraw; r++){
		// Random Length
			tryagain:
                                double RandomL = (rand()/ double(RAND_MAX));
		double RandomLength = RandomL*AxisLength;
		                bool subconeintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,RandomLength,openingangle, smoothness);
		if(subconeintpc){

		// small cones cause problems somtimes
		//if(RandomLength<15) continue;
		if(RandomLength<5) continue;
			
		for( unsigned int a = 0 ; a<nplanes; a++){
		auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir, RandomLength, openingangle, a ,smoothness);
			if(polyproj.size()==1){
			// This means things are bad this probably happens with small cones
			std::cout<<"getting out of here"<<std::endl;
			goto tryagain;}
	
		std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
		ConePolygonProjection[a] = polyproj;
		ContainedPxHitsVect[a] = contp;
		}// loop over the planes
		// First What is the total charge on each plane	
		std::vector<double> chargebyplane(nplanes);
		for(unsigned int a=0; a<nplanes ; a++){
			double charge=0;
			for(auto const h: PxHitsVect[a]) charge+=h.charge;
			chargebyplane[a] = charge;
			}
		// Now I can look at the contained fraction that get caught in the contained polygon
		std::vector<double> Contchargebyplane(nplanes);
		for(unsigned int a=0; a<nplanes ; a++){
			double charge=0;
			for(auto const h: ContainedPxHitsVect[a]) charge+=h.charge;
			Contchargebyplane[a] = charge;
			}

		std::vector<double> CRatio(nplanes);
		for(unsigned int a=0; a<nplanes ; a++){
		double ratio;
		ratio= Contchargebyplane[a]/chargebyplane[a];
		CRatio[a] = ratio;
			}
		// Now c ratio is ready to fill 
		cratio0->Fill(RandomLength,CRatio[0]);
		cratio1->Fill(RandomLength,CRatio[1]);
		cratio2->Fill(RandomLength,CRatio[2]);
			//Fill the cratio tree
		FullTree->Fill();
		_CRatio0=CRatio[0];
		_CRatio1=CRatio[1];
		_CRatio2=CRatio[2];
		_RandomLength = RandomLength;
		CRTree->Fill();
		}// if in subtpc

		}
	}// in cone

	// Let's keep playing
		//======================================================================================
		//======  LETS TRY TO LOOK Reco Angle Based on simple chargewt/dist linear fit =========
		//======================================================================================
		if(coneintpc){
			//This gets me the pair of points for start and end of the cone	in 2D
			      AxisSEpt = fconeproj.startendpt(StartConePos,StartConeDir,AxisLength);
			//This gets me the truth slope
			     ConeAxisSlopeCept = fconeproj.ConeAxisSC(StartConePos,StartConeDir,AxisLength);
			// This Gets me the Fit of the balances charge
			for(unsigned int a = 0 ; a < nplanes ; a++){
                        recofitvec[a] = fctools.ChargeDistFit(PxHitsVect[a],ConeAxisSlopeCept[a],-1.0/ConeAxisSlopeCept[a].first);
			}
			
			if((recofitvec[0].first>0.0&&ConeAxisSlopeCept[0].first>0.0)||(recofitvec[0].first<0.0&&ConeAxisSlopeCept[0].first<0.0) ){
			double dangled0 = atan(recofitvec[0].first) - atan(ConeAxisSlopeCept[0].first);
			_DangleDeg0 = dangled0*180/PI;
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}

			if(!((recofitvec[0].first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(recofitvec[0].first<0.0 && ConeAxisSlopeCept[0].first<0.0)) ){
			double dangled0 = atan(recofitvec[0].first) + atan(ConeAxisSlopeCept[0].first);
			_DangleDeg0 = dangled0*180/PI;
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}

			if((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0) || (recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0) ){
			double dangled1 = atan(recofitvec[1].first) - atan(ConeAxisSlopeCept[1].first);
			_DangleDeg1 = dangled1*180/PI;
			dangle1->Fill(dangled1);
			dangledeg1->Fill(dangled1*180/PI);
				}

			if(!((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0)||(recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0)) ){
			double dangled1 = atan(recofitvec[1].first) + atan(ConeAxisSlopeCept[1].first);
			_DangleDeg1 = dangled1*180/PI;
			dangle1->Fill(dangled1);
			dangledeg1->Fill(dangled1*180/PI);
				}

			if((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0) || (recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0) ){
			double dangled2 = atan(recofitvec[2].first) - atan(ConeAxisSlopeCept[2].first);
			_DangleDeg2 = dangled2*180/PI;
			dangle2->Fill(dangled2);
			dangledeg2->Fill(dangled2*180/PI);
				}

			if(!((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0)||(recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0)) ){
			double dangled2 = atan(recofitvec[2].first) + atan(ConeAxisSlopeCept[2].first);
			_DangleDeg2 = dangled2*180/PI;
			dangle2->Fill(dangled2);
			dangledeg2->Fill(dangled2*180/PI);
				}
			DATree->Fill();
	}// if cone in tpc

		//===============================================================================
		//======  LETS TRY TO LOOK AT Fit DCos based on contained hits  =================
		//===============================================================================







/*

// Lets push this into a cluster
   // set event variables
	AssSet_t  hit_ass_set;
	AssUnit_t hit_ass;

	if(coneintpc){
		// loop over all the planes 
	for(unsigned int a=0; a<nplanes;a++)
	{
	// if cone in tpc
	
	//
	// Save the clusters 
	//
        ::larlite::cluster lite_cluster;
		// skip of there are no hits. 
		//ContainedHitsVect
    // Clear data products
        if(!ContainedHitsVect[a][0].size()) continue;
        hit_ass.clear();
        hit_ass.reserve(ContainedHitsVect[a][0].size());
        lite_cluster.clear_data();

        // Make association
        //for( auto const& hit : ContainedHitsVect[a][0] ) hit_ass.push_back(hit);// just thinking ahead
		//if no hit


        for( auto const& hit : ContainedHitsVect[a][0]) hit_ass.push_back(hit);

        // Add association
        hit_ass_set.push_back(hit_ass);

        // Set cluster ID
        lite_cluster.set_id(Output_cluster->size());

        // Set cluster view
        lite_cluster.set_view(hits->at(hit_ass.at(0)).View());


        // Add a cluster to the output
        Output_cluster->push_back(lite_cluster);
	}// for loop over planes

		}// if cone in tpc
    Output_cluster->set_association(hits->id(),hit_ass_set);
*/



    return true;
  }

  bool ConeProject::finalize() {

	if(_fout)
	_fout->cd();

	cratio0->Write();
	cratio1->Write();
	cratio2->Write();

	dangle0->Write();
	dangle1->Write();
	dangle2->Write();

	dangledeg0->Write();
	dangledeg1->Write();
	dangledeg2->Write();
	
	FullTree->Write();
	CRTree->Write();
	DATree->Write();
	////////////////

    return true;
  }


        void ConeProject::InitializeAnaTree()
        {
        FullTree = new TTree("FullTree","FullTree");
                FullTree->Branch("Energy",&_Energy,"Energy/D");
        CRTree = new TTree("CRTree","CRTree");
                CRTree->Branch("Vertexx",&_Vertexx,"Vertexx/D");
                CRTree->Branch("Vertexy",&_Vertexy,"Vertexy/D");
                CRTree->Branch("Vertexz",&_Vertexz,"Vertexz/D");
                CRTree->Branch("Energy",&_Energy,"Energy/D");
                CRTree->Branch("CRatio0",&_CRatio0,"CRatio0/D");
                CRTree->Branch("CRatio1",&_CRatio1,"CRatio1/D");
                CRTree->Branch("CRatio2",&_CRatio2,"CRatio2/D");
                CRTree->Branch("RandomLength",&_RandomLength,"RandomLength/D");

        DATree = new TTree("DATree","DATree");
                DATree->Branch("Vertexx",&_Vertexx,"Vertexx/D");
                DATree->Branch("Vertexy",&_Vertexy,"Vertexy/D");
                DATree->Branch("Vertexz",&_Vertexz,"Vertexz/D");
                DATree->Branch("Energy",&_Energy,"Energy/D");
                DATree->Branch("DD0",&_DangleDeg0,"DD0/D");
                DATree->Branch("DD1",&_DangleDeg1,"DD1/D");
                DATree->Branch("DD2",&_DangleDeg2,"DD2/D");
	}// initialize tree






}
#endif
