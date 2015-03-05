#ifndef LARLITE_CONEPROJECT_CXX
#define LARLITE_CONEPROJECT_CXX

#include "ConeProject.h"
#include "DataFormat/cluster.h"


namespace larlite {

  bool ConeProject::initialize() {

	// Initialize tree and hists 
        InitializeAnaTree();

        cratio0 = new TH2D("cratio0","Ratio of Contained Charge / Total (Open Angle 70 deg); Length of Shower Axis; Charge: Contained/Deposit ",100,0,AxisLength,100,0,1);
        cratio1 = new TH2D("cratio1","Ratio of Contained Charge / Total (Open Angle 70 deg); Length of Shower Axis; Charge: Contained/Deposit",100,0,AxisLength,100,0,1);
        cratio2 = new TH2D("cratio2","Ratio of Contained Charge / Total (Open Angle 70 deg); Length of Shower Axis; Charge: Contained/Deposit",100,0,AxisLength,100,0,1);

        dangle0 = new TH1D("dangle0","bsaha",100,-1,1);
        dangle1 = new TH1D("dangle1","bsaha",100,-1,1);
        dangle2 = new TH1D("dangle2","bsaha",100,-1,1);

        dangledeg0 = new TH1D("dangledeg0","bsaha",100,-180,180);
        dangledeg1 = new TH1D("dangledeg1","bsaha",100,-180,180);
        dangledeg2 = new TH1D("dangledeg2","bsaha",100,-180,180);

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
	//== 2D: This had the vertex points in two 2d 
		std::vector<larutil::PxPoint> vtx2d(nplanes);// Not used
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> AxisSEpt(nplanes);
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::vector<larutil::PxPoint>> ConePolygonProjection(nplanes);
	//== Slope and Cept of the cone axis
		std::vector<std::pair<double,double>> ConeAxisSlopeCept;
	//== Slope and Cept of the surface axis
		std::vector<std::pair<double,double>> ConeSurfaceSlopeCept;
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
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<larutil::PxHit>> ContainedPxHitsVect(nplanes);
        //== Make the pxhit  vector by plane for now... using hits
	//	std::vector<std::vector<std::vector<larlite::hit>>> ContainedHitsVect(nplanes);
        //== Make the pxhit  vector by plane for now... using hits
		std::vector<std::vector<std::vector<unsigned int>>> ContainedHitsVect(nplanes);

//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define some variables-----------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================


//=========================================
//========== Bring in info  ===============
//=========================================
// - Bring in Hits 
        auto hits = storage->get_data<event_hit>("gaushit");
        if(!hits){print(msg::kERROR,__FUNCTION__,"No DBCluster or associated hits found!");
                throw std::exception();
                return false;}
        // Fill the hits into hit vector 
        for(auto const& h : *hits) hitsvect.push_back(h);
                        for(auto const& hit : *hits){ 
                              ::larutil::PxHit h;
                              h.t = (hit.PeakTime() + tick_offset )* geom->TimeToCm() ;
                              h.w = hit.Wire()     * geom->WireToCm();
                              h.charge = hit.Charge();
                              h.peak   = hit.Charge(false);
                              h.plane  = hit.View();
                                if( (int)hit.View() ==0) PxHitsVect[0].push_back(h);
                                if( (int)hit.View() ==1) PxHitsVect[1].push_back(h);
                                if( (int)hit.View() ==2) PxHitsVect[2].push_back(h);
                                  }
//  This uses the truth info to get the cone
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
			for(auto const& mcs : *mcshower){
			// use Det Profile 
                                auto ShowerDetProf =  mcs.DetProfile();
				StartConePos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				StartConeDir = ShowerDetProf.Momentum();
				auto dir = ShowerDetProf.Momentum();
				std::cout<<"Det Prof Energy "<< ShowerDetProf.E()<<std::endl;
		// Check if this cone is inside of the tpc. look at edges. 
                coneintpc = fgeoconic.ConeInTPC(pos,dir,AxisLength,openingangle, smoothness);
	}//mcshower 
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================

//=========================================
//========== Define the output  ===========
//=========================================
    // First of all create an output
//      std::cout<<"Making it here ? "<<hits->event_id() <<std::endl;
		// use NC Filter for now. as a name holder
    auto Output_cluster = storage->get_data<event_cluster>("ncfilter");
    Output_cluster->clear_data();
    Output_cluster->set_event_id(hits->event_id());
    Output_cluster->set_run(hits->run());
    Output_cluster->set_subrun(hits->subrun());

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
		auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,AxisLength , openingangle, a ,smoothness);
	// $$2  Then find which hits are inside. 
		// $$ POLYCONTAIN returns the hits that are contained in that polygon
		std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
		ConePolygonProjection[a] = polyproj;
		ContainedPxHitsVect[a] = contp;

		//  this uses the hits
		//std::vector<larlite::hit> conth = fgeoconic.PolyContainHit(hitsvect, polyproj);
		std::vector<unsigned int> conth = fgeoconic.PolyContainHit(hitsvect, polyproj, a);
		ContainedHitsVect[a].push_back(conth);
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
		if(RandomLength<15) continue;
			
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
		if(Contchargebyplane[a]/chargebyplane[a]<1.0) ratio=Contchargebyplane[a]/chargebyplane[a];
		else{ ratio=1.0;}// Sometimes we will get greater than 1 due to rounding errors 1.02
		CRatio[a] = ratio;
			}
		// Now c ratio is ready to fill 
		cratio0->Fill(RandomLength,CRatio[0]);
		cratio1->Fill(RandomLength,CRatio[1]);
		cratio2->Fill(RandomLength,CRatio[2]);
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
			if((recofitvec[0].first>0.0 &&ConeAxisSlopeCept[0].first>0.0) || (recofitvec[0].first<0.0 &&ConeAxisSlopeCept[0].first<0.0) ){
			double dangled0 = atan(recofitvec[0].first) - atan(ConeAxisSlopeCept[0].first);
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}

			if(!((recofitvec[0].first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(recofitvec[0].first<0.0 && ConeAxisSlopeCept[0].first<0.0)) ){
			double dangled0 = atan(recofitvec[0].first) + atan(ConeAxisSlopeCept[0].first);
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}

			if((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0) || (recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0) ){
			double dangled1 = atan(recofitvec[1].first) - atan(ConeAxisSlopeCept[1].first);
			dangle1->Fill(dangled1);
			dangledeg1->Fill(dangled1*180/PI);
				}

			if(!((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0)||(recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0)) ){
			double dangled1 = atan(recofitvec[1].first) + atan(ConeAxisSlopeCept[1].first);
			dangle1->Fill(dangled1);
			dangledeg1->Fill(dangled1*180/PI);
				}

			if((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0) || (recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0) ){
			double dangled2 = atan(recofitvec[2].first) - atan(ConeAxisSlopeCept[2].first);
			dangle2->Fill(dangled2);
			dangledeg2->Fill(dangled2*180/PI);
				}

			if(!((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0)||(recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0)) ){
			double dangled2 = atan(recofitvec[2].first) + atan(ConeAxisSlopeCept[2].first);
			dangle2->Fill(dangled2);
			dangledeg2->Fill(dangled2*180/PI);
				}
	}// if cone in tpc


		//===============================================================================
		//======  LETS TRY TO LOOK AT Angle Fit based on Cone Length  ===================
		//===============================================================================

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
        hit_ass.clear();
        hit_ass.reserve(ContainedHitsVect[a][0].size());
        lite_cluster.clear_data();

        // Make association
        //for( auto const& hit : ContainedHitsVect[a][0] ) hit_ass.push_back(hit);// just thinking ahead
		//if no hit
        if(!ContainedHitsVect[a][0].size()) continue;


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
	////////////////

    return true;
  }


        void ConeProject::InitializeAnaTree()
        {
        FullTree = new TTree("FullTree","FullTree");
                //FullTree->Branch("ccnc",&ccnc,"ccnc/I");
	}// initialize tree






}
#endif
