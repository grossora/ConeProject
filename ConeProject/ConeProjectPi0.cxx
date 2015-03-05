#ifndef LARLITE_CONEPROJECTPI0_CXX
#define LARLITE_CONEPROJECTPI0_CXX

#include "ConeProjectPi0.h"
#include "DataFormat/cluster.h"


namespace larlite {

  bool ConeProjectPi0::initialize() {

	// Initialize tree and hists 
        InitializeAnaTree();

        cratio0 = new TH2D("cratio0","Ratio of Contained Charge / ; Length of Shower Axis; Charge: Contained/Deposit ",100,0,AxisLength,100,0,1);
        cratio1 = new TH2D("cratio1","Ratio of Contained Charge / ; Length of Shower Axis; Charge: Contained/Deposit",100,0,AxisLength,100,0,1);
        cratio2 = new TH2D("cratio2","Ratio of Contained Charge / ; Length of Shower Axis; Charge: Contained/Deposit",100,0,AxisLength,100,0,1);

        dangle0 = new TH1D("dangle0","bsaha",100,-1,1);
        dangle1 = new TH1D("dangle1","bsaha",100,-1,1);
        dangle2 = new TH1D("dangle2","bsaha",100,-1,1);

        dangledeg0 = new TH1D("dangledeg0","bsaha",100,-180,180);
        dangledeg1 = new TH1D("dangledeg1","bsaha",100,-180,180);
        dangledeg2 = new TH1D("dangledeg2","bsaha",100,-180,180);

   return true;
  }
  
  bool ConeProjectPi0::analyze(storage_manager* storage) {
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
	//== Truth Photon Start Position  
		TLorentzVector ConePos_A;
		TLorentzVector ConePos_B;
	//== Truth Photon Start Dir  
		TLorentzVector ConeDir_A;
		TLorentzVector ConeDir_B;
	//== If the cone is contained in TPC
		bool coneintpc_A;
		bool coneintpc_B;
        //== Make the pxhit  vector by plane for now... using hits
		std::vector<std::vector<std::vector<unsigned int>>> ContainedHitsVect_A(nplanes);
		std::vector<std::vector<std::vector<unsigned int>>> ContainedHitsVect_B(nplanes);
        //== Check if there are no overlap 
		bool ConesOverlapAreGood = true; 







	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> AxisSEpt(nplanes);
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::vector<larutil::PxPoint>> ConePolygonProjection(nplanes);
	//== Slope and Cept of the cone axis
		std::vector<std::pair<double,double>> ConeAxisSlopeCept;
	//== Slope and Cept of the surface axis
		std::vector<std::pair<double,double>> ConeSurfaceSlopeCept;
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<larutil::PxHit>> PxHitsVect(nplanes);
	//== RecoFit Based on weights 
		std::vector<std::pair<double,double>> recofitvec(nplanes);
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<larutil::PxHit>> ContainedPxHitsVect(nplanes);
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
					coneintpc_B = fgeoconic.ConeInTPC(pos,dir,AxisLength,openingangle, smoothness);
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
					coneintpc_A = fgeoconic.ConeInTPC(pos,dir,AxisLength,openingangle, smoothness);
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









//=========================================
//========== Define the output  ===========
//=========================================
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




//========================================================
//==== Look to see if there are overlap in the cones =====
//========================================================
	if(coneintpc_A && coneintpc_B){
		for( unsigned int a = 0 ; a<nplanes; a++)
		{
		auto polyprojA = fgeoconic.ConicalFeatures(ConePos_A, ConeDir_A ,AxisLength , openingangle, a ,smoothness);
		auto polyprojB = fgeoconic.ConicalFeatures(ConePos_B, ConeDir_B ,AxisLength , openingangle, a ,smoothness);
		auto Overlap = fgeoconic.ConicalOverlap(polyprojA , polyprojB);
			if(Overlap) ConesOverlapAreGood = false;
		}
	}
	
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define The Output  -------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================



















//========================================================
//==== Look to see What we can do with the overlap hits ==
//========================================================
	// This needs to be thought about... Many cases the cones will overlap
	if(coneintpc_A && coneintpc_B && !ConesOverlapAreGood){
	std::cout<<" There is IS overlap ZZZZZZZZZZZZZZZZZZZZ "<<std::endl;
		}


//=================================================
//========== Find hits that are in the cone  ======
//=================================================
	if(coneintpc_A && coneintpc_B &&ConesOverlapAreGood)
	{
		for( unsigned int a = 0 ; a<nplanes; a++)
		{
		auto polyprojA = fgeoconic.ConicalFeatures(ConePos_A, ConeDir_A ,AxisLength , openingangle, a ,smoothness);
		auto polyprojB = fgeoconic.ConicalFeatures(ConePos_B, ConeDir_B ,AxisLength , openingangle, a ,smoothness);
		std::vector<unsigned int> conth_A = fgeoconic.PolyContainHit(hitsvect, polyprojA, a);
		std::vector<unsigned int> conth_B = fgeoconic.PolyContainHit(hitsvect, polyprojB, a);
		ContainedHitsVect_A[a].push_back(conth_A);
		ContainedHitsVect_B[a].push_back(conth_B);
		}
	}//coneintpc(A and B)

//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------- Find hits that are in the cone ---------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================





// Let's see what we can do... Everything abover here gets me all the info that I need. 
	// Now I am just looking at somethings

//=================================================
//========== Fill Into Clusters  ======
//=================================================
// Lets push this into a cluster
   // set event variables
	AssSet_t  hit_ass_set;
	AssUnit_t hit_ass;

	if(coneintpc_A && coneintpc_B &&ConesOverlapAreGood){
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
        hit_ass.reserve(ContainedHitsVect_A[a][0].size());
        lite_cluster.clear_data();

        // Make association
        //for( auto const& hit : ContainedHitsVect[a][0] ) hit_ass.push_back(hit);// just thinking ahead
		//if no hit
        if(!ContainedHitsVect_A[a][0].size()) continue;
        if(!ContainedHitsVect_B[a][0].size()) continue;


        for( auto const& hit : ContainedHitsVect_A[a][0]) hit_ass.push_back(hit);

        // Add association
        hit_ass_set.push_back(hit_ass);

        // Set cluster ID
        lite_cluster.set_id(Output_cluster->size());

        // Set cluster view
        lite_cluster.set_view(hits->at(hit_ass.at(0)).View());

        // Add a cluster to the output
        Output_cluster->push_back(lite_cluster);

        hit_ass.clear();
        hit_ass.reserve(ContainedHitsVect_B[a][0].size());
        lite_cluster.clear_data();
        for( auto const& hit : ContainedHitsVect_B[a][0]) hit_ass.push_back(hit);
        hit_ass_set.push_back(hit_ass);
        lite_cluster.set_id(Output_cluster->size());
        lite_cluster.set_view(hits->at(hit_ass.at(0)).View());
        Output_cluster->push_back(lite_cluster);



	}// for loop over planes

		}// if cone in tpc

    Output_cluster->set_association(hits->id(),hit_ass_set);


    return true;
  }

  bool ConeProjectPi0::finalize() {

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


        void ConeProjectPi0::InitializeAnaTree()
        {
        FullTree = new TTree("FullTree","FullTree");
                //FullTree->Branch("ccnc",&ccnc,"ccnc/I");
	}// initialize tree






}
#endif
