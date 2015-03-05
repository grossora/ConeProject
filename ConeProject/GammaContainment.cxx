#ifndef LARLITE_GAMMACONATINMENT_CXX
#define LARLITE_GAMMACONATINMENT_CXX

#include "GammaContainment.h"


namespace larlite {

  bool GammaContainment::initialize() {

        fLength = new TH1D("fLength","bsaha",100,0,200);
        fLengthA = new TH2D("fLengthA","bsaha;Length of Shower in CM; Cone Open Angle",100,0,100,90,0,90);


        C0L50 = new TH1D("C0L50","bsaha",100,0,1.0);
        C1L50 = new TH1D("C1L50","bsaha",100,0,1.0);
        C2L50 = new TH1D("C2L50","bsaha",100,0,1.0);

        C0L70 = new TH1D("C0L70","bsaha",100,0,1.0);
        C1L70 = new TH1D("C1L70","bsaha",100,0,1.0);
        C2L70 = new TH1D("C2L70","bsaha",100,0,1.0);

        C0L100 = new TH1D("C0L100","bsaha",100,0,1.0);
        C1L100 = new TH1D("C1L100","bsaha",100,0,1.0);
        C2L100 = new TH1D("C2L100","bsaha",100,0,1.0);

        C0L140 = new TH1D("C0L140","bsaha",100,0,1.0);
        C1L140 = new TH1D("C1L140","bsaha",100,0,1.0);
        C2L140 = new TH1D("C2L140","bsaha",100,0,1.0);

        CTheory0 = new TH1D("CTheory0","bsaha",100,0,1.0);
        CTheory1 = new TH1D("CTheory1","bsaha",100,0,1.0);
        CTheory2 = new TH1D("CTheory2","bsaha",100,0,1.0);
   return true;
  }
  
  bool GammaContainment::analyze(storage_manager* storage) {
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
                                auto ShowerDetProf =  mcs.DetProfile();
				StartConePos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				StartConeDir = ShowerDetProf.Momentum();
				auto dir = ShowerDetProf.Momentum();
		// Check if this cone is inside of the tpc. look at edges. 
	}//mcshower 


//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================

		//============================================================================
		//======  LETS TRY TO LOOK AT ENERGYCONTAIN  =================================
		//== Theory Guess... Params are angle = 32.94  : Based on 1 Moliere Radius  ==
		//== Theory Guess 2...Params are angle = 61.2  : Based on 2 Moliere Radius  ==
		//== Theory Guess... Params are Length = 34.16cm =============================
		//== Theory Guess... 90% or better of charge     =============================
		//============================================================================

		//double angleT = 32.94;
		double angleT = 61.2;
		bool coneintpcT = fgeoconic.ConeInTPC(StartConePos,StartConeDir,34.94, angleT, smoothness);
		std::vector<double> ratioT(nplanes);
		if(coneintpcT){
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyprojT = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,34.94, angleT, a ,smoothness);
			std::vector<larutil::PxHit> contpT = fgeoconic.PolyContain(PxHitsVect[a], polyprojT);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contpT) chargecont+=h.charge;
		
			double ratioplane = chargecont/chargeonplane;
			ratioT[a] = ratioplane;
			}
		CTheory0->Fill(ratioT[0]);
		CTheory1->Fill(ratioT[1]);
		CTheory2->Fill(ratioT[2]);
		}// if contained






//=============================================================================================================================
		//=================================================
		//======  LETS TRY TO LOOK AT ENERGYCONTAIN  ======
		//== Step through the lengths untill getting .95 ==
		//=================================================
	for(double ang=10; ang<90.0; ang+=1){
		double angle = ang;
		bool coneintpc = true;
	for(double length =10.0 ; length<100; length+=1){
		coneintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,length,angle, smoothness);
		std::vector<double> ratio(nplanes);
		if(!coneintpc){ coneintpc = false; break;}
			bool allplanes95 = true;
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,length, angle, a ,smoothness);
			std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp) chargecont+=h.charge;

			double ratioplane = chargecont/chargeonplane;
			if(ratioplane<0.90){ allplanes95 = false; break;}
			}//loop over all the planes
		
			//if(allplanes95){fLength->Fill(length); break;}
			if(allplanes95){fLengthA->Fill(length,ang); break;}
		}// for loop over length

	}//
//=============================================================================================================================
		

		//==================================================================
		//======  Look at profile of contained 90% Contained Showers  ======
		//==================================================================




		//=================================================
		//======  LETS TRY TO LOOK AT ENERGYCONTAIN  ======
		//=================================================
		//double angle50 = 2*atan(fconeprofile.ShowerRadius()/50.0)*180.0/PI;
		double angle50 = 60.0;
		bool coneintpc50 = fgeoconic.ConeInTPC(StartConePos,StartConeDir,50.0,angle50, smoothness);
		std::vector<double> ratio50(nplanes);
		if(coneintpc50){
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj50 = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,50.0, angle50, a ,smoothness);
			std::vector<larutil::PxHit> contp50 = fgeoconic.PolyContain(PxHitsVect[a], polyproj50);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp50) chargecont+=h.charge;
		
			double ratioplane = chargecont/chargeonplane;
			if(ratioplane>1) { std::cout<<" Ratio " <<chargecont/chargeonplane<<std::endl;
					   std::cout<<" PxHitsVect a  " <<PxHitsVect[a].size()<<std::endl;
					   std::cout<<" Contain   " <<contp50.size()<<std::endl;
						}//
			ratio50[a] = ratioplane;
			}
		C0L50->Fill(ratio50[0]);
		C1L50->Fill(ratio50[1]);
		C2L50->Fill(ratio50[2]);
		}// if contained



		double angle70 = 60.0;
		bool coneintpc70 = fgeoconic.ConeInTPC(StartConePos,StartConeDir,70.0,angle70, smoothness);
		std::vector<double> ratio70(nplanes);
		if(coneintpc70){
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj70 = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,70.0, angle70, a ,smoothness);
			std::vector<larutil::PxHit> contp70 = fgeoconic.PolyContain(PxHitsVect[a], polyproj70);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp70) chargecont+=h.charge;
		
			double ratioplane = chargecont/chargeonplane;
			if(ratioplane>1) { std::cout<<" Ratio " <<chargecont/chargeonplane<<std::endl;
					   std::cout<<" PxHitsVect a  " <<PxHitsVect[a].size()<<std::endl;
					   std::cout<<" Contain   " <<contp70.size()<<std::endl;
						}//
			ratio70[a] = ratioplane;
			}
		C0L70->Fill(ratio70[0]);
		C1L70->Fill(ratio70[1]);
		C2L70->Fill(ratio70[2]);
		}// if contained



		double angle100 = 60.0;
		bool coneintpc100 = fgeoconic.ConeInTPC(StartConePos,StartConeDir,100.0,angle100, smoothness);
		std::vector<double> ratio100(nplanes);
		if(coneintpc100){
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj100 = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,100.0, angle100, a ,smoothness);
			std::vector<larutil::PxHit> contp100 = fgeoconic.PolyContain(PxHitsVect[a], polyproj100);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp100) chargecont+=h.charge;
		
			double ratioplane = chargecont/chargeonplane;
			if(ratioplane>1) { std::cout<<" Ratio " <<chargecont/chargeonplane<<std::endl;
					   std::cout<<" PxHitsVect a  " <<PxHitsVect[a].size()<<std::endl;
					   std::cout<<" Contain   " <<contp100.size()<<std::endl;
						}//
			ratio100[a] = ratioplane;
			}
		C0L100->Fill(ratio100[0]);
		C1L100->Fill(ratio100[1]);
		C2L100->Fill(ratio100[2]);
		}// if contained


		double angle140 = 60.0;
		bool coneintpc140 = fgeoconic.ConeInTPC(StartConePos,StartConeDir,140.0,angle140, smoothness);
		std::vector<double> ratio140(nplanes);
		if(coneintpc140){
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj140 = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,140.0, angle140, a ,smoothness);
			std::vector<larutil::PxHit> contp140 = fgeoconic.PolyContain(PxHitsVect[a], polyproj140);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp140) chargecont+=h.charge;
		
			double ratioplane = chargecont/chargeonplane;
			if(ratioplane>1) { std::cout<<" Ratio " <<chargecont/chargeonplane<<std::endl;
					   std::cout<<" PxHitsVect a  " <<PxHitsVect[a].size()<<std::endl;
					   std::cout<<" Contain   " <<contp140.size()<<std::endl;
						}//
			ratio140[a] = ratioplane;
			}
		C0L140->Fill(ratio140[0]);
		C1L140->Fill(ratio140[1]);
		C2L140->Fill(ratio140[2]);
		}// if contained




//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------End of 50 Contain-------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================




			


    return true;
  }

  bool GammaContainment::finalize() {

	if(_fout)
	_fout->cd();
	////////////////
	fLength->Write();
	fLengthA->Write();

	C0L50->Write();
	C1L50->Write();
	C2L50->Write();

	C0L70->Write();
	C1L70->Write();
	C2L70->Write();

	C0L100->Write();
	C1L100->Write();
	C2L100->Write();

	C0L140->Write();
	C1L140->Write();
	C2L140->Write();

	CTheory0->Write();
	CTheory1->Write();
	CTheory2->Write();

    return true;
  }






}
#endif
