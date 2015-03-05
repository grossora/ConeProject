#ifndef LARLITE_GAMMAPROFILE_CXX
#define LARLITE_GAMMAPROFILE_CXX

#include "GammaProfile.h"


namespace larlite {

  bool GammaProfile::initialize() {

        fFracL0 = new TH2D("fLengthA0","DCharge/DVolume AKA DeDx ShowerProfile ;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);
        fFracL1 = new TH2D("fLengthA1","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);
        fFracL2 = new TH2D("fLengthA2","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);


        dangledeg0 = new TH1D("dangledeg0","bsaha",201,-180,180);
        dangledeg1 = new TH1D("dangledeg1","bsaha",201,-180,180);
        dangledeg2 = new TH1D("dangledeg2","bsaha",201,-180,180);


   return true;
  }
  
  bool GammaProfile::analyze(storage_manager* storage) {
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
				auto SP = mcs.Start();
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

/*
//=============================================================================================================================
		//=================================================
		//== Theory Guess... 90% or better of charge     =============================
		//======  LETS TRY TO LOOK AT ENERGYProfile  ======
		//== Theory Guess... Params are angle = 32.94  : Based on 1 Moliere Radius  ==
		//== Look at rate of changing charge ==============
		//=================================================
		double angle = 60;
		//double angle = 32.94;
		bool coneintpc = true;
		double ConeLength = -999;
		std::vector<double>  ChargeinCone(nplanes);
	for(double length =5.0 ; length<150; length+=1){
		coneintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,length,angle, smoothness);
		std::vector<double> ratio(nplanes);
		if(!coneintpc){ coneintpc = false; break;}
			bool allplanes = true;
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,length, angle, a ,smoothness);
			std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);

			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;

			double chargecont=0;
			for(auto const h: contp) chargecont+=h.charge;
			ChargeinCone[a] = chargecont;

			double ratioplane = chargecont/chargeonplane;
			if(ratioplane<0.90){ allplanes = false; break;}
			//if(ratioplane<0.95){ allplanes = false; break;}
			}//loop over all the planes
	
			if(allplanes)
			{
			ConeLength = length;	
			break;
			}
		}// for loop over length
	*/	
//=--=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
		double angle = 60;
		//double angle = 32.94;
		bool coneintpc = true;
		double ConeLength = 80;	
		std::vector<double>  ChargeinCone(nplanes);
		coneintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,ConeLength,angle, smoothness);
		std::vector<double> ratio(nplanes);
		if(coneintpc){ 
			for( unsigned int a = 0 ; a<nplanes; a++){
			auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,ConeLength, angle, a ,smoothness);
			std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
			double chargeonplane=0;
			for(auto const h: PxHitsVect[a]) chargeonplane+=h.charge;
			double chargecont=0;
			for(auto const h: contp) chargecont+=h.charge;
			ChargeinCone[a] = chargecont;
			}//loop over all the planes
		}// cone in tpc
		
	//if(coneintpc)
	if(coneintpc && energy>250)
		{
			std::vector<double> DChargeDL(nplanes,0.0);
			for(unsigned int b=0; b<nplanes; b++)
			{
			auto pb = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,7.0, angle, b ,smoothness);
			std::vector<larutil::PxHit> cp = fgeoconic.PolyContain(PxHitsVect[b], pb);
			double chargecont=0.0;
                        for(auto const h: cp) chargecont+=h.charge;
			DChargeDL[b] = chargecont; 
			}

			for(unsigned int plane=0; plane<nplanes; plane++)
			{
				//double DL = 5.0;
				double DL = 7.0;
				for(double a=7.0; a<ConeLength ; a+=DL)
				{
				// is it true that as I move the radius the volume also grows proportially? 
			auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,a, angle, plane ,smoothness);
			std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[plane], polyproj);
				double newchargecont=0.0;
				for(auto const h: contp) newchargecont+=h.charge;
				double deltacharge;
				deltacharge = newchargecont-DChargeDL[plane]; 
				DChargeDL[plane] = newchargecont;
				// fill the deltacharge and the length; 
				double frac = deltacharge/DL/ChargeinCone[plane];
				if(plane==0) fFracL0->Fill(a,frac);
				if(plane==1) fFracL1->Fill(a,frac);
				if(plane==2) fFracL2->Fill(a,frac);
				}
			}

		}

////////Fill The Angles
                //if(ConeLength>5.0){
                if(coneintpc){
                        //This gets me the pair of points for start and end of the cone in 2D
                              AxisSEpt = fconeproj.startendpt(StartConePos,StartConeDir,ConeLength);
                        //This gets me the truth slope
                             ConeAxisSlopeCept = fconeproj.ConeAxisSC(StartConePos,StartConeDir,ConeLength);
                        // This Gets me the Fit of the balances charge
                        for(unsigned int a = 0 ; a < nplanes ; a++){
			auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,ConeLength, angle, a ,smoothness);
			std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
                        recofitvec[a] = fctools.ChargeDistFit(contp,ConeAxisSlopeCept[a],-1.0/ConeAxisSlopeCept[a].first);
                        }

                        if((recofitvec[0].first>0.0&&ConeAxisSlopeCept[0].first>0.0)||(recofitvec[0].first<0.0&&ConeAxisSlopeCept[0].first<0.0) ){
                        double dangled0 = atan(recofitvec[0].first) - atan(ConeAxisSlopeCept[0].first);
                        dangledeg0->Fill(dangled0*180/PI);
                                }

                        if(!((recofitvec[0].first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(recofitvec[0].first<0.0 && ConeAxisSlopeCept[0].first<0.0)) ){
                        double dangled0 = atan(recofitvec[0].first) + atan(ConeAxisSlopeCept[0].first);
                        dangledeg0->Fill(dangled0*180/PI);
                                }

                        if((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0) || (recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0) ){
                        double dangled1 = atan(recofitvec[1].first) - atan(ConeAxisSlopeCept[1].first);
                        dangledeg1->Fill(dangled1*180/PI);
                                }

                        if(!((recofitvec[1].first>0.0 && ConeAxisSlopeCept[1].first>0.0)||(recofitvec[1].first<0.0 && ConeAxisSlopeCept[1].first<0.0)) ){
                        double dangled1 = atan(recofitvec[1].first) + atan(ConeAxisSlopeCept[1].first);
                        dangledeg1->Fill(dangled1*180/PI);
                                }

                        if((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0) || (recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0) ){
                        double dangled2 = atan(recofitvec[2].first) - atan(ConeAxisSlopeCept[2].first);
                        dangledeg2->Fill(dangled2*180/PI);
                                }

                        if(!((recofitvec[2].first>0.0 && ConeAxisSlopeCept[2].first>0.0)||(recofitvec[2].first<0.0 && ConeAxisSlopeCept[2].first<0.0)) ){
                        double dangled2 = atan(recofitvec[2].first) + atan(ConeAxisSlopeCept[2].first);
                        dangledeg2->Fill(dangled2*180/PI);
                                }
        }// if cone in tpc




//=============================================================================================================================
		

		//==================================================================
		//======  Look at profile of contained 90% Contained Showers  ======
		//==================================================================





//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------End of 50 Contain-------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================




			

    return true;
  }

  bool GammaProfile::finalize() {

	if(_fout)
	_fout->cd();
	////////////////
	fFracL0->Write();
	fFracL1->Write();
	fFracL2->Write();
	dangledeg0->Write();
	dangledeg1->Write();
	dangledeg2->Write();

    return true;
  }






}
#endif
