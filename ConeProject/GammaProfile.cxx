#ifndef LARLITE_GAMMAPROFILE_CXX
#define LARLITE_GAMMAPROFILE_CXX

#include "GammaProfile.h"


namespace larlite {

  bool GammaProfile::initialize() {
/*
        fFracL0 = new TH2D("fLengthA0","DCharge/DVolume AKA DeDx ShowerProfile ;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);
        fFracL1 = new TH2D("fLengthA1","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);
        fFracL2 = new TH2D("fLengthA2","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1);
        f3FracL0 = new TH3D("f3LengthA0","DCharge/DVolume AKA DeDx ShowerProfile ;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1,100,0,600);
        f3FracL1 = new TH3D("f3LengthA1","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1,100,0,600);
        f3FracL2 = new TH3D("f3LengthA2","DCharge/DVolume AKA DeDx ShowerProfile;Length of Shower in CM; Cone Open Angle",200,0,200,50,0,0.1,100,0,600);
	*/
        fFracL0 = new TH2D("fLengthA0","dCharge/dX :DeDx ShowerProfile ;Length of Shower in CM; Energy Deposit/Length ",200,0,200,200,0,3000);
        fFracL1 = new TH2D("fLengthA1","dCharge/dx :DeDx ShowerProfile;Length of Shower in CM; Energy Deposit/Length",200,0,200,200,0,3000);
        fFracL2 = new TH2D("fLengthA2","dCharge/dx :DeDx ShowerProfile;Length of Shower in CM; Energy Deposit/Length",200,0,200,200,0,3000);
        f3FracL0 = new TH3D("f3LengthA0","dCharge/dx :DeDx ShowerProfile;Length of Shower in CM; Energy Deposit/Length",200,0,200,50,0,10000,1000,0,1000);
        f3FracL1 = new TH3D("f3LengthA1","dCharge/dx :DeDx ShowerProfile;Length of Shower in CM;Energy Deposit/Length",200,0,200,50,0,10000,1000,0,1000);
        f3FracL2 = new TH3D("f3LengthA2","dCharge/dx :DeDx ShowerProfile;Length of Shower in CM; Energy Deposit/Length",200,0,200,50,0,10000,1000,0,1000);

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
		TLorentzVector StartShowerPos;
	//== Truth Photon Start Dir  
		TLorentzVector StartConeDir;
		TLorentzVector StartShowerDir;
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
//========== Where are the gamma  =========
//=========================================
        auto geo = larutil::Geometry::GetME();
        auto x = geo->DetHalfWidth();
        //auto y = geo->DetHalfHeight();
        auto z = geo->DetLength();
        double xlo = x-50.0;
        double xhi = x+50.0;
        double ylo = -50.00;
        double yhi = 50.00;
        double zlo = z/2.0-50.00;
        double zhi = z/2.0+50.00;

	double ax = 0.1;
	double ay = 0.1;
	double az = 0.9;
	double dxl = ax-0.1; 
	double dxh = ax+0.1; 
	double dyl = ay-0.1; 
	double dyh = ay+0.1; 
	double dzl = az-0.1; 
	double dzh = az+0.1; 

	bool GoodPosition = false;
	bool GoodDirection = false;


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
				StartConePos = SP.Position();
				StartShowerPos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				StartConeDir = SP.Momentum();
				StartShowerDir = ShowerDetProf.Momentum();
				//StartShowerDir.SetPx(286.924);
				//StartShowerDir.SetPy(466.329);
				//StartShowerDir.SetPz(669.632);
				auto dir = ShowerDetProf.Momentum();
				energy = ShowerDetProf.E();
		if(pos.X()>xlo&&pos.X()<xhi&&pos.Y()>ylo&&pos.Y()<yhi&&pos.Z()>zlo&&pos.Z()<zhi){
			 GoodPosition = true;
			}
		double pmag = sqrt(dir.Px()*dir.Px()+dir.Py()*dir.Py()+dir.Pz()*dir.Pz());	
		//if(dir.Px()/pmag>dxl&&dir.Px()/pmag<dxh&&dir.Py()/pmag>dyl&&dir.Py()/pmag<dyh&&dir.Pz()/pmag>dzl&&dir.Pz()/pmag<dzh){
		if(dir.Pz()/pmag>0.9){
			 GoodDirection = true;
			}
		// Check if this cone is inside of the tpc. look at edges. 
	}//mcshower 
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================







//=--=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
		double angle = 60;
		//double angle = 32.94;
		bool coneintpc = true;
		double ConeLength = 100;	
		//double ConeLength = 120;	
		std::vector<double>  ChargeinCone(nplanes);
		//coneintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,ConeLength,angle, smoothness);
		bool showerintpc =  false;
		try{
		showerintpc = fgeoconic.ConeInTPC(StartShowerPos,StartShowerDir,ConeLength,angle, smoothness);
		showerintpc = fgeoconic.ConeInTPC(StartShowerPos,StartShowerDir,2,angle, smoothness);
		}catch(const ::larutil::LArUtilException& e){
		return false;
		showerintpc = false;
                }
		std::vector<double> ratio(nplanes);
	//if(coneintpc ) 
	//if(coneintpc && energy<200 && GoodDirection)
	//if(coneintpc && energy<200 && energy>150)
	//if(coneintpc &&  energy>610 && energy<620)
	//if(coneintpc  && energy>100)
	//if(showerintpc  && energy>150 && energy<200 && GoodDirection)
	//if(showerintpc && energy>150 && energy<200 )
	//if(coneintpc && energy>100 && energy<200 && GoodDirection )
	if(coneintpc && energy>100 && energy<105)
	//if(showerintpc && energy>200 && energy<400 )
	//if(showerintpc && energy>100 && energy<150 )
	//if(showerintpc && energy>100 && energy<200 && GoodDirection)
	//asdf
	{
			std::vector<double> DChargeDL(nplanes,0.0);
		for(unsigned int plane=0; plane<nplanes; plane++)
		{
			double DL = 2;
			//for(double a=5.0; a<ConeLength ; a+=DL)
			for(double a=2.0; a<ConeLength ; a+=DL)
			{
			//auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,a, angle, plane ,smoothness);
			std::vector<larutil::PxHit> contp;
			try{
			auto polyproj = fgeoconic.ConicalFeatures(StartShowerPos, StartShowerDir,a, angle, plane ,smoothness);
			contp = fgeoconic.PolyContain(PxHitsVect[plane], polyproj);
			}catch(const ::larutil::LArUtilException& e){
			continue;}
			// change the dir
				double newchargecont=0.0;
				for(auto const h: contp) newchargecont+=h.charge;
				double dCdX = newchargecont/a;
				if(dCdX<100 && a<14 ) continue;	
				if(plane==0) fFracL0->Fill(a,dCdX);
				if(plane==1) fFracL1->Fill(a,dCdX);
				if(plane==2) fFracL2->Fill(a,dCdX);
				if(plane==0) f3FracL0->Fill(a,dCdX,energy);
				if(plane==1) f3FracL1->Fill(a,dCdX,energy);
				if(plane==2) f3FracL2->Fill(a,dCdX,energy);
			}
		}
	}//


/*
////////Fill The Angles
//	if(coneintpc && GoodDirection ){
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
*/



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
	f3FracL0->Write();
	f3FracL1->Write();
	f3FracL2->Write();
	dangledeg0->Write();
	dangledeg1->Write();
	dangledeg2->Write();

    return true;
  }



    // void GammaProfile::InitializeAnaTree()
      //  {
       // FullTree = new TTree("FullTree","FullTree");
                //FullTree->Branch("Energy",&_Energy,"Energy/D");
                //FullTree->Branch("Vertexx",&_Vertexx,"Vertexx/D");
                //FullTree->Branch("Vertexy",&_Vertexy,"Vertexy/D");
                //FullTree->Branch("dcdxv",&dcdxv,"dcdxv/[D]");




}
#endif
