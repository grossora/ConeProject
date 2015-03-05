#ifndef LARLITE_CONEPROJECT_CXX
#define LARLITE_CONEPROJECT_CXX

#include "ConeProject.h"

namespace larlite {

  bool ConeProject::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //
        InitializeAnaTree();


        dangle0 = new TH1D("dangle0","bsaha",100,-1,1);
        dangle1 = new TH1D("dangle1","bsaha",100,-1,1);
        dangle2 = new TH1D("dangle2","bsaha",100,-1,1);
        dangledeg0 = new TH1D("dangledeg0","bsaha",100,-180,180);
        dangledeg1 = new TH1D("dangledeg1","bsaha",100,-180,180);
        dangledeg2 = new TH1D("dangledeg2","bsaha",100,-180,180);

        chbal0 = new TH1D("chbal0","bsaha",1000,-0.3,0.3);
        chbal1 = new TH1D("chbal1","bsaha",1000,-0.3,0.3);
        chbal2 = new TH1D("chbal2","bsaha",1000,-0.3,0.3);

        DifTrueVal0 = new TH2D("Dif0","bsaha",50,0,14,500,-0.4,0.4);
        DifTrueVal1 = new TH2D("Dif1","bsaha",50,0,14,500,-0.4,0.4);
        DifTrueVal2 = new TH2D("Dif2","bsaha",50,0,14,500,-0.4,0.4);

        vertexdist0 = new TH1D("vdist0","bsaha",100,-10,10);
        vertexdist1 = new TH1D("vdist1","bsaha",100,-10,10);
        vertexdist2 = new TH1D("vdist2","bsaha",100,-10,10);

        aproj0 = new TH1D("aprodtruth0","bsaha",100,-1,1);
        aproj1 = new TH1D("aprodtruth1","bsaha",100,-1,1);
        aproj2 = new TH1D("aprodtruth2","bsaha",100,-1,1);


        Truedcosx01 = new TH1D("tdcx01","bsaha",100,-1,1);
        Truedcosy01 = new TH1D("tdcy01","bsaha",100,-1,1);
        Truedcosz01 = new TH1D("tdcz01","bsaha",100,-1,1);

        dcosx01 = new TH1D("dcx01","bsaha",100,-1,1);
        dcosy01 = new TH1D("dcy01","bsaha",100,-1,1);
        dcosz01 = new TH1D("dcz01","bsaha",100,-1,1);

        dphi01 = new TH1D("dphi01","bsaha",100,-180,180);
        dtheta01 = new TH1D("dtheta01","bsaha",100,-180,180);



    return true;
  }
  
  bool ConeProject::analyze(storage_manager* storage) {
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                        std::cout<<"Start of Event\n"<<std::endl;
//=================================================
//Define some things 
//=================================================
        auto geom = larutil::GeometryUtilities::GetME();
        auto tservice = larutil::TimeService::GetME();

	auto const& tpc_clock = tservice->TPCClock();
	double tick_offset = tservice->TriggerOffsetTPC() * tpc_clock.Frequency();

	//$$//$$ Define the Cone //
        double openingangle = 70.0; // magic number place holder for now. 
        double length = 14;// length of the cone axis
	int smoothness = 16;// would be nice if this was even
	//$$//$$ Define the Cone //

        unsigned int nplanes = larutil::Geometry::GetME()->Nplanes();
        std::vector<larlite::hit> hitsvect;
        // these things will be filled and used 
	//== 2D: This had the vertex points in two 2d 
		std::vector<larutil::PxPoint> vtx2d(nplanes);
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
		std::vector<larutil::PxHit> pxhits0;
		std::vector<larutil::PxHit> pxhits1;
		std::vector<larutil::PxHit> pxhits2;
	//== RecoFit Based on weights 
		std::pair<double,double> recofit0;
		std::pair<double,double> recofit1;
		std::pair<double,double> recofit2;
	//== Truth Based on weights 
		double TruthCB0;
		double TruthCB1;
		double TruthCB2;
	//== Truth Photon Start Position  
		TLorentzVector StartPhotPos;
		TLorentzVector StartConePos;
	//== Truth Photon Start Dir  
		TLorentzVector StartPhotDir;
		TLorentzVector StartConeDir;
	//== Truth Charge per Plane  
	std::vector<double> ChargePerPlane;
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
                                if( (int)hit.View() ==0) pxhits0.push_back(h);
                                if( (int)hit.View() ==1) pxhits1.push_back(h);
                                if( (int)hit.View() ==2) pxhits2.push_back(h);
                                  }
//		PxHitsVect.push_back(pxhits0);
//		PxHitsVect.push_back(pxhits1);
//		PxHitsVect.push_back(pxhits2);

// // This uses the truth info to get the cone
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
			for(auto const& mcs : *mcshower){
			// use Det Profile 
                                auto ShowerDetProf =  mcs.DetProfile();
				StartConePos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				StartConeDir = ShowerDetProf.Momentum();
				auto dir = ShowerDetProf.Momentum();
				ChargePerPlane = mcs.Charge();
			std::cout<<"Charge per plane 0 "<<ChargePerPlane[0]<<std::endl;
			std::cout<<"Charge per plane 1 "<<ChargePerPlane[1]<<std::endl;
			std::cout<<"Charge per plane 2 "<<ChargePerPlane[2]<<std::endl;
		// Check if this is inside of the tpc.look at edges. 
                coneintpc = fgeoconic.ConeInTPC(pos,dir,length,openingangle, smoothness);
	}// look 
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================






//=================================================
//========== Find hits that are in the cone  ======
//=================================================
	if(coneintpc){
// $$1  Need first find the cone polygon. 
// $$2  Then find which hits are inside. 
	//Loop over the planes 
	for( unsigned int a = 2 ; a<nplanes; a++){
	auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir, length, openingangle, a ,smoothness);
	//std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], polyproj);
	ConePolygonProjection[a] = polyproj;
//	ContainedPxHitsVect[a] = contp;
	}// loop over the planes
	
/*		
	for( unsigned int a = 0 ; a<nplanes; a++){
	std::vector<larutil::PxHit> contp = fgeoconic.PolyContain(PxHitsVect[a], ConePolygonProjection[a]);
	ContainedPxHitsVect[a] = contp;
	}
*/	

}// cone in tpc




			if(coneintpc){
// $$3  Need first find the cone polygon. 
			//This gets me the pair of points for start and end of the cone	in 2D
			      AxisSEpt = fconeproj.startendpt(StartConePos,StartConeDir,length);
			//This gets me the slope and cept for each cone axis in each plane	
			     //ConeAxisSlopeCept = fconeproj.ConeAxisSC(StartConePos,StartConeDir);
			     ConeAxisSlopeCept = fconeproj.ConeAxisSC(StartConePos,StartConeDir,length);
			// This Gets me the Truth Value of the balances charge
		         TruthCB0 = fctools.ChargeDistBalanceVal(pxhits0,ConeAxisSlopeCept[0],-1.0/ConeAxisSlopeCept[0].first);
		         TruthCB1 = fctools.ChargeDistBalanceVal(pxhits1,ConeAxisSlopeCept[1],-1.0/ConeAxisSlopeCept[1].first);
		         TruthCB2 = fctools.ChargeDistBalanceVal(pxhits2,ConeAxisSlopeCept[2],-1.0/ConeAxisSlopeCept[2].first);
			// This Gets me the Fit of the balances charge
                        recofit0 = fctools.ChargeDistFit(pxhits0,ConeAxisSlopeCept[0],-1.0/ConeAxisSlopeCept[0].first);
                        recofit1 = fctools.ChargeDistFit(pxhits1,ConeAxisSlopeCept[1],-1.0/ConeAxisSlopeCept[1].first);
                        recofit2 = fctools.ChargeDistFit(pxhits2,ConeAxisSlopeCept[2],-1.0/ConeAxisSlopeCept[2].first);

	// This is really ugly now because I wanted to get things to work
				chbal0->Fill(TruthCB0);
				chbal1->Fill(TruthCB1);
				chbal2->Fill(TruthCB2);

			// now look at the angle diff between the fit slope and the true slope
			if((recofit0.first>0.0 && ConeAxisSlopeCept[0].first>0.0) || (recofit0.first<0.0 && ConeAxisSlopeCept[0].first<0.0) ){
			double dangled0 = atan(recofit0.first) - atan(ConeAxisSlopeCept[0].first);
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}

			if(!((recofit0.first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(recofit0.first<0.0 && ConeAxisSlopeCept[0].first<0.0)) ){
			double dangled0 = atan(recofit0.first) + atan(ConeAxisSlopeCept[0].first);
			dangle0->Fill(dangled0);
			dangledeg0->Fill(dangled0*180/PI);
				}


			double dangled1 = atan(recofit1.first) - atan(ConeAxisSlopeCept[1].first);
			dangle1->Fill(dangled1);
			dangledeg1->Fill(dangled1*180/PI);
			double dangled2 = atan(recofit2.first) - atan(ConeAxisSlopeCept[2].first);
			dangle2->Fill(dangled2);
			dangledeg2->Fill(dangled2*180/PI);
	
double apr0 =atan(ConeAxisSlopeCept[0].first);
double apr1 =atan(ConeAxisSlopeCept[1].first);
double apr2 =atan(ConeAxisSlopeCept[2].first);
aproj0->Fill(apr0);
aproj1->Fill(apr1);
aproj2->Fill(apr2);

	}

	//////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///////////////////
//	Double_t dw0 = AxisSEpt[0].second.w - AxisSEpt[0].first.w;
//	Double_t dt0 = AxisSEpt[0].second.t - AxisSEpt[0].first.t;
//	Double_t dw1 = AxisSEpt[1].second.w - AxisSEpt[1].first.w;
//	Double_t dt1 = AxisSEpt[1].second.t - AxisSEpt[1].first.t;
	
//	Double_t omega0 = geom->Get2Dslope(dw0,dt0);
//	Double_t omega1 = geom->Get2Dslope(dw1,dt1);
	// using larutil::point
	//double omega0 = geom->Get2Dslope(AxisSEpt[0].second,AxisSEpt[0].first);
	//double omega1 = geom->Get2Dslope(AxisSEpt[1].second,AxisSEpt[1].first);
	
	
//	Double_t rdw0 = AxisSEpt[0].second.w - AxisSEpt[0].first.w;
//	Double_t rdt0 = recofit0.first*AxisSEpt[0].second.w+recofit0.second - recofit0.first*AxisSEpt[0].first.w+recofit0.second;
//	Double_t rdw1 = AxisSEpt[1].second.w - AxisSEpt[1].first.w;
//	Double_t rdt1 = recofit1.first*AxisSEpt[1].second.w+recofit1.second - recofit1.first*AxisSEpt[1].first.w+recofit1.second;
//	
//	Double_t romega0 = geom->Get2Dslope(rdw0,rdt0);
//	Double_t romega1 = geom->Get2Dslope(rdw1,rdt1);
	
/*
	// Stupid mistake.... not sure why....
    double xphi=0,xtheta=0;
    geom->Get3DaxisN(0,
                      1,
                      //omega0*TMath::Pi()/180.,
                      //omega1*TMath::Pi()/180.,
                      //atan(ConeAxisSlopeCept[0].first)*TMath::Pi()/180.,
                      //atan(ConeAxisSlopeCept[1].first)*TMath::Pi()/180.,
                      atan(ConeAxisSlopeCept[0].first),
                      atan(ConeAxisSlopeCept[1].first),
                      xphi,
                      xtheta);
		std::cout<<" phi "<< xphi <<" theta :" << xtheta <<std::endl;
    double tdirs[3]={0};
    geom->GetDirectionCosines(xphi,xtheta,tdirs);
    TVector3 truthdirs(tdirs[0],tdirs[1],tdirs[2]);
	std::cout<<"TRUTH dirsx "<<truthdirs.X()<<"dirsy "<<truthdirs.Y()<<"dirsz "<<truthdirs.Z()<<std::endl;
	
double rxphi=0,rxtheta=0;
    geom->Get3DaxisN(0,
                      1,
                      //recofit0.first*TMath::Pi()/180.,
                      //recofit1.first*TMath::Pi()/180.,
                      atan(recofit0.first)*TMath::Pi()/180.,
                      atan(recofit1.first)*TMath::Pi()/180.,
                      rxphi,
                      rxtheta);

    double dirs[3]={0};
    geom->GetDirectionCosines(rxphi,rxtheta,dirs);
    TVector3 vdirs(dirs[0],dirs[1],dirs[2]);


	std::cout<<"dirsx "<<vdirs.X()<<"dirsy "<<vdirs.Y()<<"dirsz "<<vdirs.Z()<<std::endl;
	std::cout<<" NORM  "<<vdirs.X()*vdirs.X()+vdirs.Y()*vdirs.Y()+vdirs.Z()*vdirs.Z()<<std::endl;

	// x 
	double delx01 = vdirs.X() - truthdirs.X();
	dcosx01->Fill(delx01);
	Truedcosx01->Fill(truthdirs.X());
	// y 
	double dely01 = vdirs.Y() - truthdirs.Y();
	dcosy01->Fill(dely01);
	Truedcosy01->Fill(truthdirs.Y());
	// z
	double delz01 = vdirs.Z() - truthdirs.Z();
	dcosz01->Fill(delz01);
	Truedcosz01->Fill(truthdirs.Z());


	// diff
	double dphi = xphi-rxphi;
	double dtheta = xtheta-rxtheta;
	dphi01->Fill(dphi);
	dtheta01->Fill(dtheta);

	//dphi01->Fill(xphi);
	//dtheta01->Fill(xtheta);


	//////%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///////////////////





	}// if in tpc


}// for loop


*/


/*

// Random draw away from truth vertex
for(unsigned int a =0; a<1; a++){
				double RandomLength = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomx= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomy= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomz= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double adx = randomx/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				double ady = randomy/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				double adz = randomz/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				TLorentzVector varypos; 
				varypos.SetX(StartPhotPos.X()+ ( RandomLength*PertBox*adx));
				varypos.SetY(StartPhotPos.Y()+ ( RandomLength*PertBox*ady));
				varypos.SetZ(StartPhotPos.Z()+ ( RandomLength*PertBox*adz));
				double ShiftDistance = sqrt( pow(RandomLength*PertBox*adx,2)+ pow(RandomLength*PertBox*ady,2)+pow(RandomLength*PertBox*adz,2));
                auto varyconeintpc = fgeoconic.ConeContained(varypos,StartPhotDir,openingangle,length);
			if(varyconeintpc && coneintpc ){
                           auto  varyConeAxisSlopeCept = fconeproj.ConeAxisSC(varypos,StartPhotDir);
			// the 2d fit
			auto VaryfitCB0 = fctools.ChargeDistFit(pxhits0,varyConeAxisSlopeCept[0],-1.0/varyConeAxisSlopeCept[0].first);
			auto VaryfitCB1 = fctools.ChargeDistFit(pxhits1,varyConeAxisSlopeCept[1],-1.0/varyConeAxisSlopeCept[1].first);
			auto VaryfitCB2 = fctools.ChargeDistFit(pxhits2,varyConeAxisSlopeCept[2],-1.0/varyConeAxisSlopeCept[2].first);



			if((VaryfitCB0.first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(VaryfitCB0.first<0.0 && ConeAxisSlopeCept[0].first<0.0)){
			double dangledif0 = atan(VaryfitCB0.first) - atan(ConeAxisSlopeCept[0].first);
			DifTrueVal0->Fill(ShiftDistance,dangledif0);
				}//	
			if(!((VaryfitCB0.first>0.0 && ConeAxisSlopeCept[0].first>0.0)||(VaryfitCB0.first<0.0& ConeAxisSlopeCept[0].first<0.0))){
			double dangledif0 = atan(VaryfitCB0.first) - atan(ConeAxisSlopeCept[0].first);
			DifTrueVal0->Fill(ShiftDistance,dangledif0);
				}//	


			double dangledif1 = atan(VaryfitCB1.first) - atan(ConeAxisSlopeCept[1].first);
			DifTrueVal1->Fill(ShiftDistance,dangledif1);
			double dangledif2 = atan(VaryfitCB2.first) - atan(ConeAxisSlopeCept[2].first);
			DifTrueVal2->Fill(ShiftDistance,dangledif2);
			}// if vary cone is in tpc

		}// looping over things




*/




















/*


	// Do the random draw in space
			for(auto const& mcs : *mcshower){
				auto photstart =  mcs.Start();
				auto ppos = photstart.Position();
				auto pdir = photstart.Momentum();
				auto varydir = photstart.Momentum();
				double PertBox = 14.0;// This is the box that position can wiggle in 

for(unsigned int a =0; a<20000; a++){
				double RandomLength = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomx= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomy= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomz= (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double adx = randomx/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				double ady = randomy/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				double adz = randomz/(sqrt(randomx*randomx+randomy*randomy+randomz*randomz)); 
				TLorentzVector varypos; 
				varypos.SetX(ppos.X()+ ( RandomLength*PertBox*adx));
				varypos.SetY(ppos.Y()+ ( RandomLength*PertBox*ady));
				varypos.SetZ(ppos.Z()+ ( RandomLength*PertBox*adz));
				double ShiftDistance = sqrt( pow(RandomLength*PertBox*adx,2)+ pow(RandomLength*PertBox*ady,2)+pow(RandomLength*PertBox*adz,2));
                auto varyconeintpc = fgeoconic.ConeContained(varypos,varydir,openingangle,length);
		if(varyconeintpc & coneintpc){
		//	std::cout<<" ARE WE IN THEEEEE VARYYYY"<<std::endl;
                           auto  varyConePoints = fconeproj.vertexboundsproj(varypos,varydir,openingangle,length);
                           auto  varyConeAxisSlopeCept = fconeproj.ConeAxisSC(varypos,varydir);
                           auto  varySurfaceSlope0 = fconeproj.perpslope2D(varyConePoints[0].second.first,varyConePoints[0].second.second);
                           auto  varySurfaceSlope1 = fconeproj.perpslope2D(varyConePoints[1].second.first,varyConePoints[1].second.second);
                           auto  varySurfaceSlope2 = fconeproj.perpslope2D(varyConePoints[2].second.first,varyConePoints[2].second.second);


			double VaryCB0 = fctools.ChargeDistBalanceVal(pxhits0,varyConeAxisSlopeCept[0],varySurfaceSlope0);
			double VaryCB1 = fctools.ChargeDistBalanceVal(pxhits1,varyConeAxisSlopeCept[1],varySurfaceSlope1);
			double VaryCB2 = fctools.ChargeDistBalanceVal(pxhits2,varyConeAxisSlopeCept[2],varySurfaceSlope2);
		
			// the 2d fit
			auto VaryfitCB0 = fctools.ChargeDistFit(pxhits0,varyConeAxisSlopeCept[0],varySurfaceSlope0);
			auto VaryfitCB1 = fctools.ChargeDistFit(pxhits1,varyConeAxisSlopeCept[1],varySurfaceSlope1);
			auto VaryfitCB2 = fctools.ChargeDistFit(pxhits2,varyConeAxisSlopeCept[2],varySurfaceSlope2);

				// Compare and fill 
				double DiffCB0 = TruthCB0 - VaryCB0;
				std::cout<< " Value of Diff CB0 "<<DiffCB0<<std::endl;
				std::cout<< " Value of ShiftDistance "<<ShiftDistance<<std::endl;
				double DiffCB1 = TruthCB1 - VaryCB1;
				double DiffCB2 = TruthCB2 - VaryCB2;
	
			//	DifTrueVal0->Fill(ShiftDistance,DiffCB0);
			//	DifTrueVal1->Fill(ShiftDistance,DiffCB1);
			//	DifTrueVal2->Fill(ShiftDistance,DiffCB2);
			//	dangle0->Fill(DiffCB0);
			//	dangle1->Fill(DiffCB1);
			//	dangle2->Fill(DiffCB2);


      double dangledif0 = atan(VaryfitCB0.first) - atan(ConeAxisSlopeCept[0].first);
        dangle0->Fill(dangledif0);
      double dangledif1 = atan(VaryfitCB1.first) - atan(ConeAxisSlopeCept[1].first);
        dangle1->Fill(dangledif1);
      double dangledif2 = atan(VaryfitCB2.first) - atan(ConeAxisSlopeCept[2].first);
        dangle2->Fill(dangledif2);

				DifTrueVal0->Fill(ShiftDistance,dangledif0);
				DifTrueVal1->Fill(ShiftDistance,dangledif1);
				DifTrueVal2->Fill(ShiftDistance,dangledif2);


				}// if vary incone

			}//	loop over the fill

}
*/


	/*
	//TRY BRING IN A PERTURBED VERSION
	        auto mcshower = storage->get_data<event_mcshower>("mcreco");
			for(auto const& mcs : *mcshower){
				auto photstart =  mcs.Start();
				auto ppos = photstart.Position();
				auto pdir = photstart.Momentum();
			
				TLorentzVector pos; 
				double PertBox = 20.0;// This is the box that position can wiggle in 
				double randomx = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomy = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randomz = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				//pos.SetZ(ppos.Z()+ PertBox*randomz);
				pos.SetX(ppos.X()+ +( ( (randomx/fabs(randomx))*5))+ PertBox*randomx);
				pos.SetY(ppos.Y()+ +( ( (randomy/fabs(randomy))*5))+ PertBox*randomy);
				pos.SetZ(ppos.Z()+ +( ( (randomz/fabs(randomz))*5))+ PertBox*randomz);
				TLorentzVector dir; 
				double PertPBox = 0.2;// This is the wiggle in direction direction ranges from 0-1 in cos0
				double randompx = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randompy = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				double randompz = (rand()/ double(RAND_MAX) -0.5 )*2.0;
				//double newpx = pdir.Px()+ PertPBox*randompx;
					//if(newpx>1.0||newpx<-1.0){newpx = pdir.Px()- PertPBox*randompx;}
				//double newpy = pdir.Py()+ PertPBox*randompy;	
				//	if(newpy>1.0||newpy<-1.0){newpy = pdir.Py()- PertPBox*randompy;}
				//double newpz = pdir.Pz()+ PertPBox*randompz;	
				//	if(newpz>1.0||newpz<-1.0){newpz = pdir.Pz()- PertPBox*randompz;}
				double newpx = pdir.Px() + (( (randomx/fabs(randomx))*0.2)+ PertPBox*randompx);
					if(newpx>1.0||newpx<-1.0){newpx = pdir.Px()-(( (randomx/fabs(randomx))*0.2)+ PertPBox*randompx);}
				double newpy = pdir.Py() + (( (randomy/fabs(randomy))*0.2)+ PertPBox*randompy);
					if(newpy>1.0||newpy<-1.0){newpy = pdir.Py()-(( (randomy/fabs(randomy))*0.2)+ PertPBox*randompy);}
				double newpz = pdir.Pz() + (( (randomz/fabs(randomz))*0.2)+ PertPBox*randompz);
					if(newpz>1.0||newpz<-1.0){newpz = pdir.Pz()-(( (randomz/fabs(randomz))*0.2)+ PertPBox*randompz);}
				dir.SetPx(newpx);
				dir.SetPy(newpy);
				dir.SetPz(newpz);
                coneintpc = fgeoconic.ConeContained(pos,dir,openingangle,length);
			if(coneintpc){
			     ConePoints = fconeproj.vertexboundsproj(pos,dir,openingangle,length);
			     ConeAxisSlopeCept = fconeproj.ConeAxisSC(pos,dir);
				 SurfaceSlope0 = fconeproj.perpslope2D(ConePoints[0].second.first,ConePoints[0].second.second);
				 SurfaceSlope1 = fconeproj.perpslope2D(ConePoints[1].second.first,ConePoints[1].second.second);
				 SurfaceSlope2 = fconeproj.perpslope2D(ConePoints[2].second.first,ConePoints[2].second.second);
				// BIG LOG... Make this a verbose mode
				std::cout<<" \n \t THIS IS THE START OF THE LOG::::: "<<std::endl;
				std::cout<<"========= PLANE 0============"<<std::endl;
				std::cout<<"Something I made Vertex wire"<< ConePoints[0].first.w<<std::endl;
				std::cout<<"Something I made Vertex time "<< ConePoints[0].first.t<<std::endl;
				std::cout<<"Cone Axis: time = " <<ConeAxisSlopeCept[0].first<<"* Wire "<<ConeAxisSlopeCept[0].second<<std::endl;
				std::cout<<"Surface slope : "<< SurfaceSlope0<<std::endl;
				std::cout<<"Something I made bounda w"<< ConePoints[0].second.first.w<<std::endl;
				std::cout<<"Something I made bounda t"<< ConePoints[0].second.first.t<<std::endl;
				std::cout<<"Something I made boundb w"<< ConePoints[0].second.second.w<<std::endl;
				std::cout<<"Something I made boundb t"<< ConePoints[0].second.second.t<<std::endl;

				std::cout<<"========= PLANE 1============"<<std::endl;
				std::cout<<"Something I made Vertex wire"<< ConePoints[1].first.w<<std::endl;
				std::cout<<"Something I made Vertex time "<< ConePoints[1].first.t<<std::endl;
				std::cout<<"Cone Axis: time = " <<ConeAxisSlopeCept[1].first<<"* Wire "<<ConeAxisSlopeCept[1].second<<std::endl;
				std::cout<<"Surface slope : "<< SurfaceSlope1<<std::endl;
				std::cout<<"Something I made bounda w"<< ConePoints[1].second.first.w<<std::endl;
				std::cout<<"Something I made bounda t"<< ConePoints[1].second.first.t<<std::endl;
				std::cout<<"Something I made boundb w"<< ConePoints[1].second.second.w<<std::endl;
				std::cout<<"Something I made boundb t"<< ConePoints[1].second.second.t<<std::endl;

				std::cout<<"========= PLANE 2============"<<std::endl;
				std::cout<<"Something I made Vertex wire"<< ConePoints[2].first.w<<std::endl;
				std::cout<<"Something I made Vertex time "<< ConePoints[2].first.t<<std::endl;
				std::cout<<"Cone Axis: time = " <<ConeAxisSlopeCept[2].first<<"* Wire "<<ConeAxisSlopeCept[2].second<<std::endl;
				std::cout<<"Surface slope : "<< SurfaceSlope2<<std::endl;
				std::cout<<"Something I made bounda w"<< ConePoints[2].second.first.w<<std::endl;
				std::cout<<"Something I made bounda t"<< ConePoints[2].second.first.t<<std::endl;
				std::cout<<"Something I made boundb w"<< ConePoints[2].second.second.w<<std::endl;
				std::cout<<"Something I made boundb t"<< ConePoints[2].second.second.t<<std::endl;

				std::cout<<"\t THIS IS THE END OF THE LOG::::: \n"<<std::endl;
							}// if cone in tpc
						}// Loop over the mcshower
	
*/

//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//---------------Bring in Info---------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================


//=================================================
//-- Do some stuff 
//=================================================
std::cout<<" ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl;
/*	
	
	if(coneintpc){
	recofit0 = fctools.ChargeDistFit(pxhits0,ConeAxisSlopeCept[0],SurfaceSlope0);
		std::cout<<"RecoCone Axis: time = " <<recofit0.first<<"* Wire "<<recofit0.second -249.925<<std::endl;
	recofit1 = fctools.ChargeDistFit(pxhits1,ConeAxisSlopeCept[1],SurfaceSlope1);
		std::cout<<"RecoCone Axis: time = " <<recofit1.first<<"* Wire "<<recofit1.second -249.925<<std::endl;
	recofit2 = fctools.ChargeDistFit(pxhits2,ConeAxisSlopeCept[2],SurfaceSlope2);
		std::cout<<"RecoCone Axis: time = " <<recofit2.first<<"* Wire "<<recofit2.second -249.925<<std::endl;
			}// if cone in tpc
*/
//=================================================
//-- Comare things 
//=================================================
		// look at how well we did
/*
	if(coneintpc){
      double dangledif0 = atan(recofit0.first) - atan(ConeAxisSlopeCept[0].first);
	dangle0->Fill(dangledif0);
      double dangledif1 = atan(recofit1.first) - atan(ConeAxisSlopeCept[1].first);
	dangle1->Fill(dangledif1);
      double dangledif2 = atan(recofit2.first) - atan(ConeAxisSlopeCept[2].first);
	dangle2->Fill(dangledif2);

	double cb0 = fctools.ChargeDistBalanceVal(pxhits0,ConeAxisSlopeCept[0],SurfaceSlope0);
	double cb1 = fctools.ChargeDistBalanceVal(pxhits1,ConeAxisSlopeCept[1],SurfaceSlope1);
	double cb2 = fctools.ChargeDistBalanceVal(pxhits2,ConeAxisSlopeCept[2],SurfaceSlope2);

	chbal0->Fill(cb0);
	chbal1->Fill(cb1);
	chbal2->Fill(cb2);
*/	

		//}
























  
    return true;
  }

  bool ConeProject::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //

	if(_fout)
	_fout->cd();

	DifTrueVal0->Write();
	DifTrueVal1->Write();
	DifTrueVal2->Write();
	//FullTree->Write();
	dangle0->Write();
	dangle1->Write();
	dangle2->Write();
  
	chbal0->Write();
	chbal1->Write();
	chbal2->Write();


	vertexdist0->Write();
	vertexdist1->Write();
	vertexdist2->Write();


	aproj0->Write();
	aproj1->Write();
	aproj2->Write();

	dangledeg0->Write();
	dangledeg1->Write();
	dangledeg2->Write();


	Truedcosx01->Write();
	Truedcosy01->Write();
	Truedcosz01->Write();

	dcosx01->Write();
	dcosy01->Write();
	dcosz01->Write();


	dphi01->Write();
	dtheta01->Write();

    return true;
  }


        void ConeProject::InitializeAnaTree()
        {
        FullTree = new TTree("FullTree","FullTree");
                //FullTree->Branch("ccnc",&ccnc,"ccnc/I");
	}// initialize tree






}
#endif
