#ifndef RECOTOOL_CONICALPROJECTION_CXX
#define RECOTOOL_CONICALPROJECTION_CXX


#include "conicalprojection.h"

namespace larlite {




//-----------------------------------------------------------------------------------------------------------------
        // This is the pxpoint for the vertex?
	std::vector<larutil::PxPoint> conicalprojection::vertexproj(TLorentzVector& pos){
	// Has to have good inputs... or put up a flag
        auto geom = larutil::GeometryUtilities::GetME();
	std::vector<larutil::PxPoint> vertex;
		auto p0 =geom->Get2DPointProjectionCM(&pos,0);
		vertex.push_back(p0);
		auto p1 =geom->Get2DPointProjectionCM(&pos,1);
		vertex.push_back(p1);
		auto p2 =geom->Get2DPointProjectionCM(&pos,2);
		vertex.push_back(p2);
	return vertex;
	}

//-----------------------------------------------------------------------------------------------------------------
        // This is a pair of PxPoints for the start and end point?
	std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> conicalprojection::startendpt(TLorentzVector& pos,TLorentzVector& dir, double Length){
	// Has to have good inputs... or put up a flag
        auto geom = larutil::GeometryUtilities::GetME();
	std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> Vsep;
		auto p0 =geom->Get2DPointProjectionCM(&pos,0);
		auto p1 =geom->Get2DPointProjectionCM(&pos,1);
		auto p2 =geom->Get2DPointProjectionCM(&pos,2);
	double mag =  sqrt(dir.Px()* dir.Px()+dir.Py()* dir.Py()+dir.Pz()* dir.Pz());
        double vsx = dir.Px()/mag;
        double vsy = dir.Py()/mag;
        double vsz = dir.Pz()/mag;
        TLorentzVector vv;
        double xtpc = pos.X() + Length*vsx;
        vv.SetX(xtpc);
        double ytpc = pos.Y() + Length*vsy;
        vv.SetY(ytpc);
        double ztpc = pos.Z() + Length*vsz;
        vv.SetZ(ztpc);
			// Need to know that these are inside of the TPC
		auto e0 =geom->Get2DPointProjectionCM(&vv,0);
		auto e1 =geom->Get2DPointProjectionCM(&vv,1);
		auto e2 =geom->Get2DPointProjectionCM(&vv,2);
	// Make the pairs 
	std::pair<larutil::PxPoint,larutil::PxPoint> se0(p0,e0);
	std::pair<larutil::PxPoint,larutil::PxPoint> se1(p1,e1);
	std::pair<larutil::PxPoint,larutil::PxPoint> se2(p2,e2);
		Vsep.push_back(se0);
		Vsep.push_back(se1);
		Vsep.push_back(se2);
	return Vsep;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------------------------
        // This is a pair that contains the slope and cept for the axis in cm space 
	std::vector<std::pair<double,double>> conicalprojection::ConeAxisSC(TLorentzVector& pos, TLorentzVector& dir, double length){
		std::vector<std::pair<double,double>> ret;
        auto geom = larutil::GeometryUtilities::GetME();
		auto v0 =geom->Get2DPointProjectionCM(&pos,0);
		auto v1 =geom->Get2DPointProjectionCM(&pos,1);
		auto v2 =geom->Get2DPointProjectionCM(&pos,2);
        double mag =  sqrt(dir.Px()* dir.Px()+dir.Py()* dir.Py()+dir.Pz()* dir.Pz());
        double vsx = dir.Px()/mag;
        double vsy = dir.Py()/mag;
        double vsz = dir.Pz()/mag;
        TLorentzVector vv;
	// avoid infinities 
        double xtpc = pos.X() + length*vsx;
        vv.SetX(xtpc);
        double ytpc = pos.Y() + length*vsy;
        vv.SetY(ytpc);
        double ztpc = pos.Z() + length*vsz;
        vv.SetZ(ztpc);
		auto vv0 =geom->Get2DPointProjectionCM(&vv,0);
		auto vv1 =geom->Get2DPointProjectionCM(&vv,1);
		auto vv2 =geom->Get2DPointProjectionCM(&vv,2);
                  double slope0 = (v0.t -vv0.t)/(v0.w - vv0.w);
                  double cept0 = v0.t - slope0 * v0.w ; //   is there a time issue? 
		std::pair<double,double> sc0(slope0,cept0);
		ret.push_back(sc0);	
                  double slope1 = (v1.t -vv1.t)/(v1.w - vv1.w);
                  double cept1 = v1.t - slope1 * v1.w; //   is there a time issue? 
		std::pair<double,double> sc1(slope1,cept1);
		ret.push_back(sc1);	
                  double slope2 = (v2.t -vv2.t)/(v2.w - vv2.w);
                  double cept2 = v2.t - slope2 * v2.w; //   is there a time issue? 
		std::pair<double,double> sc2(slope2,cept2);
		ret.push_back(sc2);	
	return ret;
	}



}

#endif

