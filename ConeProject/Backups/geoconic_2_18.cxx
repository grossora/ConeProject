#ifndef RECOTOOL_GEOCONIC_CXX
#define RECOTOOL_GEOCONIC_CXX


#include "geoconic.h"

namespace larlite {

	

//-----------------------------------------------------------------------------------------------------------
	// Is the 3d point in the tpc?
bool geoconic::TPCContained(const TLorentzVector& pos  ){
        auto geo = larutil::Geometry::GetME();
	auto x = geo->DetHalfWidth();
	auto y = geo->DetHalfHeight();
	auto z = geo->DetLength();
		double fid = 1; // This needs to be here because we are not doing the position box properly 
	double xlo = 0+fid;
	double xhi = 2*x-fid;
	double ylo = -y+fid;
	double yhi = y-fid;
	double zlo = 0+fid;
	double zhi = z-fid;
        if(pos.X()>xhi || pos.X()<xlo || pos.Y()>yhi || pos.Y()<ylo || pos.Z()<zlo || pos.Z()>zhi){
		return  false;}
	return true;
	}//geoconic::TPCContained

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------------------
	// Are the 3d Cone Edges Contained in the tpc?
bool geoconic::ConeInTPC(const TLorentzVector& Pos, const TLorentzVector& dir, double Length, double OpeningAngle, int smoothness){
        // Calculate the radius at a given length and opening angle 
        double Radius = Length*tan(OpeningAngle/2 *PI/180);
        auto geo = larutil::Geometry::GetME();
	// Get the proper bounds for tpc 
        auto x = geo->DetHalfWidth();
        auto y = geo->DetHalfHeight();
        auto z = geo->DetLength();
		double fid = 1; // This needs to be here because we are not doing the position box properly 
	double xlo = 0+fid;
	double xhi = 2*x-fid;
	double ylo = -y+fid;
	double yhi = y-fid;
	double zlo = 0+fid;
	double zhi = z-fid;

	// Make the points 
                TLorentzVector pos = Pos;
                TVector3 Axis;
                double axisnorm = sqrt(dir.Px()*dir.Px()+dir.Py()*dir.Py()+dir.Pz()*dir.Pz());
                Axis.SetX(dir.Px()/axisnorm);
                Axis.SetY(dir.Py()/axisnorm);
                Axis.SetZ(dir.Pz()/axisnorm);
                // Get orthogonal vector
                TVector3 avec;
                avec = Axis.Orthogonal();
                // cross these vects to get the other mut-orthoganal one
                TVector3 bvec;
                bvec = Axis.Cross(avec);

	// If vertex point is not in tpc return false
        if(Pos.X()>xhi || Pos.X()<xlo || Pos.Y()>yhi || Pos.Y()<ylo || Pos.Z()<zlo || Pos.Z()>zhi){
        return false;}
		
  for(double theta =0 ; theta<2*PI; theta+=PI/smoothness){
                double paramX = Pos.X() + Length*Axis.X()/Axis.Mag() + Radius*cos(theta)*avec.X() + Radius*sin(theta)*bvec.X();
                double paramY = Pos.Y() + Length*Axis.Y()/Axis.Mag() + Radius*cos(theta)*avec.Y() + Radius*sin(theta)*bvec.Y();
                double paramZ = Pos.Z() + Length*Axis.Z()/Axis.Mag() + Radius*cos(theta)*avec.Z() + Radius*sin(theta)*bvec.Z();
	// if any of the cone edges are not in tpc return false 
        if(paramX>xhi || paramX<xlo || paramY>yhi || paramY<ylo || paramZ<zlo || paramZ>zhi) return false;
        	}// theta loop
return true;
	}// ConeInTPC

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------------------
	// This returns Polygon edges for the cone If you already have the cone points
std::vector<larutil::PxPoint> geoconic::ConicalEdge(std::vector<larutil::PxPoint> incone){
	std::vector<larutil::PxPoint> edges; 
	int n = incone.size();
	int p = 0;
	int l ;
    for (int i =1 ; i < n; i++){
        if(incone[p].w > incone[i].w){p = i;}
	else if(incone[p].w == incone[i].w && incone[p].t > incone[i].t){p = i;}
	}// for loop over all the points
	l=p;
	int q = (p+1)%n;
	int count=0;
	std::vector<larutil::PxPoint> result(n+1);
		do
		{
			result[count]=incone[p];
			count++;
				for( int i=0; i<n; i++){
					    if(orientation(incone[p],incone[i],incone[q])==2){q=i;}
						}// for over i 
			p=q;
			q = (p+1)%n;
		}while(p!=l);
	for(unsigned int b = 0 ; b<result.size(); b++){
		if( result[b].w==0 && result[b].t==0) break;
		edges.push_back(result[b]);
		std::cout<<" The output values "<<result[b].w << " , " <<result[b].t<<std::endl;}
	
	edges = result;
	return edges;
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------------------
	// This returns Polygon edges for the cone
        // List of Conical Features that make the polygon 
        // Here we will need to give the input 
                // 1 Position of vertex 
                // 2 Direction of vertex 
                // 3 Length of cone 
                // 4 Opening Angle 
		// 5 Plane
		// 6 smoothness 
        std::vector<larutil::PxPoint> geoconic::ConicalFeatures(const TLorentzVector& Pos, const TLorentzVector& dir, double Length, double OpeningAngle, int plane, int smoothness){
	std::vector<larutil::PxPoint> ret;
	// Calculate the radius at a given length and opening angle 
 	double Radius = Length*tan(OpeningAngle/2 *PI/180);
        auto geom = larutil::GeometryUtilities::GetME();

		TLorentzVector pos = Pos;
                TVector3 Axis;
                double axisnorm = sqrt(dir.Px()*dir.Px()+dir.Py()*dir.Py()+dir.Pz()*dir.Pz());
                Axis.SetX(dir.Px()/axisnorm);
                Axis.SetY(dir.Py()/axisnorm);
                Axis.SetZ(dir.Pz()/axisnorm);
                // Get orthogonal vector
                TVector3 avec;
                avec = Axis.Orthogonal();
                // cross these vects to get the other mut-orthoganal one
                TVector3 bvec;
                bvec = Axis.Cross(avec);
		// Fill up the vertex 2d point
                 std::vector<larutil::PxPoint> pConeFullHits;
		// Fill up the Final 2d point
                 std::vector<larutil::PxPoint> CFHA;
		// Fill up the Final 2d point
                 std::vector<larutil::PxPoint> CFH;
                auto vert = geom->Get2DPointProjectionCM(&pos,plane);
                pConeFullHits.push_back(vert);
		std::cout<<"Vertex Position 3D "<< pos.X()<<" , " <<pos.Y()<<" , "<<pos.Z()<<std::endl;
		std::cout<<"Vertex Dir 3D "<< dir.Px()<<" , " <<dir.Py()<<" , "<<dir.Pz()<<std::endl;
		std::cout<<"Vertex Position 2D "<< vert.w<<" , "<<vert.t<<std::endl;
//		std::cout<<"Plane : "<<plane<<std::endl;

  for(double theta =0 ; theta<2*PI; theta+=PI/smoothness){
                double paramx = Pos.X() + Length*Axis.X()/Axis.Mag() + Radius*cos(theta)*avec.X() + Radius*sin(theta)*bvec.X();
                double paramy = Pos.Y() + Length*Axis.Y()/Axis.Mag() + Radius*cos(theta)*avec.Y() + Radius*sin(theta)*bvec.Y();
                double paramz = Pos.Z() + Length*Axis.Z()/Axis.Mag() + Radius*cos(theta)*avec.Z() + Radius*sin(theta)*bvec.Z();
		TLorentzVector FillPoint;
                FillPoint.SetX(paramx);
                FillPoint.SetY(paramy);
                FillPoint.SetZ(paramz);
                auto pt =geom->Get2DPointProjectionCM(&FillPoint,plane);
                pConeFullHits.push_back(pt);
                }

	// Check if there are collinear points... If so deal with them
	//&&&&&&&&&&&&
		std::vector<std::vector<unsigned int>> cofix;
		for(unsigned int a=0; a<pConeFullHits.size(); a++){
			for(unsigned int b=a+1; b<pConeFullHits.size(); b++){
				for(unsigned int c=b+1; c<pConeFullHits.size(); c++){
				double collinear = (pConeFullHits[b].t - pConeFullHits[a].t)*(pConeFullHits[c].w- pConeFullHits[b].w) - (pConeFullHits[c].t-pConeFullHits[b].t)*(pConeFullHits[b].w - pConeFullHits[a].w);
			if(collinear==0){
				// log these points and deal with them later 
				//std::cout<<" Look we have collinear set of three points."<<std::endl;	
				std::vector<unsigned int> copt(3);
				copt[0] = a;
				copt[1] = b;
				copt[2] = c;
				cofix.push_back(copt);
					a = c;
				goto moveon;
						}// if collinear
				}
			}
		moveon:;
		}
		

			// Fill all the hits that are not collinear
		for(unsigned int a=0; a<pConeFullHits.size(); a++){
			bool goodpt = true;
			for(unsigned int co=0; co<cofix.size();co++){
				if(a==cofix[co][0] ||a==cofix[co][1] ||a==cofix[co][2]){ goodpt = false; break;} // don't fill this point
				}
			if(goodpt){CFHA.push_back(pConeFullHits[a]);}
			}

			// Need to average the collinear points and put them in too 
			// this is fudging for now... we could write a ghram scann to put in here. 
			// this will fail when we have  5 collinear points... but fix later

			for(unsigned int co=0; co<cofix.size();co++){
			//double avgw= pConeFullHits[cofix[co][0]].w ;	
			//double avgy= (pConeFullHits[cofix[co][0]].t + pConeFullHits[cofix[co][1]].t +pConeFullHits[cofix[co][2]].t)/3.0;	
			// just pick the middle one for now
			CFHA.push_back(pConeFullHits[cofix[co][1]]);
			}
			
			// Check for duplicate points
			std::vector<unsigned int > exclude;
			for(unsigned int a=0; a< CFHA.size();a++){
				for(unsigned int b=a+1; b< CFHA.size();b++){
					if(CFHA[a].w==CFHA[b].w && CFHA[a].t==CFHA[b].t){
					// need to get rid of point b
					exclude.push_back(b);	
					}
				}
			}
			

	
			for(unsigned int a=0; a<CFHA.size() ;a++){
			bool goodpt = true;
				for(unsigned int b=0; b< exclude.size(); b++){
						if(a==exclude[b]){goodpt = false;}
				}
			if(goodpt){CFH.push_back(CFHA[a]);}
			}
			

				
			
			if(cofix.size()!=0){
				std::cout<<" Look we have collinear set of three points."<<std::endl;	
				if(cofix.size()>=1){ std::cout<< "WE HAVE GREATER COFIX "<<std::endl;
					for(unsigned int a = 0; a<cofix.size();a++){
					//std::cout<<pConeFullHits[cofix[a][2]].t<<std::endl;
std::cout<<"Collin points:"<<pConeFullHits[cofix[a][0]].w<<" , " <<pConeFullHits[cofix[a][0]].t<<" : "<<pConeFullHits[cofix[a][1]].w<<" , " <<pConeFullHits[cofix[a][1]].t<<" : "<<pConeFullHits[cofix[a][2]].w<<" , " <<pConeFullHits[cofix[a][2]].t<<std::endl;
//std::cout<<"Collin points:"<<pConeFullHits[cofix[a][0]].w<<" , " <<pConeFullHits[cofix[a][0]].t<<" : "<<pConeFullHits[cofix[a][1]].w<<" , " <<pConeFullHits[cofix[a][1]].t<<" : "<<pConeFullHits[cofix[a][2]].w<<" , " <<pConeFullHits[cofix[a][2]].t<<std::endl;
					}
					}
			}
	
	/////////
	//&&&&&&&&&&&&
	/////////

	// Now do the walk. 
	//int n = pConeFullHits.size();
	int n = CFH.size();
	int p = 0;
	int l ;
    for (int i =1 ; i < n; i++){
        //if(pConeFullHits[p].w > pConeFullHits[i].w){p = i;}
        if(CFH[p].w > CFH[i].w){p = i;}
	//else if(pConeFullHits[p].w == pConeFullHits[i].w && pConeFullHits[p].t > pConeFullHits[i].t){p = i;}
	else if(CFH[p].w == CFH[i].w && CFH[p].t > CFH[i].t){p = i;}
	}//
	l=p;
	int q = (p+1)%n;
	int count=0;
	std::vector<larutil::PxPoint> result(n+1);
		do
		{
			//result[count]=pConeFullHits[p];
			result[count]=CFH[p];
			count++;
				for( int i=0; i<n; i++){
	//			if(orientation(pConeFullHits[p],pConeFullHits[i],pConeFullHits[q])==0) std::cout<<"Collinear point"<<std::endl;
					    //if(orientation(pConeFullHits[p],pConeFullHits[i],pConeFullHits[q])==2){q=i;}
					    if(orientation(CFH[p],CFH[i],CFH[q])==2){q=i;}
						// Are there colinear points		?
						}// for over i 
			p=q;
			q = (p+1)%n;
				if(count> n+1){	
				//ret.push_back(pConeFullHits[0]);
				ret.push_back(CFH[0]);
				std::cout<< " Going to bail "<<ret.size()<<std::endl;
					//// looking at stuff before the bail 
			for(unsigned int a=0; a<CFH.size(); a++){
		std::cout<<"Circle edge w , t 2D "<< CFH[a].w<<" , "<<CFH[a].t<<std::endl;
				}//
		
		
	/*
	  for(double theta =0 ; theta<2*PI; theta+=PI/smoothness){
                double paramx = Pos.X() + Length*Axis.X()/Axis.Mag() + Radius*cos(theta)*avec.X() + Radius*sin(theta)*bvec.X();
                double paramy = Pos.Y() + Length*Axis.Y()/Axis.Mag() + Radius*cos(theta)*avec.Y() + Radius*sin(theta)*bvec.Y();
                double paramz = Pos.Z() + Length*Axis.Z()/Axis.Mag() + Radius*cos(theta)*avec.Z() + Radius*sin(theta)*bvec.Z();
		TLorentzVector FillPoint;
                FillPoint.SetX(paramx);
                FillPoint.SetY(paramy);
                FillPoint.SetZ(paramz);
		//std::cout<<" circle edge xyz : "<<paramx<<" , "<<paramy<<" , "<<paramz<<std::endl;
                auto pt =geom->Get2DPointProjectionCM(&FillPoint,plane);
		std::cout<<"Circle edge w , t 2D "<< pt.w<<" , "<<pt.t<<std::endl;
               // pConeFullHits.push_back(pt);
          //      CFH.push_back(pt);
                }
*/
				// can all be deleted at any point 

					
				return ret;}
		}while(p!=l);
	for(unsigned int b = 0 ; b<result.size(); b++){
		if( result[b].w==0 && result[b].t==0) break;
		ret.push_back(result[b]);
		}
	return ret;
		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------------------------------------
	// This returns a vector of hits that are contained in the polygon
       std::vector<larutil::PxHit> geoconic::PolyContain(std::vector<larutil::PxHit> hits, std::vector<larutil::PxPoint> polygon ){
	double InfAdd = 10000000;
	std::vector<larutil::PxHit> rethits;
	// Loop over all hits
		for( auto hit: hits){
	// Need to make an extreme point for each hit 
	larutil::PxPoint extreme;
	extreme.w = hit.w+InfAdd;// are these in cm 
	extreme.t = hit.t;// are these in cm 
	larutil::PxPoint currenthit;
	currenthit.t = hit.t ; // are these in cm 
	currenthit.w = hit.w;// are these in cm 
	// Count intersections of the above line with sides of polygon
	int count = 0, i = 0;
	int n = polygon.size();
	    do
		{
		int next = (i+1)%n;
			if (doIntersect(polygon[i], polygon[next], currenthit, extreme))
			{
				count++;
				//Fill the hit
				rethits.push_back(hit);
			}
		i = next;
		} while (i != 0);
	}// for loop over all hits to check 
return rethits;
}// end 


///////////////====Needed for the cone conatain =========================================
int geoconic::orientation(larutil::PxPoint p,larutil::PxPoint q, larutil::PxPoint r){
    double val = (q.t - p.t) * (r.w - q.w) - (q.w - p.w) * (r.t - q.t);
    if (val == 0) return 0;  // colinear
    if (val > 0) return 1 ; // clock or counterclock wise
	else return 2;
}
//----------------------------------------------------------------------------------------
bool geoconic::onSegment(larutil::PxPoint p, larutil::PxPoint q, larutil::PxPoint r){
    if (q.w <= std::max(p.w, r.w) && q.w >= std::min(p.w, r.w) &&
            q.t <= std::max(p.t, r.t) && q.t >= std::min(p.t, r.t))
        return true;
    return false;
}
//----------------------------------------------------------------------------------------
bool geoconic::doIntersect(larutil::PxPoint p1, larutil::PxPoint q1, larutil::PxPoint p2, larutil::PxPoint q2){
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    // General case
    if (o1 != o2 && o3 != o4)
        return true;
    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;

    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;

     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;

    return false; // Doesn't fall in any of the above cases
}
///////////////===== End for the Contain ===========================================


























	}

#endif

