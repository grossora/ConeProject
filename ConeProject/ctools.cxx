#ifndef RECOTOOL_CTOOLS_CXX
#define RECOTOOL_CTOOLS_CXX


#include "ctools.h"

namespace larlite {

	

//-----------------------------------------------------------------------------------------------------------
	// This is a charge/distance fit?
std::pair<double,double>  ctools::ChargeDistFit(std::vector<larutil::PxHit> hits , std::pair<double,double> coneaxis, double surfaceslope){
////////try fitting
        double na = 0;
        double wiretime = 0;
        double wire = 0;
        double time = 0;
        double wirewire = 0;

        for(auto const p : hits){
                double tildaw = (coneaxis.second - (p.t) + surfaceslope * p.w)/(surfaceslope - coneaxis.first);
                double tildat = coneaxis.first*tildaw + coneaxis.second;
                double distance= sqrt(pow(p.w - tildaw,2)+pow((p.t)- tildat,2) );
                        double weightloop = (unsigned int) p.charge/distance;
				// watch out for infinity
			if(weightloop> 1000000000) weightloop = 1000000000;
                        for(unsigned int a=0; a< weightloop; a++){
                            na +=1;
                            wiretime += p.w * p.t;
                            wire += p.w;
                            time += p.t;
                            wirewire += p.w * p.w;
                                        }//end of loop over the hits
        }//

       double slope = (na * wiretime - wire * time)/(na*wirewire -wire * wire);
       double cept = time/na - slope *(wire / na);
	std::pair<double,double> ret(slope,cept);
	return ret;
	}
	

}
#endif

