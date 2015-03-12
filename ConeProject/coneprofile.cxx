#ifndef RECOTOOL_CONEPROFILE_CXX
#define RECOTOOL_CONEPROFILE_CXX


#include "coneprofile.h"

namespace larlite {

	

//-----------------------------------------------------------------------------------------------------------

double coneprofile::Length(const double energy) const {

  /// This formula taken from Andrzej
  //double rad_length_cm = 48.; //Assumed energy independent (roughly correct)
  double rad_length_cm = 14.; //Assumed energy independent (roughly correct)
  //this should be 14cm, but studying mcshowers tells me g4 is using 
  //47cm (fit)
//  double epsilon_mev = 5715.; //Parameter fit from geant4 output apparently
  double epsilon_mev = 30.5; // from using argon
  //double shower_length = rad_length_cm * ( log( energy/epsilon_mev ) - 1 + (0.08*18) + 9.6 );
  double shower_length =  log( energy/epsilon_mev ) - 0.5 + (0.08*18) + 9.6*rad_length_cm;
  return shower_length;
}


double coneprofile::EnergyContainment(const double dist) const {

  // This formula taken from my understanding of electromagnetic showers
  // They deposit (1- 1/e) of their energy in each radiation length, X
  // So E(length of shower, d) = E_o * (1-e^(-d/X))
  // This is assuming radiation length is energy independent.
  // Ex: plug in d = 0: No energy deposited
  // Ex: plug in d = inf: All of initial energy E_o is deposited
  // Ex: plug in d = X: All but 1/e of energy E_o is deposited

  double rad_length_cm = 48.; //Assumed energy independent (roughly correct)
  //this should be 14cm, but studying mcshowers tells me g4 is using 
  //47cm (fit)

  // If you're outside of the TPC, no energy is contained
  if(dist < 0) return 0;

  return (1 - pow(2.71828,(dist*(-1)/rad_length_cm)));
}


double coneprofile::OpeningAngle(const double energy) const {

  //          /                ^ 
  //         /                 |R
  //        /)(1/2)*o_angle    |
  //      <- - - - - - - - - - - - - - - 
  //        \    <-- length -->
  //         \
  //          \

  //(1/2)*opening angle = arctan(R/L)
  return 2*atan(ShowerRadius()/GammaLength(energy));

}

double coneprofile::ShowerRadius() const {
  //Moleiere radius in LAr is 10.1cm, according to Wikipedia
  return 10.1;
}



//===========================================================================
//=================================== RG=====================================
//===========================================================================

double coneprofile::GammaLength(const double energy) const {
  double rad_length_cm = 14.; //Assumed energy independent (roughly correct)
  double epsilon_mev = 30.5; // from using argon

  double shower_length =  log( energy/epsilon_mev )*rad_length_cm - 0.5 + (0.08*18) + 9.6*rad_length_cm;
	// Since we do not know where the Photon will convert... we can add an extra 9/7*radLength
  double gshower_length = shower_length+ 9.0/7.0 *rad_length_cm;
std::cout<<"Shower of Energy: "<< energy <<" Has Length(cm of :"<<gshower_length<<std::endl;
	return gshower_length;
	}
//===========================================================================

//double coneprofile::Gammatmax(const double energy) const {
  //double rad_length_cm = 14.; //Assumed energy independent (roughly correct)
  //double epsilon_mev = 30.5; // from using argon
  //double shower_length =  log( energy/epsilon_mev ) - 0.5 + ;

}
#endif

