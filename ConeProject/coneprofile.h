/**
 * \file ConicalProjections.h
 *
 * \ingroup NCfilter
 * 
 * \brief Class def header for a class ForceRegions
 *
 * @author ryan
	-- Took these from EMShowerProfile 
	-- Not sure if these are ok for photon EM Showers
 */

///** \addtogroup Pi0Enhancer

#ifndef RECOTOOL_CONEPROFILE_H
#define RECOTOOL_CONEPROFILE_H

#include <iostream>
#include <vector>
#include "LArUtil/LArUtilManager.h"
#include "LArUtil/LArUtilBase.h"
#include "LArUtil/PxUtils.h"
#include "LArUtil/LArUtil-TypeDef.h"
#include <math.h>       /* Tan */
#define PI 3.14159265358979323846264338327950288419




namespace larlite {
  /**
     \class coneprofile 
     Describe me!
   */
  class coneprofile{

  public:

    /// Default constructor
    coneprofile(){};

    /// Default destructor
    virtual ~coneprofile(){};

  /// Given shower energy [MeV], returns the length of this shower
  double Length(const double energy) const;

  /// Given the distance [cm] from the shower start along the trunk, calculate energy containment
  double EnergyContainment(const double dist) const;

  /// Given shower energy find opening angle
  /// This is done with simple arctangent and calculation of shower length
  double OpeningAngle(const double energy) const;

  /// Given (nothing) return shower radius
  /// Shower radius doesn't depend on energy (see: Moliere radius on google)
  double ShowerRadius() const;


	};

}
#endif


