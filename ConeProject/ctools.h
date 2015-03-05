/**
 * \file ConicalProjections.h
 *
 * \ingroup NCfilter
 * 
 * \brief Class def header for a class ForceRegions
 *
 * @author ryan
 */

///** \addtogroup Pi0Enhancer

#ifndef RECOTOOL_CTOOLS_H
#define RECOTOOL_CTOOLS_H

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
     \class ctools 
     Describe me!
   */
  class ctools{

  public:

    /// Default constructor
    ctools(){};

    /// Default destructor
    virtual ~ctools(){};

        //$ Is a 3D point Contained in TPC
std::pair<double,double>  ChargeDistFit(std::vector<larutil::PxHit> hits , std::pair<double,double> coneaxis, double surfaceslope);

	};

}
#endif


