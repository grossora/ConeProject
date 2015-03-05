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

#ifndef RECOTOOL_CONICALPROJECTION_H
#define RECOTOOL_CONICALPROJECTION_H

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
     \class ForceRegions
     Describe me!
   */
  class conicalprojection{

  public:

    /// Default constructor
    conicalprojection(){};

    /// Default destructor
    virtual ~conicalprojection(){};


//-----------------------------------------------------------------------------------------------------------
        // This is the pxpoint for the vertex?
	std::vector<larutil::PxPoint> vertexproj(TLorentzVector& pos);

//-----------------------------------------------------------------------------------------------------------
        // This is a pair of PxPoints for the start and end point
	std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> startendpt(TLorentzVector& pos,TLorentzVector& dir, double Length);

//-----------------------------------------------------------------------------------------------------------
        // This is a pair that contains the slope and cept for the axis in cm space 
	std::vector<std::pair<double,double>> ConeAxisSC(TLorentzVector& pos, TLorentzVector& dir, double length);







	};

}
#endif


