/**
 * \file GammaContainment.h
 *
 * \ingroup GammaContainment
 * 
 * \brief Class def header for a class GammaContainment
 *
 * @author ryan
 */

/** \addtogroup GammaContainment

    @{*/

#ifndef LARLITE_GAMMACONTAINMENT_H
#define LARLITE_GAMMACONTAINMENT_H

#include "Analysis/ana_base.h"
#include "LArUtil/LArUtilManager.h"
#include "LArUtil/PxUtils.h"
#include "LArUtil/LArUtilBase.h"
#include "LArUtil/LArUtil-TypeDef.h"
#include "LArUtil/TimeService.h"

#include "DataFormat/DataFormat-TypeDef.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/data_base.h"

#include <math.h>       /* aTan */
#define PI 3.14159265358979323846264338327950288419
#include <stdlib.h>     /* srand, rand */

#include "TH1.h"
#include "TH2.h"


#include "geoconic.h"
#include "ctools.h"
#include "conicalprojection.h"
#include "coneprofile.h"


#include "ClusterRecoUtil/ClusterParamsAlg.h"
#include "ClusterRecoUtil/CRUHelper.h"
#include "EMShowerTools/EMShowerProfile.h"

namespace larlite {
  /**
     \class GammaContainment
     User custom analysis class made by SHELL_USER_NAME
   */
  class GammaContainment : public ana_base{
  
  public:

    /// Default constructor
    GammaContainment(){ _name="GammaContainment"; _fout=0;}

    /// Default destructor
    virtual ~GammaContainment(){}

    /** IMPLEMENT in GammaContainment.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GammaContainment.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GammaContainment.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:


        ::larlite::conicalprojection fconeproj;
        ::larlite::geoconic fgeoconic;
        ::larlite::ctools fctools;
        ::larlite::coneprofile fconeprofile;

//	double AxisLength = 200;// Setting this here... This is the length of the cone
        //double openingangle = 50.0; // magic number place holder for now. 
        double openingangle ; // magic number place holder for now. 
        int smoothness = 16;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone
	

	TH1D *fLength;
	TH2D *fLengthA;

	TH1D *C0L50;
	TH1D *C1L50;
	TH1D *C2L50;

	TH1D *C0L70;
	TH1D *C1L70;
	TH1D *C2L70;

	TH1D *C0L100;
	TH1D *C1L100;
	TH1D *C2L100;

	TH1D *C0L140;
	TH1D *C1L140;
	TH1D *C2L140;

	TH1D *CTheory0;
	TH1D *CTheory1;
	TH1D *CTheory2;
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
