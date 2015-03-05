/**
 * \file Pi0Contained.h
 *
 * \ingroup Pi0Contained
 * 
 * \brief Class def header for a class Pi0Contained
 *
 * @author ryan
 */

/** \addtogroup Pi0Contained

    @{*/

#ifndef LARLITE_PI0CONTAINED_H
#define LARLITE_PI0CONTAINED_H

#include "Analysis/ana_base.h"
#include "LArUtil/LArUtilManager.h"
#include "LArUtil/PxUtils.h"
#include "LArUtil/LArUtilBase.h"
#include "LArUtil/LArUtil-TypeDef.h"
#include "LArUtil/TimeService.h"

#include "DataFormat/DataFormat-TypeDef.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/data_base.h"

#include <math.h>       /* aTan */
#define PI 3.14159265358979323846264338327950288419
#include <stdlib.h>     /* srand, rand */

#include "TH1.h"
#include "TH2.h"


#include "geoconic.h"
#include "conicalprojection.h"
#include "coneprofile.h"


namespace larlite {
  /**
     \class Pi0Contained
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi0Contained : public ana_base{
  
  public:

    /// Default constructor
    Pi0Contained(){ _name="Pi0Contained"; _fout=0;}

    /// Default destructor
    virtual ~Pi0Contained(){}

    /** IMPLEMENT in Pi0Contained.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in Pi0Contained.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in Pi0Contained.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
        ::larlite::conicalprojection fconeproj;
        ::larlite::geoconic fgeoconic;
        ::larlite::coneprofile fconeprofile;


        int smoothness = 16;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone


	
	void InitializeAnaTree();

        TTree *ConeTree;
	double Energy_A;
	double Energy_B;
	double AxisLength_A;
	double AxisLength_B;
	double OpeningAngle_A;
	double OpeningAngle_B;
	bool coneintpc_A;
	bool coneintpc_B;
	bool overlap;
//////////////////




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
