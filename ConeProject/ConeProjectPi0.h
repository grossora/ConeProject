/**
 * \file ConeProjectPi0.h
 *
 * \ingroup ConeProjectPi0
 * 
 * \brief Class def header for a class ConeProjectPi0
 *
 * @author ryan
 */

/** \addtogroup ConeProjectPi0

    @{*/

#ifndef LARLITE_CONEPROJECTPI0_H
#define LARLITE_CONEPROJECTPI0_H

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


#include "ClusterRecoUtil/ClusterParamsAlg.h"
#include "ClusterRecoUtil/CRUHelper.h"


namespace larlite {
  /**
     \class ConeProjectPi0
     User custom analysis class made by SHELL_USER_NAME
   */
  class ConeProjectPi0 : public ana_base{
  
  public:

    /// Default constructor
    ConeProjectPi0(){ _name="ConeProjectPi0"; _fout=0;}

    /// Default destructor
    virtual ~ConeProjectPi0(){}

    /** IMPLEMENT in ConeProjectPi0.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ConeProjectPi0.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ConeProjectPi0.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
        ::larlite::conicalprojection fconeproj;
        ::larlite::geoconic fgeoconic;
        ::larlite::ctools fctools;

	double AxisLength = 14*3;// Setting this here... This is the length of the cone
        double openingangle = 50.0; // magic number place holder for now. 
        int smoothness = 16;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone
	
        unsigned int RDraw = 50;// How many times to draw... just used for playying ... not necessary


	
	void InitializeAnaTree();

        TTree *FullTree;


	TH2D *cratio0;
	TH2D *cratio1;
	TH2D *cratio2;

	TH1D *dangle0;
	TH1D *dangle1;
	TH1D *dangle2;

	TH1D *dangledeg0;
	TH1D *dangledeg1;
	TH1D *dangledeg2;

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
