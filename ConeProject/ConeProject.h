/**
 * \file ConeProject.h
 *
 * \ingroup ConeProject
 * 
 * \brief Class def header for a class ConeProject
 *
 * @author ryan
 */

/** \addtogroup ConeProject

    @{*/

#ifndef LARLITE_CONEPROJECT_H
#define LARLITE_CONEPROJECT_H

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
     \class ConeProject
     User custom analysis class made by SHELL_USER_NAME
   */
  class ConeProject : public ana_base{
  
  public:

    /// Default constructor
    ConeProject(){ _name="ConeProject"; _fout=0;}

    /// Default destructor
    virtual ~ConeProject(){}

    /** IMPLEMENT in ConeProject.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ConeProject.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ConeProject.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
        ::larlite::conicalprojection fconeproj;
        ::larlite::geoconic fgeoconic;
        ::larlite::ctools fctools;
        ::larlite::coneprofile fconeprofile;

	double AxisLength = 200;// Setting this here... This is the length of the cone
        //double openingangle = 50.0; // magic number place holder for now. 
        double openingangle ; // magic number place holder for now. 
        //int smoothness = 16;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone
        int smoothness = 64;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone
	
        unsigned int RDraw = 5;// How many times to draw... just used for playying ... not necessary
        //unsigned int RDraw = 500;// How many times to draw... just used for playying ... not necessary


	
	void InitializeAnaTree();

        TTree *FullTree;
	double _Energy;
	double _Vertexx;
	double _Vertexy;
	double _Vertexz;

        TTree *CRTree;
	double _RandomLength;
	double _CRatio0;
	double _CRatio1;
	double _CRatio2;
        TTree *DATree;
	double _DangleDeg0;
	double _DangleDeg1;
	double _DangleDeg2;

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
