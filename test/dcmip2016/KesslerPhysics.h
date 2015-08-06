//
//  KesslerPhysics.h
//  
//
//  Created by Antonin Verlet-Banide on 5/19/15.
//
//

#ifndef ____KesslerPhysics__
#define ____KesslerPhysics__

#include <iostream>


#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include "Defines.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "DataArray4D.h"


#include "WorkflowProcess.h"

///////////////////////////////////////////////////////////////////////////////

class GridPatch;
class Time;

  extern "C" {

      void tc_kessler_(
                       double * t,
                       double * qv,
                       double * qc,
                       double * qr,
                       double * rho,
                       double * pk,
                       double * z,
                       int * nz
                       );
  }


///	<summary>
///		KesslerPhysics type physics forcing.
///	</summary>
class KesslerPhysics : public WorkflowProcess {
    
public:
	///	<summary>
	///		Constructor.
	///	</summary>
	KesslerPhysics(
                      Model & model,
                      const Time & timeFrequency
                      );
    
public:
  	///	<summary>
	///		Initializer.
	///	</summary>
	virtual void Initialize(
        const Time & timeStart
    );
    
	///	<summary>
	///		Apply KesslerPhysics physics.
	///	</summary>
	virtual void Perform(
        const Time & time
    );
    
protected:
    DataArray1D<double> qv;
    DataArray1D<double> qc;
    DataArray1D<double> qr;
    DataArray1D<double> rho;
    DataArray1D<double> zc;
    DataArray1D<double> pk;
    DataArray1D<double> t;
    
    int nz;

};

///////////////////////////////////////////////////////////////////////////////



#endif /* defined(____KesslerPhysics__) */
