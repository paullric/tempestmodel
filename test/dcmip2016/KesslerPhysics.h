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
#include "DataVector.h"
#include "DataMatrix.h"
#include "GridData3D.h"
#include "GridData4D.h"


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
    DataVector<double> qv;
    DataVector<double> qc;
    DataVector<double> qr;
    DataVector<double> rho;
    DataVector<double> zc;
    DataVector<double> pk;
    DataVector<double> t;
    
    int nz;

};

///////////////////////////////////////////////////////////////////////////////



#endif /* defined(____KesslerPhysics__) */
