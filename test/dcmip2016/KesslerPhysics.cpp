//
//  KesslerPhysics.cpp
//  
//
//  Created by Antonin Verlet-Banide on 5/19/15.
//
//

#include "KesslerPhysics.h"

#include "Model.h"

#include "Announce.h" 

///////////////////////////////////////////////////////////////////////////////

extern "C" {
	void kessler(
		double * t,
		double * qv,
		double * qc,
		double * qr,
		double * rho,
		double * pk,
		double * z,
		int * nz,
		double * rainnc);
}

///////////////////////////////////////////////////////////////////////////////

KesslerPhysics::KesslerPhysics(
	Model & model,
	const Time & timeFrequency
) :
WorkflowProcess(
	model,
	timeFrequency)
{ }

//////////////////////////////////////////////////////////////////////////////

	void KesslerPhysics::Initialize(
		const Time & timeStart
) {
		// Indices of EquationSet variables
		const int UIx = 0;
		const int VIx = 1;
		const int TIx = 2;
		const int WIx = 3;
		const int RIx = 4;


		// Get a copy of the GLL grid
		Grid * pGrid = m_model.GetGrid();

		// Check position of variables
		if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
			_EXCEPTIONT("Not implemented for --vstagger CPH, use --vstagger LOR");
		}

		int ndim=pGrid->GetRElements();

		qv.Allocate(ndim);
		qc.Allocate(ndim);
		qr.Allocate(ndim);
		t.Allocate(ndim);
		rho.Allocate(ndim);
		zc.Allocate(ndim);
		pk.Allocate(ndim);
	
}

	
///////////////////////////////////////////////////////////////////////////////
	
	void KesslerPhysics::Perform(
		const Time & time
) {
		// Indices of EquationSet variables
		const int UIx = 0;
		const int VIx = 1;
		const int TIx = 2;
		const int WIx = 3;
		const int RIx = 4;

		// Get DeltaT
		double dDeltaT = m_timeFrequency.GetSeconds();
		
		// Get a copy of the GLL grid
		Grid * pGrid = m_model.GetGrid();

		// Check position of variables
		if (pGrid->GetVarLocation(TIx) == DataLocation_REdge) {
			_EXCEPTIONT("Not implemented for --vstagger CPH, use --vstagger LOR");
		}

		// Physical constants
		const PhysicalConstants & phys = m_model.GetPhysicalConstants();
	
		// Number of radial elements
		int nz = pGrid->GetRElements();

		// Perform local update
		for (int n = 0; n < pGrid->GetActivePatchCount(); n++) {
		GridPatch * pPatch = pGrid->GetActivePatch(n);
		
		const PatchBox & box = pPatch->GetPatchBox();
		
		// Get latitude
		const DataArray2D<double> & dataLatitude = pPatch->GetLatitude();
		
		// Grid data
		DataArray4D<double> & dataNode =
		pPatch->GetDataState(0, DataLocation_Node);
		
		DataArray4D<double> & dataTracer =
		pPatch->GetDataTracers(0);
		
	   const DataArray3D<double> & dataZLevels	= pPatch->GetZLevels();
		
			// Loop over all horizontal nodes in GridPatch
			for (int i = box.GetAInteriorBegin(); i < box.GetAInteriorEnd(); i++) {
			for (int j = box.GetBInteriorBegin(); j < box.GetBInteriorEnd(); j++) {
		 
				// Loop over all levels in column
				for (int k = 0; k < pGrid->GetRElements(); k++) {
		  
					// Calculate pressure
					double dPressure = phys.PressureFromRhoTheta(dataNode[RIx][k][i][j] *dataNode[TIx][k][i][j]);
		 
					// Pointwise temperature
					double dT = dPressure / (dataNode[RIx][k][i][j] * phys.GetR());
		
					qv[k]=dataTracer[0][k][i][j]/dataNode[RIx][k][i][j];

					qc[k]=dataTracer[1][k][i][j]/dataNode[RIx][k][i][j];
	
					qr[k]=dataTracer[2][k][i][j]/dataNode[RIx][k][i][j];
		
					t[k]=dataNode[TIx][k][i][j];
					
					rho[k]=dataNode[RIx][k][i][j];
 
					zc[k]=dataZLevels[k][i][j];
	
					pk[k]=dataNode[TIx][k][i][j]/dT;
				}
		
				double rainnc = 0.0;
				kessler(t,qv,qc,qr,rho,pk,zc,&nz,&rainnc);
				for (int k = 0; k < pGrid->GetRElements(); k++) {
					dataNode[TIx][k][i][j]=t[k];
					dataNode[RIx][k][i][j]=rho[k];

					dataTracer[0][k][i][j]=qv[k]*dataNode[RIx][k][i][j];
					dataTracer[1][k][i][j]=qc[k]*dataNode[RIx][k][i][j];
					dataTracer[2][k][i][j]=qr[k]*dataNode[RIx][k][i][j];
				}
		
		
		
	
			}
			}
		}
	
		// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

