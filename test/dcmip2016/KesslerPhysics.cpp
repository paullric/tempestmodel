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
		double * dt,
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

				// Total density
				double dRho = dataNode[RIx][k][i][j];

				// Dry air density
				double dRhoD =
					dRho
					- dataTracer[0][k][i][j]
					- dataTracer[1][k][i][j]
					- dataTracer[2][k][i][j];

				// Virtual potential temperature
				double dThetav = dataNode[TIx][k][i][j];

				// Calculate moist pressure
				double dPressure =
					phys.PressureFromRhoTheta(dRho * dThetav);

				// Pointwise virtual temperature
				double dTv = dPressure / (dRho * phys.GetR());
	
				// Water vapor (RhoQv / Rho)
				qv[k] = dataTracer[0][k][i][j] / dataNode[RIx][k][i][j];

				// Cloud water (RhoQc / Rho)
				qc[k] = dataTracer[1][k][i][j] / dataNode[RIx][k][i][j];

				// Rain water (RhoQr / Rho)
				qr[k] = dataTracer[2][k][i][j] / dataNode[RIx][k][i][j];

				// Potential temperature
				t[k] = dThetav / (1.0 + 0.61 * qv[k]);

				// Dry air density
				rho[k] = dRhoD;
 
				// Heights of each level
				zc[k] = dataZLevels[k][i][j];
	
				// Exner function (Tv/theta)
				pk[k] = dTv / dThetav;
			}

			double rainnc = 0.0;
			kessler(
				&(t[0]),
				&(qv[0]),
				&(qc[0]),
				&(qr[0]),
				&(rho[0]),
				&(pk[0]),
				&(dDeltaT),
				&(zc[0]),
				&nz,
				&rainnc);

			// Update variables
			for (int k = 0; k < pGrid->GetRElements(); k++) {

				// Moist density
				dataNode[RIx][k][i][j]=
					rho[k] / (1.0 - qv[k] - qc[k] - qr[k]);

				// Virtual potential temperature
				dataNode[TIx][k][i][j] =
					t[k] * (1.0 + 0.61 * qv[k]);

				dataTracer[0][k][i][j] = qv[k] * dataNode[RIx][k][i][j];
				dataTracer[1][k][i][j] = qc[k] * dataNode[RIx][k][i][j];
				dataTracer[2][k][i][j] = qr[k] * dataNode[RIx][k][i][j];
			}
		}
		}
	}
	
	// Call up the stack to update performance time
	WorkflowProcess::Perform(time);
}

