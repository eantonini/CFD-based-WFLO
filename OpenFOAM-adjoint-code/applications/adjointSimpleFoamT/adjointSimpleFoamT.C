/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createControl.H"
    #include "createFvOptions.H"
    #include "initAContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
	    // scalar Cmu = 0.09; scalar beta = 0.072; scalar alpha = 0.52; scalar alphaK = 0.5; scalar alphaOmega = 0.5;
	    scalar Cmu = 0.0333; scalar beta = 0.0277; scalar alpha = 0.42; scalar alphaK = 0.45; scalar alphaOmega = 0.45;

	    volScalarField nuEff = nu + nut;
	    volScalarField DomegaEff = nu + alphaOmega*nut;
	    volScalarField DkEff = nu + alphaK*nut;

	    volTensorField gradU = fvc::grad(U);
	    volTensorField gradU_a = fvc::grad(U_a);
	    volVectorField gradK = fvc::grad(k);
	    volVectorField gradK_a = fvc::grad(k_a);
	    volVectorField gradOmega = fvc::grad(omega);
	    volVectorField gradOmega_a = fvc::grad(omega_a);
	    volSymmTensorField S = symm(gradU);
	    volSymmTensorField S_a = symm(gradU_a);
	    volScalarField SS = magSqr(S);

	    // Momentum equation
            fvVectorMatrix UEqn_a
    	    (
	        fvm::div(-phi, U_a)
	      - fvm::laplacian(nuEff, U_a)
	      - fvc::div(nuEff*dev2(T(fvc::grad(U_a))))
	      ==
	        fvOptions(U_a)
	      - (gradU & U_a)
	      - gradK*k_a
	      - gradOmega*omega_a
	      - 2.0*(fvc::laplacian(k_a*k/omega, U) + fvc::div(k_a*k/omega*dev2(T(gradU))))
	      - 2.0*alpha*(fvc::laplacian(omega_a, U) + fvc::div(omega_a*dev2(T(gradU))))
     	    );
    	    UEqn_a.relax();
    	    solve(UEqn_a == -fvc::grad(p_a));

	    volScalarField rAU_a(1.0/UEqn_a.A());
	    volVectorField HbyA_a("HbyA_a", U_a);
	    HbyA_a = rAU_a*UEqn_a.H();
	    surfaceScalarField phiHbyA_a("phiHbyA_a", fvc::interpolate(HbyA_a) & mesh.Sf());
	    adjustPhi(phiHbyA_a, U_a, p_a);
	    // Non-orthogonal pressure corrector loop
	    while (simple.correctNonOrthogonal())
	    {
		// Pressure equation
    		fvScalarMatrix pEqn_a
   		(
    		    fvm::laplacian(rAU_a, p_a) == fvc::div(phiHbyA_a)
    		);
		pEqn_a.setReference(p_aRefCell, p_aRefValue);
        	pEqn_a.solve();
		if (simple.finalNonOrthogonalIter())
        	{
            	    phi_a = phiHbyA_a - pEqn_a.flux();
        	}
	    }
	    #include "aContinuityErrs.H"
	    // Explicitly relax pressure for momentum corrector
	    p_a.relax();
	    // Momentum corrector
	    U_a = HbyA_a - rAU_a*fvc::grad(p_a);
	    U_a.correctBoundaryConditions();

	    // Update omega and G at the wall
	    // Turbulence specific dissipation rate equation
	    fvScalarMatrix omegaEqn_a
	    (
	        fvm::div(-phi, omega_a)
	      - fvm::laplacian(DomegaEff, omega_a)
	     ==
	        k/sqr(omega)*(2.0*S && gradU_a)
	      + alphaK*k/sqr(omega)*(gradK & gradK_a)
	      + alphaOmega*k/sqr(omega)*(gradOmega & gradOmega_a)
	      - k*k_a*(2.0*SS/sqr(omega)+Cmu)
	      - 2.0*beta*omega*omega_a
	    );
	    omegaEqn_a.relax();
	    solve(omegaEqn_a);

	    // Turbulent kinetic energy equation
	    fvScalarMatrix kEqn_a
	    (
	        fvm::div(-phi, k_a)
	      - fvm::laplacian(DkEff, k_a)
	     ==
	      - 1.0/omega*(2.0*S && gradU_a)
	      - alphaK/omega*(gradK & gradK_a)
	      - alphaOmega/omega*(gradOmega & gradOmega_a)
	      + (2.0*SS/omega-Cmu*omega)*k_a
	    );
	    kEqn_a.relax();
	    solve(kEqn_a);
        }
	runTime.write();
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
