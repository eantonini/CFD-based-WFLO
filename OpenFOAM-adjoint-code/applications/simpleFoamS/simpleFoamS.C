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
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
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
    #include "initContinuityErrs.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity SIMPLE corrector
        {
	    // Momentum predictor
	    fvVectorMatrix UEqn
	    (
	        fvm::div(phi, U)
	      + turbulence->divDevReff(U)
	      ==
	        fvOptions(U)
	    );
	    UEqn.relax();
	    fvOptions.constrain(UEqn);
	    solve(UEqn == -fvc::grad(p));
	    fvOptions.correct(U);

	    volScalarField rAU(1.0/UEqn.A());
	    volVectorField HbyA("HbyA", U);
	    HbyA = rAU*UEqn.H();
	    surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(HbyA) & mesh.Sf());
	    adjustPhi(phiHbyA, U, p);
	    // Non-orthogonal pressure corrector loop
	    while (simple.correctNonOrthogonal())
	    {
	        fvScalarMatrix pEqn
	        (
	            fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
	        );
	        pEqn.setReference(pRefCell, pRefValue);
	        pEqn.solve();
	        if (simple.finalNonOrthogonalIter())
	        {
	            phi = phiHbyA - pEqn.flux();
	        }
	    }
	    #include "continuityErrs.H"
	    // Explicitly relax pressure for momentum corrector
	    p.relax();
	    // Momentum corrector
	    U = HbyA - rAU*fvc::grad(p);
	    U.correctBoundaryConditions();
	    fvOptions.correct(U);
        }
	laminarTransport.correct();
        turbulence->correct();
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
