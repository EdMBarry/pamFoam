/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    pamFoam

Description
    This is the generic solver for the photoanaerobic model extended to
    CFD. This model includes the resolution of the radiative transfer equation,
    as well as the process biokinetics and multiphase flows. Functionality
    exists to neglect any of the aforementioned modules. An additional biofilm
    model is currently in development and will be included as a modular add-on.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiphaseSystem.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "IOMRFZoneList.H"
#include "CorrectPhi.H"
//#include "oxygenTransferModel.H"
#include "photoBioModel.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "readBiokineticsProperties.H"
    #include "createPamFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    scalar slamDampCoeff
    (
        fluid.lookupOrDefault<scalar>("slamDampCoeff", 1)
    );

    dimensionedScalar maxSlamVelocity
    (
        "maxSlamVelocity",
        dimVelocity,
        fluid.lookupOrDefault<scalar>("maxSlamVelocity", GREAT)
    );

    #include "fluidPhases.H"
    //#include "createOxygenTransferModel.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- ASM rates
        #include "pamRates.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            turbulence->correct();
            fluid.solve();
            rho = fluid.rho();
            #include "zonePhaseVolumes.H"

            //#include "TEqns.H"
            #include "UEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            #include "DDtU.H"

            // --- Limit liquid phase velocities in gas phase
            liquidPhase.U() *= 1 - pos(alphaGas - 0.99);

            // --- Calculate dissipation coefficient for pure gas regions
            dissipationCoeff = pos(alphaGas - 0.99)/runTime.deltaT();

            // --- Update scalar diffusivities
            DSAC = (turbulence->nut()/ScT + DSACValue)*(1 - pos(alphaGas - 0.5));
            DSS  = (turbulence->nut()/ScT + DSSValue) *(1 - pos(alphaGas - 0.5));
            DXS  = (turbulence->nut()/ScT + DXSValue) *(1 - pos(alphaGas - 0.5));
            DXPB = (turbulence->nut()/ScT + DXPBValue)*(1 - pos(alphaGas - 0.5));
            DXI = (turbulence->nut()/ScT + DXIValue)*(1 - pos(alphaGas - 0.5));
            DSIN  = (turbulence->nut()/ScT + DSINValue) *(1 - pos(alphaGas - 0.5));
            DSIP = (turbulence->nut()/ScT + DSIPValue)*(1 - pos(alphaGas - 0.5));
            DSI  = (turbulence->nut()/ScT + DSIValue) *(1 - pos(alphaGas - 0.5));

            // --- Solve ASM scalar equations
            #include "SACEqn.H"
            #include "SSEqn.H"
            #include "XSEqn.H"
            #include "XPBEqn.H"
            #include "XIEqn.H"
            #include "SINEqn.H"
            #include "SIPEqn.H"
            #include "SIEqn.H"

            #include "solveRadiativeField.H"
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //