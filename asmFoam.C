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
    asmFoam

Description
    Solver for an activated sludge system based on ASM1 biokinetics.

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
#include "oxygenTransferModel.H"

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
    #include "createAsmFields.H"
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
    #include "createOxygenTransferModel.H"

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
        #include "asmRates.H"

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
            DSS  = (turbulence->nut()/ScT + DSSValue) *(1 - pos(alphaGas - 0.5));
            DXS  = (turbulence->nut()/ScT + DXSValue) *(1 - pos(alphaGas - 0.5));
            DXBH = (turbulence->nut()/ScT + DXBHValue)*(1 - pos(alphaGas - 0.5));
            DXBA = (turbulence->nut()/ScT + DXBAValue)*(1 - pos(alphaGas - 0.5));
            DXP  = (turbulence->nut()/ScT + DXPValue) *(1 - pos(alphaGas - 0.5));
            DSO  = (turbulence->nut()/ScT + DSOValue) *(1 - pos(alphaGas - 0.5));
            DSNO = (turbulence->nut()/ScT + DSNOValue)*(1 - pos(alphaGas - 0.5));
            DSNH = (turbulence->nut()/ScT + DSNHValue)*(1 - pos(alphaGas - 0.5));
            DSND = (turbulence->nut()/ScT + DSNDValue)*(1 - pos(alphaGas - 0.5));
            DXND = (turbulence->nut()/ScT + DXNDValue)*(1 - pos(alphaGas - 0.5));

            // --- Solve ASM scalar equations
            #include "SSEqn.H"
            #include "XSEqn.H"
            #include "XBHEqn.H"
            #include "XBAEqn.H"
            #include "XPEqn.H"
            #include "SOEqn.H"
            #include "SNOEqn.H"
            #include "SNHEqn.H"
            #include "SNDEqn.H"
            #include "XNDEqn.H"
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
