/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "dynamic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace oxygenTransferModels
{
    defineTypeNameAndDebug(dynamic, 0);
    addToRunTimeSelectionTable(oxygenTransferModel, dynamic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oxygenTransferModels::dynamic::dynamic
(
    const dictionary& oxygenTransferModelDict,
    const phaseModel& liquidPhase,
    const phaseModel& gasPhase
)
:
    oxygenTransferModel(oxygenTransferModelDict, liquidPhase, gasPhase),
    DL_(oxygenTransferModelDict.subDict("dynamicCoeffs").lookup("DL"))
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::oxygenTransferModels::dynamic::~dynamic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::oxygenTransferModels::dynamic::kLa()
{
    scalar pi = constant::mathematical::pi;
    volScalarField Vr = mag(liquidPhase_.U() - gasPhase_.U());
    volScalarField d = gasPhase_.d();
    volScalarField kL = 2*sqrt(DL_*Vr/pi/d);
    volScalarField alphaGas = gasPhase_;
    volScalarField a = 6*alphaGas/d/max((1 - alphaGas), SMALL);
    volScalarField factor = 1 - pos(alphaGas - 0.5);
    kLa_ = kL*a*factor;
    return kLa_;
}


// ************************************************************************* //
