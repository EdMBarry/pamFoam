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

#include "oxygenTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oxygenTransferModel, 0);
    defineRunTimeSelectionTable(oxygenTransferModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oxygenTransferModel::oxygenTransferModel
(
    const dictionary& oxygenTransferModelDict,
    const phaseModel& liquidPhase,
    const phaseModel& gasPhase
)
:
    oxygenTransferModelDict_(oxygenTransferModelDict),
    liquidPhase_(liquidPhase),
    gasPhase_(gasPhase),
    kLa_
    (
        IOobject
        (
            "kLa",
            liquidPhase_.U().time().timeName(),
            liquidPhase_.U().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        liquidPhase_.U().mesh(),
        dimensionedScalar("kLa", dimless/dimTime, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::oxygenTransferModel::~oxygenTransferModel()
{}


// ************************************************************************* //
