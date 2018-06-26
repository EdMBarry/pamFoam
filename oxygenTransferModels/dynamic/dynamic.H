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

Class
    Foam::oxygenTransferModels::dynamic

Description
     Oxygen transfer model with dynamic kLa, with kL based on penetration 
     theory

SourceFiles
    dynamic.C

\*---------------------------------------------------------------------------*/

#ifndef oxygenTransferModel_dynamic_H
#define oxygenTransferModel_dynamic_H

#include "oxygenTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace oxygenTransferModels
{

/*---------------------------------------------------------------------------*\
      Class oxygenTransferModels::dynamic Declaration
\*---------------------------------------------------------------------------*/

class dynamic
:
    public oxygenTransferModel
{

protected:

    // Protected data
    dimensionedScalar DL_;

public:

    //- Runtime type information
    TypeName("dynamic");

    // Constructors

        //- Default constructor
        dynamic
        (
            const dictionary& oxygenTransferModelDict,
            const phaseModel& liquidPhase,
            const phaseModel& gasPhase
        );

    // Destructor
    virtual ~dynamic();


    // Member Functions

        //- Return kLa field
        tmp<volScalarField> kLa();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace oxygenTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

