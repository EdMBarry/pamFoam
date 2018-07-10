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
    Foam::oxygenTransferModel

Description
     Abstract base class for oxygen transfer models

SourceFiles
    oxygenTransferModel.C
    newOxygenTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef oxygenTransferModel_H
#define oxygenTransferModel_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class oxygenTransferModel Declaration
\*---------------------------------------------------------------------------*/

class oxygenTransferModel
{

protected:

    // Protected data


        //- oxygenTransferModel dictionary
        dictionary oxygenTransferModelDict_;

        //- Phase models
        const phaseModel& liquidPhase_;
        const phaseModel& gasPhase_;

        //- kLa field
        volScalarField kLa_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        oxygenTransferModel(const oxygenTransferModel&);

        //- Disallow default bitwise assignment
        void operator=(const oxygenTransferModel&);

public:

    //- Runtime type information
    TypeName("oxygenTransferModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        oxygenTransferModel,
        dictionary,
        (
            const dictionary& oxygenTransferModelDict,
            const phaseModel& liquidPhase,
            const phaseModel& gasPhase
        ),
        (oxygenTransferModelDict, liquidPhase, gasPhase)
    );

    // Constructors

        //- Default constructor
        oxygenTransferModel(
            const dictionary& oxygenTransferModelDict,
            const phaseModel& liquidPhase,
            const phaseModel& gasPhase
        );

    // Destructor
    virtual ~oxygenTransferModel();

    // Selectors

        //- Return a reference to the oxygen transfer model
        static autoPtr<oxygenTransferModel> New
        (
            const dictionary& oxygenTransferModelDict,
            const phaseModel& liquidPhase,
            const phaseModel& gasPhase
        );

    // Member Functions

        //- Return kLa field
        virtual tmp<volScalarField> kLa() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

