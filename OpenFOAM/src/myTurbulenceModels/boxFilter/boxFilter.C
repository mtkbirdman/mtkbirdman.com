/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "boxFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boxFilter, 0);
    addToRunTimeSelectionTable(LESfilter, boxFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boxFilter::boxFilter
(
    const fvMesh& mesh, 
    scalar widthCoeff
)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff)

{}


Foam::boxFilter::boxFilter(const fvMesh& mesh, const dictionary& dict)
:
    LESfilter(mesh),
    widthCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "widthCoeff",
            0.5
        )
    )

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boxFilter::read(const dictionary& dict)
{
    dict.optionalSubDict(type() + "Coeffs").readIfPresent<scalar>
    (
        "widthCoeff",
        widthCoeff_
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::boxFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    volScalarField neighbors(fvc::surfaceSum(mesh().magSf()/mesh().magSf()));

    tmp<volScalarField> filteredField = 
        (unFilteredField*(scalar(1) - widthCoeff_*neighbors) + fvc::surfaceSum(fvc::interpolate(unFilteredField)))/(scalar(1) + widthCoeff_*neighbors);

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::boxFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    volScalarField neighbors(fvc::surfaceSum(mesh().magSf()/mesh().magSf()));

    tmp<volVectorField> filteredField = 
        (unFilteredField*(scalar(1) - widthCoeff_*neighbors) + fvc::surfaceSum(fvc::interpolate(unFilteredField)))/(scalar(1) + widthCoeff_*neighbors);

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::boxFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    volScalarField neighbors(fvc::surfaceSum(mesh().magSf()/mesh().magSf()));

    tmp<volSymmTensorField> filteredField = 
        (unFilteredField*(scalar(1) - widthCoeff_*neighbors) + fvc::surfaceSum(fvc::interpolate(unFilteredField)))/(scalar(1) + widthCoeff_*neighbors);

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::boxFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    volScalarField neighbors(fvc::surfaceSum(mesh().magSf()/mesh().magSf()));

    tmp<volTensorField> filteredField = 
        (unFilteredField*(scalar(1) - widthCoeff_*neighbors) + fvc::surfaceSum(fvc::interpolate(unFilteredField)))/(scalar(1) + widthCoeff_*neighbors);

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //
