/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "epsilonLowReWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonLowReWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonLowReWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    //os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    //os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    //os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonLowReWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonLowReWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonLowReWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    epsilon_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


Foam::epsilonLowReWallFunctionFvPatchScalarField&
Foam::epsilonLowReWallFunctionFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const epsilonLowReWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonLowReWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonLowReWallFunctionFvPatchScalarField&>(epf);
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonLowReWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            //epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
            epf.calculate(turbulence, w, epf.patch(), epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonLowReWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));
    
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];

        const scalar w = cornerWeights[facei];

        epsilon0[celli] += w*2.0*k[celli]*nuw[facei]/sqr(y[facei]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    //Cmu_(0.09),
    //kappa_(0.41),
    //E_(9.8),
    //yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    //Cmu_(ptf.Cmu_),
    //kappa_(ptf.kappa_),
    //E_(ptf.E_),
    //yPlusLam_(ptf.yPlusLam_),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    //Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    //kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    //E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    //yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // Apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


Foam::epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    //Cmu_(ewfpsf.Cmu_),
    //kappa_(ewfpsf.kappa_),
    //E_(ewfpsf.E_),
    //yPlusLam_(ewfpsf.yPlusLam_),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonLowReWallFunctionFvPatchScalarField::
epsilonLowReWallFunctionFvPatchScalarField
(
    const epsilonLowReWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    //Cmu_(ewfpsf.Cmu_),
    //kappa_(ewfpsf.kappa_),
    //E_(ewfpsf.E_),
    //yPlusLam_(ewfpsf.yPlusLam_),
    epsilon_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonLowReWallFunctionFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, epsilon(true));
    }

    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;
    
    FieldType& epsilon = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        epsilon[celli] = epsilon0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}

void Foam::epsilonLowReWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintEpsilon(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& epsilon
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // Anly set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintEpsilon.append(epsilon[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }
    /*
    matrix.setValues
    (
        constraintCells,
        scalarField(constraintEpsilon.xfer())
    );
    */
    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonLowReWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonLowReWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
