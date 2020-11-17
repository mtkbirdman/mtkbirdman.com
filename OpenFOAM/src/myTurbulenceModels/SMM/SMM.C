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

\*---------------------------------------------------------------------------*/

#include "SMM.H"
#include "wallDist.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SMM, 0);
addToRunTimeSelectionTable(LESModel, SMM, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> SMM::fSGS() const
{
    const volScalarField epsilon_(Ce_*k_*sqrt(k_)/delta() + 2.0*nu()*k_/sqr(y_)); //eq.(12)
    const volScalarField uEps_(pow(nu()*epsilon_,1.0/4.0)); //eq.(11)
    const volScalarField yEps_(((uEps_*y_)/nu())*sqrt((Cl_*y_)/delta())); //eq.(11)
    return scalar(1)-exp(-pow(yEps_/A0_,4.0/3.0)); //(eq.10)
}

void SMM::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}

void SMM::correctNonlinearStress(const volTensorField& gradU)
{
    nut_ = CSGS_*fSGS()*sqrt(k_)*delta(); //eq.(10)
    nut_.correctBoundaryConditions();
    
    dimensionedScalar smallValue1("smallValue1",dimensionSet(0,0,-2,0,0,0,0),SMALL);
    dimensionedScalar smallValue2("smallValue2",dimensionSet(0,2,-2,0,0,0,0),SMALL);
    
    volVectorField UHat_(filter_(U_));
    volSymmTensorField TauPrime(CB_*symm((U_ - UHat_)*(U_ - UHat_))); //eq.(4)
    volSymmTensorField TauPrimeA(dev(TauPrime));
    
    volSymmTensorField S(symm(gradU)); //eq.(2)
    volScalarField nuPrime(-(TauPrimeA && S)/(2.0*max(magSqr(S),smallValue1))); //eq.(5)

    nonlinearStress_ = 2.0*k_*(TauPrimeA - (-2.0*nuPrime*S))/max(tr(TauPrime),smallValue2); //eq.(9)
}

tmp<fvScalarMatrix> SMM::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*k_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SMM::SMM
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    nonlinearEddyViscosity<incompressible::LESModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_(dimensioned<scalar>::lookupOrAddToDict("Ck",coeffDict_,0.1)),
    Ce_(dimensioned<scalar>::lookupOrAddToDict("Ce",coeffDict_,0.835)),
    CB_(dimensioned<scalar>::lookupOrAddToDict("CB",coeffDict_,1.0)),
    CSGS_(dimensioned<scalar>::lookupOrAddToDict("CSGS",coeffDict_,0.05)),
    A0_(dimensioned<scalar>::lookupOrAddToDict("A0",coeffDict_,30.0)),
    Cl_(dimensioned<scalar>::lookupOrAddToDict("Cl",coeffDict_,4.0)),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    LESDelta_
    (
        IOobject
        (
            "LESDelta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        delta()
    ),
    
    y_(wallDist::New(mesh_).y()),
    filterPtr_(LESfilter::New(mesh_, coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, kMin_);

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool SMM::read()
{
    if (nonlinearEddyViscosity<incompressible::LESModel>::read())
    {
        Ck_.readIfPresent(coeffDict());
        Ce_.readIfPresent(coeffDict());
        CB_.readIfPresent(coeffDict());
        CSGS_.readIfPresent(coeffDict());
        A0_.readIfPresent(coeffDict());
        Cl_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

tmp<volScalarField> SMM::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Ce_*k_*sqrt(k_)/delta() + 2.0*nu()*k_/sqr(y_) //eq.(12)
        )
    );
}


void SMM::correct()
{
    if (!turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<incompressible::LESModel>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();
    volScalarField divU(fvc::div(fvc::absolute(phi_, U_)));
    
    volScalarField G(GName(), (nut_*dev(twoSymm(gradU)) - nonlinearStress_) && gradU);
    
    tmp<fvScalarMatrix> kEqn //eq.(13)
    (
        fvm::ddt(k_)                                        //Dk/dt
      + fvm::div(phi_, k_)                                  //Dk/dt
      - fvm::laplacian(nu() + Ck_*fSGS()*sqrt(k_)*delta(), k_) //d/dxj(nu+Ck*fsgs*sqrt(k)*delta)dk/dxj)
     ==
        G                                                   //tauij*dUi/dxj
      - fvm::SuSp((2.0/3.0)*divU, k_)                       //tauij*dUi/dxj
      - fvm::Sp(Ce_*sqrt(k_)/delta() + 2.0*nu()/sqr(y_), k_)  //epsilon
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, kMin_);

    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // Edn namespace incompressible
} // End namespace Foam

// ************************************************************************* //
