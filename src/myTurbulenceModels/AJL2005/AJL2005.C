/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "AJL2005.H"
#include "wallDist.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(AJL2005, 0);
addToRunTimeSelectionTable(RASModel, AJL2005, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> AJL2005::Rt() const
{
    return sqr(k_)/(nu()*epsilon_);
}

tmp<volScalarField> AJL2005::fw(const dimensioned<scalar> xi) const
{
    return exp(-sqr((pow(nu()*epsilon_,0.25)*y_/nu())/xi));
}

tmp<volScalarField> AJL2005::fMu(const volScalarField& Rt,const volScalarField& fw) const
{
    return (scalar(1) + (35.0/pow(Rt,0.75))*exp(-pow(Rt/30.0,0.75)))
            *(scalar(1) - fw);
}

tmp<volScalarField> AJL2005::f2(const volScalarField& Rt,const volScalarField& fw) const //feps
{
    return (scalar(1) - 0.3*exp(-sqr(Rt/6.5)))
            *(scalar(1) - fw);
}

void AJL2005::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}


void AJL2005::correctNonlinearStress(const volTensorField& gradU)
{
    nut_ = Cmu_*fMu(Rt(),fw(26.0))*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    volSymmTensorField S(symm(gradU));
    volTensorField W(-skew(gradU));

    volScalarField tau(nut_/k_);
    
    volScalarField fB(scalar(1) + Ceta_*(CD_*tau)*(mag(W) - mag(S)));
    volScalarField CB(1.0/(scalar(1) + (22.0/3.0)*sqr(CD_*tau)*magSqr(W) + (2.0/3.0)*sqr(CD_*tau)*(magSqr(W) - magSqr(S))*fB));
    
    volScalarField fr1((magSqr(W) - magSqr(S))/max(magSqr(W) + magSqr(S),dimensionedScalar("minWS",dimensionSet(0,0,-2,0,0,0,0),SMALL)));
    volScalarField fr2(magSqr(S)/max(magSqr(W) + magSqr(S),dimensionedScalar("minWS",dimensionSet(0,0,-2,0,0,0,0),SMALL)));
    
    volScalarField fs1(fr1*fr2*(0.15*Ceta_)*sqr(CD_*tau)*(magSqr(W) - magSqr(S)));
    volScalarField fs2(-fr1*fr2*(scalar(1)+(0.07*Ceta_)*(CD_*tau)*(mag(W) - mag(S))));
    
    volScalarField nStar((nu()*epsilon_,0.25)*y_/nu());
    volVectorField N(fvc::grad(nStar));
    volVectorField d(N/mag(N));
    volSymmTensorField dd(symm(d*d));

    volScalarField taud1((scalar(1) - fw(15.0))*(k_/epsilon_) + fw(15.0)*(1.0)*sqrt(nu()/epsilon_));
    volSymmTensorField bw1(
        fw(26.0)*(
            - (1.0)*(1.0/2.0)*dev(dd)
            + (scalar(1) - sqr(fr1))*sqr(taud1)*(
                - ((0.25*0.5)/(scalar(1) + 0.5*sqr(taud1)*mag(S)*mag(W)))*twoSymm(S&W)
                + ((1.5*0.5)/(scalar(1) + 0.5*sqr(taud1)*magSqr(S)))*dev(innerSqr(S))
            )
        )
    );
    
    volScalarField taud2((scalar(1)-fw(15.0))*(k_/epsilon_) + fw(15.0)*(3.0)*sqrt(nu()/epsilon_));
    volSymmTensorField bw2(
        fw(26.0)*(
            - (0.0)*(1.0/2.0)*dev(dd)
            + (scalar(1) - sqr(fr1))*sqr(taud2)*(
                - (((13.0/30.0)*1.0)/(scalar(1) + 1.0*sqr(taud2)*mag(S)*mag(W)))*twoSymm(S&W)
                + ((0.6*1.0)/(scalar(1) + 1.0*sqr(taud2)*magSqr(S)))*dev(innerSqr(S))
            )
        )
    );

    volScalarField ftau(exp(-3.0*(scalar(1) - fw(26.0))*(k_/epsilon_)*mag(S)));
    volSymmTensorField bw(ftau*bw1 + (scalar(1) - ftau)*bw2);
    
    nonlinearStress_ = 
        2.0*nut_*S // Cancel linear term
        +2.0*(k_/CD_)*(
            - CB*(tau*CD_)*S //b1
            + (scalar(1) - fw(26.0))*(
                2.0*CB*sqr(tau*CD_)*(-twoSymm(S&W) + dev(innerSqr(S))) //b2
                - CB*fs1*(tau*CD_)*S + 2.0*fs2*sqr(tau*CD_)*dev(innerSqr(S)) //bs
            )
            + CD_*bw //bw
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AJL2005::AJL2005
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
    nonlinearEddyViscosity<incompressible::RASModel>
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

    Ceps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps1",
            coeffDict_,
            1.45
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.83
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.2
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.5
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.12
        )
    ),
    CD_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CD",
            coeffDict_,
            0.8
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            coeffDict_,
            100.0
        )
    ),
    
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

    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    y_(wallDist::New(mesh_).y())
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool AJL2005::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        CD_.readIfPresent(coeffDict());
        Ceta_.readIfPresent(coeffDict());

        return true;
    }

    return false;
}


void AJL2005::correct()
{
    if (!turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        Ceps1_*G*epsilon_/k_
      - fvm::Sp(Ceps2_*f2(Rt(),fw(3.3))*epsilon_/k_, epsilon_)
    );

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
