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

Class
    Foam::incompressible::LESModels::SMM

Group
    grpLESTurbulence

Description
    

SourceFiles
    SMM.C

\*---------------------------------------------------------------------------*/

#ifndef SMM_H
#define SMM_H

#include "turbulentTransportModel.H"
#include "nonlinearEddyViscosity.H"
#include "LESfilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class SMM Declaration
\*---------------------------------------------------------------------------*/

class SMM
:
    public nonlinearEddyViscosity<incompressible::LESModel>
{
    // Private Member Functions

protected:

    // Protected data

        // Model constants

            dimensionedScalar Ck_;
            dimensionedScalar Ce_;
            dimensionedScalar CB_;
            dimensionedScalar CSGS_;
            dimensionedScalar A0_;
            dimensionedScalar Cl_;

        // Fields

            volScalarField k_;
            volScalarField LESDelta_;

            //- Wall distance
            const volScalarField& y_;
        
        // Filters

            //simpleFilter simpleFilter_;
            autoPtr<LESfilter> filterPtr_;
            LESfilter& filter_;



    // Protected Member Functions

        tmp<volScalarField> fSGS() const;
        virtual void correctNut();
        virtual void correctNonlinearStress(const volTensorField& gradU);
        virtual tmp<fvScalarMatrix> kSource() const;


public:

    //- Runtime type information
    TypeName("SMM");


    // Constructors

        //- Constructor from components
        SMM
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~SMM()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the (estimated) specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            volScalarField epsilon_(Ce_*k_*sqrt(k_)/delta());
            return tmp<volScalarField>::New
            (
                IOobject
                (
                    IOobject::groupName("omega", this->alphaRhoPhi_.group()),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                epsilon_/(0.09*k_)
            );
        }

        //- Correct eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // Edn namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
