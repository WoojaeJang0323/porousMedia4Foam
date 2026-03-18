/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mixingPrecipitation.H"
#include "addToRunTimeSelectionTable.H"

#include "fvMatrix.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvm.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace geochemicalModels
    {
        defineTypeNameAndDebug(mixingPrecipitation, 0);

        addToRunTimeSelectionTable
        (
            basicGeochemicalModel,
            mixingPrecipitation,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geochemicalModels::mixingPrecipitation::mixingPrecipitation
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicGeochemicalModel(mesh, dict),
    mixingPrecipitationDict_(dict.subDict(typeName)),
    transportPropertiesDict_(dict),
    mineralSubDict_( mineralList_.size() ),
    Vm_( mineralList_.size() ),
    k_( mixingPrecipitationDict_.lookup("k") ),
    Ksp_( mixingPrecipitationDict_.lookup("Ksp") ),
    Di_( mixingPrecipitationDict_.lookup("Di") ),
    speciesAName_
    (
        mixingPrecipitationDict_.lookupOrDefault<word>("speciesA", "C_A")
    ),
    speciesBName_
    (
        mixingPrecipitationDict_.lookupOrDefault<word>("speciesB", "C_B")
    ),
    As_
    (
        IOobject
        (
            "As",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    )
{
    Info << "Initialization of mixingPrecipitation model ...." << nl;

    Info << "  As       : min = " << min(As_).value()
         << ", max = " << max(As_).value() << nl;

    Y_.resize(2);

    Y_.set
    (
        0,
        new volScalarField
        (
            IOobject
            (
                speciesAName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    Y_.set
    (
        1,
        new volScalarField
        (
            IOobject
            (
                speciesBName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
    );

    readMineralProperties();

    Info << "  speciesA : " << speciesAName_ << nl
         << "  speciesB : " << speciesBName_ << nl
         << "  k        = " << k_.value() << " (intrinsic rate constant)" << nl
         << "  Ksp      = " << Ksp_.value() << nl
         << "  Di       = " << Di_.value() << nl
         << "OK" << nl << endl;
}


// -------------------------------------------------------------------------//

void Foam::geochemicalModels::mixingPrecipitation::readMineralProperties()
{
    forAll(mineralList_, s)
    {
        word currentMineral = mineralList_[s];
        Info << "  Reading properties for mineral: " << currentMineral << endl;

        mineralSubDict_.set
        (
            s,
            new dictionary
            (
                transportPropertiesDict_.subDict(currentMineral + "Properties")
            )
        );

        Vm_.set
        (
            s,
            new dimensionedScalar
            (
                mineralSubDict_[s].lookup("Vm")
            )
        );
    }
}


void Foam::geochemicalModels::mixingPrecipitation::updateFluidComposition()
{
    word divScheme = "div(phi,Yi)";
    const volTensorField& Deff = effectiveDispersionTensor();

    volScalarField& C_A = Y_[0];
    volScalarField& C_B = Y_[1];

    // Effective volumetric rate = k * As [1/s]
    // k [m/s] (intrinsic rate constant, numerically mol/m^2/s)
    // As [1/m] (specific surface area = S_react / V_cell)
    volScalarField kAs("kAs", k_ * As_);

    // Saturation ratio: sigma = C_A * C_B / Ksp
    volScalarField sigma
    (
        "sigma",
        C_A * C_B / Ksp_
    );

    // Precipitation mask: 1 where supersaturated (sigma >= 1), 0 elsewhere
    volScalarField mask
    (
        "mask",
        pos(sigma - dimensionedScalar("one", dimless, 1.0))
    );

    // --- Solve species A (using C_B from previous iteration) ---
    //
    // Source = -kAs * (C_A*C_B/Ksp - 1) * mask
    //        = (-kAs/Ksp * C_B * mask) * C_A  +  kAs * mask
    //
    // Semi-implicit linearization:
    //   Sp coefficient = -kAs/Ksp * C_B * mask  (<= 0, stable implicit sink)
    //   Su source      = kAs * mask             (>= 0, explicit source)
    {
        fvScalarMatrix CAEqn
        (
            fvm::ddt(eps_, C_A)
          + fvm::div(phi_, C_A, divScheme)
          - fvm::laplacian(Deff, C_A, "laplacian(Di,Yi)")
          ==
            fvm::Sp(-(kAs / Ksp_) * C_B * mask, C_A)
          + kAs * mask
        );
        CAEqn.solve();
    }

    // --- Solve species B (using updated C_A) ---
    //
    // Recompute sigma and mask with updated C_A
    sigma = C_A * C_B / Ksp_;
    mask = pos(sigma - dimensionedScalar("one", dimless, 1.0));

    // Recompute kAs (As_ is constant, but kept for clarity)
    kAs = k_ * As_;

    {
        fvScalarMatrix CBEqn
        (
            fvm::ddt(eps_, C_B)
          + fvm::div(phi_, C_B, divScheme)
          - fvm::laplacian(Deff, C_B, "laplacian(Di,Yi)")
          ==
            fvm::Sp(-(kAs / Ksp_) * C_A * mask, C_B)
          + kAs * mask
        );
        CBEqn.solve();
    }
}


void Foam::geochemicalModels::mixingPrecipitation::updateMineralDistribution()
{
    volScalarField& C_A = Y_[0];
    volScalarField& C_B = Y_[1];

    // Effective volumetric rate = k * As [1/s]
    volScalarField kAs("kAs", k_ * As_);

    // Rate of precipitation (only when supersaturated)
    // R_vol = k * As * max(sigma - 1, 0)  [1/s]
    volScalarField sigma("sigma", C_A * C_B / Ksp_);
    volScalarField R
    (
        "R",
        kAs * max
        (
            sigma - dimensionedScalar("one", dimless, 1.0),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    forAll(Ys_, s)
    {
        // dYs/dt = R * Vm  (volume fraction grows as mineral precipitates)
        solve
        (
            fvm::ddt(Ys_[s]) == R * Vm_[s]
        );

        // Clamp mineral volume fraction
        Ys_[s].max(0.0);
        Ys_[s].min(0.999);
    }
}


// ************************************************************************* //
