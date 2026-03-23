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

Class
    Foam::geochemicalModels::mixingPrecipitation

Description
    Mixing-induced precipitation model for two-species systems.

    Solves transport of two dissolved species (A and B) with a nonlinear
    precipitation sink term based on the saturation ratio:

        R_vol = k * As * max(sigma - 1, 0)

    where sigma = C_A * C_B / Ksp is the saturation ratio,
    k is the intrinsic rate constant [mol/m^2/s] (dimension [m/s] in code),
    and As is the specific surface area [m^2/m^3] = [1/m].

    The reaction consumes equal moles of A and B (1:1 stoichiometry):
        dC_A/dt = -R_vol
        dC_B/dt = -R_vol
        dYs/dt  = R_vol * Vm

    Parameters in transportProperties:
        k    - intrinsic rate constant [mol/m^2/s] (code dimension: [m/s])
        Ksp  - solubility product in concentration units squared (e.g. mM^2)
        Di   - molecular diffusion coefficient [m^2/s]
        speciesA - name of species A field (default: "C_A")
        speciesB - name of species B field (default: "C_B")

    The reactive surface area is obtained from the surfaceAreaModel
    framework (via porousMedia_[0].surfaceArea()), allowing dynamic
    evolution when a suitable surfaceAreaModel is selected.

SourceFiles
    mixingPrecipitation.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPrecipitation_H
#define mixingPrecipitation_H

#include "volFields.H"
#include "speciesTable.H"
#include "runTimeSelectionTables.H"
#include "basicGeochemicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace geochemicalModels
{

/*---------------------------------------------------------------------------*\
                    Class mixingPrecipitation Declaration
\*---------------------------------------------------------------------------*/

class mixingPrecipitation : public basicGeochemicalModel
{

    dictionary mixingPrecipitationDict_;

  protected:

    dictionary transportPropertiesDict_;

    PtrList<dictionary> mineralSubDict_;

    PtrList<dimensionedScalar> Vm_;

    dimensionedScalar k_;
    dimensionedScalar Ksp_;
    dimensionedScalar Di_;

    word speciesAName_;
    word speciesBName_;

  // Private Member Functions

    void readMineralProperties();

    //- Disallow copy construct
    mixingPrecipitation(const mixingPrecipitation&);


  public:

    //- Runtime type information
      TypeName("mixingPrecipitation");


  // Constructors

    //- Construct from components
    mixingPrecipitation
    (
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Destructor
    virtual ~mixingPrecipitation()
    {}


    // Member Functions

    virtual void updateFluidComposition();

    virtual void updateMineralDistribution();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace geochemicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
