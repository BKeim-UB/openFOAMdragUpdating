/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
  \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "carnahanStarlingSolidDrag.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //   

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(carnahanStarlingSolidDrag, 0);
    addToRunTimeSelectionTable(dragModel, carnahanStarlingSolidDrag, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //     


Foam::tmp<Foam::volScalarField>
Foam::dragModels::carnahanStarlingSolidDrag::KSolidSolid
(
    const phaseModel& gas,
    const phaseModel& solid1,
    const phaseModel& solid2
) const
{

    const phaseSystem& fluid = gas.fluid();
    const volScalarField& alphag = gas;
    const volScalarField& alphas1 = solid1;
    const volScalarField& alphas2 = solid2;
    const scalar Pi = constant::mathematical::pi;
    const volScalarField magURel(mag(solid1.U() - solid2.U()));
    const volScalarField alphasMax = fluid.alfasMax();
    const word&continuousPhaseName = fluid.continuousPhaseName();
    const phaseModel& continuousPhase = fluid.phases()[continuousPhaseName];

    volScalarField alphas = 1.0 - alphag;
    volScalarField g0 = 1.0/(1 - cbrt(alphas/alphasMax));

    volScalarField const_sum = 0.0/solid1.d();
    
    //----------------------------------------------
    //inserted
    // --- MFIX-style RDF calculation for a TWO-SOLID system ---                                                                                            
    const volScalarField zeta3 = alphas1 + alphas2;
    const volScalarField zeta2 = (alphas1 / solid1.d()) + (alphas2 / solid2.d());

    const volScalarField term1 = 1.0 / (1.0 - zeta3 ); // removed small here 
							    const volScalarField term2 = (3.0 * solid1.d() * solid2.d() / (solid1.d() + solid2.d())) * (zeta2 / sqr(1.0 - zeta3)); // and here
    const volScalarField term3 = (2.0 * sqr(solid1.d() * solid2.d() / (solid1.d() + solid2.d()))) * (sqr(zeta2) / pow(1.0 - zeta3, 3)); // and here

    const volScalarField g0_ls = term1 + term2 + term3;
    // ----------------------------------------------------------------------
    
    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if (&phase != &continuousPhase)
        {
            const volScalarField& alpha = phase;
            const_sum += alpha / phase.d();
        }

    }

    // --- Final drag calculation ---                                                                 
    volScalarField fractNum =
          3.0 * (1.0 + E_) * (Pi / 2.0 + Cf_ * sqr(Pi) / 8.0)
        * alphas1 * solid1.rho() * alphas2 * solid2.rho()
        * sqr(solid1.d() + solid2.d())
        * g0_ls * magURel;

    volScalarField fractDen = 2.0 * Pi * (solid1.rho() * pow(solid1.d(), 3.0) +  solid2.rho() * pow(solid2.d(), 3.0));

    return (fractNum / (fractDen));

}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //       

Foam::dragModels::carnahanStarlingSolidDrag::carnahanStarlingSolidDrag
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dragModel(dict, interface, registerObject),
    interface_(interface),
    solid1Name_(dict.lookup("solid1")),
    solid2Name_(dict.lookup("solid2")),
    gasName_(dict.lookup("gas")),
    E_("E", dimless, dict),
    Cf_("Cf", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::carnahanStarlingSolidDrag::~carnahanStarlingSolidDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
// Main K function
Foam::tmp<Foam::volScalarField>
Foam::dragModels::carnahanStarlingSolidDrag::K() const
{
  //  const phaseModel& solid1 = mesh().lookup<phaseModel>(solid1Name_);
  //  const phaseModel& solid2 = mesh().lookup<phaseModel>(solid2Name_);

    const phaseModel& gas = interface_.fluid().phases()[gasName_];
    const phaseModel& solid1 = interface_.fluid().phases()[solid1Name_];
    const phaseModel& solid2 = interface_.fluid().phases()[solid2Name_];

     if (interface_.contains(solid1) && interface_.contains(solid2))
    {
        return KSolidSolid(gas, solid1, solid2);
    }

    FatalErrorInFunction
        << "The interface " << interface_.name() << " does not contain two "
        << "out of the gas, liquid and solid phase models."
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::dragModels::carnahanStarlingSolidDrag::Kf() const
{
    return fvc::interpolate(K());
}

// ************************************************************************* //
