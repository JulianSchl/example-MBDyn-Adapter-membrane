/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = b06b3cd26b12c13c97018c1738b3a225b7a4e089
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void oscillatingLinearInlet_b06b3cd26b12c13c97018c1738b3a225b7a4e089(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    oscillatingLinearInletFixedValueFvPatchVectorField
);


const char* const oscillatingLinearInletFixedValueFvPatchVectorField::SHA1sum =
    "b06b3cd26b12c13c97018c1738b3a225b7a4e089";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089"
            " from patch/DimensionedField\n";
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089"
            " from patch/DimensionedField/mapper\n";
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089"
            " from patch/dictionary\n";
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089"
            " as copy\n";
    }
}


oscillatingLinearInletFixedValueFvPatchVectorField::
oscillatingLinearInletFixedValueFvPatchVectorField
(
    const oscillatingLinearInletFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

oscillatingLinearInletFixedValueFvPatchVectorField::
~oscillatingLinearInletFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void oscillatingLinearInletFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs oscillatingLinearInlet sha1: b06b3cd26b12c13c97018c1738b3a225b7a4e089\n";
    }

//{{{ begin code
    #line 29 "/home/julian/software/mbdyn-adapter/tutorials/cavityFSI/fluid/0/U/boundaryField/inlet"
const vector axis(0, 1, 0);
            vectorField v(this->patch().Cf().size());
            scalar t(this->db().time().value());
            forAll(this->patch().Cf(),patchid){
                vector vi((this->patch().Cf()[patchid].y()-0.875)/0.125*(1.0-cos(2.0*3.14159*t/5.0)),0,0);
                v[patchid] = vi;
            }
            operator==(v);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

