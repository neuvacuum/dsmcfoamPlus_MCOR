/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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
    dsmcKuangFreeStreamInflowPatch

Description

    Fully deletion boundary patch: atomisticMolecules get deleted instantly upon
    collision with a face on the boundary. (aka vacuum boundary condition)

    Selective deletion of one specie is possible. Other species get
    specularly reflected.

SourceFiles
    dsmcKuangFreeStreamInflowPatch.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcKuangFreeStreamInflowPatch_H
#define dsmcKuangFreeStreamInflowPatch_H

#include "dsmcGeneralBoundary.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dsmcKuangFreeStreamInflowPatch Declaration
\*---------------------------------------------------------------------------*/

class dsmcKuangFreeStreamInflowPatch
:
	public dsmcGeneralBoundary
{

private:

    // Private data
        dictionary propsDict_;

        labelList typeIds_;

        // Boundary temperatures
        scalar translationalTemperature_;
        scalar rotationalTemperature_;
        scalar vibrationalTemperature_;
        scalar electronicTemperature_;


        // list of number densities corresponding to the type Ids
        scalarField numberDensities_;

        // list of species
        List <scalarField> accumulatedParcelsToInsert_;


         scalar omega_;  // 粒子的角速度
         vector centerPoint_; // 旋转中心
         vector rotationAxis_; // 旋转轴

        // 速度参考系：
        // - local：velocity 与 omega 直接生成当前 cell 内 parcel 存储速度；
        // - absolute：velocity 与 omega 先生成实验室系速度，再按旋转域转成局部速度。
        word velocityFrame_;

        vector velocity_;

        // Private Member Functions
         void setProperties();

         // 根据边界给定的平移速度、角速度和速度参考系，计算局部 parcel 速度。
         vector reservoirVelocity(const point& p, const label cellI) const;

public:

    //- Runtime type information
    TypeName("dsmcKuangFreeStreamInflowPatch");

    // Constructors

    //- Construct from components
    dsmcKuangFreeStreamInflowPatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

// Destructor

    ~dsmcKuangFreeStreamInflowPatch();


// Member Functions

    void initialConfiguration();

    //- calculate any required properties
    void calculateProperties();

    //- control the mols
    void controlParcelsBeforeMove();
    void controlParcelsBeforeCollisions();
    void controlParcelsAfterCollisions();

    //- output data
    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    //- update properties from a modified dictionary
    void updateProperties(const dictionary&);

    //reset fields if mesh has changed
    void setNewBoundaryFields();




        // Access





};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
