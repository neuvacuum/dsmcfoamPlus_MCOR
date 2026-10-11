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
    dsmcRotationPatchBoundary

Description
    Wall interaction setting microscopic velocity to a random one
    drawn from a Maxwellian distribution corresponding to a specified
    temperature

    FIELD VERSION: require boundaryU and boundaryT fields to be defined in
    your time directory.

SourceFiles
    dsmcRotationPatchBoundary.C


\*---------------------------------------------------------------------------*/

#ifndef dsmcRotationPatchBoundary_H
#define dsmcRotationPatchBoundary_H

#include "dsmcPatchBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dsmcRotationPatchBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcRotationPatchBoundary
:
	  public virtual dsmcPatchBoundary
{

private:

    // Private data

        dictionary propsDict_;


protected:

    // Protected data

        //- 旋转中心点，默认取全局原点；壁面局部速度使用 r = x - centrePoint_。
        const vector centrePoint_;

        //- 旋转轴方向；构造函数中会单位化，后续公式把它作为单位矢量 eOmega 使用。
        vector rotationAxis_;

        //- 角速度标量 Omega [rad/s]；壁面速度公式为 Uw = Omega*(eOmega ^ r)。
        const scalar angularVelocity_;


    // Protected member functions

        //- Wall velocity induced by the rotation, to be added to U
        vector wallVelocity(const dsmcParcel& p);


public:

    //- Runtime type information
    TypeName("dsmcRotationPatchBoundary");

    // Constructors

        //- Construct from components
        dsmcRotationPatchBoundary
		    (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
		    );

    // Destructor

 		~dsmcRotationPatchBoundary();


    // Member Functions

        void initialConfiguration();

        //- calculate any required properties
        void calculateProperties();

        //- control the mols
        virtual void controlParticle(dsmcParcel& p, dsmcParcel::trackingData& td) = 0;

        //- output data
        void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        );

        //- update properties from a modified dictionary
        void updateProperties(const dictionary&);


    // Access



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
