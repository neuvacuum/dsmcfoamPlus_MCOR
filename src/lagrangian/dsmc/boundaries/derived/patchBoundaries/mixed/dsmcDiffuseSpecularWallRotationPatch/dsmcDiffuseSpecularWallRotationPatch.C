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

Description

\*---------------------------------------------------------------------------*/

#include "dsmcDiffuseSpecularWallRotationPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallRotationPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcDiffuseSpecularWallRotationPatch,
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationPatch::setProperties()
{
    dsmcDiffuseSpecularWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallRotationPatch::dsmcDiffuseSpecularWallRotationPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcRotationPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseSpecularWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallRotationPatch::~dsmcDiffuseSpecularWallRotationPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallRotationPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallRotationPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    // 先计算当前碰撞点的局部运动壁面速度：
    // Uw = Omega*(eOmega ^ (x - x0))。
    // 后续漫反射和镜面反射必须使用同一个 Uw，避免同一壁面模型中两条支路物理不一致。
    const vector localPatchVelocity =
        dsmcRotationPatchBoundary::wallVelocity(p);

    if (diffuseFraction() > cloud_.rndGen().sample01<scalar>())
    {
        // 漫反射：在壁面局部静止系中按壁温采样 Maxwellian 速度，
        // 然后由 performDiffuseReflection() 把 localPatchVelocity 加回绝对系速度。
        dsmcDiffuseWallPatch::performDiffuseReflection
        (
            p,
            0,
            localPatchVelocity
        );
    }
    else
    {
        // 镜面反射必须先转到运动壁面参考系：
        // Urel = Uabs - Uw。
        vector& U = p.U();
        vector relativeVelocity = U - localPatchVelocity;

        // 取壁面法向并单位化；该方向与原始 specular 壁面实现保持一致。
        vector nw = p.normal();
        nw /= mag(nw);

        // 在壁面参考系中只反射法向速度分量：
        // Urel,new = Urel,old - 2*(Urel,old · n)*n。
        const scalar relativeNormalVelocity = relativeVelocity & nw;

        if (relativeNormalVelocity > 0.0)
        {
            relativeVelocity -= 2.0*relativeNormalVelocity*nw;
        }

        // 转回绝对系：
        // Uabs,new = Urel,new + Uw。
        U = relativeVelocity + localPatchVelocity;

        // 保持与基础镜面壁面一致的壁面相互作用统计入口。
        cloud_.porousMeas().specularInteraction(p, nw);
    }

    measurePropertiesAfterControl(p);
}


void dsmcDiffuseSpecularWallRotationPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallRotationPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
