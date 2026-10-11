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

#include "dsmcDiffuseSpecularWallRotationFieldPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcDiffuseSpecularWallRotationFieldPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcDiffuseSpecularWallRotationFieldPatch,
    dictionary
);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationFieldPatch::setProperties()
{
    dsmcDiffuseSpecularWallPatch::setProperties();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseSpecularWallRotationFieldPatch::dsmcDiffuseSpecularWallRotationFieldPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    dsmcRotationPatchBoundary(t, mesh, cloud, dict),
    dsmcFieldPatchBoundary(t, mesh, cloud, dict),
    dsmcDiffuseSpecularWallPatch(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties"))

{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseSpecularWallRotationFieldPatch::~dsmcDiffuseSpecularWallRotationFieldPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseSpecularWallRotationFieldPatch::initialConfiguration()
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::calculateProperties()
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    measurePropertiesBeforeControl(p);

    // 先计算当前碰撞点的局部运动壁面速度：
    // Uw = Omega*(eOmega ^ (x - x0))。
    // Field 版本只改变壁温来源，不应改变旋转壁面速度与反射物理。
    const vector localPatchVelocity =
        dsmcRotationPatchBoundary::wallVelocity(p);

    if (diffuseFraction() > cloud_.rndGen().sample01<scalar>())
    {
        // 从边界场读取局部壁温；该温度只影响漫反射 Maxwellian 采样宽度。
        const scalar& localPatchTemperature =
            dsmcFieldPatchBoundary::patchLocalTemperature(p);

        // 漫反射：在局部壁面静止系采样热速度，再叠加 localPatchVelocity 回绝对系。
        dsmcDiffuseWallPatch::performDiffuseReflection
        (
            p,
            localPatchTemperature,
            localPatchVelocity
        );
    }
    else
    {
        // 镜面反射必须在运动壁面参考系中完成：
        // Urel = Uabs - Uw。
        vector& U = p.U();
        vector relativeVelocity = U - localPatchVelocity;

        // 使用当前 parcel 记录的壁面法向，并保持与基础 specular 模型一致的方向约定。
        vector nw = p.normal();
        nw /= mag(nw);

        // 只翻转相对速度的法向分量：
        // Urel,new = Urel,old - 2*(Urel,old · n)*n。
        const scalar relativeNormalVelocity = relativeVelocity & nw;

        if (relativeNormalVelocity > 0.0)
        {
            relativeVelocity -= 2.0*relativeNormalVelocity*nw;
        }

        // 转回绝对速度，供后续 DSMC 追踪和碰撞使用。
        U = relativeVelocity + localPatchVelocity;

        // 保持原始 specularInteraction() 统计入口，避免破坏壁面通量统计。
        cloud_.porousMeas().specularInteraction(p, nw);
    }

    measurePropertiesAfterControl(p);
}


void dsmcDiffuseSpecularWallRotationFieldPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcDiffuseSpecularWallRotationFieldPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}


} // End namespace Foam

// ************************************************************************* //
