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

#include "dsmcRotationPatchBoundary.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcRotationPatchBoundary, 0);


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

vector dsmcRotationPatchBoundary::wallVelocity(const dsmcParcel& p)
{
    // 取粒子碰撞点相对旋转中心的位置矢量：
    // r = x - x0，其中 x 是当前 DSMC parcel 的壁面碰撞位置，x0 是旋转中心。
    const vector relativePosition = p.position() - centrePoint_;

    // 使用刚体转动速度公式：
    // Uw = Omega_vec ^ r = Omega * (eOmega ^ r)。
    // rotationAxis_ 已在构造函数中单位化，所以这里不再除以 |rotationAxis_|。
    const vector localWallVelocity =
        angularVelocity_*(rotationAxis_ ^ relativePosition);

    // 返回绝对系下的局部壁面速度，漫反射采样后会把该速度叠加到 Maxwellian 随机速度上。
    return localWallVelocity;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcRotationPatchBoundary::dsmcRotationPatchBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    centrePoint_
    (
        // 优先读取正确拼写 centrePoint；若旧算例仍使用 cetrePoint，也兼容读取。
        propsDict_.lookupOrDefault<vector>
        (
            "centrePoint",
            propsDict_.lookupOrDefault<vector>("cetrePoint", vector::zero)
        )
    ),
    rotationAxis_(propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0))),
    angularVelocity_
    (
        // 优先使用 Holweck 算例中实际写入的 Omega [rad/s]。
        propsDict_.lookupOrDefault<scalar>
        (
            "Omega",
            // 兼容小写 omega，便于与其他旋转参考系字典保持一致。
            propsDict_.lookupOrDefault<scalar>
            (
                "omega",
                // 兼容旧字段 wallVelocity；历史代码实际把它当角速度使用。
                propsDict_.lookupOrDefault<scalar>("wallVelocity", 0.0)
            )
        )
    )
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    // 旋转轴不能为零；否则 Uw = Omega*(eOmega ^ r) 中的单位方向 eOmega 不存在。
    const scalar axisMag = mag(rotationAxis_);

    if (axisMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_ << nl
            << "Please set dsmcRotationPatchBoundaryProperties.rotationAxis."
            << exit(FatalError);
    }

    // 单位化旋转轴，使 angularVelocity_ 保持纯角速度量纲 [rad/s]。
    rotationAxis_ /= axisMag;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcRotationPatchBoundary::~dsmcRotationPatchBoundary()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcRotationPatchBoundary::initialConfiguration()
{}


void dsmcRotationPatchBoundary::calculateProperties()
{}


void dsmcRotationPatchBoundary::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{}


void dsmcRotationPatchBoundary::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);
}


} // End namespace Foam

// ************************************************************************* //
