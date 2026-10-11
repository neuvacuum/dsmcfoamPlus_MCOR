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
    dsmcAzimuthalZoneFill

Description
    在指定 cellZone 内初始化 DSMC 粒子，并为每个新粒子施加绕给定转轴的局部周向速度。

    该模型的设计目标是服务于 MCOR1_0 旋转周期 AMI 验证算例。原有 dsmcZoneFill 只能
    叠加一个常量笛卡尔速度，这会把“周期面传输”与“壁面邻近单元的初始推进方向”混在一起。
    对于绕 x 轴的旋转周期验证，更合理的初始平均速度应满足：

        U_theta = Uaz * e_theta

    其中：

        e_theta = (omegaHat x r_perp) / |omegaHat x r_perp|

        r_perp = (x - x_c) - ((x - x_c) · omegaHat) omegaHat

    这样每个粒子的平均速度方向都会跟随它在截面中的实际角位置改变，从而避免把大量粒子
    从 t = 0 开始就错误地推向 innerwall / outerwall。

SourceFiles
    dsmcAzimuthalZoneFill.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcAzimuthalZoneFill_H
#define dsmcAzimuthalZoneFill_H

#include "dsmcAllConfigurations.H"
#include "vector.H"
#include "point.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dsmcAzimuthalZoneFill Declaration
\*---------------------------------------------------------------------------*/

class dsmcAzimuthalZoneFill
:
    public dsmcConfiguration
{
private:

    // Private data

        //- 需要初始化的 cellZone 名称。
        //  这样可以把粒子限制在一个更干净的核心区域，而不是整个网格全域填充。
        word zoneName_;

        //- 旋转轴方向的单位向量。
        //  物理上它定义了局部周向基向量 e_theta 的右手定则方向。
        vector rotationAxisHat_;

        point rotationCentre_;

        //- 周向速度的标量幅值。
        //  真正施加给粒子的平均速度是 azimuthalSpeed_ * e_theta。
        scalar azimuthalSpeed_;

        //- 可选的背景平移速度。
        //  该项用于保留扩展能力，例如将来若想在周向速度上再叠加轴向平移，
        //  可以直接在字典里给出一个常量 backgroundVelocity。
        vector backgroundVelocity_;


    // Private Member Functions

        //- 根据粒子位置计算局部平均速度。
        //  公式对应：
        //      U_mean = U_background + Uaz * e_theta
        //  其中：
        //      e_theta = (omegaHat x r_perp) / |omegaHat x r_perp|
        vector meanVelocityAtPosition(const point& p) const;


public:

    //- Runtime type information
    TypeName("dsmcAzimuthalZoneFill");

    // Constructors

        //- Construct from components
        dsmcAzimuthalZoneFill
        (
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        ~dsmcAzimuthalZoneFill();


    // Member Functions

        void setInitialConfiguration();
};


} // End namespace Foam

#endif

// ************************************************************************* //
