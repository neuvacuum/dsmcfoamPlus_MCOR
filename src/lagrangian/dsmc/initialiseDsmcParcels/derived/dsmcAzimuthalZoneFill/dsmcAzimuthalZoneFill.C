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
    在指定 zone 中初始化 DSMC 粒子，并为每个新粒子构造局部周向平均速度。

\*---------------------------------------------------------------------------*/

#include "dsmcAzimuthalZoneFill.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

namespace Foam
{

defineTypeNameAndDebug(dsmcAzimuthalZoneFill, 0);
addToRunTimeSelectionTable(dsmcConfiguration, dsmcAzimuthalZoneFill, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dsmcAzimuthalZoneFill::dsmcAzimuthalZoneFill
(
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    // 基类负责保留 mesh、cloud、初始化字典和随机数发生器。
    dsmcConfiguration(cloud, dict),

    // 下面这些成员全部从当前 configuration 子字典中直接读取。
    // 这样运行时选择到该模型后，就不需要在 setInitialConfiguration()
    // 里反复查询和解释同一组控制参数。
    zoneName_(dsmcInitialiseDict_.lookup("zoneName")),
    rotationAxisHat_(dsmcInitialiseDict_.lookup("rotationAxis")),
    rotationCentre_(dsmcInitialiseDict_.lookup("rotationCentre")),
    azimuthalSpeed_(readScalar(dsmcInitialiseDict_.lookup("azimuthalSpeed"))),
    backgroundVelocity_
    (
        dsmcInitialiseDict_.found("backgroundVelocity")
      ? vector(dsmcInitialiseDict_.lookup("backgroundVelocity"))
      : vector::zero
    )
{
    // 这里必须把 rotationAxis 归一化。
    // 因为后续所有局部周向方向、轴向投影和截面半径的几何计算，
    // 都建立在“rotationAxisHat_ 是单位向量”这一前提上。
    const scalar axisMag = mag(rotationAxisHat_);

    if (axisMag <= SMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude is too small in "
            << mesh_.time().system()/"dsmcInitialiseDict" << nl
            << "The azimuthal initializer requires a valid non-zero axis."
            << exit(FatalError);
    }

    rotationAxisHat_ /= axisMag;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcAzimuthalZoneFill::~dsmcAzimuthalZoneFill()
{}


// * * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * //

vector dsmcAzimuthalZoneFill::meanVelocityAtPosition(const point& p) const
{
    // 先构造相对旋转中心的位置向量：
    //     r = x - x_c
    const vector relativePosition = p - rotationCentre_;

    // 提取沿旋转轴方向的分量：
    //     r_parallel = (r · omegaHat) omegaHat
    const vector axialComponent = (relativePosition & rotationAxisHat_)*rotationAxisHat_;

    // 得到垂直于旋转轴的截面半径向量：
    //     r_perp = r - r_parallel
    const vector radialVector = relativePosition - axialComponent;

    // 理论局部周向方向：
    //     e_theta = (omegaHat x r_perp) / |omegaHat x r_perp|
    const vector tangentialVector = rotationAxisHat_ ^ radialVector;
    const scalar tangentialMag = mag(tangentialVector);

    if (tangentialMag <= SMALL)
    {
        return backgroundVelocity_;
    }

    const vector eTheta = tangentialVector/tangentialMag;

    // 最终平均速度：
    //     U_mean = U_background + Uaz * e_theta
    return backgroundVelocity_ + azimuthalSpeed_*eTheta;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcAzimuthalZoneFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles with azimuthal zone fill" << endl;

    const scalar translationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("translationalTemperature"))
    );

    const scalar rotationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("rotationalTemperature"))
    );

    const scalar vibrationalTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("vibrationalTemperature"))
    );

    const scalar electronicTemperature
    (
        readScalar(dsmcInitialiseDict_.lookup("electronicTemperature"))
    );

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict_.subDict("numberDensities")
    );

    wordList molecules(numberDensitiesDict.toc());
    scalarField numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar(numberDensitiesDict.lookup(molecules[i]));
    }

    // 在当前 polyMesh 中查找目标 cellZone。
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const label zoneId = cellZones.findZoneID(zoneName_);

    if (zoneId == -1)
    {
        FatalErrorInFunction
            << "Cannot find region: " << zoneName_ << nl << "in: "
            << mesh_.time().system()/"dsmcInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    if (!zone.size())
    {
        FatalErrorInFunction
            << "Zone " << zoneName_ << " is empty. "
            << "The azimuthal initializer needs a non-empty seeding zone."
            << exit(FatalError);
    }

    Info << "Lattice in zone: " << zoneName_ << endl;

    forAll(zone, c)
    {
        const label cellI = zone[c];

        // DSMC 仍按 cell tet 来随机撒点，保持和原有 dsmcZoneFill 同样的体积分配方式。
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices(mesh_, cellI);

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];
            tetPointRef tet = cellTetIs.tet(mesh_);
            const scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName = molecules[i];
                const label typeId = findIndex(cloud_.typeIdList(), moleculeName);

                if (typeId == -1)
                {
                    FatalErrorIn("Foam::dsmcCloud<dsmcParcel>::initialise")
                        << "typeId " << moleculeName << " not defined." << nl
                        << abort(FatalError);
                }

                const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);
                const scalar numberDensity = numberDensities[i];

                // 需要插入的宏观粒子数仍沿用原有体积分布公式：
                //     N_required = numberDensity * tetVolume / nEquivalentParticles(cell)
                scalar particlesRequired = numberDensity*tetVolume;
                particlesRequired /= cloud_.nParticles(cellI);

                label nParticlesToInsert = label(particlesRequired);

                if ((particlesRequired - nParticlesToInsert) > rndGen_.sample01<scalar>())
                {
                    nParticlesToInsert++;
                }

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    // 在 tet 内随机取一点，保证 zone 内体积分布均匀。
                    const point p = tet.randomPoint(rndGen_);

                    // 先生成热运动部分。
                    vector U = cloud_.equipartitionLinearVelocity
                    (
                        translationalTemperature,
                        cP.mass()
                    );

                    // 再叠加局部周向平均速度。
                    U += meanVelocityAtPosition(p);

                    const scalar ERot = cloud_.equipartitionRotationalEnergy
                    (
                        rotationalTemperature,
                        cP.rotationalDegreesOfFreedom()
                    );

                    labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                    (
                        vibrationalTemperature,
                        cP.nVibrationalModes(),
                        typeId
                    );

                    const label ELevel = cloud_.equipartitionElectronicLevel
                    (
                        electronicTemperature,
                        cP.electronicDegeneracyList(),
                        cP.electronicEnergyList()
                    );

                    label newParcel = -1;
                    label classification = 0;

                    const scalar& RWF = cloud_.coordSystem().RWF(cellI);

                    cloud_.addNewParcel
                    (
                        p,
                        U,
                        RWF,
                        ERot,
                        ELevel,
                        cellI,
                        cellTetIs.face(),
                        cellTetIs.tetPt(),
                        typeId,
                        newParcel,
                        classification,
                        vibLevel
                    );

                    // 记录当前 configuration 真正新增的 parcel 数量，
                    // 便于后续若需要在 initialise 阶段做更细统计时直接使用。
                    nParcelsAdded_++;
                }
            }
        }
    }

    // sigmaTcRMax_ 仍按当前 zone 内最丰富组分的热最概速初始化。
    // 这一步的物理意义和原 dsmcZoneFill 完全一致：
    //     sigmaTcRMax = sigma_T * c_mp(T)
    const label mostAbundantType = findMax(numberDensities);
    const dsmcParcel::constantProperties& cP = cloud_.constProps(mostAbundantType);

    forAll(zone, c)
    {
        const label cellI = zone[c];
        cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] =
            cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
            (
                translationalTemperature,
                cP.mass()
            );
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


} // End namespace Foam

// ************************************************************************* //
