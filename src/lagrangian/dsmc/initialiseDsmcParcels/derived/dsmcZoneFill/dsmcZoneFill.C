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

#include "dsmcZoneFill.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcZoneFill, 0);

addToRunTimeSelectionTable(dsmcConfiguration, dsmcZoneFill, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcZoneFill::dsmcZoneFill
(
    dsmcCloud& cloud,
    const dictionary& dict
//     const word& name
)
:
    dsmcConfiguration(cloud, dict)
{

}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcZoneFill::~dsmcZoneFill()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void dsmcZoneFill::setInitialConfiguration()
{
    Info<< nl << "Initialising particles" << endl;

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

    const vector velocity(dsmcInitialiseDict_.lookup("velocity"));

    // 提前读取 zoneName，既用于后续 cellZone 查找，也用于 velocityFrame 报错信息。
    const word regionName(dsmcInitialiseDict_.lookup("zoneName"));

    // velocityFrame 是本开发阶段新增的可选项，默认 local 以保持历史算例不变。
    // local    ：velocity 直接作为 parcel 存储速度使用；
    // absolute ：velocity 表示实验室绝对系速度，若 cell 在旋转域内，则逐粒子转换为
    //            Urel = Uabs - Omega ^ (x - origin) 后再写入 parcel；
    // solidBody：velocity 表示叠加在刚体共转速度上的平移速度，即
    //            Uabs = velocity + Omega ^ (x - origin)。这样 rotor 中存储
    //            velocity，stator 中存储完整 Uabs，可一致初始化共转流体。
    const word velocityFrame =
        dsmcInitialiseDict_.lookupOrDefault<word>("velocityFrame", word("local"));

    // 初始速度参考系一旦拼错，会直接导致整片 rotor 的初始角向速度错误。
    // 因此这里不做静默回退，而是立即中止初始化。
    if
    (
        velocityFrame != "local"
     && velocityFrame != "absolute"
     && velocityFrame != "solidBody"
    )
    {
        FatalErrorIn("dsmcZoneFill::setInitialConfiguration()")
            << "velocityFrame must be local, absolute or solidBody for zone "
            << regionName
            << ", but got " << velocityFrame << nl
            << "in: " << mesh_.time().system()/"dsmcInitialiseDict"
            << exit(FatalError);
    }

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict_.subDict("numberDensities")
    );

    wordList molecules(numberDensitiesDict.toc());

    scalarField numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = readScalar
        (
            numberDensitiesDict.lookup(molecules[i])
        );
    }

    //numberDensities /= cloud_.nParticles();

    const cellZoneMesh& cellZones = mesh_.cellZones();
    label zoneId = cellZones.findZoneID(regionName);

    if(zoneId == -1)
    {
        FatalErrorIn("dsmcZoneFill::setInitialConfiguration()")
            << "Cannot find region: " << regionName << nl << "in: "
            << mesh_.time().system()/"dsmcInitialiseDict"
            << exit(FatalError);
    }

    const cellZone& zone = cellZones[zoneId];

    if (zone.size())
    {
        Info << "Lattice in zone: " << regionName << endl;

        forAll(zone, c)
        {
            const label& cellI = zone[c];

            List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
            (
                mesh_,
                cellI
            );

            forAll(cellTets, tetI)
            {
                const tetIndices& cellTetIs = cellTets[tetI];

                tetPointRef tet = cellTetIs.tet(mesh_);

                scalar tetVolume = tet.mag();

                forAll(molecules, i)
                {
                    const word& moleculeName(molecules[i]);

                    label typeId(findIndex(cloud_.typeIdList(), moleculeName));

                    if (typeId == -1)
                    {
                        FatalErrorIn("Foam::dsmcCloud<dsmcParcel>::initialise")
                            << "typeId " << moleculeName << "not defined." << nl
                            << abort(FatalError);
                    }

                    const dsmcParcel::constantProperties& cP = cloud_.constProps(typeId);

                    scalar numberDensity = numberDensities[i];

                    // Calculate the number of particles required
                    scalar particlesRequired = numberDensity*tetVolume;

                    particlesRequired /= cloud_.nParticles(cellI);

                    // Only integer numbers of particles can be inserted
                    label nParticlesToInsert = label(particlesRequired);

                    // Add another particle with a probability proportional to the
                    // remainder of taking the integer part of particlesRequired
                    if
                    (
                        (particlesRequired - nParticlesToInsert)
                            > rndGen_.sample01<scalar>()
                    )
                    {
                        nParticlesToInsert++;
                    }

                    for (label pI = 0; pI < nParticlesToInsert; pI++)
                    {
                        point p = tet.randomPoint(rndGen_);

                        vector U = cloud_.equipartitionLinearVelocity
                        (
                            translationalTemperature,
                            cP.mass()
                        );

                        scalar ERot = cloud_.equipartitionRotationalEnergy
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

                        label ELevel = cloud_.equipartitionElectronicLevel
                        (
                            electronicTemperature,
                            cP.electronicDegeneracyList(),
                            cP.electronicEnergyList()
                        );

                        // 把用户给定的 bulk velocity 转成 parcel 实际存储参考系。
                        // 在非旋转 cell 中 local 与 absolute 完全相同。
                        // 在旋转 cell 中，DSMC 运动方程使用相对速度 Urel，
                        // 所以 absolute 初值必须扣除局部刚体速度：
                        //     Urel = Uabs - Omega ^ (p - origin)
                        vector bulkVelocity = velocity;

                        if (velocityFrame == "solidBody")
                        {
                            const vector Urot =
                                cloud_.rotatingFrameOmega()
                              ^ (p - cloud_.rotatingFrameOrigin());

                            // rotor parcel 存相对速度，stator parcel 存绝对速度。
                            // 两者恢复到实验室系后都满足 Uabs=velocity+Omega^r。
                            if (!cloud_.isInRotatingRegion(cellI))
                            {
                                bulkVelocity += Urot;
                            }
                        }
                        else if
                        (
                            velocityFrame == "absolute"
                         && cloud_.isInRotatingRegion(cellI)
                        )
                        {
                            const vector r =
                                p - cloud_.rotatingFrameOrigin();
                            const vector Urot =
                                cloud_.rotatingFrameOmega() ^ r;
                            bulkVelocity = velocity - Urot;
                        }

                        U += bulkVelocity;

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
                    }
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const dsmcParcel::constantProperties& cP = cloud_.constProps
    (
        mostAbundantType
    );


    forAll(zone, c)
    {
        const label& cellI = zone[c];

        cloud_.sigmaTcRMax().primitiveFieldRef()[cellI] = cP.sigmaT()*cloud_.maxwellianMostProbableSpeed
        (
            translationalTemperature,
            cP.mass()
        );
    }

    cloud_.sigmaTcRMax().correctBoundaryConditions();
}


} // End namespace Foam

// ************************************************************************* //
