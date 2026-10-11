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

\*---------------------------------------------------------------------------*/

#include "dsmcAdaptiveCellWeightDiagnostics.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"
#include "OFstream.H"
#include "PstreamReduceOps.H"

using namespace Foam::constant;

namespace Foam
{

defineTypeNameAndDebug(dsmcAdaptiveCellWeightDiagnostics, 0);

addToRunTimeSelectionTable
(
    dsmcField,
    dsmcAdaptiveCellWeightDiagnostics,
    dictionary
);

namespace
{

scalar parcelVibrationalEnergy
(
    const dsmcCloud& cloud,
    const dsmcParcel& parcel
)
{
    // 逐模态重建当前 parcel 的振动能。
    // 这里必须显式调用 constantProperties::eVib_m()，因为后续 split/merge
    // 的守恒诊断需要基于“真实储存在量子能级里的振动能”而不是仅统计一个
    // 平均温度。
    const dsmcParcel::constantProperties& constProp =
        cloud.constProps(parcel.typeId());

    scalar energy = 0.0;

    forAll(parcel.vibLevel(), mode)
    {
        energy += constProp.eVib_m(mode, parcel.vibLevel()[mode]);
    }

    return energy;
}


scalar parcelElectronicEnergy
(
    const dsmcCloud& cloud,
    const dsmcParcel& parcel
)
{
    // 电子能这里直接按当前离散能级查表。
    // 这与后续反应/碰撞模型里对 ELevel 的定义保持一致，便于把守恒监测
    // 和真实反应状态一一对应。
    const dsmcParcel::constantProperties& constProp =
        cloud.constProps(parcel.typeId());

    const scalarList& levels = constProp.electronicEnergyList();
    const label levelId = parcel.ELevel();

    if (levelId >= 0 && levelId < levels.size())
    {
        return levels[levelId];
    }

    return 0.0;
}

} // End unnamed namespace


dsmcAdaptiveCellWeightDiagnostics::dsmcAdaptiveCellWeightDiagnostics
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    sampleInterval_(1),
    sampleCounter_(0),
    nTimeSteps_(0.0),
    fieldName_(propsDict_.lookup("fieldName")),
    speciesIds_(),
    targetParticlesPerSubCell_(10.0),
    assumedSubCellsPerCell_(1.0),
    writeDiagnosticDictionary_(true),
    particlesPerCell_
    (
        IOobject
        (
            "particlesPerCell_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    realMoleculesPerCell_
    (
        IOobject
        (
            "realMoleculesPerCell_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cellWeightFactor_
    (
        IOobject
        (
            "cellWeightFactor_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    weightedRhoN_
    (
        IOobject
        (
            "weightedRhoN_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimVolume, 0.0)
    ),
    weightedRhoM_
    (
        IOobject
        (
            "weightedRhoM_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume, 0.0)
    ),
    weightedMomentum_
    (
        IOobject
        (
            "weightedMomentum_" + fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -2, -1, 0, 0),
            vector::zero
        )
    ),
    particlesPerCellCum_(mesh_.nCells(), 0.0),
    realMoleculesPerCellCum_(mesh_.nCells(), 0.0),
    weightedMassCum_(mesh_.nCells(), 0.0),
    weightedMomentumCum_(mesh_.nCells(), vector::zero)
{}


dsmcAdaptiveCellWeightDiagnostics::~dsmcAdaptiveCellWeightDiagnostics()
{}


void dsmcAdaptiveCellWeightDiagnostics::updateConfiguration()
{
    // 先建立一组硬编码默认值。
    // 这样即便当前 case 尚未写 adaptiveProperties，本模块也能以一个明确、
    // 可解释的保守配置启动，而不是读到未初始化垃圾值。
    sampleInterval_ = 1;
    targetParticlesPerSubCell_ = 10.0;
    assumedSubCellsPerCell_ = 1.0;
    writeDiagnosticDictionary_ = true;

    // 然后读取 dsmcProperties 里的 adaptiveProperties 作为“全局默认配置”。
    // 这一层表达的是求解器/算例级意图，例如整个 case 希望采用多少目标
    // 代表粒子数、是否默认开启守恒摘要。
    if (cloud_.particleProperties().found("adaptiveProperties"))
    {
        const dictionary& adaptiveDict =
            cloud_.particleProperties().subDict("adaptiveProperties");

        targetParticlesPerSubCell_ =
            adaptiveDict.lookupOrDefault<scalar>
            (
                "targetParticlesPerSubCell",
                targetParticlesPerSubCell_
            );

        assumedSubCellsPerCell_ =
            adaptiveDict.lookupOrDefault<scalar>
            (
                "assumedSubCellsPerCell",
                assumedSubCellsPerCell_
            );

        writeDiagnosticDictionary_ =
            adaptiveDict.lookupOrDefault<bool>
            (
                "enableConservationMonitor",
                writeDiagnosticDictionary_
            );
    }

    // 最后再用 fieldPropertiesDict 中当前 field 的本地选项做覆盖。
    // 这样同一个算例就可以：
    // 1. 在 dsmcProperties 里定义全局默认 targetParticlesPerSubCell；
    // 2. 仅对某个特定 field（例如 mixture 或 rotorZone）做局部细调。
    sampleInterval_ = max
    (
        label(1),
        propsDict_.lookupOrDefault<label>("sampleInterval", sampleInterval_)
    );

    targetParticlesPerSubCell_ =
        propsDict_.lookupOrDefault<scalar>
        (
            "targetParticlesPerSubCell",
            targetParticlesPerSubCell_
        );

    assumedSubCellsPerCell_ =
        propsDict_.lookupOrDefault<scalar>
        (
            "assumedSubCellsPerCell",
            assumedSubCellsPerCell_
        );

    writeDiagnosticDictionary_ =
        propsDict_.lookupOrDefault<bool>
        (
            "writeDiagnosticDictionary",
            writeDiagnosticDictionary_
        );

    // 目标代表粒子数和假定 sub-cell 数若非正，目标权重场就失去物理意义。
    // 这里直接中止，而不是静默夹紧到一个很小正数，避免用户误以为
    // 输出字段已经可信。
    if (targetParticlesPerSubCell_ <= SMALL)
    {
        FatalErrorIn
        (
            "dsmcAdaptiveCellWeightDiagnostics::updateConfiguration()"
        )
            << "targetParticlesPerSubCell must be positive, but got "
            << targetParticlesPerSubCell_ << nl
            << "for field " << fieldName_
            << exit(FatalError);
    }

    if (assumedSubCellsPerCell_ <= SMALL)
    {
        FatalErrorIn
        (
            "dsmcAdaptiveCellWeightDiagnostics::updateConfiguration()"
        )
            << "assumedSubCellsPerCell must be positive, but got "
            << assumedSubCellsPerCell_ << nl
            << "for field " << fieldName_
            << exit(FatalError);
    }
}


void dsmcAdaptiveCellWeightDiagnostics::initialiseSpeciesIds()
{
    const List<word>& species(propsDict_.lookup("typeIds"));

    DynamicList<word> uniqueSpecies(0);

    forAll(species, i)
    {
        const word& speciesName = species[i];

        if (findIndex(uniqueSpecies, speciesName) == -1)
        {
            uniqueSpecies.append(speciesName);
        }
    }

    uniqueSpecies.shrink();
    speciesIds_.setSize(uniqueSpecies.size(), -1);

    forAll(uniqueSpecies, i)
    {
        const label typeId = findIndex(cloud_.typeIdList(), uniqueSpecies[i]);

        if (typeId == -1)
        {
            FatalErrorIn
            (
                "dsmcAdaptiveCellWeightDiagnostics::initialiseSpeciesIds()"
            )
                << "Cannot find typeId " << uniqueSpecies[i]
                << " in cloud typeIdList for field " << fieldName_ << nl
                << "in: " << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        speciesIds_[i] = typeId;
    }
}


void dsmcAdaptiveCellWeightDiagnostics::resetAccumulators()
{
    nTimeSteps_ = 0.0;
    sampleCounter_ = 0;
    particlesPerCellCum_ = 0.0;
    realMoleculesPerCellCum_ = 0.0;
    weightedMassCum_ = 0.0;
    weightedMomentumCum_ = vector::zero;
}


bool dsmcAdaptiveCellWeightDiagnostics::includeParcel
(
    const dsmcParcel& p
) const
{
    const label celli = p.cell();

    return
    (
        p.isFree()
     && celli >= 0
     && celli < mesh_.nCells()
     && findIndex(speciesIds_, p.typeId()) != -1
    );
}


scalar dsmcAdaptiveCellWeightDiagnostics::targetRepresentativeParcelsPerCell() const
{
    return targetParticlesPerSubCell_*assumedSubCellsPerCell_;
}


void dsmcAdaptiveCellWeightDiagnostics::writeConservationSummary() const
{
    scalar totalRepresentativeParcels = 0.0;
    scalar totalRealMolecules = 0.0;
    scalar totalMass = 0.0;
    scalar translationalEnergy = 0.0;
    scalar rotationalEnergy = 0.0;
    scalar vibrationalEnergy = 0.0;
    scalar electronicEnergy = 0.0;
    vector totalMomentum = vector::zero;

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& parcel = iter();

        if (!includeParcel(parcel))
        {
            continue;
        }

        const scalar representedMolecules =
            cloud_.representedMolecules(parcel);
        const dsmcParcel::constantProperties& constProp =
            cloud_.constProps(parcel.typeId());

        const scalar mass = constProp.mass();
        const vector& U = parcel.U();

        totalRepresentativeParcels += 1.0;
        totalRealMolecules += representedMolecules;
        totalMass += representedMolecules*mass;
        totalMomentum += representedMolecules*mass*U;
        translationalEnergy +=
            representedMolecules*0.5*mass*magSqr(U);
        rotationalEnergy += representedMolecules*parcel.ERot();
        vibrationalEnergy +=
            representedMolecules*parcelVibrationalEnergy(cloud_, parcel);
        electronicEnergy +=
            representedMolecules*parcelElectronicEnergy(cloud_, parcel);
    }

    if (Pstream::parRun())
    {
        reduce(totalRepresentativeParcels, sumOp<scalar>());
        reduce(totalRealMolecules, sumOp<scalar>());
        reduce(totalMass, sumOp<scalar>());
        reduce(totalMomentum, sumOp<vector>());
        reduce(translationalEnergy, sumOp<scalar>());
        reduce(rotationalEnergy, sumOp<scalar>());
        reduce(vibrationalEnergy, sumOp<scalar>());
        reduce(electronicEnergy, sumOp<scalar>());
    }

    const scalar totalEnergy =
        translationalEnergy
      + rotationalEnergy
      + vibrationalEnergy
      + electronicEnergy;

    const fileName uniformDir =
        mesh_.time().path()/mesh_.time().timeName()/"uniform";

    if (Pstream::master())
    {
        if (!isDir(uniformDir))
        {
            mkDir(uniformDir);
        }

        const fileName outputFile =
            uniformDir/("adaptiveConservation_" + fieldName_ + ".dict");

        OFstream os(outputFile);

        os  << "time " << mesh_.time().value() << ";" << nl
            << "fieldName " << fieldName_ << ";" << nl
            << "targetParticlesPerSubCell " << targetParticlesPerSubCell_
            << ";" << nl
            << "assumedSubCellsPerCell " << assumedSubCellsPerCell_
            << ";" << nl
            << "targetRepresentativeParcelsPerCell "
            << targetRepresentativeParcelsPerCell() << ";" << nl
            << "nEquivalentParticlesReference " << cloud_.nParticle()
            << ";" << nl
            << "totalRepresentativeParcels " << totalRepresentativeParcels
            << ";" << nl
            << "totalRealMolecules " << totalRealMolecules << ";" << nl
            << "totalMass " << totalMass << ";" << nl
            << "totalMomentum " << totalMomentum << ";" << nl
            << "translationalEnergy " << translationalEnergy << ";" << nl
            << "rotationalEnergy " << rotationalEnergy << ";" << nl
            << "vibrationalEnergy " << vibrationalEnergy << ";" << nl
            << "electronicEnergy " << electronicEnergy << ";" << nl
            << "totalEnergy " << totalEnergy << ";" << nl;
    }

    Info<< "Adaptive cell-weight diagnostics (" << fieldName_ << "): "
        << "totalRepresentativeParcels=" << totalRepresentativeParcels
        << " totalRealMolecules=" << totalRealMolecules
        << " totalMass=" << totalMass
        << " totalMomentum=" << totalMomentum
        << " totalEnergy=" << totalEnergy
        << endl;
}


void dsmcAdaptiveCellWeightDiagnostics::createField()
{
    Info<< "Initialising dsmcAdaptiveCellWeightDiagnostics field "
        << fieldName_ << endl;

    initialiseSpeciesIds();
    updateConfiguration();
    resetAccumulators();
}


void dsmcAdaptiveCellWeightDiagnostics::calculateField()
{
    sampleCounter_++;

    if (sampleCounter_ < sampleInterval_)
    {
        return;
    }

    sampleCounter_ = 0;
    nTimeSteps_ += 1.0;

    forAllConstIter(dsmcCloud, cloud_, iter)
    {
        const dsmcParcel& parcel = iter();

        if (!includeParcel(parcel))
        {
            continue;
        }

        const label celli = parcel.cell();

        // includeParcel() 已保证 celli 合法。
        // 这里保留显式局部变量，是为了让后续所有 cell 级累计都只经过这一个
        // 已验证索引，避免未来维护时把 parcel.cell() 原样散落到多处。
        const scalar representedMolecules =
            cloud_.representedMolecules(parcel);
        const scalar mass = cloud_.constProps(parcel.typeId()).mass();

        particlesPerCellCum_[celli] += 1.0;
        realMoleculesPerCellCum_[celli] += representedMolecules;
        weightedMassCum_[celli] += representedMolecules*mass;
        weightedMomentumCum_[celli] += representedMolecules*mass*parcel.U();
    }

    if (!time_.time().outputTime())
    {
        return;
    }

    const scalar nAvTimeSteps = max(nTimeSteps_, scalar(1.0));
    const scalar targetRepParcelsPerCell = targetRepresentativeParcelsPerCell();
    const scalar referenceEquivalentParticles = max(cloud_.nParticle(), SMALL);

    forAll(particlesPerCellCum_, celli)
    {
        const scalar cellVolume = mesh_.cellVolumes()[celli];

        particlesPerCell_[celli] =
            particlesPerCellCum_[celli]/nAvTimeSteps;

        realMoleculesPerCell_[celli] =
            realMoleculesPerCellCum_[celli]/nAvTimeSteps;

        weightedRhoN_[celli] =
            realMoleculesPerCellCum_[celli]/(nAvTimeSteps*cellVolume);

        weightedRhoM_[celli] =
            weightedMassCum_[celli]/(nAvTimeSteps*cellVolume);

        weightedMomentum_[celli] =
            weightedMomentumCum_[celli]/(nAvTimeSteps*cellVolume);

        if (targetRepParcelsPerCell > SMALL)
        {
            cellWeightFactor_[celli] =
                realMoleculesPerCell_[celli]
               /(targetRepParcelsPerCell*referenceEquivalentParticles);
        }
        else
        {
            cellWeightFactor_[celli] = 0.0;
        }
    }

    particlesPerCell_.write();
    realMoleculesPerCell_.write();
    cellWeightFactor_.write();
    weightedRhoN_.write();
    weightedRhoM_.write();
    weightedMomentum_.write();

    if (writeDiagnosticDictionary_)
    {
        writeConservationSummary();
    }

    if (time_.resetFieldsAtOutput())
    {
        resetAccumulators();
    }
}


void dsmcAdaptiveCellWeightDiagnostics::resetField()
{
    particlesPerCell_ = 0.0;
    realMoleculesPerCell_ = 0.0;
    cellWeightFactor_ = 0.0;
    weightedRhoN_ = 0.0;
    weightedRhoM_ = 0.0;
    weightedMomentum_ = vector::zero;

    resetAccumulators();
}


void dsmcAdaptiveCellWeightDiagnostics::writeField()
{}


void dsmcAdaptiveCellWeightDiagnostics::updateProperties
(
    const dictionary& newDict
)
{
    updateBasicFieldProperties(newDict);
}


scalar dsmcAdaptiveCellWeightDiagnostics::translationalT(const label)
{
    // 当前诊断模块只负责权重、守恒和动量密度，不生成温度定义。
    return -1.0;
}


scalar dsmcAdaptiveCellWeightDiagnostics::overallT(const label)
{
    // 当前诊断模块只负责权重、守恒和动量密度，不生成温度定义。
    return -1.0;
}

} // End namespace Foam

// ************************************************************************* //
