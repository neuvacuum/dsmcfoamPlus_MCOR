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

#include "dsmcCloud.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<dsmcParcel>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dsmcCloud::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id(typeIdList_[i]);

        Info<< tab << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = dsmcParcel::constantProperties(molDict);
    }
}

void Foam::dsmcCloud::initialiseRotatingReferenceFrame()
{
    rotatingReferenceFrameActive_ = false;
    // 每次重新初始化 SRF 配置时，都必须先把上一轮读取到的 cellZone
    // 名称显式清空。原因是下面使用的是 readIfPresent("cellZone", ...)，
    // 若字典热重读后删除了 cellZone 项，而这里又不先清空成员变量，
    // rotatingRegionZoneName_ 就会继续保留旧值，导致“本应回到全域 SRF”
    // 的配置仍然错误地沿用上一次的限域 zone。
    rotatingRegionZoneName_ = word::null;
    rotatingFrameSubCycles_ = 1;
    rotatingFrameMaxSubStep_ = -1.0;
    rotatingFrameMaxRotationAngle_ = Foam::degToRad(5.0);
    rotatingFrameMeshTolFraction_ = 0.25;
    rotatingRegionCellMask_.setSize(mesh_.nCells(), false);

    if (!particleProperties_.found("rotatingReferenceFrameProperties"))
    {
        return;
    }

    const dictionary& rotatingDict =
        particleProperties_.subDict("rotatingReferenceFrameProperties");

    const bool enabled =
        rotatingDict.lookupOrDefault<bool>("active", false);

    if (!enabled)
    {
        return;
    }

    rotatingDict.readIfPresent("cellZone", rotatingRegionZoneName_);
    rotatingFrameOmega_ =
        rotatingDict.lookupOrDefault<vector>("omega", vector::zero);
    rotatingFrameOrigin_ =
        rotatingDict.lookupOrDefault<vector>("origin", vector::zero);
    rotatingFrameSubCycles_ =
        max(label(1), rotatingDict.lookupOrDefault<label>("subCycles", 1));
    rotatingFrameMaxSubStep_ =
        rotatingDict.lookupOrDefault<scalar>("maxSubStep", -1.0);
    rotatingFrameMaxRotationAngle_ =
        rotatingDict.lookupOrDefault<scalar>
        (
            "maxRotationAngle",
            Foam::degToRad(5.0)
        );
    rotatingFrameMeshTolFraction_ =
        max
        (
            scalar(0),
            rotatingDict.lookupOrDefault<scalar>("meshToleranceFraction", 0.25)
        );

    if (rotatingFrameMaxSubStep_ <= ROOTVSMALL)
    {
        rotatingFrameMaxSubStep_ = -1.0;
    }

    if (rotatingFrameMaxRotationAngle_ <= ROOTVSMALL)
    {
        rotatingFrameMaxRotationAngle_ = GREAT;
    }

    if (rotatingRegionZoneName_.empty())
    {
        forAll(rotatingRegionCellMask_, celli)
        {
            rotatingRegionCellMask_[celli] = true;
        }

        Info<< "No rotating-reference-frame cellZone specified. Applying SRF "
            << "source terms to all " << rotatingRegionCellMask_.size()
            << " cells." << endl;
    }
    else
    {
        const cellZoneMesh& cellZones = mesh_.cellZones();
        const label zoneId = cellZones.findZoneID(rotatingRegionZoneName_);

        if (zoneId == -1)
        {
            // 这里不能再降级成“打印 warning 后静默关闭 SRF”。
            // 物理原因是：用户既然显式写了 active=true 且给出了 cellZone，
            // 就已经明确要求“仅在该区域内施加非惯性系源项”。
            // 若 zone 名称拼错或缺失而继续运行，求解会悄悄退回惯性系，
            // 阶段 1.1 的所有 SRF 验证都会被污染，且批处理时极易漏看 warning。
            // 因此这里改成 FatalError，使“显式限域 SRF 配置错误”立即失败。
            FatalErrorInFunction
                << "Rotating-reference-frame cellZone '"
                << rotatingRegionZoneName_ << "' was not found." << nl
                << "SRF is configured with active=true and an explicit cellZone,"
                << " so this case cannot silently fall back to the inertial frame."
                << exit(FatalError);
        }

        const labelList& zoneCells = cellZones[zoneId];

        forAll(zoneCells, i)
        {
            const label celli = zoneCells[i];

            if (celli >= 0 && celli < rotatingRegionCellMask_.size())
            {
                rotatingRegionCellMask_[celli] = true;
            }
        }
    }

    rotatingReferenceFrameActive_ = true;

    Info<< "Activated rotating reference frame";

    if (rotatingRegionZoneName_.empty())
    {
        Info<< " on all cells";
    }
    else
    {
        Info<< " for cellZone '" << rotatingRegionZoneName_ << "'";
    }

    Info<< " with omega=" << rotatingFrameOmega_
        << ", origin=" << rotatingFrameOrigin_
        << ", subCycles=" << rotatingFrameSubCycles_
        << ", maxSubStep=" << rotatingFrameMaxSubStep_
        << ", maxRotationAngle=" << rotatingFrameMaxRotationAngle_
        << ", meshToleranceFraction=" << rotatingFrameMeshTolFraction_
        << endl;
}

Foam::vector Foam::dsmcCloud::rotatingFrameAcceleration
(
    const point& position,
    const vector& Urel
) const
{
    const vector r = position - rotatingFrameOrigin_;

    return
        -2.0*(rotatingFrameOmega_ ^ Urel)
      - (rotatingFrameOmega_ ^ (rotatingFrameOmega_ ^ r));
}

Foam::scalar Foam::dsmcCloud::rotatingFrameTrackSubStep
(
    const label celli,
    const point& position,
    const vector& Urel,
    const scalar remainingTime
) const
{
    scalar dt = remainingTime;

    if (!isInRotatingRegion(celli))
    {
        return dt;
    }

    dt = min(dt, deltaTValue(celli)/scalar(rotatingFrameSubCycles_));

    if (rotatingFrameMaxSubStep_ > ROOTVSMALL)
    {
        dt = min(dt, rotatingFrameMaxSubStep_);
    }

    const scalar omegaMag = mag(rotatingFrameOmega_);

    if
    (
        omegaMag > ROOTVSMALL
     && rotatingFrameMaxRotationAngle_ < GREAT
    )
    {
        dt = min(dt, rotatingFrameMaxRotationAngle_/omegaMag);
    }

    if
    (
        rotatingFrameMeshTolFraction_ > ROOTVSMALL
     && celli >= 0
     && celli < mesh_.nCells()
    )
    {
        const scalar cellLength = cbrt(mesh_.cellVolumes()[celli]);
        const scalar accelMag = mag(rotatingFrameAcceleration(position, Urel));
        const scalar tolLength = rotatingFrameMeshTolFraction_*cellLength;

        if (accelMag > ROOTVSMALL && tolLength > ROOTVSMALL)
        {
            dt = min(dt, sqrt(2.0*tolLength/accelMag));
        }
    }

    return min(remainingTime, max(dt, ROOTVSMALL));
}

void Foam::dsmcCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, celli)
    {
        cellOccupancy_[celli].clear();
    }

    forAllIter(dsmcCloud, *this, iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }
}

void Foam::dsmcCloud::relocateStuckParcels()
{
    forAllIter(dsmcCloud, *this, iter)
    {
        if (iter().isStuck())
        {
            bool success = iter().relocateStuckParcel(mesh_);
            if (!success)
            {
                FatalErrorInFunction
                    << "Could not relocate stuck parcel!"
                    << exit(FatalError);
            }
        }
    }
}

void Foam::dsmcCloud::buildCellOccupancyFromScratch()
{
    initialiseRotatingReferenceFrame();

    cellOccupancy_.clear();
    cellOccupancy_.setSize(mesh_.nCells());

    buildCellOccupancy();
    relocateStuckParcels();
}

void Foam::dsmcCloud::buildCollisionSelectionRemainderFromScratch()
{
    collisionSelectionRemainder_.clear();
    collisionSelectionRemainder_.setSize(mesh_.nCells());

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, cO)
    {
        collisionSelectionRemainder_[cO] = rndGen_.sample01<scalar>();
    }
}

void Foam::dsmcCloud::resetBoundaries()
{
    boundaryMeas_.reset();
    boundaries_.setNewConfig();
}

void Foam::dsmcCloud::resetMeasurementTools()
{
    trackingInfo_.reset();
    fields_.resetFields();

    cellMeas_.reset();
}

void Foam::dsmcCloud::removeElectrons()
{
    forAll(cellOccupancy_, c)
    {
        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const dsmcParcel::constantProperties& constProp =
                constProps(p->typeId());

            const label& charge = constProp.charge();

            const scalar& RWF = coordSystem().RWF(c);

            const scalar mass = constProps(p->typeId()).mass();

            momentumMean_[c] += mass*RWF*p->U();
            rhoMMean_[c] += mass*RWF;

            if (charge == -1)
            {
                rhoNMeanElectron_[c] += 1.0*RWF;
                rhoMMeanElectron_[c] += mass*RWF;
                momentumMeanElectron_[c] += mass*RWF*p->U();
                linearKEMeanElectron_[c] += mass*RWF*(p->U() & p->U());

                //- found an electron
                deleteParticle(*p);
            }
        }
    }
}

void Foam::dsmcCloud::addElectrons()
{
    label electronTypeId = -1;

    //- find electron typeId
    forAll(constProps_, cP)
    {
        const label& particleCharge = constProps_[cP].charge();

        if (particleCharge == -1)
        {
            electronTypeId = cP;
            break;
        }
    }

    forAll(cellOccupancy_, c)
    {
        if (rhoMMeanElectron_[c] > VSMALL)
        {
            scalar V = mesh_.cellVolumes()[c];

            scalar rhoMMeanElectron = rhoMMeanElectron_[c]*nParticles(c)/V;
            scalar rhoNMeanElectron = rhoNMeanElectron_[c]*nParticles(c)/V;
            vector UElectron = momentumMeanElectron_[c] /(rhoMMeanElectron*V);
            scalar linearKEMeanElectron =
		            (0.5*linearKEMeanElectron_[c]*nParticles(c))/V;

            electronTemperature_[c] = 2.0/(3.0*physicoChemical::k.value()
                * rhoNMeanElectron)*(linearKEMeanElectron
                - 0.5*rhoMMeanElectron
                * (UElectron & UElectron));
        }

        const DynamicList<dsmcParcel*>& molsInCell = cellOccupancy_[c];

        forAll(molsInCell, mIC)
        {
            dsmcParcel* p = molsInCell[mIC];

            const dsmcParcel::constantProperties& constProp =
                constProps(p->typeId());

            label charge = constProp.charge();

            if (charge == 1)
            {
                const label& cellI = p->cell();

                //- found an ion, add an electron here

                //- electron temperature will be zero if there have been no
                //  electrons in the cell during the simulation

                if (electronTemperature_[cellI] < SMALL)
                {
                    electronTemperature_[cellI] = 6000.0;
                }
                if (electronTemperature_[cellI] > 8.0e4)
                {
                    electronTemperature_[cellI] = 30000.0;
                }

                vector electronVelocity = equipartitionLinearVelocity
                    (
                        electronTemperature_[cellI],
                        constProps_[electronTypeId].mass()
                    );

                if (rhoMMean_[cellI] > VSMALL)
                {
                    cellVelocity_[cellI] = momentumMean_[cellI]
                        /rhoMMean_[cellI];
                }

                labelList vibLevel;

                electronVelocity += cellVelocity_[cellI];

                addNewParcel
                (
                    p->position(),
                    electronVelocity,
                    p->RWF(),
                    0.0,
                    0,
                    p->cell(),
                    p->tetFace(),
                    p->tetPt(),
                    electronTypeId,
                    -1,
                    0,
                    vibLevel
                );
            }
        }
    }
}

Foam::label Foam::dsmcCloud::pickFromCandidateList
(
    DynamicList<label>& candidatesInCell
)
{
    label entry = -1;
    label size = candidatesInCell.size();

    if (size > 0)
    {
        // Choose a random candidate index.
        label randomIndex = randomLabel(0, size-1);
        entry = candidatesInCell[randomIndex];

        // build a new list without the chosen entry
        DynamicList<label> newCandidates(0);

        forAll(candidatesInCell, i)
        {
            if (i != randomIndex)
            {
                newCandidates.append(candidatesInCell[i]);
            }
        }

        // transfer the new list
        candidatesInCell.transfer(newCandidates);
        candidatesInCell.shrink();
    }

    return entry;
}

void Foam::dsmcCloud::updateCandidateSubList
(
    const label& candidate,
    DynamicList<label>& candidatesInSubCell
)
{
    label newIndex = findIndex(candidatesInSubCell, candidate);

    DynamicList<label> newCandidates(0);

    forAll(candidatesInSubCell, i)
    {
        if (i != newIndex)
        {
            newCandidates.append(candidatesInSubCell[i]);
        }
    }

    // transfer the new list
    candidatesInSubCell.transfer(newCandidates);
    candidatesInSubCell.shrink();

//     Info <<  " list (after) " << candidatesInSubCell << endl;
}

Foam::label Foam::dsmcCloud::pickFromCandidateSubList
(
    DynamicList<label>& candidatesInCell,
    DynamicList<label>& candidatesInSubCell
)
{
//     Info << " list (before) " << candidatesInCell << endl;
//     Info << " sub list (before) " << candidatesInSubCell << endl;

    label entry = -1;
    label subCellSize = candidatesInSubCell.size();

    if (subCellSize > 0)
    {
        label randomIndex = randomLabel(0, subCellSize-1);
        entry = candidatesInSubCell[randomIndex];

        DynamicList<label> newSubCellList(0);

        forAll(candidatesInSubCell, i)
        {
            if (i != randomIndex)
            {
                newSubCellList.append(candidatesInSubCell[i]);
            }
        }

        candidatesInSubCell.transfer(newSubCellList);
        candidatesInSubCell.shrink();

//         Info <<  " sub list (after) " << candidatesInSubCell << endl;

        label newIndex = findIndex(candidatesInCell, entry);

        DynamicList<label> newList(0);

        forAll(candidatesInCell, i)
        {
            if (i != newIndex)
            {
                newList.append(candidatesInCell[i]);
            }
        }

        candidatesInCell.transfer(newList);
        candidatesInCell.shrink();

//         Info <<  " list (after) " << candidatesInCell << endl;
    }

    return entry;
}

void Foam::dsmcCloud::collisions()
{
    collisionPartnerSelectionModel_->collide();
}

void Foam::dsmcCloud::addNewParcel
(
    const vector& position,
    const vector& U,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification,
    const labelList& vibLevel
)
{
    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification,
        vibLevel
    );

    // parcels that have been freshly injected on boundary patch faces should
    // be tracked as having crossed that boundary patch face in the time step
    // in which they have been inserted. This is justified as the trackFraction
    // is initialized to a random value in the interval [0, 1] (cf.
    // dsmcParcel::move newParcel handling).
    if (newParcel != -1)
    {
        tracker().trackFaceTransition(typeId, U, RWF, tetFaceI);
    }

    porousMeas().additionInteraction(*pPtr, newParcel);

    addParticle(pPtr);
}

void Foam::dsmcCloud::addNewStuckParcel
(
    const vector& position,
    const vector& U,
    const scalar RWF,
    const scalar ERot,
    const label ELevel,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const label newParcel,
    const label classification,
    const labelList& vibLevel,
    const scalarField& wallTemperature,
    const vectorField& wallVectors
)
{
    dsmcParcel* pPtr = new dsmcParcel
    (
        mesh_,
        position,
        U,
        RWF,
        ERot,
        ELevel,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        newParcel,
        classification,
        vibLevel
    );

    // set parcel stuck
    dsmcParcel& p = *pPtr;
    p.setStuck();
    p.stuck().wallTemperature() = wallTemperature;
    p.stuck().wallVectors() = wallVectors;

    porousMeas().additionInteraction(p, newParcel);

    addParticle(pPtr);
}

Foam::scalar Foam::dsmcCloud::energyRatio
(
    scalar ChiA,
    scalar ChiB
)
{
    scalar ChiAMinusOne = ChiA - 1;

    scalar ChiBMinusOne = ChiB - 1;

    if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
    {
        return rndGen_.sample01<scalar>();
    }

    scalar energyRatio;

    scalar P;

    do
    {
        P = 0;

        energyRatio = rndGen_.sample01<scalar>();

        if (ChiAMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiBMinusOne);
        }
        else if (ChiBMinusOne < SMALL)
        {
            P = pow((1.0 - energyRatio),ChiAMinusOne);
        }
        else
        {
            P =
                pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                    ChiAMinusOne
                )
               *pow
                (
                    (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                    /ChiBMinusOne,
                    ChiBMinusOne
                );
        }
    } while (P < rndGen_.sample01<scalar>());

    return energyRatio;
}

Foam::scalar Foam::dsmcCloud::PSIm
(
    scalar DOFm,
    scalar DOFtot
)
{
    if (DOFm == DOFtot)
    {
        return 1.0;
    }

    if (DOFm == 2.0 && DOFtot == 4.0)
    {
        return rndGen_.sample01<scalar>();
    }

    if (DOFtot < 4.0)
    {
        return (DOFm/DOFtot);
    }

    scalar rPSIm = 0.0;
    scalar prob = 0.0;

    scalar h1 = 0.5*DOFtot - 2.0;
    scalar h2 = 0.5*DOFm - 1.0 + 1.0e-5;
    scalar h3 = 0.5*(DOFtot-DOFm)-1.0 + 1.0e-5;

    do
    {
        rPSIm = rndGen_.sample01<scalar>();
        prob = pow(h1,h1)/(pow(h2,h2)*pow(h3,h3))*pow(rPSIm,h2)*pow(1.0-rPSIm,h3);
    } while (prob < rndGen_.sample01<scalar>());

    return rPSIm;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// for running dsmcFoam+
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const dynamicFvMesh& mesh,
    bool readFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    controlDict_(t.controlDict()),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    dsmcCoordinateSystem_(dsmcCoordinateSystem::New(t, mesh, *this)),
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    rotatingReferenceFrameActive_(false),
    rotatingRegionZoneName_(word::null),
    rotatingFrameOmega_(vector::zero),
    rotatingFrameOrigin_(vector::zero),
    rotatingFrameSubCycles_(1),
    rotatingFrameMaxSubStep_(-1.0),
    rotatingFrameMaxRotationAngle_(Foam::degToRad(5.0)),
    rotatingFrameMeshTolFraction_(0.25),
    rotatingRegionCellMask_(),
    cellOccupancy_(),
    rhoNMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMeanElectron_(mesh_.nCells(), 0.0),
    rhoMMean_(mesh_.nCells(), 0.0),
    momentumMeanElectron_(mesh_.nCells(), vector::zero),
    momentumMean_(mesh_.nCells(), vector::zero),
    linearKEMeanElectron_(mesh_.nCells(), 0.0),
    electronTemperature_(mesh_.nCells(), 0.0),
    cellVelocity_(mesh_.nCells(), vector::zero),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    collisionSelectionRemainder_(),
    constProps_(),
    rndGen_
    (
        particleProperties_.lookupOrDefault<label>
        (
            "seedNumber",
            label(clock::getTime())
        )
      + particleProperties_.lookupOrDefault<label>("seedStride", 7183)
       *Pstream::myProcNo()
    ), 
    controllers_(t, mesh, *this),
    dynamicLoadBalancing_(t, mesh, *this),
    boundaryMeas_(mesh, *this, true),
    fields_(t, mesh, *this),
    boundaries_(t, mesh, *this),
    trackingInfo_(mesh, *this, true),
    binaryCollisionModel_
    (
        BinaryCollisionModel::New
        (
            particleProperties_,
            *this
        )
    ),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh, *this),
    cellMeas_(mesh, *this, true)
{
    if (readFields)
    {
        dsmcParcel::readFields(*this);
    }

    buildConstProps();

    coordSystem().checkCoordinateSystemInputs();
    porousMeas().checkPorousMeasurementsInputs();

    reactions_.initialConfiguration();

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();

    collisionPartnerSelectionModel_ = autoPtr<collisionPartnerSelection>
    (
        collisionPartnerSelection::New(mesh, *this, particleProperties_)
    );

    collisionPartnerSelectionModel_->initialConfiguration();

    fields_.createFields();
    boundaryMeas_.setInitialConfig();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
}

// running dsmcInitialise+
Foam::dsmcCloud::dsmcCloud
(
    Time& t,
    const word& cloudName,
    const dynamicFvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict,
    const bool& clearFields
)
:
    Cloud<dsmcParcel>(mesh, cloudName, false),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    controlDict_(t.controlDict()),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    dsmcCoordinateSystem_(dsmcCoordinateSystem::New(t, mesh, *this)),
    porousMeasurements_(porousMeasurements::New(t, mesh, *this)),
    nTerminalOutputs_
    (
        controlDict_.lookupOrDefault<label>("nTerminalOutputs", 1)
    ),
    rotatingReferenceFrameActive_(false),
    rotatingRegionZoneName_(word::null),
    rotatingFrameOmega_(vector::zero),
    rotatingFrameOrigin_(vector::zero),
    rotatingFrameSubCycles_(1),
    rotatingFrameMaxSubStep_(-1.0),
    rotatingFrameMaxRotationAngle_(Foam::degToRad(5.0)),
    rotatingFrameMeshTolFraction_(0.25),
    rotatingRegionCellMask_(),
    cellOccupancy_(),
    rhoNMeanElectron_(),
    rhoMMeanElectron_(),
    rhoMMean_(),
    momentumMeanElectron_(),
    momentumMean_(),
    linearKEMeanElectron_(),
    electronTemperature_(),
    cellVelocity_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 3, -1, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    collisionSelectionRemainder_(),
    constProps_(),
    rndGen_
    (
        particleProperties_.lookupOrDefault<label>
        (
            "seedNumber",
            label(clock::getTime())
        )
      + particleProperties_.lookupOrDefault<label>("seedStride", 1526)
       *Pstream::myProcNo()
    ),
    controllers_(t, mesh),
    dynamicLoadBalancing_(t, mesh, *this),
    boundaryMeas_(mesh, *this),
    fields_(t, mesh),
    boundaries_(t, mesh),
    trackingInfo_(mesh, *this),
    binaryCollisionModel_(),
    collisionPartnerSelectionModel_(),
    reactions_(t, mesh),
    cellMeas_(mesh, *this)
{
    if (!clearFields)
    {
        dsmcParcel::readFields(*this);
    }

    label initialParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(initialParcels, sumOp<label>());
    }

    if (clearFields)
    {
        Info << "clearing existing field of parcels " << endl;

        clear();

        initialParcels = 0;
    }

    buildConstProps();

    coordSystem().checkCoordinateSystemInputs(true);

    dsmcAllConfigurations conf(dsmcInitialiseDict, *this);
    conf.setInitialConfig();

    label finalParcels = this->size();

    if (Pstream::parRun())
    {
        reduce(finalParcels, sumOp<label>());
    }

    Info << nl << "Initial no. of parcels: " << initialParcels
         << " added parcels: " << finalParcels - initialParcels
         << ", total no. of parcels: " << finalParcels
         << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dsmcCloud::~dsmcCloud()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dsmcCloud::evolve()
{
    evolve_moveAndCollide();
    evolve_fields();
}

void Foam::dsmcCloud::evolve_moveAndCollide()
{
    boundaries_.updateTimeInfo();
    fields_.updateTimeInfo();
    controllers_.updateTimeInfo();

    dsmcParcel::trackingData td(*this);

    if (debug)
    {
        this->dumpParticlePositions();
    }

    controllers_.controlBeforeMove();
    boundaries_.controlBeforeMove();

    //- Remove electrons
    if (findIndex(typeIdList_, "e-") != -1)
    {
        // TODO VINCENT: there is a clever way than rebuilding entire cell occ.
        removeElectrons();
        buildCellOccupancy();
    }

    // Boundary measurements use the parcel's event-time weight, not the destination cell weight.

    Cloud<dsmcParcel>::move(td, deltaTValue());

    // Rebuild cell occupancy after motion.
    buildCellOccupancy();

    //- Add electrons back after the move function
    if (findIndex(typeIdList_, "e-") != -1)
    {
        // TODO VINCENT: there is a clever way than rebuilding entire cell occ.
        addElectrons();
        buildCellOccupancy();
    }

    // Apply the selected coordinate-system weighting or Cartesian regulation.
    coordSystem().evolve();

    controllers_.controlBeforeCollisions();
    boundaries_.controlBeforeCollisions();

    // Advance molecular collisions.
    collisions();

    //- Reactions may have changed cell occupancy, update if any reaction
    if (reactions_.nReactions() != 0)
    {
        buildCellOccupancy();
    }

    // Apply post-collision regulation according to its own configured interval.
    coordSystem().evolveAfterCollisions();

    controllers_.controlAfterCollisions();
    boundaries_.controlAfterCollisions();
}

void Foam::dsmcCloud::evolve_fields()
{
    reactions_.outputData();

    fields_.calculateFields();

    fields_.writeFields();

    controllers_.calculateProps();
    controllers_.outputResults();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    boundaryMeas_.outputResults();

    trackingInfo_.clean();
    boundaryMeas_.clean();
    cellMeas_.clean();
}

Foam::label Foam::dsmcCloud::nTerminalOutputs()
{
    return nTerminalOutputs_;
}

void Foam::dsmcCloud::info()
{
    // Parcel count and reference-equivalent molecule count are distinct for unequal weights.
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    const scalarList& iM = infoMeasurements();

    scalar nMol = iM[6];
    reduce(nMol, sumOp<scalar>());

    scalar linearKineticEnergy = iM[1];
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar rotationalEnergy = iM[2];
    reduce(rotationalEnergy, sumOp<scalar>());

    scalar vibrationalEnergy = iM[3];
    reduce(vibrationalEnergy, sumOp<scalar>());

    scalar electronicEnergy = iM[4];
    reduce(electronicEnergy, sumOp<scalar>());

    scalar stuckMolecules = iM[5];
    reduce(stuckMolecules, sumOp<scalar>());

    Info<< "    Number of DSMC particles        = "
        << nDsmcParticles
        << endl;

    if (nDsmcParticles > VSMALL)
    {
        // Reference-equivalent counts divide represented molecules by nParticle().
       Info << "    Reference-equivalent stuck particles = "
            << stuckMolecules/nParticle() << nl
            << "    Reference-equivalent free particles  = "
            << nMol/nParticle() << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average rotational energy       = "
            << rotationalEnergy/nMol << nl
            << "    Average vibrational energy      = "
            << vibrationalEnergy/nMol << nl
            << "    Average electronic energy       = "
            << electronicEnergy/nMol << nl
            << "    Total energy                    = "
            << (linearKineticEnergy + rotationalEnergy
                + vibrationalEnergy + electronicEnergy)
            << endl;

        porousMeas().writePorousMeasurementsInfo();
    }
}

void Foam::dsmcCloud::loadBalanceCheck()
{
    dynamicLoadBalancing_.update();
}

void Foam::dsmcCloud::loadBalance(const label noRefinement)
{
    dynamicLoadBalancing_.perform(noRefinement);
}

void Foam::dsmcCloud::autoMap(const mapPolyMesh& mapper)
{
    dsmcParcel::trackingData td(*this);

    Cloud<dsmcParcel>::autoMap(td, mapper);

    coordSystem().dtModel().update();

    buildCellOccupancyFromScratch();
    buildCollisionSelectionRemainderFromScratch();
    resetBoundaries();
    resetMeasurementTools();
}

Foam::label Foam::dsmcCloud::randomLabel
(
    const label valOne,
    const label valTwo
)
{
    if (valOne == valTwo)
    {
        return valOne;
    }
    else
    {
        const label start = Foam::min(valOne, valTwo);
        const label end = Foam::max(valOne, valTwo);

        label val = start + label(rndGen_.sample01<scalar>()*(end - start + 1));

        // Rare case when scalar01() returns exactly 1.000 and the truncated
        // value would be out of range.
        if(val == end + 1)
        {
            val = randomLabel(start, end);
        }
        return val;
    }
}

Foam::vector Foam::dsmcCloud::equipartitionLinearVelocity
(
    const scalar temperature,
    const scalar mass
)
{
    return sqrt(physicoChemical::k.value()*temperature/mass)
        *rndGen_.GaussNormal<vector>();
}

Foam::vector Foam::dsmcCloud::chapmanEnskogVelocity
(
    const scalar temperature,
    const scalar mass,
    const vector& q,
    const tensor& tau
)
{
    const scalar B = max(mag(q), mag(tau));
    const scalar A = 1.0 + 30.0*B;

    bool repeatTry = true;

    vector CTry = vector::zero;

    while (repeatTry)
    {
        CTry = rndGen_.GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + (q & CTry)*(0.4*(CTry & CTry) - 1.0)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen_.sample01<scalar>())
        {
            repeatTry = false;
        }
    }

    return CTry*sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

void Foam::dsmcCloud::generalisedChapmanEnskog
(
    const label& typeID,
    const scalar& translationalTemperature,
    const scalar& rotationalTemperature,
    const scalar& vibrationalTemperature,
    const scalar& mass,
    const vector& D,
    const vector& qTra,
    const vector& qRot,
    const vector& qVib,
    const tensor& tau,
    scalar& ERot,
    labelList& vibLevel,
    vector& U
)
{
    const scalar kB = physicoChemical::k.value();
    const label nVibModes = constProps(typeID).nVibrationalModes();

    scalar B = max(mag(D), mag(tau));
    B = max(B, mag(qTra));
    B = max(B, mag(qRot));
    B = max(B, mag(qVib));
    const scalar A = 1.0 + 30.0*B;

    const scalarList& thetaV = constProps(typeID).thetaV();
    const scalar epsRotAv = rotationalTemperature/translationalTemperature;

    scalar epsVibAv = 0.0;
    if (nVibModes > 0 && vibrationalTemperature > 5.)
    {
        forAll(vibLevel, mode)
        {
            epsVibAv += thetaV[mode]/vibrationalTemperature
                /(exp(thetaV[mode]/vibrationalTemperature) - 1.0);
        }

    }

    vector CTry = vector::zero;

    scalar epsRot = 0.0;
    scalar epsVib = 0.0;

    bool repeatTry = true;

    while (repeatTry)
    {
        ERot =
            equipartitionRotationalEnergy
            (
                rotationalTemperature,
                constProps(typeID).rotationalDegreesOfFreedom()
            );

        epsRot = ERot/(kB*translationalTemperature);

        if (nVibModes > 0 && vibrationalTemperature > 5.)
        {
            vibLevel =
                equipartitionVibrationalEnergyLevel
                (
                    vibrationalTemperature,
                    nVibModes,
                    typeID
                );

            scalar epsVib = 0.0;
            forAll(vibLevel, mode)
            {
                epsVib += vibLevel[mode]*thetaV[mode];
            }

            epsVib /= vibrationalTemperature;
        }

        CTry = rndGen_.GaussNormal<vector>()/sqrt(2.0);

        const scalar gammaTry = 1.0 + 2.0*(D & CTry)
            + (qTra & CTry)*(0.4*(CTry & CTry) - 1.0)
            + (qRot & CTry)*(epsRot - epsRotAv)
            + (qVib & CTry)*(epsVib - epsVibAv)
            - (CTry & (tau & CTry));

        if (gammaTry >= A*rndGen_.sample01<scalar>())
        {
            repeatTry = false;
        }
    }

    U = CTry*sqrt(2.0*kB*translationalTemperature/mass);
}

Foam::scalar Foam::dsmcCloud::equipartitionRotationalEnergy
(
    const scalar temperature,
    const scalar rotationalDof
)
{
    scalar ERot = 0.0;

    if (rotationalDof < SMALL)
    {
        return ERot;
    }
    else if (rotationalDof < 2.0 + SMALL && rotationalDof > 2.0 - SMALL)
    {
        // Special case for rDof = 2, i.e. diatomics;
        ERot = -log(rndGen_.sample01<scalar>())*physicoChemical::k.value()*temperature;
    }
    else
    {
        scalar a = 0.5*rotationalDof - 1;

        scalar energyRatio;

        scalar P = -1;

        do
        {
            energyRatio = 10*rndGen_.sample01<scalar>();

            P = pow((energyRatio/a), a)*exp(a - energyRatio);

        } while (P < rndGen_.sample01<scalar>());

        ERot = energyRatio*physicoChemical::k.value()*temperature;
    }

    return ERot;
}

Foam::labelList Foam::dsmcCloud::equipartitionVibrationalEnergyLevel
(
    const scalar temperature,
    const label nVibrationalModes,
    const label typeId
)
{
    labelList vibLevel(nVibrationalModes, 0);

    if (nVibrationalModes == 0)
    {
        return vibLevel;
    }
    else
    {
        forAll(vibLevel, mode)
        {
            vibLevel[mode] = -log(rndGen_.sample01<scalar>())*temperature
                /constProps(typeId).thetaV_m(mode);
        }
    }

    return vibLevel;
}

Foam::label Foam::dsmcCloud::equipartitionElectronicLevel
(
    const scalar temperature,
    const labelList& electronicDegeneracyList,
    const scalarList& electronicEnergyList
)
{
    const scalar EMax = physicoChemical::k.value()*temperature;
    const label jMax = electronicDegeneracyList.size() - 1;

    //- Random integer between 0 and jMax
    label jDash = 0;
    //- Maximum possible electronic energy level within list based on k*TElec
    scalar EJ = 0.0;
    //- Maximum possible degeneracy level within list
    label gJ = 0;
    //- Selected intermediate integer electronic level (0 to jMax)
    label jSelect = 0;
    //- Maximum denominator value in Liechty pdf (see below)
    scalar expMax = 0.0;
    //- Summation term based on random electronic level
    scalar expSum = 0.0;
    //- Boltzmann distribution of Eq. 3.1.1 of Liechty thesis
    scalar boltz = 0.0;
    //- Distribution function Eq. 3.1.2 of Liechty thesis
    scalar func = 0.0;

    if (jMax > 0 and temperature > SMALL)
    {
        //- Calculate summation term in denominator of Eq. 3.1.1 in Liechty
        //  thesis
        forAll(electronicDegeneracyList, i)
        {
            expSum += electronicDegeneracyList[i]
                *exp(-electronicEnergyList[i]/EMax);
        }

        //- Select maximum integer energy level based on boltz value.
        //  Note that this depends on the temperature.
        scalar boltzMax = 0.0;

        forAll(electronicDegeneracyList, i)
        {
            //- Eq. 3.1.1 of Liechty thesis.
            boltz =
                electronicDegeneracyList[i]
               *exp(-electronicEnergyList[i]/EMax)
               /expSum;

            if (boltzMax < boltz)
            {
                boltzMax = boltz;
                jSelect = i;
            }
        }

        //- Max. poss energy in list: list goes from 0 to jMax
        EJ = electronicEnergyList[jSelect];
        //- Max. poss degeneracy in list: list goes from 0 to jMax
        gJ = electronicDegeneracyList[jSelect];
        //- Max. in denominator of Liechty pdf for initialisation/wall
        //  bcs/freestream EEle etc..
        expMax = gJ*exp(-EJ/EMax);

        //- Acceptance - rejection based on Eq. 3.1.2 of Liechty thesis
        do
        {
            jDash = randomLabel(0, jMax);
            func =
                electronicDegeneracyList[jDash]
               *exp(-electronicEnergyList[jDash]/EMax)
               /expMax;
        } while(func < rndGen_.sample01<scalar>());
    }

    return jDash;
}

Foam::scalar Foam::dsmcCloud::postCollisionRotationalEnergy
(
    const scalar rotationalDof,
    const scalar ChiB
)
{
    scalar energyRatio = 0.0;

    if (rotationalDof == 2.0)
    {
        energyRatio = 1.0 - pow(rndGen_.sample01<scalar>(), 1.0/ChiB);
    }
    else
    {
        const scalar ChiA = 0.5*rotationalDof;

        scalar ChiAMinusOne = ChiA - 1.;

        scalar ChiBMinusOne = ChiB - 1.;

        if (ChiAMinusOne < SMALL && ChiBMinusOne < SMALL)
        {
            return rndGen_.sample01<scalar>();
        }

        scalar P = 0.0;

        do
        {
            P = 0;

            energyRatio = rndGen_.sample01<scalar>();

            if (ChiAMinusOne < SMALL)
            {
                P = pow(1.0 - energyRatio, ChiBMinusOne);
            }
            else if (ChiBMinusOne < SMALL)
            {
                P = pow(1.0 - energyRatio, ChiAMinusOne);
            }
            else
            {
                P =
                    pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*energyRatio/ChiAMinusOne,
                        ChiAMinusOne
                    )
                *pow
                    (
                        (ChiAMinusOne + ChiBMinusOne)*(1 - energyRatio)
                        /ChiBMinusOne,
                        ChiBMinusOne
                    );
            }
        } while (P < rndGen_.sample01<scalar>());
    }

    return energyRatio;
}

Foam::label Foam::dsmcCloud::postCollisionVibrationalEnergyLevel
(
    bool postReaction,
    const label vibLevel,
    const label iMax,
    const scalar thetaV,
    const scalar thetaD,
    const scalar refTempZv,
    const scalar omega,
    const scalar Zref,
    const scalar Ec,
    const scalar fixedZv,
    const label invZvFormulation,
    const label cellI
)
{
    label iDash = vibLevel;

    if (postReaction)
    {
        // post-collision quantum number
        scalar func = 0.0;
        scalar EVib = 0.0;

        do // acceptance - rejection
        {
            iDash = randomLabel(0, iMax);
            EVib = iDash*physicoChemical::k.value()*thetaV;

            // - equation 5.61, Bird
            func = pow(1.0 - EVib/Ec, 1.5 - omega);

        } while(func < rndGen_.sample01<scalar>());
    }
    else
    {
        scalar inverseVibrationalCollisionNumber = 1.0;

        if (fixedZv == 0)
        {
            //- Temperature used to calculate Zv
            scalar T = 0;

            if (invZvFormulation == 0)
            {
                //- "Quantised collision temperature" (equation 3, Bird 2010)
                //  denominator from Bird 5.42
                T = iMax*thetaV/(3.5 - omega);
            }
            else if (invZvFormulation == 1)
            {
                //- Macroscopic (overall) temperature
                const scalar TMacro = fields().overallT(cellI);

                if (TMacro > SMALL)
                {
                    T = TMacro;
                }
                else
                {
                    //- Collision temperature used instead
                    //  the pre-2008 formulation is recovered
                    T = iMax*thetaV/(3.5 - omega);
                }

            }
            else
            {
                //- Macroscopic (translational) temperature
                
 //TODO
                // Collision temperature for the time being
                // it gives a better agreement
                T = iMax*thetaV/(3.5 - omega);
            }

            const scalar pow1 = pow(thetaD/T, 1./3.) - 1.0;

            const scalar pow2 = pow(thetaD/refTempZv, 1./3.) - 1.0;

            //- vibrational collision number (equation 2, Bird 2010)
            const scalar ZvP1 = pow(thetaD/T, omega);

            const scalar ZvP2 =
                pow
                (
                    Zref*pow(thetaD/refTempZv, -omega),
                    pow1/pow2
                );

            const scalar Zv = ZvP1*ZvP2;

            //- In order to obtain the relaxation rate corresponding to Zv with the collision
            //  energy-based procedure, the inelastic fraction should be set to about 1/(5Zv)
            //  Bird 2008 RGD "A Comparison of Collision Energy-Based and Temperature-Based..."
            if (invZvFormulation == 2)
            {
                inverseVibrationalCollisionNumber = 1.0/(5.0*Zv);
            }
            else
            {
                inverseVibrationalCollisionNumber = 1.0/Zv;
            }
        }
        else
        {
            inverseVibrationalCollisionNumber = 1.0/fixedZv;
        }

        if (inverseVibrationalCollisionNumber > rndGen_.sample01<scalar>())
        {
            // post-collision quantum number
            scalar func = 0.0;
            scalar EVib = 0.0;

            do // acceptance - rejection
            {
                iDash = randomLabel(0, iMax);

                EVib = iDash*physicoChemical::k.value()*thetaV;

                // - equation 5.61, Bird
                func = pow(1.0 - EVib/Ec, 1.5 - omega);

            } while(func < rndGen_.sample01<scalar>());
        }
    }

    return iDash;
}

Foam::label Foam::dsmcCloud::postCollisionElectronicEnergyLevel
(
    const scalar Ec,
    const label jMax,
    const scalar omega,
    const scalarList& EElist,
    const labelList& gList
)
{

    //- Maximum allowable electronic level obtainable from Ecoll
    label jSelectA = 0;
    //- Energy level maximazing expression gList[j]*pow(Ec - EElist[j], 1.5 - omega)
    label jSelectB = 0;

    scalar g = 0.0;
    scalar gMax = 0.0;

    //- Determine the maximum possible integer energy level immediately below Ec
    forAll(gList, i)
    {
       if (EElist[i] > Ec)
       {
           break;
       }

       jSelectA = i;

       //- Eq. 3.1.6 of Liechty thesis.
       g = gList[i]*pow(Ec - EElist[i], 1.5 - omega);

       if (gMax < g)
       {
           gMax = g;
           jSelectB = i;
       }
    }

    //- Liechty's procedure - step 3
    //  Minimum of jSelectA and jSelectB
    const label jSelect = min(jSelectA, jSelectB);

    //- Max. poss energy in list: list goes from 0 to jSelect
    const scalar EJ = EElist[jSelect];
    //- Max. poss degeneracy in list: list goes from 0 to jSelect
    const label gJ = gList[jSelect];
    //- Max. denominator of Liechty pdf for post-collision pdf
    const scalar denomMax = gJ*pow(Ec - EJ, 1.5 - omega);

    //- Acceptance - rejection based on Eq. 3.1.8 of Liechty thesis
    //- Post-collision electronic energy
    label jDash = 0;
    scalar prob = 0.0;

    do
    {
       jDash = randomLabel(0, jSelectA);
       prob = gList[jDash]*pow(Ec - EElist[jDash], 1.5 - omega)/denomMax;

    } while(prob < rndGen_.sample01<scalar>());

    return jDash;

}

void Foam::dsmcCloud::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    forAllConstIter(dsmcCloud, *this, iter)
    {
        const dsmcParcel& p = iter();

        pObj<< "v " << p.position().x()
            << " "  << p.position().y()
            << " "  << p.position().z()
            << nl;
    }

    pObj.flush();
}

void Foam::dsmcCloud::reBuildCellOccupancy()
{
    buildCellOccupancy();
}

Foam::dsmcParcel* Foam::dsmcCloud::cloneFreeParcelInCell
(
    const dsmcParcel& parentParcel,
    const scalar newRWF
)
{
    // Local cloning is restricted to free parcels, including uncollided residuals.
    if (!parentParcel.isFree())
    {
        FatalErrorInFunction
            << "cloneFreeParcelInCell was called for a non-free parcel."
            << abort(FatalError);
    }

    // 这里直接使用 dsmcParcel 的拷贝构造，而不是手工逐成员重建。
    // 这样可以完整保留父粒子的底层状态，包括：
    // 1. particle 基类中的 face/stepFraction/origProc/origId；
    // 2. DSMC parcel 的速度、内部态、tracking/stuck 元信息。
    // 随后我们只修改当前 clone 的 RWF，使其成为父粒子真实分子集合中的一个
    // 等比例子集。
    dsmcParcel* clonePtr = new dsmcParcel(parentParcel);

    dsmcParcel& clone = *clonePtr;
    clone.RWF() = newRWF;

    addParticle(clonePtr);
    insertParcelInCellOccupancy(clonePtr);

    return clonePtr;
}

Foam::label Foam::dsmcCloud::equaliseCellParcelWeightsForCollisions
(
    const label cellI,
    const scalar integerRatioTolerance,
    const label maxMultiplicityPerParcel,
    const label maxAutoEqualiseDenominator
)
{
    DynamicList<dsmcParcel*>& cellParcels = cellOccupancy_[cellI];

    if (cellParcels.size() <= 1)
    {
        return 0;
    }

    if (integerRatioTolerance <= SMALL)
    {
        FatalErrorInFunction
            << "integerRatioTolerance must be greater than SMALL. Got "
            << integerRatioTolerance
            << abort(FatalError);
    }

    if (maxMultiplicityPerParcel < 1)
    {
        FatalErrorInFunction
            << "maxMultiplicityPerParcel must be >= 1. Got "
            << maxMultiplicityPerParcel
            << abort(FatalError);
    }

    if (maxAutoEqualiseDenominator < 1)
    {
        FatalErrorInFunction
            << "maxAutoEqualiseDenominator must be >= 1. Got "
            << maxAutoEqualiseDenominator
            << abort(FatalError);
    }

    // 先锁定当前 cell 内最小的真实代表分子数。
    // 在同一 cell 内，representedMolecules(parcel) 与 RWF(parcel) 线性等价，
    // 因为 dtModel(cellI) 因子对该 cell 的所有 parcel 都相同。
    scalar minRepresentedMolecules = GREAT;

    forAll(cellParcels, parcelI)
    {
        const scalar represented =
            max(representedMolecules(*cellParcels[parcelI]), SMALL);

        minRepresentedMolecules = min(minRepresentedMolecules, represented);
    }

    if (minRepresentedMolecules <= SMALL)
    {
        FatalErrorInFunction
            << "Encountered non-positive minimum represented molecules in cell "
            << cellI << "." << abort(FatalError);
    }

    const label originalSize = cellParcels.size();
    label clonesAdded = 0;
    label chosenDenominator = -1;
    labelList chosenMultiplicities(originalSize, 1);

    // 这里从最小分母开始搜索，优先选择“生成 clone 最少”的共同基准权重。
    // 当 denominator = 1 时，算法就退化成上一轮已经验证通过的“最小权重
    // 整数倍”模式；只有它失败时，才继续尝试 denominator = 2, 3, ...
    for
    (
        label denominator = 1;
        denominator <= maxAutoEqualiseDenominator;
        ++denominator
    )
    {
        bool denominatorAccepted = true;
        labelList trialMultiplicities(originalSize, 1);

        for (label parcelI = 0; parcelI < originalSize; ++parcelI)
        {
            const dsmcParcel& parcel = *cellParcels[parcelI];

            const scalar represented =
                max(representedMolecules(parcel), SMALL);

            // 对应的共同基准权重为：
            // $$
            // N_{\mathrm{base}} = \frac{N_{\min}}{d}
            // $$
            // 因而每个 parcel 需要满足：
            // $$
            // m_p \approx \frac{N_{\mathrm{real}}(p)}{N_{\mathrm{base}}}
            //      = d\frac{N_{\mathrm{real}}(p)}{N_{\min}}
            // $$
            const scalar multiplicityExact =
                scalar(denominator)*represented/minRepresentedMolecules;

            const label multiplicity =
                max(label(::floor(multiplicityExact + 0.5)), 1);

            const scalar relativeRatioError =
                mag(multiplicityExact - scalar(multiplicity))
               /max(multiplicityExact, SMALL);

            if
            (
                relativeRatioError > integerRatioTolerance
             || multiplicity > maxMultiplicityPerParcel
            )
            {
                denominatorAccepted = false;
                break;
            }

            trialMultiplicities[parcelI] = multiplicity;
        }

        if (denominatorAccepted)
        {
            chosenDenominator = denominator;
            chosenMultiplicities.transfer(trialMultiplicities);
            break;
        }
    }

    if (chosenDenominator == -1)
    {
        FatalErrorInFunction
            << "Cell " << cellI
            << " could not be collision-time auto-equalised with any common "
            << "base weight up to denominator "
            << maxAutoEqualiseDenominator << "." << nl
            << "minimum represented molecules = " << minRepresentedMolecules << nl
            << "allowed ratio tolerance = " << integerRatioTolerance << nl
            << "max multiplicity per parcel = " << maxMultiplicityPerParcel << nl
            << "This means the current cell is not representable as a finite "
            << "set of equal-weight same-state parcels under the current "
            << "weighted-collision model."
            << abort(FatalError);
    }

    // 这里只遍历进入函数时已有的 parcel。
    // 新增 clone 会被 append 到 cellParcels 尾部，但它们本身已经具有目标
    // 基准权重，因此不需要在同一轮再次参与拆分。
    for (label parcelI = 0; parcelI < originalSize; ++parcelI)
    {
        dsmcParcel& parcel = *cellParcels[parcelI];
        const label multiplicity = chosenMultiplicities[parcelI];

        if (multiplicity <= 1)
        {
            continue;
        }

        // 由于 representedMolecules 与 RWF 在线性成正比，
        // 一旦 multiplicity 已确定，只需把父粒子的 RWF 平均分到每个子粒子：
        // $$
        // RWF_{child} = \frac{RWF_{parent}}{m_p}
        // $$
        // 这样每个子粒子的真实代表分子数都会回到统一的基准权重。
        const scalar childRWF = parcel.RWF()/scalar(multiplicity);

        parcel.RWF() = childRWF;

        for (label cloneI = 1; cloneI < multiplicity; ++cloneI)
        {
            cloneFreeParcelInCell(parcel, childRWF);
            clonesAdded++;
        }
    }

    return clonesAdded;
}

void Foam::dsmcCloud::insertParcelInCellOccupancy(dsmcParcel* p)
{
    // Append without shrinking so repeated residual insertion can reuse list capacity.
    cellOccupancy_[p->cell()].append(p);
}

bool Foam::dsmcCloud::performDeferredParallelTransfers()
{
    // DSMC 当前唯一需要“延迟到 Cloud::move() 并行 transfer while-loop 内
    // 再统一处理”的粒子转移，来自某些 cyclic boundary model 的自定义
    // 远端注入路径。
    //
    // 这里把调度权继续交给 boundaries_，而不是把 stage-interface 特例
    // 直接写死在 dsmcCloud 里。这样：
    // 1. Cloud 层只负责提供时序钩子；
    // 2. 具体哪一个边界模型需要远端收发，由边界对象自己决定；
    // 3. 后续若别的 DSMC 周期边界也需要类似机制，可以直接复用。
    return boundaries_.controlAfterMove();
}

void Foam::dsmcCloud::removeParcelFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<dsmcParcel*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if (c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    molsInCell.shrink();
    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}

// NEW DANIEL *****************************************************************

void Foam::dsmcCloud::resetHybridTraRotVib
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    wordList& zonesToReset
)
{
    //Info << "Deleting (" << typeOfReset << ")" << endl;
    //scalar time1_ = mesh_.time().elapsedCpuTime();
    //molsToDelete(mesh_, *this, typeOfReset);
    //scalar time2_ = mesh_.time().elapsedCpuTime();
    //Pout << "Proc " << UPstream::myProcNo() << " deletion time: " << time2_ - time1_ << "s (" << time1_ << ", " << time2_ << ")" << nl << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();
    forAll(zonesToReset, zoneToResetI)
    {
        List<scalar> massToIntroduce(TvInitial.size(), 0.0);//////////////
        List<scalar> massIntroduced(TvInitial.size(), 0.0);///////////////
        word regionName(zonesToReset[zoneToResetI]);
        label zoneId = cellZones.findZoneID(regionName);

        if (zoneId == -1)
        {
            FatalErrorIn("resetHybrid")
                << "Cannot find region: " << regionName << nl << "in: "
                << mesh_.time().constant()/"cellZones"
                << exit(FatalError);
        }

        const cellZone& zone = cellZones[zoneId];

        if (zone.size())
        {
            Info << "\nInserting particles in" << endl;

            forAll(zone, c)
            {
                const label cellI = zone[c];
                forAll(typeIdList_, i)  ////////////////////////////////
                {                       ////////////////////////////////
                    massToIntroduce[i] += numberDensitiesField[i][cellI]*mesh_.V()[cellI];
                }                       ////////////////////////////////
                List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
                (
                    mesh_,
                    cellI
                );

                forAll(cellTets, tetI)
                {
                    const tetIndices& cellTetIs = cellTets[tetI];

                    tetPointRef tet = cellTetIs.tet(mesh_);

                    const scalar tetVolume = tet.mag();

                    forAll(typeIdList_, i)
                    {
                        scalar numberDensity = numberDensitiesField[i][cellI];
                        scalar translationalTemperature = TtInitial[i][cellI];
                        scalar rotationalTemperature = TrInitial[i][cellI];
                        scalar vibrationalTemperature = TvInitial[i][cellI];
                        vector velocity = UInitial[cellI];

                        // Calculate the number of particles required
                        scalar particlesRequired = numberDensity*tetVolume;

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

                        massIntroduced[i] += nParticlesToInsert;///////////////

                        for (label pI = 0; pI < nParticlesToInsert; pI++)
                        {
                            point p = tet.randomPoint(rndGen_);

                            vector U = vector::zero;

                            scalar ERot = 0.0;

                            labelList vibLevel
                            (
                                constProps(i).thetaV().size(),
                                0
                            );

                            label ELevel = 0; // TODO by generalisedChapmanEnskog

                            generalisedChapmanEnskog
                            (
                                i,
                                translationalTemperature,
                                rotationalTemperature,
                                vibrationalTemperature,
                                constProps(i).mass(),
                                DInitial[i][cellI],
                                qtInitial[i][cellI],
                                qrInitial[i][cellI],
                                qvInitial[i][cellI],
                                tauInitial[i][cellI],
                                ERot,
                                vibLevel,
                                U
                            );

                            U += velocity;

                            label newParcel = -1;

                            label classification = 0;

                            const scalar RWF = coordSystem().RWF(cellI);

                            addNewParcel
                            (
                                p,
                                U,
                                RWF,
                                ERot,
                                ELevel,
                                cellI,
                                cellTetIs.face(),
                                cellTetIs.tetPt(),
                                i,
                                newParcel,
                                classification,
                                vibLevel
                            );
                        }
                    }
                }
            }
        }

        Info<< "      Zone: " << regionName << endl;

        forAll(typeIdList_, i)
        {
            const scalar mass = this->constProps(i).mass();
            massToIntroduce[i] *= mass * nParticle();
            massIntroduced[i] *= mass * nParticle();
            Info<< "        Specie " << typeIdList_[i]
                << ", mTI: "
                << massToIntroduce[i] << "; mI: " << massIntroduced[i]
                << " (" << 100.0 * (massIntroduced[i]
                / (massToIntroduce[i] + VSMALL) - 1.0) << "% off)" << endl;
        }

        Info<< endl;
    }

    buildCellOccupancy();
}

void Foam::dsmcCloud::resetHybridWhenUpdated
(
    volVectorField& UInitial,
    PtrList<volScalarField>& TtInitial,
    PtrList<volScalarField>& TrInitial,
    PtrList<volScalarField>& TvInitial,
    PtrList<volScalarField>& numberDensitiesField,
    PtrList<volVectorField>& DInitial,
    PtrList<volVectorField>& qtInitial,
    PtrList<volVectorField>& qrInitial,
    PtrList<volVectorField>& qvInitial,
    PtrList<volTensorField>& tauInitial,
    word& typeOfReset,
    word& zoneToReset
)
{

}

void Foam::dsmcCloud::shockReset()
{
    label nDsmcParticles = this->size();
    reduce(nDsmcParticles, sumOp<label>());

    const IOdictionary& shockDict
    (
        IOobject
        (
            "shockDict",
            mesh_.time().system(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar initialParticles = readScalar(shockDict.lookup("initialParticles"));
    scalar maxPercentageParticles = readScalar(shockDict.lookup("maxPercentageParticles"));
    scalar rhoLeft = readScalar(shockDict.lookup("rhoLeft"));
    scalar rhoRight = readScalar(shockDict.lookup("rhoRight"));
    scalar xLeft = readScalar(shockDict.lookup("xLeft"));
    scalar xRight = readScalar(shockDict.lookup("xRight"));
    scalar areaTube = readScalar(shockDict.lookup("areaTube"));

    scalar maxDeltaParticles = maxPercentageParticles * initialParticles / 100.0;

    scalar deltaParticles = nDsmcParticles - initialParticles;
    Info << "deltaParticles [%]: " << 100.0 * deltaParticles / initialParticles
        << nl << "updating:" << (mag(deltaParticles) >= maxDeltaParticles)
        << endl;

    scalar deltaXTonParticles = deltaParticles * this->constProps(0).mass()
        / areaTube / (rhoRight - rhoLeft);
    Info << "deltaX/nEqParticles: " << deltaXTonParticles << endl;

    if (mag(deltaParticles) >= maxDeltaParticles)
    {
        forAll(mesh_.cells(), cellI)
        {
            const List<dsmcParcel*>& parcelsInCell
                = cellOccupancy_[cellI];

            scalar deltaX = deltaXTonParticles*nParticles(cellI); // NEW VINCENT

            forAll(parcelsInCell, pIC)
            {
                dsmcParcel* p = parcelsInCell[pIC];

                // deltaX > 0: if position > xRight + deltaX => out of domain
                // deltaX < 0: if position < xLeft + deltaX => out of domain
                if
                (
                    (p->classification() < 1000)
                    &&
                    (
                        (p->position().x() < xLeft + deltaX)
                            || (p->position().x() > xRight + deltaX)
                    )
                )
                {
                    p->classification() += 2000;
                    dsmcParcel* pPtr = new dsmcParcel(*p);
                    addParticle(pPtr);
                    p->classification() -= 2000;
                }
            }
        }

        forAllIter(dsmcCloud, *this, iter)
        {
            dsmcParcel& p = iter();

            if (p.classification() > 1000)
            {
                p.classification() -= 2000;
            }
            else
            {
                p.position().x() += deltaXTonParticles*nParticles(p.cell()); // NEW VINCENT

                if ((p.position().x() > xRight) || (p.position().x() < xLeft))
                {
                    deleteParticle(p);
                }
            }
        }

        buildCellOccupancy();
    }
}
// END NEW DANIEL *************************************************************

bool Foam::dsmcCloud::read()
{
    if (regIOobject::read())
    {
        initialiseRotatingReferenceFrame();
        return true;
    }
    else
    {
        return false;
    }
}

// ************************************************************************* //
