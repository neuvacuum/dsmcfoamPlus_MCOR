/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dsmcTranslationalStageInterfacePatch.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(dsmcTranslationalStageInterfacePatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcTranslationalStageInterfacePatch,
    dictionary
);


namespace
{
    void patchXBounds
    (
        const polyMesh& mesh,
        const polyPatch& patch,
        scalar& minX,
        scalar& maxX
    )
    {
        minX = GREAT;
        maxX = -GREAT;

        for(label i = 0; i < patch.size(); ++i)
        {
            const label globalFaceI = patch.start() + i;
            const face& f = mesh.faces()[globalFaceI];

            forAll(f, fp)
            {
                const scalar x = mesh.points()[f[fp]].x();
                minX = min(minX, x);
                maxX = max(maxX, x);
            }
        }

        if (minX > maxX)
        {
            FatalErrorInFunction
                << "Unable to calculate x bounds for patch " << patch.name()
                << exit(FatalError);
        }
    }
}


void dsmcTranslationalStageInterfacePatch::readProperties()
{
    sourceRegion_ =
        propsDict_.lookupOrDefault<word>("sourceRegion", word("source"));
    targetRegion_ =
        propsDict_.lookupOrDefault<word>("targetRegion", word("target"));
    sourceRegionVelocity_ =
        propsDict_.lookupOrDefault<vector>
        (
            "sourceRegionVelocity",
            vector::zero
        );
    targetRegionVelocity_ =
        propsDict_.lookupOrDefault<vector>
        (
            "targetRegionVelocity",
            vector::zero
        );
    sourceRegionWidth_ =
        propsDict_.lookupOrDefault<scalar>("sourceRegionWidth", 1.0);
    targetRegionWidth_ =
        propsDict_.lookupOrDefault<scalar>("targetRegionWidth", 1.0);
    positionMapping_ =
        propsDict_.lookupOrDefault<word>
        (
            "positionMapping",
            word("normalizedRandomPhase")
        );
    randomPhaseOperator_ =
        propsDict_.lookupOrDefault<word>
        (
            "randomPhaseOperator",
            word("addModulo")
        );
    randomPhaseMode_ =
        propsDict_.lookupOrDefault<word>
        (
            "randomPhaseMode",
            word("perTransfer")
        );
    randomPhaseSeed_ =
        propsDict_.lookupOrDefault<label>("randomPhaseSeed", 13579);
    randomPhaseValue_ =
        propsDict_.lookupOrDefault<scalar>("randomPhaseValue", 0.0);
    transformPosition_ =
        propsDict_.lookupOrDefault<bool>("transformPosition", true);
    approximateFaceCentreInsertion_ =
        propsDict_.lookupOrDefault<bool>("approximateFaceCentreInsertion", true);
    approximateInsertionFraction_ =
        propsDict_.lookupOrDefault<scalar>
        (
            "approximateInsertionFraction",
            2.0e-2
        );
    transformVelocity_ =
        propsDict_.lookupOrDefault<bool>("transformVelocity", true);
    debug_ = propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);

    if (approximateInsertionFraction_ < 0.0)
    {
        FatalErrorInFunction
            << "approximateInsertionFraction must be non-negative for patch "
            << patchName_
            << ", but got " << approximateInsertionFraction_
            << exit(FatalError);
    }

    if (sourceRegionWidth_ <= ROOTVSMALL || targetRegionWidth_ <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "sourceRegionWidth and targetRegionWidth must be positive for "
            << "patch " << patchName_
            << ", got sourceRegionWidth=" << sourceRegionWidth_
            << " targetRegionWidth=" << targetRegionWidth_
            << exit(FatalError);
    }

    if (positionMapping_ != "normalizedRandomPhase")
    {
        FatalErrorInFunction
            << "Unsupported positionMapping '" << positionMapping_
            << "' for " << typeName << ". Supported: normalizedRandomPhase."
            << exit(FatalError);
    }

    if (randomPhaseOperator_ != "addModulo")
    {
        FatalErrorInFunction
            << "Unsupported randomPhaseOperator '" << randomPhaseOperator_
            << "' for " << typeName << ". Supported: addModulo."
            << exit(FatalError);
    }

    if
    (
        randomPhaseMode_ != "perTransfer"
     && randomPhaseMode_ != "perTimeStep"
     && randomPhaseMode_ != "fixed"
    )
    {
        FatalErrorInFunction
            << "Unsupported randomPhaseMode '" << randomPhaseMode_
            << "' for " << typeName
            << ". Supported: perTransfer, perTimeStep, fixed."
            << exit(FatalError);
    }

    while (randomPhaseValue_ < 0.0)
    {
        randomPhaseValue_ += 1.0;
    }
    while (randomPhaseValue_ >= 1.0)
    {
        randomPhaseValue_ -= 1.0;
    }

    const polyPatch& sourcePatch = mesh_.boundaryMesh()[patchId_];
    const polyPatch& targetPatch = mesh_.boundaryMesh()[neighbPatchId_];

    if
    (
        isA<cyclicAMIPolyPatch>(sourcePatch)
     || isA<cyclicAMIPolyPatch>(targetPatch)
    )
    {
        FatalErrorInFunction
            << "Patch pair (" << patchName_ << ", " << neighbPatchName_ << ") "
            << "contains cyclicAMI topology, but "
            << "dsmcTranslationalStageInterfacePatch is only implemented for "
            << "strictly matched cyclic/cyclic interfaces."
            << exit(FatalError);
    }

    const bool isStrictCyclicPair =
        isA<cyclicPolyPatch>(sourcePatch)
     && isA<cyclicPolyPatch>(targetPatch);
    const bool isPlainPatchPair =
        !isA<cyclicPolyPatch>(sourcePatch)
     && !isA<cyclicPolyPatch>(targetPatch)
     && !isA<cyclicAMIPolyPatch>(sourcePatch)
     && !isA<cyclicAMIPolyPatch>(targetPatch);

    if (!isStrictCyclicPair && !isPlainPatchPair)
    {
        FatalErrorInFunction
            << "Patch pair (" << patchName_ << ", " << neighbPatchName_ << ") "
            << "must be either strict cyclic/cyclic or plain patch/patch for "
            << "dsmcTranslationalStageInterfacePatch."
            << exit(FatalError);
    }
}


vector dsmcTranslationalStageInterfacePatch::transformedVelocity
(
    const vector& U,
    const vector& fromRegionVelocity,
    const vector& toRegionVelocity
) const
{
    return U + fromRegionVelocity - toRegionVelocity;
}


scalar dsmcTranslationalStageInterfacePatch::transferPhase() const
{
    if (randomPhaseMode_ == "fixed")
    {
        return randomPhaseValue_;
    }

    if (randomPhaseMode_ == "perTimeStep")
    {
        const label currentTimeIndex = time_.timeIndex();

        if (cachedRandomPhaseTimeIndex_ != currentTimeIndex)
        {
            cachedRandomPhase_ = cloud_.rndGen().sample01<scalar>();
            cachedRandomPhaseTimeIndex_ = currentTimeIndex;
        }

        return cachedRandomPhase_;
    }

    return cloud_.rndGen().sample01<scalar>();
}


point dsmcTranslationalStageInterfacePatch::transferPatchParcel
(
    dsmcParcel& p,
    const bool fromSideA
) const
{
    const point positionBeforeTransfer = p.position();
    const label sourcePatchId = fromSideA ? patchId_ : neighbPatchId_;
    const label targetPatchId = fromSideA ? neighbPatchId_ : patchId_;
    const polyPatch& sourcePatch = mesh_.boundaryMesh()[sourcePatchId];
    const polyPatch& targetPatch = mesh_.boundaryMesh()[targetPatchId];

    scalar sourceMinX = GREAT;
    scalar sourceMaxX = -GREAT;
    scalar targetMinX = GREAT;
    scalar targetMaxX = -GREAT;

    patchXBounds(mesh_, sourcePatch, sourceMinX, sourceMaxX);
    patchXBounds(mesh_, targetPatch, targetMinX, targetMaxX);

    const scalar sourceWidth = max(sourceMaxX - sourceMinX, ROOTVSMALL);
    const scalar targetWidth = max(targetMaxX - targetMinX, ROOTVSMALL);

    scalar xi = (positionBeforeTransfer.x() - sourceMinX)/sourceWidth;
    while (xi < 0.0)
    {
        xi += 1.0;
    }
    while (xi >= 1.0)
    {
        xi -= 1.0;
    }

    const scalar phase = transferPhase();
    scalar xiTarget = xi + phase;
    while (xiTarget >= 1.0)
    {
        xiTarget -= 1.0;
    }
    while (xiTarget < 0.0)
    {
        xiTarget += 1.0;
    }

    const scalar targetX = targetMinX + xiTarget*targetWidth;

    label targetLocalFaceI = -1;
    scalar bestDistance = GREAT;

    for(label i = 0; i < targetPatch.size(); ++i)
    {
        const label globalFaceI = targetPatch.start() + i;
        const face& f = mesh_.faces()[globalFaceI];

        scalar minX = GREAT;
        scalar maxX = -GREAT;

        forAll(f, fp)
        {
            const scalar x = mesh_.points()[f[fp]].x();
            minX = min(minX, x);
            maxX = max(maxX, x);
        }

        if (targetX >= minX - SMALL && targetX <= maxX + SMALL)
        {
            targetLocalFaceI = i;
            break;
        }

        const scalar centreDistance =
            mag(mesh_.faceCentres()[globalFaceI].x() - targetX);

        if (centreDistance < bestDistance)
        {
            bestDistance = centreDistance;
            targetLocalFaceI = i;
        }
    }

    if (targetLocalFaceI < 0)
    {
        FatalErrorInFunction
            << "Unable to choose target face on patch " << targetPatch.name()
            << " for targetX=" << targetX << " targetWidth=" << targetWidth
            << " targetMinX=" << targetMinX << " targetMaxX=" << targetMaxX
            << exit(FatalError);
    }

    p.face() = targetPatch.start() + targetLocalFaceI;
    p.cell() = targetPatch.faceCells()[targetLocalFaceI];
    p.tetFace() = p.face();
    p.tetPt() = mesh_.faces()[p.tetFace()].size() - 1 - p.tetPt();

    const point& targetFaceCentre = mesh_.faceCentres()[p.face()];
    const point& ownerCellCentre = mesh_.cellCentres()[p.cell()];

    point mappedCrossPosition
    (
        targetX,
        targetFaceCentre.y(),
        targetFaceCentre.z()
    );

    if (transformPosition_)
    {
        p.position() = mappedCrossPosition;

        vector inwardDirection = ownerCellCentre - p.position();

        if (mag(inwardDirection) <= ROOTVSMALL)
        {
            const vector faceNormal = mesh_.faceAreas()[p.face()];

            if (mag(faceNormal) > ROOTVSMALL)
            {
                inwardDirection = -faceNormal;
            }
        }

        const scalar inwardMag = mag(inwardDirection);

        if (inwardMag > ROOTVSMALL)
        {
            const scalar faceLengthScale =
                Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));
            const scalar insertionDistance =
                approximateFaceCentreInsertion_
              ? max(ROOTVSMALL, approximateInsertionFraction_*faceLengthScale)
              : max(ROOTVSMALL, 1.0e-6*faceLengthScale);

            p.position() +=
                (insertionDistance/inwardMag)*inwardDirection;
        }

        if (mag(p.position() - positionBeforeTransfer) > ROOTVSMALL)
        {
            p.initCellFacePt();
        }
    }

    return mappedCrossPosition;
}


point dsmcTranslationalStageInterfacePatch::transferParcel(dsmcParcel& p) const
{
    const label faceI = p.face();
    const label patchI = mesh_.boundaryMesh().whichPatch(faceI);
    const point positionBeforeTransfer = p.position();

    if (patchI == -1)
    {
        FatalErrorInFunction
            << "Face " << faceI
            << " is not on a boundary patch during stage-interface transfer."
            << exit(FatalError);

        return positionBeforeTransfer;
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchI];

    if (!isA<cyclicPolyPatch>(patch))
    {
        if (patchI == patchId_)
        {
            return transferPatchParcel(p, true);
        }
        else if (patchI == neighbPatchId_)
        {
            return transferPatchParcel(p, false);
        }
        else
        {
            FatalErrorInFunction
                << "Patch " << patch.name()
                << " is not part of the translational stage pair ("
                << patchName_ << ", " << neighbPatchName_ << ")."
                << exit(FatalError);

            return positionBeforeTransfer;
        }
    }

    const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(patch);

    p.face() = cpp.transformGlobalFace(faceI);
    p.cell() = mesh_.faceOwner()[p.face()];
    p.tetFace() = p.face();
    p.tetPt() = mesh_.faces()[p.tetFace()].size() - 1 - p.tetPt();

    const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
    const label patchFaceI = receiveCpp.whichFace(p.face());

    if (patchFaceI == -1)
    {
        FatalErrorInFunction
            << "Unable to locate transformed face " << p.face()
            << " on neighbour cyclic patch " << receiveCpp.name()
            << exit(FatalError);

        return positionBeforeTransfer;
    }

    point mappedCrossPosition = positionBeforeTransfer;
    receiveCpp.transformPosition(mappedCrossPosition, patchFaceI);

    const point& targetFaceCentre = mesh_.faceCentres()[p.face()];
    const point& ownerCellCentre = mesh_.cellCentres()[p.cell()];

    if (transformPosition_)
    {
        if (approximateFaceCentreInsertion_)
        {
            vector insertionDirection = ownerCellCentre - targetFaceCentre;

            if (mag(insertionDirection) <= ROOTVSMALL)
            {
                const vector faceNormal = mesh_.faceAreas()[p.face()];

                if (mag(faceNormal) > ROOTVSMALL)
                {
                    insertionDirection = -faceNormal;
                }
            }

            if (mag(insertionDirection) > ROOTVSMALL)
            {
                const scalar faceLengthScale =
                    Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));
                const scalar insertionDistance =
                    max
                    (
                        ROOTVSMALL,
                        approximateInsertionFraction_*faceLengthScale
                    );

                p.position() =
                    targetFaceCentre
                  + (insertionDistance/mag(insertionDirection))*insertionDirection;
            }
            else
            {
                p.position() = targetFaceCentre;
            }
        }
        else
        {
            p.position() = mappedCrossPosition;
        }

        vector inwardDirection = ownerCellCentre - p.position();

        if (mag(inwardDirection) <= ROOTVSMALL)
        {
            const vector faceNormal = mesh_.faceAreas()[p.face()];

            if (mag(faceNormal) > ROOTVSMALL)
            {
                inwardDirection = -faceNormal;
            }
        }

        const scalar inwardMag = mag(inwardDirection);

        if (inwardMag > ROOTVSMALL)
        {
            const scalar faceLengthScale =
                Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));
            const scalar pushDistance =
                max(ROOTVSMALL, 1.0e-6*faceLengthScale);

            p.position() += (pushDistance/inwardMag)*inwardDirection;
        }

        if (mag(p.position() - positionBeforeTransfer) > ROOTVSMALL)
        {
            p.initCellFacePt();
        }
    }

    return mappedCrossPosition;
}


dsmcTranslationalStageInterfacePatch::dsmcTranslationalStageInterfacePatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    sourceRegion_("source"),
    targetRegion_("target"),
    sourceRegionVelocity_(vector::zero),
    targetRegionVelocity_(vector::zero),
    sourceRegionWidth_(1.0),
    targetRegionWidth_(1.0),
    positionMapping_("normalizedRandomPhase"),
    randomPhaseOperator_("addModulo"),
    randomPhaseMode_("perTransfer"),
    randomPhaseSeed_(13579),
    randomPhaseValue_(0.0),
    cachedRandomPhase_(0.0),
    cachedRandomPhaseTimeIndex_(-1),
    transformPosition_(true),
    approximateFaceCentreInsertion_(true),
    approximateInsertionFraction_(2.0e-2),
    transformVelocity_(true),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCountA_(0),
    debugPrintCountB_(0),
    transferCountAtoB_(0),
    transferCountBtoA_(0),
    transformCountAtoB_(0),
    transformCountBtoA_(0),
    lastSummaryTimeIndex_(-1)
{
    writeInTimeDir_ = false;
    writeInCase_ = false;
    readProperties();
}


dsmcTranslationalStageInterfacePatch::~dsmcTranslationalStageInterfacePatch()
{}


void dsmcTranslationalStageInterfacePatch::calculateProperties()
{}


void dsmcTranslationalStageInterfacePatch::initialConfiguration()
{}


void dsmcTranslationalStageInterfacePatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    (void)td;

    const label faceI = p.face();
    const label oldPatchI = mesh_.boundaryMesh().whichPatch(faceI);
    const point oldPosition = p.position();
    const vector oldVelocity = p.U();
    const scalar oldRadialWeight = p.RWF();
    const label oldCell = p.cell();
    const label oldTetFace = p.tetFace();
    const label oldTetPt = p.tetPt();

    word side("none");
    bool fromSideA = false;
    word fromRegion(sourceRegion_);
    word toRegion(targetRegion_);
    vector fromRegionVelocity(sourceRegionVelocity_);
    vector toRegionVelocity(targetRegionVelocity_);
    label expectedNewPatchId = -1;

    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    if (sourceIndexA != -1)
    {
        side = "A";
        fromSideA = true;
        fromRegion = sourceRegion_;
        toRegion = targetRegion_;
        fromRegionVelocity = sourceRegionVelocity_;
        toRegionVelocity = targetRegionVelocity_;
        expectedNewPatchId = neighbPatchId_;
    }
    else if (sourceIndexB != -1)
    {
        side = "B";
        fromSideA = false;
        fromRegion = targetRegion_;
        toRegion = sourceRegion_;
        fromRegionVelocity = targetRegionVelocity_;
        toRegionVelocity = sourceRegionVelocity_;
        expectedNewPatchId = patchId_;
    }
    else
    {
        return;
    }

    const point mappedCrossPosition = transferParcel(p);
    const scalar newRadialWeight = p.RWF();

    if
    (
        mag(newRadialWeight - oldRadialWeight)
      > SMALL*max(scalar(1), mag(oldRadialWeight))
    )
    {
        FatalErrorInFunction
            << "Translational stage-interface transfer changed parcel "
            << "radialWeight on a strictly matched cyclic interface." << nl
            << "Interface patch pair: (" << patchName_ << ", "
            << neighbPatchName_ << ")" << nl
            << "Transfer face: " << faceI << nl
            << "Old radialWeight: " << oldRadialWeight << nl
            << "New radialWeight: " << newRadialWeight
            << exit(FatalError);
    }

    const label newPatchI = mesh_.boundaryMesh().whichPatch(p.face());

    if (newPatchI != expectedNewPatchId)
    {
        FatalErrorInFunction
            << "Translational stage-interface transfer sent parcel to "
            << "unexpected patch." << nl
            << "Interface patch pair: (" << patchName_ << ", "
            << neighbPatchName_ << ")" << nl
            << "Transfer side: " << side << nl
            << "Old face: " << faceI << ", old patch: "
            << (oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none")) << nl
            << "New face: " << p.face() << ", new patch: "
            << (newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none")) << nl
            << "Expected patch: " << mesh_.boundaryMesh()[expectedNewPatchId].name()
            << exit(FatalError);
    }

    if (fromSideA)
    {
        ++transferCountAtoB_;
    }
    else
    {
        ++transferCountBtoA_;
    }

    const vector expectedVelocity =
        transformVelocity_
      ? transformedVelocity(oldVelocity, fromRegionVelocity, toRegionVelocity)
      : oldVelocity;
    const vector expectedDeltaU = expectedVelocity - oldVelocity;

    if (transformVelocity_)
    {
        p.U() = expectedVelocity;

        if (fromSideA)
        {
            ++transformCountAtoB_;
        }
        else
        {
            ++transformCountBtoA_;
        }
    }

    label* sideDebugPrintCount = fromSideA ? &debugPrintCountA_ : &debugPrintCountB_;

    if (debug_ && *sideDebugPrintCount < maxDebugPrints_)
    {
        const word oldPatchName =
            oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none");
        const word newPatchName =
            newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none");

        Info<< "dsmcTranslationalStageInterfacePatch debug: patch="
            << patchName_
            << " side=" << side
            << " oldPatch=" << oldPatchName
            << " newPatch=" << newPatchName
            << " oldFace=" << faceI
            << " newFace=" << p.face()
            << " oldCell=" << oldCell
            << " newCell=" << p.cell()
            << " oldTetFace=" << oldTetFace
            << " newTetFace=" << p.tetFace()
            << " oldTetPt=" << oldTetPt
            << " newTetPt=" << p.tetPt()
            << " oldPosition=" << oldPosition
            << " newPosition=" << p.position()
            << " mappedCrossPosition=" << mappedCrossPosition
            << " oldU=" << oldVelocity
            << " newU=" << p.U()
            << " expectedDeltaU=" << expectedDeltaU
            << " actualDeltaU=" << (p.U() - oldVelocity)
            << " fromRegion=" << fromRegion
            << " toRegion=" << toRegion
            << " fromRegionVelocity=" << fromRegionVelocity
            << " toRegionVelocity=" << toRegionVelocity
            << " oldRadialWeight=" << oldRadialWeight
            << " newRadialWeight=" << newRadialWeight
            << nl;

        ++(*sideDebugPrintCount);
    }

    if (debug_ && time_.outputTime() && lastSummaryTimeIndex_ != time_.timeIndex())
    {
        Info<< "dsmcTranslationalStageInterfacePatch summary: patch="
            << patchName_
            << " time=" << time_.timeName()
            << " proc=" << Pstream::myProcNo()
            << " localAtoBTransfers=" << transferCountAtoB_
            << " localBtoATransfers=" << transferCountBtoA_
            << " localAtoBVelocityTransforms=" << transformCountAtoB_
            << " localBtoAVelocityTransforms=" << transformCountBtoA_
            << " debugPrintsA=" << debugPrintCountA_
            << " debugPrintsB=" << debugPrintCountB_
            << nl;

        lastSummaryTimeIndex_ = time_.timeIndex();
    }
}


void dsmcTranslationalStageInterfacePatch::output
(
    const fileName&,
    const fileName&
)
{}


void dsmcTranslationalStageInterfacePatch::updateProperties
(
    const dictionary& newDict
)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
}

} // End namespace Foam

// ************************************************************************* //
