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

#include "dsmcParcel.H"
#include "dsmcCloud.H"
#include "meshTools.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * * //

Foam::label Foam::dsmcParcel::TrackedParcel::nDELETED = 0;


namespace Foam
{
namespace
{

struct BoundaryEventStepDiagnostics
{
    scalar timeValue_;
    label totalEvents_;
    label cyclicEvents_;
    label wallEvents_;
    label genericEvents_;
    label deletedEvents_;
    label retainedEvents_;
    label nextDiagnosticThreshold_;
    HashTable<label> patchCounts_;

    BoundaryEventStepDiagnostics()
    :
        timeValue_(-GREAT),
        totalEvents_(0),
        cyclicEvents_(0),
        wallEvents_(0),
        genericEvents_(0),
        deletedEvents_(0),
        retainedEvents_(0),
        nextDiagnosticThreshold_(0),
        patchCounts_(16)
    {}
};


BoundaryEventStepDiagnostics& boundaryEventStepDiagnostics()
{
    // 每个 MPI rank 本地维护一个“单时间步边界事件计数器”，
    // 用于聚合 cyclic / wall / generic 三类事件数及各 patch 命中次数。
    // 这样即便最后一条 before/after trace 没有直接落在真实卡点上，
    // 也能先判断哪一类边界交互在当前时间步内主导了事件累积。
    static BoundaryEventStepDiagnostics diagnostics;

    return diagnostics;
}


label& parcelTrackSubStepRescueWarningCount()
{
    // 每个 MPI rank 本地维护 rescue 警告打印计数。
    // 该计数只用于限制日志量，不参与物理统计；真正的 rescue 总数仍通过
    // Cloud::trackingRescue() 进入 OpenFOAM 原有 tracking rescue 计数器。
    static label warningCount = 0;

    return warningCount;
}


bool boundaryEventSummaryEnabled
(
    dsmcParcel::trackingData& td
)
{
    // 聚合 summary 用于判断整步内哪一类边界事件在快速累积，
    // 与 before/after trace 的“最后一个具体事件落点”形成互补。
    const dictionary& controlDict = td.cloud().controlDict();

    if
    (
        !controlDict.lookupOrDefault<bool>
        (
            "traceBoundaryEventSummary",
            false
        )
    )
    {
        return false;
    }

    const scalar boundaryTraceStartTime =
        controlDict.lookupOrDefault<scalar>
        (
            "boundaryTraceStartTime",
            0.0
        );

    return td.cloud().pMesh().time().value() >= boundaryTraceStartTime;
}


label boundaryEventDiagnosticInterval
(
    dsmcParcel::trackingData& td
)
{
    // 允许用户控制每累计多少个 boundary event 打印一次中间 summary。
    const label interval =
        td.cloud().controlDict().lookupOrDefault<label>
        (
            "boundaryEventDiagnosticInterval",
            5000
        );

    return interval > 1 ? interval : 1;
}


void resetBoundaryEventStepDiagnostics
(
    BoundaryEventStepDiagnostics& diagnostics,
    const scalar timeValue,
    const label diagnosticInterval
)
{
    diagnostics.timeValue_ = timeValue;
    diagnostics.totalEvents_ = 0;
    diagnostics.cyclicEvents_ = 0;
    diagnostics.wallEvents_ = 0;
    diagnostics.genericEvents_ = 0;
    diagnostics.deletedEvents_ = 0;
    diagnostics.retainedEvents_ = 0;
    diagnostics.nextDiagnosticThreshold_ = diagnosticInterval;
    diagnostics.patchCounts_.clear();
}


void writeBoundaryEventStepSummary
(
    const BoundaryEventStepDiagnostics& diagnostics,
    const word& reason
)
{
    Info<< "DSMC boundary event summary: "
        << "reason=" << reason
        << " time=" << diagnostics.timeValue_
        << " totalEvents=" << diagnostics.totalEvents_
        << " cyclicEvents=" << diagnostics.cyclicEvents_
        << " wallEvents=" << diagnostics.wallEvents_
        << " genericEvents=" << diagnostics.genericEvents_
        << " retainedEvents=" << diagnostics.retainedEvents_
        << " deletedEvents=" << diagnostics.deletedEvents_
        << " patchCounts=(";

    bool firstPatch = true;

    forAllConstIter(HashTable<label>, diagnostics.patchCounts_, iter)
    {
        if (!firstPatch)
        {
            Info<< ",";
        }

        Info<< iter.key() << ":" << iter();
        firstPatch = false;
    }

    Info<< ")" << endl;
}


void ensureBoundaryEventStepDiagnosticsCurrent
(
    dsmcParcel::trackingData& td
)
{
    if (!boundaryEventSummaryEnabled(td))
    {
        return;
    }

    BoundaryEventStepDiagnostics& diagnostics = boundaryEventStepDiagnostics();
    const scalar currentTime = td.cloud().pMesh().time().value();
    const label diagnosticInterval = boundaryEventDiagnosticInterval(td);

    if (diagnostics.timeValue_ < -0.5*GREAT)
    {
        resetBoundaryEventStepDiagnostics
        (
            diagnostics,
            currentTime,
            diagnosticInterval
        );
        return;
    }

    if (mag(currentTime - diagnostics.timeValue_) > SMALL)
    {
        if (diagnostics.totalEvents_ > 0)
        {
            writeBoundaryEventStepSummary
            (
                diagnostics,
                "timeStepCompleted"
            );
        }

        resetBoundaryEventStepDiagnostics
        (
            diagnostics,
            currentTime,
            diagnosticInterval
        );
    }
}


void recordBoundaryEventEntry
(
    dsmcParcel::trackingData& td,
    const word& patchName,
    const word& category
)
{
    if (!boundaryEventSummaryEnabled(td))
    {
        return;
    }

    ensureBoundaryEventStepDiagnosticsCurrent(td);

    BoundaryEventStepDiagnostics& diagnostics = boundaryEventStepDiagnostics();

    diagnostics.totalEvents_++;

    if (category == "cyclic")
    {
        diagnostics.cyclicEvents_++;
    }
    else if (category == "wall")
    {
        diagnostics.wallEvents_++;
    }
    else
    {
        diagnostics.genericEvents_++;
    }

    label patchCount = 0;

    if (diagnostics.patchCounts_.found(patchName))
    {
        patchCount = diagnostics.patchCounts_[patchName];
    }

    diagnostics.patchCounts_.set
    (
        patchName,
        patchCount + 1
    );

    if (diagnostics.totalEvents_ >= diagnostics.nextDiagnosticThreshold_)
    {
        writeBoundaryEventStepSummary
        (
            diagnostics,
            "intervalReached"
        );

        diagnostics.nextDiagnosticThreshold_ +=
            boundaryEventDiagnosticInterval(td);
    }
}


void recordBoundaryEventExit
(
    dsmcParcel::trackingData& td
)
{
    if (!boundaryEventSummaryEnabled(td))
    {
        return;
    }

    BoundaryEventStepDiagnostics& diagnostics = boundaryEventStepDiagnostics();

    if (td.keepParticle)
    {
        diagnostics.retainedEvents_++;
    }
    else
    {
        diagnostics.deletedEvents_++;
    }
}

bool boundaryTraceEnabled
(
    dsmcParcel::trackingData& td
)
{
    // 这一层运行期开关用于定位“trackToFace() 已经命中 boundary face，
    // 但随后卡在 patch interaction 入口之后”的病态。
    // 通过记录 before/after trace，可以直接判断停在 cyclic、wall
    // 还是 generic patch 的哪一段调用路径上。
    const dictionary& controlDict = td.cloud().controlDict();

    if
    (
        !controlDict.lookupOrDefault<bool>
        (
            "traceBoundaryModelControl",
            false
        )
    )
    {
        return false;
    }

    const scalar boundaryTraceStartTime =
        controlDict.lookupOrDefault<scalar>
        (
            "boundaryTraceStartTime",
            0.0
        );

    return td.cloud().pMesh().time().value() >= boundaryTraceStartTime;
}


label nextBoundaryTraceEventId()
{
    // 每个 MPI rank 本地维护事件编号即可。
    // 这里只需要保证同一次边界交互的 before/after 共享同一个 eventId。
    static label boundaryTraceEventCounter = 0;

    return ++boundaryTraceEventCounter;
}


bool boundaryTraceEventWithinLimit
(
    dsmcParcel::trackingData& td,
    const label eventId
)
{
    // 限制 trace 事件数量，避免日志无限膨胀。
    const label maxBoundaryTraceEvents =
        td.cloud().controlDict().lookupOrDefault<label>
        (
            "maxBoundaryTraceEvents",
            2000
        );

    return eventId <= maxBoundaryTraceEvents;
}


bool boundaryTracePatchAllowed
(
    dsmcParcel::trackingData& td,
    const word& patchName
)
{
    // 可选 patch 过滤器：
    // - 未提供 traceBoundaryPatchName 时追所有 patch；
    // - 提供后只追指定 patch。
    const word tracedPatchName =
        td.cloud().controlDict().lookupOrDefault<word>
        (
            "traceBoundaryPatchName",
            word("")
        );

    return tracedPatchName.empty() || tracedPatchName == patchName;
}


bool boundaryTraceFaceAllowed
(
    dsmcParcel::trackingData& td,
    const label faceI
)
{
    // face 级可选过滤：
    // - 未提供 traceBoundaryFaces 时放行所有 face；
    // - 提供后只记录列出的全局 face 编号。
    labelList tracedFaces(0);

    if (!td.cloud().controlDict().readIfPresent("traceBoundaryFaces", tracedFaces))
    {
        return true;
    }

    return findIndex(tracedFaces, faceI) != -1;
}


bool boundaryTraceCellAllowed
(
    dsmcParcel::trackingData& td,
    const label cellI
)
{
    // cell 级可选过滤，与 face 过滤配合使用：
    // - 未提供 traceBoundaryCells 时放行所有 cell；
    // - 提供后只记录列出的 cell。
    labelList tracedCells(0);

    if (!td.cloud().controlDict().readIfPresent("traceBoundaryCells", tracedCells))
    {
        return true;
    }

    return findIndex(tracedCells, cellI) != -1;
}


void writeBoundaryTrace
(
    const dsmcParcel& p,
    dsmcParcel::trackingData& td,
    const label eventId,
    const word& stage,
    const word& patchName,
    const word& patchType,
    const label modelId
)
{
    if (!boundaryTraceEventWithinLimit(td, eventId))
    {
        return;
    }

    const polyMesh& mesh = td.cloud().pMesh();
    const label faceI = p.face();
    const label patchI =
        faceI >= 0 ? mesh.boundaryMesh().whichPatch(faceI) : -1;

    // 统一输出几何定位、时间推进、物理状态与控制流状态。
    // 若程序停在某个 boundary model 内部，最后一条 before/after trace
    // 通常就能直接指向具体 patch/model。
    Info<< "DSMC boundary trace: "
        << "event=" << eventId
        << " stage=" << stage
        << " time=" << mesh.time().value()
        << " patch=" << patchName
        << " patchType=" << patchType
        << " modelId=" << modelId
        << " face=" << faceI
        << " patchIndex=" << patchI
        << " cell=" << p.cell()
        << " tetFace=" << p.tetFace()
        << " tetPt=" << p.tetPt()
        << " stepFraction=" << p.stepFraction()
        << " position=" << p.position()
        << " U=" << p.U()
        << " keepParticle=" << td.keepParticle
        << " switchProcessor=" << td.switchProcessor
        << " isFree=" << p.isFree()
        << endl;
}

void applyPostWallInwardPush
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    // 处理“壁面反射后粒子几何上仍贴在 wall face 上”的病态。
    // 这里统一沿 owner cell 内部方向施加一个极小 inward push：
    //     x_new = x_hit + beta_wall * L_f * n_in
    // 只修正几何稳健性，不改变 wall interaction 的物理公式。
    if (!td.keepParticle || !p.isFree())
    {
        return;
    }

    const polyMesh& mesh = td.cloud().pMesh();
    const label faceI = p.face();
    const label cellI = p.cell();

    if (faceI < 0 || cellI < 0)
    {
        return;
    }

    if (mesh.boundaryMesh().whichPatch(faceI) == -1)
    {
        return;
    }

    vector inwardDirection = mesh.cellCentres()[cellI] - p.position();

    if (mag(inwardDirection) <= ROOTVSMALL)
    {
        // 若 cell centre 方向退化，则回退到 boundary face 外法向的反方向。
        const vector faceNormal = mesh.faceAreas()[faceI];

        if (mag(faceNormal) > ROOTVSMALL)
        {
            inwardDirection = -faceNormal;
        }
    }

    const scalar inwardMag = mag(inwardDirection);

    if (inwardMag <= ROOTVSMALL)
    {
        return;
    }

    const scalar wallPostCollisionPushScale =
        td.cloud().controlDict().lookupOrDefault<scalar>
        (
            "wallPostCollisionPushScale",
            1.0e-6
        );

    if (wallPostCollisionPushScale <= 0.0)
    {
        return;
    }

    const scalar faceLengthScale = Foam::sqrt(mag(mesh.faceAreas()[faceI]));
    const scalar pushDistance =
        max(ROOTVSMALL, wallPostCollisionPushScale*faceLengthScale);

    p.position() += (pushDistance/inwardMag)*inwardDirection;

    // inward push 后还需重建 tetFace/tetPt。
    // 否则粒子虽然位置离面，但仍可能保留上一张 boundary face 的局部 tet 拓扑，
    // 下一轮追踪继续从不一致的 face/tet 状态出发。
    p.initCellFacePt();
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dsmcParcel::move
(
    dsmcParcel::trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    if (isFree())
    {
        const polyMesh& mesh = td.cloud().pMesh();
        const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

        if (newParcel() != -1)
        {
            // note: this justifies that freshly inserted parcels should be
            // tracked as if they had passed the boundary face on which they
            // have been inserted in the time step in which they are inserted.
            stepFraction() = td.cloud().rndGen().sample01<scalar>();
            newParcel() = -1;
        }

        //scalar tEnd = (1.0 - stepFraction())*trackTime; // OLD FORMULATION
        label orgCell = cell(); // NEW VINCENT
        scalar tEnd = (1.0 - stepFraction())*td.cloud().deltaTValue(orgCell); // NEW VINCENT

        // For reduced-D cases, the velocity used to track needs to be
        // constrained, but the actual U_ of the parcel must not be
        // altered or used, as it is altered by patch interactions one
        // needs to retain its 3D value for collision purposes.
        vector Utracking = U_;
        const scalar totalTrackTime = td.cloud().deltaTValue(orgCell);

        // 单颗粒单时间步最大追踪子步数保护。
        // 用于把几何容差、周期面重入或极小 trackFraction 造成的静默满 CPU
        // 转成显式报错。
        const label maxParcelTrackSubSteps =
            td.cloud().controlDict().lookupOrDefault<label>
            (
                "maxParcelTrackSubSteps",
                td.cloud().rotatingReferenceFrameActive() ? 20000 : 10000
            );
        const bool maxParcelTrackSubStepsRescue =
            td.cloud().controlDict().lookupOrDefault<bool>
            (
                "maxParcelTrackSubStepsRescue",
                false
            );
        const label maxParcelTrackSubStepsRescueTrigger =
            max
            (
                td.cloud().controlDict().lookupOrDefault<label>
                (
                    "maxParcelTrackSubStepsRescueTrigger",
                    200
                ),
                label(1)
            );
        const scalar maxParcelTrackSubStepsRescueMaxTrackFraction =
            max
            (
                td.cloud().controlDict().lookupOrDefault<scalar>
                (
                    "maxParcelTrackSubStepsRescueMaxTrackFraction",
                    1e-4
                ),
                scalar(0)
            );
        const scalar rescueCellCentreFraction =
            min
            (
                max
                (
                    td.cloud().controlDict().lookupOrDefault<scalar>
                    (
                        "maxParcelTrackSubStepsRescueCellCentreFraction",
                        0.05
                    ),
                    scalar(0)
                ),
                scalar(1)
            );
        const label maxParcelTrackSubStepsRescueWarnings =
            max
            (
                td.cloud().controlDict().lookupOrDefault<label>
                (
                    "maxParcelTrackSubStepsRescueWarnings",
                    20
                ),
                label(0)
            );
        label parcelTrackSubStepCount = 0;

        // 下面这些变量记录“上一轮已经执行过的 trackToFace()”的关键状态。
        // watchdog 触发点位于下一轮 trackToFace() 之前，因此若不保存上一轮
        // 信息，就只能知道当前粒子在哪里，却不知道上一轮到底推进了多少。
        scalar lastTrackFraction = -GREAT;
        scalar lastAttemptDt = -GREAT;
        scalar lastEffectiveDt = -GREAT;
        scalar lastRemainingTrackTimeBeforeTrack = -GREAT;
        point lastStartPosition = vector::zero;
        point lastRequestedEndPosition = vector::zero;
        point lastEndPosition = vector::zero;
        label lastStartCell = -1;
        label lastEndCell = -1;
        label lastFace = -1;
        bool lastEndedOnBoundaryFace = false;

        while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
        {
            orgCell = cell();

            const bool inRotatingRegion =
                td.cloud().isInRotatingRegion(orgCell);

            const scalar baseDt = td.cloud().deltaTValue(orgCell);
            scalar dt = tEnd; // NEW VINCENT

            if (inRotatingRegion)
            {
                dt =
                    td.cloud().rotatingFrameTrackSubStep
                    (
                        orgCell,
                        position(),
                        U_,
                        tEnd
                    );

                const point startPos(position());
                const vector startU(U_);
                const vector a0 =
                    td.cloud().rotatingFrameAcceleration
                    (
                        startPos,
                        startU
                    );

                const vector Uhalf = startU + 0.5*dt*a0;
                Utracking = Uhalf;
            }
            else
            {
                Utracking = U_;
            }

            // 在调用 trackToFace() 之前做 watchdog 计数，防止外层 while
            // 因极小推进量或重复重入而无界增长。
            ++parcelTrackSubStepCount;

            const bool earlyRescueCheck =
            (
                maxParcelTrackSubStepsRescue
             && parcelTrackSubStepCount > maxParcelTrackSubStepsRescueTrigger
            );
            const bool watchdogExceeded =
                parcelTrackSubStepCount > maxParcelTrackSubSteps;

            if (earlyRescueCheck || watchdogExceeded)
            {
                const label faceI = face();
                const label patchI =
                    faceI >= 0 ? mesh.boundaryMesh().whichPatch(faceI) : -1;
                const word patchName =
                    patchI >= 0 ? pbMesh[patchI].name() : word("none");

                const bool smallProgressInternalLoop =
                (
                    earlyRescueCheck
                 && patchI == -1
                 && lastTrackFraction >= 0
                 && lastTrackFraction
                  <= maxParcelTrackSubStepsRescueMaxTrackFraction
                );

                if
                (
                    maxParcelTrackSubStepsRescue
                 && (smallProgressInternalLoop || watchdogExceeded)
                )
                {
                    const label cellI = cell();
                    const bool validCell = cellI >= 0 && cellI < mesh.nCells();
                    const point rescueStartPosition = position();
                    const scalar stepFractionBeforeRescue = stepFraction();
                    const scalar remainingTrackTimeBeforeRescue = tEnd;
                    scalar rescueDisplacement = 0.0;

                    if (validCell && rescueCellCentreFraction > SMALL)
                    {
                        // 对这类 watchdog 命中，日志已经表明 face=-1、patch=none，
                        // 即粒子没有处在一个可交给 patch 模型处理的真实边界事件上。
                        // 最小破坏性的处理是沿当前 cell 几何中心方向做极小回拉：
                        //
                        //     x_rescue = x_old + alpha (C_cell - x_old)
                        //
                        // 其中 alpha 由 controlDict 控制。这个公式保证修正方向
                        // 指向当前 owner cell 内部，而不是任意改变粒子速度或跨面
                        // 迁移拓扑；随后结束本时间步剩余追踪，避免同一颗粒子继续
                        // 在极小 trackFraction 循环中耗尽 CPU。
                        const vector toCellCentre =
                            mesh.cellCentres()[cellI] - position();
                        const vector displacement =
                            rescueCellCentreFraction*toCellCentre;

                        rescueDisplacement = mag(displacement);
                        position() += displacement;
                        face() = -1;
                        tetFace() = -1;
                        tetPt() = -1;
                        initCellFacePt();
                    }

                    tEnd = 0.0;
                    stepFraction() = 1.0;
                    td.cloud().trackingRescue();

                    label& rescueWarningCount =
                        parcelTrackSubStepRescueWarningCount();

                    if (rescueWarningCount < maxParcelTrackSubStepsRescueWarnings)
                    {
                        WarningInFunction
                            << "Rescued a DSMC parcel after exceeding "
                            << "maxParcelTrackSubSteps." << nl
                            << "This avoids a pathological tracking loop by "
                            << "ending only the remaining part of the current "
                            << "time step for this parcel." << nl
                            << "maxParcelTrackSubSteps=" << maxParcelTrackSubSteps << nl
                            << "maxParcelTrackSubStepsRescueTrigger="
                            << maxParcelTrackSubStepsRescueTrigger << nl
                            << "maxParcelTrackSubStepsRescueMaxTrackFraction="
                            << maxParcelTrackSubStepsRescueMaxTrackFraction << nl
                            << "smallProgressInternalLoop="
                            << smallProgressInternalLoop << nl
                            << "watchdogExceeded=" << watchdogExceeded << nl
                            << "parcelTrackSubStepCount="
                            << parcelTrackSubStepCount << nl
                            << "cell=" << cell() << nl
                            << "face=" << faceI << nl
                            << "patch=" << patchName << nl
                            << "stepFractionBeforeRescue="
                            << stepFractionBeforeRescue << nl
                            << "remainingTrackTimeBeforeRescue="
                            << remainingTrackTimeBeforeRescue << nl
                            << "cellDeltaT=" << baseDt << nl
                            << "positionBeforeRescue="
                            << rescueStartPosition << nl
                            << "positionAfterRescue=" << position() << nl
                            << "rescueCellCentreFraction="
                            << rescueCellCentreFraction << nl
                            << "rescueDisplacement="
                            << rescueDisplacement << nl
                            << "U=" << U_ << nl
                            << "Utracking=" << Utracking << nl
                            << "inRotatingRegion=" << inRotatingRegion << nl
                            << "newParcel=" << newParcel() << nl
                            << "classification=" << classification() << nl
                            << "lastTrackFraction=" << lastTrackFraction << nl
                            << "lastAttemptDt=" << lastAttemptDt << nl
                            << "lastEffectiveDt=" << lastEffectiveDt << nl
                            << "lastRemainingTrackTimeBeforeTrack="
                            << lastRemainingTrackTimeBeforeTrack << nl
                            << "lastStartCell=" << lastStartCell << nl
                            << "lastEndCell=" << lastEndCell << nl
                            << "lastFace=" << lastFace << nl
                            << "lastEndedOnBoundaryFace="
                            << lastEndedOnBoundaryFace << nl
                            << "lastStartPosition=" << lastStartPosition << nl
                            << "lastRequestedEndPosition="
                            << lastRequestedEndPosition << nl
                            << "lastEndPosition=" << lastEndPosition
                            << endl;

                        ++rescueWarningCount;
                    }

                    break;
                }

                if (watchdogExceeded)
                {
                    FatalErrorInFunction
                        << "Exceeded maxParcelTrackSubSteps while moving a DSMC parcel."
                        << nl
                        << "This indicates a pathological tracking loop rather than a "
                        << "normal parcel move." << nl
                        << "maxParcelTrackSubSteps=" << maxParcelTrackSubSteps << nl
                        << "parcelTrackSubStepCount=" << parcelTrackSubStepCount << nl
                        << "cell=" << cell() << nl
                        << "face=" << faceI << nl
                        << "patch=" << patchName << nl
                        << "stepFraction=" << stepFraction() << nl
                        << "remainingTrackTime=" << tEnd << nl
                        << "cellDeltaT=" << baseDt << nl
                        << "position=" << position() << nl
                        << "U=" << U_ << nl
                        << "Utracking=" << Utracking << nl
                        << "inRotatingRegion=" << inRotatingRegion << nl
                        << "newParcel=" << newParcel() << nl
                        << "classification=" << classification() << nl
                        << "lastTrackFraction=" << lastTrackFraction << nl
                        << "lastAttemptDt=" << lastAttemptDt << nl
                        << "lastEffectiveDt=" << lastEffectiveDt << nl
                        << "lastRemainingTrackTimeBeforeTrack="
                        << lastRemainingTrackTimeBeforeTrack << nl
                        << "lastStartCell=" << lastStartCell << nl
                        << "lastEndCell=" << lastEndCell << nl
                        << "lastFace=" << lastFace << nl
                        << "lastEndedOnBoundaryFace="
                        << lastEndedOnBoundaryFace << nl
                        << "lastStartPosition=" << lastStartPosition << nl
                        << "lastRequestedEndPosition="
                        << lastRequestedEndPosition << nl
                        << "lastEndPosition=" << lastEndPosition
                        << exit(FatalError);
                }
            }

            if (td.cloud().coordSystem().type() == "dsmcCartesian")
            {
                // Apply correction to position for reduced-D cases,
                // but not for axisymmetric cases
                meshTools::constrainToMeshCentre(mesh, position());

                // Apply correction to velocity to constrain tracking for
                // reduced-D cases,  but not for axisymmetric cases
                meshTools::constrainDirection(mesh, mesh.solutionD(), Utracking);
            }

            const point startPos(position());
            const vector startU(U_);
            const scalar requestedDt = dt;
            const point requestedEndPosition = position() + requestedDt*Utracking;
            const label startCell = cell();

            const scalar trackFraction =
                trackToFace(requestedEndPosition, td, true);

            lastTrackFraction = trackFraction;
            lastAttemptDt = requestedDt;
            lastEffectiveDt = requestedDt*trackFraction;
            lastRemainingTrackTimeBeforeTrack = tEnd;
            lastStartPosition = startPos;
            lastRequestedEndPosition = requestedEndPosition;
            lastEndPosition = position();
            lastStartCell = startCell;
            lastEndCell = cell();
            lastFace = face();
            lastEndedOnBoundaryFace =
            (
                face() >= 0
             && pbMesh.whichPatch(face()) != -1
            );

            if (trackFraction <= particle::minStepFractionTol)
            {
                // trackToFace() 返回极小步长有两种截然不同的语义：
                // 1. 纯几何/容差病态：粒子几乎没推进，且也没有真正以边界事件结束；
                // 2. 合法边界命中：patch interaction 已经执行，但这次命中的
                //    几何推进量非常小。
                //
                // 对第 2 类情况，wall / cyclic / stage-interface 可能已经在
                // hitPatch()/hitWallPatch() 内写回了新的 U_。若这里再沿用旧的
                // 旋转区预测速度：
                //     U_new = U_start + dt * a0
                // 就会把边界模型刚写入的速度覆盖掉。
                //
                // 因此这里必须先区分“是否已经以 boundary face 结束”。
                const bool endedOnBoundaryFace =
                (
                    face() >= 0
                 && pbMesh.whichPatch(face()) != -1
                );

                if (endedOnBoundaryFace)
                {
                    // 对已经命中边界的极小步长事件：
                    // 1. 保留边界模型刚写入的 U_；
                    // 2. 不再把整段剩余 tEnd 直接吞掉；
                    // 3. 至少按 minStepFractionTol*dt 消耗一个极小但正的时间份额，
                    //    让当前时间步还能继续推进剩余部分。
                    //
                    // 这样既避免“边界后速度被改回旧值”，也避免把大段剩余时间
                    // 在一次零长度/近零长度命中后直接丢弃。
                    const scalar consumedDt =
                        max(trackFraction*dt, particle::minStepFractionTol*dt);

                    // 正常子步路径会在后面的统一出口调用
                    // tracker().trackParcelFaceTransition(*this) 记录一次面跃迁。
                    // 但极小步长边界命中会在这里直接 continue，若不补记，
                    // 当前这次合法的 patch/face 穿越就会完全绕过 faceTracker，
                    // 使依赖 tracker() 的面通量与事件计数系统性偏小。
                    //
                    // 这里复用与常规路径相同的记账条件：只要 face() 有效，
                    // 就按照当前 parcel 已经被边界模型更新后的状态记一次跃迁。
                    if (face() != -1)
                    {
                        td.cloud().tracker().trackParcelFaceTransition(*this);
                    }

                    tEnd = max(0.0, tEnd - consumedDt);
                    stepFraction() = 1.0 - tEnd/totalTrackTime;
                    continue;
                }

                if (inRotatingRegion)
                {
                    const vector a0 =
                        td.cloud().rotatingFrameAcceleration
                        (
                            startPos,
                            startU
                        );

                    U_ = startU + dt*a0;
                }

                // 若既没有完成有效边界交互、又只得到近零推进量，则仍保留
                // 现有的保护性早退策略，把这次病态子步直接终止。
                tEnd = 0.0;
                stepFraction() = 1.0;
                break;
            }

            dt *= trackFraction;

            if (inRotatingRegion)
            {
                const vector a0 =
                    td.cloud().rotatingFrameAcceleration
                    (
                        startPos,
                        startU
                    );
                const vector Uhalf = startU + 0.5*dt*a0;
                const point xhalf = startPos + 0.5*dt*Uhalf;
                const vector ahalf =
                    td.cloud().rotatingFrameAcceleration
                    (
                        xhalf,
                        Uhalf
                    );

                // 若本子步已经以 boundary face 结束，则 patch interaction
                // 可能已经修改了 U_。此时不能再用 SRF 子步积分结果覆盖它，
                // 否则 wall 反射或 stage-interface 变换后的速度会被写回旧值。
                const bool endedOnBoundaryFace =
                (
                    face() >= 0
                 && pbMesh.whichPatch(face()) != -1
                );

                if (!endedOnBoundaryFace)
                {
                    U_ = startU + dt*ahalf;
                }
            }

            tEnd -= dt;

            stepFraction() = 1.0 - tEnd/totalTrackTime; // NEW VINCENT
            //stepFraction() = 1.0 - tEnd/trackTime; // OLD FORMULATION

            //- face tracking info
            if (face() != -1)
            {
                //- measure flux properties
                td.cloud().tracker().trackParcelFaceTransition(*this);
            }

            if (onBoundary() && td.keepParticle)
            {
                if (isA<processorPolyPatch>(pbMesh[patch(face())]))
                {
                    td.switchProcessor = true;
                }

            }
        }
    }
    else
    {
        //- The stuck particle is considered for desorption
        //  NOTE: there should be a better way to locate the patch than looping through them all
        forAll(td.cloud().boundaries().patchBoundaryModels(), c)
        {
            if (td.cloud().boundaries().patchBoundaryModels()[c]->patchId() == stuck().wallTemperature()[1])
            {
                // then this patch is the "dsmc*Sticking*WallPatch" the particle is stuck on
                td.cloud().boundaries().patchBoundaryModels()[c]->controlParticle(*this, td);
                break;
            }
        }
    }

    return td.keepParticle;
}


bool Foam::dsmcParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td,
    const label,
    const scalar,
    const tetIndices&
)
{
    // 只有物理 cyclic / cyclicAMI patch 需要在这里被 DSMC 周期边界模型提前接管。
    // processorCyclic 是并行分解产生的处理器通信边界，粒子跨 rank 时必须保留
    // OpenFOAM 标准 processor transfer 语义，由 hitProcessorPatch() 设置
    // switchProcessor。若在这里拦截，会把 processor patch-local face 当成物理
    // 周期面处理，混淆本地/远程 rank 的 face 与 cell 语义，容易导致并行追踪卡死。
    if
    (
        !isA<cyclicPolyPatch>(pp)
     && !isA<cyclicAMIPolyPatch>(pp)
    )
    {
        return false;
    }

    const label hitFace = face();

    forAll(td.cloud().boundaries().cyclicBoundaryModels(), c)
    {
        const labelList& faces =
            td.cloud().boundaries().cyclicBoundaryModels()[c]->allFaces();

        if (findIndex(faces, hitFace) != -1)
        {
            const bool traceThisBoundaryEvent =
                boundaryTraceEnabled(td)
             && boundaryTracePatchAllowed(td, pp.name());
            const bool traceThisBoundaryLocation =
                traceThisBoundaryEvent
             && boundaryTraceFaceAllowed(td, face())
             && boundaryTraceCellAllowed(td, cell());
            const label traceEventId =
                traceThisBoundaryLocation ? nextBoundaryTraceEventId() : -1;

            recordBoundaryEventEntry(td, pp.name(), "cyclic");

            if (traceThisBoundaryLocation)
            {
                writeBoundaryTrace
                (
                    *this,
                    td,
                    traceEventId,
                    "beforeCyclicBoundaryModelControl",
                    pp.name(),
                    pp.type(),
                    c
                );
            }

            td.cloud().boundaries().cyclicBoundaryModels()[c]->controlMol(*this, td);
            recordBoundaryEventExit(td);

            if (traceThisBoundaryLocation)
            {
                writeBoundaryTrace
                (
                    *this,
                    td,
                    traceEventId,
                    "afterCyclicBoundaryModelControl",
                    pp.name(),
                    pp.type(),
                    c
                );
            }

            return true;
        }
    }

    // 若该 cyclic patch 没有注册 DSMC 周期模型，则退回 OpenFOAM 默认处理。
    return false;
}


void Foam::dsmcParcel::hitCyclicPatch
(
    const cyclicPolyPatch& cpp,
    trackingData& td
)
{
    // 只有在 hitPatch(...) 没有截获该 cyclic patch 时，才会走到这里。
    // 此时直接复用 OpenFOAM 标准 cyclic 行为即可。
    particle::hitCyclicPatch(cpp, td);
}


void Foam::dsmcParcel::hitCyclicAMIPatch
(
    const cyclicAMIPolyPatch& cpp,
    trackingData& td,
    const vector& direction
)
{
    particle::hitCyclicAMIPatch(cpp, td, direction);
}


void Foam::dsmcParcel::hitProcessorPatch
(
    const processorPolyPatch&,
    trackingData& td
)
{
    td.switchProcessor = true;
}


void Foam::dsmcParcel::hitWallPatch
(
    const wallPolyPatch& wpp,
    trackingData& td,
    const tetIndices& tetIs
)
{
    //- find which patch has been hit
    label patchIndex = wpp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];
    const bool traceThisBoundaryEvent =
        boundaryTraceEnabled(td)
     && boundaryTracePatchAllowed(td, wpp.name());
    const bool traceThisBoundaryLocation =
        traceThisBoundaryEvent
     && boundaryTraceFaceAllowed(td, face())
     && boundaryTraceCellAllowed(td, cell());
    const label traceEventId =
        traceThisBoundaryLocation ? nextBoundaryTraceEventId() : -1;

    recordBoundaryEventEntry(td, wpp.name(), "wall");

    if (traceThisBoundaryLocation)
    {
        writeBoundaryTrace
        (
            *this,
            td,
            traceEventId,
            "beforeWallPatchModelControl",
            wpp.name(),
            wpp.type(),
            patchModelId
        );
    }

    //- apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);

    if (traceThisBoundaryLocation)
    {
        writeBoundaryTrace
        (
            *this,
            td,
            traceEventId,
            "afterWallPatchModelControl",
            wpp.name(),
            wpp.type(),
            patchModelId
        );
    }

    // 壁面模型返回后，对仍保持自由态的粒子统一追加一次极小 inward push，
    // 防止它几何上仍贴在 wall face 上，导致同一 DSMC 步里再次命中同一壁面。
    applyPostWallInwardPush(*this, td);
    recordBoundaryEventExit(td);

    if (traceThisBoundaryLocation)
    {
        writeBoundaryTrace
        (
            *this,
            td,
            traceEventId,
            "afterWallInwardPush",
            wpp.name(),
            wpp.type(),
            patchModelId
        );
    }
}


void Foam::dsmcParcel::hitPatch
(
    const polyPatch& pp,
    trackingData& td
)
{
    //- find which patch has been hit
    label patchIndex = pp.index();

    const label& patchModelId = td.cloud().boundaries().patchToModelIds()[patchIndex];
    const bool traceThisBoundaryEvent =
        boundaryTraceEnabled(td)
     && boundaryTracePatchAllowed(td, pp.name());
    const bool traceThisBoundaryLocation =
        traceThisBoundaryEvent
     && boundaryTraceFaceAllowed(td, face())
     && boundaryTraceCellAllowed(td, cell());
    const label traceEventId =
        traceThisBoundaryLocation ? nextBoundaryTraceEventId() : -1;

    recordBoundaryEventEntry(td, pp.name(), "generic");

    if (traceThisBoundaryLocation)
    {
        writeBoundaryTrace
        (
            *this,
            td,
            traceEventId,
            "beforeGenericPatchModelControl",
            pp.name(),
            pp.type(),
            patchModelId
        );
    }

    //- apply a boundary model when a molecule collides with this poly patch
    td.cloud().boundaries().patchBoundaryModels()[patchModelId]->controlParticle(*this, td);
    recordBoundaryEventExit(td);

    if (traceThisBoundaryLocation)
    {
        writeBoundaryTrace
        (
            *this,
            td,
            traceEventId,
            "afterGenericPatchModelControl",
            pp.name(),
            pp.type(),
            patchModelId
        );
    }
}


void Foam::dsmcParcel::transformProperties
(
    const tensor& T
)
{
   particle::transformProperties(T);
   U_ = transform(T, U_);
}


void Foam::dsmcParcel::transformProperties
(
    const vector& separation
)
{
    particle::transformProperties(separation);
}

bool Foam::dsmcParcel::relocateStuckParcel
(
    const polyMesh& mesh
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // find the closest patch and patch face indices, as this is the patch/face
    // on which this parcel is stuck.
    scalar closestFaceDistance = GREAT;
    label closestPatchi = -1;
    label closestPatchFacei = -1;

    forAll(mesh.cells()[cell()], i)
    {
        const label facei = mesh.cells()[cell()][i];

        // find corresponding boundary patch
        const label patchi = bMesh.whichPatch(facei);

        if (patchi != -1)
        {
            // this is a boundary patch, i.e. a potential candidate for the
            // patch to which parcel is stuck.
            const polyPatch& wpp = bMesh[patchi];

            // patch face index:
            const label patchFacei = wpp.whichFace(facei);

            // calculate distance between parcel position and this patch face:
            pointHit pHit
            (
                wpp[patchFacei].nearestPoint(position(), wpp.points())
            );

            if (pHit.hit())
            {
                const scalar distance = pHit.distance();
                if (distance < closestFaceDistance)
                {
                    // found new closest patch face
                    closestFaceDistance = distance;
                    closestPatchi = patchi;
                    closestPatchFacei = patchFacei;
                }
            }
        }
    }

    // did we find a closest patch face?
    if (closestPatchi != -1)
    {
        // yes, reset stuck parcel information to this patch face
        stuck().wallTemperature()[1] = closestPatchi;
        stuck().wallTemperature()[2] = closestPatchFacei;
        return true;
    }
    // no, indicates fatal error
    return false;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "dsmcParcelIO.C"


// ************************************************************************* //
