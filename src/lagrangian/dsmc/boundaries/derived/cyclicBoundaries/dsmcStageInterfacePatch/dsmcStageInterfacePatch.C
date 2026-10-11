/*---------------------------------------------------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\\\  /    A nd           |
     \\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/

#include "dsmcStageInterfacePatch.H"
#include "dsmcCloud.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(dsmcStageInterfacePatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcStageInterfacePatch,
    dictionary
);


void dsmcStageInterfacePatch::readProperties()
{
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0));
    rotationCentre_ =
        propsDict_.lookupOrDefault<vector>("rotationCentre", vector::zero);
    omega_ = propsDict_.lookupOrDefault<vector>("omega", vector::zero);
    sourceFrame_ =
        propsDict_.lookupOrDefault<word>("sourceFrame", word("absolute"));
    targetFrame_ =
        propsDict_.lookupOrDefault<word>("targetFrame", word("absolute"));
    transformPosition_ =
        propsDict_.lookupOrDefault<bool>("transformPosition", true);
    approximateFaceCentreInsertion_ =
        propsDict_.lookupOrDefault<bool>("approximateFaceCentreInsertion", true);
    approximateInsertionFraction_ =
        propsDict_.lookupOrDefault<scalar>("approximateInsertionFraction", 2.0e-2);
    transformVelocity_ =
        propsDict_.lookupOrDefault<bool>("transformVelocity", false);
    debug_ = propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);

    const scalar axisMag = mag(rotationAxis_);

    // rotationAxis_ 必须合法，后续位置/速度变换都依赖该轴定义。
    if (axisMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_
            << exit(FatalError);
    }

    rotationAxis_ /= axisMag;

    // 近似插入深度比例必须非负，否则会把粒子朝 patch 外侧推送。
    if (approximateInsertionFraction_ < 0.0)
    {
        FatalErrorInFunction
            << "approximateInsertionFraction must be non-negative for patch "
            << patchName_
            << ", but got " << approximateInsertionFraction_
            << exit(FatalError);
    }

    // frame 只允许 absolute / relative。
    if (!isValidFrame(sourceFrame_))
    {
        FatalErrorInFunction
            << "Invalid sourceFrame '" << sourceFrame_
            << "' for patch " << patchName_
            << ". Valid values are absolute or relative."
            << exit(FatalError);
    }

    if (!isValidFrame(targetFrame_))
    {
        FatalErrorInFunction
            << "Invalid targetFrame '" << targetFrame_
            << "' for patch " << patchName_
            << ". Valid values are absolute or relative."
            << exit(FatalError);
    }

    // 速度变换请求必须与 cloud 端 SRF 配置自洽。
    validateVelocityTransformConfiguration();

    // 当前实现只服务于严格匹配的 cyclic/cyclic interface。
    // 若输入的是 cyclicAMI，则会与本实现的 transformPosition/parcel 交换语义混淆，
    // 因此这里显式拒绝。
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
            << "contains cyclicAMI topology, but the current MCOR1_0 minimal "
            << "dsmcStageInterfacePatch is only validated for strictly matched "
            << "cyclic/cyclic interfaces." << nl
            << "Please keep MCOR1_0/2.1.b on a cyclic interface first, and "
            << "only introduce cyclicAMI after the minimal data-exchange path is "
            << "fully validated."
            << exit(FatalError);
    }

    if
    (
        !isA<cyclicPolyPatch>(sourcePatch)
     || !isA<cyclicPolyPatch>(targetPatch)
    )
    {
        FatalErrorInFunction
            << "Patch pair (" << patchName_ << ", " << neighbPatchName_ << ") "
            << "must be a cyclic/cyclic pair for the current MCOR1_0 minimal "
            << "stage-interface implementation." << nl
            << "Source patch type: " << sourcePatch.type()
            << ", target patch type: " << targetPatch.type()
            << exit(FatalError);
    }
}


bool dsmcStageInterfacePatch::isValidFrame(const word& frame) const
{
    return frame == "absolute" || frame == "relative";
}


bool dsmcStageInterfacePatch::isRelativeFrame(const word& frame) const
{
    return frame == "relative";
}


void dsmcStageInterfacePatch::validateVelocityTransformConfiguration() const
{
    // transformVelocity_ 关闭时，当前 patch 只承担 topology transfer，
    // 不需要检查 SRF 物理配置。
    if (!transformVelocity_)
    {
        return;
    }

    // 这里不能直接用 cloud_.rotatingReferenceFrameActive()、
    // cloud_.rotatingFrameOmega() 和 cloud_.rotatingFrameOrigin() 做校验。
    // 原因是当前校验发生在 dsmcStageInterfacePatch 构造阶段，而此时
    // dsmcCloud 还没走到 initialiseRotatingReferenceFrame()，
    // 这些缓存状态仍保持默认值。
    //
    // 因此这里必须直接读取 cloud 已经持有的原始 dsmcProperties 字典，
    // 以避免“字典其实是 active=true，但缓存尚未初始化，结果误报 false”。
    const IOdictionary& particleProperties = cloud_.particleProperties();
    const bool hasRotatingReferenceFrameDict =
        particleProperties.found("rotatingReferenceFrameProperties");

    const dictionary* rotatingDictPtr = nullptr;
    bool rotatingReferenceFrameEnabled = false;
    vector rotatingFrameOmega(vector::zero);
    vector rotatingFrameOrigin(vector::zero);
    word rotatingRegionZoneName(word::null);

    if (hasRotatingReferenceFrameDict)
    {
        rotatingDictPtr =
            &particleProperties.subDict("rotatingReferenceFrameProperties");
        rotatingReferenceFrameEnabled =
            rotatingDictPtr->lookupOrDefault<bool>("active", false);
        rotatingDictPtr->readIfPresent("cellZone", rotatingRegionZoneName);
        rotatingFrameOmega =
            rotatingDictPtr->lookupOrDefault<vector>("omega", vector::zero);
        rotatingFrameOrigin =
            rotatingDictPtr->lookupOrDefault<vector>("origin", vector::zero);
    }

    // 若字典要求在 absolute / relative 之间切换，但原始 dsmcProperties
    // 中并未启用 rotatingReferenceFrame，则边界和域内对速度的解释会立刻失配：
    // - 边界把 U 改成相对速度；
    // - 域内却仍按绝对速度推进粒子。
    if
    (
        sourceFrame_ != targetFrame_
     && !rotatingReferenceFrameEnabled
    )
    {
        FatalErrorInFunction
            << "transformVelocity=true requests a reference-frame change on "
            << "patch " << patchName_ << ", but "
            << "rotatingReferenceFrameProperties.active is false in "
            << "constant/dsmcProperties." << nl
            << "Stage-interface velocity conversion and in-domain SRF "
            << "tracking must be enabled together."
            << exit(FatalError);
    }

    if
    (
        sourceFrame_ != targetFrame_
     && rotatingRegionZoneName.empty()
    )
    {
        FatalErrorInFunction
            << "transformVelocity=true requests an absolute/relative frame "
            << "change on patch " << patchName_ << ", but "
            << "rotatingReferenceFrameProperties.cellZone is not specified in "
            << "constant/dsmcProperties." << nl
            << "For MCOR1_0 multi-domain rotor/stator cases, cellZone must be "
            << "explicit so that only the rotor-side domain is treated as the "
            << "rotating region."
            << exit(FatalError);
    }

    // 一旦启用速度变换，stage-interface 使用的角速度必须和 cloud 端
    // rotatingReferenceFrame 的角速度完全一致。
    if
    (
        rotatingReferenceFrameEnabled
     && mag(omega_ - rotatingFrameOmega) > SMALL
    )
    {
        FatalErrorInFunction
            << "Stage-interface omega " << omega_
            << " is inconsistent with rotatingReferenceFrame omega "
            << rotatingFrameOmega
            << " on patch " << patchName_
            << exit(FatalError);
    }

    // 同理，旋转中心也必须一致。
    // 否则界面上计算的 Omega x r 与域内 SRF 追踪使用的参考中心不同。
    if
    (
        rotatingReferenceFrameEnabled
     && mag(rotationCentre_ - rotatingFrameOrigin) > SMALL
    )
    {
        FatalErrorInFunction
            << "Stage-interface rotationCentre " << rotationCentre_
            << " is inconsistent with rotatingReferenceFrame origin "
            << rotatingFrameOrigin
            << " on patch " << patchName_
            << exit(FatalError);
    }

    // 配置级 SRF 参数一致以后，还要继续确认 patch 两侧 owner cell 的
    // zone 归属确实与 sourceFrame_/targetFrame_ 一致。
    // 这一步专门拦截“字典写成 relative->absolute，但 patch 实际贴在
    // stator/stator 或 rotor/rotor 两侧”这类几何-物理语义错位问题。
    if (sourceFrame_ != targetFrame_)
    {
        validatePatchFrameConfiguration(rotatingRegionZoneName);
    }
}


void dsmcStageInterfacePatch::validatePatchFrameConfiguration
(
    const word& rotatingRegionZoneName
) const
{
    if (sourceFrame_ == targetFrame_)
    {
        return;
    }

    // 走到这里时，cellZone 缺失已经在上一层被拦截。
    // 这里仍保留一道显式保护，防止后续维护时有人绕过前面的判断。
    if (rotatingRegionZoneName.empty())
    {
        FatalErrorInFunction
            << "Internal error: rotatingRegionZoneName is empty while "
            << "validating patch/frame consistency on patch " << patchName_
            << exit(FatalError);
    }

    const cellZoneMesh& cellZones = mesh_.cellZones();
    const label zoneId = cellZones.findZoneID(rotatingRegionZoneName);

    // 这里同样不能容忍“写了 cellZone 但 mesh 里没有”。
    // 因为一旦 zone 名称拼错，sourceFrame_/targetFrame_ 的物理含义
    // 就不再可判定，继续运行只会把错误推迟到更难排查的阶段。
    if (zoneId == -1)
    {
        FatalErrorInFunction
            << "Rotating-reference-frame cellZone '" << rotatingRegionZoneName
            << "' was not found while validating stage-interface patch "
            << patchName_ << "." << nl
            << "MCOR1_0 requires the interface patch sides and the SRF "
            << "cellZone definition to be mutually consistent."
            << exit(FatalError);
    }

    // 这里显式重建一个“配置级旋转区掩码”。
    // 不能依赖 cloud_.isInRotatingRegion() 的原因是：
    // dsmcStageInterfacePatch 构造阶段早于 cloud 的 SRF 缓存初始化，
    // 此时直接查 cloud 缓存会把所有 cell 都误读成 false。
    boolList rotatingRegionCellMask(mesh_.nCells(), false);
    const labelList& zoneCells = cellZones[zoneId];

    forAll(zoneCells, i)
    {
        const label celli = zoneCells[i];

        // zone 内若出现越界 cell 编号，说明 mesh / topoSet 数据本身已损坏。
        // 这种情况绝不能继续跑，因为后续 frame 校验会建立在非法拓扑上。
        if (celli < 0 || celli >= rotatingRegionCellMask.size())
        {
            FatalErrorInFunction
                << "Cell " << celli << " from rotatingReferenceFrame cellZone '"
                << rotatingRegionZoneName << "' is outside the mesh cell range "
                << "[0, " << rotatingRegionCellMask.size() - 1 << "]."
                << exit(FatalError);
        }

        rotatingRegionCellMask[celli] = true;
    }

    // patchId_ 这一侧按 sourceFrame_ 解释；
    // neighbPatchId_ 这一侧按 targetFrame_ 解释。
    // 这与 controlMol() 中 A 面命中时 sourceFrame_ -> targetFrame_ 的
    // 方向约定保持完全一致。
    validatePatchOwnerCellsAgainstFrame
    (
        mesh_.boundaryMesh()[patchId_],
        sourceFrame_,
        rotatingRegionCellMask
    );

    validatePatchOwnerCellsAgainstFrame
    (
        mesh_.boundaryMesh()[neighbPatchId_],
        targetFrame_,
        rotatingRegionCellMask
    );
}


void dsmcStageInterfacePatch::validatePatchOwnerCellsAgainstFrame
(
    const polyPatch& patch,
    const word& frame,
    const boolList& rotatingRegionCellMask
) const
{
    const bool expectedInRotatingRegion = isRelativeFrame(frame);
    const labelList& patchCells = patch.faceCells();

    forAll(patchCells, patchFaceI)
    {
        const label celli = patchCells[patchFaceI];

        // patch owner cell 本身必须先是一个合法 mesh cell。
        // 若这里越界，说明 patch 拓扑已经不可用，后续无论 frame 写什么都无意义。
        if (celli < 0 || celli >= rotatingRegionCellMask.size())
        {
            FatalErrorInFunction
                << "Patch " << patch.name()
                << " references invalid owner cell " << celli
                << " while validating stage-interface frame consistency."
                << exit(FatalError);
        }

        const bool actualInRotatingRegion = rotatingRegionCellMask[celli];

        if (actualInRotatingRegion != expectedInRotatingRegion)
        {
            FatalErrorInFunction
                << "Stage-interface patch/frame mismatch detected on patch "
                << patch.name() << "." << nl
                << "frame=" << frame
                << ", expectedOwnerCellInRotatingRegion="
                << expectedInRotatingRegion << nl
                << "Offending mesh face=" << patch.start() + patchFaceI
                << ", ownerCell=" << celli
                << ", actualOwnerCellInRotatingRegion="
                << actualInRotatingRegion << nl
                << "This usually means the patch pairing, sourceFrame/"
                << "targetFrame definition, or rotatingReferenceFrame cellZone "
                << "assignment is inconsistent."
                << exit(FatalError);
        }
    }
}


void dsmcStageInterfacePatch::validateParcelFrameConsistency
(
    const label oldCell,
    const label newCell,
    const word& fromFrame,
    const word& toFrame
) const
{
    // 这层校验只在真正启用速度变换时生效。
    if (!transformVelocity_)
    {
        return;
    }

    const bool oldCellInRotatingRegion = cloud_.isInRotatingRegion(oldCell);
    const bool newCellInRotatingRegion = cloud_.isInRotatingRegion(newCell);
    const bool expectedOldCellInRotatingRegion = isRelativeFrame(fromFrame);
    const bool expectedNewCellInRotatingRegion = isRelativeFrame(toFrame);

    // oldCell 必须和 fromFrame 对应的物理语义一致。
    if (oldCellInRotatingRegion != expectedOldCellInRotatingRegion)
    {
        FatalErrorInFunction
            << "Source-cell/frame mismatch detected on patch " << patchName_
            << "." << nl
            << "oldCell=" << oldCell
            << ", fromFrame=" << fromFrame
            << ", oldCellInRotatingRegion=" << oldCellInRotatingRegion
            << ", expectedOldCellInRotatingRegion="
            << expectedOldCellInRotatingRegion
            << exit(FatalError);
    }

    // newCell 必须和 toFrame 对应的物理语义一致。
    if (newCellInRotatingRegion != expectedNewCellInRotatingRegion)
    {
        FatalErrorInFunction
            << "Target-cell/frame mismatch detected on patch " << patchName_
            << "." << nl
            << "newCell=" << newCell
            << ", toFrame=" << toFrame
            << ", newCellInRotatingRegion=" << newCellInRotatingRegion
            << ", expectedNewCellInRotatingRegion="
            << expectedNewCellInRotatingRegion
            << exit(FatalError);
    }
}


vector dsmcStageInterfacePatch::transformedVelocity
(
    const vector& U,
    const point& position,
    const word& fromFrame,
    const word& toFrame
) const
{
    // 如果源参考系与目标参考系相同，则不需要做任何速度修正。
    if (fromFrame == toFrame)
    {
        return U;
    }

    // Stage Interface 的速度转换使用粒子真正穿越界面的空间位置，
    // 而不是 cell centre 或预设半径。
    const vector r = position - rotationCentre_;
    const vector uRot = omega_ ^ r;

    // absolute -> relative：减去局部牵连速度。
    if (fromFrame == "absolute" && toFrame == "relative")
    {
        return U - uRot;
    }

    // relative -> absolute：加回局部牵连速度。
    if (fromFrame == "relative" && toFrame == "absolute")
    {
        return U + uRot;
    }

    // 正常流程不会到这里；若到了，说明 frame 组合超出了当前实现能力。
    FatalErrorInFunction
        << "Unsupported frame transform from " << fromFrame
        << " to " << toFrame << exit(FatalError);

    return U;
}


point dsmcStageInterfacePatch::transferParcel(dsmcParcel& p) const
{
    // 这一函数只负责“合法跨面 + 合法几何插入”。
    // 速度参考系变换不在这里做，但这里必须把接收侧真实 crossing point 的
    // 几何映射位置返回给上层。原因是：
    // 1. 速度变换需要严格使用接收侧真实穿越点：
    //        U_new = U_old +/- Omega x r
    // 2. 数值稳健性的近似插入可能把粒子人为放到目标面中心附近；
    // 3. 因此“物理取点”和“数值插入点”必须分离。
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
        FatalErrorInFunction
            << "Patch " << patch.name()
            << " is not cyclic, but dsmcStageInterfacePatch::transferParcel() "
            << "was called."
            << exit(FatalError);

        return positionBeforeTransfer;
    }

    const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(patch);

    // 先完成严格的一一对应拓扑迁移：
    //     targetFace = transformGlobalFace(sourceFace)
    // 这一步不允许近似，确保 face/cell 归属始终正确。
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

    // 先把真实 crossing point 精确映射到接收侧 patch。
    // 上层速度变换必须使用这个 mappedCrossPosition，而不是最终插入点。
    point mappedCrossPosition = positionBeforeTransfer;
    receiveCpp.transformPosition(mappedCrossPosition, patchFaceI);

    const point& targetFaceCentre = mesh_.faceCentres()[p.face()];
    const point& ownerCellCentre = mesh_.cellCentres()[p.cell()];

    if (transformPosition_)
    {
        if (approximateFaceCentreInsertion_)
        {
            // 近似插入只用于几何稳健性：
            // 1. 从目标面中心出发；
            // 2. 沿 owner cell 内部方向推进一个与面尺度成比例的小距离；
            // 3. 不再让这个近似位置污染速度参考系变换取点。
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
            // 若用户关闭近似插入，则最终几何位置直接采用精确映射后的 crossing point。
            p.position() = mappedCrossPosition;
        }

        // 只有在 transformPosition=true 时，才允许再做这层 inward push。
        // 它纯粹服务于防止粒子继续贴在接收侧 patch 上。
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

        // 一旦跨面后重新写了几何位置，就必须重建 tetFace/tetPt。
        // 否则粒子虽然已在接收侧单元内部，但局部 tet 拓扑仍可能保留为
        // 命中 interface 时的旧状态，后续 trackToFace() 会从不一致状态起步。
        if (mag(p.position() - positionBeforeTransfer) > ROOTVSMALL)
        {
            p.initCellFacePt();
        }
    }

    return mappedCrossPosition;
}


dsmcStageInterfacePatch::dsmcStageInterfacePatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcCyclicBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    rotationAxis_(vector::zero),
    rotationCentre_(vector::zero),
    omega_(vector::zero),
    sourceFrame_("absolute"),
    targetFrame_("absolute"),
    transformPosition_(true),
    approximateFaceCentreInsertion_(true),
    approximateInsertionFraction_(2.0e-2),
    transformVelocity_(false),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCountA_(0),
    debugPrintCountB_(0),
    transferCountAtoB_(0),
    transferCountBtoA_(0),
    transformCountRelToAbs_(0),
    transformCountAbsToRel_(0),
    cumulativeFinalPositionVelocityDifferenceRelToAbs_(0.0),
    cumulativeFinalPositionVelocityDifferenceAbsToRel_(0.0),
    cumulativeFinalPositionDeltaUDifferenceRelToAbs_(0.0),
    cumulativeFinalPositionDeltaUDifferenceAbsToRel_(0.0),
    maxFinalPositionVelocityDifferenceRelToAbs_(0.0),
    maxFinalPositionVelocityDifferenceAbsToRel_(0.0),
    maxFinalPositionDeltaUDifferenceRelToAbs_(0.0),
    maxFinalPositionDeltaUDifferenceAbsToRel_(0.0),
    lastSummaryTimeIndex_(-1)
{
    // 当前边界不写额外统计文件，只在调试阶段通过 Info 输出检查跨界面行为。
    writeInTimeDir_ = false;
    writeInCase_ = false;
    readProperties();
}


dsmcStageInterfacePatch::~dsmcStageInterfacePatch()
{}


void dsmcStageInterfacePatch::calculateProperties()
{}


void dsmcStageInterfacePatch::initialConfiguration()
{}


void dsmcStageInterfacePatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    (void)td;

    // 粒子已命中本 stage-interface 管理的一对 cyclic faces 之一。
    // 这里依次完成侧别识别、跨面拓扑迁移和可选的速度参考系变换。

    const label faceI = p.face();
    const label oldPatchI = mesh_.boundaryMesh().whichPatch(faceI);
    const point crossPosition = p.position();
    const point oldPosition = crossPosition;
    const vector oldVelocity = p.U();
    // 严格匹配的 stage-interface 只代表两个几何/参考系区域之间的拓扑迁移。
    // 对固定网格权重自适应而言，parcel.RWF() 是该 macro parcel 的统计权重因子：
    //
    //     N_real(parcel) = nParticle * cellWeight * parcel.RWF()
    //
    // 这里没有发生粒子分裂、合并，也没有发生 AMI 面积比例重映射；
    // 因此穿越接口前后的 RWF 必须完全保持不变。
    // 若这个值被误改，质量、动量和能量统计都会按同一个错误比例漂移。
    const scalar oldRadialWeight = p.RWF();
    const label oldCell = p.cell();
    const label oldTetFace = p.tetFace();
    const label oldTetPt = p.tetPt();
    const bool oldCellInRotatingRegion = cloud_.isInRotatingRegion(oldCell);
    word side("none");
    bool fromSideA = false;
    word fromFrame(sourceFrame_);
    word toFrame(targetFrame_);
    label expectedNewPatchId = -1;

    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    // 命中 A 面时，按 sourceFrame_ -> targetFrame_ 解释。
    if (sourceIndexA != -1)
    {
        side = "A";
        fromSideA = true;
        fromFrame = sourceFrame_;
        toFrame = targetFrame_;
        expectedNewPatchId = neighbPatchId_;
    }
    // 命中 B 面时，方向自动反过来，变成 targetFrame_ -> sourceFrame_。
    else if (sourceIndexB != -1)
    {
        side = "B";
        fromSideA = false;
        fromFrame = targetFrame_;
        toFrame = sourceFrame_;
        expectedNewPatchId = patchId_;
    }
    else
    {
        // face 不属于当前接口模型管理范围时，直接返回。
        return;
    }

    // 先执行跨 domain 交换，并拿到接收侧真实 crossing point 的映射位置。
    // 后续速度参考系变换必须基于这个 mappedCrossPosition：
    //     U_new = U_old +/- Omega x r(mappedCrossPosition)
    // 而不能再退回到源侧旧的 crossPosition，也不能偷换成最终近似插入点。
    const point mappedCrossPosition = transferParcel(p);
    const scalar newRadialWeight = p.RWF();

    // 固定网格权重守恒断言：
    // 1. 对 dsmcStageInterfacePatch 这种严格匹配 cyclic interface，
    //    几何映射是一一对应，真实分子代表数不应随跨面改变；
    // 2. transformVelocity_ 只改变速度参考系表达，不改变 parcel 的统计权重；
    // 3. 如果未来某段代码在 transferParcel() 或其它接口处理里意外改写 RWF，
    //    这里必须立即报错，而不是让守恒误差传播到宏观场采样。
    //
    // 容差采用相对 SMALL，只用于屏蔽浮点往返级别的不可见误差。
    // 正常路径下该差值应为严格的 0，因为当前 transferParcel() 不写 RWF。
    if
    (
        mag(newRadialWeight - oldRadialWeight)
      > SMALL*max(scalar(1), mag(oldRadialWeight))
    )
    {
        FatalErrorInFunction
            << "Stage-interface transfer changed parcel radialWeight on a "
            << "strictly matched cyclic interface." << nl
            << "Interface patch pair: (" << patchName_ << ", "
            << neighbPatchName_ << ")" << nl
            << "Transfer face: " << faceI << nl
            << "Old radialWeight: " << oldRadialWeight << nl
            << "New radialWeight: " << newRadialWeight << nl
            << "This would violate fixed-mesh statistical-weight "
            << "conservation because no split, merge, or AMI remap is active."
            << exit(FatalError);
    }

    // 先确认粒子真的被送到了配对 patch，而不是仍留在源侧或漂到不相关 patch。
    const label newPatchI = mesh_.boundaryMesh().whichPatch(p.face());

    if (newPatchI != expectedNewPatchId)
    {
        FatalErrorInFunction
            << "Stage-interface transfer sent parcel to unexpected patch." << nl
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

    // 只有跨面与目标 patch 校验都通过后，才把本次事件记作成功交换。
    if (fromSideA)
    {
        ++transferCountAtoB_;
    }
    else
    {
        ++transferCountBtoA_;
    }

    const bool newCellInRotatingRegion = cloud_.isInRotatingRegion(p.cell());
    const vector expectedVelocity =
        transformVelocity_
      ? transformedVelocity
        (
            oldVelocity,
            mappedCrossPosition,
            fromFrame,
            toFrame
        )
      : oldVelocity;
    const vector expectedDeltaU = expectedVelocity - oldVelocity;
    const vector reconstructedVelocityFromFinalPosition =
        transformVelocity_
      ? transformedVelocity
        (
            oldVelocity,
            p.position(),
            fromFrame,
            toFrame
        )
      : oldVelocity;
    const vector reconstructedDeltaUFromFinalPosition =
        reconstructedVelocityFromFinalPosition - oldVelocity;

    // 第二步按开关决定是否做速度参考系变换。
    if (transformVelocity_)
    {
        // 写回速度之前，先确认 oldCell/newCell 与 fromFrame/toFrame
        // 的对应关系没有被 zone 划分或 patch 方向定义破坏。
        validateParcelFrameConsistency(oldCell, p.cell(), fromFrame, toFrame);
        p.U() = expectedVelocity;

        // 误差统计不能再做“赋值后自比较”：
        // - expectedVelocity / expectedDeltaU 基于 mappedCrossPosition
        //   这个真实物理取点；
        // - reconstructedVelocityFromFinalPosition /
        //   reconstructedDeltaUFromFinalPosition 基于最终插入位置 p.position()
        //   独立重构。
        //
        // 两者的差值真实反映“近似插入把 Omega x r 改变了多少”，而不是按构造
        // 永远为零。
        const scalar velocityError =
            mag(reconstructedVelocityFromFinalPosition - expectedVelocity);
        const scalar deltaUError =
            mag(reconstructedDeltaUFromFinalPosition - expectedDeltaU);

        if (fromFrame == "relative" && toFrame == "absolute")
        {
            ++transformCountRelToAbs_;
            cumulativeFinalPositionVelocityDifferenceRelToAbs_ += velocityError;
            cumulativeFinalPositionDeltaUDifferenceRelToAbs_ += deltaUError;
            maxFinalPositionVelocityDifferenceRelToAbs_ =
                max
                (
                    maxFinalPositionVelocityDifferenceRelToAbs_,
                    velocityError
                );
            maxFinalPositionDeltaUDifferenceRelToAbs_ =
                max
                (
                    maxFinalPositionDeltaUDifferenceRelToAbs_,
                    deltaUError
                );
        }
        else if (fromFrame == "absolute" && toFrame == "relative")
        {
            ++transformCountAbsToRel_;
            cumulativeFinalPositionVelocityDifferenceAbsToRel_ += velocityError;
            cumulativeFinalPositionDeltaUDifferenceAbsToRel_ += deltaUError;
            maxFinalPositionVelocityDifferenceAbsToRel_ =
                max
                (
                    maxFinalPositionVelocityDifferenceAbsToRel_,
                    velocityError
                );
            maxFinalPositionDeltaUDifferenceAbsToRel_ =
                max
                (
                    maxFinalPositionDeltaUDifferenceAbsToRel_,
                    deltaUError
                );
        }
    }

    // A/B 两侧分开计数，避免单侧高频事件耗尽另一侧的 debug 配额。
    label* sideDebugPrintCount = fromSideA ? &debugPrintCountA_ : &debugPrintCountB_;

    if (debug_ && *sideDebugPrintCount < maxDebugPrints_)
    {
        const vector r = mappedCrossPosition - rotationCentre_;
        const vector uRot = omega_ ^ r;
        const word oldPatchName =
            oldPatchI >= 0 ? mesh_.boundaryMesh()[oldPatchI].name() : word("none");
        const word newPatchName =
            newPatchI >= 0 ? mesh_.boundaryMesh()[newPatchI].name() : word("none");

        Info<< "dsmcStageInterfacePatch debug: patch=" << patchName_
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
            << " sourceCrossPosition=" << crossPosition
            << " mappedCrossPosition=" << mappedCrossPosition
            << " uRot=" << uRot
            << " oldCellInRotatingRegion=" << oldCellInRotatingRegion
            << " newCellInRotatingRegion=" << newCellInRotatingRegion
            << " oldU=" << oldVelocity
            << " expectedU=" << expectedVelocity
            << " reconstructedUFromFinalPosition="
            << reconstructedVelocityFromFinalPosition
            << " newU=" << p.U()
            << " oldRadialWeight=" << oldRadialWeight
            << " newRadialWeight=" << newRadialWeight
            << " expectedDeltaU=" << expectedDeltaU
            << " reconstructedDeltaUFromFinalPosition="
            << reconstructedDeltaUFromFinalPosition
            << " actualDeltaU=" << (p.U() - oldVelocity)
            << " fromFrame=" << fromFrame
            << " toFrame=" << toFrame
            << " transformPosition=" << transformPosition_
            << " approximateFaceCentreInsertion=" << approximateFaceCentreInsertion_
            << " approximateInsertionFraction=" << approximateInsertionFraction_
            << " transformVelocity=" << transformVelocity_
            << nl;

        ++(*sideDebugPrintCount);
    }

    // 这里只输出本 rank 本地摘要，不在事件路径里做 MPI collective。
    // controlMol() 按单颗粒命中事件触发；若在此处调用 reduce(...)，
    // 由于并非所有 rank 都同步进入该路径，会直接造成并行死锁。
    if (debug_ && time_.outputTime() && lastSummaryTimeIndex_ != time_.timeIndex())
    {
        const scalar meanFinalPositionVelocityDifferenceRelToAbs =
            transformCountRelToAbs_ > 0
          ? cumulativeFinalPositionVelocityDifferenceRelToAbs_/transformCountRelToAbs_
          : 0.0;
        const scalar meanFinalPositionVelocityDifferenceAbsToRel =
            transformCountAbsToRel_ > 0
          ? cumulativeFinalPositionVelocityDifferenceAbsToRel_/transformCountAbsToRel_
          : 0.0;
        const scalar meanFinalPositionDeltaUDifferenceRelToAbs =
            transformCountRelToAbs_ > 0
          ? cumulativeFinalPositionDeltaUDifferenceRelToAbs_/transformCountRelToAbs_
          : 0.0;
        const scalar meanFinalPositionDeltaUDifferenceAbsToRel =
            transformCountAbsToRel_ > 0
          ? cumulativeFinalPositionDeltaUDifferenceAbsToRel_/transformCountAbsToRel_
          : 0.0;

        Info<< "dsmcStageInterfacePatch summary: patch=" << patchName_
            << " time=" << time_.timeName()
            << " proc=" << Pstream::myProcNo()
            << " localAtoBTransfers=" << transferCountAtoB_
            << " localBtoATransfers=" << transferCountBtoA_
            << " relToAbsTransforms=" << transformCountRelToAbs_
            << " absToRelTransforms=" << transformCountAbsToRel_
            << " meanFinalPositionVelocityDifferenceRelToAbs="
            << meanFinalPositionVelocityDifferenceRelToAbs
            << " meanFinalPositionVelocityDifferenceAbsToRel="
            << meanFinalPositionVelocityDifferenceAbsToRel
            << " maxFinalPositionVelocityDifferenceRelToAbs="
            << maxFinalPositionVelocityDifferenceRelToAbs_
            << " maxFinalPositionVelocityDifferenceAbsToRel="
            << maxFinalPositionVelocityDifferenceAbsToRel_
            << " meanFinalPositionDeltaUDifferenceRelToAbs="
            << meanFinalPositionDeltaUDifferenceRelToAbs
            << " meanFinalPositionDeltaUDifferenceAbsToRel="
            << meanFinalPositionDeltaUDifferenceAbsToRel
            << " maxFinalPositionDeltaUDifferenceRelToAbs="
            << maxFinalPositionDeltaUDifferenceRelToAbs_
            << " maxFinalPositionDeltaUDifferenceAbsToRel="
            << maxFinalPositionDeltaUDifferenceAbsToRel_
            << " debugPrintsA=" << debugPrintCountA_
            << " debugPrintsB=" << debugPrintCountB_
            << nl;

        lastSummaryTimeIndex_ = time_.timeIndex();
    }
}


void dsmcStageInterfacePatch::output
(
    const fileName&,
    const fileName&
)
{}


void dsmcStageInterfacePatch::updateProperties(const dictionary& newDict)
{
    // 支持运行中重读字典参数。
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
}

} // End namespace Foam

// ************************************************************************* //
