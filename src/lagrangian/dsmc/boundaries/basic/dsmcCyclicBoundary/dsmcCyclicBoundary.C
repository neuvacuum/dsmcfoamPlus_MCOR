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

#include "dsmcCyclicBoundary.H"
#include "polyBoundaryMeshEntries.H"
#include "graph.H"
#include "processorCyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

namespace Foam
{

namespace
{
    vector normaliseRotationAxis(const vector& axis)
    {
        const scalar axisMag = mag(axis);

        if (axisMag <= ROOTVSMALL)
        {
            return vector(1, 0, 0);
        }

        return axis/axisMag;
    }


    scalar axialCoordinateAboutAxis
    (
        const point& position,
        const point& rotationCentre,
        const vector& rotationAxis
    )
    {
        const vector relative = position - rotationCentre;
        return relative & rotationAxis;
    }


    scalar cylindricalRadiusAboutAxis
    (
        const point& position,
        const point& rotationCentre,
        const vector& rotationAxis
    )
    {
        const vector relative = position - rotationCentre;
        const scalar axial = relative & rotationAxis;
        const vector radial = relative - axial*rotationAxis;

        return mag(radial);
    }

}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dsmcCyclicBoundary, 0);

defineRunTimeSelectionTable(dsmcCyclicBoundary, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcCyclicBoundary::dsmcCyclicBoundary
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    mesh_(refCast<const fvMesh>(mesh)),
    cloud_(cloud),
    time_(t),
//     rndGen_(cloud_.rndGen()),
    boundaryDict_(dict.subDict("cyclicBoundaryProperties")),

    patchName_(),
    patchId_(0),
    neighbPatchName_(),
    neighbPatchId_(0),

    faces_(),
    coupledFacesA_(),
    coupledFacesB_(),
//     cells_(),
    cellsA_(),
    cellsB_(),
    densities_(),
    velocities_(),
    temperatures_(),
    writeInTimeDir_(true),
    writeInCase_(true),
    parallelSearchDebug_(false),
    maxParallelSearchDebugPrints_(0),
    parallelSearchDebugPrintCount_(0)
{
    refreshCoupledPatchState();
}


void dsmcCyclicBoundary::refreshCoupledPatchState()
{
    patchName_ = word(boundaryDict_.lookup("patchName"));
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchId_ == -1)
    {
        FatalErrorInFunction
            << "Cannot find cyclic patch '" << patchName_ << "' in "
            << time_.system()/"boundariesDict"
            << exit(FatalError);
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];
    word meshNeighbourName;
    bool sourceIsCyclicAMI = false;

    if (isA<cyclicPolyPatch>(patch))
    {
        meshNeighbourName =
            refCast<const cyclicPolyPatch>(patch).neighbPatchName();
    }
    else if (isA<cyclicAMIPolyPatch>(patch))
    {
        meshNeighbourName =
            refCast<const cyclicAMIPolyPatch>(patch).neighbPatchName();
        sourceIsCyclicAMI = true;
    }
    else
    {
        FatalErrorInFunction
            << "Patch '" << patchName_
            << "' is neither cyclic nor cyclicAMI."
            << exit(FatalError);
    }

    const word explicitNeighbourName =
        boundaryDict_.lookupOrDefault<word>
        (
            "neighbourPatchName",
            word::null
        );

    if
    (
        !explicitNeighbourName.empty()
     && explicitNeighbourName != meshNeighbourName
    )
    {
        FatalErrorInFunction
            << "Configured neighbourPatchName='" << explicitNeighbourName
            << "' does not match mesh coupling for patch '" << patchName_
            << "', which names neighbour '" << meshNeighbourName << "'." << nl
            << "The DSMC model must use the same reciprocal patch pair as "
            << "polyMesh; fix boundariesDict or the mesh boundary definition."
            << exit(FatalError);
    }

    neighbPatchName_ = meshNeighbourName;
    neighbPatchId_ = mesh_.boundaryMesh().findPatchID(neighbPatchName_);

    if (neighbPatchId_ == -1 || neighbPatchId_ == patchId_)
    {
        FatalErrorInFunction
            << "Invalid neighbour patch '" << neighbPatchName_
            << "' for cyclic patch '" << patchName_ << "'."
            << exit(FatalError);
    }

    const polyPatch& neighbourPatch =
        mesh_.boundaryMesh()[neighbPatchId_];
    word reciprocalNeighbourName;

    if (sourceIsCyclicAMI && isA<cyclicAMIPolyPatch>(neighbourPatch))
    {
        reciprocalNeighbourName =
            refCast<const cyclicAMIPolyPatch>(neighbourPatch).neighbPatchName();
    }
    else if (!sourceIsCyclicAMI && isA<cyclicPolyPatch>(neighbourPatch))
    {
        reciprocalNeighbourName =
            refCast<const cyclicPolyPatch>(neighbourPatch).neighbPatchName();
    }
    else
    {
        FatalErrorInFunction
            << "Patch pair ('" << patchName_ << "', '" << neighbPatchName_
            << "') mixes cyclic and cyclicAMI patch types."
            << exit(FatalError);
    }

    if (reciprocalNeighbourName != patchName_)
    {
        FatalErrorInFunction
            << "Patch pair is not reciprocal: '" << patchName_
            << "' names '" << neighbPatchName_ << "', but '"
            << neighbPatchName_ << "' names '" << reciprocalNeighbourName
            << "'."
            << exit(FatalError);
    }

    parallelSearchDebug_ =
        boundaryDict_.lookupOrDefault<bool>("parallelSearchDebug", false);
    maxParallelSearchDebugPrints_ =
        boundaryDict_.lookupOrDefault<label>
        (
            "maxParallelSearchDebugPrints",
            0
        );
    parallelSearchDebugPrintCount_ = 0;

    getCoupledFaces(mesh_);
}

void dsmcCyclicBoundary::getCoupledFaces(const polyMesh& mesh)
{
    // in parallel, openFOAM replaces cyclic boundaries with processor boundaries
    // if processor boundaries coincide with cyclic boundaries. We therefore need
    // to do the following corrections to get the correct faces on the selected
    // cyclic patch.

    if (Pstream::parRun())
    {
        DynamicList<label> coupledFacesA(0);
        DynamicList<label> coupledFacesB(0);
        DynamicList<label> cellsA(0);
        DynamicList<label> cellsB(0);


        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        for(label i = 0; i < patch.size(); i++)
        {
            label globalFaceI = patch.start() + i;
            label cellI = patch.faceCells()[i];
            coupledFacesA.append(globalFaceI);
            cellsA.append(cellI);
        }

        const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

        for(label i = 0; i < patchN.size(); i++)
        {
            label globalFaceI = patchN.start() + i;
            label cellI = patchN.faceCells()[i];
            coupledFacesB.append(globalFaceI);
            cellsB.append(cellI);
        }

        // 额外补入 processorCyclicPolyPatch 上的局部 faces。
        // 这样在并行划分后，coupledFacesA_/B_ 仍然覆盖本 rank 上
        // 与这对 cyclic/cyclicAMI 边界相关的全部局部 faces。

        forAll( mesh_.boundaryMesh(), patchI )
        {
            const polyPatch& patch = mesh.boundaryMesh()[patchI];

            if(isA<processorCyclicPolyPatch>(patch))
            {
                word patchName = refCast<const processorCyclicPolyPatch>
                (
                    patch
                ).referPatchName();

                if(patchName == patchName_)
                {
                    for(label i = 0; i < patch.size(); i++)
                    {
                        label globalFaceI = patch.start() + i;
                        label cellI = patch.faceCells()[i];
                        coupledFacesA.append(globalFaceI);
                        cellsA.append(cellI);
                    }
                }
                else if(patchName == neighbPatchName_)
                {
                    for(label i = 0; i < patch.size(); i++)
                    {
                        label globalFaceI = patch.start() + i;
                        label cellI = patch.faceCells()[i];
                        coupledFacesB.append(globalFaceI);
                        cellsB.append(cellI);
                    }
                }
            }
        }

        coupledFacesA.shrink();
        coupledFacesB.shrink();
        cellsA.shrink();
        cellsB.shrink();

        coupledFacesA_.transfer(coupledFacesA);
        coupledFacesB_.transfer(coupledFacesB);
        cellsA_.transfer(cellsA);
        cellsB_.transfer(cellsB);
    }
    else
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

        coupledFacesA_.setSize(patch.size(), -1);
        cellsA_.setSize(patch.size(), -1);

        for(label i = 0; i < patch.size(); i++)
        {
            label globalFaceI = patch.start() + i;
            label cellI = patch.faceCells()[i];
            coupledFacesA_[i] = globalFaceI;
            cellsA_[i] = cellI;
        }

        const polyPatch& patchN = mesh_.boundaryMesh()[neighbPatchId_];

        coupledFacesB_.setSize(patchN.size(), -1);
        cellsB_.setSize(patchN.size(), -1);

        for(label i = 0; i < patchN.size(); i++)
        {
            label globalFaceI = patchN.start() + i;
            label cellI = patchN.faceCells()[i];
            coupledFacesB_[i] = globalFaceI;
            cellsB_[i] = cellI;
        }
    }

    label totalFaces = coupledFacesA_.size() + coupledFacesB_.size();

    faces_.setSize(totalFaces, -1);

    label index = 0;

    forAll(coupledFacesA_, fA)
    {
        faces_[index] = coupledFacesA_[fA];
        index++;
    }

    forAll(coupledFacesB_, fB)
    {
        faces_[index] = coupledFacesB_[fB];
        index++;
    }

//     Info << "coupledFacesA: "  << coupledFacesA_ << endl;
//     Info << "coupledFacesB: "  << coupledFacesB_ << endl;
//     Info << "totalFaces: "  << faces_ << endl;

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<dsmcCyclicBoundary> dsmcCyclicBoundary::New
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
{
    word dsmcCyclicBoundaryName
    (
        dict.lookup("boundaryModel")
    );

    Info<< "Selecting dsmcCyclicBoundaryModel "
         << dsmcCyclicBoundaryName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(dsmcCyclicBoundaryName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "dsmcCyclicBoundary::New(const dictionary&) : " << endl
            << "    unknown dsmcCyclicBoundary type "
            << dsmcCyclicBoundaryName
            << ", constructor not in hash table" << endl << endl
            << "    Valid  types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<dsmcCyclicBoundary>
	  (
		    cstrIter()(t, mesh, cloud, dict)
	  );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcCyclicBoundary::~dsmcCyclicBoundary()
{}


void dsmcCyclicBoundary::applyStandardCyclicTransfer
(
    dsmcParcel& p,
    const bool transformPosition
) const
{
    // 这里实现的是对 OpenFOAM 原生 cyclic 拓扑路径的统一复用。
    // 只要派生类已经确认粒子确实命中了它所管理的 cyclic 面，
    // 就应通过这条路径完成跨面迁移，而不应在各个模型里重复手写
    // face/cell/tet/tetPt 的更新规则。
    const label faceI = p.face();

    // 当前 face 必须属于某个真实的 boundary patch；
    // 如果 whichPatch 返回 -1，说明上游调用时机已经不合法。
    const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

    if (patchI == -1)
    {
        FatalErrorInFunction
            << "Face " << faceI
            << " is not on a boundary patch during cyclic transfer."
            << exit(FatalError);
    }

    const polyPatch& patch = mesh_.boundaryMesh()[patchI];

    if (!isA<cyclicPolyPatch>(patch))
    {
        FatalErrorInFunction
            << "Patch " << patch.name()
            << " is not cyclic, but applyStandardCyclicTransfer was called."
            << exit(FatalError);
    }

    const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(patch);
    const point positionBeforeTransfer = p.position();

    // transformGlobalFace 会返回接收侧与当前 patchFace 一一对应的全局 face。
    // 这一步是 OpenFOAM 对 cyclic 配对关系的标准定义。
    p.face() = cpp.transformGlobalFace(faceI);

    // face 切换后，粒子所属 cell 也必须切到接收侧 face owner。
    p.cell() = mesh_.faceOwner()[p.face()];

    // 对 DSMC lagrangian 跟踪来说，tetFace/tetPt 共同决定粒子在
    // 局部面四面体分解中的拓扑状态。这里保持与标准 cyclic 命中逻辑一致。
    p.tetFace() = p.face();
    p.tetPt() = mesh_.faces()[p.tetFace()].size() - 1 - p.tetPt();

    // 只有在模型显式要求时，才对粒子几何位置执行接收侧 patch 的标准变换。
    // 这样做可以把“只换拓扑不挪位置”的调试需求与正式迁移需求区分开。
    if (transformPosition)
    {
        const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
        const label patchFaceI = receiveCpp.whichFace(p.face());

        if (patchFaceI == -1)
        {
            FatalErrorInFunction
                << "Unable to locate transformed face " << p.face()
                << " on neighbour cyclic patch " << receiveCpp.name()
                << exit(FatalError);
        }

        receiveCpp.transformPosition(p.position(), patchFaceI);
    }

    const point& ownerCellCentre = mesh_.cellCentres()[p.cell()];
    vector inwardDirection = ownerCellCentre - p.position();

    if (mag(inwardDirection) <= ROOTVSMALL)
    {
        // 极端情况下如果粒子正好与 cell centre 重合，就退化为沿接收面内法向推进。
        const vector faceNormal = mesh_.faceAreas()[p.face()];

        if (mag(faceNormal) > ROOTVSMALL)
        {
            inwardDirection = -faceNormal;
        }
    }

    const scalar inwardMag = mag(inwardDirection);

    if (inwardMag > ROOTVSMALL)
    {
        // 位移尺度取接收面特征长度的百万分之一。
        // 这样足以让粒子脱离病态边界容差区，又不会在物理上造成可见扰动。
        const scalar faceLengthScale = Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));
        const scalar pushDistance = max(ROOTVSMALL, 1e-6*faceLengthScale);

        p.position() += (pushDistance/inwardMag)*inwardDirection;
    }

    // 只要跨面后真实改写了几何位置，就必须重建局部 tet 拓扑状态。
    // 否则粒子虽然已被送入接收侧单元内部，但 tetFace/tetPt 仍可能保留为
    // 命中周期面那一刻的旧局部状态。
    if (mag(p.position() - positionBeforeTransfer) > ROOTVSMALL)
    {
        p.initCellFacePt();
    }

}


bool dsmcCyclicBoundary::tryStandardCyclicAMITransfer
(
    dsmcParcel& p,
    const vector& direction,
    word* searchModeLabel,
    const point* targetReferencePosition,
    const point* targetInsertionReferencePosition,
    const vector* rotationAxis,
    const point* rotationCentre
) const
{
    // 这里统一复用 OpenFOAM 原生 AMI 粒子跨面路径的核心拓扑更新步骤。
    // 对 AMI 而言，目标面不是固定一一对应，而是要根据粒子当前几何位置
    // 和穿越方向在接收侧 patch 上搜索得到。因此辅助函数必须显式接收
    // direction，供 cyclicAMI 的 pointFace() 使用。
    const label faceI = p.face();
    const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

    if (patchI == -1)
    {
        FatalErrorInFunction
            << "Face " << faceI
            << " is not on a boundary patch during cyclicAMI transfer."
            << exit(FatalError);
    }

    const polyPatch& sourcePatch = mesh_.boundaryMesh()[patchI];

    if (!isA<cyclicAMIPolyPatch>(sourcePatch))
    {
        FatalErrorInFunction
            << "Patch " << sourcePatch.name()
            << " is not cyclicAMI, but tryStandardCyclicAMITransfer was called."
            << exit(FatalError);
    }

    const cyclicAMIPolyPatch& cpp = refCast<const cyclicAMIPolyPatch>(sourcePatch);

    // 当前 dsmcRotationalCyclicAMIPatch 只在 boundariesDict 中注册一次，
    // 但它实际上要同时接管 periodicA 和 periodicB 两张面。
    // 因此这里绝不能把“接收 patch”写死成 neighbPatchId_，否则会出现：
    // 1. A->B 时看起来正常，因为源 patch 恰好是 patchId_；
    // 2. B->A 时仍然在 periodicB 上继续搜索目标面，导致粒子被送回原角度附近。
    // 这正是前面日志中 side=B 但 dThetaAfterInsertion≈0 的根因。
    //
    // 正确做法是根据“当前真正命中的源 patch”动态决定接收 patch：
    //     若源 patch 是 patchId_       ，则接收 patch 应为 neighbPatchId_
    //     若源 patch 是 neighbPatchId_ ，则接收 patch 应为 patchId_
    // 也就是始终把粒子送到这对周期面的另一侧。
    const bool sourceOnPrimaryPatch = (patchI == patchId_);
    const bool sourceOnNeighbourPatch = (patchI == neighbPatchId_);

    if (!sourceOnPrimaryPatch && !sourceOnNeighbourPatch)
    {
        FatalErrorInFunction
            << "Patch " << sourcePatch.name()
            << " is not part of the cyclic pair ("
            << patchName_ << ", " << neighbPatchName_ << ") handled by this model."
            << exit(FatalError);
    }

    const label receivePatchId =
        sourceOnPrimaryPatch ? neighbPatchId_ : patchId_;
    const polyPatch& receivePatch = mesh_.boundaryMesh()[receivePatchId];

    label receivePatchFaceI = -1;

    point searchPosition = p.position();

    // 并行非严格匹配 AMI 下，派生类可以额外传入一个 caller-defined
    // targetReferencePosition，作为目标面评分的主导几何参考点。
    // 这个 referencePosition 不是固定只能取“源面中心”：
    // 1. 对 rotational cyclic AMI，可传“旋转后的源侧面中心”；
    // 2. 对 stage-interface AMI，也可以直接传粒子真实 crossing point；
    // 3. 若未提供，则退回到 searchPosition，保持与旧行为兼容。
    const point referencePosition =
        targetReferencePosition ? *targetReferencePosition : searchPosition;

    // 对普通非旋转 AMI，调用方可以不提供 rotationAxis / rotationCentre。
    // 这时这里仍保持旧行为，默认按 x 轴、原点定义柱坐标不变量。
    // 但对于当前 z 轴 rotational cyclic AMI，派生类会显式传入
    // 真实的 rotationAxis_ / rotationCentre_，从而把局部目标面评分中的
    // “轴向 + 半径”几何不变量真正切换到对应旋转轴上。
    const vector searchAxis =
        normaliseRotationAxis(rotationAxis ? *rotationAxis : vector(1, 0, 0));
    const point searchCentre =
        rotationCentre ? *rotationCentre : point::zero;

    if (Pstream::parRun())
    {
        // 并行时不能再无条件调用 cpp.pointFace()。
        // 原因是当前 v1706/hyStrath 路径下，局部 processor 上的 cyclicAMI
        // patch 片段一旦不满足 OpenFOAM 内部的 source/target 建址假设，
        // pointFace() 会在内部直接 FatalError，调用方无法像普通返回值那样
        // 捕获失败并继续尝试其他搜索方向。
        //
        // 因此并行分支改为：
        // 1. 利用 boundaryDict 中显式给出的 neighbourPatchName，直接拿到
        //    当前 processor 上接收侧 patch 的本地面片；
        // 2. 以“旋转后粒子位置”作为几何目标，在本地接收 patch 上做最近面搜索；
        // 3. 若本地根本没有接收侧面片，或最近面仍明显不可信，则返回 false，
        //    交回派生类继续尝试下一组搜索方向，而不是在这里直接崩溃。
        //
        // 这样做牺牲的是一部分对 OpenFOAM 原生 AMI 建址的直接依赖，
        // 换来的是并行 DSMC 粒子跨面路径的可控性和可回退性。
        const label receivePatchSize = receivePatch.size();

        if (receivePatchSize <= 0)
        {
            return false;
        }

        // 并行非严格匹配 AMI 下，局部 processor 往往只持有接收 patch 的一小片。
        // 如果仍然沿用过于严格的“局部尺度阈值”过滤，最容易出现的失败模式不是
        // “选错面”，而是“所有候选都被过滤掉，导致第一次真实跨面事件立刻失败”。
        //
        // 因此这里采用两级策略：
        // 1. strict 候选：继续保留局部尺度过滤，优先使用更可信的目标面；
        // 2. loose  候选：即使 strict 全部失败，也保留一个全局最优最近面，
        //    作为并行局部片段场景下的保底退路。
        //
        // 这样做的目的不是放弃几何判别，而是把“过于保守导致完全无解”
        // 改成“先求一个局部最合理解，再交给后续追踪继续暴露问题”。
        scalar bestStrictScore = GREAT;
        label bestStrictLocalFace = -1;
        label bestStrictGlobalFace = -1;
        scalar bestStrictCentreDistance = GREAT;
        scalar bestStrictPlaneDistance = GREAT;
        scalar bestStrictRadialDistance = GREAT;
        scalar bestStrictAxialDistance = GREAT;
        scalar bestStrictReferenceCentreDistance = GREAT;
        scalar bestStrictReferencePlaneDistance = GREAT;
        scalar bestStrictReferenceRadialDistance = GREAT;
        scalar bestStrictReferenceAxialDistance = GREAT;
        scalar bestStrictFaceRadius = ROOTVSMALL;
        scalar bestStrictDirectionalPenalty = GREAT;
        point bestStrictCentre(Zero);

        scalar bestLooseScore = GREAT;
        label bestLooseLocalFace = -1;
        label bestLooseGlobalFace = -1;
        scalar bestLooseCentreDistance = GREAT;
        scalar bestLoosePlaneDistance = GREAT;
        scalar bestLooseRadialDistance = GREAT;
        scalar bestLooseAxialDistance = GREAT;
        scalar bestLooseReferenceCentreDistance = GREAT;
        scalar bestLooseReferencePlaneDistance = GREAT;
        scalar bestLooseReferenceRadialDistance = GREAT;
        scalar bestLooseReferenceAxialDistance = GREAT;
        scalar bestLooseFaceRadius = ROOTVSMALL;
        scalar bestLooseDirectionalPenalty = GREAT;
        point bestLooseCentre(Zero);

        const scalar directionMag = mag(direction);
        const vector unitDirection =
            directionMag > ROOTVSMALL ? direction/directionMag : vector::zero;

        // 对任意旋转轴的旋转周期边界来说，理想映射必须同时满足：
        // 1. 轴向坐标 s = (x - c) · a_hat 保持不变；
        // 2. 半径   r = |(x - c) - s a_hat| 保持不变；
        // 3. 仅绕 a_hat 的极角发生 ±Δtheta 变化。
        //
        // 因此前一版只用 centreDistance / planeDistance / directionalPenalty
        // 的最近面评分还不够，它会把“整体最近但半径明显不对”的面也选进来。
        // 这里显式把 searchPosition 的“任意轴柱坐标不变量”提出来，
        // 后面对每个候选面额外计算 axialDistance 和 radialDistance。
        const scalar searchRadius =
            cylindricalRadiusAboutAxis(searchPosition, searchCentre, searchAxis);
        const scalar searchAxial =
            axialCoordinateAboutAxis(searchPosition, searchCentre, searchAxis);

        // referencePosition 是由派生类决定的“目标面评分主导点”。
        // 对不同模型它可以代表不同物理含义：
        // 1. rotational cyclic AMI 中，它常取旋转后的源侧面中心；
        // 2. stage-interface AMI 中，它也可以直接等于 crossing point；
        // 因此这里不再把它写死为某一种几何来源，只保留“由调用方定义”的语义。
        const scalar referenceRadius =
            cylindricalRadiusAboutAxis(referencePosition, searchCentre, searchAxis);
        const scalar referenceAxial =
            axialCoordinateAboutAxis(referencePosition, searchCentre, searchAxis);

        // 先做一轮 reference scan，找出“旋转后源面中心”在接收 patch 上最接近的
        // 局部 (x,r) 参考带。后续 strict / loose 都只在这条参考带附近找目标面，
        // 避免继续在整块接收 patch 上全局比距离，把粒子吸到错误的靠壁面带。
        scalar minReferenceBandScore = GREAT;
        scalar minReferenceRadialDistance = GREAT;
        scalar minReferenceAxialDistance = GREAT;

        for (label localFaceI = 0; localFaceI < receivePatchSize; ++localFaceI)
        {
            const label globalFaceI = receivePatch.start() + localFaceI;
            const point candidateCentre = mesh_.faceCentres()[globalFaceI];
            const vector candidateArea = mesh_.faceAreas()[globalFaceI];
            const scalar faceAreaMag = mag(candidateArea);

            if (faceAreaMag <= ROOTVSMALL)
            {
                continue;
            }

            const scalar candidateRadius =
                cylindricalRadiusAboutAxis(candidateCentre, searchCentre, searchAxis);
            const scalar referenceRadialDistance =
                mag(referenceRadius - candidateRadius);
            const scalar referenceAxialDistance =
                mag
                (
                    referenceAxial
                  - axialCoordinateAboutAxis(candidateCentre, searchCentre, searchAxis)
                );
            const scalar referenceBandScore =
                referenceRadialDistance + referenceAxialDistance;

            if (referenceBandScore < minReferenceBandScore)
            {
                minReferenceBandScore = referenceBandScore;
                minReferenceRadialDistance = referenceRadialDistance;
                minReferenceAxialDistance = referenceAxialDistance;
            }
        }

        if (minReferenceBandScore >= GREAT)
        {
            return false;
        }

        for (label localFaceI = 0; localFaceI < receivePatchSize; ++localFaceI)
        {
            const label globalFaceI = receivePatch.start() + localFaceI;
            const face& candidateFace = mesh_.faces()[globalFaceI];
            const point candidateCentre = mesh_.faceCentres()[globalFaceI];
            const vector candidateArea = mesh_.faceAreas()[globalFaceI];
            const scalar faceAreaMag = mag(candidateArea);

            if (faceAreaMag <= ROOTVSMALL)
            {
                continue;
            }

            // 先用接收面中心距离筛选最可能的目标面。
            const vector centreOffset = searchPosition - candidateCentre;
            const scalar centreDistance = mag(centreOffset);

            // 再叠加与面法向的偏离程度，避免在同一 patch 局部片段上
            // 仅按中心距离选到“几何上邻近但法向明显不匹配”的面。
            const vector candidateNormal = candidateArea/faceAreaMag;
            const scalar planeDistance = mag(centreOffset & candidateNormal);

            // 若搜索方向与接收面法向几乎完全平行于面内，局部目标面搜索会更敏感。
            // 这里把法向投影也计入评分，优先选择与当前跨面方向更一致的面。
            const scalar directionalPenalty =
                directionMag > ROOTVSMALL
              ? 1.0 - mag(unitDirection & candidateNormal)
              : 0.0;

            // 额外计算“候选面中心相对于 referencePosition 的偏差”。
            // 与 searchPosition 相比，referencePosition 不那么容易受到贴面推进、
            // 局部 processor 片段截断和近壁单元尺度差异的污染，更适合表达
            // “周期面几何上真正应该对应到哪里”。
            const vector referenceOffset = referencePosition - candidateCentre;
            const scalar referenceCentreDistance = mag(referenceOffset);
            const scalar referencePlaneDistance = mag(referenceOffset & candidateNormal);

            const scalar candidateRadius =
                cylindricalRadiusAboutAxis(candidateCentre, searchCentre, searchAxis);
            const scalar radialDistance = mag(searchRadius - candidateRadius);
            const scalar axialDistance =
                mag
                (
                    searchAxial
                  - axialCoordinateAboutAxis(candidateCentre, searchCentre, searchAxis)
                );
            const scalar referenceRadialDistance =
                mag(referenceRadius - candidateRadius);
            const scalar referenceAxialDistance =
                mag
                (
                    referenceAxial
                  - axialCoordinateAboutAxis(candidateCentre, searchCentre, searchAxis)
                );

            // 面尺度取该面顶点到面心最大距离，作为局部可接受搜索半径。
            scalar faceRadius = ROOTVSMALL;
            forAll(candidateFace, pointI)
            {
                const point& fp = mesh_.points()[candidateFace[pointI]];
                faceRadius = max(faceRadius, mag(fp - candidateCentre));
            }

            // 这里的评分主导量始终是 referencePosition，而 referencePosition
            // 的来源由派生类控制：
            // 1. 若派生类传入“旋转后的源侧面中心”，则目标面会按周期面区块语义选取；
            // 2. 若派生类传入真实 crossing point，则目标面会按粒子瞬时 crossing
            //    几何位置选取；
            // 3. searchPosition 仍保留在局部距离与柱坐标不变量计算里，但不再单独
            //    作为最终评分的主导量。
            //
            // 一旦目标面确定，后面直接把粒子插到目标面中心附近，再轻推入 owner cell。
            //
            // 这么做的物理含义是：
            //     先保证粒子落到 referencePosition 所代表的正确局部区块
            // 再接受面内位置只是近似，而不是追求当前非严格匹配网格上
            // 很难稳定获得的“精确面内对应点”。
            //
            // 因此这里把评分函数定义成 referencePosition 主导。
            // referenceCentreDistance / referencePlaneDistance 分别代表：
            // 1. 候选面中心与调用方定义的 referencePosition 之间的距离；
            // 2. referencePosition 到候选面平面的法向距离。
            // 只要这两个量小，说明候选面在周期面几何对应关系上更可信。
            const scalar referenceScore =
                referenceCentreDistance
              + referencePlaneDistance
              + 0.25*faceRadius*directionalPenalty;

            const scalar localReferenceRadialTol =
                minReferenceRadialDistance
              + max(2.0*faceRadius, 0.01*max(referenceRadius, ROOTVSMALL));
            const scalar localReferenceAxialTol =
                minReferenceAxialDistance
              + max(2.0*faceRadius, 0.008*max(referenceRadius, ROOTVSMALL));

            const bool inLocalReferenceBand =
                (referenceRadialDistance <= localReferenceRadialTol)
             && (referenceAxialDistance <= localReferenceAxialTol);

            // loose 仍然保留为 fallback，但只能在 local reference band 内工作。
            const scalar looseDirectionalTol = 0.995;
            const bool looseReferenceAcceptable =
                inLocalReferenceBand
             && (directionalPenalty <= looseDirectionalTol);

            if (looseReferenceAcceptable && referenceScore < bestLooseScore)
            {
                bestLooseScore = referenceScore;
                bestLooseLocalFace = localFaceI;
                bestLooseGlobalFace = globalFaceI;
                bestLooseCentreDistance = centreDistance;
                bestLoosePlaneDistance = planeDistance;
                bestLooseRadialDistance = radialDistance;
                bestLooseAxialDistance = axialDistance;
                bestLooseReferenceCentreDistance = referenceCentreDistance;
                bestLooseReferencePlaneDistance = referencePlaneDistance;
                bestLooseReferenceRadialDistance = referenceRadialDistance;
                bestLooseReferenceAxialDistance = referenceAxialDistance;
                bestLooseFaceRadius = faceRadius;
                bestLooseDirectionalPenalty = directionalPenalty;
                bestLooseCentre = candidateCentre;
            }

            // strict 在 loose 的基础上再收紧，只允许位于更窄的 local band 内。
            const scalar strictReferenceRadialTol =
                minReferenceRadialDistance
              + max(0.75*faceRadius, 0.004*max(referenceRadius, ROOTVSMALL));
            const scalar strictReferenceAxialTol =
                minReferenceAxialDistance
              + max(0.75*faceRadius, 0.003*max(referenceRadius, ROOTVSMALL));
            const scalar strictDirectionalTol = 0.99;

            const bool strictReferenceAcceptable =
                (referenceRadialDistance <= strictReferenceRadialTol)
             && (referenceAxialDistance <= strictReferenceAxialTol)
             && (directionalPenalty <= strictDirectionalTol);

            if (strictReferenceAcceptable && referenceScore < bestStrictScore)
            {
                bestStrictScore = referenceScore;
                bestStrictLocalFace = localFaceI;
                bestStrictGlobalFace = globalFaceI;
                bestStrictCentreDistance = centreDistance;
                bestStrictPlaneDistance = planeDistance;
                bestStrictRadialDistance = radialDistance;
                bestStrictAxialDistance = axialDistance;
                bestStrictReferenceCentreDistance = referenceCentreDistance;
                bestStrictReferencePlaneDistance = referencePlaneDistance;
                bestStrictReferenceRadialDistance = referenceRadialDistance;
                bestStrictReferenceAxialDistance = referenceAxialDistance;
                bestStrictFaceRadius = faceRadius;
                bestStrictDirectionalPenalty = directionalPenalty;
                bestStrictCentre = candidateCentre;
            }
        }
        if (bestStrictLocalFace >= 0)
        {
            receivePatchFaceI = bestStrictLocalFace;

            // 并行局部搜索最终成功且命中了 strict 候选。
            // 若调用方提供了 searchModeLabel，则把当前事件标记为 strict，
            // 供上层边界模型继续做 strict / loose 分模式统计。
            if (searchModeLabel)
            {
                *searchModeLabel = "strict";
            }
        }
        else if (bestLooseLocalFace >= 0)
        {
            // strict 全空时，退回到局部 patch 上的最优最近面。
            // 这是专门为并行 AMI 局部片段场景准备的保底路径。
            receivePatchFaceI = bestLooseLocalFace;

            // 若最终只能依赖保底最近面，也把模式显式返回给上层。
            // 这样后续统计就能直接回答：剩余宏观量偏差到底主要来自
            // loose 退回，还是 strict 主路径本身。
            if (searchModeLabel)
            {
                *searchModeLabel = "loose";
            }
        }
        else
        {
            if
            (
                parallelSearchDebug_
             && parallelSearchDebugPrintCount_ < maxParallelSearchDebugPrints_
            )
            {
                ++parallelSearchDebugPrintCount_;

                Info<< "dsmcCyclicBoundary parallelSearch debug: mode=failed"
                    << " proc=" << Pstream::myProcNo()
                    << " sourcePatch=" << sourcePatch.name()
                    << " receivePatch=" << receivePatch.name()
                    << " sourceFace=" << faceI
                    << " receivePatchSize=" << receivePatchSize
                    << " searchPosition=" << searchPosition
                    << " searchRadius="
                    << cylindricalRadiusAboutAxis(searchPosition, searchCentre, searchAxis)
                    << " referencePosition=" << referencePosition
                    << " referenceRadius="
                    << cylindricalRadiusAboutAxis(referencePosition, searchCentre, searchAxis)
                    << endl;
            }

            return false;
        }

        if
        (
            parallelSearchDebug_
         && parallelSearchDebugPrintCount_ < maxParallelSearchDebugPrints_
        )
        {
            ++parallelSearchDebugPrintCount_;

            const bool usedStrict = (bestStrictLocalFace >= 0);
            const label chosenGlobalFace = usedStrict ? bestStrictGlobalFace : bestLooseGlobalFace;
            const scalar chosenScore = usedStrict ? bestStrictScore : bestLooseScore;
            const scalar chosenCentreDistance = usedStrict ? bestStrictCentreDistance : bestLooseCentreDistance;
            const scalar chosenPlaneDistance = usedStrict ? bestStrictPlaneDistance : bestLoosePlaneDistance;
            const scalar chosenRadialDistance = usedStrict ? bestStrictRadialDistance : bestLooseRadialDistance;
            const scalar chosenAxialDistance = usedStrict ? bestStrictAxialDistance : bestLooseAxialDistance;
            const scalar chosenReferenceCentreDistance =
                usedStrict ? bestStrictReferenceCentreDistance : bestLooseReferenceCentreDistance;
            const scalar chosenReferencePlaneDistance =
                usedStrict ? bestStrictReferencePlaneDistance : bestLooseReferencePlaneDistance;
            const scalar chosenReferenceRadialDistance =
                usedStrict ? bestStrictReferenceRadialDistance : bestLooseReferenceRadialDistance;
            const scalar chosenReferenceAxialDistance =
                usedStrict ? bestStrictReferenceAxialDistance : bestLooseReferenceAxialDistance;
            const scalar chosenFaceRadius = usedStrict ? bestStrictFaceRadius : bestLooseFaceRadius;
            const scalar chosenDirectionalPenalty =
                usedStrict ? bestStrictDirectionalPenalty : bestLooseDirectionalPenalty;
            const point& chosenCentre = usedStrict ? bestStrictCentre : bestLooseCentre;

            Info<< "dsmcCyclicBoundary parallelSearch debug: mode="
                << (usedStrict ? "strict" : "loose")
                << " proc=" << Pstream::myProcNo()
                << " sourcePatch=" << sourcePatch.name()
                << " receivePatch=" << receivePatch.name()
                << " sourceFace=" << faceI
                << " chosenLocalFace=" << receivePatchFaceI
                << " chosenGlobalFace=" << chosenGlobalFace
                << " receivePatchSize=" << receivePatchSize
                << " searchPosition=" << searchPosition
                << " searchRadius="
                << cylindricalRadiusAboutAxis(searchPosition, searchCentre, searchAxis)
                << " referencePosition=" << referencePosition
                << " referenceRadius="
                << cylindricalRadiusAboutAxis(referencePosition, searchCentre, searchAxis)
                << " chosenCentre=" << chosenCentre
                << " chosenRadius="
                << cylindricalRadiusAboutAxis(chosenCentre, searchCentre, searchAxis)
                << " chosenScore=" << chosenScore
                << " centreDistance=" << chosenCentreDistance
                << " planeDistance=" << chosenPlaneDistance
                << " radialDistance=" << chosenRadialDistance
                << " axialDistance=" << chosenAxialDistance
                << " referenceCentreDistance=" << chosenReferenceCentreDistance
                << " referencePlaneDistance=" << chosenReferencePlaneDistance
                << " referenceRadialDistance=" << chosenReferenceRadialDistance
                << " referenceAxialDistance=" << chosenReferenceAxialDistance
                << " faceRadius=" << chosenFaceRadius
                << " directionalPenalty=" << chosenDirectionalPenalty
                << endl;
        }
    }
    else
    {
        label patchFaceI = faceI - cpp.start();
        receivePatchFaceI = cpp.pointFace(patchFaceI, direction, searchPosition);

        if (receivePatchFaceI < 0)
        {
            // 非严格匹配曲面接口下，单一搜索方向可能因为局部法向、角点容差
            // 或粒子速度近似切向而暂时找不到合法目标面。这里不立刻报错退出，
            // 而是把“本次方向失败”的结果返回给派生类，让它继续尝试其他
            // 候选搜索方向，提高 AMI 目标面搜索的稳健性。
            return false;
        }

        // 串行路径仍然走 OpenFOAM 原生 pointFace()。
        // 它没有 strict / loose 的两级回退概念，因此这里用 native 标记，
        // 让上层统计能区分“串行原生搜索”和“并行本地 strict / loose 搜索”。
        if (searchModeLabel)
        {
            *searchModeLabel = "native";
        }
    }

    // 并行局部搜索当前采用的是“局部区块近邻插入”近似策略：
    // 先选对目标面，再把粒子直接插入到目标 owner cell 内部一个有限深度。
    // 这样做不追求面内落点的高精度，而是优先保证两件事：
    // 1. 粒子进入正确目标面所属的接收单元；
    // 2. 粒子不会仍然贴在目标 patch 上，从而在同一步里重复命中同一周期面。
    const bool approximateFaceCentreInsertion = Pstream::parRun();

    const point positionBeforeTransfer = p.position();

    p.face() = receivePatchFaceI + receivePatch.start();
    p.cell() = mesh_.faceOwner()[p.face()];
    p.tetFace() = p.face();
    p.tetPt() = mesh_.faces()[p.tetFace()].size() - 1 - p.tetPt();

    const point& targetFaceCentre = mesh_.faceCentres()[p.face()];
    const point& ownerCellCentre = mesh_.cellCentres()[p.cell()];
    const vector targetFaceArea = mesh_.faceAreas()[p.face()];
    point insertionBasePosition = targetFaceCentre;

    // 对 stage-interface AMI 来说，评分参考点和最终插入参考点不应再混为一谈：
    // 1. targetReferencePosition     用于“哪一张目标面更像正确答案”的评分；
    // 2. targetInsertionReferencePosition 用于“粒子最终应该从目标面上的哪里入内”。
    //
    // 若仍然无条件从 targetFaceCentre 起步，就会把目标 face 的面内跨度
    // 直接转成一个额外的 Omega x r 偏差，从而在 summary 里形成
    // meanFinalPositionDeltaUDifference 的下限。
    //
    // 这里采用一个更稳妥的折中：
    // 1. 先把参考落点正交投影到目标面所在平面；
    // 2. 再把面内偏移限制在局部面尺度附近；
    // 3. 最后只沿入内方向推进一个很小距离。
    //
    // 这样做不追求“精确面内几何交点”，但可以显著减弱 targetFaceCentre
    // 带来的系统性面内偏置，同时保持粒子仍然稳定落在目标 owner cell 一侧。
    if (targetInsertionReferencePosition && mag(targetFaceArea) > ROOTVSMALL)
    {
        const vector unitFaceNormal = targetFaceArea/mag(targetFaceArea);
        point projectedInsertionBasePosition =
            *targetInsertionReferencePosition
          - (((*targetInsertionReferencePosition) - targetFaceCentre) & unitFaceNormal)
            * unitFaceNormal;

        vector tangentialOffset = projectedInsertionBasePosition - targetFaceCentre;
        tangentialOffset -= (tangentialOffset & unitFaceNormal)*unitFaceNormal;

        const scalar faceLengthScale = Foam::sqrt(mag(targetFaceArea));
        const scalar maxTangentialOffset = max(ROOTVSMALL, 0.45*faceLengthScale);
        const scalar tangentialOffsetMag = mag(tangentialOffset);

        if (tangentialOffsetMag > maxTangentialOffset)
        {
            projectedInsertionBasePosition =
                targetFaceCentre
              + (maxTangentialOffset/tangentialOffsetMag)*tangentialOffset;
        }

        insertionBasePosition = projectedInsertionBasePosition;
    }

    // 串行 AMI 路径也必须先把“目标侧参考落点”真正写回粒子位置。
    // 否则像 blade-row stage-interface 这种会先做几何映射的模型，虽然上面
    // 已经算出了 insertionBasePosition，下面的 inward push 仍会从旧的源侧
    // crossing position 出发，最终把 requested owner cell 与真实几何位置
    // 拉开，触发 initCellFacePt() 在目标 cell 内找不到 tetFace/tetPt。
    if (targetInsertionReferencePosition)
    {
        p.position() = insertionBasePosition;
    }

    if (approximateFaceCentreInsertion)
    {
        // 上一版把粒子直接放到 owner cell 内部 10% 深度，虽然能消掉重复撞面，
        // 但会把粒子系统性偏向 owner cell centre 一侧，从而明显扭曲跨面通量方向。
        // 这一版改成更温和的“面法向有限深度插入”：
        // 1. 仍从目标面中心出发；
        // 2. 只沿目标面的入内方向推进一个与面尺度成比例的小距离；
        // 3. 不再按 face-centre 到 owner-centre 的长连线大比例推进。
        //
        // 这样做的目标是同时满足：
        //     避免重复撞同一周期面
        // 和
        //     避免 owner-cell bias 过强地扭曲主导跨面方向
        const vector faceToOwner = ownerCellCentre - insertionBasePosition;
        vector insertionDirection = faceToOwner;

        if (mag(insertionDirection) <= ROOTVSMALL)
        {
            // 极端退化时退回到面法向的反方向，保证仍能朝单元内部推进。
            const vector faceNormal = targetFaceArea;

            if (mag(faceNormal) > ROOTVSMALL)
            {
                insertionDirection = -faceNormal;
            }
        }

        if (mag(insertionDirection) > ROOTVSMALL)
        {
            const scalar faceLengthScale = Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));

            // 插入深度取面尺度的一个更小比例：
            //     x_new = x_face + beta * L_f * n_in
            // 其中 beta 从上一版的 0.05 进一步减到 0.02。
            // 这样做的目的，是在保持“不会贴在 patch 上重复撞面”的同时，
            // 进一步减弱插入位置对 owner-cell 一侧的偏置，避免把跨面主方向拉偏。
            const scalar insertionDistance = max(ROOTVSMALL, 2e-2*faceLengthScale);
            p.position() =
                insertionBasePosition
              + (insertionDistance/mag(insertionDirection))*insertionDirection;
        }
        else
        {
            p.position() = insertionBasePosition;
        }
    }

    // AMI 跨面后同样可能出现粒子几何上仍然贴在接收侧 patch 上的病态。
    // 即便上面已经把粒子沿入内方向推进了一个有限距离，仍保留一层极小 inward push
    // 作为兜底，避免极端退化几何下再次被识别为贴面粒子。
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
        const scalar faceLengthScale = Foam::sqrt(mag(mesh_.faceAreas()[p.face()]));
        const scalar pushDistance = max(ROOTVSMALL, 1e-6*faceLengthScale);

        p.position() += (pushDistance/inwardMag)*inwardDirection;
    }

    // AMI 路径会在目标面中心近似插入和 inward push 后改写粒子几何位置。
    // 因此同样需要把 tetFace/tetPt 与最终位置重新对齐，避免后续追踪沿用
    // 过期的 boundary-face 局部拓扑。
    if (mag(p.position() - positionBeforeTransfer) > ROOTVSMALL)
    {
        p.initCellFacePt();
    }

    return true;
}


void dsmcCyclicBoundary::applyStandardCyclicAMITransfer
(
    dsmcParcel& p,
    const vector& direction,
    word* searchModeLabel,
    const point* targetReferencePosition,
    const point* targetInsertionReferencePosition,
    const vector* rotationAxis,
    const point* rotationCentre
) const
{
    if
    (
        !tryStandardCyclicAMITransfer
        (
            p,
            direction,
            searchModeLabel,
            targetReferencePosition,
            targetInsertionReferencePosition,
            rotationAxis,
            rotationCentre
        )
    )
    {
        const label faceI = p.face();
        const label patchI = mesh_.boundaryMesh().whichPatch(faceI);
        const cyclicAMIPolyPatch& cpp =
            refCast<const cyclicAMIPolyPatch>(mesh_.boundaryMesh()[patchI]);
        const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();

        // 这里的包装接口只负责在“失败即报错”的调用点补充更清晰的报错信息。
        // sourcePatch/receivePatch 是 tryStandardCyclicAMITransfer() 内部的局部变量，
        // 在当前作用域中并不存在，因此这里回退到本函数本地已经拿到的
        // cpp / receiveCpp 名称，既能表达当前跨面的周期面对，也不会再触发作用域错误。
        FatalErrorInFunction
            << "Unable to locate target AMI face while transferring parcel across "
            << cpp.name() << " and " << receiveCpp.name()
            << " for parcel position " << p.position()
            << " using direction " << direction
            << exit(FatalError);
    }
}


void dsmcCyclicBoundary::updateBoundaryProperties
(
    const dictionary& newDict
)
{
    boundaryDict_ = newDict.subDict("cyclicBoundaryProperties");
    refreshCoupledPatchState();
}


void dsmcCyclicBoundary::setNewBoundaryFields()
{
    refreshCoupledPatchState();
}


void dsmcCyclicBoundary::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const scalarField& yData
)
{
    fileName writeFile(pathName/nameFile);

    graph outputGraph("title", "x", "y", xData, yData);

    outputGraph.write(writeFile, "raw");
}


void dsmcCyclicBoundary::writeTimeData
(
    const fileName& pathName,
    const word& nameFile,
    const scalarField& xData,
    const vectorField& yData
)
{
    OFstream file(pathName/nameFile + ".xyz");

    if(file.good())
    {
        forAll(yData, n)
        {
            file
                << xData[n] << "\t"
                << yData[n].x() << "\t" << yData[n].y()
                << "\t" << yData[n].z()
                << endl;
        }
    }
    else
    {
        FatalErrorIn("void patchBoundary::writeTimeData()")
            << "Cannot open file " << file.name()
            << abort(FatalError);
    }
}

const labelList& dsmcCyclicBoundary::controlPatch() const
{
    return coupledFacesA_;
}

const labelList& dsmcCyclicBoundary::controlZone() const
{
    return cellsA_;
}

const labelList& dsmcCyclicBoundary::allFaces() const
{
    return faces_;
}

const word& dsmcCyclicBoundary::patchName() const
{
    return patchName_;
}

const label& dsmcCyclicBoundary::patchId() const
{
    return patchId_;
}

const scalarField& dsmcCyclicBoundary::densityField() const
{
    return densities_;
}

scalarField& dsmcCyclicBoundary::densityField()
{
    return densities_;
}

const vectorField& dsmcCyclicBoundary::velocityField() const
{
    return velocities_;
}
vectorField& dsmcCyclicBoundary::velocityField()
{
    return velocities_;
}

const scalarField& dsmcCyclicBoundary::temperatureField() const
{
    return temperatures_;
}

scalarField& dsmcCyclicBoundary::temperatureField()
{
    return temperatures_;
}

const tensor& dsmcCyclicBoundary::strainRate() const
{
    return strainRate_;
}

tensor& dsmcCyclicBoundary::strainRate()
{
    return strainRate_;
}


const bool& dsmcCyclicBoundary::writeInTimeDir() const
{
    return writeInTimeDir_;
}

const bool& dsmcCyclicBoundary::writeInCase() const
{
    return writeInCase_;
}


} // End namespace Foam

// ************************************************************************* //
