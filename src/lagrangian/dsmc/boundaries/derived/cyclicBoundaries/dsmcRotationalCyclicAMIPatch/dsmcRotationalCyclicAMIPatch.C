/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "dsmcRotationalCyclicAMIPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transform.H"
#include "OFstream.H"
#include "ListOps.H"

namespace Foam
{

namespace
{
    // 对调试量 oldRadius / oldTheta / dThetaAfterInsertion 来说，
    // 关键不是“固定在 y-z 平面算角度”，而是“围绕当前配置的 rotationAxis_”
    // 在其法平面内定义一致的局部极坐标。
    // 若继续把 theta 写死成 atan2(z, y)，那么把算例从绕 x 轴切到绕 z 轴后，
    // 日志会系统性地给出错误半径和错误角度，直接掩盖真正的几何 bug。
    vector azimuthReferenceDirection(const vector& rotationAxis)
    {
        // 优先用 x 轴作为平面内参考方向；若它几乎与旋转轴平行，
        // 则退回到 y 轴；若仍退化，再退回到 z 轴。
        vector seed(1, 0, 0);

        if (mag(seed - (seed & rotationAxis)*rotationAxis) <= ROOTVSMALL)
        {
            seed = vector(0, 1, 0);
        }

        if (mag(seed - (seed & rotationAxis)*rotationAxis) <= ROOTVSMALL)
        {
            seed = vector(0, 0, 1);
        }

        vector e1 = seed - (seed & rotationAxis)*rotationAxis;
        const scalar e1Mag = mag(e1);

        if (e1Mag <= ROOTVSMALL)
        {
            return vector(1, 0, 0);
        }

        return e1/e1Mag;
    }


    scalar cylindricalRadiusAboutAxis
    (
        const vector& relativePosition,
        const vector& rotationAxis
    )
    {
        const scalar axial = relativePosition & rotationAxis;
        const vector radial = relativePosition - axial*rotationAxis;

        return mag(radial);
    }


    scalar azimuthAngleAboutAxis
    (
        const vector& relativePosition,
        const vector& rotationAxis
    )
    {
        const scalar axial = relativePosition & rotationAxis;
        const vector radial = relativePosition - axial*rotationAxis;
        const scalar radialMag = mag(radial);

        if (radialMag <= ROOTVSMALL)
        {
            return 0.0;
        }

        const vector e1 = azimuthReferenceDirection(rotationAxis);
        const vector e2 = rotationAxis ^ e1;

        return Foam::atan2(radial & e2, radial & e1);
    }
}

defineTypeNameAndDebug(dsmcRotationalCyclicAMIPatch, 0);

addToRunTimeSelectionTable
(
    dsmcCyclicBoundary,
    dsmcRotationalCyclicAMIPatch,
    dictionary
);


void dsmcRotationalCyclicAMIPatch::readProperties()
{
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(1, 0, 0));
    rotationCentre_ =
        propsDict_.lookupOrDefault<vector>("rotationCentre", vector::zero);
    rotationAngle_ =
        propsDict_.lookupOrDefault<scalar>("rotationAngle", 0.0);
    transformPosition_ =
        propsDict_.lookupOrDefault<bool>("transformPosition", true);
    debug_ =
        propsDict_.lookupOrDefault<bool>("debug", false);
    maxDebugPrints_ =
        propsDict_.lookupOrDefault<label>("maxDebugPrints", 20);

    const scalar axisMag = mag(rotationAxis_);

    if (axisMag <= ROOTVSMALL)
    {
        FatalErrorInFunction
            << "rotationAxis magnitude must be non-zero for patch "
            << patchName_
            << exit(FatalError);
    }

    rotationAxis_ /= axisMag;
}


vector dsmcRotationalCyclicAMIPatch::rotateVectorRodrigues
(
    const vector& value,
    const scalar signedAngle
) const
{
    // rotationAxis_ 在 readProperties() 中已经被单位化，这里直接把它当作
    // Rodrigues 公式中的单位旋转轴 a_hat 使用。
    const vector& axis = rotationAxis_;

    // Rodrigues 公式的三个部分分别是：
    // 1. 原向量在垂直于旋转轴平面内的 cos(theta) 投影；
    // 2. 由 a_hat x v 给出的正交旋转分量；
    // 3. 原向量在旋转轴方向上的不变量分量。
    const scalar cosTheta = Foam::cos(signedAngle);
    const scalar sinTheta = Foam::sin(signedAngle);
    const scalar axisDotValue = axis & value;
    const vector axisCrossValue = axis ^ value;

    return
        value*cosTheta
      + axisCrossValue*sinTheta
      + axis*axisDotValue*(1.0 - cosTheta);
}


point dsmcRotationalCyclicAMIPatch::rotatePointRodrigues
(
    const point& position,
    const scalar signedAngle
) const
{
    // 对绝对位置点必须先平移到以旋转中心为原点的局部坐标系，
    // 再对相对向量应用 Rodrigues 旋转，最后再平移回全局坐标系：
    //     p' = c + R(p - c)
    const vector relative = position - rotationCentre_;

    return rotationCentre_ + rotateVectorRodrigues(relative, signedAngle);
}


void dsmcRotationalCyclicAMIPatch::appendSearchDirection
(
    DynamicList<vector>& directions,
    DynamicList<word>& directionNames,
    const vector& candidate,
    const word& candidateName
) const
{
    // 非严格匹配 AMI 接口的目标面搜索对 direction 较敏感。
    // 这里先过滤掉零向量或极小向量，因为它们既不能提供稳定的
    // 几何判别方向，也容易把 pointFace() 推进到数值退化分支。
    const scalar candidateMag = mag(candidate);

    if (candidateMag <= ROOTVSMALL)
    {
        return;
    }

    const vector unitCandidate = candidate/candidateMag;

    // 只去重“同向”候选，不去掉反向候选。
    // 原因是曲面非匹配接口中，源侧法向的正反方向都可能成为有用的
    // 回退搜索方向；但若两个候选只是几乎完全同向，再重复尝试没有意义。
    forAll(directions, directionI)
    {
        if ((unitCandidate & directions[directionI]) > 1.0 - SMALL)
        {
            return;
        }
    }

    directions.append(unitCandidate);
    directionNames.append(candidateName);
}


void dsmcRotationalCyclicAMIPatch::recordCrossingEvent
(
    const word& directionLabel,
    const word& searchModeLabel,
    const label targetFace,
    const scalar absUn
)
{
    // 这里记录的是“真实发生且已经完成目标面拓扑更新”的跨面事件。
    // 因此 targetFace 使用的是迁移后的 p.face()，而不是旧的源侧面号。
    // 这样后续检查 face 命中偏置时，得到的就是接收侧真正被命中的面。
    //
    // 与上一版相比，这里新增 searchModeLabel，目的是把统计拆成 strict / loose
    // 两条路径，直接回答：当前剩余的宏观量偏差到底主要来自哪一种搜索模式。
    Map<label>* totalFaceHits = NULL;
    Map<label>* modeFaceHits = NULL;
    label* totalTransfersPtr = NULL;
    scalar* totalAbsUnPtr = NULL;
    label* modeTransfersPtr = NULL;
    scalar* modeAbsUnPtr = NULL;

    if (directionLabel == "AtoB")
    {
        totalTransfersPtr = &totalTransfersAtoB_;
        totalAbsUnPtr = &totalAbsUnAtoB_;
        totalFaceHits = &targetFaceHitsAtoB_;

        if (searchModeLabel == "strict")
        {
            modeTransfersPtr = &totalTransfersAtoBStrict_;
            modeAbsUnPtr = &totalAbsUnAtoBStrict_;
            modeFaceHits = &targetFaceHitsAtoBStrict_;
        }
        else if (searchModeLabel == "loose")
        {
            modeTransfersPtr = &totalTransfersAtoBLoose_;
            modeAbsUnPtr = &totalAbsUnAtoBLoose_;
            modeFaceHits = &targetFaceHitsAtoBLoose_;
        }
    }
    else if (directionLabel == "BtoA")
    {
        totalTransfersPtr = &totalTransfersBtoA_;
        totalAbsUnPtr = &totalAbsUnBtoA_;
        totalFaceHits = &targetFaceHitsBtoA_;

        if (searchModeLabel == "strict")
        {
            modeTransfersPtr = &totalTransfersBtoAStrict_;
            modeAbsUnPtr = &totalAbsUnBtoAStrict_;
            modeFaceHits = &targetFaceHitsBtoAStrict_;
        }
        else if (searchModeLabel == "loose")
        {
            modeTransfersPtr = &totalTransfersBtoALoose_;
            modeAbsUnPtr = &totalAbsUnBtoALoose_;
            modeFaceHits = &targetFaceHitsBtoALoose_;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown crossing direction label " << directionLabel
            << " for patch " << patchName_
            << exit(FatalError);
    }

    ++(*totalTransfersPtr);
    (*totalAbsUnPtr) += absUn;

    if (totalFaceHits->found(targetFace))
    {
        (*totalFaceHits)[targetFace] += 1;
    }
    else
    {
        totalFaceHits->insert(targetFace, 1);
    }

    // strict / loose 之外的模式（例如串行 native 路径）目前只计入总累计。
    // 这能保证旧的串行结果文件格式仍然可解释，同时把当前并行研究关心的
    // strict / loose 分叉清晰拆出来。
    if (modeTransfersPtr && modeAbsUnPtr && modeFaceHits)
    {
        ++(*modeTransfersPtr);
        (*modeAbsUnPtr) += absUn;

        if (modeFaceHits->found(targetFace))
        {
            (*modeFaceHits)[targetFace] += 1;
        }
        else
        {
            modeFaceHits->insert(targetFace, 1);
        }
    }
}


void dsmcRotationalCyclicAMIPatch::writeFaceHitTable
(
    const fileName& outputFile,
    const Map<label>& faceHits,
    const word& directionLabel
) const
{
    OFstream os(outputFile);

    os  << "# targetFace,cumulativeHits,direction" << nl;

    DynamicList<label> targetFaces(faceHits.size());

    forAllConstIter(Map<label>, faceHits, iter)
    {
        targetFaces.append(iter.key());
    }

    sort(targetFaces);

    forAll(targetFaces, faceI)
    {
        const label targetFace = targetFaces[faceI];
        os  << targetFace << ","
            << faceHits[targetFace] << ","
            << directionLabel << nl;
    }
}


dsmcRotationalCyclicAMIPatch::dsmcRotationalCyclicAMIPatch
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
    rotationAngle_(0.0),
    transformPosition_(true),
    debug_(false),
    maxDebugPrints_(20),
    debugPrintCount_(0),
    totalTransfersAtoB_(0),
    totalTransfersBtoA_(0),
    totalAbsUnAtoB_(0.0),
    totalAbsUnBtoA_(0.0),
    totalTransfersAtoBStrict_(0),
    totalTransfersAtoBLoose_(0),
    totalTransfersBtoAStrict_(0),
    totalTransfersBtoALoose_(0),
    totalAbsUnAtoBStrict_(0.0),
    totalAbsUnAtoBLoose_(0.0),
    totalAbsUnBtoAStrict_(0.0),
    totalAbsUnBtoALoose_(0.0),
    writtenTransfersAtoB_(0),
    writtenTransfersBtoA_(0),
    writtenAbsUnAtoB_(0.0),
    writtenAbsUnBtoA_(0.0),
    writtenTransfersAtoBStrict_(0),
    writtenTransfersAtoBLoose_(0),
    writtenTransfersBtoAStrict_(0),
    writtenTransfersBtoALoose_(0),
    writtenAbsUnAtoBStrict_(0.0),
    writtenAbsUnAtoBLoose_(0.0),
    writtenAbsUnBtoAStrict_(0.0),
    writtenAbsUnBtoALoose_(0.0),
    targetFaceHitsAtoB_(),
    targetFaceHitsBtoA_(),
    targetFaceHitsAtoBStrict_(),
    targetFaceHitsAtoBLoose_(),
    targetFaceHitsBtoAStrict_(),
    targetFaceHitsBtoALoose_()
{
    writeInTimeDir_ = true;
    writeInCase_ = false;
    readProperties();
}


dsmcRotationalCyclicAMIPatch::~dsmcRotationalCyclicAMIPatch()
{}


void dsmcRotationalCyclicAMIPatch::calculateProperties()
{}


void dsmcRotationalCyclicAMIPatch::initialConfiguration()
{}


void dsmcRotationalCyclicAMIPatch::controlMol
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    (void)td;

    const label faceI = p.face();
    const vector oldPosition = p.position();
    const vector oldVelocity = p.U();
    const vector oldRel = oldPosition - rotationCentre_;
    word side("none");
    word crossingDirection("none");

    scalar signedAngle = 0.0;

    const label sourceIndexA = findIndex(coupledFacesA_, faceI);
    const label sourceIndexB = findIndex(coupledFacesB_, faceI);

    if (sourceIndexA != -1)
    {
        // 对 rotational cyclic patch，rotationAngle_ 的物理含义必须固定成：
        // “把 periodicA 所在扇区旋转到 periodicB 所在扇区所需的正向角度”。
        //
        // 以当前 z 轴回归算例为例：
        // 1. periodicA 位于负 theta 一侧；
        // 2. periodicB 位于正 theta 一侧；
        // 3. 因此 A -> B 的正确 Rodrigues 旋转应该是 +rotationAngle_；
        // 4. 反向 B -> A 的正确 Rodrigues 旋转才是 -rotationAngle_。
        //
        // 旧实现把这里的符号写反了，会导致：
        // 1. 位置最终被 AMI 搜索插到目标 patch 一侧；
        // 2. 但速度却按反方向角度旋转；
        // 3. 最终在宏观平均场里表现成明显的径向偏置，而不是局部切向速度。
        signedAngle = rotationAngle_;
        side = "A";
        crossingDirection = "AtoB";
    }
    else if (sourceIndexB != -1)
    {
        signedAngle = -rotationAngle_;
        side = "B";
        crossingDirection = "BtoA";
    }
    else
    {
        return;
    }

    // 这里显式使用 Rodrigues 旋转公式，而不再依赖黑盒旋转张量。
    // 这样做的目的是把“位置点旋转”和“纯向量旋转”都写成同一套清晰公式，
    // 便于后续逐项核对：
    //     p' = c + R(p - c)
    //     U' = R U
    // 同时也更容易和日志中的 dTheta、理论旋转角做一一对应。
    const vector analyticalRotatedVelocity =
        rotateVectorRodrigues(oldVelocity, signedAngle);
    const point analyticalRotatedPosition =
        rotatePointRodrigues(oldPosition, signedAngle);

    // 旋转后的源侧面中心是当前粒子“应当落到目标周期面哪一块区域”的
    // 更稳定几何参考。与粒子瞬时位置相比，它不容易受到贴面推进和
    // 局部 processor 片段截断的影响，因此后续并行目标面搜索会同时参考它。
    const point rotatedSourceFaceCentre =
        rotatePointRodrigues(mesh_.faceCentres()[faceI], signedAngle);

    p.U() = analyticalRotatedVelocity;

    if (transformPosition_)
    {
        // 对几何位置先施加解析旋转，把粒子送到目标扇区对应的几何位置附近。
        // 随后的 AMI 传递将以该旋转后位置为输入，在目标 patch 上搜索合法面元，
        // 解决曲面非匹配条件下无法直接 face-to-face 映射的问题。
        p.position() = analyticalRotatedPosition;
    }

    // 对非严格匹配的曲面 AMI 接口，单一搜索方向并不总是足够稳健。
    // 当粒子速度接近切向、局部曲率较大或 patch 交线附近存在容差退化时，
    // 仅靠“旋转后的速度方向”调用 pointFace() 可能会偶发找不到目标面。
    // 因此这里构造一组有序候选方向，按“最物理 -> 最保守”的顺序回退尝试：
    // 1. 旋转后的粒子速度；
    // 2. 旋转后的源侧面法向；
    // 3. 其反向；
    // 4. 旋转后的源侧单元向外位移方向；
    // 5. 其反向。
    // 这样做的目的不是改变物理映射公式，而是在不改变
    //     r_new = r_c + R(±Δθ)(r_old - r_c)
    //     v_new = R(±Δθ) v_old
    // 的前提下，提高 AMI 目标面搜索的数值稳健性。
    DynamicList<vector> searchDirections(0);
    DynamicList<word> searchDirectionNames(0);

    const vector rotatedVelocity = p.U();
    const vector rotatedFaceNormal =
        rotateVectorRodrigues(mesh_.faceAreas()[faceI], signedAngle);
    const vector rotatedOutwardVector =
        rotateVectorRodrigues(oldPosition - mesh_.cellCentres()[p.cell()], signedAngle);

    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        rotatedVelocity,
        "rotatedVelocity"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        rotatedFaceNormal,
        "rotatedFaceNormal"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        -rotatedFaceNormal,
        "minusRotatedFaceNormal"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        rotatedOutwardVector,
        "rotatedOutwardVector"
    );
    appendSearchDirection
    (
        searchDirections,
        searchDirectionNames,
        -rotatedOutwardVector,
        "minusRotatedOutwardVector"
    );

    word successfulDirection("none");
    word successfulSearchMode("none");
    bool transferSucceeded = false;

    forAll(searchDirections, directionI)
    {
        if
        (
            tryStandardCyclicAMITransfer
            (
                p,
                searchDirections[directionI],
                &successfulSearchMode,
                &rotatedSourceFaceCentre,
                // 参数顺序必须与基类 helper 保持一致：
                // targetReferencePosition 用于目标面评分，targetInsertionReferencePosition
                // 用于最终插入点投影，随后才是 rotationAxis / rotationCentre。
                // 若漏掉插入参考点，z 轴算例会把旋转轴和旋转中心实参错位，
                // 目标面评分退回错误的默认几何语义。
                &analyticalRotatedPosition,
                &rotationAxis_,
                &rotationCentre_
            )
        )
        {
            successfulDirection = searchDirectionNames[directionI];
            transferSucceeded = true;
            break;
        }
    }

    if (!transferSucceeded)
    {
        FatalErrorInFunction
            << "Unable to complete cyclicAMI transfer for patch " << patchName_
            << " side " << side
            << " face " << faceI
            << " after trying search directions " << searchDirectionNames
            << ". Parcel position: " << p.position()
            << ", parcel velocity: " << p.U()
            << exit(FatalError);
    }

    scalar absUn = 0.0;
    const scalar rotatedFaceNormalMag = mag(rotatedFaceNormal);

    if (rotatedFaceNormalMag > ROOTVSMALL)
    {
        const vector unitRotatedFaceNormal = rotatedFaceNormal/rotatedFaceNormalMag;
        absUn = mag(rotatedVelocity & unitRotatedFaceNormal);
    }

    recordCrossingEvent(crossingDirection, successfulSearchMode, p.face(), absUn);

    if (debug_ && debugPrintCount_ < maxDebugPrints_)
    {
        const vector newPosition = p.position();
        const vector newVelocity = p.U();
        const vector newRel = newPosition - rotationCentre_;
        const vector analyticalRel = analyticalRotatedPosition - rotationCentre_;

        // 这里必须围绕 rotationAxis_ 计算柱坐标半径和极角。
        // 前一版把调试量写死成 y-z 平面，只适用于绕 x 轴的算例。
        // 一旦 rotationAxis_ 改成 z 轴，oldRadius / dThetaAfterInsertion
        // 会直接变成错误值，无法再用日志判断“目标面插入到底有没有带偏”。
        const scalar oldRadius =
            cylindricalRadiusAboutAxis(oldRel, rotationAxis_);
        const scalar analyticalRadius =
            cylindricalRadiusAboutAxis(analyticalRel, rotationAxis_);
        const scalar newRadius =
            cylindricalRadiusAboutAxis(newRel, rotationAxis_);
        const scalar oldTheta = azimuthAngleAboutAxis(oldRel, rotationAxis_);
        const scalar analyticalTheta =
            azimuthAngleAboutAxis(analyticalRel, rotationAxis_);
        const scalar newTheta = azimuthAngleAboutAxis(newRel, rotationAxis_);
        scalar dTheta = analyticalTheta - oldTheta;
        scalar dThetaAfterInsertion = newTheta - oldTheta;

        while (dTheta > constant::mathematical::pi)
        {
            dTheta -= constant::mathematical::twoPi;
        }

        while (dTheta < -constant::mathematical::pi)
        {
            dTheta += constant::mathematical::twoPi;
        }

        while (dThetaAfterInsertion > constant::mathematical::pi)
        {
            dThetaAfterInsertion -= constant::mathematical::twoPi;
        }

        while (dThetaAfterInsertion < -constant::mathematical::pi)
        {
            dThetaAfterInsertion += constant::mathematical::twoPi;
        }

        Info<< "dsmcRotationalCyclicAMIPatch debug: patch=" << patchName_
            << " side=" << side
            << " face=" << faceI
            << " angle=" << signedAngle
            << " oldPosition=" << oldPosition
            << " analyticalRotatedPosition=" << analyticalRotatedPosition
            << " newPosition=" << newPosition
            << " oldU=" << oldVelocity
            << " analyticalRotatedU=" << analyticalRotatedVelocity
            << " newU=" << newVelocity
            << " oldRadius=" << oldRadius
            << " analyticalRadius=" << analyticalRadius
            << " newRadius=" << newRadius
            << " oldTheta=" << oldTheta
            << " analyticalTheta=" << analyticalTheta
            << " newTheta=" << newTheta
            << " dTheta=" << dTheta
            << " dThetaAfterInsertion=" << dThetaAfterInsertion
            << " transformPosition=" << transformPosition_
            << " successfulSearchDirection=" << successfulDirection
            << " successfulSearchMode=" << successfulSearchMode
            << nl;

        ++debugPrintCount_;
    }
}


void dsmcRotationalCyclicAMIPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    (void)fixedPathName;

    // dsmcBoundaries 已经为当前模型和当前 patch 创建好：
    //   case/<time>/uniform/boundaries/dsmc/cyclicBoundaryModels/<model>/<patch>
    // 这里直接在该目录下写 summary 与 face-hit CSV，方便和其他 DSMC
    // 输出保持同一目录层级，不额外引入新的分散路径。
    const fileName summaryFile(timePath/"amiCrossingSummary.dat");
    const fileName faceHitsAtoBFile(timePath/"amiCrossingFaceHits_AtoB.csv");
    const fileName faceHitsBtoAFile(timePath/"amiCrossingFaceHits_BtoA.csv");
    const fileName faceHitsAtoBStrictFile(timePath/"amiCrossingFaceHits_AtoB_strict.csv");
    const fileName faceHitsAtoBLooseFile(timePath/"amiCrossingFaceHits_AtoB_loose.csv");
    const fileName faceHitsBtoAStrictFile(timePath/"amiCrossingFaceHits_BtoA_strict.csv");
    const fileName faceHitsBtoALooseFile(timePath/"amiCrossingFaceHits_BtoA_loose.csv");

    const label intervalTransfersAtoB = totalTransfersAtoB_ - writtenTransfersAtoB_;
    const label intervalTransfersBtoA = totalTransfersBtoA_ - writtenTransfersBtoA_;
    const scalar intervalAbsUnAtoB = totalAbsUnAtoB_ - writtenAbsUnAtoB_;
    const scalar intervalAbsUnBtoA = totalAbsUnBtoA_ - writtenAbsUnBtoA_;

    const label intervalTransfersAtoBStrict =
        totalTransfersAtoBStrict_ - writtenTransfersAtoBStrict_;
    const label intervalTransfersAtoBLoose =
        totalTransfersAtoBLoose_ - writtenTransfersAtoBLoose_;
    const label intervalTransfersBtoAStrict =
        totalTransfersBtoAStrict_ - writtenTransfersBtoAStrict_;
    const label intervalTransfersBtoALoose =
        totalTransfersBtoALoose_ - writtenTransfersBtoALoose_;
    const scalar intervalAbsUnAtoBStrict =
        totalAbsUnAtoBStrict_ - writtenAbsUnAtoBStrict_;
    const scalar intervalAbsUnAtoBLoose =
        totalAbsUnAtoBLoose_ - writtenAbsUnAtoBLoose_;
    const scalar intervalAbsUnBtoAStrict =
        totalAbsUnBtoAStrict_ - writtenAbsUnBtoAStrict_;
    const scalar intervalAbsUnBtoALoose =
        totalAbsUnBtoALoose_ - writtenAbsUnBtoALoose_;

    OFstream os(summaryFile);

    os  << "# MCOR1_0 rotational cyclic AMI crossing summary" << nl
        << "# patchName cumulativeTransfersAtoB cumulativeTransfersBtoA "
        << "intervalTransfersAtoB intervalTransfersBtoA "
        << "cumulativeAbsUnAtoB cumulativeAbsUnBtoA "
        << "intervalAbsUnAtoB intervalAbsUnBtoA "
        << "cumulativeTransfersAtoBStrict cumulativeTransfersAtoBLoose "
        << "cumulativeTransfersBtoAStrict cumulativeTransfersBtoALoose "
        << "intervalTransfersAtoBStrict intervalTransfersAtoBLoose "
        << "intervalTransfersBtoAStrict intervalTransfersBtoALoose "
        << "cumulativeAbsUnAtoBStrict cumulativeAbsUnAtoBLoose "
        << "cumulativeAbsUnBtoAStrict cumulativeAbsUnBtoALoose "
        << "intervalAbsUnAtoBStrict intervalAbsUnAtoBLoose "
        << "intervalAbsUnBtoAStrict intervalAbsUnBtoALoose" << nl
        << patchName_ << " "
        << totalTransfersAtoB_ << " "
        << totalTransfersBtoA_ << " "
        << intervalTransfersAtoB << " "
        << intervalTransfersBtoA << " "
        << totalAbsUnAtoB_ << " "
        << totalAbsUnBtoA_ << " "
        << intervalAbsUnAtoB << " "
        << intervalAbsUnBtoA << " "
        << totalTransfersAtoBStrict_ << " "
        << totalTransfersAtoBLoose_ << " "
        << totalTransfersBtoAStrict_ << " "
        << totalTransfersBtoALoose_ << " "
        << intervalTransfersAtoBStrict << " "
        << intervalTransfersAtoBLoose << " "
        << intervalTransfersBtoAStrict << " "
        << intervalTransfersBtoALoose << " "
        << totalAbsUnAtoBStrict_ << " "
        << totalAbsUnAtoBLoose_ << " "
        << totalAbsUnBtoAStrict_ << " "
        << totalAbsUnBtoALoose_ << " "
        << intervalAbsUnAtoBStrict << " "
        << intervalAbsUnAtoBLoose << " "
        << intervalAbsUnBtoAStrict << " "
        << intervalAbsUnBtoALoose << nl;

    writeFaceHitTable(faceHitsAtoBFile, targetFaceHitsAtoB_, "AtoB");
    writeFaceHitTable(faceHitsBtoAFile, targetFaceHitsBtoA_, "BtoA");
    writeFaceHitTable(faceHitsAtoBStrictFile, targetFaceHitsAtoBStrict_, "AtoB_strict");
    writeFaceHitTable(faceHitsAtoBLooseFile, targetFaceHitsAtoBLoose_, "AtoB_loose");
    writeFaceHitTable(faceHitsBtoAStrictFile, targetFaceHitsBtoAStrict_, "BtoA_strict");
    writeFaceHitTable(faceHitsBtoALooseFile, targetFaceHitsBtoALoose_, "BtoA_loose");

    // 更新“上次写出快照”，保证下一次 output() 写出的 interval 值只对应
    // 新增时段，而不是重复统计从 t=0 到当前的全部累计量。
    writtenTransfersAtoB_ = totalTransfersAtoB_;
    writtenTransfersBtoA_ = totalTransfersBtoA_;
    writtenAbsUnAtoB_ = totalAbsUnAtoB_;
    writtenAbsUnBtoA_ = totalAbsUnBtoA_;
    writtenTransfersAtoBStrict_ = totalTransfersAtoBStrict_;
    writtenTransfersAtoBLoose_ = totalTransfersAtoBLoose_;
    writtenTransfersBtoAStrict_ = totalTransfersBtoAStrict_;
    writtenTransfersBtoALoose_ = totalTransfersBtoALoose_;
    writtenAbsUnAtoBStrict_ = totalAbsUnAtoBStrict_;
    writtenAbsUnAtoBLoose_ = totalAbsUnAtoBLoose_;
    writtenAbsUnBtoAStrict_ = totalAbsUnBtoAStrict_;
    writtenAbsUnBtoALoose_ = totalAbsUnBtoALoose_;
}


void dsmcRotationalCyclicAMIPatch::updateProperties(const dictionary& newDict)
{
    updateBoundaryProperties(newDict);
    propsDict_ = newDict.subDict(typeName + "Properties");
    readProperties();
}

} // End namespace Foam

// ************************************************************************* //
