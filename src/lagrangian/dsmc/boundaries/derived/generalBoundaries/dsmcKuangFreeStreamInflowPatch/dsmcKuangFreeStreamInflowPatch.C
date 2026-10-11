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

#include "dsmcKuangFreeStreamInflowPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcKuangFreeStreamInflowPatch, 0);

addToRunTimeSelectionTable
(
    dsmcGeneralBoundary, dsmcKuangFreeStreamInflowPatch, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcKuangFreeStreamInflowPatch::dsmcKuangFreeStreamInflowPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    typeIds_(),
    translationalTemperature_(),
    rotationalTemperature_(),
    vibrationalTemperature_(),
    electronicTemperature_(),
    numberDensities_(),
    accumulatedParcelsToInsert_(),
    omega_(0.0),
    centerPoint_(vector::zero),
    rotationAxis_(0, 0, 1),
    velocityFrame_("local"),
    velocity_(vector::zero)

    // angularVelocity_(dict.lookupOrDefault<vector>("angularVelocity", vector::zero))


{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcKuangFreeStreamInflowPatch::~dsmcKuangFreeStreamInflowPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcKuangFreeStreamInflowPatch::initialConfiguration()
{}

void dsmcKuangFreeStreamInflowPatch::calculateProperties()
{

}

vector dsmcKuangFreeStreamInflowPatch::reservoirVelocity
(
    const point& p,
    const label cellI
) const
{
    // 基础平移速度来自 boundariesDict 中的 velocity。
    // velocityFrame_ == local    时，它已经是当前 cell 内 parcel 存储速度；
    // velocityFrame_ == absolute 时，它是实验室绝对系速度，需要稍后转换。
    vector UModel = velocity_;

    // 当 omega 近似为零时，不需要计算旋转牵连速度。
    if (mag(omega_) > SMALL)
    {
        // 旋转轴长度用于把 rotationAxis_ 归一化为单位方向。
        const scalar axisMag = mag(rotationAxis_);

        // omega 非零时，旋转轴不能为零向量，否则 omega*axis x r 没有定义。
        if (axisMag <= SMALL)
        {
            FatalErrorInFunction
                << "Invalid rotationAxis for dsmcKuangFreeStreamInflowPatch."
                << nl
                << "patchName = " << patchName_ << nl
                << "omega = " << omega_ << nl
                << "rotationAxis = " << rotationAxis_ << nl
                << "A non-zero omega requires a non-zero rotationAxis because "
                << "the reservoir velocity is U = velocity + omega*axisHat x "
                << "(x - centerPoint)." << abort(FatalError);
        }

        // 归一化旋转轴，得到 axisHat。
        const vector axisHat = rotationAxis_/axisMag;

        // 相对旋转中心的位置矢量 r = x - x0。
        const vector r = p - centerPoint_;

        // 牵连速度公式：U_rot = omega * axisHat x r。
        UModel += omega_*(axisHat ^ r);
    }

    // local 模式保持历史 dsmcKuang 行为：UModel 直接作为 parcel 局部速度。
    if (velocityFrame_ == "local")
    {
        return UModel;
    }

    // absolute 模式参考 dsmcWangPressureInlet：
    // UModel 是实验室绝对系速度，进入旋转 cell 前必须转换成局部相对速度。
    return velocityToLocalFrame(p, UModel, cellI);
}

void dsmcKuangFreeStreamInflowPatch::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    const scalar sqrtPi = sqrt(pi);

    // compute parcels to insert
    forAll(accumulatedParcelsToInsert_, i)
    {
        const label typeId = typeIds_[i];
        scalar mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedParcelsToInsert_[i], f)
        {
            const label faceI = faces_[f];
            const label cellI = cells_[f];
            const vector sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT = cloud_.deltaTValue(cellI);

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    translationalTemperature_,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            // 当前模型允许 reservoir 具有位置相关的速度：
            //
            //     U_res(x) = velocity + omega*axisHat x (x - centerPoint)
            //
            // 若 velocityFrame=absolute，该速度先按实验室系构造，再参考
            // dsmcWangPressureInlet 的做法通过 velocityToLocalFrame() 转成
            // 当前旋转 cell 内 parcel 实际使用的局部速度。
            //
            // Bird 的半空间注入通量公式需要使用该面上的法向漂移速度。
            // 这里用面心速度代表整个边界面，保持注入粒子数计算和后续
            // 速度采样使用同一个局部 reservoir 速度。
            const vector faceReservoirVelocity =
                reservoirVelocity(mesh_.faceCentres()[faceI], cellI);

            // n = -sF/|sF| 指向计算域内部；sCosTheta 是无量纲法向漂移速度。
            scalar sCosTheta =
                (faceReservoirVelocity & -sF/fA)/mostProbableSpeed;

            // 固定网格自适应粒子数打开后，新注入 parcel 不再一定使用 RWF=1。
            // 单个新 parcel 实际代表的真实分子数为：
            //
            //     N_real_per_parcel = nParticles(patch, face) * pRWF(patch, face)
            //
            // 如果这里只除以 nParticles，而 addNewParcel 时又给新粒子赋
            // pRWF，则真实注入质量通量会被额外乘以 pRWF，形成
            // pRWF 越大 -> 注入真实质量越大 -> 压力/密度越高 -> pRWF 继续变大
            // 的正反馈。这正是 testTurboSingle 出口压力和
            // Reference-equivalent free particles 爆炸的直接原因。
            const scalar boundaryRWF =
                cloud_.coordSystem().pRWF(patchId_, f);

            // pRWF 必须为正，否则无法把真实通量折算为计算 parcel 数。
            if (boundaryRWF <= SMALL)
            {
                FatalErrorInFunction
                    << "Invalid boundary RWF for dsmcKuangFreeStreamInflowPatch."
                    << nl
                    << "patchName = " << patchName_ << nl
                    << "patchId = " << patchId_ << nl
                    << "local face index = " << f << nl
                    << "boundaryRWF = " << boundaryRWF << nl
                    << "Boundary injection requires pRWF > SMALL because "
                    << "the parcel insertion count is divided by "
                    << "nParticles(patch,face)*pRWF(patch,face)."
                    << abort(FatalError);
            }

            // 这是一个新计算 parcel 所代表的真实分子数。
            const scalar representedMoleculesPerInsertedParcel =
                cloud_.nParticles(patchId_, f)*boundaryRWF;

            // From Bird eqn 4.22
            accumulatedParcelsToInsert_[i][f] +=
                (
                    fA*numberDensities_[i]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi*representedMoleculesPerInsertedParcel);
        }
    }

    const scalar faceTranslationalTemperature = translationalTemperature_;
    const scalar faceRotationalTemperature = rotationalTemperature_;
    const scalar faceVibrationalTemperature = vibrationalTemperature_;
    const scalar faceElectronicTemperature = electronicTemperature_;

    // insert parcels
    forAll(faces_, f)
    {
        const label faceI = faces_[f];
        const label cellI = cells_[f];
        const vector fC = mesh_.faceCentres()[faceI];
        const vector sF = mesh_.faceAreas()[faces_[f]];

        scalar fA = mag(sF);

        List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
        (
            mesh_,
            faceI,
            cellI
        );

        // Cumulative triangle area fractions
        List<scalar> cTriAFracs(faceTets.size(), 0.0);

        scalar previousCumulativeSum = 0.0;

        forAll(faceTets, triI)
        {
            const tetIndices faceTetIs = faceTets[triI];

            cTriAFracs[triI] =
                faceTetIs.faceTri(mesh_).mag()/fA
                + previousCumulativeSum;

            previousCumulativeSum = cTriAFracs[triI];
        }

        // Force the last area fraction value to 1.0 to avoid any
        // rounding/non-flat face errors giving a value < 1.0
        cTriAFracs.last() = 1.0;

        // Normal unit vector *negative* so normal is pointing into the
        // domain
        vector n = sF;
        n /= -mag(n);

        //- Wall tangential unit vector. Use the direction between the
        //  face centre and the first vertex in the list
        vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
        t1 /= mag(t1);

        //- Other tangential unit vector.  Rescaling in case face is not
        //  flat and n and t1 aren't perfectly orthogonal
        vector t2 = n^t1;
        t2 /= mag(t2);

        forAll(typeIds_, m)
        {
            const label typeId = typeIds_[m];

            scalar& faceAccumulator = accumulatedParcelsToInsert_[m][f];

            // Number of whole particles to insert
            label nI = max(label(faceAccumulator), 0);

            // Add another particle with a probability proportional to the
            // remainder of taking the integer part of faceAccumulator
            if ((faceAccumulator - nI) > rndGen.sample01<scalar>())
            {
                nI++;
            }

            faceAccumulator -= nI;

            const scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nI; i++)
            {
                // Choose a triangle to insert on, based on their relative
                // area

                scalar triSelection = rndGen.sample01<scalar>();

                // Selected triangle
                label selectedTriI = -1;

                forAll(cTriAFracs, triI)
                {
                    selectedTriI = triI;

                    if (cTriAFracs[triI] >= triSelection)
                    {
                        break;
                    }
                }

                // Randomly distribute the points on the triangle.

                const tetIndices& faceTetIs = faceTets[selectedTriI];


                // 随机生成位置
                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // 使用新 parcel 实际插入点处的 reservoir 速度。
                // 对 velocityFrame=absolute，该函数会把实验室系速度转成局部速度。
                const vector localReservoirVelocity =
                    reservoirVelocity(p, cellI);

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTranslationalTemperature,
                        mass
                    )
                );

                scalar sCosTheta =
                    (localReservoirVelocity & n)/mostProbableSpeed;

                // Coefficients required for Bird eqn 12.5
                scalar uNormProbCoeffA = sCosTheta + sqrt(sqr(sCosTheta) + 2.0);

                scalar uNormProbCoeffB =
                    0.5*
                    (
                        1.0
                        + sCosTheta*(sCosTheta - sqrt(sqr(sCosTheta) + 2.0))
                    );

                // Equivalent to the QA value in Bird's DSMC3.FOR
                scalar randomScaling = 3.0;

                if (sCosTheta < -3)
                {
                    randomScaling = mag(sCosTheta) + 1;
                }

                scalar P = -1;

                // Normalised candidates for the normal direction velocity
                // component
                scalar uNormal;
                scalar uNormalThermal;

                if(abs(localReservoirVelocity & n) > VSMALL)
                {
                    // Select a velocity using Bird eqn 12.5
                    do
                    {
                        uNormalThermal =
                            randomScaling*(2.0*rndGen.sample01<scalar>() - 1);

                        uNormal = uNormalThermal + sCosTheta;

                        if (uNormal < 0.0)
                        {
                            P = -1;
                        }
                        else
                        {
                            P = 2.0*uNormal/uNormProbCoeffA
                                *exp(uNormProbCoeffB - sqr(uNormalThermal));
                        }

                    } while (P < rndGen.sample01<scalar>());
                }
                else
                {
                    uNormal = sqrt(-log(rndGen.sample01<scalar>()));
                }



                vector U =
                    sqrt(physicoChemical::k.value()*faceTranslationalTemperature/mass)
                    *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                    + (t1 & localReservoirVelocity)*t1
                    + (t2 & localReservoirVelocity)*t2
                    + mostProbableSpeed*uNormal*n
                    ;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceRotationalTemperature,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceVibrationalTemperature,
                    cloud_.constProps(typeId).nVibrationalModes(),
                    typeId
                );

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceElectronicTemperature,
                    cloud_.constProps(typeId).electronicDegeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList()
                );

                label newParcel = patchId();

                // 新 parcel 必须携带与上面通量折算一致的边界 RWF。
                // 这样插入 parcel 的数量减少了 pRWF 倍，但每个 parcel 的
                // 代表分子数增加了 pRWF 倍，边界真实通量保持 Bird 公式
                // 给定的物理值不变。
                const scalar RWF = cloud_.coordSystem().pRWF(patchId_, f);

                // 防御性检查：防止非法 RWF 被写入 parcel，进而污染质量守恒。
                if (RWF <= SMALL)
                {
                    FatalErrorInFunction
                        << "Invalid injection RWF for "
                        << "dsmcKuangFreeStreamInflowPatch." << nl
                        << "patchName = " << patchName_ << nl
                        << "patchId = " << patchId_ << nl
                        << "local face index = " << f << nl
                        << "cellI = " << cellI << nl
                        << "RWF = " << RWF << nl
                        << "Injected parcels must carry the same positive "
                        << "pRWF used in the boundary flux-to-parcel-count "
                        << "conversion."
                        << abort(FatalError);
                }

                cloud_.addNewParcel
                (
                    p,
                    U,
                    RWF,
                    ERot,
                    ELevel,
                    cellI,
                    faces_[f],
                    faceTetIs.tetPt(),
                    typeId,
                    newParcel,
                    0,
                    vibLevel
                );
            }
        }
    }
}

void dsmcKuangFreeStreamInflowPatch::controlParcelsBeforeCollisions()
{

}

void dsmcKuangFreeStreamInflowPatch::controlParcelsAfterCollisions()
{
}

void dsmcKuangFreeStreamInflowPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcKuangFreeStreamInflowPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    // 刷新本派生模型的属性字典，避免继续使用构造时的旧 subDict。
    propsDict_ = newDict.subDict(typeName + "Properties");

    // 参考 dsmcWangPressureInlet：运行期字典更新后必须重新解析本模型属性。
    // 这样 velocityFrame、velocity、omega、温度、数密度和 typeIds 的修改
    // 才会真正进入后续注入通量计算。
    setProperties();
}



void dsmcKuangFreeStreamInflowPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");
    translationalTemperature_ = readScalar(propsDict_.lookup("translationalTemperature"));
    rotationalTemperature_ = propsDict_.lookupOrDefault<scalar>("rotationalTemperature", 0.0);
    vibrationalTemperature_ = propsDict_.lookupOrDefault<scalar>("vibrationalTemperature", 0.0);
    electronicTemperature_ = propsDict_.lookupOrDefault<scalar>("electronicTemperature", 0.0);

    // 参考 dsmcWangPressureInlet 的 velocityFrame 设计：
    // local    保持旧 dsmcKuang 行为，velocity/omega 直接表示 parcel 局部速度；
    // absolute 表示 velocity/omega 给出实验室系 reservoir 速度，注入旋转 cell
    //          前用 dsmcGeneralBoundary::velocityToLocalFrame 转成局部速度。
    velocityFrame_ =
        propsDict_.lookupOrDefault<word>("velocityFrame", word("local"));

    // 拼写错误会直接改变旋转参考系中的边界通量，因此必须立即报错。
    if (velocityFrame_ != "local" && velocityFrame_ != "absolute")
    {
        FatalErrorIn("dsmcKuangFreeStreamInflowPatch::setProperties()")
            << "velocityFrame must be either 'local' or 'absolute' for patch "
            << patchName_ << ", but got " << velocityFrame_ << nl
            << "in: " << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    omega_ = propsDict_.lookupOrDefault<scalar>("omega", 0.0);

    // 兼容旧 dsmcKuang 字典：如果没有显式给旋转轴和中心，则默认绕 z 轴、
    // 原点旋转。omega=0 时这两个量不会影响结果。
    rotationAxis_ =
        propsDict_.lookupOrDefault<vector>("rotationAxis", vector(0, 0, 1));

    centerPoint_ =
        propsDict_.lookupOrDefault<vector>("centerPoint", vector::zero);

    // 只要 omega 非零，旋转轴就必须合法；提前检查比运行中插入粒子时报错更清楚。
    if (mag(omega_) > SMALL && mag(rotationAxis_) <= SMALL)
    {
        FatalErrorIn("dsmcKuangFreeStreamInflowPatch::setProperties()")
            << "A non-zero omega requires a non-zero rotationAxis for patch "
            << patchName_ << nl
            << "omega = " << omega_ << nl
            << "rotationAxis = " << rotationAxis_ << nl
            << "in: " << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    //  read in the type ids

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcKuangFreeStreamInflowPatch::dsmcKuangFreeStreamInflowPatch()")
            << "Cannot have zero typeIds being inserted." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word moleculeName(molecules[i]);

        if(findIndex(moleculesReduced, moleculeName) == -1)
        {
            moleculesReduced.append(moleculeName);
        }
    }

    moleculesReduced.shrink();

    //  set the type ids

    typeIds_.setSize(moleculesReduced.size(), -1);

    forAll(moleculesReduced, i)
    {
        const word moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcKuangFreeStreamInflowPatch::dsmcKuangFreeStreamInflowPatch()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // read in the mass density per specie

    const dictionary& numberDensitiesDict
    (
        propsDict_.subDict("numberDensities")
    );

    numberDensities_.clear();

    numberDensities_.setSize(typeIds_.size(), 0.0);

    forAll(numberDensities_, i)
    {
        numberDensities_[i] = readScalar
        (
            numberDensitiesDict.lookup(moleculesReduced[i])
        );
    }

    // set the accumulator

    accumulatedParcelsToInsert_.setSize(typeIds_.size());

    forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}


void dsmcKuangFreeStreamInflowPatch::setNewBoundaryFields()
{
    patchId_ = mesh_.boundaryMesh().findPatchID(patchName_);

    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.setSize(patch.size());
    cells_.setSize(patch.size());

    //- loop through all faces and set the boundary cells
    //- no conflict with parallelisation because the faces are unique

    nFaces_ = 0;
    patchSurfaceArea_ = 0.0;

    for(label i = 0; i < patch.size(); i++)
    {
        label globalFaceI = patch.start() + i;

        faces_[i] = globalFaceI;
        cells_[i] = patch.faceCells()[i];
        nFaces_++;
        patchSurfaceArea_ += mag(mesh_.faceAreas()[globalFaceI]);
    }

    if(Pstream::parRun())
    {
        reduce(patchSurfaceArea_, sumOp<scalar>());
    }

   forAll(accumulatedParcelsToInsert_, m)
    {
        accumulatedParcelsToInsert_[m].setSize(nFaces_, 0.0);
    }
}




} // End namespace Foam

// ************************************************************************* //
