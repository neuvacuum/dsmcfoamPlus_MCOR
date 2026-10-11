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

Measures DSMC macroscopic fields for a single species or a gas mixture and
writes the results to a volume field that can be viewed in Paraview.

Translational, rotatational and vibrational temperature fields will also be
written automatically.

Boundary fields are measured in conjunction with the boundaryMeasurements class
and are also written.

\*---------------------------------------------------------------------------*/

#include "dsmcVolFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(dsmcVolFields, 0);

addToRunTimeSelectionTable(dsmcField, dsmcVolFields, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace
{

bool validVolFieldOutputProperty(const Foam::word& propertyName)
{
    return
       propertyName == "density"
    || propertyName == "none"
    || propertyName == "numberDensity"
    || propertyName == "massDensity"
    || propertyName == "parcelCount"
    || propertyName == "pressure"
    || propertyName == "temperature"
    || propertyName == "velocity"
    || propertyName == "absoluteVelocity"
    || propertyName == "thetaVelocity"
    || propertyName == "srfDiagnostics"
    || propertyName == "mach"
    || propertyName == "wallDiagnostics"
    || propertyName == "meanFreePath"
    || propertyName == "collisionDiagnostics"
    || propertyName == "classification"
    || propertyName == "errors"
    || propertyName == "heatFluxShearStress"
    || propertyName == "dsmcN"
    || propertyName == "dsmcNMean"
    || propertyName == "rhoN"
    || propertyName == "rhoM"
    || propertyName == "p"
    || propertyName == "Ttra"
    || propertyName == "Trot"
    || propertyName == "Tvib"
    || propertyName == "Telec"
    || propertyName == "Tov"
    || propertyName == "U"
    || propertyName == "UMag"
    || propertyName == "UAbs"
    || propertyName == "UAbsMag"
    || propertyName == "UThetaRel"
    || propertyName == "UThetaAbs"
    || propertyName == "rotatingRegion"
    || propertyName == "aCoriolis"
    || propertyName == "aCentrifugal"
    || propertyName == "aSRF"
    || propertyName == "Ma"
    || propertyName == "q"
    || propertyName == "fD"
    || propertyName == "tau"
    || propertyName == "wallHeatFlux"
    || propertyName == "wallShearStress"
    || propertyName == "mfp"
    || propertyName == "mfpToDx"
    || propertyName == "DxToMfp"
    || propertyName == "meanCollisionRate"
    || propertyName == "measuredCollisionRate"
    || propertyName == "meanCollisionSeparation"
    || propertyName == "mct"
    || propertyName == "mctToDt"
    || propertyName == "SOFP"
    || propertyName == "classIDistribution"
    || propertyName == "classIIDistribution"
    || propertyName == "classIIIDistribution"
    || propertyName == "densityError"
    || propertyName == "velocityError"
    || propertyName == "temperatureError"
    || propertyName == "pressureError"
    || propertyName == "heatFluxVector"
    || propertyName == "pressureTensor"
    || propertyName == "shearStressTensor";
}

Foam::vector absoluteVelocityFromMean
(
    const Foam::dsmcCloud& cloud,
    const Foam::label celli,
    const Foam::point& samplePosition,
    const Foam::vector& meanVelocity
)
{
    if (cloud.isInRotatingRegion(celli))
    {
        return
            meanVelocity
          + (
                cloud.rotatingFrameOmega()
              ^ (samplePosition - cloud.rotatingFrameOrigin())
            );
    }

    return meanVelocity;
}

Foam::vector coriolisAccelerationFromMean
(
    const Foam::dsmcCloud& cloud,
    const Foam::vector& relativeVelocity
)
{
    // 这里直接复用 SRF 物理公式中的科里奥利项：
    //     a_cor = -2 * (Omega x U_rel)
    // 之所以单独拆成辅助函数，是为了让后面：
    // 1. 内部单元诊断；
    // 2. boundary 采样诊断；
    // 3. 总 SRF 加速度重构；
    // 共用完全同一套表达式，避免多处手写后逐渐漂移。
    return -2.0*(cloud.rotatingFrameOmega() ^ relativeVelocity);
}

Foam::vector centrifugalAccelerationFromPosition
(
    const Foam::dsmcCloud& cloud,
    const Foam::point& samplePosition
)
{
    // 离心项只依赖旋转角速度与空间位置：
    //     a_cen = -(Omega x (Omega x r))
    // 其中
    //     r = x - origin
    const Foam::vector relPosition =
        samplePosition - cloud.rotatingFrameOrigin();
    const Foam::vector& omega = cloud.rotatingFrameOmega();

    return -(omega ^ (omega ^ relPosition));
}

Foam::scalar tangentialVelocityComponent
(
    const Foam::dsmcCloud& cloud,
    const Foam::point& samplePosition,
    const Foam::vector& velocity
)
{
    const scalar omegaMag = mag(cloud.rotatingFrameOmega());

    // 这里把 rotatingFrameOmega 仅作为“定义柱坐标轴向方向”的全局几何轴。
    // 因而即便某个单元本身不在 SRF 作用域内，它的绝对速度 U_abs
    // 仍然可以投影到同一套 e_theta 上，用于比较 rotor/stator 两侧
    // 的周向速度方向是否连续。
    //
    // 若角速度本身为 0，则柱坐标基底退化，周向速度定义也失去意义。
    if (omegaMag < VSMALL)
    {
        return 0.0;
    }

    // 柱坐标轴向单位向量 e_axis = Omega/|Omega|。
    const vector axisHat = cloud.rotatingFrameOmega()/omegaMag;
    const vector relPosition = samplePosition - cloud.rotatingFrameOrigin();

    // 径向向量 r_perp = r - (r·e_axis)e_axis，对应旋转轴垂直平面内的投影。
    const vector radialVector =
        relPosition - (axisHat & relPosition)*axisHat;
    const scalar radialMag = mag(radialVector);

    // 若采样点恰好落在旋转轴上，则 e_r 与 e_theta 都退化；
    // 这种情况下把周向速度定义为 0，避免数值上除以极小半径。
    if (radialMag < VSMALL)
    {
        return 0.0;
    }

    const vector radialHat = radialVector/radialMag;

    // 右手系周向单位向量 e_theta = e_axis x e_r。
    // 因而 u_theta = U · e_theta。
    // 这里故意不把速度限定为 U_abs 或 U_rel：
    // - 传入 UAbsMean 时，得到绝对参考系周向速度；
    // - 传入 UMean  时，得到相对参考系周向速度。
    // 这样可用完全同一套柱坐标基底比较 rotor/stator 两侧
    // “绝对速度看起来像什么”与“rotor 相对系里真实推进速度像什么”。
    const vector thetaHat = axisHat ^ radialHat;

    return velocity & thetaHat;
}

Foam::scalar absoluteTangentialVelocityComponent
(
    const Foam::dsmcCloud& cloud,
    const Foam::point& samplePosition,
    const Foam::vector& absoluteVelocity
)
{
    // 绝对速度周向分量是上面通用投影公式的一个特例。
    // 单独保留这个包装函数，可以避免下面大量调用点重命名，
    // 也能让“绝对速度”和“相对速度”两种后处理语义保持清晰可读。
    return tangentialVelocityComponent
    (
        cloud,
        samplePosition,
        absoluteVelocity
    );
}

Foam::scalar relativeTangentialVelocityComponent
(
    const Foam::dsmcCloud& cloud,
    const Foam::point& samplePosition,
    const Foam::vector& relativeVelocity
)
{
    // 相对速度周向分量与绝对速度使用相同的柱坐标基底，
    // 仅把被投影的速度改为 rotor 区域真实推进所用的 U_rel。
    return tangentialVelocityComponent
    (
        cloud,
        samplePosition,
        relativeVelocity
    );
}

}

void dsmcVolFields::calculateWallUnitVectors()
{
    forAll(n_, patchi)
    {
        const polyPatch& pPatch = mesh_.boundaryMesh()[patchi];

        if (isA<wallPolyPatch>(pPatch))
        {
            // faceCentres() 由临时场提供数据，必须在整个 patch 循环中持有副本，
            // 否则多个 dsmcVolFields 实例会反复读取已失效的临时存储。
            const vectorField fC(pPatch.faceCentres());

            forAll(n_[patchi], facei)
            {
                n_[patchi][facei] = pPatch.faceAreas()[facei]
                    /mag(pPatch.faceAreas()[facei]);

                //- Wall tangential unit vector. Use the direction between the
                // face centre and the first vertex in the list
                t1_[patchi][facei] = fC[facei]
                    - mesh_.points()[mesh_.faces()[pPatch.start() + facei][0]];
                t1_[patchi][facei] /= mag(t1_[patchi][facei]);

                //- Other tangential unit vector.  Rescaling in case face is not
                //  flat and n and t1 aren't perfectly orthogonal
                t2_[patchi][facei] = n_[patchi][facei]^t1_[patchi][facei];
                t2_[patchi][facei] /= mag(t2_[patchi][facei]);
            }
        }
    }
}

void dsmcVolFields::configureOutputProperties()
{
    outputProperties_.clear();

    if (!propsDict_.found("outputProperties"))
    {
        return;
    }

    const List<word> propertyNames(propsDict_.lookup("outputProperties"));

    if (propertyNames.size() == 0)
    {
        FatalErrorIn("dsmcVolFields::configureOutputProperties()")
            << "outputProperties is present but empty for field "
            << fieldName_ << nl << "in: "
            << time_.time().system()/"fieldPropertiesDict"
            << exit(FatalError);
    }

    outputProperties_.resize(2*propertyNames.size());

    forAll(propertyNames, i)
    {
        const word& propertyName = propertyNames[i];

        if (!validVolFieldOutputProperty(propertyName))
        {
            FatalErrorIn("dsmcVolFields::configureOutputProperties()")
                << "Cannot find dsmcVolFields output property: "
                << propertyName << nl
                << "Valid high-level groups include: "
                << "(none pressure velocity absoluteVelocity thetaVelocity "
                << "numberDensity massDensity parcelCount temperature "
                << "meanFreePath collisionDiagnostics mach srfDiagnostics "
                << "wallDiagnostics)"
                << nl << "in: "
                << time_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        if (propertyName == "density")
        {
            outputProperties_.insert("numberDensity");
            outputProperties_.insert("massDensity");
        }
        else if (propertyName == "q")
        {
            outputProperties_.insert("wallHeatFlux");
        }
        else if (propertyName == "tau")
        {
            outputProperties_.insert("wallShearStress");
        }
        else
        {
            outputProperties_.insert(propertyName);
        }
    }

    // 只在显式配置 outputProperties 时关闭 AUTO_WRITE。
    // 这样旧算例保持原来的全量写出行为；新算例只由下面的
    // writeSelected() 白名单负责落盘，避免 runTime.write() 又把全部
    // 宏观量自动写出。
    setAllFieldWriteOptions(IOobject::NO_WRITE);

    Info<< "dsmcVolFields " << fieldName_
        << " selective outputProperties " << propertyNames << endl;
}

void dsmcVolFields::setAllFieldWriteOptions
(
    const IOobject::writeOption writeOpt
)
{
    dsmcN_.writeOpt() = writeOpt;
    dsmcNMean_.writeOpt() = writeOpt;
    rhoN_.writeOpt() = writeOpt;
    rhoM_.writeOpt() = writeOpt;
    p_.writeOpt() = writeOpt;
    Ttra_.writeOpt() = writeOpt;
    Trot_.writeOpt() = writeOpt;
    Tvib_.writeOpt() = writeOpt;
    Telec_.writeOpt() = writeOpt;
    Tov_.writeOpt() = writeOpt;
    q_.writeOpt() = writeOpt;
    tau_.writeOpt() = writeOpt;
    mfp_.writeOpt() = writeOpt;
    mfpToDx_.writeOpt() = writeOpt;
    DxToMfp_.writeOpt() = writeOpt;
    meanCollisionRate_.writeOpt() = writeOpt;
    measuredCollisionRate_.writeOpt() = writeOpt;
    meanCollisionTime_.writeOpt() = writeOpt;
    mctToDt_.writeOpt() = writeOpt;
    meanCollisionSeparation_.writeOpt() = writeOpt;
    SOF_.writeOpt() = writeOpt;
    Ma_.writeOpt() = writeOpt;
    classIDistribution_.writeOpt() = writeOpt;
    classIIDistribution_.writeOpt() = writeOpt;
    classIIIDistribution_.writeOpt() = writeOpt;
    densityError_.writeOpt() = writeOpt;
    velocityError_.writeOpt() = writeOpt;
    temperatureError_.writeOpt() = writeOpt;
    pressureError_.writeOpt() = writeOpt;
    UMean_.writeOpt() = writeOpt;
    UMagMean_.writeOpt() = writeOpt;
    UAbsMean_.writeOpt() = writeOpt;
    UAbsMagMean_.writeOpt() = writeOpt;
    UThetaRelMean_.writeOpt() = writeOpt;
    UThetaAbsMean_.writeOpt() = writeOpt;
    rotatingRegionMean_.writeOpt() = writeOpt;
    aCoriolisMean_.writeOpt() = writeOpt;
    aCentrifugalMean_.writeOpt() = writeOpt;
    aSRFMean_.writeOpt() = writeOpt;
    fD_.writeOpt() = writeOpt;
    heatFluxVector_.writeOpt() = writeOpt;
    pressureTensor_.writeOpt() = writeOpt;
    shearStressTensor_.writeOpt() = writeOpt;
}

bool dsmcVolFields::outputProperty(const word& propertyName) const
{
    return outputProperties_.size() == 0 || outputProperties_.found(propertyName);
}

bool dsmcVolFields::outputProperty
(
    const word& propertyName,
    const word& propertyGroup
) const
{
    return
       outputProperty(propertyName)
    || (
           outputProperties_.size() != 0
        && propertyGroup != word::null
        && outputProperties_.found(propertyGroup)
       );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcVolFields::dsmcVolFields
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
    mfpTref_(273.0),
    fieldName_(propsDict_.lookup("fieldName")),
    speciesIds_(),
    dsmcN_
    (
        IOobject
        (
            "dsmcN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    dsmcNMean_
    (
        IOobject
        (
            "dsmcNMean_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless/dimVolume, 0.0)
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimMass/dimVolume, 0.0)
    ),
    p_
    (
        IOobject
        (
            "p_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimPressure, 0.0)
    ),
    Ttra_
    (
        IOobject
        (
            "Ttra_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Trot_
    (
        IOobject
        (
            "Trot_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Tvib_
    (
        IOobject
        (
            "Tvib_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Telec_
    (
        IOobject
        (
            "Telec_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    Tov_
    (
        IOobject
        (
            "Tov_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimTemperature, 0.0)
    ),
    q_
    (
        IOobject
        (
            "wallHeatFlux_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(1, 0, -3, 0, 0), 0.0)
    ),
    tau_
    (
        IOobject
        (
            "wallShearStress_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimPressure, 0.0)
    ),
    mfp_
    (
        IOobject
        (
            "mfp_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    mfpToDx_
    (
        IOobject
        (
            "mfpToDx_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    DxToMfp_
    (
        IOobject
        (
            "DxToMfp_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    meanCollisionRate_
    (
        IOobject
        (
            "meanCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    measuredCollisionRate_
    (
        IOobject
        (
            "measuredCollisionRate_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, -1, 0, 0), 0.0)
    ),
    meanCollisionTime_
    (
        IOobject
        (
            "mct_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimensionSet(0, 0, 1, 0, 0), 0.0)
    ),
    mctToDt_
    (
        IOobject
        (
            "mctToDt_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero",  dimless, 0.0)
    ),
    meanCollisionSeparation_
    (
        IOobject
        (
            "mcs_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    SOF_
    (
        IOobject
        (
            "SOFP_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    Ma_
    (
        IOobject
        (
            "Ma_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIDistribution_
    (
        IOobject
        (
            "classIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIDistribution_
    (
        IOobject
        (
            "classIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    classIIIDistribution_
    (
        IOobject
        (
            "classIIIDistribution_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    densityError_
    (
        IOobject
        (
            "rhoMError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    velocityError_
    (
        IOobject
        (
            "UError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    temperatureError_
    (
        IOobject
        (
            "TError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pressureError_
    (
        IOobject
        (
            "pError_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    UMean_
    (
        IOobject
        (
            "U_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    UMagMean_
    (
        IOobject
        (
            "UMag_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimLength/dimTime, 0.0)
    ),
    UAbsMean_
    (
        IOobject
        (
            "UAbs_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/dimTime, vector::zero)
    ),
    UAbsMagMean_
    (
        IOobject
        (
            "UAbsMag_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimLength/dimTime, 0.0)
    ),
    UThetaRelMean_
    (
        IOobject
        (
            "UThetaRel_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimLength/dimTime, 0.0)
    ),
    UThetaAbsMean_
    (
        IOobject
        (
            "UThetaAbs_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimLength/dimTime, 0.0)
    ),
    rotatingRegionMean_
    (
        IOobject
        (
            "rotatingRegion_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0.0", dimless, 0.0)
    ),
    aCoriolisMean_
    (
        IOobject
        (
            "aCoriolis_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/sqr(dimTime), vector::zero)
    ),
    aCentrifugalMean_
    (
        IOobject
        (
            "aCentrifugal_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/sqr(dimTime), vector::zero)
    ),
    aSRFMean_
    (
        IOobject
        (
            "aSRF_"+ fieldName_,
            time_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0.0", dimLength/sqr(dimTime), vector::zero)
    ),
    fD_
    (
        IOobject
        (
            "fD_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, -1, -2, 0, 0),
            vector::zero
        )
    ),
    heatFluxVector_
    (
        IOobject
        (
            "heatFluxVector_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            dimensionSet(1, 0, -3, 0, 0),
            vector::zero
        )
    ),
    pressureTensor_
    (
        IOobject
        (
            "pressureTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    shearStressTensor_
    (
        IOobject
        (
            "shearStressTensor_"+ fieldName_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor
        (
            "zero",
            dimPressure,
            tensor::zero
        )
    ),
    dsmcNCum_(mesh_.nCells(), 0.0),
    nCum_(mesh_.nCells(), 0.0),
    dsmcNElecLvlCum_(mesh_.nCells(), 0.0),
    dsmcMCum_(mesh_.nCells(), 0.0),
    mCum_(mesh_.nCells(), 0.0),
    dsmcLinearKECum_(mesh_.nCells(), 0.0),
    linearKECum_(mesh_.nCells(), 0.0),
    dsmcErotCum_(mesh_.nCells(), 0.0),
    dsmcZetaRotCum_(mesh_.nCells(), 0.0),
    dsmcMuuCum_(mesh_.nCells(), 0.0),
    dsmcMuvCum_(mesh_.nCells(), 0.0),
    dsmcMuwCum_(mesh_.nCells(), 0.0),
    dsmcMvvCum_(mesh_.nCells(), 0.0),
    dsmcMvwCum_(mesh_.nCells(), 0.0),
    dsmcMwwCum_(mesh_.nCells(), 0.0),
    dsmcMccCum_(mesh_.nCells(), 0.0),
    dsmcMccuCum_(mesh_.nCells(), 0.0),
    dsmcMccvCum_(mesh_.nCells(), 0.0),
    dsmcMccwCum_(mesh_.nCells(), 0.0),
    dsmcEuCum_(mesh_.nCells(), 0.0),
    dsmcEvCum_(mesh_.nCells(), 0.0),
    dsmcEwCum_(mesh_.nCells(), 0.0),
    dsmcECum_(mesh_.nCells(), 0.0),
    zetaVib_(mesh_.nCells(), 0.0),
    dsmcNClassICum_(mesh_.nCells(), 0.0),
    dsmcNClassIICum_(mesh_.nCells(), 0.0),
    dsmcNClassIIICum_(mesh_.nCells(), 0.0),
    collisionSeparation_(mesh_.nCells(), 0.0),
    dsmcNCollsCum_(mesh_.nCells(), 0.0),
    dsmcMomentumCum_(mesh.nCells(), vector::zero),
    momentumCum_(mesh.nCells(), vector::zero),
    boundaryCells_(),
    dsmcSpeciesEvibModCum_(),
    dsmcSpeciesEelecCum_(),
    dsmcNSpeciesCum_(),
    nSpeciesCum_(),
    dsmcMccSpeciesCum_(),
    speciesTvib_(),
    dsmcNGrndElecLvlSpeciesCum_(),
    dsmcN1stElecLvlSpeciesCum_(),
    speciesMfp_(),
    speciesMcr_(),
    rhoNBF_(),
    rhoMBF_(),
    linearKEBF_(),
    ErotBF_(),
    zetaRotBF_(),
    qBF_(),
    zetaVibBF_(),
    rhoNIntBF_(),
    rhoNElecBF_(),
    momentumBF_(),
    fDBF_(),
    speciesEvibBF_(),
    speciesEelecBF_(),
    speciesRhoNBF_(),
    speciesMccBF_(),
    speciesTvibBF_(),
    speciesZetaVibBF_(),
    speciesEvibModBF_(),
    n_(),
    t1_(),
    t2_(),
    averagingAcrossManyRuns_(false),
    measureClassifications_(false),
    measureMeanFreePath_(false),
    measureErrors_(false),
    densityOnly_(false),
    measureHeatFluxShearStress_(false),
    writeRotationalTemperature_(false),
    writeVibrationalTemperature_(false),
    writeElectronicTemperature_(false),
    outputProperties_()
{
    configureOutputProperties();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcVolFields::~dsmcVolFields()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcVolFields::readIn()
{
    IOdictionary resumeSampling
    (
        IOobject
        (
            "resumeSampling_" + fieldName_,
            time_.time().timeName(),
            "uniform",
            time_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (resumeSampling.size() > 0)
    {
        dictionary dict
        (
            resumeSampling.readStream(resumeSampling.filePath())
        );

        dict.readIfPresent("nTimeSteps", nTimeSteps_);

        if
        (
            nTimeSteps_ > 0.0
         && dict.lookupOrDefault<label>("physicalMomentWeightingVersion", 0) != 1
        )
        {
            FatalErrorInFunction
                << "Stored sampling moments use an incompatible weighting "
                << "convention. Start a fresh sampling window; do not mix "
                << "unweighted and represented-molecule accumulators."
                << exit(FatalError);
        }

        // DSMC parcel related cumulative values
        dict.readIfPresent("dsmcNCum", dsmcNCum_);
        dict.readIfPresent("dsmcMCum", dsmcMCum_);

        dict.readIfPresent("dsmcLinearKECum", dsmcLinearKECum_);
        dict.readIfPresent("dsmcMomentumCum", dsmcMomentumCum_);
        dict.readIfPresent("dsmcErotCum", dsmcErotCum_);
        dict.readIfPresent("dsmcZetaRotCum", dsmcZetaRotCum_);
        dict.readIfPresent("dsmcSpeciesEelecCum", dsmcSpeciesEelecCum_);
        dict.readIfPresent("dsmcNSpeciesCum", dsmcNSpeciesCum_);
        dict.readIfPresent("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
        dict.readIfPresent("dsmcMuuCum", dsmcMuuCum_);
        dict.readIfPresent("dsmcMuvCum", dsmcMuvCum_);
        dict.readIfPresent("dsmcMuwCum", dsmcMuwCum_);
        dict.readIfPresent("dsmcMvvCum", dsmcMvvCum_);
        dict.readIfPresent("dsmcMvwCum", dsmcMvwCum_);
        dict.readIfPresent("dsmcMwwCum", dsmcMwwCum_);
        dict.readIfPresent("dsmcMccCum", dsmcMccCum_);
        dict.readIfPresent("dsmcMccuCum", dsmcMccuCum_);
        dict.readIfPresent("dsmcMccvCum", dsmcMccvCum_);
        dict.readIfPresent("dsmcMccwCum", dsmcMccwCum_);
        dict.readIfPresent("dsmcEuCum", dsmcEuCum_);
        dict.readIfPresent("dsmcEvCum", dsmcEvCum_);
        dict.readIfPresent("dsmcEwCum", dsmcEwCum_);
        dict.readIfPresent("dsmcECum", dsmcECum_);
        dict.readIfPresent("dsmcNElecLvlCum", dsmcNElecLvlCum_);
        dict.readIfPresent
        (
            "dsmcNGrndElecLvlSpeciesCum",
            dsmcNGrndElecLvlSpeciesCum_
        );
        dict.readIfPresent
        (
            "dsmcN1stElecLvlSpeciesCum",
            dsmcN1stElecLvlSpeciesCum_
        );
        if (measureClassifications_)
        {
          dict.readIfPresent("dsmcNClassICum", dsmcNClassICum_);
          dict.readIfPresent("dsmcNClassIICum", dsmcNClassIICum_);
          dict.readIfPresent("dsmcNClassIIICum", dsmcNClassIIICum_);
        }
        dict.readIfPresent("dsmcSpeciesEvibModCum", dsmcSpeciesEvibModCum_);
        dict.readIfPresent("dsmcNCollsCum", dsmcNCollsCum_);

        // boundary measurements
        dict.readIfPresent("rhoNBF", rhoNBF_);
        dict.readIfPresent("rhoMBF", rhoMBF_);
        dict.readIfPresent("linearKEBF", linearKEBF_);
        dict.readIfPresent("momentumBF", momentumBF_);
        dict.readIfPresent("ErotBF", ErotBF_);
        dict.readIfPresent("zetaRotBF", zetaRotBF_);
        dict.readIfPresent("rhoNIntBF", rhoNIntBF_);
        dict.readIfPresent("rhoNElecBF", rhoNElecBF_);
        dict.readIfPresent("qBF", qBF_);
        dict.readIfPresent("fDBF", fDBF_);
        dict.readIfPresent("speciesRhoNBF", speciesRhoNBF_);
        dict.readIfPresent("speciesEvibBF", speciesEvibBF_);
        dict.readIfPresent("speciesEelecBF",  speciesEelecBF_);
        dict.readIfPresent("speciesMccBF", speciesMccBF_);
        dict.readIfPresent("speciesEvibModBF", speciesEvibModBF_);

        // cumulative values
        dict.readIfPresent("nCum", nCum_);
        dict.readIfPresent("mCum", mCum_);
        dict.readIfPresent("nSpeciesCum", nSpeciesCum_);
        dict.readIfPresent("momentumCum", momentumCum_);
        dict.readIfPresent("linearKECum", linearKECum_);
        dict.readIfPresent("collisionSeparation", collisionSeparation_);
    }
}

void dsmcVolFields::writeOut()
{
    if (time_.time().outputTime())
    {
        IOdictionary dict
        (
            IOobject
            (
                "resumeSampling_" + fieldName_,
                time_.time().timeName(),
                "uniform",
                time_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        dict.add("nTimeSteps", nTimeSteps_);
        dict.add("physicalMomentWeightingVersion", 1);

        // DSMC parcel related cumulative values
        dict.add("dsmcNCum", dsmcNCum_);
        dict.add("dsmcMCum", dsmcMCum_);

        dict.add("dsmcLinearKECum", dsmcLinearKECum_);
        dict.add("dsmcMomentumCum", dsmcMomentumCum_);
        dict.add("dsmcErotCum", dsmcErotCum_);
        dict.add("dsmcZetaRotCum", dsmcZetaRotCum_);
        dict.add("dsmcSpeciesEelecCum", dsmcSpeciesEelecCum_);
        dict.add("dsmcNSpeciesCum", dsmcNSpeciesCum_);
        dict.add("dsmcMccSpeciesCum", dsmcMccSpeciesCum_);
        dict.add("dsmcMuuCum", dsmcMuuCum_);
        dict.add("dsmcMuvCum", dsmcMuvCum_);
        dict.add("dsmcMuwCum", dsmcMuwCum_);
        dict.add("dsmcMvvCum", dsmcMvvCum_);
        dict.add("dsmcMvwCum", dsmcMvwCum_);
        dict.add("dsmcMwwCum", dsmcMwwCum_);
        dict.add("dsmcMccCum", dsmcMccCum_);
        dict.add("dsmcMccuCum", dsmcMccuCum_);
        dict.add("dsmcMccvCum", dsmcMccvCum_);
        dict.add("dsmcMccwCum", dsmcMccwCum_);
        dict.add("dsmcEuCum", dsmcEuCum_);
        dict.add("dsmcEvCum", dsmcEvCum_);
        dict.add("dsmcEwCum", dsmcEwCum_);
        dict.add("dsmcECum", dsmcECum_);
        dict.add("dsmcNElecLvlCum", dsmcNElecLvlCum_);
        dict.add("dsmcNGrndElecLvlSpeciesCum", dsmcNGrndElecLvlSpeciesCum_);
        dict.add("dsmcN1stElecLvlSpeciesCum", dsmcN1stElecLvlSpeciesCum_);
        if (measureClassifications_)
        {
            dict.add("dsmcNClassICum", dsmcNClassICum_);
            dict.add("dsmcNClassIICum", dsmcNClassIICum_);
            dict.add("dsmcNClassIIICum", dsmcNClassIIICum_);
        }
        dict.add("dsmcSpeciesEvibModCum", dsmcSpeciesEvibModCum_);
        dict.add("dsmcNCollsCum", dsmcNCollsCum_);

        // cumulative values
        dict.add("nCum", nCum_);
        dict.add("mCum", mCum_);
        dict.add("nSpeciesCum", nSpeciesCum_);
        dict.add("momentumCum", momentumCum_);
        dict.add("linearKECum", linearKECum_);
        dict.add("collisionSeparation", collisionSeparation_);

        // boundary measurements
        dict.add("rhoNBF", rhoNBF_);
        dict.add("rhoMBF", rhoMBF_);
        dict.add("linearKEBF", linearKEBF_);
        dict.add("momentumBF", momentumBF_);
        dict.add("ErotBF", ErotBF_);
        dict.add("zetaRotBF", zetaRotBF_);
        dict.add("rhoNIntBF", rhoNIntBF_);
        dict.add("rhoNElecBF", rhoNElecBF_);
        dict.add("qBF", qBF_);
        dict.add("fDBF", fDBF_);
        dict.add("speciesRhoNBF", speciesRhoNBF_);
        dict.add("speciesEvibBF", speciesEvibBF_);
        dict.add("speciesEelecBF",  speciesEelecBF_);
        dict.add("speciesMccBF", speciesMccBF_);
        dict.add("speciesEvibModBF", speciesEvibModBF_);

        IOstream::streamFormat fmt = time_.time().writeFormat();
        IOstream::versionNumber ver = time_.time().writeVersion();
        IOstream::compressionType cmp = time_.time().writeCompression();

        dict.regIOobject::writeObject(fmt, ver, cmp);
    }
}

//- Initial configuration
void dsmcVolFields::createField()
{
    Info << "Initialising dsmcVolFields field" << endl;

    const List<word>& species (propsDict_.lookup("typeIds"));

    DynamicList<word> speciesReduced(0);

    forAll(species, i)
    {
        const word& speciesName(species[i]);

        if (findIndex(speciesReduced, speciesName) == -1)
        {
            speciesReduced.append(speciesName);
        }
    }

    speciesReduced.shrink();

    speciesIds_.setSize(speciesReduced.size(), -1);

    forAll(speciesReduced, i)
    {
        const word& speciesName = speciesReduced[i];

        const label spId = findIndex(cloud_.typeIdList(), speciesName);

        if (spId == -1)
        {
            FatalErrorIn("dsmcVolFields::dsmcVolFields()")
                << "Cannot find typeId: " << speciesName << nl << "in: "
                << mesh_.time().system()/"fieldPropertiesDict"
                << exit(FatalError);
        }

        speciesIds_[i] = spId;
    }

    const label nCells = mesh_.nCells();
    const label nSpecies = speciesIds_.size();
    const label nPatches = mesh_.boundaryMesh().size();

    //- Volume fields initialisation
    dsmcNSpeciesCum_.setSize(nSpecies);
    nSpeciesCum_.setSize(nSpecies);
    dsmcMccSpeciesCum_.setSize(nSpecies);
    speciesMfp_.setSize(nSpecies);
    speciesMcr_.setSize(nSpecies);
    speciesTvib_.setSize(nSpecies);
    dsmcSpeciesEvibModCum_.setSize(nSpecies);
    dsmcSpeciesEelecCum_.setSize(nSpecies);
    dsmcNGrndElecLvlSpeciesCum_.setSize(nSpecies);
    dsmcN1stElecLvlSpeciesCum_.setSize(nSpecies);

    forAll(speciesIds_, i)
    {
        const label spId = speciesIds_[i];
        const label zetaRot =
            cloud_.constProps(spId).rotationalDegreesOfFreedom();
        const label nVibMod = cloud_.constProps(spId).nVibrationalModes();
        const label nElecLevels = cloud_.constProps(spId).nElectronicLevels();

        if (zetaRot > 0 and (not writeRotationalTemperature_))
        {
            writeRotationalTemperature_ = true;
        }

        if (nVibMod > 0 and (not writeVibrationalTemperature_))
        {
            writeVibrationalTemperature_ = true;
        }

        if (nElecLevels > 1 and (not writeElectronicTemperature_))
        {
            writeElectronicTemperature_ = true;
        }

        dsmcNSpeciesCum_[i].setSize(nCells, 0.0);
        nSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcMccSpeciesCum_[i].setSize(nCells, 0.0);
        speciesMfp_[i].setSize(nCells, 0.0);
        speciesMcr_[i].setSize(nCells, 0.0);
        speciesTvib_[i].setSize(nCells);
        dsmcSpeciesEvibModCum_[i].setSize(nVibMod);

        forAll(dsmcSpeciesEvibModCum_[i], j)
        {
            dsmcSpeciesEvibModCum_[i][j].setSize(nCells, 0.0);
        }

        dsmcSpeciesEelecCum_[i].setSize(nCells, 0.0);
        dsmcNGrndElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcN1stElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
    }

    //- Boundary fields initialisation
    boundaryCells_.setSize(nPatches);
    rhoNBF_.setSize(nPatches);
    rhoMBF_.setSize(nPatches);
    linearKEBF_.setSize(nPatches);
    momentumBF_.setSize(nPatches);
    ErotBF_.setSize(nPatches);
    zetaRotBF_.setSize(nPatches);
    qBF_.setSize(nPatches);
    fDBF_.setSize(nPatches);
    zetaVibBF_.setSize(nPatches);
    rhoNIntBF_.setSize(nPatches);
    rhoNElecBF_.setSize(nPatches);

    n_.setSize(nPatches);
    t1_.setSize(nPatches);
    t2_.setSize(nPatches);

    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        boundaryCells_[j].setSize(nFaces);
        rhoNBF_[j].setSize(nFaces, 0.0);
        rhoMBF_[j].setSize(nFaces, 0.0);
        linearKEBF_[j].setSize(nFaces, 0.0);
        momentumBF_[j].setSize(nFaces, vector::zero);
        ErotBF_[j].setSize(nFaces, 0.0);
        zetaRotBF_[j].setSize(nFaces, 0.0);
        qBF_[j].setSize(nFaces, 0.0);
        fDBF_[j].setSize(nFaces, vector::zero);
        zetaVibBF_[j].setSize(nFaces, 0.0);
        rhoNIntBF_[j].setSize(nFaces, 0.0);
        rhoNElecBF_[j].setSize(nFaces, 0.0);

        n_[j].setSize(nFaces, vector::zero);
        t1_[j].setSize(nFaces, vector::zero);
        t2_[j].setSize(nFaces, vector::zero);

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }

    calculateWallUnitVectors();

    speciesRhoNBF_.setSize(nSpecies);
    speciesMccBF_.setSize(nSpecies);
    speciesEvibBF_.setSize(nSpecies);
    speciesTvibBF_.setSize(nSpecies);
    speciesZetaVibBF_.setSize(nSpecies);
    speciesEvibModBF_.setSize(nSpecies);
    speciesEelecBF_.setSize(nSpecies);

    forAll(speciesIds_, i)
    {
        const label spId = speciesIds_[i];

        speciesRhoNBF_[i].setSize(nPatches);
        speciesMccBF_[i].setSize(nPatches);

        speciesEvibBF_[i].setSize(nPatches);
        speciesTvibBF_[i].setSize(nPatches);
        speciesZetaVibBF_[i].setSize(nPatches);
        speciesEvibModBF_[i].setSize
        (
            cloud_.constProps(spId).nVibrationalModes()
        );
        speciesEelecBF_[i].setSize(nPatches);

        forAll(speciesEvibBF_[i], j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            const label nFaces = patch.size();

            speciesRhoNBF_[i][j].setSize(nFaces, 0.0);
            speciesMccBF_[i][j].setSize(nFaces, 0.0);
            speciesEvibBF_[i][j].setSize(nFaces, 0.0);
            speciesTvibBF_[i][j].setSize(nFaces, 0.0);
            speciesZetaVibBF_[i][j].setSize(nFaces, 0.0);
            speciesEelecBF_[i][j].setSize(nFaces, 0.0);
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            speciesEvibModBF_[i][mod].setSize(nPatches);
            forAll(speciesEvibModBF_[i][mod], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                const label nFaces = patch.size();
                speciesEvibModBF_[i][mod][j].setSize(nFaces, 0.0);
            }
        }
    }

    sampleInterval_ = propsDict_.lookupOrDefault("sampleInterval", 1);

    measureClassifications_ =
        propsDict_.lookupOrDefault<bool>("measureClassifications", false);

    measureErrors_ = propsDict_.lookupOrDefault<bool>("measureErrors", false);

    densityOnly_ = propsDict_.lookupOrDefault<bool>("densityOnly", false);

    measureHeatFluxShearStress_ =
        propsDict_.lookupOrDefault<bool>("measureHeatFluxShearStress", false);

    measureMeanFreePath_ =
        propsDict_.lookupOrDefault<bool>("measureMeanFreePath", false);

    mfpTref_ =
        propsDict_.lookupOrDefault<scalar>("mfpReferenceTemperature", 273.0);

    averagingAcrossManyRuns_ =
        propsDict_.lookupOrDefault<bool>("averagingAcrossManyRuns", false);

    //- read in stored data from dictionary
    if (averagingAcrossManyRuns_)
    {
        if (!time_.resetFieldsAtOutput())
        {
            Info<< "Averaging across many runs for field " << fieldName_
                << " is enabled. Sampled data will be read from file."
                << endl;
            readIn();
        }
        else
        {
            Info<< "Averaging across many runs for field " << fieldName_
                << " will be enabled as soon as resetAtOutput is turned off."
                << endl;
        }
    }
}

void dsmcVolFields::calculateField()
{
    sampleCounter_++;

    const scalar kB = physicoChemical::k.value();
    const scalar NAvo = physicoChemical::NA.value();

    //- Reset instantaneous number of DSMC parcels
    dsmcN_ = 0.0;

    if (sampleInterval_ <= sampleCounter_)
    {
        nTimeSteps_ += 1.0;
        const scalar nAvTimeSteps = nTimeSteps_;

        if (densityOnly_)
        {
            //- Loop over the the entire parcel cloud
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label spId = findIndex(speciesIds_, p.typeId());

                //- Do not consider adsorbed parcels
                if (spId != -1 && p.isFree())
                {
                    const label cell = p.cell();
                    // 固定网格粒子数自适应启用后，同一个 cell 内的 parcel
                    // 不再必然代表相同数量的真实分子。
                    // 因此这里不能再使用 cell 级 cloud_.nParticles(cell)，
                    // 而必须使用 parcel 级
                    //
                    //     N_real(parcel) = N_cell * RWF(parcel)
                    //
                    // 否则 split 后产生的大量低权重克隆 parcel 会被当成
                    // 全权重 parcel 统计，rhoN/rhoM/p 会被系统性放大。
                    const scalar nParticles = cloud_.representedMolecules(p);
                    const scalar mass = cloud_.constProps(p.typeId()).mass();

                    // cumulative number of DSMC parcels
                    dsmcNCum_[cell] += 1.0;
                    // instantaneous number of DSMC parcels in this time step
                    dsmcN_[cell] += 1.0;
                    // cumulative number of real particles
                    nCum_[cell] += nParticles;
                    // cumulative mass of real particles
                    mCum_[cell] += mass*nParticles;
                }
            }
        }
        else
        {
            //- Loop over the the entire parcel cloud
            forAllConstIter(dsmcCloud, cloud_, iter)
            {
                const dsmcParcel& p = iter();
                const label spId = findIndex(speciesIds_, p.typeId());

                //- Do not consider adsorbed parcels
                if (spId != -1 && p.isFree())
                {
                    const dsmcParcel::constantProperties& cP =
                        cloud_.constProps(p.typeId());

                    //- Local copies
                    const label cell = p.cell();
                    // Physical moments use represented molecules; simulator counts remain unweighted.
                    const scalar nParticles = cloud_.representedMolecules(p);
                    const scalar mp = cP.mass();
                    const vector& Up = p.U();
                    const scalar linearKE = mp*(Up & Up);
                    const scalar Erotp = p.ERot();
                    const scalar zetaRotp = cP.rotationalDegreesOfFreedom();

                    scalar Evibp = 0.0;
                    forAll(cP.thetaV(), mod)
                    {
                        const label vibLvl = p.vibLevel()[mod];
                        const scalar EvibMod = cP.eVib_m(mod, vibLvl);
                        dsmcSpeciesEvibModCum_[spId][mod][cell] +=
                            nParticles*EvibMod;
                        Evibp += EvibMod;
                    }

                    const label nElecLevels = cP.nElectronicLevels();
                    const scalarList& electronicEnergies =
                        cP.electronicEnergyList();
                    const scalar Eelecp = 0.0; // TODO
                    const scalar Eintp = Erotp + Evibp + Eelecp;

                    // cumulative number of DSMC parcels
                    dsmcNCum_[cell] += 1.0;
                    // instantaneous number of DSMC parcels in this time step
                    dsmcN_[cell] += 1.0;
                    // Physical moments use represented molecules, while raw
                    // parcel counts above remain occupancy diagnostics.
                    dsmcMCum_[cell] += nParticles*mp;
                    dsmcLinearKECum_[cell] += nParticles*linearKE;
                    dsmcMomentumCum_[cell] += nParticles*mp*Up;
                    dsmcErotCum_[cell] += nParticles*Erotp;
                    dsmcZetaRotCum_[cell] += nParticles*zetaRotp;
                    // cumulative electronic energy of the DSMC parcels of
                    // each species
                    dsmcSpeciesEelecCum_[spId][cell] +=
                        electronicEnergies[p.ELevel()];
                    // cumulative number of DSMC parcels of each species
                    dsmcNSpeciesCum_[spId][cell] += 1.0;
                    dsmcMccSpeciesCum_[spId][cell] += nParticles*linearKE;
                    // cumulative number of real particles of each species
                    nSpeciesCum_[spId][cell] += nParticles;
                    // cumulative number of real particles
                    nCum_[cell] += nParticles;
                    // cumulative mass of real particles
                    mCum_[cell] += mp*nParticles;
                    // cumulative momentum of real particles
                    momentumCum_[cell] += mp*Up*nParticles;
                    // cumulative linear kinetic energy of the real particles
                    linearKECum_[cell] += linearKE*nParticles;

                    // cumulative counters for mass times squared velocity
                    // components of the DSMC parcel
                    // notation: Up = c = (u v w)
                    dsmcMuuCum_[cell] += nParticles*mp*sqr(Up.x());
                    dsmcMuvCum_[cell] += nParticles*mp*Up.x()*Up.y();
                    dsmcMuwCum_[cell] += nParticles*mp*Up.x()*Up.z();
                    dsmcMvvCum_[cell] += nParticles*mp*sqr(Up.y());
                    dsmcMvwCum_[cell] += nParticles*mp*Up.y()*Up.z();
                    dsmcMwwCum_[cell] += nParticles*mp*sqr(Up.z());

                    dsmcMccCum_[cell] += nParticles*linearKE;
                    dsmcMccuCum_[cell] += nParticles*linearKE*Up.x();
                    dsmcMccvCum_[cell] += nParticles*linearKE*Up.y();
                    dsmcMccwCum_[cell] += nParticles*linearKE*Up.z();

                    // cumulative internal energy times velocity components of
                    // the DSMC parcel
                    dsmcEuCum_[cell] += nParticles*Eintp*Up.x();
                    dsmcEvCum_[cell] += nParticles*Eintp*Up.y();
                    dsmcEwCum_[cell] += nParticles*Eintp*Up.z();
                    dsmcECum_[cell] += nParticles*Eintp;

                    if (nElecLevels > 1)
                    {
                        dsmcNElecLvlCum_[cell] += 1.0;

                        if (p.ELevel() == 0)
                        {
                            dsmcNGrndElecLvlSpeciesCum_[spId][cell]++;
                        }
                        if (p.ELevel() == 1)
                        {
                            dsmcN1stElecLvlSpeciesCum_[spId][cell]++;
                        }
                    }

                    if (measureClassifications_)
                    {
                        const label classification = p.classification();

                        if (classification == 0)
                        {
                            dsmcNClassICum_[cell] += 1.0;
                        }

                        if (classification == 1)
                        {
                            dsmcNClassIICum_[cell] += 1.0;
                        }

                        if (classification == 2)
                        {
                            dsmcNClassIIICum_[cell] += 1.0;
                        }
                    }
                }
            }

            //- Loop over all cells
            forAll(dsmcNCum_, celli)
            {
                collisionSeparation_[celli] +=
                    cloud_.cellPropMeasurements().collisionSeparation()[celli];

                dsmcNCollsCum_[celli] +=
                    cloud_.cellPropMeasurements().nColls()[celli];

                if (dsmcNCum_[celli] > 1e-3)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;

                    const scalar rhoNMean = nCum_[celli]
                        /(nAvTimeSteps*cellVolume);
                    const scalar rhoMMean = mCum_[celli]
                        /(nAvTimeSteps*cellVolume);

                    rhoN_[celli] = rhoNMean;
                    rhoM_[celli] = rhoMMean;

                    UMean_[celli] = momentumCum_[celli]/mCum_[celli];
                    UAbsMean_[celli] = absoluteVelocityFromMean
                    (
                        cloud_,
                        celli,
                        mesh_.cellCentres()[celli],
                        UMean_[celli]
                    );
                    // 速度模是后处理最常用的标量诊断。
                    // UMag 对应当前求解变量 UMean_ 的模：
                    // - rotor/SRF 区域：|U_rel|；
                    // - stator/惯性区域：|U_abs|。
                    // UAbsMag 始终对应绝对速度模 |U_abs|，用于跨旋转区和
                    // 静止区进行统一物理比较。
                    UMagMean_[celli] = mag(UMean_[celli]);
                    UAbsMagMean_[celli] = mag(UAbsMean_[celli]);
                    // UThetaRel 只对真正处于旋转参考系定义内的单元有物理意义。
                    // 原因是：
                    // 1. rotorZone 内的 parcel 速度 UMean_[celli] 被解释为 U_rel；
                    // 2. statorZone 内的 parcel 速度 UMean_[celli] 被解释为 U_abs；
                    // 3. 若在 statorZone 也直接做 “UMean · e_theta”，
                    //    就会把一个绝对速度投影伪装成“相对系周向速度”，
                    //    使 UThetaRel 这个字段名与实际物理含义脱节。
                    //
                    // 因此这里显式限定：
                    // - 旋转区：输出真实的 u_theta,rel = U_rel · e_theta；
                    // - 非旋转区：写成 0，把跨 domain 对比统一交给 UThetaAbs。
                    if (cloud_.isInRotatingRegion(celli))
                    {
                        UThetaRelMean_[celli] =
                            relativeTangentialVelocityComponent
                            (
                                cloud_,
                                mesh_.cellCentres()[celli],
                                UMean_[celli]
                            );
                    }
                    else
                    {
                        UThetaRelMean_[celli] = 0.0;
                    }

                    // 直接输出绝对速度在柱坐标 e_theta 方向上的分量，
                    // 便于按不同 z 层检查周向速度是否出现反号。
                    UThetaAbsMean_[celli] = absoluteTangentialVelocityComponent
                    (
                        cloud_,
                        mesh_.cellCentres()[celli],
                        UAbsMean_[celli]
                    );

                    // 下面这组诊断场不再回答“速度长什么样”，而是直接回答
                    // “当前 cell 是否真的挂上了 SRF 力项”。
                    // 物理解释如下：
                    // 1. rotatingRegionMean_ 直接反映 cellZone 限域；
                    // 2. aCentrifugalMean_ 只由几何位置与 Omega 决定，
                    //    是判断 rotor 区是否真的启用旋转参考系最稳健的空间量；
                    // 3. aCoriolisMean_ 与 U_rel 相关，因此它更像“当前统计到的
                    //    宏观相对速度在 SRF 下会受到多大的速度耦合修正”；
                    // 4. aSRFMean_ 则把两者重新拼回 parcel 推进时实际用的总加速度。
                    const bool inRotatingRegion =
                        cloud_.isInRotatingRegion(celli);
                    rotatingRegionMean_[celli] =
                        inRotatingRegion ? 1.0 : 0.0;

                    if (inRotatingRegion)
                    {
                        aCoriolisMean_[celli] = coriolisAccelerationFromMean
                        (
                            cloud_,
                            UMean_[celli]
                        );
                        aCentrifugalMean_[celli] =
                            centrifugalAccelerationFromPosition
                            (
                                cloud_,
                                mesh_.cellCentres()[celli]
                            );
                        aSRFMean_[celli] =
                            aCoriolisMean_[celli] + aCentrifugalMean_[celli];
                    }
                    else
                    {
                        aCoriolisMean_[celli] = vector::zero;
                        aCentrifugalMean_[celli] = vector::zero;
                        aSRFMean_[celli] = vector::zero;
                    }

                    const scalar linearKEMean = 0.5*linearKECum_[celli]
                        /(cellVolume*nAvTimeSteps);

                    Ttra_[celli] =
                        2.0/(3.0*kB*rhoNMean)
                       *(
                            linearKEMean - 0.5*rhoMMean
                           *(
                                UMean_[celli] & UMean_[celli]
                            )
                        );

                    p_[celli] = rhoNMean*kB*Ttra_[celli];
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[celli] = 0.001;
                    rhoN_[celli] = 0.0;
                    rhoM_[celli] = 0.0;
                    UMean_[celli] = vector::zero;
                    UAbsMean_[celli] = vector::zero;
                    UMagMean_[celli] = 0.0;
                    UAbsMagMean_[celli] = 0.0;
                    UThetaRelMean_[celli] = 0.0;
                    UThetaAbsMean_[celli] = 0.0;
                    rotatingRegionMean_[celli] =
                        cloud_.isInRotatingRegion(celli) ? 1.0 : 0.0;
                    if (cloud_.isInRotatingRegion(celli))
                    {
                        // 即使该 cell 当前几乎没有统计样本，离心项仍然可以
                        // 直接作为“这里是否挂上 SRF 几何场”的指示量输出；
                        // 而科里奥利项因缺乏可靠 U_rel 统计，保守置零。
                        aCoriolisMean_[celli] = vector::zero;
                        aCentrifugalMean_[celli] =
                            centrifugalAccelerationFromPosition
                            (
                                cloud_,
                                mesh_.cellCentres()[celli]
                            );
                        aSRFMean_[celli] = aCentrifugalMean_[celli];
                    }
                    else
                    {
                        aCoriolisMean_[celli] = vector::zero;
                        aCentrifugalMean_[celli] = vector::zero;
                        aSRFMean_[celli] = vector::zero;
                    }
                    Ttra_[celli] = 0.0;
                    p_[celli] = 0.0;
                }
            }

            //- Obtain boundary measurements
            forAll(speciesIds_, i)
            {
                const label spId = speciesIds_[i];

                forAll(mesh_.boundaryMesh(), j)
                {
                    forAll(mesh_.boundaryMesh()[j], k)
                    {
                        rhoNBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNBF(spId, j, k);
                        rhoMBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoMBF(spId, j, k);
                        linearKEBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesLinearKEBF(spId, j, k);
                        momentumBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesMomentumBF(spId, j, k);
                        ErotBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesErotBF(spId, j, k);
                        zetaRotBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesZetaRotBF(spId, j, k);
                        rhoNIntBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNIntBF(spId, j, k);
                        rhoNElecBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNElecBF(spId, j, k);
                        qBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesqBF(spId, j, k);
                        fDBF_[j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesfDBF(spId, j, k);

                        speciesRhoNBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesRhoNBF(spId, j, k);
                        speciesEvibBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesEvibBF(spId, j, k);
                        speciesEelecBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesEelecBF(spId, j, k);
                        speciesMccBF_[i][j][k] +=
                            cloud_.boundaryFluxMeasurements()
                              .speciesMccBF(spId, j, k);
                    }
                }

                forAll(speciesEvibModBF_[i], mod)
                {
                    forAll(mesh_.boundaryMesh(), j)
                    {
                        forAll(mesh_.boundaryMesh()[j], k)
                        {
                            speciesEvibModBF_[i][mod][j][k] +=
                                cloud_.boundaryFluxMeasurements()
                                  .speciesEvibModBF(spId, mod, j, k);
                        }
                    }
                }
            }
        }

        sampleCounter_ = 0;
    }

    if (time_.time().outputTime())
    {
        const scalar nAvTimeSteps = nTimeSteps_;

        if (densityOnly_)
        {
            forAll(dsmcNCum_, celli)
            {
                if (dsmcNCum_[celli] > SMALL)
                {
                    const scalar cellVolume = mesh_.cellVolumes()[celli];

                    dsmcNMean_[celli] = dsmcNCum_[celli]/nAvTimeSteps;

                    rhoN_[celli] = nCum_[celli]/(nAvTimeSteps*cellVolume);
                    rhoM_[celli] = mCum_[celli]/(nAvTimeSteps*cellVolume);
                }
                else
                {
                    // not zero so that weighted decomposition still works
                    dsmcNMean_[celli] = 0.001;
                    rhoN_[celli] = 0.0;
                    rhoM_[celli] = 0.0;
                }

                if (dsmcN_[celli] < SMALL)
                {
                    // not zero so that weighted decomposition still works
                    dsmcN_[celli] = 0.001;
                }
            }
        }
        else
        {
            const label nSpecies = speciesIds_.size();

            forAll(dsmcNCum_, celli)
            {
                //- Fields initialisation
                scalar moleculesRhoN = 0.0;
                Tvib_[celli] = 0.0;
                scalarList speciesTvib(nSpecies, 0.0);
                List<scalarList> speciesTvibMod(nSpecies);
                scalarList speciesZetaVib(nSpecies, 0.0);
                List<scalarList> speciesZetaVibMod(nSpecies);

                scalar molarCv_trarot = 0.0;
                scalar molarCp_trarot = 0.0;
                scalar molecularMass = 0.0;
                scalar particleCv = 0.0;
                scalar gamma = 0.0;

                //- Rotational energy mode
                const scalar zetaRotTot
                (
                    nCum_[celli] > 0.0
                  ? dsmcZetaRotCum_[celli]/nCum_[celli]
                  : 0.0
                );

                Trot_[celli] =
                (
                    dsmcZetaRotCum_[celli] > 0.0
                  ? 2.0*dsmcErotCum_[celli]/(kB*dsmcZetaRotCum_[celli])
                  : 0.0
                );

                //- Vibrational energy mode
                forAll(speciesIds_, i)
                {
                    const label spId = speciesIds_[i];
                    const label nVibMod =
                        cloud_.constProps(spId).nVibrationalModes();

                    speciesZetaVibMod[i].setSize(nVibMod, 0.0);
                    speciesTvibMod[i].setSize(nVibMod, 0.0);
                    scalar zetaByTvibMod = 0.0;

                    forAll(dsmcSpeciesEvibModCum_[i], mod)
                    {
                        if
                        (
                            dsmcSpeciesEvibModCum_[i][mod][celli] > VSMALL
                         && nSpeciesCum_[i][celli] > 0.0
                         && speciesZetaVibMod.size() > SMALL
                        )
                        {
                            const scalar thetaV =
                                cloud_.constProps(spId).thetaV_m(mod);

                            const scalar iMean =
                                dsmcSpeciesEvibModCum_[i][mod][celli]
                               /(kB*thetaV*nSpeciesCum_[i][celli]);

                            if (iMean > SMALL)
                            {
                                const scalar logFactor = log(1.0 + 1.0/iMean);

                                speciesTvibMod[i][mod] = thetaV/logFactor;

                                speciesZetaVibMod[i][mod] = 2.0*iMean*logFactor;

                                speciesZetaVib[i] += speciesZetaVibMod[i][mod];

                                zetaByTvibMod = speciesZetaVibMod[i][mod]
                                    *speciesTvibMod[i][mod];
                            }
                        }
                    }

                    if (speciesZetaVib[i] > SMALL)
                    {
                        moleculesRhoN += nSpeciesCum_[i][celli];

                        speciesTvib[i] = zetaByTvibMod/speciesZetaVib[i];

                        Tvib_[celli] += nSpeciesCum_[i][celli]*speciesTvib[i];

                        zetaVib_[celli] += nSpeciesCum_[i][celli]
                            *speciesZetaVib[i];
                    }

                } //- end species loop

                if (moleculesRhoN > SMALL)
                {
                    Tvib_[celli] /= moleculesRhoN;
                    zetaVib_[celli] /= moleculesRhoN;
                }

                //- Electronic energy mode // TODO Vincent
                //  To reintroduce - I do not trust this part
                scalar zetaElecTot = 0.0;
                Telec_[celli] = 0.0;

                //- Overall temperature
                Tov_[celli] =
                    (
                        3.0*Ttra_[celli]
                      + zetaRotTot*Trot_[celli]
                      + zetaVib_[celli]*Tvib_[celli]
                      + zetaElecTot*Telec_[celli]
                    ) /
                    (3.0 + zetaRotTot + zetaVib_[celli] + zetaElecTot);

                if (measureHeatFluxShearStress_)
                {
                    if (nCum_[celli] > 0.0)
                    {
                        pressureTensor_[celli].xx() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMuuCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].x())
                            );
                        pressureTensor_[celli].xy() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMuvCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].x()
                              * UMean_[celli].y()
                            );
                        pressureTensor_[celli].xz() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMuwCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].x()
                              * UMean_[celli].z()
                            );

                        pressureTensor_[celli].yx() =
                            pressureTensor_[celli].xy();
                        pressureTensor_[celli].yy() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMvvCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].y())
                            );
                        pressureTensor_[celli].yz() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMvwCum_[celli]
                              - dsmcMCum_[celli]*UMean_[celli].y()
                              * UMean_[celli].z()
                            );

                        pressureTensor_[celli].zx() =
                            pressureTensor_[celli].xz();
                        pressureTensor_[celli].zy() =
                            pressureTensor_[celli].yz();
                        pressureTensor_[celli].zz() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                dsmcMwwCum_[celli]
                              - dsmcMCum_[celli]*sqr(UMean_[celli].z())
                            );

                        const scalar scalarPressure =
                            1.0/3.0
                           *(
                                pressureTensor_[celli].xx()
                              + pressureTensor_[celli].yy()
                              + pressureTensor_[celli].zz()
                            );

                        shearStressTensor_[celli] = -pressureTensor_[celli];
                        shearStressTensor_[celli].xx() += scalarPressure;
                        shearStressTensor_[celli].yy() += scalarPressure;
                        shearStressTensor_[celli].zz() += scalarPressure;

                        //- terms involving pressure tensor should not be
                        //  multiplied by the number density
                        //  (see Bird corrigendum)

                        heatFluxVector_[celli].x() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                0.5*dsmcMccuCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].x()
                              + dsmcEuCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].x()
                            )
                          - pressureTensor_[celli].xx()*UMean_[celli].x()
                          - pressureTensor_[celli].xy()*UMean_[celli].y()
                          - pressureTensor_[celli].xz()*UMean_[celli].z();

                        heatFluxVector_[celli].y() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                0.5*dsmcMccvCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].y()
                              + dsmcEvCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].y()
                            )
                          - pressureTensor_[celli].yx()*UMean_[celli].x()
                          - pressureTensor_[celli].yy()*UMean_[celli].y()
                          - pressureTensor_[celli].yz()*UMean_[celli].z();

                        heatFluxVector_[celli].z() =
                            rhoN_[celli]/nCum_[celli]
                           *(
                                0.5*dsmcMccwCum_[celli]
                              - 0.5*dsmcMccCum_[celli]*UMean_[celli].z()
                              + dsmcEwCum_[celli]
                              - dsmcECum_[celli]*UMean_[celli].z()
                            )
                          - pressureTensor_[celli].zx()*UMean_[celli].x()
                          - pressureTensor_[celli].zy()*UMean_[celli].y()
                          - pressureTensor_[celli].zz()*UMean_[celli].z();
                    }
                    else
                    {
                        pressureTensor_[celli] = tensor::zero;
                        shearStressTensor_[celli] = tensor::zero;
                        heatFluxVector_[celli] = vector::zero;
                    }
                }

                if (dsmcNCum_[celli] > SMALL and Ttra_[celli] > SMALL)
                {
                    forAll(speciesIds_, i)
                    {
                        const label spId = speciesIds_[i];
                        const scalar speciesZetaRot =
                            cloud_.constProps(spId)
                              .rotationalDegreesOfFreedom();

                        const scalar Xs = nSpeciesCum_[i][celli]
                            /nCum_[celli];

                        molecularMass += Xs*cloud_.constProps(spId).mass();

                         //- Heat capacity at constant volume/(0.5*kB)
                         //  trans-rotational
                         molarCv_trarot += Xs*(3.0 + speciesZetaRot);

                        //- Heat capacity at constant pressure/(0.5*kB)
                        //  trans-rotational
                        molarCp_trarot += Xs*(5.0 + speciesZetaRot);
                    }

                    particleCv = molarCv_trarot/NAvo;

                    gamma = molarCp_trarot/molarCv_trarot;

                    const scalar speedOfSound = sqrt
                        (
                            gamma*kB/molecularMass*Ttra_[celli]
                        );

                    Ma_[celli] = mag(UMean_[celli])/speedOfSound;
                }
                else
                {
                    Ma_[celli] = 0.0;
                }

                if (measureMeanFreePath_ && Ttra_[celli] > 1.0)
                {
                    const scalar deltaT = cloud_.deltaTValue(celli);

                    mfp_[celli] = 0.0;
                    meanCollisionRate_[celli] = 0.0;

                    forAll(speciesIds_, s)
                    {
                        const label spIdp = speciesIds_[s];

                        speciesMfp_[s][celli] = 0.0;
                        speciesMcr_[s][celli] = 0.0;

                        forAll(speciesIds_, r)
                        {
                            const label spIdq = speciesIds_[r];

                            const scalar dPQ =
                                0.5*
                                (
                                    cloud_.constProps(spIdp).d()
                                  + cloud_.constProps(spIdq).d()
                                );
                            const scalar omegaPQ =
                                0.5*
                                (
                                    cloud_.constProps(spIdp).omega()
                                  + cloud_.constProps(spIdq).omega()
                                );
                            const scalar massRatio =
                                cloud_.constProps(spIdp).mass()
                               /cloud_.constProps(spIdq).mass();

                            if
                            (
                                dsmcNSpeciesCum_[r][celli] > SMALL
                             && Ttra_[celli] > SMALL
                            )
                            {
                                const scalar nDensQ =
                                    nSpeciesCum_[r][celli]
                                   /(mesh_.cellVolumes()[celli]*nAvTimeSteps);
                                const scalar reducedMass =
                                    cloud_.constProps(spIdp).mass()
                                   *cloud_.constProps(spIdq).mass()
                                   /
                                    (
                                       cloud_.constProps(spIdp).mass()
                                     + cloud_.constProps(spIdq).mass()
                                    );

                                // Bird 1994, eq (4.76)
                                speciesMfp_[s][celli] += pi*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        mfpTref_/Ttra_[celli], omegaPQ - 0.5
                                    )*sqrt(1.0+massRatio);

                                // Bird 1994, eq (4.74)
                                speciesMcr_[s][celli] +=
                                    2.0*sqrt(pi)*sqr(dPQ)*nDensQ
                                   *pow
                                    (
                                        Ttra_[celli]/mfpTref_, 1.0 - omegaPQ
                                    )
                                   *sqrt
                                    (
                                        2.0*kB*mfpTref_/reducedMass
                                    );
                            }
                        }

                        if (speciesMfp_[s][celli] > SMALL)
                        {
                            speciesMfp_[s][celli] = 1.0/speciesMfp_[s][celli];
                        }
                    }

                    meanCollisionSeparation_[celli] =
                    (
                        dsmcNCollsCum_[celli] > SMALL
                      ? collisionSeparation_[celli]/dsmcNCollsCum_[celli]
                      : GREAT
                    );

                    if (nCum_[celli] > SMALL)
                    {
                        // const scalar symmFactor = 2.0;
                        // TODO (s == r ? 1.0 : 2.0);
                        measuredCollisionRate_[celli] = dsmcNCollsCum_[celli]
                            *cloud_.nParticles(celli)/(nCum_[celli]*deltaT);
                    }

                    if (rhoN_[celli] > SMALL)
                    {
                        forAll(speciesIds_, i)
                        {
                            const scalar rhoNi = nSpeciesCum_[i][celli];

                            // Bird 1994, eq (4.77)
                            mfp_[celli] += speciesMfp_[i][celli]
                                *rhoNi/nCum_[celli];

                            // Bird 1994, eq (1.38)
                            meanCollisionRate_[celli] +=
                                speciesMcr_[i][celli]*rhoNi/nCum_[celli];
                        }
                    }

                    if (mfp_[celli] < SMALL)
                    {
                        mfp_[celli] = GREAT;
                    }

                    if (meanCollisionRate_[celli] > SMALL)
                    {
                        meanCollisionTime_[celli] =
                            1.0/meanCollisionRate_[celli];
                        mctToDt_[celli] = meanCollisionTime_[celli]/deltaT;
                    }
                    else
                    {
                        meanCollisionTime_[celli] = GREAT;
                        mctToDt_[celli] = GREAT;
                    }

                    if (mfp_[celli] != GREAT)
                    {
                        scalar maxCellDx = 0.0;
                        scalarField cellDx(3, 0.0);

                        const labelList& pLabels
                        (
                            mesh_.cells()[celli].labels(mesh_.faces())
                        );
                        pointField pLocal(pLabels.size(), vector::zero);

                        forAll (pLabels, pointi)
                        {
                            pLocal[pointi] = mesh_.points()[pLabels[pointi]];
                        }

                        cellDx[0] = Foam::max(pLocal & vector(1,0,0))
                            - Foam::min(pLocal & vector(1,0,0));
                        cellDx[1] = Foam::max(pLocal & vector(0,1,0))
                            - Foam::min(pLocal & vector(0,1,0));
                        cellDx[2] = Foam::max(pLocal & vector(0,0,1))
                            - Foam::min(pLocal & vector(0,0,1));

                        maxCellDx = cellDx[0];

                        forAll(cellDx, dim)
                        {
                            if (cellDx[dim] > maxCellDx)
                            {
                                maxCellDx = cellDx[dim];
                            }
                        }

                        mfpToDx_[celli] = mfp_[celli]/maxCellDx;

                        SOF_[celli] =
                        (
                            mfp_[celli] > SMALL
                          ? meanCollisionSeparation_[celli]/mfp_[celli]
                          : 0.0
                        );
                    }
                    else
                    {
                        mfpToDx_[celli] = GREAT;
                        SOF_[celli] = GREAT;
                    }

                    // when few particles in cell, undesired refinement
                    // this condition should eliminates this problem
                    if (dsmcN_[celli] >= 4.0)
                    {
                        DxToMfp_[celli] = 1.0/mfpToDx_[celli];
                    }
                }

                if (measureClassifications_)
                {
                    if (dsmcNCum_[celli] > SMALL)
                    {
                        classIDistribution_[celli] = dsmcNClassICum_[celli]
                            /dsmcNCum_[celli];
                        classIIDistribution_[celli] = dsmcNClassIICum_[celli]
                            /dsmcNCum_[celli];
                        classIIIDistribution_[celli] = dsmcNClassIIICum_[celli]
                            /dsmcNCum_[celli];
                    }
                }

                if (measureErrors_)
                {
                    if
                    (
                         dsmcNMean_[celli] > SMALL && Ma_[celli] > SMALL
                      && gamma > SMALL && particleCv > SMALL
                    )
                    {
                        const scalar deno = sqrt(dsmcNMean_[celli]*nAvTimeSteps);

                        densityError_[celli] = 1.0/deno;
                        velocityError_[celli] = 1.0/(deno*Ma_[celli]*sqrt(gamma));
                        temperatureError_[celli] = sqrt(kB/particleCv)/deno;
                        pressureError_[celli] = sqrt(gamma)/deno;
                    }

                }
            } //- end loop over cells

            //- Computing boundary measurements: loop over all boundary patches
            forAll(rhoNBF_, j)
            {
                //- Determine of the type of patch: patch, wall, cyclic, ...
                const polyPatch& patch = mesh_.boundaryMesh()[j];

                const bool isWall = isA<wallPolyPatch>(patch);

                const bool isNonEmptyNonCyclic = isA<polyPatch>(patch)
                    && !isA<emptyPolyPatch>(patch)
                    && !isA<cyclicPolyPatch>(patch);

                if (isWall)
                {
                    //- Loop over all wall boundary faces
                    forAll(patch, k)
                    {
                        const label celli = boundaryCells_[j][k];

                        //- Initialise face fields
                        Tvib_.boundaryFieldRef()[j][k] = 0.0;
                        zetaVibBF_[j][k] = 0.0;
                        scalar molecularMassBF = 0.0;
                        scalar molarCvBF_trarot = 0.0;
                        scalar molarCpBF_trarot = 0.0;

                        // Note: do not use the nParticles value that includes
                        // the RWF here. This is wrong because boundary
                        // measurements are performend during move steps. Hence
                        // radial weighting (if simulation is axisymmetric or
                        // spherical) is performed after the measurement. That
                        // is why the parcel RWFs are already included during
                        // the measurement step and we only need the FNUM value
                        // here.
                        const scalar nParticles = cloud_.coordSystem().dtModel()
                            .nParticles(j, k);

                        const scalar rhoNMean =
                            rhoNBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar rhoMMean =
                            rhoMBF_[j][k]*nParticles/nAvTimeSteps;
                        const scalar linearKEMean =
                            linearKEBF_[j][k]*nParticles/nAvTimeSteps;

                        rhoN_.boundaryFieldRef()[j][k] = rhoNMean;
                        rhoM_.boundaryFieldRef()[j][k] = rhoMMean;

                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_.boundaryFieldRef()[j][k] = dsmcN_[celli];
                        dsmcNMean_.boundaryFieldRef()[j][k] = dsmcNMean_[celli];

                        //- Translational energy mode and velocity
                        if (rhoMMean > VSMALL)
                        {
                            UMean_.boundaryFieldRef()[j][k] = momentumBF_[j][k]
                                /rhoMBF_[j][k];

                            // 只有当当前 wall face 确实积累到了有效动量样本时，
                            // 才允许把相对速度转回绝对速度。
                            // 否则若 UMean 只是“因为无样本被置零”，再机械地加上
                            // Omega x r，就会凭空在壁面制造出一个刚体转动速度，
                            // 误导后处理把“无统计样本”看成“有真实周向流”。
                            UAbsMean_.boundaryFieldRef()[j][k] =
                                absoluteVelocityFromMean
                                (
                                    cloud_,
                                    celli,
                                    mesh_.boundaryMesh()[j].faceCentres()[k],
                                    UMean_.boundaryField()[j][k]
                                );
                            // 边界速度模与内部 cell 使用同一口径：
                            // UMag 是当前 UMean 的模，UAbsMag 是绝对速度模。
                            UMagMean_.boundaryFieldRef()[j][k] =
                                mag(UMean_.boundaryField()[j][k]);
                            UAbsMagMean_.boundaryFieldRef()[j][k] =
                                mag(UAbsMean_.boundaryField()[j][k]);
                            // 边界统计与内部单元保持同一物理口径：
                            // 只有当 owner cell 位于旋转区时，当前 wall/patch
                            // 采样到的 UMean 才能被解释成相对速度 U_rel。
                            // 非旋转区边界若继续输出 “UMean · e_theta”，
                            // 会把 stator 绝对速度错误包装成 UThetaRel。
                            if (cloud_.isInRotatingRegion(celli))
                            {
                                UThetaRelMean_.boundaryFieldRef()[j][k] =
                                    relativeTangentialVelocityComponent
                                    (
                                        cloud_,
                                        mesh_.boundaryMesh()[j].faceCentres()[k],
                                        UMean_.boundaryField()[j][k]
                                    );
                            }
                            else
                            {
                                UThetaRelMean_.boundaryFieldRef()[j][k] = 0.0;
                            }

                            UThetaAbsMean_.boundaryFieldRef()[j][k] =
                                absoluteTangentialVelocityComponent
                                (
                                    cloud_,
                                    mesh_.boundaryMesh()[j].faceCentres()[k],
                                    UAbsMean_.boundaryField()[j][k]
                                );

                            // 对 wall/采样边界，同样输出 SRF 诊断量。
                            // 这里沿用 owner cell 的旋转区语义：
                            // - 若 owner cell 在 rotorZone，则这张边界 face
                            //   被视为“属于旋转区一侧”的采样位置；
                            // - 若 owner cell 在 statorZone，则相关 SRF 场
                            //   全部显式写零。
                            const bool boundaryInRotatingRegion =
                                cloud_.isInRotatingRegion(celli);
                            rotatingRegionMean_.boundaryFieldRef()[j][k] =
                                boundaryInRotatingRegion ? 1.0 : 0.0;

                            if (boundaryInRotatingRegion)
                            {
                                aCoriolisMean_.boundaryFieldRef()[j][k] =
                                    coriolisAccelerationFromMean
                                    (
                                        cloud_,
                                        UMean_.boundaryField()[j][k]
                                    );
                                aCentrifugalMean_.boundaryFieldRef()[j][k] =
                                    centrifugalAccelerationFromPosition
                                    (
                                        cloud_,
                                        mesh_.boundaryMesh()[j].faceCentres()[k]
                                    );
                                aSRFMean_.boundaryFieldRef()[j][k] =
                                    aCoriolisMean_.boundaryField()[j][k]
                                  + aCentrifugalMean_.boundaryField()[j][k];
                            }
                            else
                            {
                                aCoriolisMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                                aCentrifugalMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                                aSRFMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                            }

                            Ttra_.boundaryFieldRef()[j][k] =
                                2.0/(3.0*kB*rhoNMean)*
                                (
                                    linearKEMean - 0.5*rhoMMean*
                                    (
                                        UMean_.boundaryField()[j][k]
                                      & UMean_.boundaryField()[j][k]
                                    )
                                );
                        }
                        else
                        {
                            UMean_.boundaryFieldRef()[j][k] = vector::zero;
                            UAbsMean_.boundaryFieldRef()[j][k] = vector::zero;
                            UMagMean_.boundaryFieldRef()[j][k] = 0.0;
                            UAbsMagMean_.boundaryFieldRef()[j][k] = 0.0;
                            UThetaRelMean_.boundaryFieldRef()[j][k] = 0.0;
                            UThetaAbsMean_.boundaryFieldRef()[j][k] = 0.0;
                            rotatingRegionMean_.boundaryFieldRef()[j][k] =
                                cloud_.isInRotatingRegion(celli) ? 1.0 : 0.0;
                            if (cloud_.isInRotatingRegion(celli))
                            {
                                aCoriolisMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                                aCentrifugalMean_.boundaryFieldRef()[j][k] =
                                    centrifugalAccelerationFromPosition
                                    (
                                        cloud_,
                                        mesh_.boundaryMesh()[j].faceCentres()[k]
                                    );
                                aSRFMean_.boundaryFieldRef()[j][k] =
                                    aCentrifugalMean_.boundaryField()[j][k];
                            }
                            else
                            {
                                aCoriolisMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                                aCentrifugalMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                                aSRFMean_.boundaryFieldRef()[j][k] =
                                    vector::zero;
                            }
                            Ttra_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Rotational energy mode
                        const scalar zetaRotTot =
                        (
                            rhoNBF_[j][k] > SMALL
                          ? zetaRotBF_[j][k]/rhoNBF_[j][k]
                          : 0.0
                        );

                        Trot_.boundaryFieldRef()[j][k] =
                        (
                            zetaRotBF_[j][k] > SMALL
                          ? 2.0*ErotBF_[j][k]/(kB*zetaRotBF_[j][k])
                          : 0.0
                        );

                        //- Vibrational energy mode: loop over all species
                        scalar moleculesRhoN = 0.0;

                        forAll(speciesIds_, i)
                        {
                            const label spId = speciesIds_[i];
                            const label nVibMod =
                                cloud_.constProps(spId).nVibrationalModes();

                            speciesZetaVibBF_[i][j][k] = 0.0;
                            speciesTvibBF_[i][j][k] = 0.0;

                            scalar zetaByTvibMod = 0.0;
                            scalarList speciesZetaVibMod(nVibMod, 0.0);
                            scalarList speciesTvibMod(nVibMod, 0.0);

                            if (speciesRhoNBF_[i][j][k] > SMALL)
                            {
                                forAll(speciesZetaVibMod, mod)
                                {
                                    const scalar thetaV =
                                        cloud_.constProps(spId).thetaV()[mod];

                                    const scalar iMean =
                                        speciesEvibModBF_[i][mod][j][k]
                                       /(kB*thetaV*speciesRhoNBF_[i][j][k]);

                                    if (iMean > SMALL)
                                    {
                                        const scalar logFactor =
                                            log(1.0 + 1.0/iMean);

                                        speciesTvibMod[mod] = thetaV/logFactor;

                                        speciesZetaVibMod[mod] =
                                            2.0*iMean*logFactor;

                                        speciesZetaVibBF_[i][j][k] +=
                                            speciesZetaVibMod[mod];

                                        zetaByTvibMod += speciesZetaVibMod[mod]
                                            *speciesTvibMod[mod];
                                    }
                                }
                            }

                            if (speciesZetaVibBF_[i][j][k] > SMALL)
                            {
                                moleculesRhoN += speciesRhoNBF_[i][j][k];

                                speciesTvibBF_[i][j][k] = zetaByTvibMod
                                    /speciesZetaVibBF_[i][j][k];

                                Tvib_.boundaryFieldRef()[j][k] +=
                                    speciesRhoNBF_[i][j][k]
                                   *speciesTvibBF_[i][j][k];

                                zetaVibBF_[j][k] +=
                                    speciesRhoNBF_[i][j][k]
                                   *speciesZetaVibBF_[i][j][k];
                            }
                        }

                        if (moleculesRhoN > SMALL)
                        {
                            Tvib_.boundaryFieldRef()[j][k] /= moleculesRhoN;
                            zetaVibBF_[j][k] /= moleculesRhoN;
                        }

                        //- Electronic energy mode // TODO Vincent
                        //  Removed temporarily - I don't trust this part
                        scalar zetaElecTot = 0.0;
                        Telec_.boundaryFieldRef()[j][k] = 0.0;

                        Tov_.boundaryFieldRef()[j][k] =
                            (
                                (3.0*Ttra_.boundaryField()[j][k])
                              + (zetaRotTot*Trot_.boundaryField()[j][k])
                              + (
                                    zetaVibBF_[j][k]
                                   *Tvib_.boundaryField()[j][k]
                                )
                              + (
                                    zetaElecTot
                                   *Telec_.boundaryFieldRef()[j][k]
                                )
                            )
                           /(
                                3.0 + zetaRotTot + zetaVibBF_[j][k]
                              + zetaElecTot
                            );

                        if (rhoNBF_[j][k] > SMALL)
                        {
                            //- Loop over all species
                            forAll(speciesIds_, i)
                            {
                                const label spId = speciesIds_[i];

                                const scalar speciesZetaRotBF =
                                    cloud_.constProps(spId)
                                      .rotationalDegreesOfFreedom();

                                const scalar Xs =
                                    speciesRhoNBF_[i][j][k]/rhoNBF_[j][k];

                                molecularMassBF += Xs
                                    *cloud_.constProps(spId).mass();

                                //- Heat capacity at constant volume/(0.5*kB)
                                //  trans-rotational energy mode
                                molarCvBF_trarot += Xs*(3.0 + speciesZetaRotBF);

                                //- Heat capacity at constant pressure/(0.5*kB)
                                //  trans-rotational energy mode
                                molarCpBF_trarot += Xs*(5.0 + speciesZetaRotBF);
                            }

                            //- Mach number calculation
                            const scalar gasConstant =
                                (
                                    rhoNBF_[j][k] > SMALL
                                 && Ttra_.boundaryFieldRef()[j][k] > SMALL
                                  ? kB/molecularMassBF
                                  : 0.0
                                );

                            const scalar gamma = molarCpBF_trarot/molarCvBF_trarot;

                            const scalar speedOfSound =
                                sqrt
                                (
                                    gamma*gasConstant*Ttra_.boundaryField()[j][k]
                                );

                            Ma_.boundaryFieldRef()[j][k] =
                                mag(UMean_.boundaryField()[j][k])/speedOfSound;
                        }
                        else
                        {
                            Ma_.boundaryFieldRef()[j][k] = 0.0;
                        }

                        //- Force density
                        fD_.boundaryFieldRef()[j][k] = fDBF_[j][k]/nAvTimeSteps;

                        //- Surface pressure
                        p_.boundaryFieldRef()[j][k] =
                            fD_.boundaryField()[j][k] & n_[j][k];

                        //- Wall shear stress
                        tau_.boundaryFieldRef()[j][k] =
                            sqrt
                            (
                                sqr(fD_.boundaryField()[j][k] & t1_[j][k])
                              + sqr(fD_.boundaryField()[j][k] & t2_[j][k])
                            );

                        //- Heat flux
                        q_.boundaryFieldRef()[j][k] = qBF_[j][k]/nAvTimeSteps;

                        //- ZeroGradient condition assumed for Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_.boundaryFieldRef()[j][k] = mfp_[celli];
                            SOF_.boundaryFieldRef()[j][k] = SOF_[celli];
                            mfpToDx_.boundaryFieldRef()[j][k] = mfpToDx_[celli];
                            meanCollisionRate_.boundaryFieldRef()[j][k] =
                                meanCollisionRate_[celli];
                            meanCollisionTime_.boundaryFieldRef()[j][k] =
                                meanCollisionTime_[celli];
                            mctToDt_.boundaryFieldRef()[j][k] = mctToDt_[celli];
                        }
                    }
                }
                else if (isNonEmptyNonCyclic)
                {
                    //- Loop over all boundary faces and set zeroGradient
                    //  conditions
                    forAll(boundaryCells_[j], k)
                    {
                        const label celli = boundaryCells_[j][k];

                        //- Instantaneous and sampled numbers of DSMC parcels
                        //  are that of the neighbouring cell
                        dsmcN_.boundaryFieldRef()[j][k] = dsmcN_[celli];
                        dsmcNMean_.boundaryFieldRef()[j][k] =
                            dsmcNMean_[celli];

                        //- Number density and mass density fields
                        rhoN_.boundaryFieldRef()[j][k] = rhoN_[celli];
                        rhoM_.boundaryFieldRef()[j][k] = rhoM_[celli];

                        //- Temperature fields
                        Ttra_.boundaryFieldRef()[j][k] = Ttra_[celli];
                        Trot_.boundaryFieldRef()[j][k] = Trot_[celli];
                        Tvib_.boundaryFieldRef()[j][k] = Tvib_[celli];
                        Tov_.boundaryFieldRef()[j][k] = Tov_[celli];

                        //- Pressure, Mach and velocity fields
                        p_.boundaryFieldRef()[j][k] = p_[celli];
                        Ma_.boundaryFieldRef()[j][k] = Ma_[celli];
                        UMean_.boundaryFieldRef()[j][k] = UMean_[celli];
                        // 对 zeroGradient/普通 patch，这里直接沿用单元内部已经算好的
                        // UAbsMean_，不再按 face centre 二次重算。
                        // 这样可以同时避免两类问题：
                        // 1. 当内部单元根本没有有效样本、UAbsMean_[celli]=0 时，
                        //    不会在边界再凭空叠加一个 Omega x r；
                        // 2. 边界值与内部单元值保持完全一致，符合这一路径的
                        //    zeroGradient 语义。
                        UAbsMean_.boundaryFieldRef()[j][k] = UAbsMean_[celli];
                        UMagMean_.boundaryFieldRef()[j][k] =
                            UMagMean_[celli];
                        UAbsMagMean_.boundaryFieldRef()[j][k] =
                            UAbsMagMean_[celli];
                        UThetaRelMean_.boundaryFieldRef()[j][k] =
                            UThetaRelMean_[celli];
                        UThetaAbsMean_.boundaryFieldRef()[j][k] =
                            UThetaAbsMean_[celli];
                        rotatingRegionMean_.boundaryFieldRef()[j][k] =
                            rotatingRegionMean_[celli];
                        aCoriolisMean_.boundaryFieldRef()[j][k] =
                            aCoriolisMean_[celli];
                        aCentrifugalMean_.boundaryFieldRef()[j][k] =
                            aCentrifugalMean_[celli];
                        aSRFMean_.boundaryFieldRef()[j][k] =
                            aSRFMean_[celli];

                        //- Optional fields
                        if (measureMeanFreePath_)
                        {
                            mfp_.boundaryFieldRef()[j][k] = mfp_[celli];
                            SOF_.boundaryFieldRef()[j][k] = SOF_[celli];
                            mfpToDx_.boundaryFieldRef()[j][k] = mfpToDx_[celli];
                            meanCollisionRate_.boundaryFieldRef()[j][k] =
                                meanCollisionRate_[celli];
                            meanCollisionTime_.boundaryFieldRef()[j][k] =
                                meanCollisionTime_[celli];
                            mctToDt_.boundaryFieldRef()[j][k] = mctToDt_[celli];
                        }

                        if (measureHeatFluxShearStress_)
                        {
                            shearStressTensor_.boundaryFieldRef()[j][k] =
                                shearStressTensor_[celli];
                            heatFluxVector_.boundaryFieldRef()[j][k] =
                                heatFluxVector_[celli];
                            pressureTensor_.boundaryFieldRef()[j][k] =
                                pressureTensor_[celli];
                        }

                        if (measureClassifications_)
                        {
                            classIDistribution_.boundaryFieldRef()[j][k] =
                                classIDistribution_[celli];
                            classIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIDistribution_[celli];
                            classIIIDistribution_.boundaryFieldRef()[j][k] =
                                classIIIDistribution_[celli];
                        }
                    }
                }
            }

            //- Write solution fields
            // outputProperty() 在未配置 outputProperties 时总是返回 true，
            // 因此同一条写出路径即可同时覆盖旧算例的全量输出和新算例的白名单输出。
            {
                writeSelected(dsmcN_, "dsmcN", "parcelCount");
                writeSelected(dsmcNMean_, "dsmcNMean", "parcelCount");
                writeSelected(rhoN_, "rhoN", "numberDensity");
                writeSelected(rhoM_, "rhoM", "massDensity");

                writeSelected(p_, "p", "pressure");
                writeSelected(Ttra_, "Ttra", "temperature");
                writeSelected(UMean_, "U", "velocity");
                writeSelected(UMagMean_, "UMag", "velocity");
                writeSelected(UAbsMean_, "UAbs", "absoluteVelocity");
                writeSelected(UAbsMagMean_, "UAbsMag", "absoluteVelocity");
                writeSelected(UThetaRelMean_, "UThetaRel", "thetaVelocity");
                writeSelected(UThetaAbsMean_, "UThetaAbs", "thetaVelocity");
                writeSelected(rotatingRegionMean_, "rotatingRegion", "srfDiagnostics");
                writeSelected(aCoriolisMean_, "aCoriolis", "srfDiagnostics");
                writeSelected(aCentrifugalMean_, "aCentrifugal", "srfDiagnostics");
                writeSelected(aSRFMean_, "aSRF", "srfDiagnostics");
                writeSelected(Ma_, "Ma", "mach");
                if
                (
                    outputProperty("wallHeatFlux", "wallDiagnostics")
                )
                {
                    q_.write();
                }
                writeSelected(fD_, "fD", "wallDiagnostics");
                if
                (
                    outputProperty("wallShearStress", "wallDiagnostics")
                )
                {
                    tau_.write();
                }

                if (writeRotationalTemperature_)
                {
                    writeSelected(Trot_, "Trot", "temperature");
                }
                if (writeVibrationalTemperature_)
                {
                    writeSelected(Tvib_, "Tvib", "temperature");
                }
                if (writeElectronicTemperature_)
                {
                    writeSelected(Telec_, "Telec", "temperature");
                }
                if
                (
                      writeRotationalTemperature_
                   or writeVibrationalTemperature_
                   or writeElectronicTemperature_
                )
                {
                    writeSelected(Tov_, "Tov", "temperature");
                }

                if (measureMeanFreePath_)
                {
                    writeSelected(mfp_, "mfp", "meanFreePath");
                    writeSelected(mfpToDx_, "mfpToDx", "meanFreePath");
                    writeSelected(DxToMfp_, "DxToMfp", "meanFreePath");
                    writeSelected
                    (
                        meanCollisionRate_,
                        "meanCollisionRate",
                        "collisionDiagnostics"
                    );
                    writeSelected
                    (
                        measuredCollisionRate_,
                        "measuredCollisionRate",
                        "collisionDiagnostics"
                    );
                    writeSelected
                    (
                        meanCollisionSeparation_,
                        "meanCollisionSeparation",
                        "collisionDiagnostics"
                    );
                    writeSelected(meanCollisionTime_, "mct", "meanFreePath");
                    writeSelected(mctToDt_, "mctToDt", "meanFreePath");
                    writeSelected(SOF_, "SOFP", "meanFreePath");
                }

                if (measureClassifications_)
                {
                    writeSelected
                    (
                        classIDistribution_,
                        "classIDistribution",
                        "classification"
                    );
                    writeSelected
                    (
                        classIIDistribution_,
                        "classIIDistribution",
                        "classification"
                    );
                    writeSelected
                    (
                        classIIIDistribution_,
                        "classIIIDistribution",
                        "classification"
                    );
                }

                if (measureErrors_)
                {
                    writeSelected(densityError_, "densityError", "errors");
                    writeSelected(velocityError_, "velocityError", "errors");
                    writeSelected
                    (
                        temperatureError_,
                        "temperatureError",
                        "errors"
                    );
                    writeSelected(pressureError_, "pressureError", "errors");
                }

                if (measureHeatFluxShearStress_)
                {
                    writeSelected
                    (
                        heatFluxVector_,
                        "heatFluxVector",
                        "heatFluxShearStress"
                    );
                    writeSelected
                    (
                        pressureTensor_,
                        "pressureTensor",
                        "heatFluxShearStress"
                    );
                    writeSelected
                    (
                        shearStressTensor_,
                        "shearStressTensor",
                        "heatFluxShearStress"
                    );
                }
            }
        }

        //- Reset fields after printing the instantaneous solution ... or
        //  continue sampling
        if (time_.resetFieldsAtOutput())
        {
            nTimeSteps_ = 0.0;

            forAll(dsmcNCum_, celli)
            {
                dsmcNCum_[celli] = 0.0;
                dsmcMCum_[celli] = 0.0;
                dsmcLinearKECum_[celli] = 0.0;
                dsmcMomentumCum_[celli] = vector::zero;
                dsmcErotCum_[celli] = 0.0;
                dsmcZetaRotCum_[celli] = 0.0;
                dsmcNElecLvlCum_[celli] = 0.0,
                dsmcNClassICum_[celli] = 0.0;
                dsmcNClassIICum_[celli] = 0.0;
                dsmcNClassIIICum_[celli] = 0.0;
                collisionSeparation_[celli] = 0.0;
                dsmcNCollsCum_[celli] = 0.0;
                measuredCollisionRate_[celli] = 0.0;
                dsmcMuuCum_[celli] = 0.0;
                dsmcMuvCum_[celli] = 0.0;
                dsmcMuwCum_[celli] = 0.0;
                dsmcMvvCum_[celli] = 0.0;
                dsmcMvwCum_[celli] = 0.0;
                dsmcMwwCum_[celli] = 0.0;
                dsmcMccCum_[celli] = 0.0;
                dsmcMccuCum_[celli] = 0.0;
                dsmcMccvCum_[celli] = 0.0;
                dsmcMccwCum_[celli] = 0.0;
                dsmcEuCum_[celli] = 0.0;
                dsmcEvCum_[celli] = 0.0;
                dsmcEwCum_[celli] = 0.0;
                dsmcECum_[celli] = 0.0;
                zetaVib_[celli] = 0.0;
                nCum_[celli] = 0.0;
                mCum_[celli] = 0.0;
                momentumCum_[celli] = vector::zero;
                linearKECum_[celli] = 0.0;
            }

            forAll(speciesIds_, i)
            {
                forAll(speciesTvib_[i], celli)
                {
                    dsmcNSpeciesCum_[i][celli] = 0.0;
                    nSpeciesCum_[i][celli] = 0.0;
                    dsmcMccSpeciesCum_[i][celli] = 0.0;
                    speciesMfp_[i][celli] = 0.0;
                    speciesMcr_[i][celli] = 0.0;
                    speciesTvib_[i][celli] = 0.0;
                    dsmcSpeciesEelecCum_[i][celli] = 0.0;
                    dsmcNGrndElecLvlSpeciesCum_[i][celli] = 0.0;
                    dsmcN1stElecLvlSpeciesCum_[i][celli] = 0.0;
                }

                forAll(dsmcSpeciesEvibModCum_[i], mod)
                {
                    forAll(dsmcSpeciesEvibModCum_[i][mod], celli)
                    {
                        dsmcSpeciesEvibModCum_[i][mod][celli] = 0.0;
                    }
                }
            }

            //- Reset boundary information
            forAll(rhoNBF_, j)
            {
                rhoNBF_[j] = 0.0;
                rhoMBF_[j] = 0.0;
                linearKEBF_[j] = 0.0;
                rhoNIntBF_[j] = 0.0;

                ErotBF_[j] = 0.0;
                zetaRotBF_[j] = 0.0;
                zetaVibBF_[j] = 0.0;
                rhoNElecBF_[j] = 0.0;

                qBF_[j] = 0.0;
                fDBF_[j] = vector::zero;
                momentumBF_[j] = vector::zero;
            }

            forAll(speciesIds_, i)
            {
                forAll(speciesTvibBF_[i], j)
                {
                    speciesRhoNBF_[i][j] = 0.0;
                    speciesMccBF_[i][j] = 0.0;
                    speciesEvibBF_[i][j] = 0.0;
                    speciesTvibBF_[i][j] = 0.0;
                    speciesZetaVibBF_[i][j] = 0.0;
                    speciesEelecBF_[i][j] = 0.0;
                }

                forAll(speciesEvibModBF_[i], mod)
                {
                    forAll(speciesEvibModBF_[i][mod], j)
                    {
                        speciesEvibModBF_[i][mod][j] = 0.0;
                    }
                }
            }
        }

        if (averagingAcrossManyRuns_ && !time_.resetFieldsAtOutput())
        {
            writeOut();
        }
    }
}

//- reset fields when mesh is edited
void dsmcVolFields::resetField()
{
    const label nCells = mesh_.nCells();

    nTimeSteps_ = 0.0;

    //- Reset volume information
    dsmcNCum_.clear();
    dsmcMCum_.clear();
    dsmcLinearKECum_.clear();
    dsmcMomentumCum_.clear();
    dsmcErotCum_.clear();
    dsmcZetaRotCum_.clear();
    dsmcNElecLvlCum_.clear();
    dsmcNClassICum_.clear();
    dsmcNClassIICum_.clear();
    dsmcNClassIIICum_.clear();
    collisionSeparation_.clear();
    dsmcNCollsCum_.clear();
    dsmcMuuCum_.clear();
    dsmcMuvCum_.clear();
    dsmcMuwCum_.clear();
    dsmcMvvCum_.clear();
    dsmcMvwCum_.clear();
    dsmcMwwCum_.clear();
    dsmcMccCum_.clear();
    dsmcMccuCum_.clear();
    dsmcMccvCum_.clear();
    dsmcMccwCum_.clear();
    dsmcEuCum_.clear();
    dsmcEvCum_.clear();
    dsmcEwCum_.clear();
    dsmcECum_.clear();
    zetaVib_.clear();
    nCum_.clear();
    mCum_.clear();
    momentumCum_.clear();
    linearKECum_.clear();

    dsmcNCum_.setSize(nCells, 0.0);
    dsmcMCum_.setSize(nCells, 0.0);
    dsmcLinearKECum_.setSize(nCells, 0.0);
    dsmcMomentumCum_.setSize(nCells, vector::zero);
    dsmcErotCum_.setSize(nCells, 0.0);
    dsmcZetaRotCum_.setSize(nCells, 0.0);
    dsmcNElecLvlCum_.setSize(nCells, 0.0);
    dsmcNClassICum_.setSize(nCells, 0.0);
    dsmcNClassIICum_.setSize(nCells, 0.0);
    dsmcNClassIIICum_.setSize(nCells, 0.0);
    collisionSeparation_.setSize(nCells, 0.0);
    dsmcNCollsCum_.setSize(nCells, 0.0);
    measuredCollisionRate_.setSize(nCells, 0.0);
    dsmcMuuCum_.setSize(nCells, 0.0);
    dsmcMuvCum_.setSize(nCells, 0.0);
    dsmcMuwCum_.setSize(nCells, 0.0);
    dsmcMvvCum_.setSize(nCells, 0.0);
    dsmcMvwCum_.setSize(nCells, 0.0);
    dsmcMwwCum_.setSize(nCells, 0.0);
    dsmcMccCum_.setSize(nCells, 0.0);
    dsmcMccuCum_.setSize(nCells, 0.0);
    dsmcMccvCum_.setSize(nCells, 0.0);
    dsmcMccwCum_.setSize(nCells, 0.0);
    dsmcEuCum_.setSize(nCells, 0.0);
    dsmcEvCum_.setSize(nCells, 0.0);
    dsmcEwCum_.setSize(nCells, 0.0);
    dsmcECum_.setSize(nCells, 0.0);
    zetaVib_.setSize(nCells, 0.0);
    nCum_.setSize(nCells, 0.0);
    mCum_.setSize(nCells, 0.0);
    momentumCum_.setSize(nCells, vector::zero);
    linearKECum_.setSize(nCells, 0.0);

    forAll(speciesIds_, i)
    {
        dsmcNSpeciesCum_[i].clear();
        nSpeciesCum_[i].clear();
        dsmcMccSpeciesCum_[i].clear();
        speciesMfp_[i].clear();
        speciesMcr_[i].clear();
        speciesTvib_[i].clear();
        dsmcSpeciesEelecCum_[i].clear();
        dsmcNGrndElecLvlSpeciesCum_[i].clear();
        dsmcN1stElecLvlSpeciesCum_[i].clear();

        dsmcNSpeciesCum_[i].setSize(nCells, 0.0);
        nSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcMccSpeciesCum_[i].setSize(nCells, 0.0);
        speciesMfp_[i].setSize(nCells, 0.0);
        speciesMcr_[i].setSize(nCells, 0.0);
        speciesTvib_[i].setSize(nCells, 0.0);
        dsmcSpeciesEelecCum_[i].setSize(nCells, 0.0);
        dsmcNGrndElecLvlSpeciesCum_[i].setSize(nCells, 0.0);
        dsmcN1stElecLvlSpeciesCum_[i].setSize(nCells, 0.0);

        forAll(dsmcSpeciesEvibModCum_[i], mod)
        {
           dsmcSpeciesEvibModCum_[i][mod].clear();
           dsmcSpeciesEvibModCum_[i][mod].setSize(nCells, 0.0);
        }
    }

    //- Reset boundary information
    forAll(mesh_.boundaryMesh(), j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        rhoNBF_[j].clear();
        rhoMBF_[j].clear();
        linearKEBF_[j].clear();
        momentumBF_[j].clear();
        ErotBF_[j].clear();
        zetaRotBF_[j].clear();
        qBF_[j].clear();
        fDBF_[j].clear();
        zetaVibBF_[j].clear();
        rhoNIntBF_[j].clear();
        rhoNElecBF_[j].clear();

        n_[j].clear();
        t1_[j].clear();
        t2_[j].clear();

        rhoNBF_[j].setSize(nFaces, 0.0);
        rhoMBF_[j].setSize(nFaces, 0.0);
        linearKEBF_[j].setSize(nFaces, 0.0);
        momentumBF_[j].setSize(nFaces, vector::zero);
        ErotBF_[j].setSize(nFaces, 0.0);
        zetaRotBF_[j].setSize(nFaces, 0.0);
        qBF_[j].setSize(nFaces, 0.0);
        fDBF_[j].setSize(nFaces, vector::zero);
        zetaVibBF_[j].setSize(nFaces, 0.0);
        rhoNIntBF_[j].setSize(nFaces, 0.0);
        rhoNElecBF_[j].setSize(nFaces, 0.0);

        n_[j].setSize(nFaces, vector::zero);
        t1_[j].setSize(nFaces, vector::zero);
        t2_[j].setSize(nFaces, vector::zero);
    }

    forAll(speciesIds_, i)
    {
        const label nPatches = mesh_.boundaryMesh().size();

        speciesEvibBF_[i].clear();
        speciesEelecBF_[i].clear();
        speciesRhoNBF_[i].clear();
        speciesMccBF_[i].clear();
        speciesTvibBF_[i].clear();
        speciesZetaVibBF_[i].clear();
        speciesEvibModBF_[i].clear();

        speciesRhoNBF_[i].setSize(nPatches);
        speciesMccBF_[i].setSize(nPatches);
        speciesEvibBF_[i].setSize(nPatches);
        speciesTvibBF_[i].setSize(nPatches);
        speciesZetaVibBF_[i].setSize(nPatches);
        speciesEelecBF_[i].setSize(nPatches);

        forAll(mesh_.boundaryMesh(), j)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[j];
            const label nFaces = patch.size();

            speciesRhoNBF_[i][j].clear();
            speciesMccBF_[i][j].clear();
            speciesEvibBF_[i][j].clear();
            speciesTvibBF_[i][j].clear();
            speciesZetaVibBF_[i][j].clear();
            speciesEelecBF_[i][j].clear();

            speciesRhoNBF_[i][j].setSize(nFaces, 0.0);
            speciesMccBF_[i][j].setSize(nFaces, 0.0);
            speciesEvibBF_[i][j].setSize(nFaces, 0.0);
            speciesTvibBF_[i][j].setSize(nFaces, 0.0);
            speciesZetaVibBF_[i][j].setSize(nFaces, 0.0);
            speciesEelecBF_[i][j].setSize(nFaces, 0.0);
        }

        forAll(speciesEvibModBF_[i], mod)
        {
            const label nPatches = mesh_.boundaryMesh().size();
            speciesEvibModBF_[i][mod].setSize(nPatches);

            forAll(speciesEvibModBF_[i][mod], j)
            {
                const polyPatch& patch = mesh_.boundaryMesh()[j];
                const label nFaces = patch.size();
                speciesEvibModBF_[i][mod][j].setSize(nFaces, 0.0);
            }
        }
    }

    forAll(boundaryCells_, j)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[j];
        const label nFaces = patch.size();

        boundaryCells_[j].clear();
        boundaryCells_[j].setSize(nFaces);

        forAll(boundaryCells_[j], k)
        {
            boundaryCells_[j][k] = patch.faceCells()[k];
        }
    }
}

void dsmcVolFields::writeField()
{}

void dsmcVolFields::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBasicFieldProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");
    configureOutputProperties();
}

} // End namespace Foam

// ************************************************************************** //
