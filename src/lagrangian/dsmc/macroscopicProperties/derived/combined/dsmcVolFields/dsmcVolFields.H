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

Class
    dsmcVolFields

Description
    Geometric field for temperature, density, velocity, pressure.

SourceFiles
    dsmcVolFields.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcVolFields_H
#define dsmcVolFields_H

#include "dsmcField.H"
#include "vector.H"
#include "polyMesh.H"
#include "volFields.H"
#include "dsmcCloud.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dsmcVolFields Declaration
\*---------------------------------------------------------------------------*/

class dsmcVolFields
:
	public dsmcField
{

protected:

    // Protected data

        dictionary propsDict_;

        label sampleInterval_;
        label sampleCounter_;
        scalar nTimeSteps_;
        
        scalar mfpTref_;
        word fieldName_;
        
        labelList speciesIds_;
        
        volScalarField dsmcN_;
        volScalarField dsmcNMean_;
        volScalarField rhoN_;
        volScalarField rhoM_;
        volScalarField p_;
        volScalarField Ttra_;
        volScalarField Trot_;
        volScalarField Tvib_;
        volScalarField Telec_;
        volScalarField Tov_;
        volScalarField q_;
        volScalarField tau_;
        volScalarField mfp_;
        volScalarField mfpToDx_;
        volScalarField DxToMfp_;
        volScalarField meanCollisionRate_;
        volScalarField measuredCollisionRate_;
        volScalarField meanCollisionTime_;
        volScalarField mctToDt_;
        volScalarField meanCollisionSeparation_;
        volScalarField SOF_;
        volScalarField Ma_;
        volScalarField classIDistribution_;
        volScalarField classIIDistribution_;
        volScalarField classIIIDistribution_;
        volScalarField densityError_;
        volScalarField velocityError_;
        volScalarField temperatureError_;
        volScalarField pressureError_;
        volVectorField UMean_;
        //- UMean_ 的速度大小。
        //  注意：在旋转区 UMean_ 表示相对速度 U_rel，在非旋转区表示绝对速度；
        //  因此该场用于查看“当前求解变量自身”的速度模。
        volScalarField UMagMean_;
        volVectorField UAbsMean_;
        //- 绝对速度 UAbsMean_ 的速度大小。
        //  跨 rotor/stator 统一比较物理流速时，应优先查看该场。
        volScalarField UAbsMagMean_;
        //- 相对参考系下的周向速度分量。
        //  对 rotor 区域，它对应 parcel 在 SRF 方程中实际推进所使用的
        //  相对速度 U_rel 在柱坐标 e_theta 方向上的投影：
        //      u_theta,rel = U_rel · e_theta
        //  对非旋转区域（例如 statorZone），当前实现会显式写成 0。
        //  这样做不是说 stator 没有周向绝对速度，而是为了避免把
        //  “本来没有相对参考系定义的区域”误读成也存在一个可比较的
        //  u_theta,rel 场。跨 domain 的统一比较应改用 UAbs/UThetaAbs。
        //  该量专门用于区分：
        //  1. “绝对速度图里看不到明显旋转，但相对系里其实已建立周向趋势”
        //  2. “相对系里也没有建立周向趋势，说明问题不只是后处理展示方式”
        volScalarField UThetaRelMean_;
        volScalarField UThetaAbsMean_;
        //- 旋转区掩码。用于把“哪里允许施加 SRF 力”直接输出成一个
        //  可视化体场：
        //      rotatingRegion = 1  -> 当前 cell 属于 rotor 的 SRF 区
        //      rotatingRegion = 0  -> 当前 cell 属于 stator / inertial 区
        //  这个字段本身不依赖采样统计，因此是验证 cellZone 限域是否正确的
        //  最直接几何判据。
        volScalarField rotatingRegionMean_;
        //- 科里奥利加速度诊断场：
        //      a_cor = -2 * (Omega x U_rel)
        //  这里只在旋转区输出，非旋转区显式写零。
        volVectorField aCoriolisMean_;
        //- 离心加速度诊断场：
        //      a_cen = -(Omega x (Omega x r))
        //  该项只依赖几何位置与 Omega，本身就是判断 rotor 区是否真的挂上
        //  旋转参考系最干净的空间分布量之一。
        volVectorField aCentrifugalMean_;
        //- SRF 总表观加速度诊断场：
        //      a_SRF = a_cor + a_cen
        //  这里与 dsmcParcel::move() 中真正使用的 SRF 力项保持同一公式，
        //  用于直接回答“多流域 rotor 内是否真的施加了旋转力”。
        volVectorField aSRFMean_;
        volVectorField fD_;
        volVectorField heatFluxVector_;
        volTensorField pressureTensor_;
        volTensorField shearStressTensor_;

        scalarField dsmcNCum_;
        scalarField nCum_;
        scalarField dsmcNElecLvlCum_;
        scalarField dsmcMCum_;
        scalarField mCum_;
        scalarField dsmcLinearKECum_;
        scalarField linearKECum_;
        scalarField dsmcErotCum_;
        scalarField dsmcZetaRotCum_;
        scalarField dsmcMuuCum_;
        scalarField dsmcMuvCum_;
        scalarField dsmcMuwCum_;
        scalarField dsmcMvvCum_;
        scalarField dsmcMvwCum_;
        scalarField dsmcMwwCum_;
        scalarField dsmcMccCum_;
        scalarField dsmcMccuCum_;
        scalarField dsmcMccvCum_;
        scalarField dsmcMccwCum_;
        scalarField dsmcEuCum_;
        scalarField dsmcEvCum_;
        scalarField dsmcEwCum_;
        scalarField dsmcECum_;
        scalarField zetaVib_;
        scalarField dsmcNClassICum_;
        scalarField dsmcNClassIICum_;
        scalarField dsmcNClassIIICum_;
        scalarField collisionSeparation_;
        scalarField dsmcNCollsCum_;
        vectorField dsmcMomentumCum_;
        vectorField momentumCum_;
        
        List<labelList> boundaryCells_;
        List<List<scalarField>> dsmcSpeciesEvibModCum_;
        List<scalarField> dsmcSpeciesEelecCum_;
        List<scalarField> dsmcNSpeciesCum_;
        List<scalarField> nSpeciesCum_;
        List<scalarField> dsmcMccSpeciesCum_;
        List<scalarField> speciesTvib_;
        List<scalarField> dsmcNGrndElecLvlSpeciesCum_;
        List<scalarField> dsmcN1stElecLvlSpeciesCum_;
        List<scalarField> speciesMfp_;
        List<scalarField> speciesMcr_;


      //- Boundary measurements
      
        List<scalarField> rhoNBF_;
        List<scalarField> rhoMBF_;
        List<scalarField> linearKEBF_;
        List<scalarField> ErotBF_;
        List<scalarField> zetaRotBF_;
        List<scalarField> qBF_;
        List<scalarField> zetaVibBF_;
        List<scalarField> rhoNIntBF_;
        List<scalarField> rhoNElecBF_;
        List<vectorField> momentumBF_;
        List<vectorField> fDBF_;

        List<List<scalarField>> speciesEvibBF_;
        List<List<scalarField>> speciesEelecBF_;
        List<List<scalarField>> speciesRhoNBF_;
        List<List<scalarField>> speciesMccBF_;
        List<List<scalarField>> speciesTvibBF_;
        List<List<scalarField>> speciesZetaVibBF_;
        List<List<List<scalarField>>> speciesEvibModBF_;

        //- Wall unit vectors
        //  temporary storage, should be done in boundaries/patch
        List<vectorField> n_;
        List<vectorField> t1_;
        List<vectorField> t2_;

        bool averagingAcrossManyRuns_;
        bool measureClassifications_;
        bool measureMeanFreePath_;
        bool measureErrors_;
        bool densityOnly_;
        bool measureHeatFluxShearStress_;
        bool writeRotationalTemperature_;
        bool writeVibrationalTemperature_;
        bool writeElectronicTemperature_;
        wordHashSet outputProperties_;

    // Protected Member Functions

        void readIn();

        void writeOut();

        void calculateWallUnitVectors();

        void configureOutputProperties();

        void setAllFieldWriteOptions(const IOobject::writeOption writeOpt);

        bool outputProperty(const word& propertyName) const;

        bool outputProperty
        (
            const word& propertyName,
            const word& propertyGroup
        ) const;

        template<class FieldType>
        void writeSelected
        (
            FieldType& field,
            const word& propertyName,
            const word& propertyGroup = word::null
        ) const
        {
            if (outputProperty(propertyName, propertyGroup))
            {
                field.write();
            }
        }

public:

    //- Runtime type information
    TypeName("dsmcVolFields");

    // Constructors

        //- Construct from components
        dsmcVolFields
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );

    // Destructor

    ~dsmcVolFields();


    // Member Functions

      //- create initial field
      virtual void createField();

      //- calculate field
      virtual void calculateField();

      //- reset field if mesh has changed
      virtual void resetField();

      //- write field
      virtual void writeField();

      virtual void updateProperties(const dictionary& dict);

    // Access

      virtual scalar translationalT(const label cellI)
      {
          return Ttra_[cellI];
      }

      virtual scalar overallT(const label cellI)
      {
          return Tov_[cellI];
      }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
