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
    dsmcCyclicBoundary

Description

    Used for applying instantaneous actions on molecules colliding with a cyclic
    boundary.

SourceFiles
    dsmcCyclicBoundaryI.H
    dsmcCyclicBoundary.C
    dsmcCyclicBoundaryIO.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCyclicBoundary_H
#define dsmcCyclicBoundary_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "vector.H"
#include "volFields.H"
#include "polyMesh.H"
#include "dsmcParcel.H"
#include "Time.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class dsmcCyclicBoundary Declaration
\*---------------------------------------------------------------------------*/

class dsmcCyclicBoundary
{

protected:

    // Protected data

        const fvMesh& mesh_;

        dsmcCloud& cloud_;

        Time& time_;

//         Random& rndGen_;

        //- subDictionary containing the properties
        dictionary boundaryDict_;

        //- control on cyclic patch
        word patchName_;
        label patchId_;

        //- corresponding/neighbouring coupled patch
        word neighbPatchName_;
        label neighbPatchId_;


        //- faces on the dsmcCyclicBoundary
        labelList faces_;

        labelList coupledFacesA_;
        labelList coupledFacesB_;

        //- cells corresponding to the dsmcCyclicBoundary faces
//         labelList cells_;

        labelList cellsA_;
        labelList cellsB_;

        //- set all the properties below from model if required
        scalarField densities_;
        vectorField velocities_;
        scalarField temperatures_;
        tensor strainRate_;

        //- writing out properties
        bool writeInTimeDir_;
        bool writeInCase_;

        //- 并行 AMI 局部目标面搜索诊断开关。
        // 仅在需要排查并行目标面评分与候选选择时打开；
        // 默认保持关闭，不影响正常跨面路径。
        bool parallelSearchDebug_;

        //- 并行局部 AMI 搜索诊断的最大打印次数。
        label maxParallelSearchDebugPrints_;

        //- 已打印的并行局部搜索诊断计数。
        // tryStandardCyclicAMITransfer() 为 const，因此这里保持 mutable。
        mutable label parallelSearchDebugPrintCount_;


    // Private Member Functions

        //- 统一解析并校验 cyclic/cyclicAMI 成对 patch，再刷新控制面。
        //  构造、运行时字典重读和网格更新都必须走同一路径，避免
        //  patchName/neighbourPatchName 与缓存的 patchId/face 列表漂移。
        void refreshCoupledPatchState();

        void getCoupledFaces(const polyMesh& mesh);

        //- 统一执行 OpenFOAM 标准 cyclic 拓扑迁移。
        // 这个辅助函数的职责是把粒子从当前命中的 cyclic 面
        // 合法送到耦合侧，并同步更新后续 trackToFace() 继续需要的
        // face/cell/tetFace/tetPt 与几何位置。
        // 之所以把这段逻辑收敛到基类，是为了避免各个派生 cyclic
        // 边界各自复制一份实现后逐渐漂移。
        void applyStandardCyclicTransfer
        (
            dsmcParcel& p,
            const bool transformPosition
        ) const;

        bool tryStandardCyclicAMITransfer
        (
            dsmcParcel& p,
            const vector& direction,
            word* searchModeLabel = NULL,
            const point* targetReferencePosition = NULL,
            const point* targetInsertionReferencePosition = NULL,
            const vector* rotationAxis = NULL,
            const point* rotationCentre = NULL
        ) const;

        //- 统一执行 OpenFOAM 标准 cyclicAMI 拓扑迁移。
        // 这是给已经确定好搜索方向、且希望“失败即报错”的调用点使用的包装接口。
        // 当前内部实现会复用 tryStandardCyclicAMITransfer()，若其返回 false，
        // 再统一抛出 FatalError。
        // searchModeLabel 会在成功时返回：
        // - "strict"：通过严格几何一致性判据命中的目标面；
        // - "loose" ：strict 全空后退回到保底最近面；
        // - "native"：串行路径仍由原生 pointFace() 直接搜索得到。
        void applyStandardCyclicAMITransfer
        (
            dsmcParcel& p,
            const vector& direction,
            word* searchModeLabel = NULL,
            const point* targetReferencePosition = NULL,
            const point* targetInsertionReferencePosition = NULL,
            const vector* rotationAxis = NULL,
            const point* rotationCentre = NULL
        ) const;

public:

    //- Runtime type information
    TypeName("dsmcCyclicBoundary");

    // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            dsmcCyclicBoundary,
            dictionary,
            (
                Time& t,
                const polyMesh& mesh,
                dsmcCloud& cloud,
                const dictionary& dict
            ),
            (t, mesh, cloud, dict)
        );

    // Constructors

        //- Construct from components
        dsmcCyclicBoundary
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<dsmcCyclicBoundary> New
        (
            Time& t,
            const polyMesh& mesh,
            dsmcCloud& cloud,
            const dictionary& dict
        );


    // Destructor

        virtual ~dsmcCyclicBoundary();


    // Member Functions

        virtual void calculateProperties() = 0;

        //- create an initial configuration
        virtual void initialConfiguration() = 0;

        //- control the molecule
        virtual void controlMol
        (
            dsmcParcel& p,
            dsmcParcel::trackingData& td
        ) = 0;

        //- 在 Cloud::move() 的并行 transfer while-loop 内，
        //  处理该 cyclic boundary model 挂起的延迟并行粒子转移。
        //
        // 默认返回 false，表示当前模型没有任何额外远端注入/收发工作。
        // 只有像 dsmcStageInterfaceAMIPatch 这类需要任意-rank parcel 注入的
        // 模型，才会覆盖这个接口。
        virtual bool controlAfterMove()
        {
            return false;
        }

        //- output data
        virtual void output
        (
            const fileName& fixedPathName,
            const fileName& timePath
        ) = 0;

        //- update properties from a modified dictionary
        virtual void updateProperties(const dictionary&) = 0;

        void updateBoundaryProperties
        (
            const dictionary& newDict
        );

        //- scalar field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const scalarField& yData
        );

        void setNewBoundaryFields();

        //- vector field
        void writeTimeData
        (
            const fileName& pathName,
            const word& nameFile,
            const scalarField& xData,
            const vectorField& yData
        );

        // Access

            //- retrun faces on patch
            const labelList& controlPatch() const;


            //- return the adjacent cells
            const labelList& controlZone() const;

            const labelList& allFaces() const;

            //- return the name of the control patch
            const word& patchName() const;
            const label& patchId() const;

            //- return the targeted fields
            const scalarField& densityField() const;
            scalarField& densityField();

            const vectorField& velocityField() const;
            vectorField& velocityField();

            const scalarField& temperatureField() const;
            scalarField& temperatureField();

            const tensor& strainRate() const;
            tensor& strainRate();

            const bool& writeInTimeDir() const;
            const bool& writeInCase() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
