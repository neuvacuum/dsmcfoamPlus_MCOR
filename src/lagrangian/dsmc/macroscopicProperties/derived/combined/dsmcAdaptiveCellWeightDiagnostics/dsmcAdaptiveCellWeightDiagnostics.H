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
    dsmcAdaptiveCellWeightDiagnostics

Description
    Fixed-mesh parcel occupancy and represented-molecule diagnostics.
    Accumulates weighted fields, computes target cell weights and writes
    global conservation summaries.

SourceFiles
    dsmcAdaptiveCellWeightDiagnostics.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcAdaptiveCellWeightDiagnostics_H
#define dsmcAdaptiveCellWeightDiagnostics_H

#include "dsmcField.H"
#include "volFields.H"

namespace Foam
{

class dsmcAdaptiveCellWeightDiagnostics
:
    public dsmcField
{
protected:

    // 诊断字段自身的配置子字典。
    // 读取顺序是：
    // 1. 代码内硬编码默认值；
    // 2. dsmcProperties/adaptiveProperties 的全局默认值；
    // 3. fieldPropertiesDict 当前 field 条目的局部覆盖值。
    // 这样既能保证算例级统一默认配置，又允许单个 field 局部细调。
    dictionary propsDict_;

    // 采样步长控制。
    // 只有当 sampleCounter_ 累计到 sampleInterval_ 时，才把当前 parcel
    // 状态纳入一次有效统计样本；这样可以避免每个 DSMC 微步都把同一批
    // 高频波动直接写进目标权重场。
    label sampleInterval_;
    label sampleCounter_;

    // 已累计的有效样本数。
    // 所有 cell 级输出场最终都按这个样本数做时间平均。
    scalar nTimeSteps_;

    // 当前诊断字段的名字后缀，例如 mixture / H2 / rotor。
    word fieldName_;

    // 参与本次诊断的 species 列表。
    labelList speciesIds_;

    scalar targetParticlesPerSubCell_;

    scalar assumedSubCellsPerCell_;

    // 是否在每个输出时刻写出一份独立的全局守恒摘要字典。
    bool writeDiagnosticDictionary_;

    // 每个 cell 的代表粒子数时间平均。
    volScalarField particlesPerCell_;

    // 每个 cell 的真实代表分子数时间平均。
    volScalarField realMoleculesPerCell_;

    volScalarField cellWeightFactor_;

    // 按当前固定统计权重语义计算得到的加权数密度。
    volScalarField weightedRhoN_;

    // 按当前固定统计权重语义计算得到的加权质量密度。
    volScalarField weightedRhoM_;

    // 按当前固定统计权重语义计算得到的加权动量密度。
    volVectorField weightedMomentum_;

    // 以下累计量用于在一个输出窗口内做时间平均。
    scalarField particlesPerCellCum_;
    scalarField realMoleculesPerCellCum_;
    scalarField weightedMassCum_;
    vectorField weightedMomentumCum_;

    // 从 fieldPropertiesDict 与 dsmcProperties/adaptiveProperties 读取并合并
    // 当前诊断模块所需的参数。
    void updateConfiguration();

    // 把字符串 species 名解析成 cloud 内部 typeId。
    void initialiseSpeciesIds();

    // 清空一个输出窗口内的累计量。
    void resetAccumulators();

    // 判断一个 parcel 是否属于当前诊断对象。
    bool includeParcel(const dsmcParcel& p) const;

    // 计算目标 cellWeightFactor 时采用的目标代表粒子数。
    scalar targetRepresentativeParcelsPerCell() const;

    // 在每个输出时刻写出全局守恒摘要。
    void writeConservationSummary() const;


public:

    TypeName("dsmcAdaptiveCellWeightDiagnostics");

    dsmcAdaptiveCellWeightDiagnostics
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcAdaptiveCellWeightDiagnostics();

    virtual void createField();

    virtual void calculateField();

    virtual void resetField();

    virtual void writeField();

    virtual void updateProperties(const dictionary& newDict);

    virtual scalar translationalT(const label cellI);

    virtual scalar overallT(const label cellI);
};

} // End namespace Foam

#endif

// ************************************************************************* //
