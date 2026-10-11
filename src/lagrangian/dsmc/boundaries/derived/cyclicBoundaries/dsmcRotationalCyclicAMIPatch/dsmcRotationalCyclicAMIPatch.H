/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef dsmcRotationalCyclicAMIPatch_H
#define dsmcRotationalCyclicAMIPatch_H

#include "dsmcCyclicBoundary.H"
#include "tensor.H"
#include "DynamicList.H"
#include "Map.H"

namespace Foam
{

class dsmcRotationalCyclicAMIPatch
:
    public dsmcCyclicBoundary
{
private:

    dictionary propsDict_;
    vector rotationAxis_;
    vector rotationCentre_;
    scalar rotationAngle_;
    bool transformPosition_;
    bool debug_;
    label maxDebugPrints_;
    label debugPrintCount_;

    //- 从字典读取旋转周期 AMI 的几何与调试参数。
    // 这里把 rotationAxis / rotationCentre / rotationAngle / debug 等
    // 集中读取，避免构造函数、updateProperties()、运行期重读字典时
    // 各自维护一套分散逻辑。
    void readProperties();

    //- 使用 Rodrigues 旋转公式旋转一个“纯向量”。
    // 这里的 value 可以是速度、面法向面积向量、以及任何不带平移项的几何量。
    // 对应公式为：
    //     v' = v cos(theta)
    //        + (a_hat x v) sin(theta)
    //        + a_hat (a_hat · v) (1 - cos(theta))
    // 其中 a_hat 是已经单位化后的 rotationAxis_。
    vector rotateVectorRodrigues
    (
        const vector& value,
        const scalar signedAngle
    ) const;

    //- 使用 Rodrigues 旋转公式旋转一个“绝对位置点”。
    // 对位置点必须先减去旋转中心，再对相对向量做 Rodrigues 旋转，
    // 最后再把旋转中心加回去。对应公式为：
    //     p' = c + R(p - c)
    // 这里故意把“位置点旋转”和“纯向量旋转”拆成两个函数，避免后续维护时
    // 把带平移项的位置变换与不带平移项的速度/法向变换混用。
    point rotatePointRodrigues
    (
        const point& position,
        const scalar signedAngle
    ) const;

    //- 只把非退化且方向不重复的候选搜索方向压入列表。
    // 对曲面非严格匹配的 AMI 接口而言，单一搜索方向并不总是可靠，
    // 因此这里需要维护一组有序回退方向；但若多个方向几乎共线，
    // 再重复尝试只会制造噪声而不会提高成功率。
    void appendSearchDirection
    (
        DynamicList<vector>& directions,
        DynamicList<word>& directionNames,
        const vector& candidate,
        const word& candidateName
    ) const;

    //- 记录完整跨面事件统计。
    // 这些量的目标不是替代宏观量场，而是回答：
    // 1. A->B 与 B->A 的事件数是否平衡；
    // 2. 哪些目标 face 被命中过多；
    // 3. 法向速度加权后的“事件通量”是否存在系统偏置。
    //
    // 这里把统计放在边界模型内部，而不是继续依赖前 2000 条 debug 行，
    // 是为了让并行运行时也能拿到完整 run 的累计统计，而不是一个被
    // maxDebugPrints 截断的早期样本。
    label totalTransfersAtoB_;
    label totalTransfersBtoA_;
    scalar totalAbsUnAtoB_;
    scalar totalAbsUnBtoA_;

    label totalTransfersAtoBStrict_;
    label totalTransfersAtoBLoose_;
    label totalTransfersBtoAStrict_;
    label totalTransfersBtoALoose_;
    scalar totalAbsUnAtoBStrict_;
    scalar totalAbsUnAtoBLoose_;
    scalar totalAbsUnBtoAStrict_;
    scalar totalAbsUnBtoALoose_;

    //- 上一次 output() 时刻已经写出的累计值。
    // 有了这组“上次写出快照”，就能在每个输出时刻同时给出：
    // 1. 从计算开始到当前时刻的累计事件数/累计通量；
    // 2. 当前写出区间内新增的事件数/新增通量。
    label writtenTransfersAtoB_;
    label writtenTransfersBtoA_;
    scalar writtenAbsUnAtoB_;
    scalar writtenAbsUnBtoA_;

    //- 分模式写出快照。
    // 与总累计值一样，这些量用于给每个 output() 时刻同时输出：
    // 1. strict / loose 的累计贡献；
    // 2. 当前写出区间内 strict / loose 的新增贡献。
    label writtenTransfersAtoBStrict_;
    label writtenTransfersAtoBLoose_;
    label writtenTransfersBtoAStrict_;
    label writtenTransfersBtoALoose_;
    scalar writtenAbsUnAtoBStrict_;
    scalar writtenAbsUnAtoBLoose_;
    scalar writtenAbsUnBtoAStrict_;
    scalar writtenAbsUnBtoALoose_;

    //- 目标 face 命中计数。
    // key 是迁移后粒子所在的目标 face label，value 是累计命中次数。
    // 分开保存 A->B 和 B->A，便于后续检查是否出现目标面高度集中。
    Map<label> targetFaceHitsAtoB_;
    Map<label> targetFaceHitsBtoA_;

    //- 分模式目标 face 命中计数。
    // 若下一轮统计显示剩余偏差主要由 loose 路径贡献，则这些表会直接告诉我们
    // 到底是哪些目标面在 loose 路径下被反复命中。
    Map<label> targetFaceHitsAtoBStrict_;
    Map<label> targetFaceHitsAtoBLoose_;
    Map<label> targetFaceHitsBtoAStrict_;
    Map<label> targetFaceHitsBtoALoose_;

    void recordCrossingEvent
    (
        const word& directionLabel,
        const word& searchModeLabel,
        const label targetFace,
        const scalar absUn
    );

    //- 把 face 命中统计按 face label 排序写出到 CSV。
    // 这里故意输出完整列表而不是只写 top-N，避免后续做空间分桶时
    // 又需要回头改代码重新导出。
    void writeFaceHitTable
    (
        const fileName& outputFile,
        const Map<label>& faceHits,
        const word& directionLabel
    ) const;

public:

    TypeName("dsmcRotationalCyclicAMIPatch");

    dsmcRotationalCyclicAMIPatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcRotationalCyclicAMIPatch();

    void calculateProperties();
    void initialConfiguration();

    void controlMol
    (
        dsmcParcel& p,
        dsmcParcel::trackingData& td
    );

    void output
    (
        const fileName& fixedPathName,
        const fileName& timePath
    );

    void updateProperties(const dictionary&);
};

} // End namespace Foam

#endif

// ************************************************************************* //
