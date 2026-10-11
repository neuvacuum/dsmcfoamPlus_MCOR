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

#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::dsmcCloud::cloudName() const
{
    return cloudName_;
}

inline const Foam::dynamicFvMesh& Foam::dsmcCloud::mesh() const
{
    return mesh_;
}

inline const Foam::IOdictionary& Foam::dsmcCloud::particleProperties() const
{
    return particleProperties_;
}

inline const Foam::dictionary& Foam::dsmcCloud::controlDict() const
{
    return controlDict_;
}

inline const Foam::List<Foam::word>& Foam::dsmcCloud::typeIdList() const
{
    return typeIdList_;
}

inline Foam::scalar Foam::dsmcCloud::nParticle() const
{
    return coordSystem().nParticlesOrg();
}

inline const Foam::volScalarField& Foam::dsmcCloud::nParticles() const
{
    return coordSystem().nParticles();
}

inline Foam::scalar Foam::dsmcCloud::nParticles(const label celli) const
{
    return coordSystem().nParticles(celli);
}

inline Foam::scalar Foam::dsmcCloud::nParticles
(
    const label patchi,
    const label facei
) const
{
    return coordSystem().nParticles(patchi, facei);
}

inline Foam::scalar Foam::dsmcCloud::representedMolecules
(
    const dsmcParcel& parcel
) const
{
    // 这里显式把“时间步模型权重”和“parcel 级 RWF”相乘。
    // 这样：
    // 1. 对普通 Cartesian + 常时间步模型，结果退化为 nEquivalentParticles；
    // 2. 对 variableTimeStepModel，自动继承 cell 级 nParticles(celli)；
    // 3. 对 axisymmetric/spherical 或固定网格 split/merge，自然叠加
    //    parcel 自身的 RWF。
    //
    // 若 parcel 处于极少见的临时无效 cell 状态（例如某些异常调试路径），
    // 这里退回到全局 nParticle()*RWF，避免在诊断输出阶段因为越界而再次
    // 放大问题。
    const label celli = parcel.cell();

    if (celli >= 0 && celli < mesh_.nCells())
    {
        return coordSystem().dtModel().nParticles(celli)*parcel.RWF();
    }

    return nParticle()*parcel.RWF();
}

inline const Foam::scalarField& Foam::dsmcCloud::deltaTField() const
{
    return coordSystem().dtModel().deltaTField();
}

inline Foam::scalar Foam::dsmcCloud::deltaTValue(const label celli) const
{
    return coordSystem().dtModel().deltaTValue(celli);
}

inline bool Foam::dsmcCloud::rotatingReferenceFrameActive() const
{
    return rotatingReferenceFrameActive_;
}

inline bool Foam::dsmcCloud::isInRotatingRegion(const label celli) const
{
    return
    (
        rotatingReferenceFrameActive_
     && celli >= 0
     && celli < rotatingRegionCellMask_.size()
     && rotatingRegionCellMask_[celli]
    );
}

inline const Foam::vector& Foam::dsmcCloud::rotatingFrameOmega() const
{
    return rotatingFrameOmega_;
}

inline const Foam::vector& Foam::dsmcCloud::rotatingFrameOrigin() const
{
    return rotatingFrameOrigin_;
}

inline Foam::label Foam::dsmcCloud::rotatingFrameSubCycles() const
{
    return rotatingFrameSubCycles_;
}

inline Foam::scalar Foam::dsmcCloud::rotatingFrameMaxSubStep() const
{
    return rotatingFrameMaxSubStep_;
}

inline Foam::label Foam::dsmcCloud::nTerminalOutputs() const
{
    return nTerminalOutputs_;
}

inline const Foam::DynamicList<Foam::DynamicList<Foam::dsmcParcel*>>&
Foam::dsmcCloud::cellOccupancy() const
{
    return cellOccupancy_;
}

inline Foam::volScalarField& Foam::dsmcCloud::sigmaTcRMax()
{
    return sigmaTcRMax_;
}

inline Foam::scalarField& Foam::dsmcCloud::collisionSelectionRemainder()
{
    return collisionSelectionRemainder_;
}

inline const Foam::List<Foam::dsmcParcel::constantProperties>&
Foam::dsmcCloud::constProps() const
{
    return constProps_;
}

inline const Foam::dsmcParcel::constantProperties& Foam::dsmcCloud::constProps
(
    label typeId
) const
{
    if (typeId < 0 || typeId >= constProps_.size())
    {
        FatalErrorIn("Foam::dsmcCloud::constProps(label typeId)")
            << "constantProperties for requested typeId index "
            << typeId << " do not exist" << nl
            << abort(FatalError);
    }

    return constProps_[typeId];
}

inline Foam::Random& Foam::dsmcCloud::rndGen()
{
    return rndGen_;
}

inline const Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision() const
{
    return binaryCollisionModel_;
}

inline Foam::BinaryCollisionModel& Foam::dsmcCloud::binaryCollision()
{
    return binaryCollisionModel_();
}

inline const Foam::dsmcBoundaries& Foam::dsmcCloud::boundaries() const
{
    return  boundaries_;
}

inline Foam::dsmcBoundaries& Foam::dsmcCloud::boundaries()
{
    return  boundaries_;
}

inline const Foam::dsmcFieldProperties& Foam::dsmcCloud::fields() const
{
    return fields_;
}

inline Foam::dsmcFieldProperties& Foam::dsmcCloud::fields()
{
    return fields_;
}

inline const Foam::dsmcControllers& Foam::dsmcCloud::controllers() const
{
    return controllers_;
}

inline Foam::dsmcControllers& Foam::dsmcCloud::controllers()
{
    return controllers_;
}

inline const Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing() const
{
    return dynamicLoadBalancing_;
}

inline Foam::dsmcDynamicLoadBalancing& Foam::dsmcCloud::dynamicLoadBalancing()
{
    return dynamicLoadBalancing_;
}

inline const Foam::dsmcCoordinateSystem& Foam::dsmcCloud::coordSystem() const
{
    return dsmcCoordinateSystem_;
}

inline Foam::dsmcCoordinateSystem& Foam::dsmcCloud::coordSystem()
{
    return dsmcCoordinateSystem_();
}

inline const Foam::porousMeasurements&
Foam::dsmcCloud::porousMeas() const
{
    return porousMeasurements_;
}

inline Foam::porousMeasurements& Foam::dsmcCloud::porousMeas()
{
    return porousMeasurements_();
}

inline const Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker() const
{
    return trackingInfo_;
}

inline Foam::dsmcFaceTracker& Foam::dsmcCloud::tracker()
{
    return trackingInfo_;
}

inline const Foam::dsmcReactions& Foam::dsmcCloud::reactions() const
{
    return reactions_;
}

inline Foam::dsmcReactions& Foam::dsmcCloud::reactions()
{
    return reactions_;
}

inline Foam::scalarList Foam::dsmcCloud::infoMeasurements()
{
    scalarList infoMeasurements(7, 0.0);

    forAllIter(dsmcCloud, *this, iter)
    {
        dsmcParcel& p = iter();

        const dsmcParcel::constantProperties& cP = constProps(p.typeId());

        // Accumulate represented quantities using each parcel's own weight.
        const scalar nParticles = representedMolecules(p);

        infoMeasurements[0] += cP.mass()*nParticles;
        infoMeasurements[1] += 0.5*cP.mass()*(p.U() & p.U())*nParticles;
        infoMeasurements[2] += p.ERot()*nParticles;
        infoMeasurements[3] += cP.eVib_tot(p.vibLevel())*nParticles;
        infoMeasurements[4] += cP.electronicEnergyList()[p.ELevel()]*nParticles;

        if(p.isStuck())
        {
            infoMeasurements[5] += nParticles;
        }
        else
        {
            infoMeasurements[6] += nParticles;
        }
    }

    return infoMeasurements;
}

inline Foam::scalar Foam::dsmcCloud::maxwellianAverageSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

inline Foam::scalarField Foam::dsmcCloud::maxwellianAverageSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return 2.0*sqrt(2.0*physicoChemical::k.value()*temperature/(pi*mass));
}

inline Foam::scalar Foam::dsmcCloud::maxwellianRMSSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}

inline Foam::scalarField Foam::dsmcCloud::maxwellianRMSSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(3.0*physicoChemical::k.value()*temperature/mass);
}

inline Foam::scalar Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    const scalar temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

inline Foam::scalarField Foam::dsmcCloud::maxwellianMostProbableSpeed
(
    const scalarField& temperature,
    const scalar mass
) const
{
    return sqrt(2.0*physicoChemical::k.value()*temperature/mass);
}

inline void Foam::dsmcCloud::clear()
{
    return IDLList<dsmcParcel>::clear();
}

inline Foam::boundaryMeasurements& Foam::dsmcCloud::boundaryFluxMeasurements()
{
    return boundaryMeas_;
}

inline const Foam::boundaryMeasurements&
Foam::dsmcCloud::boundaryFluxMeasurements() const
{
    return boundaryMeas_;
}

inline Foam::cellMeasurements& Foam::dsmcCloud::cellPropMeasurements()
{
    return cellMeas_;
}

inline const Foam::cellMeasurements&
Foam::dsmcCloud::cellPropMeasurements() const
{
    return cellMeas_;
}

// ************************************************************************* //
