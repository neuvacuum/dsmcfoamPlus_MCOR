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


\*--------------------------------------------------------------------------*/

#include "dsmcLiouFangPressureOutletSpecifiedMolarFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcLiouFangPressureOutletSpecifiedMolarFraction, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcLiouFangPressureOutletSpecifiedMolarFraction, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcLiouFangPressureOutletSpecifiedMolarFraction::dsmcLiouFangPressureOutletSpecifiedMolarFraction
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcGeneralBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    outletPressure_(),
    referenceTemperature_(300.0),
    feedbackStartTimeSteps_(100),
    nTimeSteps_(scalar(0.0)),
    velocityFrame_("local"),
    typeIds_(),
    UMean_(faces_.size(), vector::zero),
    outletVelocity_(faces_.size(), vector::zero),
    UCollected_(faces_.size(), vector::zero),
    accumulatedMoleculesToInsert_(),
    totalVibrationalEnergy_(),
    vibT_(),
    vDof_(),
    nTotalParcelsSpecies_(),
    moleFractions_(),
    cellVolume_(faces_.size(), scalar(0.0)),
    outletNumberDensity_(faces_.size(), scalar(0.0)),
    outletMassDensity_(faces_.size(), scalar(0.0)),
    outletTemperature_(faces_.size(), scalar(300.0)),
    totalMomentum_(faces_.size(), vector::zero),
    totalMass_(faces_.size(), scalar(0.0)),
    totalMassDensity_(faces_.size(), scalar(0.0)),
    totalNumberDensity_(faces_.size(), scalar(0.0)),
    totalRotationalEnergy_(faces_.size(), scalar(0.0)),
    totalRotationalDof_(faces_.size(), scalar(0.0)),
    nTotalParcels_(faces_.size(), scalar(0.0)),
    nTotalParcelsInt_(faces_.size(), scalar(0.0)),
    mcc_(faces_.size(), scalar(0.0))
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];  //get volume of each boundary cell
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcLiouFangPressureOutletSpecifiedMolarFraction::~dsmcLiouFangPressureOutletSpecifiedMolarFraction()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcLiouFangPressureOutletSpecifiedMolarFraction::initialConfiguration()
{}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::calculateProperties()
{

}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    //loop over all species
    forAll(accumulatedMoleculesToInsert_, iD)
    {
        // loop over all faces of the patch
        forAll(accumulatedMoleculesToInsert_[iD], f)
        {
            scalar& accumulatedMolecules = accumulatedMoleculesToInsert_[iD][f];

            // outletVelocity_ 的存储参考系由 velocityFrame_ 决定。
            // 若它是实验室绝对系速度，真正回注到旋转 cell 前必须转换为局部相对速度。
            const vector faceVelocityModel = outletVelocity_[f];
            scalar faceTemperature = outletTemperature_[f];
            const label& faceI = faces_[f];
            const label& cellI = cells_[f];
            const scalar representedMoleculesPerParcel =
                cloud_.nParticles(patchId_, f)
               *max(cloud_.coordSystem().RWF(cellI), SMALL);

            if
            (
                !(representedMoleculesPerParcel > 0.0)
             || !(representedMoleculesPerParcel < VGREAT)
            )
            {
                FatalErrorInFunction
                    << "Invalid outlet molecular weight on patch "
                    << patchName_ << ", face " << f << ": "
                    << representedMoleculesPerParcel
                    << exit(FatalError);
            }

            // 余数始终以真实分子数保存；runtime RWF 改变时仅重新计算
            // 当前一个 parcel 代表的分子数，不能重新解释旧的 parcel 余数。
            const scalar moleculeTol =
                100.0*SMALL*representedMoleculesPerParcel;

            if
            (
                accumulatedMolecules < 0.0
             && accumulatedMolecules > -moleculeTol
            )
            {
                accumulatedMolecules = 0.0;
            }

            const scalar accumulatedParcelCount =
                accumulatedMolecules/representedMoleculesPerParcel;

            if
            (
                !(accumulatedMolecules >= 0.0)
             || !(accumulatedMolecules < VGREAT)
             || !(accumulatedParcelCount >= 0.0)
             || accumulatedParcelCount > scalar(labelMax)
            )
            {
                FatalErrorInFunction
                    << "Invalid accumulated outlet molecular flux on patch "
                    << patchName_ << ", face " << f << ": molecules="
                    << accumulatedMolecules << ", equivalentParcels="
                    << accumulatedParcelCount << nl
                    << "Check pressure, temperature, velocity and RWF feedback."
                    << exit(FatalError);
            }

            const vector& fC = mesh_.faceCentres()[faceI];
            const vector& sF = mesh_.faceAreas()[faces_[f]];
            scalar fA = mag(sF);

            List<tetIndices> faceTets = polyMeshTetDecomposition::faceTetIndices
            (
                mesh_,
                faceI,
                cellI
            );

            //Cumulative triangle area fractions
            List<scalar> cTriAFracs(faceTets.size(), 0.0);

            scalar previousCummulativeSum = 0.0;

            forAll(faceTets, triI)
            {
                const tetIndices& faceTetIs = faceTets[triI];

                cTriAFracs[triI] =
                    faceTetIs.faceTri(mesh_).mag()/fA
                    + previousCummulativeSum;

                previousCummulativeSum = cTriAFracs[triI];
            }

            //Force the last area fraction value to 1.0 to avoid any
            //rounding/non-flat face errors giving a value < 1.0
            cTriAFracs.last() = 1.0;

            //Normal unit vector *negative* so normal is pointing into the
            // domain
            vector n = sF;
            n /= -mag(n);

            //Wall tangential unit vector. Use the direction between the
            //face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            // flat and n and t1 aren't perfectly orthogonal
            vector t2 = n^t1;
            t2 /= mag(t2);

            const label nParcelsToInsert = label(accumulatedParcelCount);

            accumulatedMolecules -=
                scalar(nParcelsToInsert)*representedMoleculesPerParcel;

            if
            (
                accumulatedMolecules < 0.0
             && accumulatedMolecules > -moleculeTol
            )
            {
                accumulatedMolecules = 0.0;
            }

            const label& typeId = typeIds_[iD];
            scalar mass = cloud_.constProps(typeId).mass();

            for (label i = 0; i < nParcelsToInsert; i++)
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

                point p = faceTetIs.faceTri(mesh_).randomPoint(rndGen);

                // Liou-Fang 出口的回注粒子最终仍由 dsmcParcel::move() 按当前 cell
                // 的局部参考系推进。因此 absolute 出口速度必须在插入点 p 上转换为：
                //     Urel = Uabs - Omega ^ (p - origin)
                const vector faceVelocity =
                    velocityFrame_ == "absolute"
                  ? velocityToLocalFrame(p, faceVelocityModel, cellI)
                  : faceVelocityModel;

                // Velocity generation
                scalar mostProbableSpeed
                (
                    cloud_.maxwellianMostProbableSpeed
                    (
                        faceTemperature,
                        mass
                    )
                );

                scalar sCosTheta = (faceVelocity & n)/mostProbableSpeed;

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

                if(abs(faceVelocity & n) > VSMALL)
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
                    sqrt(physicoChemical::k.value()*faceTemperature/mass)
                    *(
                        rndGen.GaussNormal<scalar>()*t1
                        + rndGen.GaussNormal<scalar>()*t2
                    )
                    + (t1 & faceVelocity)*t1
                    + (t2 & faceVelocity)*t2
                    + mostProbableSpeed*uNormal*n;

                scalar ERot = cloud_.equipartitionRotationalEnergy
                (
                    faceTemperature,
                    cloud_.constProps(typeId).rotationalDegreesOfFreedom()
                );

                labelList vibLevel = cloud_.equipartitionVibrationalEnergyLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).nVibrationalModes(),
                    typeId
                );

                label newParcel = patchId();

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).electronicDegeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList()
                );

                const scalar& RWF = cloud_.coordSystem().RWF(cellI);

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

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsBeforeCollisions()
{

}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    const scalar sqrtPi = sqrt(pi);

    scalar molecularMass = 0.0;
    scalar molarcontantPressureSpecificHeat = 0.0;
    scalar molarcontantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label& typeId = typeIds_[iD];

        molecularMass += cloud_.constProps(typeId).mass()*moleFractions_[iD];
        molarcontantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[iD];
        molarcontantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId).rotationalDegreesOfFreedom())*moleFractions_[iD];
    }

    // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma = molarcontantPressureSpecificHeat/molarcontantVolumeSpecificHeat;

    // calculate properties in cells attached to each boundary face

    vectorField momentum(faces_.size(), vector::zero);
    vectorField UCollected(faces_.size(), vector::zero);
    scalarField mass(faces_.size(), scalar(0.0));
    scalarField nParcels(faces_.size(), scalar(0.0));
    scalarField nParcelsInt(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanX(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanY(faces_.size(), scalar(0.0));
//     scalarField velSqrMeanZ(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrX(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrY(faces_.size(), scalar(0.0));
//     scalarField velMeanSqrZ(faces_.size(), scalar(0.0));
    scalarField rotationalEnergy(faces_.size(), scalar(0.0));
    scalarField rotationalDof(faces_.size(), scalar(0.0));
    scalarField mcc(faces_.size(), scalar(0.0));
    scalarField vDoF(faces_.size(), scalar(0.0));
    scalarField translationalTemperature(faces_.size(), scalar(0.0));
    scalarField rotationalTemperature(faces_.size(), scalar(0.0));
    scalarField vibrationalTemperature(faces_.size(), scalar(0.0));
    scalarField overallTemperature(faces_.size(), scalar(0.0));
    scalarField numberDensity(faces_.size(), scalar(0.0));
    scalarField massDensity(faces_.size(), scalar(0.0));
    scalarField pressure(faces_.size(), scalar(0.0));
    scalarField speedOfSound(faces_.size(), scalar(0.0));
    scalarField velocityCorrection(faces_.size(), scalar(0.0));
    scalarField massDensityCorrection(faces_.size(), scalar(0.0));

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];

        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];

            // 出口宏观反馈量必须与 velocityFrame_ 处于同一参考系。
            // absolute 模式下，旋转 cell 内先把 parcel 相对速度恢复成实验室速度：
            //     Uabs = Urel + Omega ^ (x - origin)
            const vector UForBoundary =
                velocityFrame_ == "absolute"
              ? velocityToAbsoluteFrame(p->position(), p->U(), celli)
              : p->U();

            const scalar parcelMassWeight =
                cloud_.nParticles(celli)
               *p->RWF()
               *cloud_.constProps(p->typeId()).mass();

            momentum[c] += parcelMassWeight*UForBoundary;
            mass[c] += parcelMassWeight;
            mcc[c] +=
                parcelMassWeight*mag(UForBoundary)*mag(UForBoundary);
            nParcels[c] += p->RWF();
            UCollected[c] += p->RWF()*UForBoundary;
            // 内部能统计必须与质量、动量使用同一个 parcel RWF 权重。
            // 否则变权重算例中少量轻 parcel 会被当成全权重样本，污染
            // 出口温度反馈（当前单原子 smoke 中该项为零，但通用边界不能
            // 依赖这一特例）。
            rotationalEnergy[c] += p->RWF()*p->ERot();
            rotationalDof[c] +=
                p->RWF()
               *cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom();
//             velSqrMeanX[c] += sqr(p->U().x());
//             velSqrMeanY[c] += sqr(p->U().y());
//             velSqrMeanZ[c] += sqr(p->U().z());
//             velMeanSqrX[c] += p->U().x();
//             velMeanSqrY[c] += p->U().y();
//             velMeanSqrZ[c] += p->U().z();

            if(cloud_.constProps(p->typeId()).rotationalDegreesOfFreedom() > VSMALL)
            {
                nParcelsInt[c] += p->RWF();
            }
        }

        nTotalParcels_[c] += nParcels[c];

        nTotalParcelsInt_[c] += nParcelsInt[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        mcc_[c] += mcc[c];

        UCollected_[c] += UCollected[c];

        totalRotationalEnergy_[c] += rotationalEnergy[c];

        totalRotationalDof_[c] += rotationalDof[c];

        // 在累计样本足够前，稀疏出口 cell 的瞬时密度不能用于按
        // p=rho*R*T 反推温度。否则一个低占据 cell 会产生极高回注温度，
        // 并持续注入非物理高速 parcel。启动阶段以及短时空 cell 均使用
        // 明确给定的参考压力/温度状态；达到反馈门槛后才启用局部统计。
        if
        (
            totalMass_[c] <= VSMALL
         || nTotalParcels_[c] <= VSMALL
         || nTimeSteps_ <= feedbackStartTimeSteps_
        )
        {
            outletMassDensity_[c] =
                outletPressure_/(gasConstant*referenceTemperature_);
            outletNumberDensity_[c] =
                outletMassDensity_[c]/molecularMass;
            outletTemperature_[c] = referenceTemperature_;
            outletVelocity_[c] = vector::zero;
            continue;
        }

        massDensity[c] = totalMass_[c] / (cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if(nTotalParcels_[c] > VSMALL)
        {

//             velSqrMeanX_[c] += velSqrMeanX[c];
//             velSqrMeanY_[c] += velSqrMeanY[c];
//             velSqrMeanZ_[c] += velSqrMeanZ[c];
//             velMeanSqrX_[c] += velMeanSqrX[c];
//             velMeanSqrY_[c] += velMeanSqrY[c];
//             velMeanSqrZ_[c] += velMeanSqrZ[c];
//
//             totalVelSqrMeanX_[c] = velSqrMeanX_[c]/nTotalParcels_[c];
//             totalVelSqrMeanY_[c] = velSqrMeanY_[c]/nTotalParcels_[c];
//             totalVelSqrMeanZ_[c] = velSqrMeanZ_[c]/nTotalParcels_[c];
//             totalVelMeanSqrX_[c] = sqr(velMeanSqrX_[c]/nTotalParcels_[c]);
//             totalVelMeanSqrY_[c] = sqr(velMeanSqrY_[c]/nTotalParcels_[c]);
//             totalVelMeanSqrZ_[c] = sqr(velMeanSqrZ_[c]/nTotalParcels_[c]);
//
//             translationalTemperature[c] = (0.5*molecularMass)*(2.0/(3.0*physicoChemical::k.value()))
//                         *(
//                                 totalVelSqrMeanX_[c]+totalVelSqrMeanY_[c]+totalVelSqrMeanZ_[c]
//                                 -totalVelMeanSqrX_[c]-totalVelMeanSqrY_[c]-totalVelMeanSqrZ_[c]
//                             );

            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            translationalTemperature[c] = (1.0/(3.0*physicoChemical::k.value()))
                  *(
                      ((mcc_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))))
                      - (
                          (totalMass_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))
                        )*mag(UMean_[c])*mag(UMean_[c]))
                  );

//             Info << "translationalTemperature[c] = "  << translationalTemperature << endl;


            if
            (
                !(translationalTemperature[c] > VSMALL)
             || translationalTemperature[c] >= GREAT
            )
            {
                translationalTemperature[c] = 300.00;
            }

//             if(totalRotationalDof_[c] > VSMALL)
//             {
//                 rotationalTemperature[c] = (2.0/physicoChemical::k.value())*(totalRotationalEnergy_[c]/totalRotationalDof_[c]);
//             }
//             else
//             {
//                 rotationalTemperature[c] = 0.0;
//             }

//             forAll(totalVibrationalEnergy_, iD)
//             {
//                 const List<dsmcParcel*>& parcelsInCell = cellOccupancy[cells_[c]];
//
//                 forAll(parcelsInCell, pIC)
//                 {
//                     dsmcParcel* p = parcelsInCell[pIC];
//
//                     if(p->typeId() == typeIds_[iD])
//                     {
//                         totalVibrationalEnergy_[iD][c] += p->vibLevel()*physicoChemical::k.value()*cloud_.constProps(p->typeId()).thetaV();
//                         nTotalParcelsSpecies_[iD][c] += 1.0;
//                     }
//                 }
//             }

//             forAll(totalVibrationalEnergy_, iD)
//             {
//                 if(totalVibrationalEnergy_[iD][c] > VSMALL)
//                 {
//                     const scalar& thetaV = cloud_.constProps(typeIds_[iD]).thetaV();
//
//                     scalar vibrationalEMean = (totalVibrationalEnergy_[iD][c]/nTotalParcelsSpecies_[iD][c]);
//                     scalar iMean = vibrationalEMean/(physicoChemical::k.value()*thetaV);
//
//                     vibT_[iD][c] = thetaV / log(1.0 + (1.0/iMean));
//                     vDof_[iD][c] = (2.0*thetaV/vibT_[iD][c]) / (exp(thetaV/vibT_[iD][c]) - 1.0);
//
//                     scalar fraction = nTotalParcelsSpecies_[iD][c]/nTotalParcelsInt_[c];
//
//                     vDoF[c] += fraction*vDof_[iD][c];
//
//                     vibrationalTemperature[c] += fraction*vibT_[iD][c];
//                 }
//             }

            scalar nRotDof = 0.0;

            if(nTotalParcels_[c] > VSMALL)
            {
                nRotDof = totalRotationalDof_[c]/nTotalParcels_[c];
            }


//             overallTemperature[c] = (
//                                         (3.0*translationalTemperature[c])
//                                         + (nRotDof*rotationalTemperature[c])
//                                         + (vDoF[c]*vibrationalTemperature[c])
//                                     )
//                         / (3.0 + nRotDof + vDoF[c]);

            pressure[c] = numberDensity[c]*physicoChemical::k.value()*translationalTemperature[c];

            speedOfSound[c] = sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n /= mag(n);

//             if(faceNormalVelocity < VSMALL)
//             {
            massDensityCorrection[c] =
                (outletPressure_ - pressure[c])
               /(speedOfSound[c]*speedOfSound[c]);
            outletMassDensity_[c] =
                max
                (
                    massDensity[c] + massDensityCorrection[c],
                    VSMALL
                );
//                 }
//             else
//             {
//                 massDensityCorrection[c] = (pressure[c] - outletPressure_) / (speedOfSound[c]*speedOfSound[c]);
//                 outletMassDensity_[c] = massDensity[c] - massDensityCorrection[c];
//             }

            // Liou and Fang, 2000, equation 26 STEP 1
            outletNumberDensity_[c] = outletMassDensity_[c] / molecularMass;

            outletTemperature_[c] = outletPressure_ / (gasConstant*outletMassDensity_[c]);


            outletVelocity_[c] = totalMomentum_[c]/totalMass_[c];

            // 判断出口流向时应使用 parcel 实际穿过网格面的局部速度。
            const vector outletVelocityForFlux =
                velocityFrame_ == "absolute"
              ? velocityToLocalFrame
                (
                    mesh_.faceCentres()[faceI],
                    outletVelocity_[c],
                    celli
                )
              : outletVelocity_[c];

            const scalar faceNormalVelocity = n & outletVelocityForFlux;

            if(faceNormalVelocity < VSMALL)
            {
                velocityCorrection[c] =
                    (pressure[c] - outletPressure_)
                   /(massDensity[c]*speedOfSound[c]);
                outletVelocity_[c] += velocityCorrection[c]*n;
            }
            else
            {
                velocityCorrection[c] =
                    (outletPressure_ - pressure[c])
                   /(massDensity[c]*speedOfSound[c]);
                outletVelocity_[c] -= velocityCorrection[c]*n;
            }

//             Info << "Pressure = " << pressure[c] << endl;
//             Info << "massDensity = " << massDensity[c] << endl;
//             Info << "Temperature = " << translationalTemperature[c] << endl;
        }
    }

//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcLiouFangPressureOutletSpecifiedMolarFraction outlet velocity correction = " << velocityCorrection[(faces_.size()/2)] << endl;
//     }
//
//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcLiouFangPressureOutletSpecifiedMolarFraction outlet velocity = " << outletVelocity_[(faces_.size()/2)] << endl;
//     }

    forAll(accumulatedMoleculesToInsert_, iD)
    {
        const label& typeId = typeIds_[iD];

        forAll(accumulatedMoleculesToInsert_[iD], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT = cloud_.deltaTValue(mesh_.boundaryMesh()[patchId_].faceCells()[f]);
            scalar mass = cloud_.constProps(typeId).mass();

            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    outletTemperature_[f],
                    mass
                )
            );

                // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            // 回注粒子通量由局部推进速度控制；absolute 出口速度只在宏观反馈中保存，
            // 进入通量公式前必须按该 face owner cell 转回 local/relative 速度。
            const vector outletVelocityForFlux =
                velocityFrame_ == "absolute"
              ? velocityToLocalFrame
                (
                    mesh_.faceCentres()[faceI],
                    outletVelocity_[f],
                    mesh_.boundaryMesh()[patchId_].faceCells()[f]
                )
              : outletVelocity_[f];

            scalar sCosTheta = (outletVelocityForFlux & -sF/fA )/mostProbableSpeed;

            // Bird 4.22 给出本时间步穿过该 face 的真实分子数。
            // parcel 离散和当前 RWF 的换算统一留到 beforeMove，避免
            // RWF 更新前后的分数 parcel 余数混用不同单位。
            accumulatedMoleculesToInsert_[iD][f] +=
                moleFractions_[iD]*
                (
                    fA*outletNumberDensity_[f]*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi);
        }
    }
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::setProperties()
{
    outletPressure_ = readScalar(propsDict_.lookup("outletPressure"));

    referenceTemperature_ =
        propsDict_.lookupOrDefault<scalar>("referenceTemperature", 300.0);
    feedbackStartTimeSteps_ =
        propsDict_.lookupOrDefault<label>("feedbackStartTimeSteps", 100);

    if (referenceTemperature_ <= VSMALL || feedbackStartTimeSteps_ < 0)
    {
        FatalErrorInFunction
            << "referenceTemperature must be positive and "
            << "feedbackStartTimeSteps must be non-negative on patch "
            << patchName_ << ". Got referenceTemperature="
            << referenceTemperature_ << ", feedbackStartTimeSteps="
            << feedbackStartTimeSteps_ << exit(FatalError);
    }

    // 读取出口速度参考系。默认 local 维持旧算例行为；absolute 用于 SRF
    // 算例中希望按实验室系压力出口统计和回注速度的场景。
    velocityFrame_ =
        propsDict_.lookupOrDefault<word>("velocityFrame", word("local"));

    // velocityFrame_ 拼写错误会直接改变出口压力修正方向，因此必须立即报错。
    if (velocityFrame_ != "local" && velocityFrame_ != "absolute")
    {
        FatalErrorIn
        (
            "dsmcLiouFangPressureOutletSpecifiedMolarFraction::setProperties()"
        )
            << "velocityFrame must be either 'local' or 'absolute' for patch "
            << patchName_ << ", but got " << velocityFrame_ << nl
            << "in: " << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcLiouFangPressureOutletSpecifiedMolarFraction::setProperties()")
            << "Cannot have zero typeIds being inserd." << nl << "in: "
            << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    DynamicList<word> moleculesReduced(0);

    forAll(molecules, i)
    {
        const word& moleculeName(molecules[i]);

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
        const word& moleculeName(moleculesReduced[i]);

        label typeId(findIndex(cloud_.typeIdList(), moleculeName));

        if(typeId == -1)
        {
            FatalErrorIn("dsmcLiouFangPressureOutletSpecifiedMolarFraction::dsmcLiouFangPressureOutletSpecifiedMolarFraction()")
                << "Cannot find typeId: " << moleculeName << nl << "in: "
                << mesh_.time().system()/"boundariesDict"
                << exit(FatalError);
        }

        typeIds_[i] = typeId;
    }

    // read in the mole fraction per specie

    const dictionary& moleFractionsDict
    (
        propsDict_.subDict("moleFractions")
    );

    moleFractions_.clear();

    moleFractions_.setSize(typeIds_.size(), 0.0);

    forAll(moleFractions_, i)
    {
        moleFractions_[i] = readScalar
        (
            moleFractionsDict.lookup(moleculesReduced[i])
        );
    }

    // set the accumulator

    accumulatedMoleculesToInsert_.setSize(typeIds_.size());

    forAll(accumulatedMoleculesToInsert_, m)
    {
        accumulatedMoleculesToInsert_[m].setSize(nFaces_, 0.0);
    }

    vibT_.setSize(typeIds_.size());

    forAll(vibT_, m)
    {
        vibT_[m].setSize(nFaces_, 0.0);
    }

    vDof_.setSize(typeIds_.size());

    forAll(vDof_, m)
    {
        vDof_[m].setSize(nFaces_, 0.0);
    }

//     totalVibrationalEnergy_.setSize(typeIds_.size());
//
//     forAll(totalVibrationalEnergy_, m)
//     {
//         totalVibrationalEnergy_[m].setSize(nFaces_, 0.0);
//     }

    nTotalParcelsSpecies_.setSize(typeIds_.size());

    forAll(nTotalParcelsSpecies_, m)
    {
        nTotalParcelsSpecies_[m].setSize(nFaces_, 0.0);
    }
}

void dsmcLiouFangPressureOutletSpecifiedMolarFraction::setNewBoundaryFields()
{
    const polyPatch& patch = mesh_.boundaryMesh()[patchId_];

    //- initialise data members
    faces_.clear();
    cells_.clear();

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

    /***/

    nTimeSteps_ = 0;

    forAll(vibT_, iD)
    {
        vibT_[iD].clear();
        vibT_[iD].setSize(nFaces_, 0.0);
    }
    forAll(accumulatedMoleculesToInsert_, iD)
    {
        accumulatedMoleculesToInsert_[iD].clear();
        accumulatedMoleculesToInsert_[iD].setSize(nFaces_, 0.0);
    }

    forAll(totalVibrationalEnergy_, iD)
    {
        totalVibrationalEnergy_[iD].clear();
        totalVibrationalEnergy_[iD].setSize(nFaces_, 0.0);
    }

    forAll(vDof_, iD)
    {
        vDof_[iD].clear();
        vDof_[iD].setSize(nFaces_, 0.0);
    }

    forAll(nTotalParcelsSpecies_, iD)
    {
        nTotalParcelsSpecies_[iD].clear();
        nTotalParcelsSpecies_[iD].setSize(nFaces_, 0.0);
    }

    // recalculate volume of each boundary cell
    cellVolume_.clear();
    cellVolume_.setSize(nFaces_, 0.0);

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

    // reset the other fields
    UMean_.clear();
    outletVelocity_.clear();
    UCollected_.clear();
    outletNumberDensity_.clear();
    outletMassDensity_.clear();
    outletTemperature_.clear();
    totalMomentum_.clear();
    totalMass_.clear();
    totalMassDensity_.clear();
    totalNumberDensity_.clear();
    totalRotationalEnergy_.clear();
    totalRotationalDof_.clear();
    nTotalParcels_.clear();
    nTotalParcelsInt_.clear();
    mcc_.clear();

    UMean_.setSize(nFaces_, vector::zero);
    outletVelocity_.setSize(nFaces_, vector::zero);
    UCollected_.setSize(nFaces_, vector::zero);
    outletNumberDensity_.setSize(nFaces_, 0.0);
    outletMassDensity_.setSize(nFaces_, 0.0);
    outletTemperature_.setSize(nFaces_, referenceTemperature_);
    totalMomentum_.setSize(nFaces_, vector::zero);
    totalMass_.setSize(nFaces_, 0.0);
    totalMassDensity_.setSize(nFaces_, 0.0);
    totalNumberDensity_.setSize(nFaces_, 0.0);
    totalRotationalEnergy_.setSize(nFaces_, 0.0);
    totalRotationalDof_.setSize(nFaces_, 0.0);
    nTotalParcels_.setSize(nFaces_, 0.0);
    nTotalParcelsInt_.setSize(nFaces_, 0.0);
    mcc_.setSize(nFaces_, 0.0);
}



} // End namespace Foam

// ************************************************************************* //
