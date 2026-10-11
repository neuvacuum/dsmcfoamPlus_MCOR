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

#include "dsmcWangPressureInlet.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

using namespace Foam::constant::mathematical;

namespace Foam
{

defineTypeNameAndDebug(dsmcWangPressureInlet, 0);

addToRunTimeSelectionTable(dsmcGeneralBoundary, dsmcWangPressureInlet, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcWangPressureInlet::dsmcWangPressureInlet
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
    moleFractions_(),
    velocityFrame_("local"),
    inletPressure_(),
    inletTemperature_(),
    n_(),
    cellVolume_(faces_.size(), scalar(0.0)),
    accumulatedMoleculesToInsert_(),
    inletVelocity_(faces_.size(), vector::zero),
    totalMomentum_(faces_.size(), vector::zero),
    totalMass_(faces_.size(), scalar(0.0)),
    nTotalParcels_(faces_.size(), scalar(0.0)),
    nTimeSteps_(scalar(0.0)),
    mcc_(faces_.size(), scalar(0.0)),
    UMean_(faces_.size(), vector::zero),
    UCollected_(faces_.size(), vector::zero)
{
    writeInTimeDir_ = false;
    writeInCase_ = true;

    setProperties();

    // calculate required number density at inlet boundary
    // equation 32, Liou and Fang (2000)
    n_ = inletPressure_ / (physicoChemical::k.value()*inletTemperature_);

    //get volume of each boundary cell
    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcWangPressureInlet::~dsmcWangPressureInlet()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void dsmcWangPressureInlet::initialConfiguration()
{}

void dsmcWangPressureInlet::calculateProperties()
{

}

void dsmcWangPressureInlet::controlParcelsBeforeMove()
{
    Random& rndGen = cloud_.rndGen();

    //loop over all species
    forAll(accumulatedMoleculesToInsert_, iD)
    {
        // loop over all faces of the patch
        forAll(accumulatedMoleculesToInsert_[iD], f)
        {
            scalar& accumulatedMolecules = accumulatedMoleculesToInsert_[iD][f];

            // inletVelocity_ 的存储参考系由 velocityFrame_ 决定。
            // local    ：历史行为，后续直接作为 parcel 存储速度使用；
            // absolute ：实验室系入口速度，真正注入旋转 cell 前需按粒子位置转成相对速度。
            const vector faceVelocityModel = inletVelocity_[f];
            scalar faceTemperature = inletTemperature_;
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
                    << "Invalid inlet molecular weight on patch "
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
                    << "Invalid accumulated inlet molecular flux on patch "
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
            //domain
            vector n = sF;
            n /= -mag(n);

            // Wall tangential unit vector. Use the direction between the
            // face centre and the first vertex in the list
            vector t1 = fC - mesh_.points()[mesh_.faces()[faceI][0]];
            t1 /= mag(t1);

            //Other tangential unit vector.  Rescaling in case face is not
            //flat and n and t1 aren't perfectly orthogonal
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

                // 对旋转域入口，Bird 入口采样公式必须使用 parcel 实际存储/推进的速度。
                // 若用户把边界速度声明为 absolute，则这里在随机插入点 p 上执行
                //     Urel = Uabs - Omega ^ (p - origin)
                // 从而保留径向位置导致的局部刚体速度差异。
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

                label ELevel = cloud_.equipartitionElectronicLevel
                (
                    faceTemperature,
                    cloud_.constProps(typeId).electronicDegeneracyList(),
                    cloud_.constProps(typeId).electronicEnergyList()
                );

                label newParcel = patchId();

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

void dsmcWangPressureInlet::controlParcelsBeforeCollisions()
{

}

void dsmcWangPressureInlet::controlParcelsAfterCollisions()
{
    nTimeSteps_ += 1.0;

    const scalar sqrtPi = sqrt(pi);

    scalar molecularMass = 0.0;
    scalar molarconstantPressureSpecificHeat = 0.0;
    scalar molarconstantVolumeSpecificHeat = 0.0;

    forAll(moleFractions_, iD)
    {
        const label& typeId_ = typeIds_[iD];

        molecularMass += cloud_.constProps(typeId_).mass()*moleFractions_[iD];
        molarconstantPressureSpecificHeat += (5.0 + cloud_.constProps(typeId_).rotationalDegreesOfFreedom())*moleFractions_[iD];
        molarconstantVolumeSpecificHeat += (3.0 + cloud_.constProps(typeId_).rotationalDegreesOfFreedom())*moleFractions_[iD];
    }

     // R = k/m
    const scalar gasConstant = physicoChemical::k.value()/molecularMass;

    const scalar gamma = molarconstantPressureSpecificHeat/molarconstantVolumeSpecificHeat;

    vectorField momentum(faces_.size(), vector::zero);
    vectorField UCollected(faces_.size(), vector::zero);
    scalarField mass(faces_.size(), scalar(0.0));
    scalarField nParcels(faces_.size(), scalar(0.0));
    scalarField mcc(faces_.size(), scalar(0.0));
    scalarField translationalTemperature(faces_.size(), scalar(0.0));
    scalarField numberDensity(faces_.size(), scalar(0.0));
    scalarField massDensity(faces_.size(), scalar(0.0));
    scalarField pressure(faces_.size(), scalar(0.0));
    scalarField speedOfSound(faces_.size(), scalar(0.0));
    scalarField velocityCorrection(faces_.size(), scalar(0.0));

    const List<DynamicList<dsmcParcel*> >& cellOccupancy = cloud_.cellOccupancy();

    forAll(cells_, c)
    {
        const label celli = cells_[c];

        const List<dsmcParcel*>& parcelsInCell = cellOccupancy[celli];

        forAll(parcelsInCell, pIC)
        {
            dsmcParcel* p = parcelsInCell[pIC];
            label iD = findIndex(typeIds_, p->typeId());

            if(iD != -1)
            {
                // 压力入口的宏观反馈量应与 velocityFrame_ 保持同一语义。
                // local    ：沿用历史行为，直接统计 parcel 存储速度；
                // absolute ：旋转 cell 内先恢复 Uabs = Urel + Omega ^ r，再参与压力修正。
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
            }
        }

        nTotalParcels_[c] += nParcels[c];

        totalMomentum_[c] += momentum[c];

        totalMass_[c] += mass[c];

        mcc_[c] += mcc[c];

        UCollected_[c] += UCollected[c];

        // 变权重重映射或开边界流动可能让局部 cell 暂时为空；没有有效
        // 统计量时沿用入口给定状态，避免压力反馈中的除零。
        if
        (
            totalMass_[c] <= VSMALL
         || nTotalParcels_[c] <= VSMALL
        )
        {
            continue;
        }

        massDensity[c] = totalMass_[c] / (cellVolume_[c]*nTimeSteps_);

        numberDensity[c] = massDensity[c]/molecularMass;

        if(nTotalParcels_[c] > VSMALL)
        {
            UMean_[c] = UCollected_[c]/nTotalParcels_[c];

            translationalTemperature[c] = (1.0/(3.0*physicoChemical::k.value()))
                *(
                    ((mcc_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))))
                    - (
                        (totalMass_[c]/(nTotalParcels_[c]*cloud_.nParticles(celli))
                        )*mag(UMean_[c])*mag(UMean_[c]))
                );

            if
            (
                !(translationalTemperature[c] > VSMALL)
             || translationalTemperature[c] >= GREAT
            )
            {
                translationalTemperature[c] = 300.00;
            }

            pressure[c] = numberDensity[c]*physicoChemical::k.value()*translationalTemperature[c];


            speedOfSound[c] = sqrt(gamma*gasConstant*translationalTemperature[c]);

            label faceI = faces_[c];
            vector n = mesh_.faceAreas()[faceI];
            n /= mag(n);

//             scalar faceNormalVelocity = (n & inletVelocity_[c]);

            inletVelocity_[c] = totalMomentum_[c]/totalMass_[c];

            if(nTimeSteps_ > 100)
            {
            //velocity correction for each boundary cellI
//             if(faceNormalVelocity < VSMALL)
//             {
                velocityCorrection[c] = (pressure[c] - inletPressure_) / (massDensity[c]*speedOfSound[c]);
                inletVelocity_[c] += velocityCorrection[c]*n;
//             }
//             else
//             {
//                 velocityCorrection[c] = (inletPressure_ - pressure[c]) / (massDensity[c]*speedOfSound[c]);
//                 inletVelocity_[c] -= velocityCorrection[c]*n;
//             }
            }
        }
    }

//     if(faces_.size() > VSMALL)
//     {
//         Pout << "dsmcWangPressureInlet inlet velocity correction = " << velocityCorrection[(faces_.size()/2)] << endl;
//         Pout << "dsmcWangPressureInlet inlet velocity = " << inletVelocity_[(faces_.size()/2)] << endl;
//     }

    // compute number of parcels to insert
    forAll(accumulatedMoleculesToInsert_, iD)
    {
        const label& typeId = typeIds_[iD];
        scalar mass = cloud_.constProps(typeId).mass();

        forAll(accumulatedMoleculesToInsert_[iD], f)
        {
            const label& faceI = faces_[f];
            const vector& sF = mesh_.faceAreas()[faceI];
            const scalar fA = mag(sF);

            const scalar deltaT = cloud_.deltaTValue(mesh_.boundaryMesh()[patchId_].faceCells()[f]);
            scalar mostProbableSpeed
            (
                cloud_.maxwellianMostProbableSpeed
                (
                    inletTemperature_,
                    mass
                )
            );

            // Dotting boundary velocity with the face unit normal
            // (which points out of the domain, so it must be
            // negated), dividing by the most probable speed to form
            // molecularSpeedRatio * cosTheta

            // 注入通量公式使用的是 parcel 穿过当前网格面的局部推进速度。
            // 当入口速度以 absolute 给出且 owner cell 处于旋转域时，需要在面心处
            // 先转换成 local/relative 速度，再取法向分量。
            const vector inletVelocityForFlux =
                velocityFrame_ == "absolute"
              ? velocityToLocalFrame
                (
                    mesh_.faceCentres()[faceI],
                    inletVelocity_[f],
                    mesh_.boundaryMesh()[patchId_].faceCells()[f]
                )
              : inletVelocity_[f];

            scalar sCosTheta = (inletVelocityForFlux & -sF/fA )/mostProbableSpeed;

            // Bird 4.22 给出本时间步穿过该 face 的真实分子数。
            // parcel 离散和当前 RWF 的换算统一留到 beforeMove，避免
            // RWF 更新前后的分数 parcel 余数混用不同单位。
            accumulatedMoleculesToInsert_[iD][f] +=
                moleFractions_[iD]*
                (
                    fA*n_*deltaT*mostProbableSpeed
                    *
                    (
                        exp(-sqr(sCosTheta)) + sqrtPi*sCosTheta*(1 + erf(sCosTheta))
                    )
                )
                /(2.0*sqrtPi);
        }
    }
}

void dsmcWangPressureInlet::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
}

void dsmcWangPressureInlet::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    setProperties();
}

void dsmcWangPressureInlet::setProperties()
{
    inletPressure_ = readScalar(propsDict_.lookup("inletPressure"));

    inletTemperature_ = readScalar(propsDict_.lookup("inletTemperature"));

    // 读取边界速度参考系。默认保持 local，确保历史算例不因新增 SRF 支持而改变。
    velocityFrame_ =
        propsDict_.lookupOrDefault<word>("velocityFrame", word("local"));

    // 只接受 local / absolute 两种明确语义，避免拼写错误静默污染入口速度。
    if (velocityFrame_ != "local" && velocityFrame_ != "absolute")
    {
        FatalErrorIn("dsmcWangPressureInlet::setProperties()")
            << "velocityFrame must be either 'local' or 'absolute' for patch "
            << patchName_ << ", but got " << velocityFrame_ << nl
            << "in: " << mesh_.time().system()/"boundariesDict"
            << exit(FatalError);
    }

    const List<word> molecules (propsDict_.lookup("typeIds"));

    if(molecules.size() == 0)
    {
        FatalErrorIn("dsmcWangPressureInlet::setProperties()")
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
            FatalErrorIn("dsmcWangPressureInlet::dsmcWangPressureInlet()")
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

    moleFractions_.clear();  //TODO VINCENT molefractions is not a list of scalarFields!

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
}


void dsmcWangPressureInlet::setNewBoundaryFields()
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

    forAll(accumulatedMoleculesToInsert_, iD)
    {
        accumulatedMoleculesToInsert_[iD].clear();

        accumulatedMoleculesToInsert_[iD].setSize(nFaces_, 0.0);
    }

    // recalculate volume of each boundary cell
    cellVolume_.clear();
    cellVolume_.setSize(nFaces_, 0.0);

    forAll(cellVolume_, c)
    {
        cellVolume_[c] = mesh_.cellVolumes()[cells_[c]];
    }

    // reset other fields
    inletVelocity_.clear();
    totalMomentum_.clear();
    totalMass_.clear();
    nTotalParcels_.clear();
    mcc_.clear();
    UMean_.clear();
    UCollected_.clear();

    inletVelocity_.setSize(nFaces_, vector::zero);
    totalMomentum_.setSize(nFaces_, vector::zero);
    totalMass_.setSize(nFaces_, 0.0);
    nTotalParcels_.setSize(nFaces_, 0.0);
    mcc_.setSize(nFaces_, 0.0);
    UMean_.setSize(nFaces_, vector::zero);
    UCollected_.setSize(nFaces_, vector::zero);
}


} // End namespace Foam

// ************************************************************************* //
