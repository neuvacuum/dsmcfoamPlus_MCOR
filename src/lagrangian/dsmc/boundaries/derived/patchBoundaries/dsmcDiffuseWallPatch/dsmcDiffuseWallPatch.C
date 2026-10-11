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

#include "dsmcDiffuseWallPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

namespace
{

label histogramBinCount(const dictionary& propsDict)
{
    const scalar binWidth =
        propsDict.lookupOrDefault<scalar>("velocityHistogramBinWidth", 25.0);
    const scalar maxSpeed =
        propsDict.lookupOrDefault<scalar>("velocityHistogramMaxSpeed", 10000.0);

    if (binWidth <= 0 || maxSpeed <= 0)
    {
        FatalErrorIn("histogramBinCount()")
            << "velocityHistogramBinWidth and velocityHistogramMaxSpeed "
            << "must both be positive" << exit(FatalError);
    }

    return max(label(1), label(ceil(maxSpeed/binWidth)));
}

}

defineTypeNameAndDebug(dsmcDiffuseWallPatch, 0);

addToRunTimeSelectionTable
(
    dsmcPatchBoundary,
    dsmcDiffuseWallPatch,
    dictionary
);


void dsmcDiffuseWallPatch::addHistogramSample
(
    const scalar normalSpeed,
    const scalar tangentialSpeed,
    const scalar weight,
    scalarField& normalCounts,
    scalarField& tangentialCounts,
    scalar& normalOverflow,
    scalar& tangentialOverflow
)
{
    if (normalSpeed < histogramMaxSpeed_)
    {
        const label bin = min
        (
            label(normalSpeed/histogramBinWidth_),
            nHistogramBins_ - 1
        );
        normalCounts[bin] += weight;
    }
    else
    {
        normalOverflow += weight;
    }

    if (tangentialSpeed < histogramMaxSpeed_)
    {
        const label bin = min
        (
            label(tangentialSpeed/histogramBinWidth_),
            nHistogramBins_ - 1
        );
        tangentialCounts[bin] += weight;
    }
    else
    {
        tangentialOverflow += weight;
    }
}


scalarField dsmcDiffuseWallPatch::probabilityDensity
(
    const scalarField& counts
) const
{
    scalarField pdf(counts.size(), 0.0);
    const scalar inRangeWeight = sum(counts);

    if (inRangeWeight > VSMALL)
    {
        pdf = counts/(inRangeWeight*histogramBinWidth_);
    }

    return pdf;
}


void dsmcDiffuseWallPatch::resetHistograms()
{
    incidentNormalCounts_ = 0.0;
    reflectedNormalCounts_ = 0.0;
    incidentTangentialCounts_ = 0.0;
    reflectedTangentialCounts_ = 0.0;
    incidentOverflow_ = 0.0;
    reflectedOverflow_ = 0.0;
    incidentTangentialOverflow_ = 0.0;
    reflectedTangentialOverflow_ = 0.0;
}


// * * * * * * * * * * * Protected member functions  * * * * * * * * * * * * //

void dsmcDiffuseWallPatch::setProperties()
{
    velocity_ = propsDict_.lookup("velocity");

    if(propsDict_.found("groundLevelTemperature"))
    {
        temperature_ = readScalar(propsDict_.lookup("groundLevelTemperature"));
    }
    else
    {
        temperature_ = readScalar(propsDict_.lookup("temperature"));
    }

    formationLevelTemperature_ =
        propsDict_.lookupOrDefault<scalar>("formationLevelTemperature", temperature_);
}


void dsmcDiffuseWallPatch::performDiffuseReflection
(
    dsmcParcel& p,
    const scalar& localTemperature,
    const vector& localVelocity
)
{
    scalar T = getLocalTemperature(p.position()[depthAxis_]);

    if(localTemperature != 0)
    {
        T = localTemperature;
    }

    vector& U = p.U();

    scalar& ERot = p.ERot();

    labelList& vibLevel = p.vibLevel();

    label& ELevel = p.ELevel();

    const label& typeId = p.typeId();

    //- Wall unit normal vector and wall unit tangential vectors
    vector nw, tw1, tw2 = vector::zero;

    dsmcPatchBoundary::calculateWallUnitVectors(p, nw, tw1, tw2);

    const scalar mass = cloud_.constProps(typeId).mass();

    scalar rotationalDof =
        cloud_.constProps(typeId).rotationalDegreesOfFreedom();

    scalar vibrationalDof =
        cloud_.constProps(typeId).nVibrationalModes();

    Random& rndGen = cloud_.rndGen();

    U = sqrt(physicoChemical::k.value()*T/mass)
       *(
            rndGen.GaussNormal<scalar>()*tw1
          + rndGen.GaussNormal<scalar>()*tw2
          - sqrt(-2.0*log(max(1 - rndGen.sample01<scalar>(), VSMALL)))*nw
        );


    ERot = cloud_.equipartitionRotationalEnergy(T, rotationalDof);


    vibLevel =
        cloud_.equipartitionVibrationalEnergyLevel(T, vibrationalDof, typeId);


    ELevel = cloud_.equipartitionElectronicLevel
        (
            T,
            cloud_.constProps(typeId).electronicDegeneracyList(),
            cloud_.constProps(typeId).electronicEnergyList()
        );

    if (localVelocity != vector::zero)
    {
        U += localVelocity;
    }
    else
    {
        U += velocity_;
    }

    cloud_.porousMeas().diffuseInteraction(p);
}


scalar dsmcDiffuseWallPatch::getLocalTemperature
(
    const scalar& depthPosition
) const
{
    return temperature_ + (depthPosition - maxDepth_)
        *(temperature_ - formationLevelTemperature_)/lengthPatch_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
dsmcDiffuseWallPatch::dsmcDiffuseWallPatch
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcPatchBoundary(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    writeVelocityHistograms_
    (
        propsDict_.lookupOrDefault<Switch>("writeVelocityHistograms", false)
    ),
    resetHistogramsAtOutput_
    (
        propsDict_.lookupOrDefault<Switch>("resetHistogramsAtOutput", false)
    ),
    histogramBinWidth_
    (
        propsDict_.lookupOrDefault<scalar>("velocityHistogramBinWidth", 25.0)
    ),
    histogramMaxSpeed_
    (
        propsDict_.lookupOrDefault<scalar>("velocityHistogramMaxSpeed", 10000.0)
    ),
    histogramStartTime_
    (
        propsDict_.lookupOrDefault<scalar>("velocityHistogramStartTime", 0.0)
    ),
    nHistogramBins_
    (
        histogramBinCount(propsDict_)
    ),
    incidentNormalCounts_(nHistogramBins_, 0.0),
    reflectedNormalCounts_(nHistogramBins_, 0.0),
    incidentTangentialCounts_(nHistogramBins_, 0.0),
    reflectedTangentialCounts_(nHistogramBins_, 0.0),
    incidentOverflow_(0.0),
    reflectedOverflow_(0.0),
    incidentTangentialOverflow_(0.0),
    reflectedTangentialOverflow_(0.0)
{
    writeInTimeDir_ = writeVelocityHistograms_;
    writeInCase_ = false;
    measurePropertiesAtWall_ = true;

    depthAxis_ = 1;

    word depthAxisChar = propsDict_.lookupOrDefault<word>("depthAxis", "y");

    if (depthAxisChar == "x")
    {
        depthAxis_ = 0;
    }
    else if (depthAxisChar == "z")
    {
        depthAxis_ = 2;
    }

    // Here "mesh" refers to the polyMesh (patch) so the min and max positions
    // should be that of the patch bounding box
    maxDepth_ = mesh.bounds().max().component(depthAxis_);
    lengthPatch_ = maxDepth_ - mesh.bounds().min().component(depthAxis_);

    setProperties();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dsmcDiffuseWallPatch::~dsmcDiffuseWallPatch()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dsmcDiffuseWallPatch::initialConfiguration()
{}


void dsmcDiffuseWallPatch::calculateProperties()
{}


void dsmcDiffuseWallPatch::controlParticle
(
    dsmcParcel& p,
    dsmcParcel::trackingData& td
)
{
    vector nw, tw1, tw2 = vector::zero;
    dsmcPatchBoundary::calculateWallUnitVectors(p, nw, tw1, tw2);

    const vector incidentVelocity = p.U() - velocity_;

    if (writeVelocityHistograms_ && time_.value() >= histogramStartTime_)
    {
        addHistogramSample
        (
            max(incidentVelocity & nw, scalar(0)),
            sqrt
            (
                sqr(incidentVelocity & tw1)
              + sqr(incidentVelocity & tw2)
            ),
            p.RWF(),
            incidentNormalCounts_,
            incidentTangentialCounts_,
            incidentOverflow_,
            incidentTangentialOverflow_
        );
    }

    measurePropertiesBeforeControl(p);

    performDiffuseReflection(p);

    measurePropertiesAfterControl(p);

    if (writeVelocityHistograms_ && time_.value() >= histogramStartTime_)
    {
        const vector reflectedVelocity = p.U() - velocity_;

        addHistogramSample
        (
            max(-(reflectedVelocity & nw), scalar(0)),
            sqrt
            (
                sqr(reflectedVelocity & tw1)
              + sqr(reflectedVelocity & tw2)
            ),
            p.RWF(),
            reflectedNormalCounts_,
            reflectedTangentialCounts_,
            reflectedOverflow_,
            reflectedTangentialOverflow_
        );
    }
}


void dsmcDiffuseWallPatch::output
(
    const fileName& fixedPathName,
    const fileName& timePath
)
{
    if (!writeVelocityHistograms_)
    {
        return;
    }

    scalarField incidentNormal(incidentNormalCounts_);
    scalarField reflectedNormal(reflectedNormalCounts_);
    scalarField incidentTangential(incidentTangentialCounts_);
    scalarField reflectedTangential(reflectedTangentialCounts_);
    scalar incidentOverflow = incidentOverflow_;
    scalar reflectedOverflow = reflectedOverflow_;
    scalar incidentTangentialOverflow = incidentTangentialOverflow_;
    scalar reflectedTangentialOverflow = reflectedTangentialOverflow_;

    if (Pstream::parRun())
    {
        reduce(incidentNormal, sumOp<scalarField>());
        reduce(reflectedNormal, sumOp<scalarField>());
        reduce(incidentTangential, sumOp<scalarField>());
        reduce(reflectedTangential, sumOp<scalarField>());
        reduce(incidentOverflow, sumOp<scalar>());
        reduce(reflectedOverflow, sumOp<scalar>());
        reduce(incidentTangentialOverflow, sumOp<scalar>());
        reduce(reflectedTangentialOverflow, sumOp<scalar>());
    }

    if (Pstream::master())
    {
        scalarField binCentres(nHistogramBins_, 0.0);
        forAll(binCentres, bin)
        {
            binCentres[bin] = (scalar(bin) + 0.5)*histogramBinWidth_;
        }

        Foam::writeTimeData
        (
            timePath,
            "incidentNormalVelocityPdf",
            binCentres,
            probabilityDensity(incidentNormal)
        );
        Foam::writeTimeData
        (
            timePath,
            "reflectedNormalVelocityPdf",
            binCentres,
            probabilityDensity(reflectedNormal)
        );
        Foam::writeTimeData
        (
            timePath,
            "incidentTangentialSpeedPdf",
            binCentres,
            probabilityDensity(incidentTangential)
        );
        Foam::writeTimeData
        (
            timePath,
            "reflectedTangentialSpeedPdf",
            binCentres,
            probabilityDensity(reflectedTangential)
        );

        OFstream summary(timePath/"wallVelocityHistogramSummary.dat");
        summary
            << "# in-range incident normal weight, overflow weight, "
            << "in-range reflected normal weight, overflow weight, "
            << "incident tangential overflow weight, "
            << "reflected tangential overflow weight" << nl
            << sum(incidentNormal) << ' ' << incidentOverflow << ' '
            << sum(reflectedNormal) << ' ' << reflectedOverflow << ' '
            << incidentTangentialOverflow << ' '
            << reflectedTangentialOverflow << nl;
    }

    if (resetHistogramsAtOutput_)
    {
        resetHistograms();
    }
}


void dsmcDiffuseWallPatch::updateProperties(const dictionary& newDict)
{
    //- the main properties should be updated first
    updateBoundaryProperties(newDict);

    propsDict_ = newDict.subDict(typeName + "Properties");

    setProperties();
}

} // End namespace Foam

// ************************************************************************* //
