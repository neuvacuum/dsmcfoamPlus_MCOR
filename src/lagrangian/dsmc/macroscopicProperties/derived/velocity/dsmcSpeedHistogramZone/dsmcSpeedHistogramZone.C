/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "dsmcSpeedHistogramZone.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"
#include "writeTimeData.H"

namespace Foam
{

defineTypeNameAndDebug(dsmcSpeedHistogramZone, 0);
addToRunTimeSelectionTable(dsmcField, dsmcSpeedHistogramZone, dictionary);


void dsmcSpeedHistogramZone::selectCells()
{
    const word selectionMode =
        propsDict_.lookupOrDefault<word>("selectionMode", "cellZone");

    if (selectionMode == "cellZone")
    {
        const word zoneName(propsDict_.lookup("zoneName"));
        const label zoneId = mesh_.cellZones().findZoneID(zoneName);

        if (zoneId == -1)
        {
            FatalErrorIn("dsmcSpeedHistogramZone::selectCells()")
                << "Cannot find cell zone " << zoneName << exit(FatalError);
        }

        selectedCells_ = mesh_.cellZones()[zoneId];
    }
    else if (selectionMode == "boundaryAdjacentCells")
    {
        const wordList patchNames(propsDict_.lookup("patchNames"));
        boolList selected(mesh_.nCells(), false);

        forAll(patchNames, patchI)
        {
            const label patchId = mesh_.boundaryMesh().findPatchID
            (
                patchNames[patchI]
            );

            if (patchId == -1)
            {
                FatalErrorIn("dsmcSpeedHistogramZone::selectCells()")
                    << "Cannot find boundary patch " << patchNames[patchI]
                    << exit(FatalError);
            }

            const labelUList& faceCells =
                mesh_.boundaryMesh()[patchId].faceCells();

            forAll(faceCells, faceI)
            {
                selected[faceCells[faceI]] = true;
            }
        }

        DynamicList<label> cells;
        forAll(selected, cellI)
        {
            if (selected[cellI])
            {
                cells.append(cellI);
            }
        }
        cells.shrink();
        selectedCells_.transfer(cells);
    }
    else
    {
        FatalErrorIn("dsmcSpeedHistogramZone::selectCells()")
            << "Unknown selectionMode " << selectionMode
            << ". Valid values are cellZone and boundaryAdjacentCells"
            << exit(FatalError);
    }
}


void dsmcSpeedHistogramZone::addSample
(
    const scalar speed,
    const scalar weight
)
{
    if (speed < maxSpeed_)
    {
        const label bin = min(label(speed/binWidth_), nBins_ - 1);
        counts_[bin] += weight;
    }
    else
    {
        overflowWeight_ += weight;
    }
}


scalarField dsmcSpeedHistogramZone::probabilityDensity
(
    const scalarField& counts
) const
{
    scalarField pdf(counts.size(), 0.0);
    const scalar inRangeWeight = sum(counts);

    if (inRangeWeight > VSMALL)
    {
        pdf = counts/(inRangeWeight*binWidth_);
    }

    return pdf;
}


dsmcSpeedHistogramZone::dsmcSpeedHistogramZone
(
    Time& t,
    const polyMesh& mesh,
    dsmcCloud& cloud,
    const dictionary& dict
)
:
    dsmcField(t, mesh, cloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
    fieldName_(propsDict_.lookup("fieldName")),
    typeIds_(),
    selectedCells_(),
    absoluteVelocity_
    (
        propsDict_.lookupOrDefault<word>("velocityFrame", "relative")
     == "absolute"
    ),
    binWidth_(readScalar(propsDict_.lookup("binWidth"))),
    maxSpeed_(readScalar(propsDict_.lookup("maxSpeed"))),
    samplingStartTime_
    (
        propsDict_.lookupOrDefault<scalar>("samplingStartTime", 0.0)
    ),
    nBins_(0),
    counts_(),
    overflowWeight_(0.0)
{
    if (binWidth_ <= 0 || maxSpeed_ <= 0)
    {
        FatalErrorIn("dsmcSpeedHistogramZone::dsmcSpeedHistogramZone()")
            << "binWidth and maxSpeed must both be positive"
            << exit(FatalError);
    }

    const word velocityFrame =
        propsDict_.lookupOrDefault<word>("velocityFrame", "relative");
    if (velocityFrame != "relative" && velocityFrame != "absolute")
    {
        FatalErrorIn("dsmcSpeedHistogramZone::dsmcSpeedHistogramZone()")
            << "velocityFrame must be relative or absolute"
            << exit(FatalError);
    }

    const wordList molecules(propsDict_.lookup("typeIds"));
    typeIds_.setSize(molecules.size(), -1);
    forAll(molecules, i)
    {
        const label typeId = findIndex(cloud_.typeIdList(), molecules[i]);
        if (typeId == -1)
        {
            FatalErrorIn("dsmcSpeedHistogramZone::dsmcSpeedHistogramZone()")
                << "Cannot find typeId " << molecules[i] << exit(FatalError);
        }
        typeIds_[i] = typeId;
    }

    nBins_ = max(label(1), label(ceil(maxSpeed_/binWidth_)));
    counts_.setSize(nBins_, 0.0);
    selectCells();
}


dsmcSpeedHistogramZone::~dsmcSpeedHistogramZone()
{}


void dsmcSpeedHistogramZone::createField()
{
    Info<< "Initialising dsmcSpeedHistogramZone " << fieldName_
        << " with " << selectedCells_.size() << " local cells" << endl;
}


void dsmcSpeedHistogramZone::calculateField()
{
    if
    (
        !time_.samplingTime()
     || time_.time().value() < samplingStartTime_
    )
    {
        return;
    }

    const List<DynamicList<dsmcParcel*> >& occupancy = cloud_.cellOccupancy();

    forAll(selectedCells_, selectedI)
    {
        const label cellI = selectedCells_[selectedI];
        const List<dsmcParcel*>& parcels = occupancy[cellI];

        forAll(parcels, parcelI)
        {
            const dsmcParcel& parcel = *parcels[parcelI];
            if (findIndex(typeIds_, parcel.typeId()) == -1 || !parcel.isFree())
            {
                continue;
            }

            vector velocity = parcel.U();
            if (absoluteVelocity_ && cloud_.isInRotatingRegion(cellI))
            {
                velocity +=
                    cloud_.rotatingFrameOmega()
                  ^ (parcel.position() - cloud_.rotatingFrameOrigin());
            }

            addSample(mag(velocity), cloud_.representedMolecules(parcel));
        }
    }
}


void dsmcSpeedHistogramZone::resetField()
{
    counts_ = 0.0;
    overflowWeight_ = 0.0;
}


void dsmcSpeedHistogramZone::writeField()
{
    if (!time_.time().outputTime())
    {
        return;
    }

    scalarField counts(counts_);
    scalar overflowWeight = overflowWeight_;

    if (Pstream::parRun())
    {
        reduce(counts, sumOp<scalarField>());
        reduce(overflowWeight, sumOp<scalar>());
    }

    if (Pstream::master())
    {
        scalarField binCentres(nBins_, 0.0);
        forAll(binCentres, bin)
        {
            binCentres[bin] = (scalar(bin) + 0.5)*binWidth_;
        }

        writeTimeData
        (
            timePath_,
            "speedProbabilityDensity_" + fieldName_,
            binCentres,
            probabilityDensity(counts)
        );

        scalarField summaryX(3, 0.0);
        scalarField summaryY(3, 0.0);
        summaryX[0] = 0;
        summaryX[1] = 1;
        summaryX[2] = 2;
        summaryY[0] = sum(counts);
        summaryY[1] = overflowWeight;
        summaryY[2] = sum(counts) + overflowWeight;
        writeTimeData
        (
            timePath_,
            "speedHistogramSummary_" + fieldName_,
            summaryX,
            summaryY
        );
    }

    if (time_.resetFieldsAtOutput())
    {
        resetField();
    }
}


void dsmcSpeedHistogramZone::updateProperties(const dictionary& newDict)
{
    updateBasicFieldProperties(newDict);
}


scalar dsmcSpeedHistogramZone::translationalT(const label)
{
    return -1.0;
}


scalar dsmcSpeedHistogramZone::overallT(const label)
{
    return -1.0;
}

} // End namespace Foam

// ************************************************************************* //
