/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

Class
    dsmcSpeedHistogramZone

Description
    Accumulates a represented-molecule-weighted speed probability density in
    either a cell zone or the cells adjacent to selected boundary patches.

\*---------------------------------------------------------------------------*/

#ifndef dsmcSpeedHistogramZone_H
#define dsmcSpeedHistogramZone_H

#include "dsmcField.H"
#include "dsmcCloud.H"
#include "scalarField.H"

namespace Foam
{

class dsmcSpeedHistogramZone
:
    public dsmcField
{
    dictionary propsDict_;
    word fieldName_;
    labelField typeIds_;
    labelList selectedCells_;
    bool absoluteVelocity_;
    scalar binWidth_;
    scalar maxSpeed_;
    scalar samplingStartTime_;
    label nBins_;
    scalarField counts_;
    scalar overflowWeight_;

    void selectCells();
    void addSample(const scalar speed, const scalar weight);
    scalarField probabilityDensity(const scalarField& counts) const;

public:

    TypeName("dsmcSpeedHistogramZone");

    dsmcSpeedHistogramZone
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    ~dsmcSpeedHistogramZone();

    virtual void createField();
    virtual void calculateField();
    virtual void resetField();
    virtual void writeField();
    virtual void updateProperties(const dictionary& dict);
    virtual scalar translationalT(const label cellI);
    virtual scalar overallT(const label cellI);
};

} // End namespace Foam

#endif

// ************************************************************************* //
