/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Copyright (C) 2016-2021 hyStrath
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef dsmcTranslationalStageInterfacePatch_H
#define dsmcTranslationalStageInterfacePatch_H

#include "dsmcCyclicBoundary.H"

namespace Foam
{

class dsmcTranslationalStageInterfacePatch
:
    public dsmcCyclicBoundary
{
private:

    dictionary propsDict_;

    word sourceRegion_;
    word targetRegion_;
    vector sourceRegionVelocity_;
    vector targetRegionVelocity_;
    scalar sourceRegionWidth_;
    scalar targetRegionWidth_;
    word positionMapping_;
    word randomPhaseOperator_;
    word randomPhaseMode_;
    label randomPhaseSeed_;
    scalar randomPhaseValue_;
    mutable scalar cachedRandomPhase_;
    mutable label cachedRandomPhaseTimeIndex_;

    bool transformPosition_;
    bool approximateFaceCentreInsertion_;
    scalar approximateInsertionFraction_;
    bool transformVelocity_;
    bool debug_;
    label maxDebugPrints_;

    label debugPrintCountA_;
    label debugPrintCountB_;
    label transferCountAtoB_;
    label transferCountBtoA_;
    label transformCountAtoB_;
    label transformCountBtoA_;
    label lastSummaryTimeIndex_;

    void readProperties();

    point transferParcel(dsmcParcel& p) const;

    point transferPatchParcel
    (
        dsmcParcel& p,
        const bool fromSideA
    ) const;

    vector transformedVelocity
    (
        const vector& U,
        const vector& fromRegionVelocity,
        const vector& toRegionVelocity
    ) const;

    scalar transferPhase() const;

public:

    TypeName("dsmcTranslationalStageInterfacePatch");

    dsmcTranslationalStageInterfacePatch
    (
        Time& t,
        const polyMesh& mesh,
        dsmcCloud& cloud,
        const dictionary& dict
    );

    virtual ~dsmcTranslationalStageInterfacePatch();

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

    void updateProperties(const dictionary& newDict);
};

} // End namespace Foam

#endif

// ************************************************************************* //
