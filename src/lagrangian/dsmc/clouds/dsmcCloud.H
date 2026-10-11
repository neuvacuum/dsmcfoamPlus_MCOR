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
    Foam::dsmcCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    dsmcCloudI.H
    dsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef dsmcCloud_H
#define dsmcCloud_H

#include "Cloud.H"

#include "dsmcParcel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "dynamicFvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"
#include "BinaryCollisionModel.H"

#include "dsmcBoundaries.H"
#include "dsmcAllConfigurations.H"
#include "dsmcFieldProperties.H"
#include "dsmcFaceTracker.H"
#include "dsmcControllers.H"
#include "dsmcCoordinateSystem.H"
#include "porousMeasurements.H"
#include "dsmcDynamicLoadBalancing.H"

#include "dsmcReactions.H"
#include "collisionPartnerSelection.H"
#include "boundaryMeasurements.H"
#include "cellMeasurements.H"

#include "molsToDelete.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dsmcCloud Declaration
\*---------------------------------------------------------------------------*/

class dsmcCloud
:
    public Cloud<dsmcParcel>
{

private:

    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word& cloudName_;

        //- References to the mesh and time databases
        const dynamicFvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Dictionary of simulation properties
        dictionary controlDict_;

        //- A list of unique instances of molecule types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Coordinate system model: Cartesian, axisymmetric, spherical
        autoPtr<dsmcCoordinateSystem> dsmcCoordinateSystem_;

        //- Porous measurements model: no, instantaneous, global
        autoPtr<porousMeasurements> porousMeasurements_;

        //- Write out the information to terminal every nTerminalOutputs
        label nTerminalOutputs_;

        //- Optional rotating-reference-frame support
        bool rotatingReferenceFrameActive_;
        word rotatingRegionZoneName_;
        vector rotatingFrameOmega_;
        vector rotatingFrameOrigin_;
        label rotatingFrameSubCycles_;
        scalar rotatingFrameMaxSubStep_;
        scalar rotatingFrameMaxRotationAngle_;
        scalar rotatingFrameMeshTolFraction_;
        boolList rotatingRegionCellMask_;

        //- A data structure holding which particles are in which cell
        //List<DynamicList<dsmcParcel*> > cellOccupancy_;
        DynamicList<DynamicList<dsmcParcel*>> cellOccupancy_; // EDITED VINCENT

        //- A data structure holding electron number density
        scalarList rhoNMeanElectron_;

        //- A data structure holding electron mass density
        scalarList rhoMMeanElectron_;

        //- A data structure holding all species mass density
        scalarList rhoMMean_;

        //- A data structure holding electron momentum
        List<vector> momentumMeanElectron_;

        //- A data structure holding all species momentum
        List<vector> momentumMean_;

        //- A data structure holding electron linearKE
        scalarList linearKEMeanElectron_;

        //- A data structure holding electron temperature
        scalarList electronTemperature_;

        //- A data structure holding cell velocity
        List<vector> cellVelocity_;

        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Parcel constant properties - one for each type
        List<dsmcParcel::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- Controllers
        dsmcControllers controllers_;

        //- Load balancing
        dsmcDynamicLoadBalancing dynamicLoadBalancing_;

        //- Boundary measurements
        boundaryMeasurements boundaryMeas_;

        //- Measurement fields
        //dsmcStandardFields standardFields_;

        dsmcFieldProperties fields_;

        //- Boundaries
        dsmcBoundaries boundaries_;

        //- Flux-face tracking
        dsmcFaceTracker trackingInfo_;

      // References to the cloud sub-models

        //- Binary collision model
        autoPtr<BinaryCollisionModel> binaryCollisionModel_;

        //- Collision selection routine model
        autoPtr<collisionPartnerSelection> collisionPartnerSelectionModel_;

        //- Chemical reactions
        dsmcReactions reactions_;

        //- Calculation of the mean collision separation
        cellMeasurements cellMeas_;

    // Private Member Functions

        label pickFromCandidateList
        (
            DynamicList<label>& candidatesInCell
        );

        // Whenever you pick from the candidate list
        // always update the candidate sublist
        void updateCandidateSubList
        (
            const label& candidate,
            DynamicList<label>& candidatesInSubCell
        );

        label pickFromCandidateSubList
        (
            DynamicList<label>& candidatesInCell,
            DynamicList<label>& candidatesInSubCell
        );

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Read optional rotating reference frame configuration
        void initialiseRotatingReferenceFrame();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- When mesh changed, relocate the patch/face of stuck parcels
        void relocateStuckParcels();

        //- When mesh has changed, changed the size of cellOccupancy_
        void buildCellOccupancyFromScratch();

        //- When mesh has changed, change the size of collisionSelectionRemainder_
        void buildCollisionSelectionRemainderFromScratch();

        //- When mesh has changed, update the boundaries
        void resetBoundaries();

        //- When mesh has changed, update the measurement tools
        void resetMeasurementTools();

        //- Remove electrons after adding their velocities to a DynamicList
        void removeElectrons();

        //- Add an electron where an ion is found, after the move function
        void addElectrons();

        //- Calculate collisions between molecules
        void collisions();

        //- Disallow default bitwise copy construct
        dsmcCloud(const dsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const dsmcCloud&);

public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //  file
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const dynamicFvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        dsmcCloud
        (
            Time& t,
            const word& cloudName,
            const dynamicFvMesh& mesh,
            const IOdictionary& dsmcInitialiseDict,
            const bool& clearFields
        );

    //- Destructor
     virtual ~dsmcCloud();

    //- Type of parcel the cloud was instantiated for
//     typedef dsmcParcel parcelType;

    // Member Functions

        // Access

    // References to the mesh and databases

        //- Return the cloud type
        inline const word& cloudName() const;

        //- Return refernce to the mesh
        inline const dynamicFvMesh& mesh() const;

    // References to the dsmc specific data

        //- Return particle properties dictionary
        inline const IOdictionary& particleProperties() const;

        //- Return particle properties dictionary
        inline const dictionary& controlDict() const;

        //- Return the idList
        inline const List<word>& typeIdList() const;

        //- Return the number of real particles represented by one parcel
        inline scalar nParticle() const;

        //- Return the number of real particles represented by one parcel
        inline const volScalarField& nParticles() const;

        //- Return the number of real particles represented by one parcel for
        //  a given cell
        inline scalar nParticles(const label celli) const;

        //- Return the number of real particles represented by one parcel for
        //  a given patch/face couple
        inline scalar nParticles(const label patchi, const label facei) const;

        // Return the time-step-model reference weight multiplied by the parcel RWF.
        inline scalar representedMolecules(const dsmcParcel& parcel) const;

        //- Return the DSMC time-step field
        inline const scalarField& deltaTField() const;

        //- Return the DSMC time-step for a given cell
        inline scalar deltaTValue(const label celli = 0) const;

        //- Return true if the optional rotating-reference-frame support is active
        inline bool rotatingReferenceFrameActive() const;

        //- Return true if the cell belongs to the configured rotating region
        inline bool isInRotatingRegion(const label celli) const;

        //- Return the configured angular velocity vector
        inline const vector& rotatingFrameOmega() const;

        //- Return the configured rotation origin
        inline const vector& rotatingFrameOrigin() const;

        //- Return the configured minimum number of sub-cycles per move step
        inline label rotatingFrameSubCycles() const;

        //- Return the configured absolute upper bound for SRF sub-step size
        inline scalar rotatingFrameMaxSubStep() const;

        //- Apparent acceleration in the configured rotating reference frame
        vector rotatingFrameAcceleration
        (
            const point& position,
            const vector& Urel
        ) const;

        //- Time sub-step used for SRF parcel tracking in a cell
        scalar rotatingFrameTrackSubStep
        (
            const label celli,
            const point& position,
            const vector& Urel,
            const scalar remainingTime
        ) const;

        //- Return the number of timesteps to write out to terminal at
        inline label nTerminalOutputs() const;

        //- Return the cell occupancy addressing
        inline const DynamicList<DynamicList<dsmcParcel*> >& //EDITED VINCENT
        cellOccupancy() const;

        //- Return the sigmaTcRMax field.  non-const access to allow
        // updating.
        inline volScalarField& sigmaTcRMax();

        //- Return the collision selection remainder field.  non-const
        // access to allow updating.
        inline scalarField& collisionSelectionRemainder();

        //- Return all of the constant properties
        inline const List<dsmcParcel::constantProperties>&
        constProps() const;

        //- Return the constant properties of the given typeId
        inline const dsmcParcel::constantProperties&
        constProps(label typeId) const;

        //- Return refernce to the random object
        inline Random& rndGen();

    // Random implementation corrections
        //Generates a random integer between valMin and valMax
        label randomLabel
        (
            const label valOne,
            const label valTwo
        );

    // Kinetic theory helper functions

        //- Generate a random velocity sampled from the Maxwellian speed
        // distribution
        vector equipartitionLinearVelocity
        (
            const scalar temperature,
            const scalar mass
        );

        //- Generate a random velocity sampled from the Chapman-Enskog speed
        // distribution
        vector chapmanEnskogVelocity
        (
            const scalar temperature,
            const scalar mass,
            const vector& q,
            const tensor& tau
        );

        //- Generate a random velocity sampled from the Generalised
        //  Chapman-Enskog distribution speed
        void generalisedChapmanEnskog
        (
            const label& typeID,
            const scalar& translationalTemperature,
            const scalar& rotationalTemperature,
            const scalar& vibrationalTemperature,
            const scalar& mass,
            const vector& D,
            const vector& qTra,
            const vector& qRot,
            const vector& qVib,
            const tensor& tau,
            scalar& ERot,
            labelList& vibLevel,
            vector& U
        );

        //- Generate a random rotational energy, sampled from the
        // equilibrium distribution (Bird eqn 11.22 and 11.23 and
        // adapting code from DSMC3.FOR)
        scalar equipartitionRotationalEnergy
        (
            const scalar temperature,
            const scalar rotationalDegreesOfFreedom
        );

        //- Generate a random vibrational energy
        labelList equipartitionVibrationalEnergyLevel
        (
            const scalar temperature,
            const label nVibrationalModes,
            const label typeId
        );

        //- Generate a random electronic energy level
        label equipartitionElectronicLevel
        (
            const scalar temperature,
            const labelList& electronicDegeneracyList,
            const scalarList& electronicEnergyList
        );

        // - Generate a LB post-collision rotational energy
        scalar postCollisionRotationalEnergy
        (
            const scalar rotationalDof,
            const scalar ChiB
        );

        // - Generate a LB post-collision vibrational energy
        label postCollisionVibrationalEnergyLevel
        (
            bool postReaction,
            const label vibLevel,
            const label iMax,
            const scalar thetaV,
            const scalar thetaD,
            const scalar refTempZv,
            const scalar omega,
            const scalar Zref,
            const scalar Ec,
            const scalar fixedZv = 0,
            const label invZvFormulation = 2,
            const label cellI = 0
        );

        // - Generate a LB post-collision electronic energy level
        label postCollisionElectronicEnergyLevel
        (
            const scalar Ec,
            const label jMax,
            const scalar omega,
            const scalarList& EElist,
            const labelList& gList
        );

        // From the Maxwellian distribution:
        //- Average particle speed
        inline scalar maxwellianAverageSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline scalarField maxwellianAverageSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;

        //- RMS particle speed
        inline scalar maxwellianRMSSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline scalarField maxwellianRMSSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;

        //- Most probable speed
        inline scalar maxwellianMostProbableSpeed
        (
            const scalar temperature,
            const scalar mass
        ) const;

        inline scalarField maxwellianMostProbableSpeed
        (
            const scalarField& temperature,
            const scalar mass
        ) const;

    // Sub-models

        //- Return reference to binary elastic collision model
        inline const BinaryCollisionModel&  binaryCollision() const;

        //- Return non-const reference to binary elastic collision model
        inline BinaryCollisionModel& binaryCollision();

        inline const dsmcBoundaries& boundaries() const;
        inline dsmcBoundaries& boundaries();

        inline const dsmcFieldProperties& fields() const;
        inline dsmcFieldProperties& fields();

        inline const dsmcControllers& controllers() const;
        inline dsmcControllers& controllers();

        //- Return reference to coordinate system model
        inline const dsmcCoordinateSystem& coordSystem() const;

        //- Return non-const reference to coordinate system model
        inline dsmcCoordinateSystem& coordSystem();

        //- Return reference to porous measurements model
        inline const porousMeasurements& porousMeas() const;

        //- Return non-const reference to porous measurements model
        inline porousMeasurements& porousMeas();

        inline const dsmcDynamicLoadBalancing& dynamicLoadBalancing() const;
        inline dsmcDynamicLoadBalancing& dynamicLoadBalancing();

        inline const dsmcFaceTracker& tracker() const;
        inline dsmcFaceTracker& tracker();

        inline const dsmcReactions& reactions() const;
        inline dsmcReactions& reactions();

        inline const boundaryMeasurements& boundaryFluxMeasurements() const;
        inline boundaryMeasurements& boundaryFluxMeasurements();

        inline const cellMeasurements& cellPropMeasurements() const;
        inline cellMeasurements& cellPropMeasurements();

      // Check

            //- Calculate properties of system
            inline scalarList infoMeasurements();

            //- Print cloud information
            void info();

             //- Check if load balance is required
            void loadBalanceCheck();

            //- Perform load balance
            void loadBalance(const label noRefinement = 0);

	    // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            void autoMap(const mapPolyMesh&);

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;

            void insertParcelInCellOccupancy(dsmcParcel* p);

        void removeParcelFromCellOccupancy
        (
            const label& cellMolId,
            const label& cell
        );

        // Clone a free parcel at the requested weight and insert it into the cloud and cell occupancy.
        dsmcParcel* cloneFreeParcelInCell
        (
            const dsmcParcel& parentParcel,
            const scalar newRWF
        );

        //- 在碰撞开始前，对某个 cell 内的 parcel 权重做“共同基准权重精确拉平”。
        // 物理思路：
        // 1. 先找出该 cell 内最小的 parcel 代表分子数；
        // 2. 再尝试把其它 parcel 的代表分子数表示成
        //    $$
        //    N_{\mathrm{real}}(p) \approx m_p \, N_{\mathrm{base}}
        //    $$
        //    其中 $m_p$ 为整数，$N_{\mathrm{base}} = N_{\min}/d$，
        //    而 $d$ 是不超过 maxAutoEqualiseDenominator 的正整数；
        // 3. 一旦找到这样的共同基准权重，就把每个 parcel 精确拆成若干同态
        //    等权重子粒子；
        // 4. 这样后续碰撞步骤就只会看到等权重 parcel，从而仍可复用现有
        //    BinaryCollisionModel::collide(parcelP, parcelQ) 接口而不丢守恒。
        //
        // 返回值：
        // - 本次为了拉平该 cell 权重而新增的 clone 数。
        label equaliseCellParcelWeightsForCollisions
        (
            const label cellI,
            const scalar integerRatioTolerance,
            const label maxMultiplicityPerParcel,
            const label maxAutoEqualiseDenominator
        );

            void reBuildCellOccupancy();

            //- 处理 DSMC 边界模型在 move() 内挂起的延迟并行粒子转移。
            // 该接口会被 Cloud<dsmcParcel>::move() 的 transfer while-loop
            // 调用，用于让远端 stage-interface 注入在同一时间步内完成收发。
            virtual bool performDeferredParallelTransfers();

        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar RWF,
                const scalar ERot,
                const label ELevel,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const label typeId,
                const label newParcel,
                const label classification,
                const labelList& vibLevel
            );

            //- Add new stuck parcel
            // used to clone stuck particles during evolution of coordinate
            // system (axisymmetric, spherical)
            void addNewStuckParcel
            (
                const vector& position,
                const vector& U,
                const scalar RWF,
                const scalar ERot,
                const label ELevel,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const label typeId,
                const label newParcel,
                const label classification,
                const labelList& vibLevel,
                // stuck properties of the parent parcel:
                // wallTemperature:
                //   [0]: Temperature
                //   [1]: Patch index
                //   [2]: Face index
                //   [3]: Pre-interaction energy (note: this is split across all
                //        splitted parcels)
                const scalarField& wallTemperature,
                // wallVectors:
                //   [0]: Wall unit tangential vector 1
                //   [1]: Wall unit tangential vector 2
                //   [2]: Wall unit normal vector
                //   [3]: Pre-interaction momentum (note: this is split across
                //        all splitted parcels)
                const vectorField& wallVectors
            );

            //- redistribute vibrational energy after a reaction
            scalar energyRatio
            (
                scalar ChiA,
                scalar ChiB
            );

            //- redistribute vibrational energy after a reaction
            scalar PSIm
            (
                scalar DOFm,
                scalar DOFtot
            );

            //- Evolve the cloud (move, collide, compute)
            void evolve();

            //- Evolve the cloud (move, collide)
            void evolve_moveAndCollide();

            //- Evolve the cloud (compute)
            void evolve_fields();

            // Return how often to write to terminal
            // (used for the solver to access this)
            label nTerminalOutputs();

            // NEW DANIEL

            void resetHybridTraRotVib
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& TrInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& DInitial,
                PtrList<volVectorField>& qtInitial,
                PtrList<volVectorField>& qrInitial,
                PtrList<volVectorField>& qvInitial,
                PtrList<volTensorField>& tauInitial,
                word& typeOfReset,
                wordList& zonesToReset
            );

            void resetHybridWhenUpdated
            (
                volVectorField& UInitial,
                PtrList<volScalarField>& TtInitial,
                PtrList<volScalarField>& TrInitial,
                PtrList<volScalarField>& TvInitial,
                PtrList<volScalarField>& numberDensitiesField,
                PtrList<volVectorField>& DInitial,
                PtrList<volVectorField>& qtInitial,
                PtrList<volVectorField>& qrInitial,
                PtrList<volVectorField>& qvInitial,
                PtrList<volTensorField>& tauInitial,
                word& typeOfReset,
                word& zoneToReset
            );

            //- Function to help with steady normal shock wave simulation
            void shockReset();
            // END NEW DANIEL

            //- Clear the Cloud
            void clear();

          // I-O

            //- Update properties from given dictionary
            bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "dsmcCloud.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
