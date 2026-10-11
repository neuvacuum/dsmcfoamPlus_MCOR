# Extension Features and Source Mapping

English | [简体中文](implementation.zh-CN.md)

## 1. Research Configuration

MCOR1_0 provides the following extensions within dsmcFoam+. The research configuration uses a nonrotating fixed mesh and a neutral, nonreacting gas. Other model combinations require targeted validation.

| Feature | Implementation relative to `src/lagrangian/dsmc` | Behavior |
| --- | --- | --- |
| Moment-conservative occupancy regulation | `momentConservativeOccupancy` in `coordinateSystem/derived/Cartesian/dsmcCartesian.C` | Splits underpopulated groups and prunes overpopulated groups subject to nonnegative weights |
| Weighted collisions | `collisionPartnerSelection/derived/noTimeCounter/noTimeCounter.C` | Handles unequal-weight collisions, minimum-weight participation, and residual clones using parcel event weights |
| Weighted macroscopic statistics | `macroscopicProperties/derived/combined/dsmcVolFields/dsmcVolFields.C` | Accumulates statistics using the number of real molecules represented by each parcel |
| Free-stream injection | `boundaries/derived/generalBoundaries/dsmcFreeStreamInflowPatch/dsmcFreeStreamInflowPatch.C` | Uses a common reference weight for Cartesian injection, independently of local target weights |
| Integration and data support | `clouds/dsmcCloud.C`, `parcels/dsmcParcel.C`, and related files | Provides movement/collision-stage calls, occupancy lists, and weight data |

Moment-conservative pruning groups parcels by species and classification, with at most 128 parcels per batch. Let K denote the number of internal-state categories. Batches containing at most K+13 parcels are skipped before singular value decomposition (SVD). Constraints cover internal-state populations, three momentum components, six second-order velocity moments, rotational energy, and three total heat-flux components. Pruning retains the states of existing parcels, adjusts nonnegative weights, and deletes zero-weight parcels.

After movement, the sequence is splitting underpopulated groups, pruning, and rechecking underpopulation. After collisions, only pruning is performed. The two stages have independent invocation intervals. Regulation acts on eligible free, untracked parcels in valid cells; the research configuration uses a nonrotating fixed mesh and a neutral, nonreacting gas.

`dynamicMajorant` updates the candidate-event rate between events using the current parcel count, allowing residual clones to enter the same-step candidate pool immediately. Cell weight and collision-kernel upper bounds remain fixed throughout the cell collision stage. Newly recorded kernel maxima are used in subsequent time steps. A raw acceptance probability above 1 causes the dynamic branch to terminate with an error.

## 2. Configuration

Set `collisionPartnerSelectionModel` to `noTimeCounter` in `constant/dsmcProperties`. The following excerpt contains the key research settings. Use `examples/closedBox` for a complete runnable dictionary.

```foam
collisionPartnerSelectionModel noTimeCounter;
adaptiveProperties
{
    active true;
    timeStepAdaptation false;
    subCellAdaptation false;
    cellWeightAdaptation true;
    cellWeightAdaptationMode runtimeMomentConservative;
    targetParticlesPerSubCell 30;
    assumedSubCellsPerCell 1;
    occupancyTargetTolerance 0.05;
    lightweightCellWeightingInterval 1;
    postCollisionOccupancyCleanup true;
    occupancyLimitMergeInterval 1;
    maxParticlesPerCellFactor 1;
    collisionWeighting true;
    collisionCandidateTimeMode dynamicMajorant;
    collisionCandidateSelectionMode uniformCell;
    collisionUnequalWeightHandlingMode partialSplitToMinWeight;
    collisionResidualCloneHandlingMode exactClone;
    collisionResidualClonesParticipateSameStep true;
    collisionAcceptanceSafetyFactor 1;
    frozenGlobalMassCorrection false;
    runtimeRemapMomentumEnergyCorrection false;
    collisionWeightAutoEqualise false;
    runtimeEmptyCellRefill false;
}
```

- The product of `targetParticlesPerSubCell` and `assumedSubCellsPerCell` defines the cell target. The latter is a target conversion factor, distinct from a geometric collision-subcell count.
- `occupancyTargetTolerance` defines a tolerance band around the target; `maxParticlesPerCellFactor` also affects the pruning trigger. Attainable occupancy depends on constraint feasibility. Empty cells remain empty until populated through the gas-transport or boundary mechanisms.
- `uniformCell` selects pairs across the entire cell.
- Applying MCOR to adaptive mesh refinement, rotating frames, reactions, adsorption, or other model combinations requires separate validation.
- The two-species, two-cell example is a short serial installation check. Set research parameters according to the physical problem and accuracy requirements.

## 3. Fixed-Weight Control

The fixed-weight control in the paper uses this package's modified solver. It shares the collision kernel, dynamic candidate clock, within-cell pairing, and weighted statistics with MCOR; occupancy regulation is the distinguishing option. Start the control from a consistent fixed-weight initial cloud, prepared separately from a regulated variable-weight cloud.

To disable regulation, explicitly set `cellWeightAdaptation false` and `postCollisionOccupancyCleanup false`, while retaining common collision options such as `collisionWeighting true`. Check independent options individually when changing `active`, and verify the full configuration against each case's inputs.

Runtime and accuracy comparisons within this shared framework assess the effect of enabling MCOR. Mesh-generation dictionaries, physical and numerical inputs, and execution instructions for the three cases are provided in the [case guide](../cases/README.md).
