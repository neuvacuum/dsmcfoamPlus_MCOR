# MCOR1_0 Tests

English | [简体中文](README.zh-CN.md)

From the package root, load OpenFOAM and this package's environment, then run:

```bash
./tests/check-installation.sh
./tests/run-smoke.sh
```

The installation check confirms that both applications start and that all three extension libraries load from this package. The closed-box smoke test checks 10 time steps, 20 regulation-stage records, nonzero collisions and pruning operations, and a maximum scaled moment residual below `1e-9`.

Use the commands in the [case guide](../cases/README.md) for the three short case tests. The tested execution coverage is:

| Case | Configurations | Mesh cells | Processes | Steps |
| --- | --- | ---: | ---: | ---: |
| Couette | Fixed weight, MCOR exact collisions, MCOR stochastic updates | 500 | 1 | 10 |
| Plate | Fixed weight, MCOR exact collisions | 6000 | 2 | 5 |
| Orifice | Fixed weight, MCOR exact collisions | 135168 | 2 | 1 |

These short tests completed mesh checks, initialization, and solving; the MPI tests also completed final-time reconstruction. Acceptance-clipping, invalid-candidate, and candidate-cap counts were zero, and moment residuals were below `1e-9`. These tests check basic execution paths. Statistical accuracy, discretization parameters, and physical models require assessment for each research problem.
