# neutrinoPpackages

Some neutrino packages written in Python and Mathematica.

## Build docs

To build the documentation, sphinx-doc is required.

```
sphinx-build -b html source build
```

will do the work, where `source` is the path for the source of the documentation and `build` is for the target path of the generated html files.

## Code of Conduct

1. For Mathematica code, parameters should be arranged in the following way
   ```
   list of n's, list of perturbation k's, list of perturbation amplitudes a's, list of perturbation phases phi's, mixing angle in matter thetam, endpoint of numerical calculation
   ```
