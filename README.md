# C-MCLP

Research code for experimenting with a capacitated or coverage-constrained variant of the Maximum Covering Location Problem (MCLP) on geospatial data. The repository combines three main pieces:

- geospatial preprocessing of regions and census-style demand data
- candidate facility-site generation over polygons
- staged optimization with Gurobi for MCLP and C-MCLP formulations

This is not a packaged library or polished application. It is a legacy prototype codebase with notebooks, solver logs, shapefiles, and Python modules used for experiments.

## What Is In The Repo

### Core modules

- `alg/cmclp_solver.py`
  - Builds and solves staged `MCLP` and `C-MCLP` models with `gurobipy`
  - Generates demand-to-facility covering sets `N`
  - Includes helper routines for evaluating actual coverage across stages

- `alg/census_parser.py`
  - Loads place, face, and block-level shapefiles with `fiona`
  - Filters geometries by place and county identifiers
  - Transforms geometries to a projected CRS for meter-based spatial calculations
  - Generates randomized demand points inside populated block groups

- `alg/data_generator.py`
  - Utilities for random demand-point generation and population assignment
  - Generates triangular candidate-site grids over a bounding box

- `alg/marching_army.py`
  - Alternative candidate-site generation approach based on marching lines through a polygon

- `alg/result_analyzer.py`
  - Post-processing helpers for coverage circles, overlap calculations, and plotting

- `alg/distance_buffer.py`
  - Low-level geometric helpers, plotting helpers, and distance calculations

- `alg/border_generators.py`
  - Region/boundary loading helpers, including an `osmnx`-based place fetcher

### Data and experiments

- `data/`
  - Contains sample shapefiles and dummy region inputs
  - `data/YukonAlaska/` includes region shapefiles used by the prototype
  - `data/DummyData1/` contains simple text-based boundary inputs

- `dev/ipynbs/`
  - Exploratory notebooks for buffer-based, marching-army, and solving-method experiments
  - Also contains historical solver logs and generated figures

- `dev/PolyScipDev/`
  - Older modeling experiments related to PolySCIP/ZIMPL workflows

## Expected Workflow

At a high level, the code supports a workflow like this:

1. Load a region and associated census/block demand data.
2. Generate demand points and weights from block populations.
3. Generate candidate facility points with either a triangular grid or marching-army method.
4. Build coverage sets based on service radius.
5. Solve staged `MCLP` or `C-MCLP` models with Gurobi.
6. Analyze coverage, overlap, and stage-by-stage results.

There is no single CLI entry point. The repository is organized around direct module use plus notebooks in `dev/ipynbs/`.

## Dependencies

The code imports the following Python packages:

- `gurobipy`
- `shapely`
- `geopandas`
- `fiona`
- `pyproj`
- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `rtree`
- `descartes`
- `osmnx`

System-level geospatial dependencies may also be required for packages such as `fiona`, `geopandas`, and `rtree`.

## Environment Notes

- The codebase uses legacy Python syntax in multiple files, including Python 2 style `print` statements.
- Some APIs used here are also older geospatial-library interfaces.
- `alg/cmclp_solver.py` requires a working Gurobi installation and license.
- Import paths are managed with `sys.path.append(...)` in several modules rather than package-relative imports.
- The repository currently has no `requirements.txt`, `pyproject.toml`, or automated test suite.

In practice, this project should be treated as an experimental codebase that may need modernization before use in a clean Python 3 environment.

## Minimal Example

There is no guaranteed end-to-end script in the repository, but the solver modules are structured for notebook-driven usage similar to:

```python
from alg import cmclp_solver as solver

# demandpts: list[(x, y)]
# demandwts: list[int | float]
# facilitypts: list[(x, y)]

stage_results = solver.cmclp_stages(
    demandpts=demandpts,
    demandwts=demandwts,
    facilitypts=facilitypts,
    stageP=5,
    stageS=1000.0,
)
```

For region and demand generation, inspect:

- `alg/census_parser.py`
- `alg/data_generator.py`
- notebooks under `dev/ipynbs/`

## Repository Layout

```text
C-MCLP/
├── alg/                 # Core geospatial and optimization code
├── data/                # Sample shapefiles and dummy inputs
├── dev/ipynbs/          # Research notebooks, logs, and figures
├── dev/PolyScipDev/     # Older modeling experiments
└── README.md
```

## Limitations

- No installation or packaging metadata
- No formal tests or reproducible run scripts
- Legacy code style and mixed experimentation artifacts
- Heavy reliance on notebook exploration
- Solver behavior depends on commercial Gurobi tooling

## Recommended Next Steps

If you plan to keep developing this repository, the highest-value improvements would be:

1. Add a reproducible Python environment definition.
2. Port the code cleanly to Python 3.
3. Add a small end-to-end example script.
4. Separate experimental notebooks from stable library code.
5. Add tests around geometry generation and solver-stage logic.
