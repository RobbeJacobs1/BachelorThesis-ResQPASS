# Bachelor Thesis – ResQPASS

This repository contains the code accompanying a bachelor thesis on the **ResQPASS algorithm**, focusing on the effect of the characteristic mesh length on convergence behaviour.

---

## Repository structure

The main MATLAB code is located in the `MatLab/` directory.

- `MatLab/`  
  Core MATLAB implementation and experiment scripts.

- `MatLab/Create_2D_meshes/`  
  Scripts for generating and processing 2D mesh-based problems.

- `MatLab/Create_3D_meshes/`  
  Scripts for generating and processing 3D mesh-based problems.

---

## Mesh preprocessing

The folders `Create_2D_meshes` and `Create_3D_meshes` contain scripts to:
- generate meshes
- parse `.msh` files using Python
- prepare data for MATLAB simulations

Python scripts used:
- `gmsh_parser2D.py`
- `gmsh_parser3D.py`

### Requirements
To use these scripts, **Gmsh must be installed and available in your system PATH**.

---

## Figures and plotting

All figures in the thesis are generated from MATLAB scripts.

- Many scripts include a `saveFigure` flag to export figures.
- Formatting is handled via `setFigParameters.m`.
- Figures are stored in `MatLab/Images/`.

---

## Note on data

Large generated data files (e.g. `.mat`, `.msh`, simulation results) are not included in this repository.

These files are required for:
- `MatLab/subsection_5_2_2_preConditioning.m`
- `MatLab/section_6_2_springProblem2D.m`
- `MatLab/section_6_3_springProblem3D.m`

### Data setup

To reconstruct the required data, either:

- Regenerate it using the provided mesh generation and preprocessing scripts:
   - `MatLab/Create_2D_meshes/`
   - `MatLab/Create_3D_meshes/`
   - Uncomment the solvers in `MatLab/subsection_5_2_2_preConditioning.m` and comment the load command.

or

- Download the precomputed `Results/` directory from ... and place it in the repository root    

    Afterwards, run:
    ```bash
    python setup_project.py

---

## Dependencies

- MATLAB (tested version: [R2025b])
- Python (for mesh parsing scripts, tested on 3.12.3)
- Gmsh