# Protein–ligand notebook: Colab environment failure

This note is for people running `Protein_ligand.ipynb` on Google Colab who hit an error in **Install Conda Colab** or **Install dependencies**.

## What you see

- `FileNotFoundError: Expected Python binary at /content/.pixi/envs/default/bin/python` while installing Conda Colab. The same failure is reported in [issue 160](https://github.com/pablo-arantes/making-it-rain/issues/160).
- `ModuleNotFoundError: No module named 'pytraj'`, or the same error for `openmm` or `rdkit`, after the dependency cell. See [issue 159](https://github.com/pablo-arantes/making-it-rain/issues/159).
- A kernel crash, or `numpy.core.multiarray failed to import`, when trajectory analysis starts.
- On a Jupyter session outside Colab, including Windows: `ModuleNotFoundError: No module named 'condacolab'`. CondaColab only starts inside Google Colab.

## Why the published notebook fails

Current Google Colab (runtime 2026.07) is Python 3.12 with NumPy 2.0.2. The notebook installs CondaColab 0.2, which creates one Pixi environment and restarts the kernel. The next cell then mixes `pixi add` with system `pip` and does not stop when a command fails.

- `pip install numpy==2.0.0`, `rdkit==2025.3.1`, `prolif==1.1.0`, `openbabel-wheel`, and `MDAnalysis==2.8.0` install into whichever `pip` is first on `PATH`. That is often not the Pixi kernel that later imports those packages.
- AmberTools, OpenMM, PDBFixer, ParmEd, and MDTraj are added with separate `pixi add` calls, so shared libraries are not resolved in one step.
- `prolif==1.1.0` (November 2022) predates Python 3.12 and NumPy 2. The LigPlot cells call `LigNetwork.from_ifp` and `to_dataframe(return_atoms=True)`. ProLIF 2 removed both.
- conda-forge AmberTools 26 already provides `pytraj`, `tleap`, `antechamber`, `parmchk2`, `pdb4amber`, and ParmEd, and its Python 3.12 builds use NumPy 2. Installing another NumPy with pip, or pinning NumPy 1.26, conflicts with that build.

## What the updated notebook does

- One CondaColab install from conda-forge for AmberTools, OpenMM, PDBFixer, RDKit, OpenBabel, MDTraj, MDAnalysis, and seaborn. It uses Colab's own Python and does not pin NumPy.
- If Pixi writes its Python binary somewhere other than the path CondaColab expects, the setup cell links that binary into place before the kernel restarts.
- After the restart, biopandas, py3Dmol, and current ProLIF are installed with that kernel's pip. A failed install stops the cell.
- The two LigPlot cells call `LigNetwork.from_fingerprint`.

## How to run it

1. Open the notebook on a Colab GPU runtime.
2. Run **Install Conda Colab** by itself and wait until the session reconnects. Do not use Run all for that cell.
3. Run **Install dependencies**.
4. Before a production simulation, confirm `import pytraj`, `tleap -h`, and `from openbabel import pybel`.
