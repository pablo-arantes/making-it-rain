## Making it rain
Cloud-based molecular simulations for everyone
## Hello there!

Welcome to Making it rain page, here you can find a Jupyter notebook scripts for running Molecular Dynamics (MD) simulations using OpenMM engine and AMBER and CHARMM force fields files on Google Colab. This site is a supplementary material of the paper "***Making it rain: Cloud-based molecular simulations for everyone***" and we encourage you to read it before using this pipeline.

![alt text](GraphAbs.png)

The main goal of this work is to demonstrate how to harness the power of cloud-computing to run microsecond-long MD simulations in a cheap and yet feasible fashion.

**Important**: We've updated the notebooks to [CondaColab](https://github.com/conda-incubator/condacolab). Now, all the dependencies will be installed faster than before (less than half of the previous time). You will see a CondaColab cell, just run and wait a few seconds, the session will restart and this is **normal and expected**. After that, you can continue running the cells like normal. Do not use the `Run all` option. Run the `condacolab` cell _individually_ and wait for the kernel to restart. **Only then**, you can run all cells if you want. Thank you for your support.

1. **AMBER** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/Making-it-rain/blob/main/Amber.ipynb)  - `Using AMBER to generate topology and to build the simulation box`
2. **CHARMM** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/Making-it-rain/blob/main/CHARMM_GUI.ipynb) - `Using inputs from CHARMM-GUI solution builder`
3. **AlphaFold2+MD** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/Making-it-rain/blob/main/AlphaFold2%2BMD.ipynb) - `Using AlphaFold2_mmseqs2 to generate protein model + MD simulation using AMBER to generate topology and to build simulation box`


**UPDATE (October 2021)**

4. **Protein-Ligand simulations** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Protein_ligand.ipynb)  - `Using AMBER to generate topology and to build the simulation box and for the ligand using GAFF2 force field`
5. **Using AMBER Inputs** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Amber_inputs.ipynb)  - `Using inputs from AMBER suite of biomolecular simulation program`
6. **Using GROMACS Inputs** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Gromacs_inputs.ipynb)  - `Using inputs from GROMACS biomolecular simulation package (AMBER, CHARMM and OPLS force fields are compatible)`

**UPDATE (March 2022)**

7. **RESP Partial Charges** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Partial_Charges.ipynb)  - `Using a SMILES as input and outputs a mol2 file with RESP derived partial charges. Options for setting method (HF, B3LYP, ...), basis set (3-21G, 6-31G*) and singlepoint or geometry optimization are available`
8. **Small Molecules MD** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/MD_Small_Molecules.ipynb)  - `Using a SMILES as a input, geometry optimization with TorchANI and topology with AMBER (GAFF2 force field)`
9. **GLYCAM** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Glycam.ipynb)  - `Using inputs from GLYCAM server`

**UPDATE (August 2022)**

10. **DRUDE** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Drude.ipynb)  - `Using inputs from CHARMM-GUI Drude Prepper`

**UPDATE (March 2024)**

11. **Protein-Membrane simulations** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Protein%2BMembranes.ipynb)  - `Using OpenFF to generate the topology and build the simulation box for protein-membrane systems with AMBER force fields.`
12. **Martini+cg2all** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Martini%2Bcg2all.ipynb)  - `Utilizing Vermouth, the Python library that powers Martinize2, to generate the topology and build the simulation box for protein systems using Martini force fields. Additionally, employing cg2at to predict all-atom trajectories from coarse-grained (CG) representations.`
13. **AMBER Mutations** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/Amber_mutations.ipynb)  - `Performing mutations on protein/nucleic acid systems and utilizing AMBER to generate the topology and build the simulation box.`

**UPDATE (March 2025)**

First, we made MD simulations rain down from the cloud. Now, we’re bringing deep learning into the mix! 

14. **Subsampled AlphaFold2** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GMdSilva/gms_natcomms_1705932980_data/blob/main/AlphaFold2_Traj_v1.ipynb) - `Colab notebook for running the subsampled AlphaFold2 approach for predicting protein conformational ensembles.`
15. **Biomolecular Emulator** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pablo-arantes/making-it-rain/blob/main/BioEmu.ipynb) - `Biomolecular Emulator (BioEmu), a model that samples from the approximated equilibrium distribution of structures for a protein monomer, given its amino acid sequence.`


Want to try MD simulation on NAMD using Google Colab? [**Mostafa Sayed**](https://github.com/mabdelmaksoud53) create a colab notebook for running MD simulations using NAMD and CHARMM-GUI inputs. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mabdelmaksoud53/Colab_NAMD_Suite/blob/main/MD_NAMD.ipynb)

## Bugs
- If you encounter any bugs, please report the issue [here](https://github.com/pablo-arantes/Making-it-rain/issues).

## Acknowledgments

- We would like to thank the [Psi4](https://psicode.org/) team for developing an excellent and open source suite of ab initio quantum chemistry.
- We would like to thank the [Roitberg](https://roitberg.chem.ufl.edu/) team for developing the fantastic [TorchANI](https://github.com/aiqm/torchani).
- We would like to thank the OpenMM team for developing an excellent and open source engine. 
- We would like to thank the AlphaFold team for developing an excellent model and open sourcing the software. 
- We would like to thank the ChemosimLab ([@ChemosimLab](https://twitter.com/ChemosimLab)) team for their incredible [ProLIF](https://prolif.readthedocs.io/en/latest/index.html#) (Protein-Ligand Interaction Fingerprints) tool.
- Credit to Sergey Ovchinnikov ([@sokrypton](https://twitter.com/sokrypton)), Milot Mirdita ([@milot_mirdita](https://twitter.com/milot_mirdita)) and Martin Steinegger ([@thesteinegger](https://twitter.com/thesteinegger)) for their fantastic [ColabFold](https://github.com/sokrypton/ColabFold).
- Making it rain by **Pablo R. Arantes** ([@pablitoarantes](https://twitter.com/pablitoarantes)), **Marcelo D. Polêto** ([@mdpoleto](https://twitter.com/mdpoleto)), **Conrado Pedebos** ([@ConradoPedebos](https://twitter.com/ConradoPedebos)) and **Rodrigo Ligabue-Braun** ([@ligabue_braun](https://twitter.com/ligabue_braun)).
- Also, credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin.
- Finally, we would like to thank Professor [Giulia Palermo](https://palermolab.com/) for her support and thoughtful comments in the development of the present work.

## Do you want to cite this work?


Arantes P.R., Depólo Polêto M., Pedebos C., Ligabue-Braun R. Making it rain: cloud-based molecular simulations for everyone. 
Journal of Chemical Information and Modeling 2021. DOI: [10.1021/acs.jcim.1c00998](https://doi.org/10.1021/acs.jcim.1c00998).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5196783.svg)](https://doi.org/10.5281/zenodo.5196783)

