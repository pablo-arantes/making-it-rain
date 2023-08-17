'''
Original file: https://colab.research.google.com/drive/1OAF63N47PNpxuVR12RSqxEuznS6IjDM9#scrollTo=H-oBjHKEBPJY
Conda env: conda create -c conda-forge -n mmgbsa python=3.7 openbabel ambertools parmed rdkit openmmforcefields mdanalysis
Installation:
'''

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw


from openmm import *
from openmm.app import *
from openmm.unit import *
import pytraj as pt

from openbabel import pybel

import os, fnmatch
import sys
from sys import stdout
import subprocess
from copy import copy
import argparse
import locale


AllChem.SetPreferCoordGen(True)

parser=argparse.ArgumentParser()
# Cell 1
#@markdown **Important:** The protonation of your ligand is crucial for the correct parameterization of the molecule.

parser.add_argument('--workDir',type=str,
                    default='./Output/',
                    help='')
parser.add_argument('--protein_file',type=str,
                    default='./Sample_Data/protein.pdb',
                    help='pdb file')
parser.add_argument('--ligand_file',type=str,
                    default='./Sample_Data/ligand.pdb',
                    help='pdb file')
parser.add_argument('--remove_waters',choices=['yes', 'no'],
                    default='yes',
                    help='')
parser.add_argument('--Add_ligand_hydrogens',choices=['Yes', 'No'],
                    default='Yes',
                    help='')

# Cell 3
#@title **Parameters to generate the topology:**
parser.add_argument('--Force_field',choices=["ff19SB", "ff14SB"],
                    default='ff19SB',
                    help='')
parser.add_argument('--Water_type',choices=["TIP3P", "OPC"],
                    default='TIP3P',
                    help='')
parser.add_argument('--size_box',type=int,
                    default=12,
                    help='Angstrons; type:"slider", min:10, max:20, step:1')
parser.add_argument('--Ions',choices=["NaCl", "KCl" ],
                    default="NaCl",
                    help='')
parser.add_argument('--Concentration',type=str,
                    default="0.15",
                    help='type:"str"')
parser.add_argument('--Ligand_Force_field',type=str,
                    default="GAFF2",
                    help='type:"str"')
parser.add_argument('--Ligand_isomer',type=str,
                    default="1",
                    help='type:"string", min:1, max:10, step:100')

# Cell 4
#@title ### **Parameters for MD Equilibration protocol:**
parser.add_argument('--Jobname',type=str,
                    default="prot_lig_equil",
                    help='type:"string"')
parser.add_argument('--Minimization_steps',type=str,
                    default="1000",
                    help='type:"string"; #@param "1000", "5000", "10000", "20000", "50000", "100000"')
parser.add_argument('--Time',type=str,
                    default="5",
                    help='Simulation time (in nanoseconds)')
parser.add_argument('--Integration_timestep',type=str,
                    default="0.5",
                    help='integration time (in femtoseconds)"; #@param ["0.5", "1", "2", "3", "4"]')
parser.add_argument('--Temperature',type=int,
                    default=298,
                    help='Temperature (in Kelvin)')
parser.add_argument('--Pressure',type=int,
                    default=1,
                    help='Pressure (in bar)')
parser.add_argument('--Force_constant',type=int,
                    default=700,
                    help='#@param {type:"slider", min:0, max:2000, step:100}')
parser.add_argument('--Write_the_trajectory',type=str,
                    default="10",
                    help='#@param ["10", "100", "200", "500", "1000"]')
parser.add_argument('--Write_the_log',type=str,
                    default="10",
                    help='#@param ["10", "100", "200", "500", "1000"]')

# Cell 6
#@markdown ### **Parameters for MD Production protocol:**
parser.add_argument('--Equilibrated_PDB',type=str,
                    default="prot_lig_equil.pdb",
                    help='')
parser.add_argument('--State_file',type=str,
                    default="prot_lig_equil.rst",
                    help='')
parser.add_argument('--Stride_Time',type=str,
                    default="10",
                    help='')
parser.add_argument('--Number_of_strides',type=str,
                    default="1",
                    help='')

# Cell 8
parser.add_argument('--Skip',type=str,
                    default="1",
                    help='#@param ["1", "2", "5", "10", "20", "50"]')
parser.add_argument('--Output_format',type=str,
                    default="dcd",
                    help='#@param ["dcd", "pdb", "trr", "xtc"]')
parser.add_argument('--first_stride',type=str,
                    default="1",
                    help='')
parser.add_argument('--trajectory_saved_frequency',type=str,
                    default="10",
                    help='#@param ["10", "100", "200", "500", "1000"]')
# Cell 10
parser.add_argument('--igb',choices=["1", "2", "5", "7", "8"],
                    default="2",
                    help='the "OBC" models (igb=2 and 5) are newer, but appear to give significant improvements and are recommended for most projects (For more information check the Section 4.1 of the [Amber Manual](https://ambermd.org/doc12/Amber20.pdf))')

options = parser.parse_args()

# Cell 1
#@markdown **Important:** The protonation of your ligand is crucial for the correct parameterization of the molecule.

Protein_PDB_file_name = options.protein_file
remove_waters = options.remove_waters
if remove_waters == "yes":
  no_waters = "nowat"
else:
  no_waters = ''

Add_ligand_hydrogens = options.Add_ligand_hydrogens
ligand_name = options.ligand_file
workDir = options.workDir

initial_pdb = os.path.join(workDir, str(Protein_PDB_file_name))
prepareforleap = os.path.join(workDir, "prepareforleap.in")
ligand_pdb = os.path.join(workDir, str(ligand_name))
ligand_pdb2 = os.path.join(workDir, "ligand_H.pdb")
starting = os.path.join(workDir, "starting1.pdb")
starting2 = os.path.join(workDir, "starting2.pdb")
starting_end = os.path.join(workDir, "starting_end.pdb")


def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def remove_lines(filename):
    with open(filename, 'r') as file:
        ter_count = 0
        for line in file:
            if line.startswith('TER'):
                ter_count += 1
                if ter_count >= 1:
                    yield line
                    for i in range(3):
                        line = next(file, None)
                        if line is not None and line.startswith('ATOM') and line.split()[2] in ['P', 'OP1', 'OP2']:
                            continue
                        else:
                            yield line
                else:
                    yield line
            else:
                yield line

if Add_ligand_hydrogens == "Yes":
  mol= [m for m in pybel.readfile(filename=ligand_pdb, format='pdb')][0]
  out=pybel.Outputfile(filename="temp.mol",format='mol',overwrite=True)
  out.write(mol)
  out.close()

  mol = Chem.MolFromMolFile('temp.mol', removeHs=True)
  hmol = Chem.AddHs(mol)
  mp = AllChem.MMFFGetMoleculeProperties(hmol)
  ff = AllChem.MMFFGetMoleculeForceField(hmol, mp)
  for a in hmol.GetAtoms():
    if (a.GetAtomicNum() > 1):
      ff.MMFFAddPositionConstraint(a.GetIdx(), 0, 1.e4)
  ff.Minimize(maxIts=1000)
  charge_mol = Chem.rdPartialCharges.ComputeGasteigerCharges(hmol)
  charge = Chem.GetFormalCharge(hmol)
  print("Charge = " + str(charge))
  # AllChem.MolToMolFile(hmol, (os.path.join(workDir, f"start_min.mol")))
  AllChem.MolToPDBFile(hmol, ligand_pdb2)
  mol_end = mol_with_atom_index(hmol)
else:
  mol= [m for m in pybel.readfile(filename=ligand_pdb, format='pdb')][0]
  out=pybel.Outputfile(filename="temp.mol",format='mol',overwrite=True)
  out.write(mol)
  out.close()

  hmol = Chem.MolFromMolFile('temp.mol', removeHs=False)
  mp = AllChem.MMFFGetMoleculeProperties(hmol)
  ff = AllChem.MMFFGetMoleculeForceField(hmol, mp)
  for a in hmol.GetAtoms():
    if (a.GetAtomicNum() > 1):
      ff.MMFFAddPositionConstraint(a.GetIdx(), 0, 1.e4)
  ff.Minimize(maxIts=1000)
  charge_mol = Chem.rdPartialCharges.ComputeGasteigerCharges(hmol)
  charge = Chem.GetFormalCharge(hmol)
  print("Charge = " + str(charge))
  # AllChem.MolToMolFile(hmol, (os.path.join(workDir, f"start_min.mol")))
  AllChem.MolToPDBFile(hmol, ligand_pdb2)
  mol_end = mol_with_atom_index(hmol)
  IPythonConsole.drawMol3D(hmol)

#Fix protein
f = open(prepareforleap, "w")
f.write("""parm """ + str(initial_pdb) + "\n"
"""loadcrd """ + str(initial_pdb) + """ name edited""" + "\n"
"""prepareforleap crdset edited name from-prepareforleap \ """ + "\n"
"""pdbout """ + str(starting) + " " + str(no_waters) + """ noh""" + "\n"
"""go """)
f.close()

prepareforleap_command = "cpptraj -i " + str(prepareforleap)
original_stdout = sys.stdout # Save a reference to the original standard output
with open('prepareforleap.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(prepareforleap_command)
    sys.stdout = original_stdout # Reset the standard output to its original value

subprocess.run(["chmod 700 prepareforleap.sh"], shell=True)
subprocess.run(["./prepareforleap.sh"], shell=True,)


pdb4amber_cmd = "pdb4amber -i " + str(starting) + " -o " + str(starting_end) + " -a"
original_stdout = sys.stdout # Save a reference to the original standard output

with open('pdb4amber.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(pdb4amber_cmd)
    sys.stdout = original_stdout # Reset the standard output to its original value

subprocess.run(["chmod 700 pdb4amber.sh"], shell=True)
subprocess.run(["./pdb4amber.sh"], shell=True)


protein_check = os.path.exists(starting_end)
ligand_check = os.path.exists(ligand_pdb2)

if protein_check == True and ligand_check == True:
  print("Successfully generated protein and ligand files! :-)")
else:
  print("ERROR: Check your inputs! ")

#Cell 2
#@title **Enumerate Stereoisomers to generate ligand topology:**
##@markdown **You can find the smiles for your lingad at: https://pubchem.ncbi.nlm.nih.gov/**

mol= [m for m in pybel.readfile(filename=ligand_pdb2, format='pdb')][0]
mol.calccharges
mol.addh()
out=pybel.Outputfile(filename="temp2.smi",format='smiles',overwrite=True)
out.write(mol)
out.close()

fileObj = open("temp2.smi", "r",) #opens the file in read mode
for aRow in fileObj:
    smi = aRow.split('\t')
fileObj.close()

Ligand_smiles = smi[0]
os.system('rm temp2.smi >/dev/null 2>&1')

mol = Chem.MolFromSmiles(Ligand_smiles)

def spam(n):
    out=[]
    for perm in getPerms(n):
        elem = [ int(i) for i in list(perm) ]
        out.append(elem)
    return out

def getPerms(n):
    from itertools import permutations
    for i in getCandidates(n):
        for perm in set(permutations(i)):
            yield ''.join(perm)

def getCandidates(n):
    for i in range(0, n+1):
        res = "1" * i + "0" * (n - i)
        yield res

def GetStereoIsomers(mol):
    out = []

    chiralCentres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    #return the molecule object when no chiral centres where identified
    if chiralCentres == []:
        return [mol]

      #All bit permutations with number of bits equals number of chiralCentres
    elements = spam(len(chiralCentres))
    os.system('rm smiles.txt temp2.smi >/dev/null 2>&1')
    for isoId,element in enumerate(elements):
        for centreId,i in enumerate(element):
            atomId = chiralCentres[centreId][0]
            if i == 0:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif i == 1:
                mol.GetAtomWithIdx(atomId).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
        outmol = copy(mol)
        out.append(outmol)
        print(Chem.MolToSmiles(mol,isomericSmiles=True), file=open("smiles.txt", "a",))
    return out

Draw.MolsToGridImage(GetStereoIsomers(mol), subImgSize=(500,200), molsPerRow=1)
chiralCentres = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
if chiralCentres != []:
  print("Follow the stereoisomers for your ligand: \n")
  fileObj = open("smiles.txt", "r",) #opens the file in read mode
  smiles = fileObj.read().splitlines() #puts the file into an array
  fileObj.close()
  x = len(smiles[:-1])
  for a in range(x+1):
    y = smiles[0+a:(a+1)]
    globals()[f"isomer{a+1}"] = str(y[0])
    print("Isomer " + str(a+1) + " = " + str(y[0]) + "\n")
else:
  isomer1 = Ligand_smiles
  print("No chiral centres were identified! \nIsomer 1 = " + str(isomer1)  )
Draw.MolsToGridImage(GetStereoIsomers(mol), subImgSize=(700,200), molsPerRow=1, returnPNG=True)

# Cell 3
#@title **Parameters to generate the topology:**

Force_field = options.Force_field
if Force_field == "ff19SB":
  ff = "leaprc.protein.ff19SB"
else:
  ff = "leaprc.protein.ff14SB"

Water_type = options.Water_type
if Water_type == "TIP3P":
  water = "leaprc.water.tip3p"
  water_box = "TIP3PBOX"
else:
  water = "leaprc.water.opc"
  water_box = "OPCBOX"

size_box = options.size_box

Ions = options.Ions

Concentration = options.Concentration

Ligand_Force_field = options.Ligand_Force_field

Ligand_isomer = options.Ligand_isomer
if chiralCentres == []:
  isomer_end = isomer1
else:
  isomer_end = globals()[f"isomer{Ligand_isomer}"]

Ligand_net_charges = charge

tleap = os.path.join(workDir, "tleap.in")
top_nw = os.path.join(workDir, "SYS_nw.prmtop")
crd_nw = os.path.join(workDir, "SYS_nw.crd")
pdb_nw = os.path.join(workDir, "SYS_nw.pdb")
top = os.path.join(workDir, "SYS_gaff2.prmtop")
crd = os.path.join(workDir, "SYS_gaff2.crd")
pdb = os.path.join(workDir, "SYS.pdb")
ligand_noh = os.path.join(workDir, "ligand_noh.pdb")
ligand_h = os.path.join(workDir, "ligand_h.pdb")
ligand_mol2 = os.path.join(workDir, "ligand.mol2")
ligand_frcmod = os.path.join(workDir, "ligand.frcmod")
lig_new = os.path.join(workDir, "ligand_gaff.pdb")
protein_ligand = os.path.join(workDir, "protein_ligand.pdb")
lib = os.path.join(workDir, "lig.lib")

#gaff_command1 = "pdb4amber -i " + str(ligand_pdb2) + " -o " + str(ligand_h)
gaff_command1 = "pdb4amber -i " + str(ligand_pdb2) + " -o " + str(ligand_h)
gaff_command3 = "antechamber -i " + str(ligand_h) + " -fi pdb -o " + str(ligand_mol2) + " -fo mol2 -c bcc -nc " + str(Ligand_net_charges) + " -rn LIG -at gaff2"
gaff_command4 = "parmchk2 -i " + str(ligand_mol2) + " -f mol2 -o " + str(ligand_frcmod) + " -s gaff2"

original_stdout = sys.stdout # Save a reference to the original standard output

with open('gaff.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(gaff_command1)
    print(gaff_command3)
    print(gaff_command4)
    sys.stdout = original_stdout # Reset the standard output to its original value

os.system('chmod 700 gaff.sh 2>&1 1>/dev/null')
os.system('bash gaff.sh >/dev/null 2>&1')

f = open(tleap, "w")
f.write("""source """ + str(ff) + "\n"
"""source leaprc.gaff2
LIG = loadmol2 """ + str(ligand_mol2) + "\n"
"""loadamberparams """ + str(ligand_frcmod) + "\n"
"""saveoff LIG """ + str(lib) + "\n"
"""savepdb LIG """ + str(lig_new) + "\n"
"""quit""")
f.close()

tleap_command = "tleap -f " + str(tleap)
cat_command = "cat " + str(starting_end) + " " + str(lig_new) + str(" > ") + str(protein_ligand)

original_stdout = sys.stdout # Save a reference to the original standard output

with open('run_tleap.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(tleap_command)
    print(cat_command)
    sys.stdout = original_stdout # Reset the standard output to its original value

os.system('chmod 700 run_tleap.sh 2>&1 1>/dev/null')
os.system('bash run_tleap.sh 2>&1 1>/dev/null')

ppdb = PandasPdb().read_pdb(protein_ligand)
ppdb.df['ATOM'] = ppdb.df['ATOM']
ppdb.df['OTHERS'] = [ppdb.df['OTHERS'] != 'OTHERS']
ppdb.to_pdb(path=protein_ligand, records=['ATOM', 'HETATM'], gz=False, append_newline=True)

f = open(tleap, "w")
f.write("""source """ + str(ff) + "\n"
"""source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.GLYCAM_06j-1
source leaprc.gaff2
source """  + str(water) + "\n"
"""loadamberparams """ + str(ligand_frcmod) + "\n"
"""loadoff """ + str(lib) + "\n"
"""SYS = loadpdb """ + str(protein_ligand) + "\n"
"""alignaxes SYS
savepdb SYS """ + str(pdb_nw) + "\n"
"""saveamberparm SYS """ + str(top_nw) + " " + str(crd_nw) + "\n"
"""solvatebox SYS """ + str(water_box) + " " + str(size_box) +  """ 0.7
saveamberparm SYS """ + str(top) + " " + str(crd) + "\n"
"""savepdb SYS """ + str(pdb) + "\n"
"""quit""")
f.close()

tleap_command = "tleap -f " + str(tleap)

original_stdout = sys.stdout # Save a reference to the original standard output

with open('run_tleap.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(tleap_command)
    sys.stdout = original_stdout # Reset the standard output to its original value

SYS = os.path.join(workDir, "SYS*")
rm_sys = "rm " + SYS

original_stdout = sys.stdout # Save a reference to the original standard output

with open('rm_sys.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(rm_sys)
    sys.stdout = original_stdout # Reset the standard output to its original value

os.system('chmod 700 rm_sys.sh 2>&1 1>/dev/null')
os.system('bash rm_sys.sh 2> /dev/null')

os.system('chmod 700 run_tleap.sh 2>&1 1>/dev/null')
os.system('bash run_tleap.sh 2>&1 1>/dev/null')


os.system('grep "Volume:" leap.log > temp.txt')
with open("temp.txt", 'r') as f:
  for line in f:
        vol = float(line.split()[1])

vol_lit  = vol * pow(10, -27)
atom_lit = 9.03 * pow(10, 22)
conc = float(Concentration)
num_ion = int(vol_lit * (conc/0.15) * atom_lit)

if Ions == "NaCl":
  pos_neut = "Na+ 0"
  pos_num = "Na+ " + str(num_ion)
  Cl_num = num_ion
else:
  pos_neut = "K+ 0"
  pos_num = "K+ " + str(num_ion)
  Cl_num = num_ion

f = open(tleap, "w")
f.write("""source """ + str(ff) + "\n"
"""source leaprc.DNA.OL15
source leaprc.RNA.OL3
source leaprc.GLYCAM_06j-1
source leaprc.gaff2
source """  + str(water) + "\n"
"""loadamberparams """ + str(ligand_frcmod) + "\n"
"""loadoff """ + str(lib) + "\n"
"""SYS = loadpdb """ + str(protein_ligand) + "\n"
"""alignaxes SYS
check SYS
charge SYS
addions SYS """ + str(pos_neut) + "\n"
"""addions SYS Cl- 0
check SYS
charge SYS
savepdb SYS """ + str(pdb_nw) + "\n"
"""saveamberparm SYS """ + str(top_nw) + " " + str(crd_nw) + "\n"
"""solvatebox SYS """ + str(water_box) + " " + str(size_box) +  """ 0.7 """ + "\n"
"""addIonsRand SYS """ + str(pos_num) + """ Cl- """ + str(Cl_num) + "\n"
"""saveamberparm SYS """ + str(top) + " " + str(crd) + "\n"
"""savepdb SYS """ + str(pdb) + "\n"
"""quit""")
f.close()


os.system('chmod 700 run_tleap.sh 2>&1 1>/dev/null')
os.system('bash run_tleap.sh 2>&1 1>/dev/null')

pdb_amber = os.path.exists(pdb)
top_amber = os.path.exists(top)
crd_amber = os.path.exists(crd)

if pdb_amber == True and top_amber == True and crd_amber == True:
  print("Successfully generated topology! :-)")
else:
  print("ERROR: Check your inputs! ")
os.system('rm *.sh  ANTECHAMBER* ATOMTYPE* temp.txt >/dev/null 2>&1')

# Cell 4
#@title ### **Parameters for MD Equilibration protocol:**

Jobname = options.Jobname

top = os.path.join(workDir, "SYS_gaff2.prmtop")
crd = os.path.join(workDir, "SYS_gaff2.crd")
pdb = os.path.join(workDir, "SYS.pdb")


Minimization_steps = options.Minimization_steps  

Time = options.Time 
stride_time_eq = Time
Integration_timestep = options.Integration_timestep 
dt_eq = Integration_timestep

Temperature = options.Temperature
temperature_eq = Temperature
Pressure = options.Pressure
pressure_eq = Pressure

Force_constant = options.Force_constant

Write_the_trajectory = options.Write_the_trajectory 
write_the_trajectory_eq = Write_the_trajectory

Write_the_log = options.Write_the_log  
write_the_log_eq = Write_the_log

# Cell 5
#@title **Runs an Equilibration MD simulation (NPT ensemble)**

# Defining MD simulation parameters

jobname = os.path.join(workDir, Jobname)
coordinatefile = crd
pdbfile = pdb
topologyfile = top

time_ps = float(Time)*1000
simulation_time = float(time_ps)*picosecond		# in ps
dt = int(dt_eq)*femtosecond
temperature = float(temperature_eq)*kelvin
savcrd_freq = int(write_the_trajectory_eq)*picosecond
print_freq  = int(write_the_log_eq)*picosecond

pressure	= float(pressure_eq)*bar

restraint_fc = int(Force_constant) # kJ/mol

nsteps  = int(simulation_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

#############################################
# Defining functions to use below:
def backup_old_log(pattern, string):
	result = []
	for root, dirs, files in os.walk("./"):
		for name in files:
			if fnmatch.fnmatch(name, pattern):

				try:
					number = int(name[-2])
					avail = isinstance(number, int)
					#print(name,avail)
					if avail == True:
						result.append(number)
				except:
					pass

	if len(result) > 0:
		maxnumber = max(result)
	else:
		maxnumber = 0

	backup_file = "\#" + string + "." + str(maxnumber + 1) + "#"
	os.system("mv " + string + " " + backup_file)
	return backup_file

def restraints(system, crd, fc, restraint_array):

	boxlx = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
	boxly = system.getDefaultPeriodicBoxVectors()[1][1].value_in_unit(nanometers)
	boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometers)

	if fc > 0:
		# positional restraints for all heavy-atoms
		posresPROT = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
		posresPROT.addPerParticleParameter('k')
		posresPROT.addPerParticleParameter('x0')
		posresPROT.addPerParticleParameter('y0')
		posresPROT.addPerParticleParameter('z0')

		for atom1 in restraint_array:
			atom1 = int(atom1)

			xpos  = crd.positions[atom1].value_in_unit(nanometers)[0]
			ypos  = crd.positions[atom1].value_in_unit(nanometers)[1]
			zpos  = crd.positions[atom1].value_in_unit(nanometers)[2]

			posresPROT.addParticle(atom1, [fc, xpos, ypos, zpos])

		system.addForce(posresPROT)

	return system
##############################################

#############################################
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCoordinate file = " + str(coordinatefile))
print("\tPDB file = " + str(pdbfile))
print("\tTopology file = " + str(topologyfile))

print("\n\tSimulation_time = " + str(simulation_time))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps))

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure))
#############################################

print("\n> Setting the system:\n")

if Ligand_Force_field == "OpenFF 2.0.0 (Sage)":
  print("\t- Reading topology and structure file...")
  prmtop = pmd.load_file(topologyfile)
  inpcrd = AmberInpcrdFile(coordinatefile)

  print("\t- Creating system and setting parameters...")
  nonbondedMethod = PME
  nonbondedCutoff = 1.0*nanometers
  ewaldErrorTolerance = 0.0005
  constraints = HBonds
  rigidWater = True
  constraintTolerance = 0.000001
  friction = 1.0
  system = complex_structure.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
                                          constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
else:
  print("\t- Reading topology and structure file...")
  prmtop = AmberPrmtopFile(topologyfile)
  inpcrd = AmberInpcrdFile(coordinatefile)

  print("\t- Creating system and setting parameters...")
  nonbondedMethod = PME
  nonbondedCutoff = 1.0*nanometers
  ewaldErrorTolerance = 0.0005
  constraints = HBonds
  rigidWater = True
  constraintTolerance = 0.000001
  friction = 1.0
  system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
                                          constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)


print("\t- Applying restraints. Force Constant = " + str(Force_constant) + "kJ/mol")
pt_system = pt.iterload(coordinatefile, topologyfile)
pt_topology = pt_system.top
restraint_array = pt.select_atoms('!(:H*) & !(:WAT) & !(:Na+) & !(:Cl-) & !(:Mg+) & !(:K+)', pt_topology)

system = restraints(system, inpcrd, restraint_fc, restraint_array)

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

print("\t- Setting integrator...")
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print("\t- Energy minimization: " + str(Minimization_steps) + " steps")
simulation.minimizeEnergy(tolerance=10*kilojoule/mole, maxIterations=int(Minimization_steps))
print("\t-> Potential Energy = " + str(simulation.context.getState(getEnergy=True).getPotentialEnergy()))

print("\t- Setting initial velocities...")
simulation.context.setVelocitiesToTemperature(temperature)

#############################################
# Running Equilibration on NPT ensemble

dcd_file = jobname + ".dcd"
log_file = jobname + ".log"
rst_file = jobname + ".rst"
prv_rst_file = jobname + ".rst"
pdb_file = jobname + ".pdb"

# Creating a trajectory file and reporters
dcd = DCDReporter(dcd_file, nsavcrd)
firstdcdstep = (nsteps) + nsavcrd
dcd._dcd = DCDFile(dcd._out, simulation.topology, simulation.integrator.getStepSize(), firstdcdstep, nsavcrd) # charmm doesn't like first step to be 0

simulation.reporters.append(dcd)
simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, speed=True, progress=True, totalSteps=nsteps, remainingTime=True, separator='\t\t'))
simulation.reporters.append(StateDataReporter(log_file, nprint, step=True, kineticEnergy=True, potentialEnergy=True, totalEnergy=True, temperature=True, volume=True, speed=True))

print("\n> Simulating " + str(nsteps) + " steps...")
simulation.step(nsteps)

simulation.reporters.clear() # remove all reporters so the next iteration don't trigger them.


##################################
# Writing last frame information of stride
print("\n> Writing state file (" + str(rst_file) + ")...")
state = simulation.context.getState( getPositions=True, getVelocities=True )
with open(rst_file, 'w') as f:
	f.write(XmlSerializer.serialize(state))

last_frame = int(nsteps/nsavcrd)
print("> Writing coordinate file (" + str(pdb_file) + ", frame = " + str(last_frame) + ")...")
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open(pdb_file, 'w'))

print("\n> Finished!\n")

# Cell 6
#@markdown ### **Provide input file names below:**

# Equilibrated_PDB = options.Equilibrated_PDB #@param {type:"string"}
# State_file = options.State_file  #@param {type:"string"}

# #@markdown ---
# #@markdown ### **Parameters for MD Production protocol:**

# Jobname = options.Jobname

# top = os.path.join(workDir, "SYS_gaff2.prmtop")
# crd = os.path.join(workDir, "SYS_gaff2.crd")
# pdb = os.path.join(workDir, "SYS.pdb")

# #@markdown Simulation time (in nanoseconds), number of strides (integers) and integration timestep (in femtoseconds):
# Stride_Time = options.Stride_Time #@param {type:"string"}
# stride_time_prod = Stride_Time
# Number_of_strides = options.Number_of_strides #@param {type:"string"}
# nstride = Number_of_strides
# Integration_timestep = options.Integration_timestep #@param ["0.5", "1", "2", "3", "4"]
# dt_prod = Integration_timestep

# #@markdown Temperature (in Kelvin) and Pressure (in bar)
# Temperature = options.Temperature #@param {type:"string"}
# temperature_prod = Temperature
# Pressure =  options.Pressure #@param {type:"string"}
# pressure_prod = Pressure

# #@markdown Frequency to write the trajectory file (in picoseconds):
# Write_the_trajectory = options.Write_the_trajectory #@param ["10", "100", "200", "500", "1000"]
# write_the_trajectory_prod = Write_the_trajectory

# #@markdown Frequency to write the log file (in picoseconds):
# Write_the_log = options.Write_the_log #@param ["10", "100", "200", "500", "1000"]
# write_the_log_prod = Write_the_log

# #@markdown ---

# Cell 8
#@title **Concatenate and align the trajectory**
#@markdown **Important**: The **Google Drive Path**, **Jobname**, **Number of strides**, **stride time** and **trajectory saved frequency** should be the same you have been used to run your simulation in the previous steps.

Equilibrated_PDB = options.Equilibrated_PDB
Skip = options.Skip #@param ["1", "2", "5", "10", "20", "50"]
stride_traj = Skip
Output_format = options.Output_format  #@param ["dcd", "pdb", "trr", "xtc"]
first_stride = options.first_stride #@param {type:"string"}
Number_of_strides = options.Number_of_strides #@param {type:"string"}
nstride = int(Number_of_strides)
stride_time = options.Stride_time #@param {type:"string"}
trajectory_saved_frequency = options.trajectory_saved_frequency #@param ["10", "100", "200", "500", "1000"]
traj_save_freq = trajectory_saved_frequency
Remove_waters = options.remove_waters #@param ["yes", "no"]
# stride_id_as_ref_for_alignment = "1" #@param {type: "string"}
output_prefix = first_stride+"-"+str(int(first_stride)+nstride-1)

stride_time_ps = float(stride_time)*1000
simulation_time_analysis = stride_time_ps*nstride
simulation_ns = float(stride_time)*int(Number_of_strides)
number_frames = int(simulation_time_analysis)/int(traj_save_freq)
number_frames_analysis = number_frames/int(Skip)


nw_dcd = os.path.join(workDir, str(Jobname) + output_prefix + "_nw." + str(Output_format))
whole_dcd = os.path.join(workDir, str(Jobname) + output_prefix + "_whole." + str(Output_format))
template =  os.path.join(workDir, str(Jobname) + '_%s.dcd')
pdb = os.path.join(workDir, Equilibrated_PDB)

flist = [template % str(i) for i in range(int(first_stride), int(first_stride) + nstride)]

if Remove_waters == "yes":
  #Save topology without waters
  gaff_top = pt.load_topology(os.path.join(workDir, "SYS_gaff2.prmtop"))
  gaff_nw = gaff_top['!:WAT']
  gaff_nw.save(os.path.join(workDir, "SYS_gaff2_nw.prmtop"))
  # Save trajectory without waters
  trajlist = pt.load(flist, os.path.join(workDir, "SYS_gaff2.prmtop"), stride=Skip)
  t0 = trajlist.strip(':WAT')
  traj_image = t0.iterframe(autoimage=True, rmsfit=0)
  traj_nw = pt.write_traj(nw_dcd, traj_image, overwrite=True, options=Output_format)
  traj_dcd_check = os.path.exists(nw_dcd)
  traj = nw_dcd
  pdb_ref = os.path.join(workDir, "SYS_gaff2_nw.prmtop")
else:
  trajlist = pt.load(flist, os.path.join(workDir, "SYS_gaff2.prmtop"), stride=Skip)
  traj_image = trajlist.iterframe(autoimage=True, rmsfit=0)
  traj = pt.write_traj(whole_dcd, traj_image, overwrite=True, options=Output_format)
  traj_dcd_check = os.path.exists(whole_dcd)
  traj = whole_dcd
  pdb_ref = os.path.join(workDir, "SYS_gaff2.prmtop")

traj_load = pt.load(traj, pdb_ref)
print(traj_load)

if traj_dcd_check == True:
  print("Trajectory concatenated successfully! :-)")
else:
  print("ERROR: Check your inputs! ")

# Cell 10
#@title **MM-PBSA method to calculate the binding free energy**
#@markdown **Important:** We will now calculate the interaction energy and solvation free energy for the complex, receptor and ligand and average the results to obtain an estimate of the binding free energy. Please note that we will not perform a calculation of the entropy contribution to binding and so strictly speaking our result will not be a true free energy but could be used to compare against similar systems. We will carry out the binding energy calculation using both the MM-GBSA method and the MM-PBSA method for comparison.

#@markdown Select the GB/SA input parameters,  the "OBC" models (igb=2 and 5) are newer, but appear to give significant improvements and are recommended for most projects (For more information check the Section 4.1 of the [Amber Manual](https://ambermd.org/doc12/Amber20.pdf)):
igb = options.igb #@param ["1", "2", "5", "7", "8"]

def getpreferredencoding(do_setlocale = True):
    return "UTF-8"
locale.getpreferredencoding = getpreferredencoding

if igb == "1":
  mbondi = 'mbondi'
elif igb == "2" or igb == "5":
  mbondi = 'mbondi2'
elif igb == "7":
  mbondi = 'bondi'
elif igb == "8":
  mbondi = 'mbondi3'
else:
  pass

Salt_concentration = '0.15' #@param {type:"string"}
fold_MMPBSA = "MMPBSA_igb_" + igb
#@markdown **Provide output file names below:**
Output_name = 'FINAL_RESULTS_MMPBSA' #@param {type:"string"}

final_mmpbsa = os.path.join(workDir, Output_name)

if number_frames_analysis > 10:
  stride = number_frames_analysis/10
else:
  stride = 1

f = open("mmpbsa.in", "w")
f.write("""&general """  + "\n"
"""  endframe=""" + str(int(number_frames_analysis)) + """,  interval=""" + str(int(stride)) + """, strip_mask=:WAT:Na+:Cl-:Mg+:K+, """ + "\n"
"""/ """ + "\n"
"""&gb """ + "\n"
""" igb=""" + str(igb) +  """, saltcon=""" + str(Salt_concentration) +  """, """ + "\n"
"""/ """ + "\n"
"""&pb """ + "\n"
""" istrng=""" + str(Salt_concentration) +  """, inp=2, radiopt=0, prbrad=1.4, """ + "\n"
"""/""")
f.close()

amberhome = "source /usr/local/amber.sh"
ante_MMPBSA = "ante-MMPBSA.py  -p " + str(pdb_ref) + " -c com.prmtop -r rec.prmtop -l ligand.prmtop -s :WAT:Na+:Cl-:Mg+:K+ -n :LIG --radii " + str(mbondi)
MMPBSA = "MMPBSA.py -O -i mmpbsa.in -o " + str(final_mmpbsa) +  ".dat -sp " + str(pdb_ref) + " -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y "  + str(traj)
mkdir = "mkdir " + os.path.join(workDir, fold_MMPBSA)
mv = "mv _MMPBSA* com.prmtop rec.prmtop ligand.prmtop reference.frc mmpbsa.in " + os.path.join(workDir, fold_MMPBSA)

original_stdout = sys.stdout # Save a reference to the original standard output

with open('run_MMPBSA.sh', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    print(amberhome)
    print(ante_MMPBSA)
    print(MMPBSA)
    print(mkdir)
    print(mv)
    sys.stdout = original_stdout # Reset the standard output to its original value

os.system('chmod 700 run_MMPBSA.sh 2>&1 1>/dev/null')
os.system('bash run_MMPBSA.sh 2>&1 1>/dev/null')

f_mmpbsa = open(final_mmpbsa + '.dat', 'r')
file_contents = f_mmpbsa.read()
print(file_contents)
f_mmpbsa.close()