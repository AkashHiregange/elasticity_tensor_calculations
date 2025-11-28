#IMPORTS
from ase.io import read, write
from ase.build import make_supercell
import sys
import numpy as np
# from carmm.phonon.pre_process import make_displaced_supercells, get_charges_and_moments, creating_files_and_directories
sys.path.append(r"C:\Users\\akash\OneDrive - Cardiff University\Desktop\elasticity_tensor_calculations")
from getting_deformed_geometries import generate_deformed_strutures, create_files_and_directories


ase_atoms_eq = read('CoO_RS_NM_pri_Opt_rSCAN.traj')
t = np.array([[-1, 1, 1],[1, -1, 1],[1, 1, -1]])
ase_atoms_eq = make_supercell(ase_atoms_eq, t)


eq_structure, structure, deformations = generate_deformed_strutures(atoms_object=ase_atoms_eq)
create_files_and_directories(eq_structure, structure, deformations)

