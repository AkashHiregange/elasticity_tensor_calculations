#IMPORTS
from ase.io import read, write
from ase.build import make_supercell
import sys
import numpy as np
# from carmm.phonon.pre_process import make_displaced_supercells, get_charges_and_moments, creating_files_and_directories
sys.path.append(r"C:\Users\\akash\OneDrive - Cardiff University\Desktop\elasticity_tensor_calculations")
from elastic_tensor_calculation import read_strain_tensor_from_pkl, read_stress_from_outputs, compute_elasticity_tensor


strain_tensor = read_strain_tensor_from_pkl('strain_tensor.pkl')
stress_tensor = read_stress_from_outputs(aims_out_file=True)
elasticity_tensor = compute_elasticity_tensor(strain_tensor, stress_tensor=stress_tensor)

print('This is the computed elasticity tensor')
print(elasticity_tensor)