#IMPORTS
from ase.io import read, write
from ase.build import make_supercell
import sys
import numpy as np
import os
from mace.calculators import mace_mp

# from carmm.phonon.pre_process import make_displaced_supercells, get_charges_and_moments, creating_files_and_directories
sys.path.append(r"C:\Users\\akash\OneDrive - Cardiff University\Desktop\elasticity_tensor_calculations")
from getting_deformed_geometries import generate_deformed_strutures, create_files_and_directories


ase_atoms_eq = read('CoO_Opt_mace_mp.traj')
t = np.array([[-1, 1, 1],[1, -1, 1],[1, 1, -1]])
ase_atoms_eq = make_supercell(ase_atoms_eq, t)


eq_structure, structure, deformations = generate_deformed_strutures(atoms_object=ase_atoms_eq)
create_files_and_directories(eq_structure, structure, deformations)

home = os.getcwd()
for defor in os.listdir():
    if os.path.isdir(f'{home}/{defor}'):
        for file in os.listdir(f'{home}/{defor}'):
            atoms = read(f'{home}/{defor}/geometry.in')
            calc = mace_mp(model='large', device='cpu', default_dtype='float64')
            atoms.set_calculator(calc)
            print(f'{defor}, {atoms.get_potential_energy()}')
            atoms.write(f'{home}/{defor}/{defor}.xyz')
            print(atoms.get_stress(voigt=False))


from elastic_tensor_calculation import read_strain_tensor_from_pkl, read_stress_from_outputs, compute_elasticity_tensor


strain_tensor = read_strain_tensor_from_pkl('strain_tensor.pkl')
stress_tensor = read_stress_from_outputs(output_file_type='.xyz', aims_out_file=False)
elasticity_tensor = compute_elasticity_tensor(strain_tensor, stress_tensor=stress_tensor)

print('This is the computed elasticity tensor')
print(elasticity_tensor)