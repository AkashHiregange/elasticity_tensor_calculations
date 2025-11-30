import os
import shutil
from ase.io import read, write
import pickle
import numpy as np
from pymatgen.analysis.elasticity import Strain, Deformation
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from ase.build import make_supercell

def generate_deformed_strutures(atoms_object, norm_strains = [0.01, 0.025], shear_strains = [0.032, 0.02], write_strain=True):


    '''
    :param atoms_object:
    :param norm_strains: # TODOs: I remember that there was a certain minimum requirement on number of deformed
    structures to compute the elasticity tensor. I will have to go through the documentation to see why is that!!
    :param shear_strains:
    :param write_strain:
    :return:
    '''

    eq_structure = atoms_object
    structure = AseAtomsAdaptor.get_structure(eq_structure)

    deformations: list[Deformation] = []

    strain_list = []
    for ind in [(0, 0), (1, 1), (2, 2)]:
        for amount in norm_strains:
            strain = Strain.from_index_amount(ind, amount)
            strain_list.append(strain)
            deformations.append(strain.get_deformation_matrix())
    for ind in [(0, 1), (0, 2), (1, 2)]:
        for amount in shear_strains:
            strain = Strain.from_index_amount(ind, amount)
            strain_list.append(strain)
            deformations.append(strain.get_deformation_matrix())

    if write_strain:
        strain_tensor = np.array(strain_list)
        with open('strain_tensor.pkl', 'wb') as fp:
            pickle.dump(strain_tensor, fp)

    return eq_structure, structure, deformations

def create_files_and_directories(eq_structure, structure, deformations):
    deformed_struct = [defo.apply_to_structure(structure) for defo in deformations]

    for i, def_struc in enumerate(deformed_struct):
        dir_no = i + 1
        eq_structure_copy = eq_structure.copy()
        directory = f'defor_{dir_no}'
        parent_dir = os.getcwd()
        path_final = os.path.join(parent_dir, directory)
        if os.path.exists(path_final):
            shutil.rmtree(path_final)
        os.mkdir(path_final)
        atoms = AseAtomsAdaptor.get_atoms(def_struc)
        eq_structure_copy.set_cell(atoms.get_cell(), scale_atoms=True)
        eq_structure_copy.write(path_final + '/geometry.in')
        #shutil.copy(parent_dir + '/input.py', path_final + '/input.py')
        #shutil.copy(parent_dir + '/submission.script', path_final + '/submission.script')

ase_atoms_eq = read('CoO_RS_NM_pri_Opt_rSCAN.traj')
t = np.array([[-1, 1, 1],[1, -1, 1],[1, 1, -1]])
ase_atoms_eq = make_supercell(ase_atoms_eq, t)


eq_structure, structure, deformations = generate_deformed_strutures(atoms_object=ase_atoms_eq)
create_files_and_directories(eq_structure, structure, deformations)



''' -------------------------------------------------------------------------------------------------------'''

'''

This section is for getting multiple deformed structures with strain tensor having zeros in some places in the array i.e
normal and shear stress are being applied separately and the deformed structures are obtained.
If you want to have a single deformed structure with all normal and shear stress being applied at the same time, comment
out the section below (until the next section starts) and use the section below this.

'''

# home = os.getcwd()
#
# ase_atoms_eq = read('CoO_RS_NM_pri_Opt_rSCAN.traj')
# t = np.array([[-1, 1, 1],[1, -1, 1],[1, 1, -1]])
# ase_atoms_eq = make_supercell(ase_atoms_eq, t)
#
# structure = AseAtomsAdaptor.get_structure(ase_atoms_eq)
# structure.add_oxidation_state_by_element({"Co": 0, "O": 0})
#
# norm_strains = [0.01, 0.025]
# shear_strains = [0.032, 0.02]
# deformations: list[Deformation] = []
# deformed_struct: list[Structure] = []
#
# strain_list = []
# for ind in [(0, 0), (1, 1), (2, 2)]:
#     for amount in norm_strains:
#         strain = Strain.from_index_amount(ind, amount)
#         # print(strain)
#         strain_list.append(strain)
#         deformations.append(strain.get_deformation_matrix())
# for ind in [(0, 1), (0, 2), (1, 2)]:
#     for amount in shear_strains:
#         strain = Strain.from_index_amount(ind, amount)
#         strain_list.append(strain)
#         # print(strain)
#         deformations.append(strain.get_deformation_matrix())
#
# deformed_struct = [defo.apply_to_structure(structure) for defo in deformations]
# # print(deformed_struct)
# strain_tensor = np.array(strain_list)
# with open('strain_tensor.pkl', 'wb') as fp:
#     pickle.dump(strain_tensor, fp)
#
# for def_struc in deformed_struct:
#     dir_no = deformed_struct.index(def_struc) + 1
#     try:
#         os.mkdir(f'defor_{dir_no}')
#     except Exception as e:
#         print(e)
#         shutil.rmtree(f'{os.getcwd()}/defor_{dir_no}')
#         os.mkdir(f'defor_{dir_no}')
#     os.chdir(f'defor_{dir_no}')
#     atoms = AseAtomsAdaptor.get_atoms(def_struc)
#     ase_atoms_eq.set_cell(atoms.get_cell(), scale_atoms=True)
#     from ase.visualize import view
#     ase_atoms_eq.write('geometry.in')
#     os.chdir(home)
#     shutil.copyfile(f'{home}/input.py', f'defor_{dir_no}/input.py')
#     shutil.copyfile(f'{home}/submission.script', f'defor_{dir_no}/submission.script')

    

''' ---------------------------------------------------------------------------------------------------------------- '''

'''

This section is for getting a single deformed structure with all normal and shear stress being applied at the same time.
If you want to have a multiple deformed structures with all normal and shear stress being applied separately, comment
out the section below and use the section above this.

DON'T USE THIS. FOR CALCULATING ELASTIC CONSTANTS, WE NEED MORE THAN ONE STRESS AND STRAIN TENSOR

'''

# strain = np.zeros((3,3))
# print(strain)
# for ind in [(0,0), (1,1), (2,2)]:
#     for amount in norm_strains:
#         strain = strain + Strain.from_index_amount(ind, amount)
#         # print(strain)
#         # deformations.append(strain.get_deformation_matrix())
# for ind in [(0, 1), (0, 2), (1, 2)]:
#     for amount in shear_strains:
#         strain = strain + Strain.from_index_amount(ind, amount)
#         # print(strain)
#         # deformations.append(strain.get_deformation_matrix())
#
# deformations.append(strain.get_deformation_matrix())
# deformed_struct = [defo.apply_to_structure(structure) for defo in deformations]
# print(deformed_struct)
#
# for def_struc in deformed_struct:
#     os.mkdir(f'all_strain_in_one')
#     os.chdir(f'all_strain_in_one')
#     atoms = AseAtomsAdaptor.get_atoms(def_struc)
#     ase_atoms_eq.set_cell(atoms.get_cell(),scale_atoms=True)
#     ase_atoms_eq.write('geometry.in')
#     shutil.copyfile(home + f'/input.py', home + f'/all_strain_in_one/input.py')
#     shutil.copyfile(home + f'/submission.script', home + f'/all_strain_in_one/submission.script')
#     os.chdir(home)

''' ---------------------------------------------------------------------------------------------------------------- '''





