from ase import Atoms
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from ase.io import read, write
import os
from ase.io.trajectory import Trajectory
from ase.constraints import ExpCellFilter

from carmm.run.aims_path import set_aims_command
set_aims_command(hpc="isambard", basis_set="intermediate", defaults=2020)

crys = read('geometry.in')

from carmm.run.aims_calculator import get_aims_calculator
fhi_calc = get_aims_calculator(dimensions=3, xc='libxc MGGA_X_RSCAN+MGGA_C_RSCAN')

from carmm.run.aims_calculator import get_k_grid
k_grid = get_k_grid(crys, sampling_density=0.0234)

fhi_calc.set(
         spin='none',
         compute_forces=True,
         compute_numerical_stress=True,
         k_grid=k_grid,
         relativistic=('atomic_zora','scalar'),
         occupation_type=('gaussian','0.01'),
         n_max_pulay=20,
         charge_mix_param=0.05,
         sc_accuracy_etot=1e-6,
         sc_accuracy_eev=1e-3,
         sc_accuracy_rho=1e-6,
         sc_accuracy_forces=1e-4,
        )
crys.set_calculator(fhi_calc)

energy = crys.get_potential_energy()
print(f" Total energy: {crys.get_potential_energy()}")
print(f" Stress: {crys.get_stress()}")
#print(f" Cell parameters: {crys.get_cell_lengths_and_angles()}")

os.chdir(home)


