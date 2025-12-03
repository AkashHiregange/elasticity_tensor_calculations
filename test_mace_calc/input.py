from ase import Atoms
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from ase.io import read, write
import os
from ase.io.trajectory import Trajectory
from ase.constraints import ExpCellFilter
from mace.calculators import mace_mp

crys = read('Co_metal_Opt_rSCAN_12x12x12.traj')

calc = mace_mp(model='large', device='cpu', default_dtype='float64')

crys.set_calculator(calc)

ecf = ExpCellFilter(crys)

# Setup optimisation
dynamics = BFGS(ecf)
traj = Trajectory(f'Co_Opt_mace_mp.traj', 'w', crys)
dynamics.attach(traj)

# Run optimisation
dynamics.run(fmax=0.01)


