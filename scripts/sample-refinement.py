from rosetta import *
rosetta.init()
import random
from random import *

# Creating a pose object which stores all the residues and their phi and psi angles
p = Pose()
start = pose_from_pdb("../models/1LGL.pdb")
p.assign(start)

#read fasta file for sequence
with open("../fasta/1LGL_fasta.txt",'r') as fasta:
    sequence = fasta.readlines()[1].strip()

#making an initial pose with A chain of CueR
make_pose_from_sequence(p, sequence, "fa_standard")

#Checking initial phi, psi angle
for i in range(1, p.total_residue() + 1):
	print i, " phi = ", p.phi(i), "psi = ", p.psi(i)

#setup score function and weights
scorefxn = ScoreFunction()
scorefxn.set_weight(fa_atr, 0.8)
scorefxn.set_weight(fa_rep, 0.44)
scorefxn.set_weight(fa_intra_rep, 0.004)
scorefxn.set_weight(fa_sol, 0.75)
scorefxn.set_weight(fa_elec, 0.07)
scorefxn.set_weight(hbond_sr_bb, 1.7)
scorefxn.set_weight(hbond_lr_bb, 1.7)
scorefxn.set_weight(hbond_bb_sc, 1.7)
scorefxn.set_weight(hbond_sc, 1.1)
scorefxn.set_weight(rama, 0.2)
scorefxn.set_weight(ref, 1)
scorefxn.set_weight(fa_dun, 0.56)

#set up simulation parameters
ncycles = 6000
kT = 1.0
mc = MonteCarlo(p, scorefxn, kT)
n_moves = 10

#set up Movemap
movemap = MoveMap()
movemap.set_bb(True)

#smallmover, shearmover and minmover
small_mover = SmallMover(movemap, kT, n_moves)
shear_mover = ShearMover(movemap, kT, n_moves)
min_mover = MinMover()

# sequence mover mover
seq_mover = SequenceMover()
seq_mover.add_mover(small_mover)
seq_mover.add_mover(shear_mover)
seq_mover.add_mover(min_mover)


#runing simulation through monte carlo
for i in range(1, ncycles):
    print i
    seq_mover.apply(p)
    mc.boltzmann(p)
    mc.show_scores()
    mc.show_counters()
    mc.show_state()
    mc.recover_low(p)

dump_pdb(p, "../models/1LGL-refined.pdb")

