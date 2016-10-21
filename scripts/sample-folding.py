from rosetta import *
rosetta.init()
import random
from random import *

# Creating a pose object which stores all the residues and their phi and psi angles
p = Pose()

#read fasta file for sequence
with open("../fasta/1LGL_fasta.txt",'r') as fasta:
    sequence = fasta.readlines()[1].strip()

#making an initial pose with A chain of CueR
make_pose_from_sequence(p, sequence, "centroid")

#Checking initial phi, psi angle
for i in range(1, p.total_residue() + 1):
	print i, " phi = ", p.phi(i), "psi = ", p.psi(i)

#setup score function and weights
scorefxn = ScoreFunction()
scorefxn.set_weight(hbond_lr_bb, 1.17)
scorefxn.set_weight(hbond_sr_bb, 1.17)
scorefxn.set_weight(vdw, 1.0)
#scorefxn.set_weight(env, 1.0)
#scorefxn.set_weight(pair, 1.0)
#scorefxn.set_weight(cbeta, 1.0)

#set up simulation parameters
ncycles = 100000
kT = 1.0
mc = MonteCarlo(p, scorefxn, kT)

#set up Mover
movemap = MoveMap()
movemap.set_bb(True)

# fragment mover
fragset = ConstantLengthFragSet(9)
fragset.read_fragment_file("../fraglibs/1LGL_9.txt")
cost = GunnCost()

frag_mover = SmoothFragmentMover(fragset, movemap, cost)

#runing simulation through monte carlo
for i in range(1, ncycles):
    print i
    frag_mover.apply(p)
    mc.boltzmann(p)
    mc.show_scores()
    mc.show_counters()
    mc.show_state()
    mc.recover_low(p)

dump_pdb(p, "../models/1LGL.pdb")

