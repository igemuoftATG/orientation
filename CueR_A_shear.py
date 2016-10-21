#Cite paper
#PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta
import random
from random import *
import rosetta
rosetta.init()


model_object = Pose()

#making an initial pose with A chain of CueR
#make_pose_from_sequence(init_pose, "MNISDVAKITGLTSKAIRFYEEKGLVTPPMRSENGYRTYTQQHLNELTLLRQARQVGFNLEESGELVNLFNDPQRHSADVKRRTLEKVAEIERHIEELQSMRDQLLALANACPGDDSADCPIIENLSGCCHHRAG","centroid_rot")
seq = pose_from_sequence(model_object, "MNISDVAKITGLTSKAIRFYEEKGLVTPPMRSENGYRTYTQQHLNELTLLRQARQVGFNLEESGELVNLFNDPQRHSADVKRRTLEKVAEIERHIEELQSMRDQLLALANACPGDDSADCPIIENLSGCCHHRAG","centroid")

model_object.assign(seq)

#p.assign(init_pose)

#Checking initial phi, psi angle
for i in range(1, model_object.total_residue() + 1):
	print (i, " phi = ", model_object.phi(i), "psi = ", model_object.psi(i))

#setup score function
scorefxn = ScoreFunction()
scorefxn.set_weight(hbond_lr_bb, 1.0)
scorefxn.set_weight(vdw, 1.0)
scorefxn.set_weight(env, 1.0)
scorefxn.set_weight(pair, 1.0)
scorefxn.set_weight(cbeta, 1.0)

#set up simulation parameters
ncycles = 500
kT = 1.0
mc = MonteCarlo(model_object, scorefxn, kT)

#set up Mover
movemap = MoveMap()
movemap.set_bb(True)
shear_mover = ShearMover(movemap, kT, 5)

#run simulation
for i in range(1, ncycles):
	print (i)
	shear_mover.apply(model_object)
	mc.boltzmann(model_object)
	mc.show_scores()
	mc.show_counters()
	mc.show_state()

#dump into pdb file
mc.recover_low(model_object)
dump_pdb(model_object, "./output.pdb")
