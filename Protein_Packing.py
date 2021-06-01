from __future__ import print_function
from rosetta import *
from pyrosetta import *
init()
import os; os.chdir('.test.output')

print('Packing and Design ----------------------------------------------')

print('mover: PackRotamersMover')
pose = core.import_pose.pose_from_file("../test/data/FoldedBackbone.pdb")
#New configurations evaluated using the Rosetta all-atom energy function
scorefxn = get_fa_scorefxn()
task_pack = standard_packer_task(pose)
#Restrict residues 1, 11-21 to swapping rotamers of SAME residue
#Packing is unrestricted on residues 12-20, so packers facilitate sequence design
#Each swap is in the simulated annealing search is tested against the Metropolis criterion
task_pack.nonconst_residue_task(1).restrict_to_repacking()
for i in range (11, 21):
    task_pack.nonconst_residue_task(i).restrict_to_repacking()
    print(i,"--------------s-------")
 
print(task_pack)
#Ramdom sampling from a rotamer library to optimize side chain conformations in the pose
pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover( scorefxn, task_pack )
pack.apply(pose)
#Faster version of standard packer (PackRotamersMover) finds local minima in side-chain conformation space
rotamer_trials = pyrosetta.rosetta.protocols.minimization_packing.RotamerTrialsMover(scorefxn, task_pack)
rotamer_trials.apply(pose)
print('rotamer trials--------------------------------------------------------------------------')
print(pack)
print('printed pack--------------------------------------------------------------------------')


task_design = core.pack.task.TaskFactory.create_packer_task(pose)
task_design.nonconst_residue_task(1).restrict_to_repacking()
for i in range (11, 21):
    task_design.nonconst_residue_task(i).restrict_to_repacking()
    print(i,"--------------s-------")
print(task_design)
print('printed task design--------------------------------------------------------------------------')
pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover( scorefxn, task_design )
pack.apply(pose)

tf = standard_task_factory

print('restrict to repacking--------------------------------------------------------------------------')
pr = core.pack.task.operation.RestrictResidueToRepacking()
pr.include_residue(1)
for i in range (11, 21):
    pr.include_residue(i)
    tf.push_back(pr)
    
print('almost ending--------------------------------------------------------------------------')

pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover( scorefxn )
pack.task_factory(tf)
pack.apply(pose)
dump_pdb(pose,'PackedBackbone.pdb')


