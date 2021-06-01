#STEP 1: PROTEIN FOLDING
#Create an optimally folded backbone by sampling torsion angles which maximize hydrogen bonding and van der Waals interactions between residues
from pyrosetta import *
from random import *
from pyrosetta.rosetta.core.scoring import *
init()
p = Pose()
p = pose_from_pdb("UnfoldedBackbone.pdb")
#Score function is computed using two equally weighted terms representing short-range hydrogen bonds and van der Waals interactions
scorefxn = ScoreFunction()
scorefxn.set_weight(hbond_sr_bb, 1.0)
scorefxn.set_weight(vdw, 1.0)
#Monte Carlo algorithm accepts or rejects new protein conformations based on the Metropolis criterion
mc = MonteCarlo(p, scorefxn, 1.0)
#Subroutine selects a random residue, preturbs its phi or psi torsion angle by a number chosen from a Gaussian distribution
def perturb_bb(pose):
    resnum = randint(1, pose.total_residue())
    pose.set_phi(resnum, pose.phi(resnum)-25+random()*50)
    pose.set_psi(resnum, pose.psi(resnum)-25+random()*50)
    return pose
#Loop iterates through subroutine 100000 times
def my_protocol(pose):
    mc.reset(pose)
    for i in range(1,100000):
        perturb_bb(p)
        print (scorefxn(p))
        mc.boltzmann(p) 

        if (i%1000 == 0):
            mc.recover_low(p)
    #Output lowest-energy structure
    mc.recover_low(p)
    return pose
my_protocol(p)
dump_pdb(p, "FoldedBackbone.pdb") 
mc.show_scores()
mc.show_state()
