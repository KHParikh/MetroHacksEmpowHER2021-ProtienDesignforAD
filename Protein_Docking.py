from __future__ import print_function

import optparse

from rosetta import *
from pyrosetta import *
from rosetta.protocols.rigid import *


init(extra_options = "-constant_seed")
import os; os.chdir('.test.output')


def sample_docking(pdb_filename, partners,
        translation = 3.0, rotation = 8.0,
        jobs = 1, job_output = 'SF2021DeNovooDocked2'):
    # Create a pose from the desired PDB file
    pose = Pose()
    pose_from_file(pose, pdb_filename)
    print ("pose created")

    # Setup the docking FoldTree
    # using this method, the jump number 1 is automatically set to be the
    #    inter-body jump
    dock_jump = 1
    pyrosetta.rosetta.protocols.docking.setup_foldtree(pose, partners, Vector1([dock_jump]))

    # create centroid <--> fullatom conversion Movers
    to_centroid = SwitchResidueTypeSetMover('centroid')
    to_fullatom = SwitchResidueTypeSetMover('fa_standard')
    recover_sidechains = pyrosetta.rosetta.protocols.simple_moves.ReturnSidechainMover(pose)

    # Convert to centroid (general regions instead of specific side chains)
    to_centroid.apply(pose)

    # Create a (centroid) test pose
    test_pose = Pose()
    test_pose.assign(pose)

    # Create ScoreFunctions for centroid and fullatom docking
    scorefxn_low = create_score_function('interchain_cen')
    scorefxn_high = create_score_function('docking')

    scorefxn_high_min = create_score_function('docking', 'docking_min')

    # Create Movers for producing an initial perturbation of the structure
    randomize_upstream = RigidBodyRandomizeMover(pose, dock_jump,
        partner_upstream)
    randomize_downstream = RigidBodyRandomizeMover(pose, dock_jump,
        partner_downstream)
    dock_pert = RigidBodyPerturbMover(dock_jump, translation, rotation)
    spin = RigidBodySpinMover(dock_jump)
    slide_into_contact = pyrosetta.rosetta.protocols.docking.DockingSlideIntoContact(dock_jump)

    # Setup the MinMover
    # the MoveMap can set jumps (by jump number) as degrees of freedom
    movemap = MoveMap()
    movemap.set_jump(dock_jump, True)
    # the MinMover can minimize score based on a jump degree of freedom, this
    #    will find the distance between the docking partners which minimizes
    #    the score
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn_high_min)

    # Create a SequenceMover for the perturbation step
    perturb = pyrosetta.rosetta.protocols.moves.SequenceMover()
    perturb.add_mover(randomize_upstream)
    perturb.add_mover(randomize_downstream)
    perturb.add_mover(dock_pert)
    perturb.add_mover(spin)
    perturb.add_mover(slide_into_contact)
    perturb.add_mover(to_fullatom)
    perturb.add_mover(recover_sidechains)
    perturb.add_mover(minmover)

    # Setup the DockingProtocol
    dock_prot = pyrosetta.rosetta.protocols.docking.DockingProtocol() 
    dock_prot.set_movable_jumps(Vector1([1]))
    dock_prot.set_lowres_scorefxn(scorefxn_low)
    dock_prot.set_highres_scorefxn(scorefxn_high_min)
 
    # Setup the PyJobDistributor
    jd =  PyJobDistributor(job_output, jobs, scorefxn_high)
    temp_pose = Pose()    # a temporary pose to export to PyMOL
    temp_pose.assign(pose)
    to_fullatom.apply(temp_pose)    # the original pose was fullatom
    recover_sidechains.apply(temp_pose)    # with these sidechains
    jd.native_pose = temp_pose    # for RMSD comparison

    # Perform protein-protein docking
    counter = 0
    while not jd.job_complete:
        test_pose.assign(pose)
        counter += 1
        test_pose.pdb_info().name(job_output + '_' + str(counter))
        perturb.apply(test_pose)
        dock_prot.apply(test_pose)
        to_fullatom.apply(test_pose)
        test_pose.pdb_info().name(job_output + '_' + str( counter ) + '_fa')
        jd.output_decoy(test_pose)
