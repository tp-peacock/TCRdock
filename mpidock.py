from pyrosetta import *
from pyrosetta.rosetta import *
pyrosetta.init(extra_options="-extrachi_cutoff 12 -ex1 -ex2 -ex3")
from pyrosetta.toolbox import mutate_residue
import pyrosetta.rosetta.protocols.grafting
from rosetta.core.pack.dunbrack import *
from pyrosetta.toolbox import pose_from_rcsb
import rosetta.protocols.rigid
import rosetta.protocols.rigid as rigid_moves
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta.protocols.loops.loop_closure.ccd import *
import pyrosetta.mpi
from pyrosetta import logger

import time
from IPython import embed

class Protocol:
  def __init__(self):
    self.data = {}

  def run():
    p = Pose()
    p.assign(centroid_complex)
    randomize1 = rigid_moves.RigidBodyRandomizeMover(p,1,rigid_moves.partner_upstream)
    randomize2 = rigid_moves.RigidBodyRandomizeMover(p,1,rigid_moves.partner_downstream)
    randomize1.apply(p)
    randomize2.apply(p)
    slide = rosetta.protocols.docking.DockingSlideIntoContact(1)
    slide.apply(p)
    dock_lowres = rosetta.protocols.docking.DockingLowRes(scorefxn_low,1)
    dock_lowres.apply(p)
    fa_switch.apply(p)
    recover_sidechains.apply(p)
    # print 'Pre packing high score:',  scorefxn_high(p)
    # print 'Pre packing low score:',  scorefxn_low(p)
    # print "starting pack mover:"
    pack_mover_start_time = time.time()
    pack_mover.apply(p)
    # print "Time for pack mover:", time.time() - pack_mover_start_time
    # print 'Post packing high score:',  scorefxn_high(p)
    # print 'Post packing low score:',  scorefxn_low(p)
    return p
    

def dock(protocol, decoy_num):
 # filename = "/home/ucbptpe/Scratch/docking_26022019/decoys/decoy_" + str(decoy_num) + ".pdb"
  filename = "decoy_" + str(decoy_num) + ".pdb"
  pose = protocol.run()
  pose.dump_pdb(filename)
  score = scorefxn_high(pose)
  protocol.data[decoy_num] = [score,filename]


start_time = time.time()

dockingpdb = "/home/ucbptpe/RosettaTest/1ao7transrot.pdb"
pose = Pose()
ready_for_docking = pose_from_pdb(dockingpdb)

scorefxn_low = create_score_function('interchain_cen')
recover_sidechains = rosetta.protocols.simple_moves.ReturnSidechainMover(ready_for_docking)
fa_switch = SwitchResidueTypeSetMover('fa_standard')
censwitch = SwitchResidueTypeSetMover('centroid')

centroid_complex = Pose()
centroid_complex.assign(ready_for_docking)
censwitch.apply(centroid_complex)

jump = Vector1([1])
rosetta.protocols.docking.setup_foldtree(centroid_complex, "AC_DE", jump)
scorefxn_high = get_fa_scorefxn()

full_atom_pose = Pose()
full_atom_pose.assign(pose_from_pdb("/home/ucbptpe/RosettaTest/1AO7.pdb"))

task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(full_atom_pose)
task.or_include_current(True)
task.restrict_to_repacking()
pack_mover = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn_high, task)

#jd.native_pose = full_atom_pose

pyrosetta.mpi.mpi_init()
protocol = Protocol()
pyrosetta.mpi.MPIJobDistributor(protocol,10, dock)

total_time = time.time() - start_time
print "Total time taken:", total_time

print protocol.data
embed() 