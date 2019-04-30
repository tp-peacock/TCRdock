import time
start_time = time.time()


from IPython import embed
import dockPrepare
import selectTemplate

tcr_in = "/home/tom/phd/pdbTestFiles/1ao7transrot_tcr_only.pdb"
alpha = "D"
beta = "E"
pmhc_in = "/home/tom/phd/pdbTestFiles/1ao7transrot_pmhc_only.pdb"
dockingpdb = dockPrepare.prepare(tcr_in, alpha, beta, pmhc_in)
templates = "/home/tom/phd/templateBuilder/templates.csv"
templatepdb = selectTemplate.select(templates, "D", "E")

#templatepdb = "/home/ucbptpe/TCRdock/atlas_3QDJ.pdb"



embed()



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
from mpi4py import MPI



OUTDATA = []
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

tcr_in = "/home/tom/phd/pdbTestFiles/1ao7transrot_tcr_only.pdb"
alpha = "D"
beta = "E"
pmhc_in = "/home/tom/phd/pdbTestFiles/1ao7transrot_pmhc_only.pdb"
dockingpdb = dockPrepare.prepare(tcr_in, alpha, beta, pmhc_in)



templatepdb = "/home/ucbptpe/TCRdock/atlas_3QDJ.pdb"

#dockingpdb = "/home/ucbptpe/RosettaTest/1ao7transrot.pdb"
scorefxn_low = create_score_function('interchain_cen')
pose = Pose()
ready_for_docking = pose_from_pdb(dockingpdb)

recover_sidechains = rosetta.protocols.simple_moves.ReturnSidechainMover(ready_for_docking)

fa_switch = SwitchResidueTypeSetMover('fa_standard')
censwitch = SwitchResidueTypeSetMover('centroid')
centroid_complex = Pose()
centroid_complex.assign(ready_for_docking)
censwitch.apply(centroid_complex)
jump = Vector1([1])
rosetta.protocols.docking.setup_foldtree(centroid_complex, "AC_DE", jump)
scorefxn_high = get_fa_scorefxn()
#jd = pyrosetta.PyJobDistributor('/home/ucbptpe/Scratch/docking_lowres/rmsd_test2/rmsd_dock', 100, scorefxn_high)

native_pose = Pose()
native_pose.assign(pose_from_pdb("/home/ucbptpe/TCRdock/1AO7.pdb"))

template_pose = Pose()
template_pose.assign(pose_from_pdb(templatepdb))

task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(ready_for_docking)
task.or_include_current(True)
task.restrict_to_repacking()
pack_mover = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn_high, task)
print pack_mover

#jd.native_pose = full_atom_pose

import pyrosetta.mpi
from pyrosetta import logger

pyrosetta.mpi.mpi_init()

def my_function(decoy_num):
  filename = "/home/ucbptpe/Scratch/docking_20032019/decoys/decoy_" + str(decoy_num)

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
  print 'Pre packing high score:',  scorefxn_high(p)
  print 'Pre packing low score:',  scorefxn_low(p)
  print "starting pack mover:"
  pack_mover_start_time = time.time()
  pack_mover.apply(p)
  print "Time for pack mover:", time.time() - pack_mover_start_time
  print 'Post packing high score:',  scorefxn_high(p)
  print 'Post packing low score:',  scorefxn_low(p)
  rmsd_template = core.scoring.CA_rmsd(p, template_pose)
  rmsd_native = core.scoring.CA_rmsd(p,native_pose)
  OUTDATA.append([filename,scorefxn_low(p),rmsd_template, rmsd_native])
  p.dump_pdb(filename)

#while not jd.job_complete:
#  p = Pose()
#  p.assign(centroid_complex)
#  randomize1 = rigid_moves.RigidBodyRandomizeMover(p,1,rigid_moves.partner_upstream)
#  randomize2 = rigid_moves.RigidBodyRandomizeMover(p,1,rigid_moves.partner_downstream)
#  randomize1.apply(p)
#  randomize2.apply(p)
#  slide = rosetta.protocols.docking.DockingSlideIntoContact(1)
#  slide.apply(p)
#  dock_lowres = rosetta.protocols.docking.DockingLowRes(scorefxn_low,1)
#  dock_lowres.apply(p)
#  fa_switch.apply(p)
#  recover_sidechains.apply(p)
#  print 'Pre packing high score:',  scorefxn_high(p)
#  print 'Pre packing low score:',  scorefxn_low(p) 
#  print "starting pack mover:"
#  pack_mover_start_time = time.time()
#  pack_mover.apply(p)
#  print "Time for pack mover:", time.time() - pack_mover_start_time
#  print 'Post packing high score:',  scorefxn_high(p)
#  print 'Post packing low score:',  scorefxn_low(p)  
#  jd.output_decoy(p)

pyrosetta.mpi.MPIJobDistributor(10000, my_function)

#total_time = time.time() - start_time
#print "Total time taken:", total_time
#print ""
data = MPI.COMM_WORLD.gather(OUTDATA, root=0)
if rank == 0:
  outdata = open("/home/ucbptpe/Scratch/docking_20032019/data_out.csv","w")
  print "------------------------------summary---------------------------------"
  print "docking input:", dockingpdb
  for i in data:
    for j in i:
      print j
      for k in j:
        outdata.write(str(k)+", ")
      outdata.write("\n")
  outdata.close()

  total_time = time.time() - start_time
  print "Total run time:", total_time
  print "----------------------------------------------------------------------"
  print ""



