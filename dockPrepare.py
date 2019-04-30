#simple script (employing other scripts) to dock a TCR with pMHC

# read in separate TCR structure, pMHC structure (store as files, they will
# get set as poses in Rosetta script)

# select a template based on sequence similarity of TCR (maybe also pMHC). Requires
# labelled chains - default alpha "A", beta "B"

# submit tcr, pmhc, template for docking

# clear up decoy data? Or keep in Rosetta script? Argument-based (optional) - default
# is delete data

# argument for how many cores available

# Note: This script conforms to the labelling of chains used in the ATLAS database.
# the docking file created will relabel chains as:
# TCR alpha chain - "D"
# TCR beta chain - "E"
# peptide chain  - "C"
# MHC chain - "A"

import os
import sys
import argparse
from Bio.PDB import *

from IPython import embed

def args():
	parser = argparse.ArgumentParser(description='Dock a TCR with pMHC.')
	parser.add_argument('-tcr', '--tcr', type=str, help='Input pdb file for TCR', required=True)
	parser.add_argument('-a', '--achain', type=str, help='TCR alpha chain label',required=False, default="D")
	parser.add_argument('-b', '--bchain', type=str, help='TCR beta chain label',required=False, default="E")
	parser.add_argument('-pmhc', '--pmhc', type=str, help='Input pdb file for pMHC', required=True)
	return parser.parse_args()

def relabel(model, old_chain, new_chain):
	model[old_chain].id = new_chain

def detectAntigenChains(pmhc):
	# probably this should be performed at the beginning of the pipeline
	# to avoid performing every time
	# probably extract all these prep stages out of the docking scripts
	if len(pmhc.child_list) > 2:
		print "Error: pmhc has more than two chains"
		sys.exit()
	mhc_chain = max(pmhc.child_list, key=len)
	p_chain = min(pmhc.child_list, key=len)
	return p_chain.id, mhc_chain.id

def prepare(tcr, alpha, beta, pmhc):
	# maybe extract all this to another initial script to prep for docking...
	# with relabelling of chains...
	parser = PDBParser()
	io = PDBIO()
	tcr_structure = parser.get_structure('tcr', tcr)
	pmhc_structure = parser.get_structure('pmhc', pmhc)
	relabel(tcr_structure[0], alpha, "D")
	relabel(tcr_structure[0], beta, "E")

	p_chain, mhc_chain = detectAntigenChains(pmhc_structure[0]) 
	relabel(pmhc_structure[0], p_chain, "C")
	relabel(pmhc_structure[0], mhc_chain, "A")
	
	# add chains from tcr structure to pmhc structure
	for i in tcr_structure[0].child_dict:
		pmhc_structure[0].child_dict[i] = tcr_structure[0].child_dict[i]
	pmhc_structure[0].child_list += tcr_structure[0].child_list
	# set and save structure as same model
	io.set_structure(pmhc_structure)
	output_pdb = "ready_for_docking.pdb"
	io.save(output_pdb)
	return output_pdb




if __name__ == "__main__":
	args = args()
	output_pdb = prepare(args.tcr, args.achain, args.bchain, args.pmhc)
	print "Docking peparation file written to ready_for_docking.pdb"