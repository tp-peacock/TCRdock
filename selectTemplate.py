import pandas as pd
from Bio import pairwise2
import argparse
from IPython import embed

def args():
	parser = argparse.ArgumentParser(description='Select TCR template by sequence alignment')
	parser.add_argument('-a', '--achain', type=str, help='alpha chain sequence', required=True)
	parser.add_argument('-b', '--bchain', type=str, help='beta chain sequence', required=True)
	parser.add_argument('-e', '--exclude_list',  nargs='+', help='Templates pdb codes to exclude', required=False)
	return parser.parse_args()

def takeAlignScore(elem):
    return elem[1][2]

def takeSecond(elem):
    return elem[1]

def getAlignments(tdata, in_seq, chain, exclude = []):
	alignments = {}

	if chain in ["a","A","alpha"]:
		for index, row in tdata.iterrows():
			if index in exclude:
				continue
			a = pairwise2.align.globalxx(row.alpha_seq, in_seq) # consider only top alignment
			alignments[index] = a[0]
	
	elif chain in ["b", "B", "beta"]:
		for index, row in tdata.iterrows():
			if index in exclude:
				continue
			a = pairwise2.align.globalxx(row.beta_seq, in_seq) # consider only top alignment
			alignments[index] = a[0]
	
	return alignments

def getBestAlignment(a_alignments, b_alignments):
	# at the moment, take only based on sequence similarity.
	# in the future, probably also consider loop length
	totals = []
	for pdb in a_alignments:
		totals.append([pdb, a_alignments[pdb][2] + b_alignments[pdb][2]])
	totals.sort(key=takeSecond, reverse=True)	
	# maybe need to re-rank based on pMHC for close scores...
	return totals[0][0] # highest score is tempate


if __name__ == "__main__":

	tdata = pd.DataFrame.from_csv('templates.csv')

	args = args()

#	exclude_list = ["1kj2", "3vxm"]

	#test case 3tjh
#	in_seq_a = "AQSVTQPDARVTVSEGASLQLRCKYSYSATPYLFWYVQYPRQGLQMLLKYYSGDPVVQGVNGFEAEFSKSDSSFHLRKASVHWSDSAVYFCAVSAKGTGSKLSFGKGAKLTVS"
#	in_seq_b = "AAVTQSPRNKVTVTGGNVTLSCRQTNSHNYMYWYRQDTGHGLRLIHYSYGAGNLQIGDVPDGYKATRTTQEDFFLLLELASPSQTSLYFCASSDAPGQLYFGEGSKLTVL"

	#test case 1ao7
#	in_seq_a = "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVT"
#	in_seq_b = "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVT"
	
	#3vxm
#	in_seq_a = "AQSVTQPDIHITVSEGASLELRCNYSYGATPYLFWYVQSPGQGLQLLLKYFSGDTLVQGIKGFEAEFKRSQSSFNLRKPSVHWSDAAEYFCAVGAPSGAGSYQLTFGKGTKLSVI"
	#1kj2
#	in_seq_b = "VTLLEQNPRWRLVPRGQAVNLRCILKNSQYPWMSWYQQDLQKQLQWLFTLRSPGDKEVKSLPGADYLATRVTDTELRLQVANMSQGRTLYCTCSAAPDWGASAETLYFGSGTRLTVL"

	a_alignments = getAlignments(tdata, args.achain, "alpha", exclude = args.exclude_list)
	b_alignments = getAlignments(tdata, args.bchain, "beta")

	template = getBestAlignment(a_alignments, b_alignments)
	print "template is", template
	embed()