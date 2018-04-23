

import argparse			# for command-line arguments
from tesselate_and_infer import tesselate_and_infer
from calculate import calculate

from constants import version, output_folder, bl_produce_maps

def main(arg_str):
	"""
	Main analysis file that uses TRamWAy
	"""
	
	## Define arguments
	arg_parser = argparse.ArgumentParser(description = 'TRamWAy wrapper for analysis of a random walk trajectory')
	arg_parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s ' + str(version))
	arg_parser.add_argument('-f', '--file', required = True, action = 'store', type = str, 
		help = 'Path to a trajectory')


	## Analyze arguments
	input_args = arg_parser.parse_args(arg_str.split())


	## Use the analyzed arguments
	file = input_args.file

	## Tesselate and perform inference
	tesselate_and_infer(file)

	## Calculate Bayes factors and output results
	calculate(file, output_folder, bl_output_map = False)

