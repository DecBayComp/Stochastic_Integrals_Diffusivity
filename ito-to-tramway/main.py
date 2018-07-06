"""
List of commands launched by the job manager for tessellation, inference, calculation of Bayes factors and results output with the help of TRamWAy
"""


import argparse			# for command-line arguments
from tesselate_and_infer import tesselate_and_infer
from calculate import calculate

from constants import version, folders, bl_produce_maps, dt, localization_error


def main(arg_str):
    """
    Main analysis file that uses TRamWAy
    """

    # Define arguments
    arg_parser = argparse.ArgumentParser(
        description='TRamWAy wrapper for analysis of a random walk trajectory')
    arg_parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s ' + str(version))
    arg_parser.add_argument('-f', '--file', required=True, action='store', type=str,
                            help='Path to a trajectory')

    # Analyze arguments
    input_args = arg_parser.parse_args(arg_str.split())

    # Use the analyzed arguments
    file = input_args.file

    print(file)
    # Tesselate and perform inference
    tesselate_and_infer(file, localization_error=localization_error)

    # Calculate Bayes factors and output results
    _, output_folder = folders()
    calculate(file, output_folder, bl_produce_maps=False, dt=dt,
              snr_label='snr', localization_error=localization_error)
