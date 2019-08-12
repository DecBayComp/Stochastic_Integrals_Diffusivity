"""
List of commands launched by the job manager for tessellation, inference, calculation of Bayes factors and results output with the help of TRamWAy
"""


import argparse  # for command-line arguments

from calculate import calculate
from constants import bl_produce_maps, dt, folders, version
from tesselate_and_infer import tesselate_and_infer

snr_label = 'snr'
extension = '.csv'
sigma = 0
kwargs = {
    'method': 'gwr',
    'force': True,
    'max_iter': 50,
    'ref_distance': 0.1,
    'sep': ';',
    'dt': dt,
    #     'method': 'hexagon', 'avg_distance': 0.1, 'min_location_count': 0,
    # 'method': 'kohonen', 'min_probability': 20,  # 'avg_probability': 20,
    # 'avg_distance': 0.1,
    #     'reset_origin': False,
    #
    # 'avg_location_count': 20,
    # 'ref_distance': 0,
    # 'prune': False,
    # 'xlims': np.array([0, 2]),
    # 'ylims': np.array([0, 2]),
    # 'clip_grad': 0.9,
    # 'clip_alpha': 0.95,
}


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

    # Calculate Bayes factors and output results
    _, output_folder = folders()
    calculate(file, output_folder, bl_produce_maps=False,
              snr_label=snr_label, sigma=sigma, extension=extension, **kwargs)
