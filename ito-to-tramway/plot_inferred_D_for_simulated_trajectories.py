
try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2
    has_run = 1
else:
    print("Graphic interface NOT re-initialized")

# import numpy as np
# from scipy.integrate import dblquad, quad
# from scipy.special import factorial, gammainc, gammaln

from calculate import calculate
from constants import dt
# import this
from tesselate_and_infer import tesselate_and_infer

# from tramway.inference.bayes_factors.calculate_bayes_factors import \
#     calculate_bayes_factors
# from tramway.inference.bayes_factors.calculate_marginalized_integral import (calculate_integral_ratio,
#                                                                              calculate_marginalized_integral)
# from tramway.inference.bayes_factors.calculate_posteriors import (calculate_one_1D_posterior_in_2D,
#                                                                   calculate_one_1D_prior_in_2D,
#                                                                   calculate_one_2D_posterior)
# from tramway.inference.bayes_factors.convenience_functions import n_pi_func, p

# %% Produce a diffusivity map for one of the simulated trajectories for the article
file = r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory\sim_data_000000050.csv'
output_folder = r'D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\input\diffusivity_map_for_article_sim_trajectory\result'

tesselate_and_infer(file, localization_error=0, load=False)
calculate(file, output_folder, bl_produce_maps=True,
          snr_label='snr', localization_error=0, ticks=True)
