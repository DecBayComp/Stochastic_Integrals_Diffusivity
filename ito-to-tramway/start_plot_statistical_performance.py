"""
This interactve python file combines data for the article statistical performance plot and then makes the plot.
It doesn't perform the inference itself. It should have already been performed by launching the `start_me.py`.
"""

# try:
#     has_run
# except NameError:
# %matplotlib
#     %load_ext autoreload
#     %autoreload 2
#
#     has_run = 1
# else:
#     print("Graphic interface NOT re-initialized")

# import matplotlib
import numpy as np

from combine_results import combine_results
from estimate_theoretical_performance import estimate_theoretical_performance
from get_expected_B import get_expected_B
from plot_for_article import plot_for_article

# %% combine results
data, data_lgB, ksis_unique, avg_data, expect_mean_n, trials_number = combine_results(
    bl_force_reload=False)

# %% >>> Theoretical estimates <<<

# %% >> > Plot << <
# print(exp_zeta_ts_over_zeta_sps)
plot_for_article(ksis_unique, avg_data,
                 data_lgB, expect_mean_n, trials_number)
