"""
This interactve python file combines data for the article statistical performance plot and then makes the plot.
It doesn't perform the inference itself. It should have already been performed by launching the `start_me.py`.
"""

try:
    has_run
except NameError:
    %matplotlib
    %load_ext autoreload
    %autoreload 2

    has_run = 1
else:
    print("Graphic interface NOT re-initialized")


from combine_results import combine_results
from estimate_theoretical_performance import estimate_theoretical_performance
import numpy as np
from plot_for_article import plot_for_article

# %% combine results
data, ksis_unique, avg_data, expect_mean_n, trials_number = combine_results(
    bl_force_reload=True)

# %% >>> Theoretical estimates <<<
exp_zeta_ts_over_zeta_sps = estimate_theoretical_performance(expect_mean_n)

# %% >> > Plot << <
# print(exp_zeta_ts_over_zeta_sps)
plot_for_article(ksis_unique, avg_data,
                 exp_zeta_ts_over_zeta_sps, expect_mean_n, trials_number)
