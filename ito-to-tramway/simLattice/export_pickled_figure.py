
import glob
import logging
import multiprocessing
import os
import pickle
from multiprocessing import Pool, Process

import matplotlib
from tqdm import tqdm

# %matplotlib


def export_pickled_figure(pickle_file, show=False, pdf=True, png=False, **kwargs):
    """
    Load and save a pickled figure in the same place
    """
    # sys.stdout = open(str(os.getpid()) + ".out", "a", buffering=0)
    # sys.stderr = open(str(os.getpid()) + "_error.out", "a", buffering=0)

    print('>>> Plotting ', pickle_file)
    # print(png)
    with open(pickle_file, 'rb') as file:
        fig = pickle.load(file)
    if show:
        fig.show()

    basename = os.path.splitext(pickle_file)[0]
    if pdf:
        pdf_file = basename + '.pdf'
        fig.savefig(pdf_file, bbox_inches='tight', pad_inches=0)

    if png:
        # print('png', png)
        factor = 5
        figsize = fig.get_size_inches()
        fig.set_figwidth(figsize[0] * factor)
        fig.set_figheight(figsize[1] * factor)
        # print('A1')
        try:
            # print(basename + '.png')
            fig.savefig(basename + '.png', pad_inches=0, bbox_inches='tight')  # )

        except:
            logging.warn('Unable to save figure. The file might be open')

    fig = None
    print('<<< Finished plotting ', pickle_file)


def export_all_pickled_figures_in_folder(folder, pdf=False, png=False, **kwargs):
    try:
        cpus = multiprocessing.cpu_count() - 1
    except NotImplementedError:
        cpus = 1

    matplotlib.use('Agg')  # enable for console runs with no displays

    file_list = [f for f in glob.iglob(os.path.join(folder, r"*.pickle"), recursive=False)]
    # print(file_list)

    processes = []
    args_list = []
    print('Using {cpus} threads'.format(cpus=cpus))
    with Pool(cpus) as pool:
        if pdf:
            for f in file_list:
                pool.apply_async(export_pickled_figure, args=(f, ),
                                 kwds={'pdf': True, 'png': False})
        if png:
            for f in file_list:
                pool.apply_async(export_pickled_figure, args=(f, ),
                                 kwds={'png': True, 'pdf': False})

        pool.close()
        pool.join()
        # print('Prosessing ' + f)
        # args_list.apply_async((f,))
        # p = Process(target=export_pickled_figure, args=(f, **kwargs))
        # processes.append(p)
        # p.start()

    # Launch processes
    # with Pool(cpus) as pool:
        # pool.map(export_pickled_figure, args_list)

        # Collect results of the processes
        # for p in tqdm(processes):
        # p.join()


if __name__ == '__main__':
    pickle_file = r"D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\bayes_factors\lattice_diffusivity_obstacles_N=1000_all_trajectories_hexagon_D.pickle"
    # folder = r"D:\Google Drive\git\Stochastic_Integrals_Diffusivity\ito-to-tramway\bayes_factors"
    # folder = r"D:\calculated_data\sim_performance_2D_with_perp_1000_internal_steps\dat"
    folder = r"D:\\"

    # export_pickled_figure(pickle_file)
    export_all_pickled_figures_in_folder(folder, png=True, pdf=True)
