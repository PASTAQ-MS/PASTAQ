import math

from .tapp import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# TODO(alex): Write documentation.
def tic(raw_data, min_rt = -math.inf, max_rt = math.inf):
    rt = []
    intensity = []
    for i in range(0, len(raw_data.scans)):
        if (raw_data.scans[i].retention_time < min_rt or
            raw_data.scans[i].retention_time > max_rt):
            continue
        sum = 0
        for j in range(0, raw_data.scans[i].num_points):
            sum = sum + raw_data.scans[i].intensity[j]
        intensity = intensity + [sum]
        rt = rt + [raw_data.scans[i].retention_time]

    return (rt, intensity)

def load_example_data():
     raw_data = read_mzxml(
            '/data/toydata/toy_data.mzXML',
            instrument_type = 'orbitrap',
            resolution_ms1 = 75000,
            resolution_msn = 30000,
            reference_mz = 200,
            polarity = 'pos',
            min_mz = 801,
            max_mz = 803,
            min_rt = 2808,
            max_rt = 2928,
        )
     # raw_data = read_mzxml(
            # '/data/Spike-In QEx UKE/HD5YD15ED_DDA_R1.mzXML',
            # instrument_type = 'orbitrap',
            # resolution_ms1 = 75000,
            # resolution_msn = 30000,
            # reference_mz = 200,
            # polarity = 'pos',
            # min_mz = 801,
            # max_mz = 803,
            # min_rt = 2808,
            # max_rt = 2928,
        # )
     # raw_data = read_mzxml(
            # '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
            # instrument_type = 'orbitrap',
            # resolution_ms1 = 70000,
            # resolution_msn = 30000,
            # reference_mz = 200,
            # polarity = 'pos',
        # )
     return raw_data

def example_pipeline():
    print("Loading data")
    raw_data = load_example_data()
    n, m = calculate_dimensions(raw_data, 9, 10, 10)
    print("Estimated memory consumption of the [{0}x{1}] grid: {2:.2f} (MB)".format(n, m, n * m /1024/1024))
    mesh = resample(raw_data, 9, 10, 10)
    return (raw_data, mesh)

# NOTE: This is not the best design for this function and could be greatly improved.
def plot_grid(mesh, transform='none', figure=None):
    plt.style.use('dark_background')
    plt.ion()
    plt.show()

    if figure is None:
        figure = plt.figure()

    grid = mesh.matrix
    grid = np.reshape(grid, (mesh.m, mesh.n))
    bins_rt = mesh.bins_rt
    bins_mz = mesh.bins_mz
    num_bins_mz = len(bins_mz)
    num_bins_rt = len(bins_rt)
    min_mz = np.array(bins_mz).min()
    max_mz = np.array(bins_mz).max()
    min_rt = np.array(bins_rt).min()
    max_rt = np.array(bins_rt).max()
    if transform == 'sqrt':
        grid = np.sqrt(grid)
    elif transform == 'log':
        grid = np.log(grid + 0.00001)

    plt.figure(figure.number)
    plt.clf()
    gs = gridspec.GridSpec(5,5)
    mz_plot = plt.subplot(gs[0, :-1])
    mz_plot.clear()
    mz_plot.plot(grid.sum(axis=0))
    mz_plot.margins(x=0)
    mz_plot.set_xticks([]) 
    mz_plot.set_ylabel("Intensity")

    rt_plot = plt.subplot(gs[1:, -1])
    rt_plot.plot(grid.sum(axis=1), bins_rt)  
    rt_plot.margins(y=0)
    rt_plot.set_yticks([]) 
    rt_plot.set_xlabel("Intensity")

    grid_plot = plt.subplot(gs[1:, :-1])
    grid_plot.pcolormesh(mesh.bins_mz, mesh.bins_rt, grid)

    grid_plot.set_xlabel("m/z")
    grid_plot.set_ylabel("retention time (s)")

    return({
        "grid_plot": grid_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })

RawData.tic = tic
