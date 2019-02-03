import math

from .tapp import *
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import time

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
            # '/data/ftp_data/150210_11_01.mzXML',
            # instrument_type = 'orbitrap',
            # resolution_ms1 = 75000,
            # resolution_msn = 30000,
            # reference_mz = 200,
            # polarity = 'pos',
            # min_mz = 400,
            # max_mz = 1000,
            # min_rt = 2000,
            # max_rt = 4000,
        # )
     raw_data = read_mzxml(
            '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
            instrument_type = 'orbitrap',
            resolution_ms1 = 70000,
            resolution_msn = 30000,
            reference_mz = 200,
            polarity = 'pos',
        )
     return raw_data

# NOTE: This is not the best design for this function and could be greatly improved.
def plot_mesh(mesh, transform='none', figure=None):
    plt.style.use('dark_background')
    plt.ion()
    plt.show()

    if figure is None:
        figure = plt.figure()

    img = mesh.matrix
    img = np.reshape(img, (mesh.m, mesh.n))
    bins_rt = mesh.bins_rt
    bins_mz = mesh.bins_mz
    num_bins_mz = len(bins_mz)
    num_bins_rt = len(bins_rt)
    min_mz = np.array(bins_mz).min()
    max_mz = np.array(bins_mz).max()
    min_rt = np.array(bins_rt).min()
    max_rt = np.array(bins_rt).max()
    if transform == 'sqrt':
        img = np.sqrt(img)
    elif transform == 'log':
        img = np.log(img + 0.00001)

    plt.figure(figure.number)
    plt.clf()
    gs = gridspec.GridSpec(5,5)
    mz_plot = plt.subplot(gs[0, :-1])
    mz_plot.clear()
    mz_plot.plot(img.sum(axis=0))
    mz_plot.margins(x=0)
    mz_plot.set_xticks([]) 
    mz_plot.set_ylabel("Intensity")

    rt_plot = plt.subplot(gs[1:, -1])
    rt_plot.plot(img.sum(axis=1), bins_rt)  
    rt_plot.margins(y=0)
    rt_plot.set_yticks([]) 
    rt_plot.set_xlabel("Intensity")

    img_plot = plt.subplot(gs[1:, :-1])
    # img_plot.pcolormesh(mesh.bins_mz, mesh.bins_rt, img)
    img_plot.imshow(img, aspect='auto', origin="lower")
    
    # This only approximates thew ticks assuming we are on a small region, on a
    # warped grid this doesn't work for a large mz range.
    delta_mz = bins_mz[1] - bins_mz[0]
    delta_rt = bins_rt[1] - bins_rt[0]
    locs_x = img_plot.get_xticks()
    labels_x = []
    for loc in locs_x:
        if loc < 0 or loc >= num_bins_mz:
            labels_x = labels_x + [""]
        else:
            labels_x = labels_x + ["{0:.2f}".format(min_mz + delta_mz * loc)]

    locs_y = img_plot.get_yticks()
    labels_y = []
    for loc in locs_y:
        if loc < 0 or loc >= num_bins_rt:
            labels_y = labels_y + [""]
        else:
            labels_y = labels_y + ["{0:.2f}".format(min_rt + delta_rt * loc)]
    
    img_plot.set_xticklabels(labels_x)
    img_plot.set_yticklabels(labels_y)
    img_plot.set_xlabel("m/z")
    img_plot.set_ylabel("retention time (s)")

    return({
        "img_plot": img_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })

def example_pipeline(show_mesh_plot=False, show_plot_fit=False, silent=True):
    start = time.time()
    print("Loading data...")
    raw_data = load_example_data()
    end = time.time()
    print("> Finished in: {}".format(end - start))

    print("Resampling...")
    start = time.time()
    n, m = calculate_dimensions(raw_data, 9, 10, 10)
    print("Estimated memory consumption of the [{0}x{1}] grid: {2:.2f} (MB)".format(n, m, n * m /1024/1024 * 8))
    mesh = resample(raw_data, 9, 10, 10)
    end = time.time()
    print("> Finished in: {}".format(end - start))

    print("Saving mesh to disk...")
    start = time.time()
    mesh.save("mesh.dat")
    end = time.time()
    print("> Finished in: {}".format(end - start))

    print("Finding local maxima in mesh...")
    start = time.time()
    local_max = find_local_max(mesh)
    local_max = pd.DataFrame(local_max)
    local_max.columns = ['i', 'j', 'mz', 'rt', 'intensity']
    local_max = local_max.sort_values('intensity', ascending=False)
    end = time.time()
    print("> Finished in: {}".format(end - start))

    if show_mesh_plot:
        print("Plotting mesh...")
        mesh_plot = plot_mesh(mesh)

        print("Plotting local maxima...")
        mesh_plot['img_plot'].scatter(local_max['i'], local_max['j'], color='aqua', s=5, marker="s", alpha=0.9)

    # print("Fitting the top 10 peaks...")
    print("Fitting peaks...")
    start = time.time()
    def sigma_at_mz(mz, fwhm_ref, mz_ref):
        return fwhm_ref * (mz/mz_ref) ** 1.5  # NOTE: Orbitrap only

    # FIXME: The plotting should be independant of the fitting loop. This it is
    # terrible design.
    if show_plot_fit:
        fig_2 = plt.figure()
        fig_3 = plt.figure()

    retention_times = [scan.retention_time for scan in raw_data.scans]
    fitted_peaks = []
    for i in range(0, len(local_max)):
    # for i in range(0, 10):
        # Show log message each 10%
        # if i % (len(local_max)/10) == 0 and not silent:
        # if i % (len(local_max)/10) == 0:
        print("peak {0} out of {1}".format(i, local_max.shape[0]))
        selected_peak = local_max.iloc[i]
        # print(selected_peak)
        theoretical_sigma_mz = fwhm_at(raw_data, selected_peak['mz']) / 2.355 # FIXME: This is a rough approximation.
        theoretical_sigma_rt = 6
        tolerance_mz = 3 * theoretical_sigma_mz
        # print(theoretical_sigma_mz)
        # print(tolerance_mz)
        # print("mz_min: {}".format(selected_peak['mz'] - tolerance_mz))
        # print("mz_max: {}".format(selected_peak['mz'] + tolerance_mz))
        min_mz = selected_peak['mz'] - tolerance_mz
        max_mz = selected_peak['mz'] + tolerance_mz

        print("> Finding the closest scans:")
        start_internal = time.time()
        # Find the required scans.
        closest_scan_index = np.abs(retention_times - selected_peak['rt']).argmin()

        # The closest N scans are used for fitting.
        scan_indices = [idx for idx in list(range(closest_scan_index - 5, closest_scan_index + 5)) if idx >= 0 and idx < len(raw_data.scans)]
        end_internal = time.time()
        print("> Finished in: {}".format(end_internal - start_internal))

        if len(scan_indices) < 3:
            continue

        # NOTE: Should we store the total number of non zero scans? this is the thing that we should use for filtering.
        # print(scan_indices)
        print("> Finding the mz points:")
        start_internal = time.time()
        mzs = []
        rts = []
        intensities = []
        if show_plot_fit:
            xic_x = []
            xic_y_total = []
            xic_y_max = []
        for idx in scan_indices:
            scan = raw_data.scans[idx]

            if show_plot_fit:
                total_intensity = 0
                max_intensity = 0
            # for i in range(0, scan.num_points):
                # # Find the mz values within  the desired range.
                # if scan.mz[i] > max_mz:
                    # break
                # if scan.mz[i] < min_mz:
                    # continue

                # mzs = mzs + [scan.mz[i]]
                # rts = rts + [scan.retention_time]
                # intensities = intensities + [scan.intensity[i]]
                # if show_plot_fit:
                    # total_intensity = total_intensity + scan.intensity[i]
                    # if scan.intensity[i] > max_intensity:
                        # max_intensity = scan.intensity[i]

            # DEBUG: is this more efficient?
            mz_index = np.argwhere(np.array((scan.mz > min_mz) & (scan.mz < max_mz))).flatten()
            # mzs = np.concatenate([np.array([1,2,3]), np.array([4,5,6]]))
            mzs = np.concatenate([mzs, np.array(scan.mz)[mz_index]])
            intensities = np.concatenate([intensities, np.array(scan.intensity)[mz_index]])
            rts  = np.concatenate([rts, np.repeat(scan.retention_time, len(mz_index))])
            # print(mzs)
            # print(intensities)
            # print(rts)

            if show_plot_fit:
                xic_x = xic_x + [scan.retention_time]
                xic_y_total = xic_y_total + [total_intensity]
                xic_y_max = xic_y_max + [max_intensity]
        end_internal = time.time()
        print("> Finished in: {}".format(end_internal - start_internal))

        def gaus(x, a, x0, sigma):
                return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        def gaus2d(X, a, x_0, sigma_x, y_0, sigma_y):
            x = X[0]
            y = X[1]
            return a * np.exp(-0.5 * ((x - x_0) /sigma_x) ** 2 ) * np.exp(-0.5 * ((y - y_0)/sigma_y) ** 2)

        print("> Fitting the 2D Gaussian:")
        start_internal = time.time()
        try:
            X = np.array([mzs, rts])
            mean_x = sum(X[0] * intensities) / sum(intensities)
            sigma_x = np.sqrt(sum(intensities * (X[0] - mean_x)**2) / sum(intensities))
            mean_y = sum(X[1] * intensities) / sum(intensities)
            sigma_y = np.sqrt(sum(intensities * (X[1] - mean_y)**2) / sum(intensities))
            popt_2d, pcov_2d = curve_fit(gaus2d, X, intensities, p0=[max(intensities), mean_x, sigma_x, mean_y, sigma_y])
            # print(popt_2d)
        except:
            if not silent:
                print("error when fitting 2D gaussian on peak: {}".format(i))
            continue
        end_internal = time.time()
        print("> Finished in: {}".format(end_internal - start_internal))

        # # # Discard peaks where the fit is not conforming with the theoretical distributions.
        # if popt_2d[2] == 0 or popt_2d[4] == 0:
            # continue

        # if popt_2d[2] > 3 * theoretical_sigma_mz or popt_2d[4] > 3 * theoretical_sigma_rt:
            # continue

        # Store the fitted parameters.
        fitted_peaks = fitted_peaks + [{
                'smoothed_mz': selected_peak['mz'],
                'smoothed_rt': selected_peak['rt'],
                'estimated_mz': mean_x,
                'estimated_rt': mean_y,
                'estimated_sigma_mz': sigma_x,
                'estimated_sigma_rt': sigma_y,
                'fitted_mz': popt_2d[1],
                'fitted_rt': popt_2d[3],
                'fitted_height': popt_2d[0],
                'fitted_sigma_mz': popt_2d[2],
                'fitted_sigma_rt': popt_2d[4],
                'fitted_total_intensity': np.array(intensities).sum(),
                'roi_mz_min': min_mz,
                'roi_mz_max': max_mz,
                'roi_rt_min': rts[0],
                'roi_rt_max': rts[-1],
            }]

        if show_plot_fit:
            # PLOTTING
            color = np.random.rand(3,1).flatten()

            # MZ fit plot.
            sort_idx_mz = np.argsort(mzs)
            fitted_intensity_2d_mz = gaus(np.array(mzs)[sort_idx_mz], popt_2d[0],popt_2d[1],popt_2d[2])
            plt.figure(fig_2.number)
            markerline, stemlines, baseline = plt.stem(np.array(mzs)[sort_idx_mz], np.array(intensities)[sort_idx_mz], label='intensities', markerfmt=' ')
            plt.setp(baseline, color=color, alpha=0.5)
            plt.setp(stemlines, color=color, alpha=0.5)
            plt.plot(np.array(mzs)[sort_idx_mz], fitted_intensity_2d_mz, linestyle='--', color=color, label='2d_fitting')
            plt.xlabel('m/z')
            plt.ylabel('Intensity')

            # RT fit plot.
            sort_idx_rt = np.argsort(xic_x)
            fitted_intensity_2d_rt = gaus(np.array(xic_x)[sort_idx_rt], popt_2d[0],popt_2d[3],popt_2d[4])
            plt.figure(fig_3.number)
            plt.plot(np.array(xic_x)[sort_idx_rt], fitted_intensity_2d_rt, color=color, linestyle='--')
            plt.plot(xic_x, xic_y_max, label=str(i), linestyle='-', color=color, alpha=0.5)
            plt.xlabel('retention time (s)')
            plt.ylabel('Intensity')

            # # # Contour plot
            # # x = grid['bins_mz']
            # # y = grid['bins_rt']
            # # x, y = np.meshgrid(x, y)
            # # Z = gaus2d((x,y), *popt_2d) 
            # # x = range(0, len(grid['bins_mz']))
            # # y = range(0, len(grid['bins_rt']))
            # # x, y = np.meshgrid(x, y)
            # # grid_plot['grid_plot'].contour(x, y, Z, levels=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10], colors=[color], alpha=0.5)
            # # grid_plot['grid_plot'].scatter(selected_peak['i'], selected_peak['j'], color=color, marker='x')
            # # grid_plot['grid_plot'].scatter(popt_2d[1]/(np.array(grid['bins_mz']).max() - np.array(grid['bins_mz']).min()), popt_2d[3] - np.array(grid['bins_rt']).min(), color=color, marker='P')
            # # print(np.abs(popt_2d[1] - selected_peak['mz']))

    fitted_peaks = pd.DataFrame(fitted_peaks)
    end = time.time()
    print("> Finished in: {}".format(end - start))

    print("Saving fitted peaks to disk...")
    start = time.time()
    pd.DataFrame(fitted_peaks).to_csv('fitted_peaks.csv')
    end = time.time()
    print("> Finished in: {}".format(end - start))

    return (raw_data, mesh, local_max, fitted_peaks)

RawData.tic = tic
