import math

from .tapp import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

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
            fwhm_rt = 9,
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
     # raw_data = read_mzxml(
            # '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
            # instrument_type = 'orbitrap',
            # resolution_ms1 = 70000,
            # resolution_msn = 30000,
            # reference_mz = 200,
            # fwhm_rt = 9,
            # polarity = 'pos',
        # )
     return raw_data

# NOTE: This is not the best design for this function and could be greatly improved.
def plot_mesh(mesh, transform='none', figure=None):
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

def find_scan_indexes(raw_data, peak_candidate):
    # Find min/max scans.
    rts = np.array([scan.retention_time for scan in raw_data.scans])
    scan_idx = np.where((rts >= peak_candidate['roi_min_rt']) & (rts <= peak_candidate['roi_max_rt']))[0]
    return scan_idx

def find_mz_indexes(raw_data, peak_candidate, scan_idx):
    mz_idx = []
    for j in scan_idx:
        scan = raw_data.scans[j]
        mz_i = np.where(np.array(
            (scan.mz >= peak_candidate['roi_min_mz']) &
            (scan.mz <= peak_candidate['roi_max_mz'])))[0]
        mz_idx = mz_idx + [mz_i]
    return mz_idx

def find_raw_points_py(raw_data, scan_idx, mz_idx):
    mzs = []
    rts = []
    intensities = []
    for i in range(0, len(scan_idx)):
        scan = raw_data.scans[scan_idx[i]]
        mzs = mzs + [scan.mz[j] for j in mz_idx[i]]
        intensities = intensities + [scan.intensity[j] for j in mz_idx[i]]
        rts  = np.concatenate([rts, np.repeat(scan.retention_time, len(mz_idx[i]))])
    return (np.array(mzs), np.array(intensities), np.array(rts))

def gaus(x, a, x0, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def gaus2d(X, a, x_0, sigma_x, y_0, sigma_y):
    x = X[0]
    y = X[1]
    return a * np.exp(-0.5 * ((x - x_0) /sigma_x) ** 2 ) * np.exp(-0.5 * ((y - y_0)/sigma_y) ** 2)

def fit_raw_points(mzs, intensities, rts):
    X = np.array([mzs, rts])
    mean_x = sum(X[0] * intensities) / sum(intensities)
    sigma_x = np.sqrt(sum(intensities * (X[0] - mean_x)**2) / sum(intensities))
    mean_y = sum(X[1] * intensities) / sum(intensities)
    sigma_y = np.sqrt(sum(intensities * (X[1] - mean_y)**2) / sum(intensities))
    fitted_parameters, pcov_2d = curve_fit(gaus2d, X, intensities, p0=[max(intensities), mean_x, sigma_x, mean_y, sigma_y])
    return fitted_parameters

def fit(raw_data, peak_candidate):
    scan_idx = find_scan_indexes(raw_data, peak_candidate)
    mz_idx = find_mz_indexes(raw_data, peak_candidate, scan_idx)
    data_points = find_raw_points_py(raw_data, scan_idx, mz_idx)
    fitted_parameters = fit_raw_points(data_points[0], data_points[1], data_points[2])
    return fitted_parameters

def fit2(raw_data, peak_candidate):
    data_points = find_raw_points(
            raw_data,
            peak_candidate['roi_min_mz'],
            peak_candidate['roi_max_mz'],
            peak_candidate['roi_min_rt'],
            peak_candidate['roi_max_rt']
        )
    fitted_parameters = fit_raw_points(data_points.mz, data_points.intensity, data_points.rt)
    return fitted_parameters

def plot_peak_fit(raw_data, peak, fig_mz, fig_rt):
    # PLOTTING
    color = np.random.rand(3,1).flatten()

    data_points = find_raw_points(
            raw_data,
            peak['roi_min_mz'],
            peak['roi_max_mz'],
            peak['roi_min_rt'],
            peak['roi_max_rt']
        )

    rts = data_points.rt
    mzs = data_points.mz
    intensities = data_points.intensity

    # MZ fit plot.
    sort_idx_mz = np.argsort(mzs)
    fitted_intensity_2d_mz = gaus(
            np.array(mzs)[sort_idx_mz],
            peak['fitted_height'],
            peak['fitted_mz'],
            peak['fitted_sigma_mz'],
        )
    plt.figure(fig_mz.number)
    markerline, stemlines, baseline = plt.stem(np.array(mzs)[sort_idx_mz], np.array(intensities)[sort_idx_mz], label='intensities', markerfmt=' ')
    plt.setp(baseline, color=color, alpha=0.5)
    plt.setp(stemlines, color=color, alpha=0.5)
    plt.plot(np.array(mzs)[sort_idx_mz], fitted_intensity_2d_mz, linestyle='--', color=color, label='2d_fitting')
    plt.xlabel('m/z')
    plt.ylabel('Intensity')

    # RT fit plot.
    xic_x = np.unique(rts)
    # xic_y_max = []
    # for x,y in zip(rts, intensities):
        # pass
    sort_idx_rt = np.argsort(xic_x)
    fitted_intensity_2d_rt = gaus(
            np.array(xic_x)[sort_idx_rt],
            peak['fitted_height'],
            peak['fitted_rt'],
            peak['fitted_sigma_rt'],
        )
    plt.figure(fig_rt.number)
    plt.plot(np.array(xic_x)[sort_idx_rt], fitted_intensity_2d_rt, color=color, linestyle='--')
    # plt.plot(xic_x, xic_y_max, label=str(i), linestyle='-', color=color, alpha=0.5)
    plt.xlabel('retention time (s)')
    plt.ylabel('Intensity')

    return

def fit_and_plot(raw_data, peak_candidate):
    data_points = find_raw_points(
            raw_data,
            peak_candidate['roi_min_mz'],
            peak_candidate['roi_max_mz'],
            peak_candidate['roi_min_rt'],
            peak_candidate['roi_max_rt']
        )
    fitted_parameters = fit_raw_points(data_points.mz, data_points.intensity, data_points.rt)
    plot_peak_candidate(data_points, fitted_parameters)
    return fitted_parameters

def find_roi(raw_data, local_max, avg_rt_fwhm=10):
    peak_candidates = []
    for i in range(0, len(local_max)):
        selected_peak = local_max.iloc[i]
        theoretical_sigma_mz = fwhm_at(raw_data, selected_peak['mz']) / (2 * math.sqrt(2 * math.log(2)))
        theoretical_sigma_rt = avg_rt_fwhm / (2 * math.sqrt(2 * math.log(2)))
        tolerance_mz = 3 * theoretical_sigma_mz
        tolerance_rt = 3 * theoretical_sigma_rt
        min_mz = selected_peak['mz'] - tolerance_mz
        max_mz = selected_peak['mz'] + tolerance_mz
        min_rt = selected_peak['rt'] - tolerance_rt
        max_rt = selected_peak['rt'] + tolerance_rt

        peak_candidates = peak_candidates + [{
                'i': selected_peak['i'],
                'j': selected_peak['j'],
                'estimated_mz': selected_peak['mz'],
                'estimated_rt': selected_peak['rt'],
                'estimated_height': selected_peak['intensity'],
                'roi_min_mz': min_mz,
                'roi_max_mz': max_mz,
                'roi_min_rt': min_rt,
                'roi_max_rt': max_rt,
            }]

    return peak_candidates

def profile_peak_fitting(max_peaks=20):
    print("Loading data...")
    raw_data = load_example_data()

    print("Resampling...")
    mesh = resample(raw_data, 5, 5, 0.5, 0.5)

    print("Saving mesh to disk...")
    mesh.save("mesh.dat")

    print("Finding local maxima in mesh...")
    local_max = find_local_max(mesh)
    local_max = pd.DataFrame(local_max)
    local_max.columns = ['i', 'j', 'mz', 'rt', 'intensity']
    local_max = local_max.sort_values('intensity', ascending=False)
    if max_peaks != math.inf:
        local_max = local_max[0:max_peaks]

    peak_candidates = find_roi(raw_data, local_max)
    fitted_parameters = []
    fitted_peaks = []
    for peak_candidate in peak_candidates:
        try:
            # fitted_parameters = fitted_parameters + [fit(raw_data, peak_candidate)]
            fitted_parameters = fitted_parameters + [fit2(raw_data, peak_candidate)]
            peak = peak_candidate
            peak['fitted_height'] = fitted_parameters[0]
            peak['fitted_mz'] = fitted_parameters[1]
            peak['fitted_sigma_mz'] = fitted_parameters[2]
            peak['fitted_rt'] = fitted_parameters[3]
            peak['fitted_sigma_rt'] = fitted_parameters[4]
            fitted_peaks = fitted_peaks + [peak]
        except:
            # print("Couldn't fit peak candidate: {}".format(peak_candidate))
            pass

    return fitted_peaks

def example_pipeline(show_mesh_plot=False, show_plot_fit=True, silent=True, max_peaks=15):
    if show_plot_fit or show_mesh_plot:
        plt.style.use('dark_background')
        plt.ion()
        plt.show()

    print("Loading data...")
    raw_data = read_mzxml(
        '/data/toydata/toy_data.mzXML',
        instrument_type = 'orbitrap',
        resolution_ms1 = 75000,
        resolution_msn = 30000,
        reference_mz = 200,
        fwhm_rt = 9,
        polarity = 'pos',
        min_mz = 801,
        max_mz = 803,
        min_rt = 2808,
        max_rt = 2928,
    )

    print("Resampling...")
    mesh = resample(raw_data, 10, 10, 0.5, 0.5)

    print("Saving mesh to disk...")
    mesh.save("mesh.dat")

    print("Finding local maxima in mesh...")
    local_max = find_local_max(mesh)
    local_max = pd.DataFrame(local_max)
    local_max.columns = ['i', 'j', 'mz', 'rt', 'intensity']
    local_max = local_max.sort_values('intensity', ascending=False)
    if max_peaks != math.inf:
        local_max = local_max[0:max_peaks]

    if show_mesh_plot:
        print("Plotting mesh...")
        mesh_plot = plot_mesh(mesh, transform='sqrt')

        print("Plotting local maxima...")
        mesh_plot['img_plot'].scatter(local_max['i'], local_max['j'], color='aqua', s=5, marker="s", alpha=0.9)

    print("Fitting peaks...")
    peak_candidates = find_roi(raw_data, local_max)
    fitted_peaks = []
    if show_plot_fit:
        fig_mz = plt.figure()
        fig_rt = plt.figure()
    for peak_candidate in peak_candidates:
        try:
            # fitted_parameters = fitted_parameters + [fit(raw_data, peak_candidate)]
            fitted_parameters = fit2(raw_data, peak_candidate)
            peak = peak_candidate
            peak['fitted_height'] = fitted_parameters[0]
            peak['fitted_mz'] = fitted_parameters[1]
            peak['fitted_sigma_mz'] = fitted_parameters[2]
            peak['fitted_rt'] = fitted_parameters[3]
            peak['fitted_sigma_rt'] = fitted_parameters[4]
            fitted_peaks = fitted_peaks + [peak]
            if show_plot_fit:
                plot_peak_fit(raw_data, peak, fig_mz, fig_rt)
        except Exception as e:
            print(e)
            pass

    # fitted_peaks = fit_peaks(raw_data, local_max, show_plot_fit=show_plot_fit)
    # fitted_peaks_tuple = [tuple(fitted_peaks.iloc[row]) for row in range(0, fitted_peaks.shape[0])]
    # print("Saving fitted peaks to disk...")
    # tapp.save_fitted_peaks(list(fitted_peaks_tuple), "fitted_peaks.bpks")
    # pd.DataFrame(fitted_peaks).to_csv('fitted_peaks.csv')

    return (raw_data, mesh, local_max, fitted_peaks)

RawData.tic = tic
