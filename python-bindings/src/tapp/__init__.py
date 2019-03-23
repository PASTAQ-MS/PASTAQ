import math

from .tapp import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

# TODO(alex): Write documentation.


def tic(raw_data, min_rt=-math.inf, max_rt=math.inf):
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
        instrument_type='orbitrap',
        resolution_ms1=75000,
        resolution_msn=30000,
        reference_mz=200,
        fwhm_rt=9,
        polarity='pos',
        min_mz=801,
        max_mz=803,
        min_rt=2808,
        max_rt=2928,
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
    gs = gridspec.GridSpec(5, 5)
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
    scan_idx = np.where((rts >= peak_candidate['roi_min_rt']) & (
        rts <= peak_candidate['roi_max_rt']))[0]
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
        rts = np.concatenate(
            [rts, np.repeat(scan.retention_time, len(mz_idx[i]))])
    return (np.array(mzs), np.array(intensities), np.array(rts))


def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def gaus2d(X, a, x_0, sigma_x, y_0, sigma_y):
    x = X[0]
    y = X[1]
    return a * np.exp(-0.5 * ((x - x_0) / sigma_x) ** 2) * np.exp(-0.5 * ((y - y_0)/sigma_y) ** 2)


def fit_curvefit(mzs, intensities, rts):
    X = np.array([mzs, rts])
    mean_x = sum(X[0] * intensities) / sum(intensities)
    sigma_x = np.sqrt(sum(intensities * (X[0] - mean_x)**2) / sum(intensities))
    mean_y = sum(X[1] * intensities) / sum(intensities)
    sigma_y = np.sqrt(sum(intensities * (X[1] - mean_y)**2) / sum(intensities))
    fitted_parameters, pcov_2d = curve_fit(gaus2d, X, intensities, p0=[
                                           max(intensities), mean_x, sigma_x, mean_y, sigma_y])
    return fitted_parameters

# NOTE: Testing different fitting methods


def generate_gaussian():
    x = np.linspace(801.38, 801.42, 100)
    y = gaus(x, 1, 801.40, 0.0025)
    # x = x - x.mean()
    return (x, y)

# NOTE: Testing different fitting methods


def fit_inle(x, y):
    mean_x = np.sum(x * y) / np.sum(y)
    sigma_x = np.sqrt(sum(y * (x - mean_x)**2) / sum(y))
    fitted_parameters, cov = curve_fit(gaus, x, y, p0=[max(y), mean_x, sigma_x])
    return fitted_parameters


def fit_caruana(x, y):
    x_mean = x.mean()
    x = x - x_mean
    X = np.array(
        [
            [len(x), np.array(x).sum(), np.power(x, 2).sum()],
            [np.array(x).sum(), np.power(x, 2).sum(), np.power(x, 3).sum()],
            [np.power(x, 2).sum(), np.power(x, 3).sum(), np.power(x, 4).sum()]
        ],
    )
    Y = np.array([
        np.log(y).sum(),
        (x * np.log(y)).sum(),
        (np.power(x, 2) * np.log(y)).sum()
    ])
    a, b, c = np.linalg.solve(X, Y)
    mean = -b / (2 * c) + x_mean
    sigma = np.sqrt(-1 / (2 * c))
    height = np.exp(a - (b ** 2) / (4 * c))
    # print(np.allclose(np.dot(X, A), Y))
    return np.array([height, mean, sigma])


def fit_guos(x, y):
    x_mean = x.mean()
    x = x - x_mean
    X = np.array(
        [
            [
                np.power(y, 2).sum(),
                (x * np.power(y, 2)).sum(),
                (np.power(x, 2) * np.power(y, 2)).sum(),
            ],
            [
                (x * np.power(y, 2)).sum(),
                (np.power(x, 2) * np.power(y, 2)).sum(),
                (np.power(x, 3) * np.power(y, 2)).sum(),
            ],
            [
                (np.power(x, 2) * np.power(y, 2)).sum(),
                (np.power(x, 3) * np.power(y, 2)).sum(),
                (np.power(x, 4) * np.power(y, 2)).sum(),
            ],
        ],
    )
    Y = np.array([
        (np.power(y, 2) * np.log(y)).sum(),
        (np.power(y, 2) * x * np.log(y)).sum(),
        (np.power(y, 2) * np.power(x, 2) * np.log(y)).sum()
    ])
    a, b, c = np.linalg.solve(X, Y)
    mean = -b / (2 * c) + x_mean
    sigma = np.sqrt(-1 / (2 * c))
    height = np.exp(a - (b ** 2) / (4 * c))
    # print(np.allclose(np.dot(X, A), Y))
    return np.array([height, mean, sigma])


def fit_guos_2d(x, y, z):
    x_mean = x.mean()
    x = x - x_mean
    y_mean = y.mean()
    y = y - y_mean

    z_2 = np.power(z, 2)
    x_2 = np.power(x, 2)
    x_3 = np.power(x, 3)
    x_4 = np.power(x, 4)
    y_2 = np.power(y, 2)
    y_3 = np.power(y, 3)
    y_4 = np.power(y, 4)

    X = np.array(
        [
            [
                z_2.sum(),
                (x * z_2).sum(),
                (x_2 * z_2).sum(),
                (y * z_2).sum(),
                (y_2 * z_2).sum(),
            ],
            [
                (x * z_2).sum(),
                (x_2 * z_2).sum(),
                (x_3 * z_2).sum(),
                (x * y * z_2).sum(),
                (x * y_2 * z_2).sum(),
            ],
            [
                (x_2 * z_2).sum(),
                (x_3 * z_2).sum(),
                (x_4 * z_2).sum(),
                (x_2 * y * z_2).sum(),
                (x_2 * y_2 * z_2).sum(),
            ],
            [
                (y * z_2).sum(),
                (x * y * z_2).sum(),
                (x_2 * y * z_2).sum(),
                (y_2 * z_2).sum(),
                (y_3 * z_2).sum(),
            ],
            [
                (y_2 * z_2).sum(),
                (x * y_2 * z_2).sum(),
                (x_2 * y_2 * z_2).sum(),
                (y_3 * z_2).sum(),
                (y_4 * z_2).sum(),
            ],
        ],
    )
    Y = np.array([
        (z_2 * np.log(z)).sum(),
        (z_2 * x * np.log(z)).sum(),
        (z_2 * x_2 * np.log(z)).sum(),
        (z_2 * y * np.log(z)).sum(),
        (z_2 * y_2 * np.log(z)).sum(),
    ])
    # a, b, c, d, e = np.linalg.solve(X, Y)
    beta = np.linalg.lstsq(X, Y)
    a, b, c, d, e = beta[0]

    sigma_mz = np.sqrt(1/(-2 * c))
    mz = b / (-2 * c) + x_mean
    sigma_rt = np.sqrt(1/(-2 * e))
    rt = d / (-2 * e) + y_mean
    height = np.exp(a - ((b ** 2) / (4 * c)) - ((d ** 2) / (4 * e)))

    # print(np.allclose(np.dot(X, A), Y))
    return np.array([height, mz, sigma_mz, rt, sigma_rt])


def fit_guos_2d_from_peak(peak):
    x_mean = peak.raw_roi_mz
    y_mean = peak.raw_roi_rt
    X = np.array(
        [
            [
                peak.a_0_0(),
                peak.a_0_1(),
                peak.a_0_2(),
                peak.a_0_3(),
                peak.a_0_4(),
            ],
            [
                peak.a_1_0(),
                peak.a_1_1(),
                peak.a_1_2(),
                peak.a_1_3(),
                peak.a_1_4(),
            ],
            [
                peak.a_2_0(),
                peak.a_2_1(),
                peak.a_2_2(),
                peak.a_2_3(),
                peak.a_2_4(),
            ],
            [
                peak.a_3_0(),
                peak.a_3_1(),
                peak.a_3_2(),
                peak.a_3_3(),
                peak.a_3_4(),
            ],
            [
                peak.a_4_0(),
                peak.a_4_1(),
                peak.a_4_2(),
                peak.a_4_3(),
                peak.a_4_4(),
            ],
        ],
    )
    Y = np.array([
        peak.c_0(),
        peak.c_1(),
        peak.c_2(),
        peak.c_3(),
        peak.c_4(),
    ])
    # a, b, c, d, e = np.linalg.solve(X, Y)
    beta = np.linalg.lstsq(X, Y, rcond=None)
    a, b, c, d, e = beta[0]

    # if c < 0:
        # sigma_mz = np.sqrt(1/(-2 * c))
        # mz = b / (-2 * c) + x_mean
    # else:
        # sigma_mz = np.nan
        # mz = np.nan
    # if e < 0:
        # sigma_rt = np.sqrt(1/(-2 * e))
        # rt = d / (-2 * e) + y_mean
    # else:
        # sigma_rt = np.nan
        # rt = np.nan

    # if c != 0 and e != 0:
        # height = np.exp(a - ((b ** 2) / (4 * c)) - ((d ** 2) / (4 * e)))
    # else:
        # height = np.nan

    sigma_mz = np.sqrt(1/(-2 * c))
    mz = b / (-2 * c) + x_mean
    sigma_rt = np.sqrt(1/(-2 * e))
    rt = d / (-2 * e) + y_mean
    height = np.exp(a - ((b ** 2) / (4 * c)) - ((d ** 2) / (4 * e)))

    # print(np.allclose(np.dot(X, A), Y))
    return np.array([height, mz, sigma_mz, rt, sigma_rt])


def test_gaus_fit():
    x, y = generate_gaussian()
    parameters_inle = fit_inle(x, y)
    parameters_guos = fit_guos(x, y)
    parameters_caruana = fit_caruana(x, y)
    print("parameters_inle:", parameters_inle)
    print("parameters_guos:", parameters_guos)
    print("parameters_caruana:", parameters_caruana)
    plt.style.use('dark_background')
    plt.ion()
    plt.show()
    fig = plt.figure()
    plt.scatter(x, y, label='Raw data')
    plt.plot(x, gaus(x, *parameters_inle), label='curve_fit',
             linestyle='--', color='crimson')
    plt.plot(x, gaus(x, *parameters_guos), label='guos')
    plt.plot(x, gaus(x, *parameters_caruana),
             label='caruana', linestyle=':', color='aqua')
    plt.legend(loc='upper left')


def fit(raw_data, peak_candidate):
    scan_idx = find_scan_indexes(raw_data, peak_candidate)
    mz_idx = find_mz_indexes(raw_data, peak_candidate, scan_idx)
    data_points = find_raw_points_py(raw_data, scan_idx, mz_idx)
    fitted_parameters = fit_curvefit(
        data_points[0], data_points[1], data_points[2])
    return fitted_parameters


def fit2(raw_data, peak_candidate):
    data_points = find_raw_points(
        raw_data,
        peak_candidate['roi_min_mz'],
        peak_candidate['roi_max_mz'],
        peak_candidate['roi_min_rt'],
        peak_candidate['roi_max_rt']
    )
    fitted_parameters = fit_curvefit(
        data_points.mz, data_points.intensity, data_points.rt)
    return fitted_parameters


def fit3(raw_data, peak_candidate):
    data_points = find_raw_points(
        raw_data,
        peak_candidate['roi_min_mz'],
        peak_candidate['roi_max_mz'],
        peak_candidate['roi_min_rt'],
        peak_candidate['roi_max_rt']
    )
    return fit_guos_2d(
        np.array(data_points.mz),
        np.array(data_points.rt),
        np.array(data_points.intensity))


def fit_raw_weighted_estimate(raw_data, peak_candidate):
    data_points = find_raw_points(
        raw_data,
        peak_candidate['roi_min_mz'],
        peak_candidate['roi_max_mz'],
        peak_candidate['roi_min_rt'],
        peak_candidate['roi_max_rt']
    )
    mzs = np.array(data_points.mz)
    rts = np.array(data_points.rt)
    intensities = np.array(data_points.intensity)

    mean_x = sum(mzs * intensities) / sum(intensities)
    sigma_x = np.sqrt(sum(intensities * (mzs - mean_x)**2) / sum(intensities))
    mean_y = sum(rts * intensities) / sum(intensities)
    sigma_y = np.sqrt(sum(intensities * (rts - mean_y)**2) / sum(intensities))

    return np.array([intensities.max(), mean_x, sigma_x, mean_y, sigma_y])


def plot_peak_fit(raw_data, peak, fig_mz, fig_rt):
    # PLOTTING
    color = np.random.rand(3, 1).flatten()

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
    markerline, stemlines, baseline = plt.stem(np.array(mzs)[sort_idx_mz], np.array(
        intensities)[sort_idx_mz], label='intensities', markerfmt=' ')
    plt.setp(baseline, color=color, alpha=0.5)
    plt.setp(stemlines, color=color, alpha=0.5)
    plt.plot(np.array(mzs)[sort_idx_mz], fitted_intensity_2d_mz,
             linestyle='--', color=color, label='2d_fitting')
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
    plt.plot(np.array(xic_x)[sort_idx_rt],
             fitted_intensity_2d_rt, color=color, linestyle='--')
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
    fitted_parameters = fit_curvefit(
        data_points.mz, data_points.intensity, data_points.rt)
    plot_peak_candidate(data_points, fitted_parameters)
    return fitted_parameters


def find_roi(raw_data, local_max, avg_rt_fwhm=10):
    peak_candidates = []
    for i in range(0, len(local_max)):
        selected_peak = local_max.iloc[i]
        theoretical_sigma_mz = fwhm_at(
            raw_data, selected_peak['mz']) / (2 * math.sqrt(2 * math.log(2)))
        theoretical_sigma_rt = avg_rt_fwhm / (2 * math.sqrt(2 * math.log(2)))
        tolerance_mz = 2.5 * theoretical_sigma_mz
        tolerance_rt = 2.5 * theoretical_sigma_rt
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


def profile_resample():
    # raw_data = read_mzxml(
        # '/data/toydata/toy_data.mzXML',
        # instrument_type = 'orbitrap',
        # resolution_ms1 = 75000,
        # resolution_msn = 30000,
        # reference_mz = 200,
        # fwhm_rt = 9,
        # polarity = 'pos',
        # min_mz = 801,
        # max_mz = 803,
        # min_rt = 2808,
        # max_rt = 2928,
    # )
    raw_data = read_mzxml(
        '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
        instrument_type='orbitrap',
        resolution_ms1=75000,
        resolution_msn=30000,
        reference_mz=200,
        fwhm_rt=9,
        polarity='pos',
        # min_mz = 200,
        # max_mz = 800,
        # min_rt = 0,
        # max_rt = 1000,
    )

    mesh = resample(raw_data, 5, 5, 0.5, 0.5)


def profile_peak_fitting(max_peaks=20):
    print("Loading data...")
    # raw_data = read_mzxml(
    # '/data/qatar/17122018/mzXML/Acute2U_3001.mzXML',
    # instrument_type = 'orbitrap',
    # resolution_ms1 = 70000,
    # resolution_msn = 30000,
    # reference_mz = 200,
    # fwhm_rt = 9,
    # polarity = 'pos',
    # )
    raw_data = read_mzxml(
        '/data/toydata/toy_data.mzXML',
        instrument_type='orbitrap',
        resolution_ms1=75000,
        resolution_msn=30000,
        reference_mz=200,
        fwhm_rt=9,
        polarity='pos',
        min_mz=801,
        max_mz=803,
        min_rt=2808,
        max_rt=2928,
    )

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

    peak_candidates = find_roi(raw_data, local_max, 9)
    fitted_parameters = []
    fitted_peaks = []
    for peak_candidate in peak_candidates:
        try:
            # fitted_parameters = fitted_parameters + [fit(raw_data, peak_candidate)]
            # fitted_parameters = fitted_parameters + [fit2(raw_data, peak_candidate)]
            fitted_parameters = fitted_parameters + \
                [fit3(raw_data, peak_candidate)]
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
    # raw_data = read_mzxml(
    # '/data/toydata/toy_data.mzXML',
    # instrument_type = 'orbitrap',
    # resolution_ms1 = 75000,
    # resolution_msn = 30000,
    # reference_mz = 200,
    # fwhm_rt = 9,
    # polarity = 'pos',
    # min_mz = 801,
    # max_mz = 803,
    # min_rt = 2808,
    # max_rt = 2928,
    # )
    raw_data = read_mzxml(
        '/data/toydata/toy_data_tof.mzXML',
        # FIXME: This is not correct, should be TOF, but not currently available.
        instrument_type='tof',
        resolution_ms1=30000,
        resolution_msn=30000,
        reference_mz=200,
        fwhm_rt=9,
        polarity='pos',
        min_mz=510,
        max_mz=531,
        min_rt=2390,
        max_rt=2510,
    )

    print("Resampling...")
    mesh = resample(raw_data, 10, 10, 0.5, 0.5)

    print("Saving mesh to disk...")
    mesh.save("mesh.dat")

    # print("Finding local maxima in mesh...")
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
        mesh_plot['img_plot'].scatter(
            local_max['i'], local_max['j'], color='aqua', s=5, marker="s", alpha=0.9)

    print("Fitting peaks...")
    peak_candidates = find_roi(raw_data, local_max)
    fitted_peaks = []
    if show_plot_fit:
        fig_mz = plt.figure()
        fig_rt = plt.figure()
    for peak_candidate in peak_candidates:
        try:
            # fitted_parameters = fitted_parameters + [fit(raw_data, peak_candidate)]
            # fitted_parameters = fit2(raw_data, peak_candidate)
            fitted_parameters = fit3(raw_data, peak_candidate)
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


def debugging_qatar():
    file_name = "/data/qatar/17122018/mzXML/Acute2U_3001.mzXML"
    tapp_parameters = {
        'instrument_type': 'orbitrap',
        'resolution_ms1': 75000,
        'resolution_msn': 30000,
        'reference_mz': 200,
        'avg_fwhm_rt': 9,
        'num_samples_mz': 5,
        'num_samples_rt': 5,
        'max_peaks': 1000000,
        # 'max_peaks': 20,
    }
    print("Reading raw data")
    raw_data = tapp.read_mzxml(
        file_name,
        instrument_type=tapp_parameters['instrument_type'],
        resolution_ms1=tapp_parameters['resolution_ms1'],
        resolution_msn=tapp_parameters['resolution_msn'],
        reference_mz=tapp_parameters['reference_mz'],
        # NOTE: For testing purposes
        fwhm_rt=tapp_parameters['avg_fwhm_rt'],
        # min_mz=313.06909,
        # max_mz=313.07223,
        min_mz=313.05,
        max_mz=313.08,
        min_rt=6.5 * 60,
        max_rt=14 * 60,
        # min_rt=417,
        # max_rt=440,
        polarity='pos',
    )
    print("Resampling")
    mesh = resample(raw_data, 10, 10, 0.5, 0.5)
    plt.style.use('dark_background')
    plt.ion()
    plt.show()
    # plot_mesh(mesh)

    # Testing internal peak finding routine.
    print("Finding peaks")
    peaks = find_peaks(raw_data, mesh)
    print("Fitting peaks via least_squares")
    fitted_peaks = []
    for peak_candidate in peaks:
        fitted_peak = fit_guos_2d_from_peak(peak_candidate)
        fitted_peaks = fitted_peaks + [fitted_peak]

    fitted_peaks = pd.DataFrame(fitted_peaks)
    fitted_peaks.columns = ['fitted_height', 'fitted_mz',
                            'fitted_sigma_mz', 'fitted_rt', 'fitted_sigma_rt']

    # peak = peaks[0]
    # print(fit_guos_2d_from_peak(peak))
    # data_points = find_raw_points(
    # raw_data,
    # peak.roi_min_mz,
    # peak.roi_max_mz,
    # peak.roi_min_rt,
    # peak.roi_max_rt
    # )
    # print(fit_guos_2d(
    # np.array(data_points.mz),
    # np.array(data_points.rt),
    # np.array(data_points.intensity)))
    peaks = pd.DataFrame(
        {
            'local_max_mz': np.array([peak.local_max_mz for peak in peaks]),
            'local_max_rt': np.array([peak.local_max_rt for peak in peaks]),
            'local_max_height': np.array([peak.local_max_height for peak in peaks]),
            'slope_descent_mz': np.array([peak.mesh_boundary_mz for peak in peaks]),
            'slope_descent_rt': np.array([peak.mesh_boundary_rt for peak in peaks]),
            'slope_descent_sigma_mz': np.array([peak.mesh_boundary_sigma_mz for peak in peaks]),
            'slope_descent_sigma_rt': np.array([peak.mesh_boundary_sigma_rt for peak in peaks]),
            'slope_descent_total_intensity': np.array([peak.mesh_boundary_total_intensity for peak in peaks]),
            'slope_descent_border_background': np.array([peak.mesh_boundary_border_background for peak in peaks]),
            'mesh_roi_mz': np.array([peak.mesh_roi_mz for peak in peaks]),
            'mesh_roi_rt': np.array([peak.mesh_roi_rt for peak in peaks]),
            'mesh_roi_sigma_mz': np.array([peak.mesh_roi_sigma_mz for peak in peaks]),
            'mesh_roi_sigma_rt': np.array([peak.mesh_roi_sigma_rt for peak in peaks]),
            'mesh_roi_total_intensity': np.array([peak.mesh_roi_total_intensity for peak in peaks]),
            'raw_roi_mz': np.array([peak.raw_roi_mz for peak in peaks]),
            'raw_roi_rt': np.array([peak.raw_roi_rt for peak in peaks]),
            'raw_roi_sigma_mz': np.array([peak.raw_roi_sigma_mz for peak in peaks]),
            'raw_roi_sigma_rt': np.array([peak.raw_roi_sigma_rt for peak in peaks]),
            'raw_roi_total_intensity': np.array([peak.raw_roi_total_intensity for peak in peaks]),
            'raw_roi_max_height': np.array([peak.raw_roi_max_height for peak in peaks]),
            'raw_roi_num_points': np.array([peak.raw_roi_num_points for peak in peaks]),
            'raw_roi_num_scans': np.array([peak.raw_roi_num_scans for peak in peaks]),
        })
    peaks = pd.concat([peaks, fitted_peaks], axis=1)

    # local_max = find_local_max(mesh)
    # local_max = pd.DataFrame(local_max)
    # local_max.columns = ['i', 'j', 'mz', 'rt', 'intensity']
    # local_max = local_max.sort_values('intensity', ascending=False)
    # # if max_peaks != math.inf:
    # # local_max = local_max[0:max_peaks]

    # if True:
    # print("Plotting mesh...")
    # mesh_plot = plot_mesh(mesh, transform='sqrt')

    # print("Plotting local maxima...")
    # mesh_plot['img_plot'].scatter(
    # local_max['i'], local_max['j'], color='aqua', s=5, marker="s", alpha=0.9)

    # print("Preparing peak candidates")
    # peak_candidates = find_roi(raw_data, local_max)
    # print("Fitting peaks...")
    # fitted_peaks = []
    # if True:
    # fig_mz = plt.figure()
    # fig_rt = plt.figure()
    # for peak_candidate in peak_candidates[0:10]:
    # # for peak_candidate in peak_candidates:
    # try:
    # # fitted_parameters = fitted_parameters + [fit(raw_data, peak_candidate)]
    # # fitted_parameters = fit2(raw_data, peak_candidate)
    # # fitted_parameters = fit3(raw_data, peak_candidate)
    # fitted_parameters = fit_raw_weighted_estimate(
    # raw_data, peak_candidate)
    # peak = peak_candidate
    # peak['fitted_height'] = fitted_parameters[0]
    # peak['fitted_mz'] = fitted_parameters[1]
    # peak['fitted_sigma_mz'] = fitted_parameters[2]
    # peak['fitted_rt'] = fitted_parameters[3]
    # peak['fitted_sigma_rt'] = fitted_parameters[4]
    # if (
    # peak['fitted_sigma_mz'] <= 0 or
    # peak['fitted_sigma_rt'] <= 0 or
    # peak['fitted_height'] <= 0 or
    # peak['fitted_mz'] > raw_data.max_mz or
    # peak['fitted_mz'] < raw_data.min_mz or
    # peak['fitted_rt'] > raw_data.max_rt or
    # peak['fitted_rt'] < raw_data.min_rt or
    # np.isnan(peak['fitted_sigma_mz']).any() or
    # np.isnan(peak['fitted_sigma_rt']).any() or
    # np.isnan(peak['fitted_height']).any() or
    # np.isnan(peak['fitted_mz']).any() or
    # np.isnan(peak['fitted_rt']).any() or
    # np.isinf(peak['fitted_sigma_mz']).any() or
    # np.isinf(peak['fitted_sigma_rt']).any() or
    # np.isinf(peak['fitted_height']).any() or
    # np.isinf(peak['fitted_mz']).any() or
    # np.isinf(peak['fitted_rt']).any()
    # ):
    # continue
    # fitted_peaks = fitted_peaks + [peak]
    # if True:
    # plot_peak_fit(raw_data, peak, fig_mz, fig_rt)
    # except Exception as e:
    # print(e)
    # pass

    # return raw_data, mesh, local_max, fitted_peaks
    return raw_data, mesh, peaks


RawData.tic = tic
