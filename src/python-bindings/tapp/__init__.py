# TODO(alex): Write documentation.

from scipy.special import erfc
import math
import os
# TODO: Use pathlib instead of os?
# from pathlib import Path

import tapp
from .tapp import *
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from scipy.optimize import curve_fit


def plot_mesh(mesh, transform='sqrt', figure=None):
    plt.style.use('dark_background')
    plt.ion()
    plt.show()

    if figure is None:
        figure = plt.figure()

    img = mesh.data
    img = np.reshape(img, (mesh.m, mesh.n))
    bins_rt = mesh.bins_rt
    bins_mz = mesh.bins_mz

    plt.figure(figure.number)
    plt.clf()
    gs = gridspec.GridSpec(5, 5)
    mz_plot = plt.subplot(gs[0, :-1])
    mz_plot.clear()
    mz_plot.plot(bins_mz, img.sum(axis=0))
    mz_plot.margins(x=0)
    mz_plot.set_xticks([])
    mz_plot.set_ylabel("Intensity")

    rt_plot = plt.subplot(gs[1:, -1])
    rt_plot.plot(img.sum(axis=1), bins_rt)
    rt_plot.margins(y=0)
    rt_plot.set_yticks([])
    rt_plot.set_xlabel("Intensity")

    img_plot = plt.subplot(gs[1:, :-1])
    offset_rt = (np.array(mesh.bins_rt).max() -
                 np.array(mesh.bins_rt).min())/mesh.m / 2
    offset_mz = (np.array(mesh.bins_mz).max() -
                 np.array(mesh.bins_mz).min())/mesh.n / 2
    if transform == 'sqrt':
        img_plot.pcolormesh(
            np.array(mesh.bins_mz) - offset_mz,
            np.array(mesh.bins_rt) - offset_rt,
            img,
            snap=True,
            norm=colors.PowerNorm(gamma=1./2.))
    elif transform == 'cubic':
        img_plot.pcolormesh(
            np.array(mesh.bins_mz) - offset_mz,
            np.array(mesh.bins_rt) - offset_rt,
            img,
            norm=colors.PowerNorm(gamma=1./3.))
    elif transform == 'log':
        img_plot.pcolormesh(
            np.array(mesh.bins_mz) - offset_mz,
            np.array(mesh.bins_rt) - offset_rt,
            img,
            norm=colors.LogNorm(vmin=img.min()+1e-8, vmax=img.max()))
    else:
        img_plot.pcolormesh(mesh.bins_mz, mesh.bins_rt, img)

    img_plot.set_xlim([np.array(mesh.bins_mz).min() - offset_mz,
                       np.array(mesh.bins_mz).max() - offset_mz])
    img_plot.set_ylim([np.array(mesh.bins_rt).min() - offset_rt,
                       np.array(mesh.bins_rt).max() - offset_rt])

    img_plot.set_xlabel("m/z")
    img_plot.set_ylabel("retention time (s)")

    return({
        "img_plot": img_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })


def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def gaus2d(X, a, x_0, sigma_x, y_0, sigma_y):
    x = X[0]
    y = X[1]
    b = np.exp(-0.5 * ((x - x_0) / sigma_x) ** 2)
    c = np.exp(-0.5 * ((y - y_0)/sigma_y) ** 2)
    return a * b * c


def fit_curvefit(mzs, intensities, rts):
    X = np.array([mzs, rts])
    mean_x = sum(X[0] * intensities) / sum(intensities)
    sigma_x = np.sqrt(sum(intensities * (X[0] - mean_x)**2) / sum(intensities))
    mean_y = sum(X[1] * intensities) / sum(intensities)
    sigma_y = np.sqrt(sum(intensities * (X[1] - mean_y)**2) / sum(intensities))
    fitted_parameters, pcov_2d = curve_fit(gaus2d, X, intensities, p0=[
        max(intensities), mean_x, sigma_x, mean_y, sigma_y])
    return fitted_parameters


def fit_inle(x, y):
    mean_x = np.sum(x * y) / np.sum(y)
    sigma_x = np.sqrt(sum(y * (x - mean_x)**2) / sum(y))
    fitted, cov = curve_fit(gaus, x, y, p0=[max(y), mean_x, sigma_x])
    return fitted


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
    beta = np.linalg.lstsq(X, Y)
    a, b, c, d, e = beta[0]

    sigma_mz = np.sqrt(1/(-2 * c))
    mz = b / (-2 * c) + x_mean
    sigma_rt = np.sqrt(1/(-2 * e))
    rt = d / (-2 * e) + y_mean
    height = np.exp(a - ((b ** 2) / (4 * c)) - ((d ** 2) / (4 * e)))

    return np.array([height, mz, sigma_mz, rt, sigma_rt])


def fit_weighted_guos_2d(x, y, z, x_center, y_center, sig_mz, sig_rt):
    x = x - x_center
    y = y - y_center

    a_0_0 = 0
    a_0_1 = 0
    a_0_2 = 0
    a_0_3 = 0
    a_0_4 = 0
    a_1_0 = 0
    a_1_1 = 0
    a_1_2 = 0
    a_1_3 = 0
    a_1_4 = 0
    a_2_0 = 0
    a_2_1 = 0
    a_2_2 = 0
    a_2_3 = 0
    a_2_4 = 0
    a_3_0 = 0
    a_3_1 = 0
    a_3_2 = 0
    a_3_3 = 0
    a_3_4 = 0
    a_4_0 = 0
    a_4_1 = 0
    a_4_2 = 0
    a_4_3 = 0
    a_4_4 = 0
    c_0 = 0
    c_1 = 0
    c_2 = 0
    c_3 = 0
    c_4 = 0
    for i, intensity in enumerate(z):
        mz = x[i]
        rt = y[i]
        w = gaus2d([mz, rt], 1, 0, sig_mz, 0, sig_rt)
        w_2 = w * w

        a_0_0 += w_2
        a_0_1 += w_2 * mz
        a_0_2 += w_2 * mz * mz
        a_0_3 += w_2 * rt
        a_0_4 += w_2 * rt * rt

        a_1_0 += w_2 * mz
        a_1_1 += w_2 * mz * mz
        a_1_2 += w_2 * mz * mz * mz
        a_1_3 += w_2 * rt * mz
        a_1_4 += w_2 * rt * rt * mz

        a_2_0 += w_2 * mz * mz
        a_2_1 += w_2 * mz * mz * mz
        a_2_2 += w_2 * mz * mz * mz * mz
        a_2_3 += w_2 * rt * mz * mz
        a_2_4 += w_2 * rt * rt * mz * mz

        a_3_0 += w_2 * rt
        a_3_1 += w_2 * mz * rt
        a_3_2 += w_2 * mz * mz * rt
        a_3_3 += w_2 * rt * rt
        a_3_4 += w_2 * rt * rt * rt

        a_4_0 += w_2 * rt * rt
        a_4_1 += w_2 * mz * rt * rt
        a_4_2 += w_2 * mz * mz * rt * rt
        a_4_3 += w_2 * rt * rt * rt
        a_4_4 += w_2 * rt * rt * rt * rt

        c_0 += w_2 * np.log(intensity)
        c_1 += w_2 * np.log(intensity) * mz
        c_2 += w_2 * np.log(intensity) * mz * mz
        c_3 += w_2 * np.log(intensity) * rt
        c_4 += w_2 * np.log(intensity) * rt * rt

    X = np.array(
        [
            [
                a_0_0,
                a_0_1,
                a_0_2,
                a_0_3,
                a_0_4,
            ],
            [
                a_1_0,
                a_1_1,
                a_1_2,
                a_1_3,
                a_1_4,
            ],
            [
                a_2_0,
                a_2_1,
                a_2_2,
                a_2_3,
                a_2_4,
            ],
            [
                a_3_0,
                a_3_1,
                a_3_2,
                a_3_3,
                a_3_4,
            ],
            [
                a_4_0,
                a_4_1,
                a_4_2,
                a_4_3,
                a_4_4,
            ],
        ],
    )
    Y = np.array([
        c_0,
        c_1,
        c_2,
        c_3,
        c_4,
    ])
    a, b, c, d, e = np.linalg.lstsq(X, Y, rcond=1)[0]

    if c >= 0 or e >= 0:
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])

    sigma_mz = np.sqrt(1/(-2 * c))
    mz = b / (-2 * c) + x_center
    sigma_rt = np.sqrt(1/(-2 * e))
    rt = d / (-2 * e) + y_center
    height = np.exp(a - ((b ** 2) / (4 * c)) - ((d ** 2) / (4 * e)))

    return np.array([height, mz, sigma_mz, rt, sigma_rt])


def fit_weighted_guos_2d_const(x, y, z, x_center, y_center, sig_mz, sig_rt):
    x = x - x_center
    y = y - y_center

    a_0_0 = 0
    a_0_2 = 0
    a_0_4 = 0
    a_2_0 = 0
    a_2_2 = 0
    a_2_4 = 0
    a_4_0 = 0
    a_4_2 = 0
    a_4_4 = 0
    c_0 = 0
    c_2 = 0
    c_4 = 0
    for i, intensity in enumerate(z):
        mz = x[i]
        rt = y[i]
        w = gaus2d([mz, rt], 1, 0, sig_mz, 0, sig_rt)
        w_2 = w * w

        a_0_0 += w_2
        a_0_2 += w_2 * mz * mz
        a_0_4 += w_2 * rt * rt

        a_2_0 += w_2 * mz * mz
        a_2_2 += w_2 * mz * mz * mz * mz
        a_2_4 += w_2 * rt * rt * mz * mz

        a_4_0 += w_2 * rt * rt
        a_4_2 += w_2 * mz * mz * rt * rt
        a_4_4 += w_2 * rt * rt * rt * rt

        c_0 += w_2 * np.log(intensity)
        c_2 += w_2 * np.log(intensity) * mz * mz
        c_4 += w_2 * np.log(intensity) * rt * rt

    X = np.array(
        [
            [
                a_0_0,
                a_0_2,
                a_0_4,
            ],
            [
                a_2_0,
                a_2_2,
                a_2_4,
            ],
            [
                a_4_0,
                a_4_2,
                a_4_4,
            ],
        ],
    )
    Y = np.array([
        c_0,
        c_2,
        c_4,
    ])
    a, c, e = np.linalg.lstsq(X, Y, rcond=1)[0]

    if c >= 0 or e >= 0:
        return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])

    sigma_mz = np.sqrt(1/(-2 * c))
    mz = x_center
    sigma_rt = np.sqrt(1/(-2 * e))
    rt = y_center
    height = np.exp(a)

    return np.array([height, mz, sigma_mz, rt, sigma_rt])


def fit2(raw_data, peak_candidate):
    data_points = tapp.find_raw_points(
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
    data_points = tapp.find_raw_points(
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
    data_points = tapp.find_raw_points(
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

    data_points = tapp.find_raw_points(
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
    markerline, stemlines, baseline = plt.stem(
        np.array(mzs)[sort_idx_mz],
        np.array(intensities)[sort_idx_mz],
        label='intensities',
        markerfmt=' '
    )
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
    plt.xlabel('retention time (s)')
    plt.ylabel('Intensity')

    return fig_mz, fig_rt


def fit_and_plot(raw_data, peak_candidate):
    data_points = tapp.find_raw_points(
        raw_data,
        peak_candidate['roi_min_mz'],
        peak_candidate['roi_max_mz'],
        peak_candidate['roi_min_rt'],
        peak_candidate['roi_max_rt']
    )
    fitted_parameters = fit_curvefit(
        data_points.mz, data_points.intensity, data_points.rt)
    tapp.plot_peak_candidate(data_points, fitted_parameters)
    return fitted_parameters


def find_roi(raw_data, local_max, avg_rt_fwhm=10):
    peak_candidates = []
    for i in range(0, len(local_max)):
        selected_peak = local_max.iloc[i]
        theoretical_sigma_mz = tapp.theoretical_fwhm(
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


def peak_extraction(file_name, tapp_parameters):
    print("Reading raw data")
    raw_data = tapp.read_mzxml(
        file_name,
        min_mz=tapp_parameters['min_mz'],
        max_mz=tapp_parameters['max_mz'],
        min_rt=tapp_parameters['min_rt'],
        max_rt=tapp_parameters['max_rt'],
        instrument_type=tapp_parameters['instrument_type'],
        resolution_ms1=tapp_parameters['resolution_ms1'],
        resolution_msn=tapp_parameters['resolution_msn'],
        reference_mz=tapp_parameters['reference_mz'],
        fwhm_rt=tapp_parameters['avg_fwhm_rt'],
        polarity=tapp_parameters['polarity'],
    )
    print("Resampling")
    mesh = tapp.resample(
        raw_data,
        tapp_parameters['num_samples_mz'],
        tapp_parameters['num_samples_rt'],
        tapp_parameters['smoothing_coefficient_mz'],
        tapp_parameters['smoothing_coefficient_rt'],
    )

    print("Finding peaks")
    peaks = tapp.find_peaks(raw_data, mesh, tapp_parameters['max_peaks'])
    print("Found {} peaks".format(len(peaks)))

    return raw_data, mesh, peaks


def plot_xic(peak, raw_data, figure=None, method="max"):
    xic = peak.xic(raw_data, method=method)
    x = xic.retention_time
    y = xic.intensity
    plt.style.use('dark_background')
    if not figure:
        figure = plt.figure()

    plt.ion()
    plt.show()
    plt.plot(x, y, label='peak_id = {}'.format(peak.id))
    plt.xlabel('Retention time (s)')
    plt.ylabel('Intensity')
    plt.legend()

    return figure


def fit_sigmas(peak):
    X = np.array(
        [
            [
                peak.a_2_2(),
                peak.a_2_4(),
            ],
            [
                peak.a_4_2(),
                peak.a_4_4(),
            ],
        ],
    )
    Y = np.array([
        peak.c_2(),
        peak.c_4(),
    ])

    print(X, Y)
    beta = np.linalg.lstsq(X, Y, rcond=None)
    beta_2, beta_4 = beta[0]
    var_x = 1/(2 * beta_2)
    var_y = 1/(2 * beta_4)

    print(var_x, var_y)
    if var_x <= 0 or var_y <= 0:
        return np.nan, np.nan

    sigma_mz = np.sqrt(var_x)
    sigma_rt = np.sqrt(var_y)

    return sigma_mz, sigma_rt


def fit_height_and_sigmas(peak):
    X = np.array(
        [
            [
                peak.a_0_0(),
                peak.a_0_2(),
                peak.a_0_4(),
            ],
            [
                peak.a_2_0(),
                peak.a_2_2(),
                peak.a_2_4(),
            ],
            [
                peak.a_4_0(),
                peak.a_4_2(),
                peak.a_4_4(),
            ],
        ],
    )
    Y = np.array([
        peak.c_0(),
        peak.c_2(),
        peak.c_4(),
    ])

    print(X, Y)
    beta = np.linalg.lstsq(X, Y, rcond=None)
    beta_0, beta_2, beta_4 = beta[0]
    var_x = -1/(2 * beta_2)
    var_y = -1/(2 * beta_4)

    print(var_x, var_y)
    # if var_x <= 0 or var_y <= 0:
    # return np.nan, np.nan

    sigma_mz = np.sqrt(var_x)
    sigma_rt = np.sqrt(var_y)
    height = np.exp(beta_0)

    return height, sigma_mz, sigma_rt


def plot_raw_points(
        peak,
        raw_data,
        img_plot=None,
        rt_plot=None,
        mz_plot=None,
        xic_method="max"):
    data_points = tapp.find_raw_points(
        raw_data,
        peak.roi_min_mz,
        peak.roi_max_mz,
        peak.roi_min_rt,
        peak.roi_max_rt,
    )
    rts = data_points.rt
    mzs = data_points.mz
    intensities = data_points.intensity

    # Calculate min/max values for the given peak.
    min_mz = peak.roi_min_mz
    max_mz = peak.roi_max_mz
    min_rt = peak.roi_min_rt
    max_rt = peak.roi_max_rt

    if not img_plot and not rt_plot and not mz_plot:
        plt.style.use('dark_background')
        plt.ion()
        plt.show()
        plt.figure()
        plt.clf()
        gs = gridspec.GridSpec(5, 5)
        mz_plot = plt.subplot(gs[0, :-1])
        mz_plot.margins(x=0)
        mz_plot.set_xticks([])
        mz_plot.set_ylabel("Intensity")
        rt_plot = plt.subplot(gs[1:, -1])
        rt_plot.margins(y=0)
        rt_plot.set_yticks([])
        rt_plot.set_xlabel("Intensity")
        img_plot = plt.subplot(gs[1:, :-1])

        # Set the min/max limits for mz/rt.
        mz_plot.set_xlim([min_mz, max_mz])
        rt_plot.set_ylim([min_rt, max_rt])
        img_plot.set_xlim([min_mz, max_mz])
        img_plot.set_ylim([min_rt, max_rt])

    # NOTE: Adding 200 for a more pleasant color map on the first peaks, found
    # this number by trial and error, dont @ me.
    np.random.seed(peak.id + 200)
    color = np.append(np.random.rand(3, 1).flatten(), 0.5)
    np.random.seed(None)

    if img_plot:
        img_plot.scatter(
            mzs, rts,
            c=np.sqrt(intensities),
            edgecolor=color,
        )
    if rt_plot:
        xic = peak.xic(raw_data, method=xic_method)
        x = xic.retention_time
        y = xic.intensity
        rt_plot.plot(y, x, color=color)
    if mz_plot:
        sort_idx_mz = np.argsort(mzs)
        markerline, stemlines, baseline = mz_plot.stem(
            np.array(mzs)[sort_idx_mz],
            np.array(intensities)[sort_idx_mz],
            markerfmt=' ',
        )
        plt.setp(baseline, color=color, alpha=0.5)
        plt.setp(stemlines, color=color, alpha=0.5)

    # Set x/y limits if necessary.
    lim_min_mz, lim_max_mz = img_plot.get_xlim()
    lim_min_rt, lim_max_rt = img_plot.get_ylim()
    if min_mz < lim_min_mz:
        lim_min_mz = min_mz
    if min_rt < lim_min_rt:
        lim_min_rt = min_rt
    if max_mz > lim_max_mz:
        lim_max_mz = max_mz
    if max_rt > lim_max_rt:
        lim_max_rt = max_rt
    mz_plot.set_xlim([lim_min_mz, lim_max_mz])
    rt_plot.set_ylim([lim_min_rt, lim_max_rt])
    img_plot.set_xlim([lim_min_mz, lim_max_mz])
    img_plot.set_ylim([lim_min_rt, lim_max_rt])

    return({
        "img_plot": img_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })


def plot_sigma(
        peak,
        height, mz, rt,
        sigma_mz, sigma_rt,
        img_plot=None, rt_plot=None, mz_plot=None,
        linestyle='--',
        label=None,
        marker='.',
):
    # Calculate min/max values for the given peak.
    min_mz = mz - 3 * sigma_mz
    max_mz = mz + 3 * sigma_mz
    min_rt = rt - 3 * sigma_rt
    max_rt = rt + 3 * sigma_rt

    if not img_plot and not rt_plot and not mz_plot:
        plt.style.use('dark_background')
        plt.ion()
        plt.show()
        plt.figure()
        plt.clf()
        gs = gridspec.GridSpec(5, 5)
        mz_plot = plt.subplot(gs[0, :-1])
        mz_plot.margins(x=0)
        mz_plot.set_xticks([])
        mz_plot.set_ylabel("Intensity")
        rt_plot = plt.subplot(gs[1:, -1])
        rt_plot.margins(y=0)
        rt_plot.set_yticks([])
        rt_plot.set_xlabel("Intensity")
        img_plot = plt.subplot(gs[1:, :-1])

        # Set the min/max limits for mz/rt.
        mz_plot.set_xlim([min_mz, max_mz])
        rt_plot.set_ylim([min_rt, max_rt])
        img_plot.set_xlim([min_mz, max_mz])
        img_plot.set_ylim([min_rt, max_rt])

    # NOTE: Adding 200 for a more pleasant color map on the first peaks, found
    # this number by trial and error, dont @ me.
    np.random.seed(peak.id + 200)
    base_color = np.random.rand(3, 1).flatten()
    np.random.seed(None)

    lim_min_mz, lim_max_mz = img_plot.get_xlim()
    lim_min_rt, lim_max_rt = img_plot.get_ylim()
    if img_plot:
        # Set the limits for the img_plot
        if min_mz < lim_min_mz:
            lim_min_mz = min_mz
        if min_rt < lim_min_rt:
            lim_min_rt = min_rt
        if max_mz > lim_max_mz:
            lim_max_mz = max_mz
        if max_rt > lim_max_rt:
            lim_max_rt = max_rt
        img_plot.set_xlim([lim_min_mz, lim_max_mz])
        img_plot.set_ylim([lim_min_rt, lim_max_rt])

        # Plotting the center of the peak.
        color_0 = np.append(base_color, 1)
        img_plot.scatter(
            mz, rt,
            marker=marker,
            label=label,
            color=color_0, facecolors='none', edgecolors=color_0,
        )

        color_1 = np.append(base_color, 0.9)
        elip_1 = Ellipse(
            (mz, rt),
            2 * sigma_mz,
            2 * sigma_rt,
            fill=False,
            color=color_1,
            linestyle=linestyle,
        )
        color_2 = np.append(base_color, 0.6)
        elip_2 = Ellipse(
            (mz, rt),
            2 * 2 * sigma_mz,
            2 * 2 * sigma_rt,
            fill=False,
            color=color_2,
            linestyle=linestyle,
        )
        color_3 = np.append(base_color, 0.4)
        elip_3 = Ellipse(
            (mz, rt),
            3 * 2 * sigma_mz,
            3 * 2 * sigma_rt,
            fill=False,
            color=color_3,
            linestyle=linestyle,
        )
        img_plot.add_artist(elip_1)
        img_plot.add_artist(elip_2)
        img_plot.add_artist(elip_3)
    if rt_plot:
        # Set the limits for mz_plot.
        if min_rt < lim_min_rt:
            lim_min_rt = min_rt
        if max_rt > lim_max_rt:
            lim_max_rt = max_rt
        img_plot.set_xlim([lim_min_mz, lim_max_mz])
        img_plot.set_ylim([lim_min_rt, lim_max_rt])
        rt_plot.set_ylim([lim_min_rt, lim_max_rt])
        x = np.linspace(min_rt, max_rt, 100)
        y = gaus(x, height, rt, sigma_rt)
        rt_plot.plot(
            y, x,
            linestyle=linestyle,
            color=base_color,
            label=label,
        )
    if mz_plot:
        # Set the limits for rt_plot.
        if min_mz < lim_min_mz:
            lim_min_mz = min_mz
        if max_mz > lim_max_mz:
            lim_max_mz = max_mz
        mz_plot.set_xlim([lim_min_mz, lim_max_mz])
        x = np.linspace(min_mz, max_mz, 100)
        y = gaus(x, height, mz, sigma_mz)
        mz_plot.plot(
            x, y,
            linestyle=linestyle,
            color=base_color,
            label=label,
        )

    return({
        "img_plot": img_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })


def plot_raw_roi_sigma(peak, img_plot=None, rt_plot=None, mz_plot=None):
    return plot_sigma(
        peak,
        peak.local_max_height,
        peak.local_max_mz,
        peak.local_max_rt,
        peak.raw_roi_sigma_mz,
        peak.raw_roi_sigma_rt,
        img_plot,
        rt_plot,
        mz_plot,
        label='raw_roi',
        marker='s',
    )


def plot_raw_roi_fitted_sigma(
        peak, raw_data, img_plot=None, rt_plot=None, mz_plot=None):
    data_points = tapp.find_raw_points(
        raw_data,
        peak.roi_min_mz,
        peak.roi_max_mz,
        peak.roi_min_rt,
        peak.roi_max_rt
    )
    mzs = np.array(data_points.mz)
    rts = np.array(data_points.rt)
    intensities = np.array(data_points.intensity)
    X = np.array([mzs, rts])
    def f(x, h, mz, sigma_mz, rt, sigma_rt): return gaus2d(
        x, h, peak.local_max_mz, sigma_mz, peak.local_max_rt, sigma_rt)
    fitted_parameters, pcov_2d = curve_fit(
        f, X, intensities, p0=[
            peak.raw_roi_max_height,
            peak.local_max_mz,
            peak.raw_roi_sigma_mz,
            peak.local_max_rt,
            peak.raw_roi_sigma_rt
        ])

    h, mz, sigma_mz, rt, sigma_rt = fitted_parameters

    return plot_sigma(
        peak,
        h,
        mz,
        rt,
        sigma_mz,
        sigma_rt,
        img_plot,
        rt_plot,
        mz_plot,
        linestyle='-',
        label='fitted_raw_roi',
        marker='.',
    )


def plot_raw_roi_fitted_sigma_fast(
        peak, raw_data, img_plot=None, rt_plot=None, mz_plot=None):
    data_points = tapp.find_raw_points(
        raw_data,
        peak.roi_min_mz,
        peak.roi_max_mz,
        peak.roi_min_rt,
        peak.roi_max_rt
    )
    mzs = np.array(data_points.mz)
    rts = np.array(data_points.rt)
    intensities = np.array(data_points.intensity)
    fitted_parameters = fit_guos_2d(mzs, rts, intensities)

    height, mz, sigma_mz, rt, sigma_rt = fitted_parameters

    return plot_sigma(
        peak,
        height,
        mz,
        rt,
        sigma_mz,
        sigma_rt,
        img_plot,
        rt_plot,
        mz_plot,
        linestyle=':',
        label='fitted_raw_roi (fast)',
        marker='P',
    )


def plot_raw_roi_fitted_sigma_weighted(
        peak, raw_data, img_plot=None, rt_plot=None, mz_plot=None):
    data_points = tapp.find_raw_points(
        raw_data,
        peak.roi_min_mz,
        peak.roi_max_mz,
        peak.roi_min_rt,
        peak.roi_max_rt
    )
    mzs = np.array(data_points.mz)
    rts = np.array(data_points.rt)
    intensities = np.array(data_points.intensity)
    fwhm_mz = raw_data.theoretical_fwhm(peak.local_max_mz)
    theoretical_sigma_rt = raw_data.fwhm_rt/(2 * np.sqrt(2 * np.log(2)))
    theoretical_sigma_mz = fwhm_mz/(2 * np.sqrt(2 * np.log(2)))

    # IMPORTANT: Since multiple peaks might appear within the 3 * sigma ROI of
    # a peak, the R2 calculated from this number can be skewed. For this
    # reason, we are only using +-1 * sigma for the estimation of R2.
    min_mz = peak.local_max_mz - theoretical_sigma_mz
    max_mz = peak.local_max_mz + theoretical_sigma_mz
    min_rt = peak.local_max_rt - theoretical_sigma_rt
    max_rt = peak.local_max_rt + theoretical_sigma_rt
    idx = (mzs > min_mz) & (mzs < max_mz) & (rts > min_rt) & (rts < max_rt)
    mzs = np.copy(mzs[idx])
    rts = np.copy(rts[idx])
    intensities = np.copy(intensities[idx])
    fitted_parameters = fit_weighted_guos_2d(
        mzs, rts, intensities,
        peak.local_max_mz, peak.local_max_rt,
        theoretical_sigma_mz, theoretical_sigma_rt)

    height, mz, sigma_mz, rt, sigma_rt = fitted_parameters

    return plot_sigma(
        peak,
        height,
        mz,
        rt,
        sigma_mz,
        sigma_rt,
        img_plot,
        rt_plot,
        mz_plot,
        linestyle='-',
        label='fitted_raw_roi (weighted)',
        marker='P',
    )


def plot_theoretical_sigma(
        peak, raw_data, img_plot=None, rt_plot=None, mz_plot=None):
    fwhm_mz = raw_data.theoretical_fwhm(peak.local_max_mz)
    theoretical_sigma_rt = raw_data.fwhm_rt/(2 * np.sqrt(2 * np.log(2)))
    theoretical_sigma_mz = fwhm_mz/(2 * np.sqrt(2 * np.log(2)))
    return plot_sigma(
        peak,
        peak.local_max_height,
        peak.local_max_mz,
        peak.local_max_rt,
        theoretical_sigma_mz,
        theoretical_sigma_rt,
        img_plot,
        rt_plot,
        mz_plot,
        linestyle=':',
        label='fitted_raw_roi (fast)',
        marker='P',
    )


def testing_different_sigmas(peaks, raw_data, i=0):
    plots = peaks[i].plot_raw_points(raw_data)
    plots = peaks[i].plot_raw_roi_sigma(
        plots['img_plot'], plots['rt_plot'], plots['mz_plot'])
    plots = peaks[i].plot_raw_roi_fitted_sigma_weighted(
        raw_data, plots['img_plot'], plots['rt_plot'], plots['mz_plot'])
    plt.legend()
    return plots


def emg(t, h, tg, sigma, tau):
    a = 1/(2 * tau)
    b = 1/2 * np.power(sigma/tau, 2) - (t - tg)/tau
    z = 1 / np.sqrt(2) * ((t - tg)/sigma - sigma/tau)
    c = erfc(-z)
    return h * a * np.exp(b) * c


def gauss_mz_emg_rt(X, h, mz_0, sigma_mz, rt_0, sigma_rt, tau):
    mz = X[0]
    rt = X[1]
    a = 1/(2 * tau)
    b = 1/2 * np.power(sigma_rt/tau, 2) - (rt - rt_0)/tau
    z = 1 / np.sqrt(2) * ((rt - rt_0)/sigma_rt - sigma_rt/tau)
    c = erfc(-z)
    d = np.exp(-0.5 * np.power((mz - mz_0) / sigma_mz, 2))
    return h * a * np.exp(b) * c * d


def calculate_r2(x, y, z, h, mz, sigma_mz, rt, sigma_rt):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    ss_tot = (np.power(z - z.mean(), 2)).sum()
    ss_res = (np.power(z - gaus2d([x, y], h, mz,
                                  sigma_mz, rt, sigma_rt), 2)).sum()
    r2 = 1 - ss_res / ss_tot
    return r2


def fit_and_evaluate_r2(peak, raw_data):
    data_points = tapp.find_raw_points(
        raw_data,
        peak.roi_min_mz,
        peak.roi_max_mz,
        peak.roi_min_rt,
        peak.roi_max_rt
    )
    mzs = np.array(data_points.mz)
    rts = np.array(data_points.rt)
    intensities = np.array(data_points.intensity)
    fwhm_mz = raw_data.theoretical_fwhm(peak.local_max_mz)
    theoretical_sigma_rt = raw_data.fwhm_rt/(2 * np.sqrt(2 * np.log(2)))
    theoretical_sigma_mz = fwhm_mz/(2 * np.sqrt(2 * np.log(2)))

    # IMPORTANT: Since multiple peaks might appear within the 3 * sigma ROI of
    # a peak, the R2 calculated from this number can be skewed. For this
    # reason, we are only using +-1 * sigma for the estimation of R2.
    min_mz = peak.local_max_mz - theoretical_sigma_mz
    max_mz = peak.local_max_mz + theoretical_sigma_mz
    min_rt = peak.local_max_rt - theoretical_sigma_rt
    max_rt = peak.local_max_rt + theoretical_sigma_rt
    idx = (mzs > min_mz) & (mzs < max_mz) & (rts > min_rt) & (rts < max_rt)
    mzs = np.copy(mzs[idx])
    rts = np.copy(rts[idx])
    intensities = np.copy(intensities[idx])
    if idx.sum() == 0:
        print(peak.id)

    # Fit 0: Theoretical.
    theoretical_r2 = calculate_r2(
        mzs, rts, intensities,
        peak.local_max_height,
        peak.local_max_mz,
        theoretical_sigma_mz,
        peak.local_max_rt,
        theoretical_sigma_rt,
    )
    # Fit 1: Estimated.
    estimated_r2 = calculate_r2(
        mzs, rts, intensities,
        peak.local_max_height,
        peak.local_max_mz,
        peak.raw_roi_sigma_mz,
        peak.local_max_rt,
        peak.raw_roi_sigma_rt,
    )
    # Fit 2: Weighted Least Square Fitting.
    fitted_parameters = fit_weighted_guos_2d_const(
        mzs, rts, intensities,
        peak.local_max_mz, peak.local_max_rt,
        theoretical_sigma_mz, theoretical_sigma_rt)
    height, mz, sigma_mz, rt, sigma_rt = fitted_parameters
    fitted_r2 = calculate_r2(
        mzs, rts, intensities,
        height,
        mz,
        sigma_mz,
        rt,
        sigma_rt,
    )
    return (peak.id, theoretical_r2, estimated_r2, fitted_r2)


def calculate_r2_all_peaks(peaks, raw_data, plot_density):
    r2_values = [fit_and_evaluate_r2(
        peak, raw_data, verbose=False) for peak in peaks]
    r2_values = pd.DataFrame(
        r2_values,
        columns=['id', 'theoretical_r2', 'estimated_r2', 'fitted_r2'])
    if plot_density:
        import seaborn as sns
        plt.style.use('dark_background')
        plt.ion()
        plt.show()
        plt.figure()
        sns.distplot(r2_values['theoretical_r2'].dropna(),
                     hist=False, label='theoretical_r2')
        sns.distplot(r2_values['estimated_r2'].dropna(),
                     hist=False, label='estimated_r2')
        sns.distplot(r2_values['fitted_r2'].dropna(),
                     hist=False, label='fitted_r2')
    return r2_values


def to_table(peaks):
    peaks_df = pd.DataFrame(
        {
            'id': np.array([peak.id for peak in peaks]),
            'local_max_mz': np.array([peak.local_max_mz for peak in peaks]),
            'local_max_rt': np.array([peak.local_max_rt for peak in peaks]),
            'local_max_height': np.array([peak.local_max_height for peak in peaks]),
            'raw_roi_mean_mz': np.array([peak.raw_roi_mean_mz for peak in peaks]),
            'raw_roi_mean_rt': np.array([peak.raw_roi_mean_rt for peak in peaks]),
            'raw_roi_sigma_mz': np.array([peak.raw_roi_sigma_mz for peak in peaks]),
            'raw_roi_sigma_rt': np.array([peak.raw_roi_sigma_rt for peak in peaks]),
            'raw_roi_skewness_mz': np.array([peak.raw_roi_skewness_mz for peak in peaks]),
            'raw_roi_skewness_rt': np.array([peak.raw_roi_skewness_rt for peak in peaks]),
            'raw_roi_kurtosis_mz': np.array([peak.raw_roi_kurtosis_mz for peak in peaks]),
            'raw_roi_kurtosis_rt': np.array([peak.raw_roi_kurtosis_rt for peak in peaks]),
            'raw_roi_max_height': np.array([peak.raw_roi_max_height for peak in peaks]),
            'raw_roi_total_intensity': np.array([peak.raw_roi_total_intensity for peak in peaks]),
            'raw_roi_num_points': np.array([peak.raw_roi_num_points for peak in peaks]),
            'raw_roi_num_scans': np.array([peak.raw_roi_num_scans for peak in peaks]),
        })
    return peaks_df


def linked_peptides_to_table(linked_peptides):
    linked_peptides_df = pd.DataFrame(
        {
            'sequence': np.array([linked_peptide.sequence for linked_peptide in linked_peptides]),
            'charge_state': np.array([linked_peptide.charge_state for linked_peptide in linked_peptides]),
            'ident_rt': np.array([linked_peptide.ident_rt for linked_peptide in linked_peptides]),
            'ident_mz': np.array([linked_peptide.ident_mz for linked_peptide in linked_peptides]),
            'number_of_isotopes': np.array([len(linked_peptide.linked_isotopes) for linked_peptide in linked_peptides]),
            'monoisotopic_height': np.array([linked_peptide.monoisotopic_height for linked_peptide in linked_peptides]),
            'monoisotopic_intensity': np.array([linked_peptide.monoisotopic_intensity for linked_peptide in linked_peptides]),
            'total_height': np.array([linked_peptide.total_height for linked_peptide in linked_peptides]),
            'total_intensity': np.array([linked_peptide.total_intensity for linked_peptide in linked_peptides]),
            'weighted_error': np.array([linked_peptide.weighted_error for linked_peptide in linked_peptides]),
            'psm_id': np.array([linked_peptide.psm_id for linked_peptide in linked_peptides]),
        })
    return linked_peptides_df


def create_psm_protein_graph(ident_data):
    unique_proteins = pd.Series(
        [
            protein_hypothesis.db_sequence_id
            for protein_hypothesis in ident_data.protein_hypotheses
        ]).unique()
    unique_psm = np.unique(np.concatenate(
        [
            protein_hypothesis.spectrum_ids
            for protein_hypothesis in ident_data.protein_hypotheses
        ]))
    # NOTE: Probably we want to use adjancency list instead of matrix, as the matrix should be quite sparse.
    inc_mat = np.zeros([len(unique_psm), len(unique_proteins)])
    for protein_hypothesis in ident_data.protein_hypotheses:
        db_sequence = protein_hypothesis.db_sequence_id
        i = np.where(unique_proteins == db_sequence)[0][0]
        for spectrum_id in protein_hypothesis.spectrum_ids:
            j = np.where(unique_psm == spectrum_id)[0][0]
            inc_mat[j, i] = 1
    return (unique_proteins, unique_psm, inc_mat)


def razor_proteins(unique_proteins, unique_psm, inc_mat):
    # Resolve shared peptides by the Occam's Razor approach.
    # 1.- Sort proteins by number of associated PSM (Descendent).
    number_of_psm_per_protein = inc_mat.sum(axis=0)
    sort_index = np.argsort(number_of_psm_per_protein)[::-1]
    unique_proteins = unique_proteins[sort_index]
    inc_mat = inc_mat[:, sort_index]
    for i in range(0, len(unique_proteins) - 1):

        # FIXME: If we were to be correct, we should reorder the matrix after
        # each iteration. This is computationally very expensive for this
        # prototype function. A better approach should be used for the C++
        # version.

        # 2.- Greedyly assign PSMs to the first protein they occur and remove
        #   PSM from the incidence matrix for the rest of proteins.
        inc_mat[np.where(inc_mat[:, i] == 1)[0], (i + 1):] = 0

    return (unique_proteins, unique_psm, inc_mat)


def psm_db_sequences(ident_data):
    unique_proteins, unique_psm, inc_mat = create_psm_protein_graph(
        ident_data)
    unique_proteins, unique_psm, inc_mat = razor_proteins(
        unique_proteins, unique_psm, inc_mat)
    db_sequences = []
    for psm in ident_data.spectrum_ids:
        unique_psm_index = np.where(psm.id == unique_psm)[0]
        if len(unique_psm_index) == 0:
            db_sequences += [""]
        else:
            unique_psm_index = unique_psm_index[0]
            db_sequence_id = unique_proteins[inc_mat[unique_psm_index, :] == 1][0]
            db_sequences += [db_sequence_id]
    db_sequences_df = pd.DataFrame(
        {
            "protein_id":
                [db_sequence.id for db_sequence in ident_data.db_sequences],
            "protein_name":
                [db_sequence.value for db_sequence in ident_data.db_sequences],
        })
    db_sequences = pd.DataFrame({"protein_id": db_sequences})
    db_sequences_df = pd.merge(db_sequences, db_sequences_df, how='left')
    db_sequences_df['psm_id'] = [psm.id for psm in ident_data.spectrum_ids]

    return db_sequences_df


def default_parameters(instrument, avg_fwhm_rt):
    if instrument == 'orbitrap':
        tapp_parameters = {
            'instrument_type': 'orbitrap',
            'resolution_ms1': 70000,
            'resolution_msn': 30000,
            'reference_mz': 200,
            'avg_fwhm_rt': avg_fwhm_rt,
            # Meshing.
            'num_samples_mz': 5,
            'num_samples_rt': 5,
            'smoothing_coefficient_mz': 0.4,
            'smoothing_coefficient_rt': 0.4,
            # Warp2D.
            'warp2d_slack': 30,
            'warp2d_window_size': 50,
            'warp2d_num_points': 2000,
            'warp2d_rt_expand_factor': 0.2,
            'warp2d_peaks_per_window': 100,
            # MetaMatch.
            'metamatch_radius_mz': 0.005,
            'metamatch_radius_rt': avg_fwhm_rt/2,
            'metamatch_fraction': 0.7,
            # Other
            'max_peaks': 100000,
            'polarity': 'both',
            'min_mz': 0,
            'max_mz': 100000,
            'min_rt': 0,
            'max_rt': 100000,
            # Quality.
            'similarity_num_peaks': 2000,
        }
        return tapp_parameters


def dda_pipeline(
        tapp_parameters,
        input_files,
        output_dir="TAPP",
        override_existing=False,
        save_mesh=False):
    # TODO: Logger should have different levels and user can configure the
    # verbosity of output.
    # TODO: Sanitize parameters.
    # TODO: Sanitize input/outputs.
    # TODO:     - Check if file names exist.
    # TODO:     - Check if there are name conflicts.
    # TODO:     - Check that input extension is valid.
    # TODO:     - Check if we have permission to write on output directory.
    # Create lists of files, and groups.
    input_raw_files = []
    input_stems = []
    input_ident_files = []
    groups = []
    for key in (input_files.keys()):
        input_raw_files += [key]
        base_name = os.path.basename(key)
        base_name = os.path.splitext(base_name)
        stem = base_name[0]
        input_stems += [stem]

        # TODO: Check that all files contain a group, if not, assign a default
        # group distinct from the rest.
        groups += [input_files[key]['group']]

        # Check that all files contain a ident_path, if not, assign 'none'.
        if input_files[key]['ident_path']:
            input_ident_files += [input_files[key]['ident_path']]
        else:
            input_ident_files += ['none']

    # Sort input files by groups and stems.
    groups, input_stems, input_ident_files, input_raw_files = list(
        zip(*sorted(zip(
            groups, input_stems, input_ident_files, input_raw_files))))

    # Create output directory and subdirectoreis if necessary.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir, 'raw')):
        os.makedirs(os.path.join(output_dir, 'raw'))
    if not os.path.exists(os.path.join(output_dir, 'quality')):
        os.makedirs(os.path.join(output_dir, 'quality'))
    if not os.path.exists(os.path.join(output_dir, 'mesh') and save_mesh):
        os.makedirs(os.path.join(output_dir, 'mesh'))
    if not os.path.exists(os.path.join(output_dir, 'peaks')):
        os.makedirs(os.path.join(output_dir, 'peaks'))
    if not os.path.exists(os.path.join(output_dir, 'warped_peaks')):
        os.makedirs(os.path.join(output_dir, 'warped_peaks'))
    if not os.path.exists(os.path.join(output_dir, 'metamatch')):
        os.makedirs(os.path.join(output_dir, 'metamatch'))
    if not os.path.exists(os.path.join(output_dir, 'linking')):
        os.makedirs(os.path.join(output_dir, 'linking'))
    if not os.path.exists(os.path.join(output_dir, 'ident')):
        os.makedirs(os.path.join(output_dir, 'ident'))
    if not os.path.exists(os.path.join(output_dir, 'features')):
        os.makedirs(os.path.join(output_dir, 'features'))
    if not os.path.exists(os.path.join(output_dir, 'quant')):
        os.makedirs(os.path.join(output_dir, 'quant'))

    # TODO: Initialize summary, log and parameters files.
    import logging
    import time
    import datetime

    class DeltaTimeFilter(logging.Filter):
        def filter(self, record):
            current_time = time.time()
            record.delta_time = datetime.timedelta(
                seconds=current_time - self.prev_time)
            self.prev_time = current_time
            return True

        def __init__(self):
            self.prev_time = time.time()

    logger = logging.getLogger('pipeline')
    logger.addFilter(DeltaTimeFilter())
    logger.setLevel(logging.INFO)
    logger_fh = logging.FileHandler(os.path.join(output_dir, 'info.log'))
    logger_fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s | %(delta_time)s | %(message)s')
    logger_fh.setFormatter(formatter)
    logger.addHandler(logger_fh)

    time_pipeline_start = time.time()

    # Raw data to binary conversion.
    logger.info('Starting raw data conversion')
    time_start = time.time()
    for i, file_name in enumerate(input_raw_files):
        # Check if file has already been processed.
        stem = input_stems[i]
        out_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        # Read raw files (MS1).
        logger.info('Reading MS1: {}'.format(file_name))
        raw_data = tapp.read_mzxml(
            file_name,
            min_mz=tapp_parameters['min_mz'],
            max_mz=tapp_parameters['max_mz'],
            min_rt=tapp_parameters['min_rt'],
            max_rt=tapp_parameters['max_rt'],
            instrument_type=tapp_parameters['instrument_type'],
            resolution_ms1=tapp_parameters['resolution_ms1'],
            resolution_msn=tapp_parameters['resolution_msn'],
            reference_mz=tapp_parameters['reference_mz'],
            fwhm_rt=tapp_parameters['avg_fwhm_rt'],
            polarity=tapp_parameters['polarity'],
            ms_level=1,
        )

        # Write raw_data to disk (MS1).
        logger.info('Writing MS1: {}'.format(out_path))
        raw_data.dump(out_path)

    for i, file_name in enumerate(input_raw_files):
        # Check if file has already been processed.
        stem = input_stems[i]
        out_path = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        # Read raw files (MS2).
        logger.info('Reading MS2: {}'.format(file_name))
        raw_data = tapp.read_mzxml(
            file_name,
            min_mz=tapp_parameters['min_mz'],
            max_mz=tapp_parameters['max_mz'],
            min_rt=tapp_parameters['min_rt'],
            max_rt=tapp_parameters['max_rt'],
            instrument_type=tapp_parameters['instrument_type'],
            resolution_ms1=tapp_parameters['resolution_ms1'],
            resolution_msn=tapp_parameters['resolution_msn'],
            reference_mz=tapp_parameters['reference_mz'],
            fwhm_rt=tapp_parameters['avg_fwhm_rt'],
            polarity=tapp_parameters['polarity'],
            ms_level=2,
        )

        # Write raw_data to disk (MS2).
        logger.info('Writing MS2: {}'.format(out_path))
        raw_data.dump(out_path)

    logger.info('Finished raw data conversion in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform resampling/smoothing and peak detection and save results to disk.
    logger.info('Starting peak detection')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        out_path = os.path.join(output_dir, 'peaks', "{}.bpks".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading raw_data from disk: {}".format(stem))
        raw_data = tapp.read_raw_data(in_path)

        logger.info("Resampling: {}".format(stem))
        mesh = tapp.resample(
            raw_data,
            tapp_parameters['num_samples_mz'],
            tapp_parameters['num_samples_rt'],
            tapp_parameters['smoothing_coefficient_mz'],
            tapp_parameters['smoothing_coefficient_rt'],
        )

        if save_mesh:
            logger.info('Writing mesh: {}'.format(out_path))
            mesh.dump(out_path)

        logger.info("Finding peaks: {}".format(stem))
        peaks = tapp.find_peaks(raw_data, mesh, tapp_parameters['max_peaks'])
        logger.info('Writing peaks:'.format(out_path))
        tapp.write_peaks(peaks, out_path)

    logger.info('Finished peak detection in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Calculate similarity matrix before alignment, generate heatmap and save
    # to disk.
    out_path = os.path.join(output_dir, 'quality', 'unwarped_similarity')
    logger.info("Starting unwarped similarity matrix calculation")
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or override_existing:
        similarity_matrix = np.zeros(
            len(input_stems) ** 2).reshape(len(input_stems), len(input_stems))
        for i in range(0, len(input_stems)):
            stem_a = input_stems[i]
            logger.info("Reading peaks_a from disk: {}".format(stem_a))
            peaks_a = tapp.read_peaks(os.path.join(
                output_dir, 'peaks', '{}.bpks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'peaks', '{}.bpks'.format(stem_b)))
                logger.info(
                    "Calculating similarity of {} vs {}".format(
                        stem_a, stem_b))
                similarity_matrix[j, i] = tapp.find_similarity(
                    peaks_a, peaks_b,
                    tapp_parameters['similarity_num_peaks']).geometric_ratio
                similarity_matrix[i, j] = similarity_matrix[j, i]
        similarity_matrix = pd.DataFrame(similarity_matrix)
        similarity_matrix_names = [input_stem.split(
            '.')[0] for input_stem in input_stems]
        similarity_matrix.columns = similarity_matrix_names
        similarity_matrix.rename(index=dict(zip(
            range(0, len(similarity_matrix_names), 1),
            similarity_matrix_names)), inplace=True)
        plt.ioff()
        fig = plt.figure()
        sns.heatmap(similarity_matrix, xticklabels=True,
                    yticklabels=True, square=True, vmin=0, vmax=1)
        # Save similarity matrix and figure to disk.
        logger.info("Saving similarity matrix: {}.csv".format(out_path))
        similarity_matrix.to_csv("{}.csv".format(out_path))
        # TODO: Use plot saving from utilities library.
        fig.set_size_inches(7.5 * 16/9, 7.5)
        plt.savefig("{}.png".format(out_path), dpi=100)
        plt.close(fig)
        logger.info("Saving similarity matrix plot: {}.png".format(out_path))
    logger.info('Finished unwarped similarity matrix calculation in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Correct retention time. If a reference sample is selected it will be
    # used, otherwise, exhaustive warping will be performed.
    # TODO: Allow usage of reference sample.
    out_path = os.path.join(output_dir, 'quality',
                            'exhaustive_warping_similarity')
    reference_index = 0
    similarity_matrix = np.zeros(
        len(input_stems) ** 2).reshape(len(input_stems), len(input_stems))
    logger.info("Starting exhaustive warping similarity calculation")
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or override_existing:
        for i in range(0, len(input_stems)):
            stem_a = input_stems[i]
            logger.info("Reading peaks_a from disk: {}".format(stem_a))
            peaks_a = tapp.read_peaks(os.path.join(
                output_dir, 'peaks', '{}.bpks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'peaks', '{}.bpks'.format(stem_b)))
                logger.info("Warping {} peaks to {}".format(stem_b, stem_a))
                peaks_b = tapp.warp_peaks(
                    [peaks_a, peaks_b],
                    0,
                    tapp_parameters['warp2d_slack'],
                    tapp_parameters['warp2d_window_size'],
                    tapp_parameters['warp2d_num_points'],
                    tapp_parameters['warp2d_rt_expand_factor'],
                    tapp_parameters['warp2d_peaks_per_window'])[1]
                logger.info(
                    "Calculating similarity of {} vs {} (warped)".format(
                        stem_a, stem_b))
                similarity_matrix[j, i] = tapp.find_similarity(
                    peaks_a, peaks_b,
                    tapp_parameters['similarity_num_peaks']).geometric_ratio
                similarity_matrix[i, j] = similarity_matrix[j, i]
        similarity_matrix = pd.DataFrame(similarity_matrix)
        similarity_matrix_names = [input_stem.split(
            '.')[0] for input_stem in input_stems]
        similarity_matrix.columns = similarity_matrix_names
        similarity_matrix.rename(index=dict(zip(
            range(0, len(similarity_matrix_names), 1),
            similarity_matrix_names)), inplace=True)
        plt.ioff()
        fig = plt.figure()
        sns.heatmap(similarity_matrix, xticklabels=True,
                    yticklabels=True, square=True, vmin=0, vmax=1)
        # Save similarity matrix and figure to disk.
        logger.info("Saving similarity matrix: {}.csv".format(out_path))
        similarity_matrix.to_csv("{}.csv".format(out_path))
        # TODO: Use plot saving from utilities library.
        fig.set_size_inches(7.5 * 16/9, 7.5)
        plt.savefig("{}.png".format(out_path), dpi=100)
        plt.close(fig)
        logger.info("Saving similarity matrix plot: {}.png".format(out_path))
    else:
        # Load exhaustive_warping_similarity to calculate the reference idx.
        similarity_matrix = pd.read_csv("{}.csv".format(out_path), index_col=0)
        reference_index = similarity_matrix.sum(axis=0).values.argmax()
    logger.info(
        'Finished exhaustive warping similarity calculation in {}'.format(
            datetime.timedelta(seconds=time.time()-time_start)))

    # Warp all peaks to the reference file.
    reference_stem = input_stems[reference_index]
    logger.info("Starting peak warping to reference ({})".format(
        reference_stem))
    time_start = time.time()
    logger.info("Reading reference peaks")
    reference_peaks = tapp.read_peaks(os.path.join(
        output_dir, 'peaks', '{}.bpks'.format(reference_stem)))
    for stem in input_stems:
        # Check if file has already been processed.
        in_path = os.path.join(output_dir, 'peaks', "{}.bpks".format(stem))
        out_path = os.path.join(output_dir, 'warped_peaks',
                                "{}.bpks".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Warping peaks: {}".format(stem))
        if stem == reference_stem:
            tapp.write_peaks(reference_peaks, out_path)
        else:
            logger.info("Reading peaks from disk: {}".format(stem))
            peaks = tapp.read_peaks(in_path)
            peaks = tapp.warp_peaks(
                [reference_peaks, peaks],
                0,
                tapp_parameters['warp2d_slack'],
                tapp_parameters['warp2d_window_size'],
                tapp_parameters['warp2d_num_points'],
                tapp_parameters['warp2d_rt_expand_factor'],
                tapp_parameters['warp2d_peaks_per_window'])[1]
            tapp.write_peaks(peaks, out_path)
    logger.info('Finished peak warping to reference ({}) in {}'.format(
        reference_stem, datetime.timedelta(seconds=time.time()-time_start)))

    # Calculate similarity matrix after alignment, generate heatmap and save to
    # disk.
    out_path = os.path.join(output_dir, 'quality', 'warped_similarity')
    logger.info("Starting warped similarity matrix calculation")
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or override_existing:
        similarity_matrix = np.zeros(
            len(input_stems) ** 2).reshape(len(input_stems), len(input_stems))
        for i in range(0, len(input_stems)):
            stem_a = input_stems[i]
            logger.info("Reading peaks_a from disk: {}".format(stem_a))
            peaks_a = tapp.read_peaks(os.path.join(
                output_dir, 'warped_peaks', '{}.bpks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'warped_peaks', '{}.bpks'.format(stem_b)))
                logger.info(
                    "Calculating similarity of {} vs {}".format(
                        stem_a, stem_b))
                similarity_matrix[j, i] = tapp.find_similarity(
                    peaks_a, peaks_b,
                    tapp_parameters['similarity_num_peaks']).geometric_ratio
                similarity_matrix[i, j] = similarity_matrix[j, i]
        similarity_matrix = pd.DataFrame(similarity_matrix)
        similarity_matrix_names = [input_stem.split(
            '.')[0] for input_stem in input_stems]
        similarity_matrix.columns = similarity_matrix_names
        similarity_matrix.rename(index=dict(zip(
            range(0, len(similarity_matrix_names), 1),
            similarity_matrix_names)), inplace=True)
        plt.ioff()
        fig = plt.figure()
        sns.heatmap(similarity_matrix, xticklabels=True,
                    yticklabels=True, square=True, vmin=0, vmax=1)
        # Save similarity matrix and figure to disk.
        logger.info("Saving similarity matrix: {}.csv".format(out_path))
        similarity_matrix.to_csv("{}.csv".format(out_path))
        # TODO: Use plot saving from utilities library.
        fig.set_size_inches(7.5 * 16/9, 7.5)
        plt.savefig("{}.png".format(out_path), dpi=100)
        plt.close(fig)
        logger.info("Saving similarity matrix plot: {}.png".format(out_path))
    logger.info('Finished warped similarity matrix calculation in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Use metamatch to match warped peaks.
    logger.info("Starting metamatch")
    time_start = time.time()
    out_path = os.path.join(output_dir, 'metamatch')
    if (not os.path.exists(os.path.join(out_path, "metamatch.clusters"))
            or override_existing):
        metamatch_input = []
        for i, stem in enumerate(input_stems):
            in_path = os.path.join(
                output_dir, 'warped_peaks', "{}.bpks".format(stem))
            logger.info("Reading peaks from disk: {}".format(stem))
            metamatch_input += [(groups[i], tapp.read_peaks(in_path))]

        metamatch_results = tapp.perform_metamatch(
            metamatch_input,
            tapp_parameters['metamatch_radius_mz'],
            tapp_parameters['metamatch_radius_rt'],
            tapp_parameters['metamatch_fraction'])

        # Save metamatch results to disk.
        logger.info("Writing metamatch results to disk")
        tapp.write_metamatch_clusters(
            metamatch_results.clusters,
            os.path.join(out_path, "metamatch.clusters"))
        tapp.write_metamatch_peaks(
            metamatch_results.orphans,
            os.path.join(out_path, "metamatch.orphans"))

    logger.info('Finished metamatch in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Match ms2 events with corresponding detected peaks.
    logger.info('Starting peaks/msms linkage')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.bpks".format(stem))
        out_path = os.path.join(output_dir, 'linking',
                                "{}.ms2_peaks.link".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading raw_data from disk (MS2): {}".format(stem))
        raw_data = tapp.read_raw_data(in_path_raw)

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Performing linkage: {}".format(stem))
        linked_msms = tapp.link_peaks(peaks, raw_data)
        logger.info('Writing linked_msms: {}'.format(out_path))
        tapp.write_linked_msms(linked_msms, out_path)
    logger.info('Finished peaks/msms linkage in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # TODO: (If there is ident information)
    # TODO:     - Perform Occam's razor protein inference in linked peptides.
    # Read mzidentdata and save binary data to disk.
    logger.info('Starting mzIdentML parsing')
    time_start = time.time()
    for i, stem in enumerate(input_stems):
        in_path = input_ident_files[i]
        out_path = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        if in_path == 'none' or (os.path.exists(out_path) and not override_existing):
            continue
        logger.info('Reading mzIdentML: {}'.format(in_path))
        ident_data = tapp.read_mzidentml(in_path)
        logger.info('Writing ident data: {}'.format(out_path))
        tapp.write_ident_data(ident_data, out_path)
    logger.info('Finished mzIdentML parsing in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform protein inference
    logger.info('Starting protein inference')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_idents = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        out_path = os.path.join(output_dir, 'ident',
                                "{}.inferred_proteins".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading ident from disk: {}".format(stem))
        ident_data = tapp.read_ident_data(in_path_idents)

        logger.info("Performing inference: {}".format(stem))
        inferred_proteins = tapp.perform_protein_inference(ident_data)
        logger.info('Writing inferred_proteins: {}'.format(out_path))
        tapp.write_inferred_proteins(inferred_proteins, out_path)
    logger.info('Finished protein inference in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Link ms2 events with ident information.
    logger.info('Starting ident/msms linkage')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_idents = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        out_path = os.path.join(output_dir, 'linking',
                                "{}.ms2_idents.link".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading raw_data from disk (MS2): {}".format(stem))
        raw_data = tapp.read_raw_data(in_path_raw)

        logger.info("Reading ident from disk: {}".format(stem))
        ident_data = tapp.read_ident_data(in_path_idents)

        logger.info("Performing linkage: {}".format(stem))
        linked_idents = tapp.link_idents(ident_data, raw_data)
        logger.info('Writing linked_msms: {}'.format(out_path))
        tapp.write_linked_msms(linked_idents, out_path)
    logger.info('Finished ident/msms linkage in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform feature detection using averagine or linked identification
    # if available in ms2 linked peaks.
    logger.info('Starting feature detection')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.bpks".format(stem))
        in_path_idents = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        in_path_ident_link = os.path.join(output_dir, 'linking',
                                          "{}.ms2_idents.link".format(stem))
        in_path_peaks_link = os.path.join(output_dir, 'linking',
                                          "{}.ms2_peaks.link".format(stem))
        out_path = os.path.join(output_dir, 'features',
                                "{}.features".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading raw_data from disk (MS2): {}".format(stem))
        raw_data = tapp.read_raw_data(in_path_raw)

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Reading ident_data from disk: {}".format(stem))
        ident_data = tapp.read_ident_data(in_path_idents)

        logger.info("Reading linked peaks from disk: {}".format(stem))
        linked_peaks = tapp.read_linked_msms(in_path_peaks_link)

        logger.info("Reading linked idents from disk: {}".format(stem))
        linked_idents = tapp.read_linked_msms(in_path_ident_link)

        logger.info("Performing feature_detection: {}".format(stem))
        features = tapp.feature_detection(
            peaks, raw_data, ident_data, linked_peaks, linked_idents)
        logger.info('Writing features: {}'.format(out_path))
        tapp.write_features(features, out_path)
    logger.info('Finished feature detection in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform metamatch on detected features.
    logger.info('Starting metamatch on features')
    time_start = time.time()
    in_path_features = os.path.join(output_dir, 'features')
    in_path_peaks = os.path.join(output_dir, 'warped_peaks')
    out_path = os.path.join(output_dir, 'metamatch', "features.clusters")
    if (not os.path.exists(out_path) or override_existing):
        logger.info("Reading features from disk")
        features = [
            tapp.read_features(
                os.path.join(in_path_features, "{}.features".format(input_stem)))
            for input_stem in input_stems]

        logger.info("Reading peaks from disk")
        peaks = [
            tapp.read_peaks(
                os.path.join(in_path_peaks, "{}.bpks".format(input_stem)))
            for input_stem in input_stems]

        logger.info("Finding feature clusters")
        feature_clusters = tapp.find_feature_clusters(groups, peaks, features)

        logger.info("Writing feature clusters to disk")
        tapp.write_feature_clusters(feature_clusters, out_path)

    logger.info('Finished metamatch on features in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Create final quantitative tables.
    logger.info('Starting creation of quantitative tables')
    time_start = time.time()
    for stem in input_stems:
        # Peaks
        # =====
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.bpks".format(stem))
        out_path_peaks = os.path.join(output_dir, 'quant',
                                      "{}_peaks.csv".format(stem))
        if os.path.exists(out_path_peaks) and not override_existing:
            continue

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Generating peaks quantitative table")
        peaks_df = pd.DataFrame({
            'peak_id': [peak.id for peak in peaks],
            'mz': [peak.local_max_mz for peak in peaks],
            'rt': [peak.local_max_rt for peak in peaks],
            'rt_delta': [peak.rt_delta for peak in peaks],
            'height': [peak.local_max_height for peak in peaks],
            'total_intensity': [peak.raw_roi_total_intensity for peak in peaks],
            'roi_min_mz': [peak.roi_min_mz for peak in peaks],
            'roi_max_mz': [peak.roi_max_mz for peak in peaks],
            'roi_min_rt': [peak.roi_min_rt for peak in peaks],
            'roi_max_rt': [peak.roi_max_rt for peak in peaks],
            'sigma_mz': [peak.raw_roi_sigma_mz for peak in peaks],
            'sigma_rt': [peak.raw_roi_sigma_rt for peak in peaks],
            'skewness_mz': [peak.raw_roi_skewness_mz for peak in peaks],
            'skewness_rt': [peak.raw_roi_skewness_rt for peak in peaks],
            'kurtosis_mz': [peak.raw_roi_kurtosis_mz for peak in peaks],
            'kurtosis_rt': [peak.raw_roi_kurtosis_rt for peak in peaks],
            'num_points': [peak.raw_roi_num_points for peak in peaks],
            'num_scans': [peak.raw_roi_num_scans for peak in peaks],
        })
        peaks_df.to_csv(out_path_peaks, index=False)

    for stem in input_stems:
        # Features
        # ========
        in_path_features = os.path.join(
            output_dir, 'features', "{}.features".format(stem))
        out_path_features = os.path.join(output_dir, 'quant',
                                         "{}_features.csv".format(stem))
        if os.path.exists(out_path_features) and not override_existing:
            continue

        logger.info("Reading features from disk: {}".format(stem))
        features = tapp.read_features(in_path_features)

        logger.info("Generating features quantitative table")
        features_df = pd.DataFrame({
            'feature_id': [feature.id for feature in features],
            'average_mz': [feature.average_mz for feature in features],
            'average_mz_sigma': [feature.average_mz_sigma for feature in features],
            'average_rt': [feature.average_rt for feature in features],
            'average_rt_sigma': [feature.average_rt_sigma for feature in features],
            'average_rt_delta': [feature.average_rt_delta for feature in features],
            'total_height': [feature.total_height for feature in features],
            'monoisotopic_mz': [feature.monoisotopic_mz for feature in features],
            'monoisotopic_height': [feature.monoisotopic_height for feature in features],
            'peak_ids': [str(feature.peak_ids) for feature in features],
        })
        features_df.to_csv(out_path_features, index=False)

    # Matched Peaks
    # =============
    logger.info("Reading peak clusters from disk")
    in_path_peak_clusters = os.path.join(
        output_dir, 'metamatch', 'metamatch.clusters')
    out_path_peak_clusters = os.path.join(output_dir, 'quant',
                                          "peak_clusters.csv")
    if (not os.path.exists(out_path_peak_clusters) or override_existing):
        peak_clusters = tapp.read_metamatch_clusters(in_path_peak_clusters)
        logger.info("Generating peak clusters quantitative table")
        peak_clusters_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in peak_clusters],
            'mz': [cluster.mz for cluster in peak_clusters],
            'rt': [cluster.rt for cluster in peak_clusters],
            'avg_height': [cluster.avg_height for cluster in peak_clusters],
        })
        for i, stem in enumerate(input_stems):
            peak_clusters_df[stem] = [cluster.file_heights[i]
                                      for cluster in peak_clusters]
        peak_clusters_df.to_csv(out_path_peak_clusters, index=False)

    # Matched Features
    # ================
    logger.info("Reading feature clusters from disk")
    in_path_feature_clusters = os.path.join(
        output_dir, 'metamatch', 'features.clusters')
    out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                             "feature_clusters.csv")
    if (not os.path.exists(out_path_feature_clusters) or override_existing):
        feature_clusters = tapp.read_feature_clusters(
            in_path_feature_clusters)
        logger.info("Generating feature clusters quantitative table")
        feature_clusters_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in feature_clusters],
            'mz': [cluster.mz for cluster in feature_clusters],
            'rt': [cluster.rt for cluster in feature_clusters],
            'avg_height': [cluster.avg_height for cluster in feature_clusters],
        })
        for i, stem in enumerate(input_stems):
            feature_clusters_df[stem] = [cluster.file_heights[i]
                                         for cluster in feature_clusters]
        feature_clusters_df.to_csv(out_path_feature_clusters, index=False)

    # Find protein information for the identified features.
    logger.info("Finding protein information on identified features")
    for stem in input_stems:
        in_path_features = os.path.join(
            output_dir, 'features', "{}.features".format(stem))
        in_path_ident_link = os.path.join(output_dir, 'linking',
                                          "{}.ms2_idents.link".format(stem))
        in_path_ident_data = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        in_path_inferred_proteins = os.path.join(
            output_dir, 'ident', "{}.inferred_proteins".format(stem))
        out_path_identified_features = os.path.join(
            output_dir, 'quant', "{}_identified_features.csv".format(stem))
        if os.path.exists(out_path_identified_features) and not override_existing:
            continue

        logger.info("Reading features from disk: {}".format(stem))
        features = tapp.read_features(in_path_features)
        features_df = pd.DataFrame({
            'feature_id': [feature.id for feature in features],
            'msms_index': [feature.msms_id for feature in features],
        })

        logger.info("Reading ident_data from disk: {}".format(stem))
        ident_data = tapp.read_ident_data(in_path_ident_data)

        logger.info("Reading inferred_proteins from disk: {}".format(stem))
        inferred_proteins = tapp.read_inferred_proteins(
            in_path_inferred_proteins)
        inferred_proteins_df = pd.DataFrame({
            'protein_id': [inferred_protein.protein_id for inferred_protein in inferred_proteins],
            'psm_id': [inferred_protein.psm_id for inferred_protein in inferred_proteins],
        })

        logger.info("Reading linked idents from disk: {}".format(stem))
        linked_idents = tapp.read_linked_msms(in_path_ident_link)
        linked_idents_df = pd.DataFrame({
            'psm_index': [linked_ident.entity_id for linked_ident in linked_idents],
            'msms_index': [linked_ident.msms_id for linked_ident in linked_idents],
        })

        logger.info("Creating linked psm protein table: {}".format(stem))
        linked_features_df = pd.merge(
            features_df, linked_idents_df, on="msms_index", how="left")
        linked_spectrum_ids = [ident_data.spectrum_ids[i]
                               for i in linked_features_df['psm_index']]
        linked_spectrum_ids_df = pd.DataFrame({
            'psm_id': [linked_spectrum_id.id for linked_spectrum_id in linked_spectrum_ids],
            'sequence': [linked_spectrum_id.sequence for linked_spectrum_id in linked_spectrum_ids],
            'charge_state': [linked_spectrum_id.charge_state for linked_spectrum_id in linked_spectrum_ids],
            'retention_time': [linked_spectrum_id.retention_time for linked_spectrum_id in linked_spectrum_ids],
        })
        linked_spectrum_ids_df = pd.merge(
            linked_spectrum_ids_df, inferred_proteins_df, on="psm_id", how="left")

        protein_names_df = pd.DataFrame({
            'protein_id': [db_seq.id for db_seq in ident_data.db_sequences],
            'protein_name': [db_seq.value for db_seq in ident_data.db_sequences],
        })
        linked_spectrum_ids_df = pd.merge(
            linked_spectrum_ids_df, protein_names_df, on="protein_id", how="left")

        linked_features_df = pd.concat(
            [linked_features_df, linked_spectrum_ids_df], axis=1)
        linked_features_df.to_csv(out_path_identified_features, index=False)

    # Link metamatch clusters and corresponding features with identification
    # information of peptides and proteins.
    in_path_feature_clusters = os.path.join(
        output_dir, 'metamatch', 'features.clusters')
    out_path_feature_clusters_info_feature_ids = os.path.join(output_dir, 'quant',
                                                              "feature_clusters_info_feature_ids.csv")
    out_path_feature_clusters_info_sequence = os.path.join(output_dir, 'quant',
                                                           "feature_clusters_info_sequence.csv")
    out_path_feature_clusters_info_protein_id = os.path.join(output_dir, 'quant',
                                                             "feature_clusters_info_protein_id.csv")
    out_path_feature_clusters_info_protein_name = os.path.join(output_dir, 'quant',
                                                               "feature_clusters_info_protein_name.csv")
    if (not os.path.exists(out_path_feature_clusters_info_feature_ids) or override_existing):
        logger.info("Reading feature clusters from disk")
        feature_clusters = tapp.read_feature_clusters(
            in_path_feature_clusters)

        logger.info("Generating feature clusters identification table")
        peptide_sequence_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in feature_clusters],
        })
        protein_id_df = peptide_sequence_df.copy()
        protein_name_df = peptide_sequence_df.copy()

        logger.info("Generating cluster to feature_id table")
        cluster_features_df = []
        for cluster in feature_clusters:
            row = np.full(len(input_stems), np.nan)
            for feature_id in cluster.feature_ids:
                row[feature_id.file_id] = feature_id.feature_id
            cluster_features_df += [row]
        cluster_features_df = pd.DataFrame(cluster_features_df, dtype='Int64')
        cluster_features_df.columns = input_stems
        cluster_features_df = pd.concat(
            [peptide_sequence_df.copy(), cluster_features_df], axis=1)

        for i, stem in enumerate(input_stems):
            logger.info(
                "Reading identified_features from disk: {}".format(stem))
            in_path_identified_features = os.path.join(
                output_dir, 'quant', "{}_identified_features.csv".format(stem))
            identified_features = pd.read_csv(in_path_identified_features)
            identified_features['feature_id'] = identified_features['feature_id'].astype(
                'Int64')
            cluster_features_info = pd.merge(
                pd.DataFrame({'feature_id': cluster_features_df[stem]}),
                identified_features,
                how='left'
            )
            peptide_sequence_df[stem] = cluster_features_info['sequence']
            protein_id_df[stem] = cluster_features_info['protein_id']
            protein_name_df[stem] = cluster_features_info['protein_name']

        peptide_sequence_df.to_csv(
            out_path_feature_clusters_info_sequence, index=False)
        protein_id_df.to_csv(
            out_path_feature_clusters_info_protein_id, index=False)
        protein_name_df.to_csv(
            out_path_feature_clusters_info_protein_name, index=False)
        cluster_features_df.to_csv(
            out_path_feature_clusters_info_feature_ids, index=False)

    logger.info('Finished creation of quantitative tables in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    logger.info('Total time elapsed: {}'.format(
        datetime.timedelta(seconds=time.time()-time_pipeline_start)))

    # Stop logger.
    logger.removeHandler(logger_fh)
    logger_fh.close()
    return


def full_dda_pipeline_test():
    input_files = {
        '/data/HYE_DDA_Orbitrap/mzXML/subset/3_1.mzXML': {'group': 3, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/3_1.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/3_2.mzXML': {'group': 3, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/3_2.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/1_1.mzXML': {'group': 1, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/1_1.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/1_2.mzXML': {'group': 1, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/1_2.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/1_3.mzXML': {'group': 1, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/1_3.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/1_4.mzXML': {'group': 1, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/1_4.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/1_5.mzXML': {'group': 1, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/1_5.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/3_3.mzXML': {'group': 3, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/3_3.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/3_4.mzXML': {'group': 3, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/3_4.mzid'},
        '/data/HYE_DDA_Orbitrap/mzXML/subset/3_5.mzXML': {'group': 3, 'ident_path': '/data/HYE_DDA_Orbitrap/mzXML/subset/identification/mzid/3_5.mzid'},
    }
    tapp_parameters = default_parameters('orbitrap', 9)
    tapp_parameters['max_peaks'] = 1000
    tapp_parameters['polarity'] = 'pos'
    # TODO: Can we optimize the default parameters for warp2d based on
    # avg_fwhm_rt and min/max rt?
    tapp_parameters['warp2d_slack'] = 100
    tapp_parameters['warp2d_window_size'] = 100
    tapp_parameters['warp2d_peaks_per_window'] = 50

    dda_pipeline(tapp_parameters, input_files, 'tapp_pipeline_test')


def testing_feature_detection():
    peaks = tapp.read_peaks('tapp_pipeline_test/warped_peaks/1_3.bpks')
    ms2_data = tapp.read_raw_data('tapp_pipeline_test/raw/1_3.ms2')
    linked_peaks = tapp.read_linked_msms(
        'tapp_pipeline_test/linking/1_3.ms2_peaks.link')
    linked_idents = tapp.read_linked_msms(
        'tapp_pipeline_test/linking/1_3.ms2_idents.link')
    ident_data = tapp.read_ident_data('tapp_pipeline_test/ident/1_3.ident')
    results = tapp.feature_detection(
        peaks, ms2_data, ident_data, linked_peaks, linked_idents)
    return (peaks, ms2_data, ident_data, linked_peaks, linked_idents, results)


def testing_feature_matching():
    file_names = ["1_1", "1_2", "1_3", "1_4",
                  "1_5", "3_1", "3_2", "3_3", "3_4", "3_5"]
    features = [
        tapp.read_features(
            "tapp_pipeline_test/features/{}.features".format(file_name))
        for file_name in file_names]
    peaks = [
        tapp.read_peaks(
            "tapp_pipeline_test/warped_peaks/{}.bpks".format(file_name))
        for file_name in file_names]
    groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    return tapp.find_feature_clusters(groups, peaks, features)


def peak_xic(peak, raw_data, method="sum"):
    return tapp.xic(raw_data, peak.roi_min_mz, peak.roi_max_mz,
                    peak.roi_min_rt, peak.roi_max_rt, method)


Peak.xic = peak_xic
Peak.plot_xic = plot_xic
Peak.plot_raw_points = plot_raw_points
Peak.plot_raw_roi_sigma = plot_raw_roi_sigma
Peak.plot_theoretical_sigma = plot_theoretical_sigma
Peak.plot_fitted_sigma = plot_raw_roi_fitted_sigma
Peak.plot_fitted_sigma_fast = plot_raw_roi_fitted_sigma_fast
Peak.plot_fitted_sigma_weighted = plot_raw_roi_fitted_sigma_weighted
Peak.plot_sigma = plot_sigma
Peak.fit_height_and_sigmas = fit_height_and_sigmas
Peak.fit_sigmas = fit_sigmas
