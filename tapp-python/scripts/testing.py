import tapp
import math
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

plt.style.use('dark_background')
plt.ion()
plt.show()
try:
    fig
except NameError:
    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5, forward=True)
try:
    ax
except NameError:
    ax = plt.axes()

# raw_data = tapp.raw_data_load_dump("150210_11_01.rawdump")
# tic = raw_data.tic()

# # raw_data = tapp.raw_data("/data/Spike-In QEx UKE/HD5YD15ED_DDA_R1.mzXML", 799.036, 809.107, 2768.83, 3171.63)
# raw_data_big = tapp.raw_data("/data/Spike-In QEx UKE/HD5YD15ED_DDA_R1.mzXML", 0, 2000, 1000, 5000)
# # raw_data_tof_small = tapp.raw_data("/data/180724_TripleTOF_Test_Data/Christoph_Krisp_Marcel_Kwiatkowski/171214_DDA/mzXML/HD+5YD+15ED_DDA.mzXML", 762, 774, 3038, 3210)
# # raw_data.dump("HD5YD15ED_DDA_R1_small.rawdata")
# # raw_data_big.dump("HD5YD15ED_DDA_R1_big.rawdata")
# raw_data_tof_small.dump("HD5YD15ED_DDA_R1_tof_small.rawdata")
# # raw_data = tapp.raw_data_load_dump("HD5YD15ED_DDA_R1_small.rawdata")
# raw_data_big = tapp.raw_data_load_dump("HD5YD15ED_DDA_R1_big.rawdata")
# # raw_data_tof_small = tapp.raw_data_load_dump("")

# raw_data_df = pd.DataFrame(
    # {
        # "mz": raw_data.mz(),
        # "rt": raw_data.rt(),
        # "intensity": raw_data.intensity()
    # })
# raw_data_big_df = pd.DataFrame(
    # {
        # "mz": raw_data_big.mz(),
        # "rt": raw_data_big.rt(),
        # "intensity": raw_data_big.intensity()
    # })
# raw_data_tof_small_df = pd.DataFrame(
    # {
        # "mz": raw_data_tof_small.mz(),
        # "rt": raw_data_tof_small.rt(),
        # "intensity": raw_data_tof_small.intensity()
    # })

# # img = tapp.splat_points(
    # # raw_data.data,
    # # min_mz,
    # # max_mz,
    # # min_rt,
    # # max_rt,
    # # delta_rt=1.5,
    # # delta_mz=0.01,
    # # sigma_mz=(0.0025, 200),
    # # sigma_rt=15
# # )

def sigma_at_mz(mz, fwhm_ref, mz_ref):
    # return fwhm_ref * (mz/mz_ref) ** 1  # NOTE: TOF only
    return fwhm_ref * (mz/mz_ref) ** 1.5  # NOTE: Orbitrap only

def create_heatmap(
        df,
        delta_rt=1.5,
        delta_mz=0.01,
        transform='none',
        type='binning',
        sigma_mz=(0.0025, 200),
        sigma_rt=15,
    ):
    unique_rt = df.sort_values("rt")["rt"].unique()
    min_rt = unique_rt.min()
    max_rt = unique_rt.max()
    num_bins_rt = math.floor((max_rt - min_rt)/ delta_rt) + 1
    bins_rt =  [min_rt + delta_rt * x for x in range(0, num_bins_rt)] 

    unique_mz = df.sort_values("mz")["mz"].unique()
    min_mz = unique_mz.min()
    max_mz = unique_mz.max()
    num_bins_mz = math.floor((max_mz - min_mz)/ delta_mz) + 1
    bins_mz =  [min_mz + delta_mz * x for x in range(0, num_bins_mz)] 

    img = np.zeros(num_bins_rt * num_bins_mz)
    img = np.reshape(img, (num_bins_rt, num_bins_mz))
    
    if type == 'binning':
        # Splat points in matrix
        for row in range(0, df.shape[0]):
            # Find y index
            rt = df.iloc[row, 1]
            x = math.floor((rt - min_rt) / delta_rt)
            # Find x index
            mz = df.iloc[row, 0]
            y = math.floor((mz - min_mz) / delta_mz)
            # Update matrix value
            img[x,y] = img[x,y] + df.iloc[row, 2]
    elif type == 'tapp_smoothing':
        # Splat points in matrix
        for row in range(0, df.shape[0]):
            rt = df.iloc[row, 1]
            mz = df.iloc[row, 0]
            intensity = df.iloc[row, 2]
            current_sigma_mz = sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])

            kernel_min_rt = rt - 2 * sigma_rt
            kernel_max_rt = rt + 2 * sigma_rt
            kernel_min_mz = mz - 2 * current_sigma_mz
            kernel_max_mz = mz + 2 * current_sigma_mz

            i_min = 0 if kernel_min_mz < min_mz else math.floor((kernel_min_mz - min_mz) / delta_mz)
            j_min = 0 if kernel_min_rt < min_rt else math.floor((kernel_min_rt - min_rt) / delta_rt)
            i_max = num_bins_mz - 1 if kernel_max_mz > max_mz else math.floor((kernel_max_mz - min_mz) / delta_mz)
            j_max = num_bins_rt - 1 if kernel_max_rt > max_rt else math.floor((kernel_max_rt - min_rt) / delta_rt)

            for j in range(j_min, j_max + 1):
                for i in range(i_min, i_max + 1):
                    x = bins_mz[i]
                    y = bins_rt[j]

                    a = (x - mz) / current_sigma_mz
                    b = (y - rt) / sigma_rt
                    weight = math.exp(-0.5 * (a * a + b * b))

                    img[j, i] = img[j, i] + intensity * weight
    # elif type == 'tapp_internal_smoothing':

    if transform == 'sqrt':
        img = np.sqrt(img)
    elif transform == 'log':
        img = np.log(img + 0.0001)

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
    img_plot.imshow(img, aspect='auto', origin="lower")

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

    fig.suptitle("min_rt: {0:.2f}, max_rt: {1:.2f}, min_mz: {2:.2f}, max_mz: {3:.2f}, delta_rt: {4:.2f}, delta_mz: {5:.2f}".format(
        min_rt, max_rt, min_mz, max_mz, delta_rt, delta_mz))

    return({
        "img": img,
        "bins_mz": bins_mz,
        "bins_rt": bins_rt,
        "rt_plot": rt_plot,
        "mz_plot": mz_plot,
        "img_plot": img_plot
    })

# # Binned
# # heatmap = create_heatmap(raw_data_df, transform='sqrt')
# # peptide_focus_heatmap = create_heatmap(
    # # raw_data_df[
        # # (raw_data_df['mz'] < 803)
        # # & (raw_data_df['rt'] > 2970)
        # # & (raw_data_df['rt'] < 3020)]
    # # )
# # peak_focus_heatmap = create_heatmap(
    # # raw_data_df[
        # # (raw_data_df['mz'] > 799.54)
        # # & (raw_data_df['mz'] < 800.04)
        # # & (raw_data_df['rt'] > 2970)
        # # & (raw_data_df['rt'] < 3020)]
    # # )

# # TAPP Smoothing
# # peptide_focus_heatmap = create_heatmap(
    # # raw_data_df[
        # # (raw_data_df['mz'] < 803)
        # # & (raw_data_df['rt'] > 2970)
        # # & (raw_data_df['rt'] < 3020)],
    # # type='tapp_smoothing',
    # # delta_rt=1,
    # # delta_mz=0.01,
    # # sigma_rt=5,
    # # sigma_mz=(0.0025 * 2, 200)
    # # )
# peptide_focus_heatmap = create_heatmap(
        # # raw_data_df[
            # # (raw_data_df['mz'] < 803)
            # # & (raw_data_df['rt'] > 2970)
            # # & (raw_data_df['rt'] < 3020)],
        # raw_data_big_df,
        # # raw_data_df,
        # type='tapp_smoothing',
        # # type='binning',
        # # transform='sqrt',
        # # transform='log',
        # delta_rt=1,
        # delta_mz=0.01,
        # sigma_rt=3,
        # sigma_mz=(0.0025, 200)
    # )
# # full_heatmap = create_heatmap(
        # # raw_data_df,
        # # type='tapp_smoothing',
        # # transform='none',
        # # delta_rt=1,
        # # delta_mz=0.01,
        # # sigma_rt=5,
        # # sigma_mz=(0.0025, 200)
    # # )
# peptide_focus_heatmap = create_heatmap(
        # # raw_data_df[
            # # (raw_data_df['mz'] < 803)
            # # & (raw_data_df['rt'] > 2970)
            # # & (raw_data_df['rt'] < 3020)],
        # # raw_data_big_df,
        # raw_data_tof_small_df,
        # type='tapp_smoothing',
        # # type='binning',
        # # transform='sqrt',
        # # transform='log',
        # delta_rt=3,
        # delta_mz=0.01,
        # sigma_rt=7,
        # sigma_mz=(0.0072/2, 200)
    # )
# local_max = pd.DataFrame(tapp.find_local_maxima(
    # peptide_focus_heatmap['img'].shape[1],
    # peptide_focus_heatmap['img'].shape[0],
    # peptide_focus_heatmap['img'].flatten()))
# local_max.columns = ['i', 'j', 'intensity']
# local_max['mz'] = local_max['i'] * 0.01 + np.min(peptide_focus_heatmap['bins_mz'])
# local_max['rt'] = local_max['j'] * 3 + np.min(peptide_focus_heatmap['bins_rt'])
# peaks = pd.DataFrame(tapp.find_peaks(
    # peptide_focus_heatmap['img'].shape[1],
    # peptide_focus_heatmap['img'].shape[0],
    # np.array(peptide_focus_heatmap['bins_mz']).min(),
    # np.array(peptide_focus_heatmap['bins_mz']).max(),
    # np.array(peptide_focus_heatmap['bins_rt']).min(),
    # np.array(peptide_focus_heatmap['bins_rt']).max(),
    # peptide_focus_heatmap['img'].flatten()))
# peaks.columns = ['i', 'j', 'mz', 'rt', 'height', 'total_intensity', 'sigma_mz', 'sigma_rt', 'border_background']
# peptide_focus_heatmap['img_plot'].scatter(local_max['i'], local_max['j'], s=3, alpha=0.7)

# Plot boundaries of the detected peaks based on the reported sigma.
# import matplotlib.patches as patches
# for i in range(0, peaks.shape[0]):
    # rect = patches.Rectangle(
            # (
                # peaks['i'][i] - peaks['sigma_mz'][i]/0.01/2,
                # peaks['j'][i] - peaks['sigma_rt'][i]/3/2
            # ),
            # peaks['sigma_mz'][i]/0.01,
            # peaks['sigma_rt'][i]/3,
            # linewidth=1,
            # edgecolor=np.random.rand(3,1).flatten(),
            # alpha=0.9,
            # facecolor='none'
        # )
    # peptide_focus_heatmap['img_plot'].add_patch(rect)

def estimate_sigma(raw_points):
        height_sum = 0
        x_sum = 0
        y_sum = 0
        x_sig = 0
        y_sig = 0
        mzs = np.array(raw_points['mz'])
        rts = np.array(raw_points['rt'])
        intensities = np.array(raw_points['intensity'])
        for row in range(0, raw_points.shape[0]):
            mz = mzs[row]
            rt = rts[row]

            height_sum += intensities[row]
            x_sum += intensities[row] * mz
            y_sum += intensities[row] * rt
            x_sig += intensities[row] * mz * mz
            y_sig += intensities[row] * rt * rt
        sigma_mz = math.sqrt((x_sig / height_sum) - (x_sum / height_sum) ** 2)
        sigma_rt = 0
        # sigma_rt = math.sqrt((y_sig / height_sum) - (y_sum / height_sum) ** 2)
        return({'sigma_mz': sigma_mz, 'sigma_rt': sigma_rt, 'total_intensity': height_sum})

detected_peaks_orbitrap = pd.read_csv("~/Projects/tapp-testbed/feature_finder_exploration/tapp_out_2/centroid/HD5YD15ED_DDA_R1.csv", sep=" ")
detected_peaks_orbitrap['mz'] = detected_peaks_orbitrap['X']
detected_peaks_orbitrap['rt'] = detected_peaks_orbitrap['Y']

def estimate_sigma_alg_1(peaks, sigma_mz, sigma_rt):
    target = pd.DataFrame(peaks)
    mzs = []
    rts = []
    estimated_sigmas_mz_fit = []
    estimated_sigmas_rt = []
    total_intensities = []
    for row in range(0, target.shape[0]):
        print("now processing row {0} of {1}".format(row, target.shape[0]))
        mz = target['mz'][row]
        rt = target['rt'][row]

        # Get all points within the tolerance range.
        min_mz = mz - sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        max_mz = mz + sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        min_rt = rt - sigma_rt
        max_rt = rt + sigma_rt
        raw_points_target = raw_data_big_df[
                (raw_data_big_df['mz'] > min_mz) &
                (raw_data_big_df['mz'] < max_mz) &
                (raw_data_big_df['rt'] > min_rt) &
                (raw_data_big_df['rt'] < max_rt)
            ]
        if raw_points_target.empty:
            continue

        estimation = estimate_sigma(raw_points_target)
        estimated_sigmas_mz = estimated_sigmas_mz + [estimation['sigma_mz']]
        estimated_sigmas_rt = estimated_sigmas_rt + [estimation['sigma_rt']]
        total_intensities = total_intensities + [estimation['total_intensity']]
        mzs = mzs + [mz]
        rts = rts + [rt]
    return({"mzs": mzs, "rts": rts, "sigmas_mz": estimated_sigmas_mz, "sigmas_rt": estimated_sigmas_rt})

# Using the central scan only for the sigma estimation.
def estimate_sigma_alg_2(peaks, sigma_mz, sigma_rt):
    target = pd.DataFrame(peaks)
    mzs = []
    rts = []
    estimated_sigmas_mz = []
    estimated_sigmas_rt = []
    total_intensities = []
    for row in range(0, target.shape[0]):
        print("now processing row {0} of {1}".format(row, target.shape[0]))
        mz = target['mz'][row]
        rt = target['rt'][row]

        # Get all points within the tolerance range.
        min_mz = mz - sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        max_mz = mz + sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        min_rt = rt - sigma_rt
        max_rt = rt + sigma_rt
        raw_points_target = raw_data_big_df[
                (raw_data_big_df['mz'] > min_mz) &
                (raw_data_big_df['mz'] < max_mz) &
                (raw_data_big_df['rt'] > min_rt) &
                (raw_data_big_df['rt'] < max_rt)
            ]
        raw_points_target = raw_points_target[
                    np.abs(raw_points_target['rt'] - rt) == np.abs(raw_points_target['rt'] - rt).min()
                ]
        if raw_points_target.empty:
            continue

        estimation = estimate_sigma(raw_points_target)
        estimated_sigmas_mz = estimated_sigmas_mz + [estimation['sigma_mz']]
        estimated_sigmas_rt = estimated_sigmas_rt + [estimation['sigma_rt']]
        total_intensities = total_intensities + [estimation['total_intensity']]
        mzs = mzs + [mz]
        rts = rts + [rt]
    return({
            "mzs": mzs,
            "rts": rts,
            "sigmas_mz": estimated_sigmas_mz,
            "sigmas_rt": estimated_sigmas_rt
        })

# Using the central scan only for the sigma estimation.
def estimate_sigma_alg_3(peaks, sigma_mz, sigma_rt):
    target = pd.DataFrame(peaks)
    mzs = []
    rts = []
    estimated_sigmas_mz_est = []
    estimated_sigmas_mz_fit = []
    estimated_sigmas_rt = []
    total_intensities = []
    for row in range(0, target.shape[0]):
        print("now processing row {0} of {1}".format(row, target.shape[0]))
        mz = target['mz'][row]
        rt = target['rt'][row]

        # Get all points within the tolerance range.
        min_mz = mz - sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        max_mz = mz + sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        min_rt = rt - sigma_rt
        max_rt = rt + sigma_rt
        raw_points_target = raw_data_big_df[
                (raw_data_big_df['mz'] > min_mz) &
                (raw_data_big_df['mz'] < max_mz) &
                (raw_data_big_df['rt'] > min_rt) &
                (raw_data_big_df['rt'] < max_rt)
            ]
        raw_points_target = raw_points_target[
                    np.abs(raw_points_target['rt'] - rt) == np.abs(raw_points_target['rt'] - rt).min()
                ]
        if raw_points_target.empty:
            continue

        try:
            # Gaussian fit of the central scan
            from scipy.optimize import curve_fit
            x = np.array(raw_points_target['mz'])
            y = np.array(raw_points_target['intensity'])

            n = len(x)                          #the number of data
            mean = np.sum(x * y)/ sum(y)
            sigma = np.sqrt(np.sum(y*(x-mean)**2)/sum(y))
            mean = sum(x * y) / sum(y)
            sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

            def gaus(x, a, x0, sigma):
                    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

            popt, pcov = curve_fit(gaus, x, y, p0=[max(y), mean, sigma])

            estimated_sigmas_mz_est = estimated_sigmas_mz_est + [sigma]
            estimated_sigmas_mz_fit = estimated_sigmas_mz_fit + [popt[2]]
            total_intensities = total_intensities + [popt[0]]
            mzs = mzs + [popt[1]]
            rts = rts + [rt]
        except:
            print("Error on row: {0} ".format(row))
    return({
            "mzs": mzs,
            "rts": rts,
            "sigmas_mz_est": estimated_sigmas_mz_est,
            "sigmas_mz_fit": estimated_sigmas_mz_fit,
            "sigmas_rt": estimated_sigmas_rt
        })

# Using the central scan and surrounding scans for the sigma estimation.
def estimate_sigma_alg_4(peaks, sigma_mz, sigma_rt):
    target = pd.DataFrame(peaks)
    mzs = []
    rts = []
    estimated_sigmas_mz_est = []
    estimated_sigmas_mz_fit = []
    estimated_sigmas_rt = []
    total_intensities = []
    for row in range(0, target.shape[0]):
        print("now processing row {0} of {1}".format(row, target.shape[0]))
        mz = target['mz'][row]
        rt = target['rt'][row]

        # Get all points within the tolerance range.
        min_mz = mz - sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        max_mz = mz + sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])
        min_rt = rt - sigma_rt
        max_rt = rt + sigma_rt
        raw_points_target = raw_data_big_df[
                (raw_data_big_df['mz'] > min_mz) &
                (raw_data_big_df['mz'] < max_mz) &
                (raw_data_big_df['rt'] > min_rt) &
                (raw_data_big_df['rt'] < max_rt)
            ]
        if raw_points_target.empty:
            continue

        try:
            # Gaussian fit of the central scan
            x = np.array(raw_points_target['mz'])
            y = np.array(raw_points_target['intensity'])

            n = len(x)                          #the number of data
            mean = np.sum(x * y)/ sum(y)
            sigma = np.sqrt(np.sum(y*(x-mean)**2)/sum(y))
            mean = sum(x * y) / sum(y)
            sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

            def gaus(x, a, x0, sigma):
                    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

            popt, pcov = curve_fit(gaus, x, y, p0=[max(y), mean, sigma])

            estimated_sigmas_mz_est = estimated_sigmas_mz_est + [sigma]
            estimated_sigmas_mz_fit = estimated_sigmas_mz_fit + [popt[2]]
            total_intensities = total_intensities + [popt[0]]
            mzs = mzs + [popt[1]]
            rts = rts + [rt]
        except:
            print("Error on row: {0} ".format(row))
    return({
            "mzs": mzs,
            "rts": rts,
            "sigmas_mz_est": estimated_sigmas_mz_est,
            "sigmas_mz_fit": estimated_sigmas_mz_fit,
            "sigmas_rt": estimated_sigmas_rt
        })

# # Exploring a single scan and a small number of peaks to see what is the
# # estimated sigma and magnitude and explore gaussian fitting.
# scan = raw_data_big_df[raw_data_big_df["rt"] == 2554.95]
# scan_peaks = pd.DataFrame(scan[(scan['mz'] > 510.24) & (scan['mz'] < 510.30)])
# scan_peaks['norm_intensity'] = scan_peaks['intensity']/scan_peaks['intensity'].max() 
# plt.cla(); plt.plot(scan_peaks['mz'], scan_peaks['norm_intensity']) 

# # First derivative.
# dydx = np.diff(scan_peaks['intensity'])/np.diff(scan_peaks['mz'])

# sigma_estimation_1 = estimate_sigma_alg_1(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025, 200), 3)
# sigma_estimation_2 = estimate_sigma_alg_2(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025, 200), 3)
# sigma_estimation_4 = estimate_sigma_alg_4(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025, 200), 3)

# Plot the results of the different estimations.
# plt.clf()
# plt.subplot(3,1,1)
# plt.scatter(sigma_estimation_1['mzs'], sigma_estimation_1['sigmas_mz'], s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_1 [sigma_mz: (0.0025, 200), sigma_rt: 3]")
# plt.subplot(3,1,2)
# plt.scatter(sigma_estimation_2['mzs'], sigma_estimation_2['sigmas_mz'], s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_2 [sigma_mz: (0.0025, 200), sigma_rt: 3]")
# plt.subplot(3,1,3)
# plt.scatter(np.array(sigma_estimation_3['mzs'])[np.array(sigma_estimation_3['sigmas_mz_fit']) < 0.03], np.array(sigma_estimation_3['sigmas_mz_fit'])[np.array(sigma_estimation_3['sigmas_mz_fit']) < 0.03], s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_3 [sigma_mz: (0.0025, 200), sigma_rt: 3]")

# sigma_estimation_3 = estimate_sigma_alg_3(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025/2, 200), 3)
# plt.clf()
# plt.subplot(2,1,1)
# plt.scatter(np.array(sigma_estimation_3['mzs']), np.array(sigma_estimation_3['sigmas_mz_fit']), s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_3 (fitted sigma) [sigma_mz: (0.0025/2, 200), sigma_rt: 3]")
# plt.subplot(2,1,2)
# plt.scatter(np.array(sigma_estimation_3['mzs']), np.array(sigma_estimation_3['sigmas_mz_est']), s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_3 (estimated sigma) [sigma_mz: (0.0025/2, 200), sigma_rt: 3]")
# fig.savefig("alg_3_sigma_mz_estimation_fitted_0p00125at200_3sec.png", dpi=150)
# sigma_estimation_3 = estimate_sigma_alg_3(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025, 200), 3)
# sigma_estimation_3 = estimate_sigma_alg_3(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025 * 3/2, 200), 3)
# sigma_estimation_3 = estimate_sigma_alg_3(detected_peaks_orbitrap.iloc[0:10000,:], (0.0025 * 2, 200), 3)
# outliers_index = [i for i in range(0, len(sigma_estimation_3['sigmas_mz_fit'])) if sigma_estimation_3['sigmas_mz_fit'][i] > 0.02]
# not_outliers_index = [i for i in range(0, len(sigma_estimation_3['sigmas_mz_fit'])) if sigma_estimation_3['sigmas_mz_fit'][i] <= 0.02]
# plt.clf()
# plt.subplot(2,1,1)
# plt.scatter(np.array(sigma_estimation_3['mzs']), np.array(sigma_estimation_3['sigmas_mz_fit']), s=4, alpha=0.4, color='crimson')
# # plt.scatter(np.array(sigma_estimation_3['mzs'])[not_outliers_index], np.array(sigma_estimation_3['sigmas_mz_fit'])[not_outliers_index], s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_3 (fitted sigma) [sigma_mz: (0.0025 * 3/2, 200), sigma_rt: 3]")
# plt.subplot(2,1,2)
# plt.scatter(np.array(sigma_estimation_3['mzs']), np.array(sigma_estimation_3['sigmas_mz_est']), s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_3 (estimated sigma) [sigma_mz: (0.0025 * 3/2, 200), sigma_rt: 3]")
# # fig.savefig("alg_3_sigma_mz_estimation_fitted_0p0025at200_3sec.png", dpi=150)
# # row 5650, 5699?

# # sigma_estimation_4 = estimate_sigma_alg_4(detected_peaks_orbitrap.iloc[0:1000,:], (0.0025 * 3/2, 200), 5)
# outliers_index_4 = [i for i in range(0, len(sigma_estimation_4['sigmas_mz_fit'])) if sigma_estimation_4['sigmas_mz_fit'][i] > 0.02]
# not_outliers_index_4 = [i for i in range(0, len(sigma_estimation_4['sigmas_mz_fit'])) if sigma_estimation_4['sigmas_mz_fit'][i] <= 0.02]
# plt.clf()
# plt.subplot(2,1,1)
# # plt.scatter(np.array(sigma_estimation_4['mzs']), np.array(sigma_estimation_4['sigmas_mz_fit']), s=4, alpha=0.4, color='crimson')
# plt.scatter(np.array(sigma_estimation_4['mzs'])[not_outliers_index_4], np.array(sigma_estimation_4['sigmas_mz_fit'])[not_outliers_index_4], s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_4 (fitted sigma) [sigma_mz: (0.0025, 200), sigma_rt: 5]")
# plt.subplot(2,1,2)
# plt.scatter(np.array(sigma_estimation_4['mzs']), np.array(sigma_estimation_4['sigmas_mz_est']), s=4, alpha=0.4, color='crimson')
# plt.title("algorithm_4 (estimated sigma) [sigma_mz: (0.0025, 200), sigma_rt: 5]")

def resample(
        df,
        delta_rt=1.5,
        delta_mz=0.01,
        method='binning',
        sigma_mz=(0.0025, 200),
        sigma_rt=15,
    ):
    print("Resampling grid with parameters:\n> delta_mz = {0}\n> delta_rt = {1}\n> sigma_mz = {2}\n> sigma_rt = {3}\n> method = {4}".format(delta_mz, delta_rt, sigma_mz, sigma_rt, method))
    unique_rt = df.sort_values("rt")["rt"].unique()
    min_rt = unique_rt.min()
    max_rt = unique_rt.max()
    num_bins_rt = math.floor((max_rt - min_rt)/ delta_rt) + 1
    bins_rt =  [min_rt + delta_rt * x for x in range(0, num_bins_rt)] 

    unique_mz = df.sort_values("mz")["mz"].unique()
    min_mz = unique_mz.min()
    max_mz = unique_mz.max()
    num_bins_mz = math.floor((max_mz - min_mz)/ delta_mz) + 1
    bins_mz =  [min_mz + delta_mz * x for x in range(0, num_bins_mz)] 
    grid_shape = (num_bins_rt, num_bins_mz)
    print("> min_mz = {0}".format(min_mz))
    print("> max_mz = {0}".format(max_mz))
    print("> min_rt = {0}".format(min_rt))
    print("> max_rt = {0}".format(max_rt))
    print("> grid_shape (rt,mz) = {0}".format(grid_shape))
    print("> grid_size (Bytes) = {0}".format(grid_shape[0] * grid_shape[1] * 8))

    grid = np.zeros(num_bins_rt * num_bins_mz)
    grid = np.reshape(grid, grid_shape)
    
    if method == 'binning':
        # Splat points in matrix
        for row in range(0, df.shape[0]):
            # Find y index
            rt = df.iloc[row, 1]
            x = math.floor((rt - min_rt) / delta_rt)
            # Find x index
            mz = df.iloc[row, 0]
            y = math.floor((mz - min_mz) / delta_mz)
            # Update matrix value
            grid[x,y] = grid[x,y] + df.iloc[row, 2]
    elif method == 'tapp_smoothing':
        # Splat points in matrix
        for row in range(0, df.shape[0]):
            rt = df.iloc[row, 1]
            mz = df.iloc[row, 0]
            intensity = df.iloc[row, 2]
            current_sigma_mz = sigma_at_mz(mz, sigma_mz[0], sigma_mz[1])

            kernel_min_rt = rt - 2 * sigma_rt
            kernel_max_rt = rt + 2 * sigma_rt
            kernel_min_mz = mz - 2 * current_sigma_mz
            kernel_max_mz = mz + 2 * current_sigma_mz

            i_min = 0 if kernel_min_mz < min_mz else math.floor((kernel_min_mz - min_mz) / delta_mz)
            j_min = 0 if kernel_min_rt < min_rt else math.floor((kernel_min_rt - min_rt) / delta_rt)
            i_max = num_bins_mz - 1 if kernel_max_mz > max_mz else math.floor((kernel_max_mz - min_mz) / delta_mz)
            j_max = num_bins_rt - 1 if kernel_max_rt > max_rt else math.floor((kernel_max_rt - min_rt) / delta_rt)

            for j in range(j_min, j_max + 1):
                for i in range(i_min, i_max + 1):
                    x = bins_mz[i]
                    y = bins_rt[j]

                    a = (x - mz) / current_sigma_mz
                    b = (y - rt) / sigma_rt
                    weight = math.exp(-0.5 * (a * a + b * b))

                    grid[j, i] = grid[j, i] + intensity * weight
    # elif method == 'tapp_internal_smoothing':

    return({
        "grid": grid,
        "bins_mz": bins_mz,
        "bins_rt": bins_rt,
        "delta_mz": delta_mz,
        "delta_rt": delta_rt,
    })

def plot_grid(
        figure,
        grid_obj,
        peaks=None,
        transform='none',
        ):
    grid = grid_obj["grid"]
    bins_rt = grid_obj["bins_rt"]
    bins_mz = grid_obj["bins_mz"]
    num_bins_mz = len(bins_mz)
    num_bins_rt = len(bins_rt)
    delta_mz = grid_obj["delta_mz"]
    delta_rt = grid_obj["delta_rt"]
    min_mz = np.array(bins_mz).min()
    max_mz = np.array(bins_mz).max()
    min_rt = np.array(bins_rt).min()
    max_rt = np.array(bins_rt).max()
    if transform == 'sqrt':
        grid = np.sqrt(grid)
    elif transform == 'log':
        grid = np.log(grid + 0.0001)

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
    grid_plot.imshow(grid, aspect='auto', origin="lower")

    locs_x = grid_plot.get_xticks()
    labels_x = []
    for loc in locs_x:
        if loc < 0 or loc >= num_bins_mz:
            labels_x = labels_x + [""]
        else:
            labels_x = labels_x + ["{0:.2f}".format(min_mz + delta_mz * loc)]

    locs_y = grid_plot.get_yticks()
    labels_y = []
    for loc in locs_y:
        if loc < 0 or loc >= num_bins_rt:
            labels_y = labels_y + [""]
        else:
            labels_y = labels_y + ["{0:.2f}".format(min_rt + delta_rt * loc)]
    
    grid_plot.set_xticklabels(labels_x)
    grid_plot.set_yticklabels(labels_y)
    grid_plot.set_xlabel("m/z")
    grid_plot.set_ylabel("retention time (s)")

    fig.suptitle("min_rt: {0:.3f}, max_rt: {1:.3f}, min_mz: {2:.3f}, max_mz: {3:.3f}, delta_rt: {4:.3f}, delta_mz: {5:.3f}".format(
        min_rt, max_rt, min_mz, max_mz, delta_rt, delta_mz))

    # Display the peak boundary.
    # FIXME: Not working.
    if False:
        # Plot boundaries of the detected peaks based on the reported sigma.
        import matplotlib.patches as patches
        for i in range(0, peaks.shape[0]):
            rect = patches.Rectangle(
                    (
                        peaks['i'][i] - peaks['sigma_mz'][i] * 100,
                        peaks['j'][i] - peaks['sigma_rt'][i]
                    ),
                    peaks['sigma_mz'][i] * 100 * 2,
                    peaks['sigma_rt'][i],
                    linewidth=1,
                    edgecolor=np.random.rand(3,1).flatten(),
                    alpha=0.9,
                    facecolor='none'
                )
            grid_plot.add_patch(rect)
    # Display the peak centers.
    if True:
        grid_plot.scatter(peaks['i'], peaks['j'], color='aqua', s=3)

    return({
        "grid_plot": grid_plot,
        "mz_plot": mz_plot,
        "rt_plot": rt_plot,
    })

# grid = resample(
    # raw_data_big_df[
        # (raw_data_big_df['mz'] > 801)
        # & (raw_data_big_df['mz'] < 803)
        # & (raw_data_big_df['rt'] > 2808)
        # & (raw_data_big_df['rt'] < 2928)],
    # delta_mz=0.002445,
    # delta_rt=2.9/3,
    # sigma_mz=(0.0025/4, 200),
    # sigma_rt=3,
    # method='tapp_smoothing'
    # )

# print("Finding peaks.")
# peaks = pd.DataFrame(tapp.find_peaks(
    # grid['grid'].shape[1],
    # grid['grid'].shape[0],
    # np.array(grid['bins_mz']).min(),
    # np.array(grid['bins_mz']).max(),
    # np.array(grid['bins_rt']).min(),
    # np.array(grid['bins_rt']).max(),
    # grid['grid'].flatten()))
# peaks.columns = ['i', 'j', 'mz', 'rt', 'height', 'total_intensity', 'sigma_mz', 'sigma_rt', 'border_background']

# # Find the sigma

# print("Plotting.")
grid_plot = plot_grid(
    fig,
    grid,
    peaks,
    transform='sqrt'
    )

# grid_plot['grid_plot'].clear()
# Plot the raw data for the first detected peak.
unique_rt = raw_data_big_df['rt'].unique() 
plt.figure(fig_2.number)
plt.clf()
plt.figure(fig_3.number)
plt.clf()
fitted_peaks = []
# for i in range(0, 100):
# for i in [19]:
# for i in range(0, 5):
# for i in range(100, 105):
for i in range(0, peaks.shape[0]):
# for i in range(57, 58):
# for i in range(85,86):
    selected_peak = peaks.iloc[i]
    theoretical_sigma_mz = sigma_at_mz(selected_peak['mz'], 0.0028/2.335, 200)
    theoretical_sigma_rt = 6
    tolerance_mz = 3 * theoretical_sigma_mz

    # Find the closest scan.
    unique_rt_idx = np.abs(unique_rt - selected_peak['rt']).argmin()
    closest_scan = raw_data_big_df[
            (raw_data_big_df['rt'] >= selected_peak['rt'] - theoretical_sigma_rt) &
            (raw_data_big_df['rt'] <= selected_peak['rt'] + theoretical_sigma_rt)
        ]

    # Select the range to display.
    selected_range = closest_scan[
            (closest_scan['mz'] >= selected_peak['mz'] - tolerance_mz) &
            (closest_scan['mz'] <= selected_peak['mz'] + tolerance_mz)
        ]

    # NOTE: This is a performance improvement parameter. Could be set to 0 to not
    # filter by minimum number of scans per peak.
    if len(selected_range['rt'].unique()) < 3:
        continue

    selected_range = selected_range.sort_values(['mz'])
    selected_range = selected_range.groupby('mz').apply(lambda x: pd.Series({'mz':x['mz'].mean(), 'rt': x['rt'].mean(), 'intensity': x['intensity'].mean()}))
    intensity = selected_range['intensity']
    norm_intensity = intensity/intensity.max()
    mzs = selected_range['mz']

    # Fit a gaussian on the selected peak.
    x = np.array(mzs)
    y = np.array(intensity)

    mean = selected_peak['mz']
    sigma = theoretical_sigma_mz

    def gaus(x, a, x_0, sigma_x):
        return a * np.exp(-0.5 * ((x - x_0)/sigma_x) ** 2)

    def gaus2d(X, a, x_0, sigma_x, y_0, sigma_y):
        x = X[0]
        y = X[1]
        return a * np.exp(-0.5 * ((x - x_0) /sigma_x) ** 2 ) * np.exp(-0.5 * ((y - y_0)/sigma_y) ** 2)
    try:
        popt, pcov = curve_fit(gaus, x, y, p0=[max(y), mean, sigma])
        fitted_intensity = gaus(x, *popt)
        norm_fitted_intensity = fitted_intensity/fitted_intensity.max()
    except:
        print("error when fitting 1D gaussian on peak: {}".format(i))
        continue

    try:
        X = np.array([selected_range['mz'], selected_range['rt']])
        mean_x = sum(X[0] * y) / sum(y)
        sigma_x = np.sqrt(sum(y * (X[0] - mean_x)**2) / sum(y))
        mean_y = sum(X[1] * y) / sum(y)
        sigma_y = np.sqrt(sum(y * (X[1] - mean_y)**2) / sum(y))
        popt_2d, pcov_2d = curve_fit(gaus2d, X, y, p0=[max(y), mean_x, sigma_x, mean_y, sigma_y])

    except:
        print("error when fitting 2D gaussian on peak: {}".format(i))
        continue

    # # Discard peaks where the fit is not conforming with the theoretical distributions.
    if popt_2d[2] == 0 or popt_2d[4] == 0:
        continue

    if popt_2d[2] > 5 * theoretical_sigma_mz or popt_2d[4] > 5 * theoretical_sigma_rt:
        continue

    # Save the fitted parameters.
    fitted_peaks = fitted_peaks + [{
            'smoothed_mz': selected_peak['mz'],
            'smoothed_rt': selected_peak['rt'],
            'smoothed_height': selected_peak['height'],
            'smoothed_sigma_mz': selected_peak['sigma_mz'],
            'smoothed_sigma_rt': selected_peak['sigma_rt'],
            'smoothed_total_intensity': selected_peak['total_intensity'],
            'fitted_mz': popt_2d[1],
            'fitted_rt': popt_2d[3],
            'fitted_height': popt_2d[0],
            'fitted_sigma_mz': popt_2d[2],
            'fitted_sigma_rt': popt_2d[4],
            'fitted_total_intensity': selected_range['intensity'].sum(),
            'roi_min_mz': selected_range['mz'].min(),
            'roi_max_mz': selected_range['mz'].max(),
            'roi_min_rt': selected_range['rt'].min(),
            'roi_max_rt': selected_range['rt'].max(),
        }]

    color = np.random.rand(3,1).flatten()
    # MZ fit plot.
    fitted_intensity_2d_mz = gaus(selected_range['mz'], popt_2d[0],popt_2d[1],popt_2d[2])
    norm_fitted_intensity_2d_mz = fitted_intensity_2d_mz/fitted_intensity_2d_mz.max()
    plt.figure(fig_2.number)
    plt.plot(selected_range['mz'], norm_intensity, linestyle='-', color=color, alpha=0.5, label='norm_intensity')
    plt.plot(selected_range['mz'], norm_fitted_intensity, linestyle=':', color=color, label='1d_fitting')
    # plt.scatter(selected_peak['mz'] - tolerance_mz, 0, color="blue")
    # plt.scatter(selected_peak['mz'] + tolerance_mz, 0, color="red")
    plt.plot(selected_range['mz'], norm_fitted_intensity_2d_mz, linestyle='--', color=color, label='2d_fitting')
    # ax = plt.axes()
    # lines = ax.get_lines()
    # legend = plt.legend([lines[i] for i in [0,3]], ["a"], loc=1)
    # axes.add_artist(legend)

    # RT fit plot.
    selected_range = selected_range.sort_values(['rt'])
    fitted_intensity_2d_rt = gaus(selected_range['rt'], popt_2d[0],popt_2d[3],popt_2d[4])
    norm_fitted_intensity_2d_rt = fitted_intensity_2d_rt/fitted_intensity_2d_rt.max()
    plt.figure(fig_3.number)
    plt.plot(selected_range['rt'], norm_fitted_intensity_2d_rt, color=color, linestyle='--')
    tics = selected_range.groupby('rt').apply(lambda x: x['intensity'].sum())
    plt.plot(tics/tics.max(), label=str(i), linestyle='-', color=color, alpha=0.5)

    # Contour plot
    x = grid['bins_mz']
    y = grid['bins_rt']
    x, y = np.meshgrid(x, y)
    Z = gaus2d((x,y), *popt_2d) 
    x = range(0, len(grid['bins_mz']))
    y = range(0, len(grid['bins_rt']))
    x, y = np.meshgrid(x, y)
    grid_plot['grid_plot'].contour(x, y, Z, levels=[1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10], colors=[color], alpha=0.5)
    grid_plot['grid_plot'].scatter(selected_peak['i'], selected_peak['j'], color=color, marker='x')
    # grid_plot['grid_plot'].scatter(popt_2d[1]/(np.array(grid['bins_mz']).max() - np.array(grid['bins_mz']).min()), popt_2d[3] - np.array(grid['bins_rt']).min(), color=color, marker='P')
    # print(np.abs(popt_2d[1] - selected_peak['mz']))

fitted_peaks = pd.DataFrame(fitted_peaks)
