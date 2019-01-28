import tapp
import math
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

# raw_data = tapp.raw_data("/data/Spike-In QEx UKE/HD5YD15ED_DDA_R1.mzXML", 799.036, 809.107, 2768.83, 3171.63)
raw_data = tapp.raw_data_load_dump("HD5YD15ED_DDA_R1_small.rawdata")
# raw_data.dump("HD5YD15ED_DDA_R1_small.rawdata")

raw_data_df = pd.DataFrame(
    {
        "mz": raw_data.mz(),
        "rt": raw_data.rt(),
        "intensity": raw_data.intensity()
    })

def create_heatmap(df, delta_rt=1.5, delta_mz=0.01, transform='none'):
    unique_rt = df.sort_values("rt")["rt"].unique()
    min_rt = unique_rt.min()
    max_rt = unique_rt.max()
    # min_rt_diff = math.inf
    # for i in range(1, len(unique_rt)):
        # rt_diff = unique_rt[i] - unique_rt[i-1]
        # if rt_diff < min_rt_diff:
            # min_rt_diff = rt_diff
    # delta_rt = 1.4
    # delta_rt = 0.07
    num_bins_rt = math.floor((max_rt - min_rt)/ delta_rt) + 1
    bins_rt =  [min_rt + delta_rt * x for x in range(0, num_bins_rt)] 

    unique_mz = df.sort_values("mz")["mz"].unique()
    min_mz = unique_mz.min()
    max_mz = unique_mz.max()
    # min_mz_diff = math.inf
    # for i in range(1, len(unique_mz)):
        # rt_diff = unique_mz[i] - unique_mz[i-1]
        # if rt_diff < min_mz_diff:
            # min_mz_diff = rt_diff
    # delta_mz = 0.007
    # delta_mz = min_mz_diff / 2
    num_bins_mz = math.floor((max_mz - min_mz)/ delta_mz) + 1
    bins_mz =  [min_mz + delta_mz * x for x in range(0, num_bins_mz)] 

    # print("min_rt_diff: {0}, delta_rt: {1}, num_bins_rt: {2}, min_rt: {3}, max_rt: {4}".format(min_rt_diff, delta_rt, num_bins_rt, min_rt, max_rt))
    # print("min_mz_diff: {0}, delta_mz: {1}, num_bins_mz: {2}".format(min_mz_diff, delta_mz, num_bins_mz))

    # img = np.random.random(num_bins_rt * num_bins_mz)
    img = np.zeros(num_bins_rt * num_bins_mz)

    img = np.reshape(img, (num_bins_rt, num_bins_mz))
    # img = pd.DataFrame(img)
    
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

heatmap = create_heatmap(raw_data_df, transform='sqrt')
peptide_focus_heatmap = create_heatmap(
    raw_data_df[
        (raw_data_df['mz'] < 803)
        & (raw_data_df['rt'] > 2970)
        & (raw_data_df['rt'] < 3020)]
    )
# peak_focus_heatmap = create_heatmap(
    # raw_data_df[
        # (raw_data_df['mz'] > 799.54)
        # & (raw_data_df['mz'] < 800.04)
        # & (raw_data_df['rt'] > 2970)
        # & (raw_data_df['rt'] < 3020)]
    # )
