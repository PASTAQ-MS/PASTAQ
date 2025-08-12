# helper functions for PASTAQ demonstration
import math as math
import numpy as np
import os
import requests
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pastaq as pq
import plotly.graph_objects as go
import plotly.tools as tls

# helper functions
# Function to calculate the Euclidean distance
def euclidean_distance(peak1, peak2):
    peak1x = peak1.fitted_rt
    peak1y = peak1.fitted_mz
    peak2x = peak2.fitted_rt + peak2.rt_delta
    peak2y = peak2.fitted_mz
    return math.sqrt(((peak1x - peak2x)/(peak1.fitted_sigma_rt + peak2.fitted_sigma_rt)) ** 2 + ((peak1y - peak2y)/(peak1.fitted_sigma_mz + peak2.fitted_sigma_mz)) ** 2)

# Function to find the closest peak
def find_closest_peak(reference_peak, peaks, mz_tolerance = 1, rt_tolerance = 1):
    closest_peak = None
    min_distance = float('inf')
    
    for peak in peaks:
        distance = euclidean_distance(reference_peak, peak)
        if distance < min_distance:
            min_distance = distance
            closest_peak = peak
    
    if reference_peak.fitted_mz - mz_tolerance <= closest_peak.fitted_mz <= reference_peak.fitted_mz + mz_tolerance and reference_peak.fitted_rt - rt_tolerance <= closest_peak.fitted_rt + closest_peak.rt_delta <= reference_peak.fitted_rt + rt_tolerance:
        return closest_peak
    else:
        print("Reference peak not found in the peaks list.")
        return None

# get segment indices for a value
def find_segments_indices(values, starts, ends):
    """
    Finds the segment index for each value in the list.
    
    Args:
        values: A list of values to find the segments for.
        starts: A list of segment start values.
        ends: A list of segment end values.
    
    Returns:
        A list of indices corresponding to the segment each value lies in, or -1 if a value is not in any segment.
    """
    # Ensure the lists are of the same length
    if len(starts) != len(ends):
        raise ValueError("The 'starts' and 'ends' lists must have the same length.")
    
    indices = []
    
    for value in values:
        found = False
        for i in range(len(starts)):
            if starts[i] <= value <= ends[i]:
                indices.append(i)
                found = True
                break
        if not found:
            indices.append(-1)
    
    return indices

# interpolate retention time to warped retention time using aligned segments limits
def interpolate_values(values, starts, ends, new_starts, new_ends):
    """
    Interpolates the input values based on corresponding segment index values.
    
    Args:
        values: A list of values to interpolate.
        starts: A list of segment start values.
        ends: A list of segment end values.
        new_starts: A list of new segment start values for interpolation.
        new_ends: A list of new segment end values for interpolation.
    
    Returns:
        A list of interpolated values.
    """
    indices = find_segments_indices(values, starts, ends)
    interpolated_values = []
    
    for value, index in zip(values, indices):
        if index == -1:
            interpolated_values.append(None)  # Value is not in any segment
        else:
            # Perform linear interpolation
            old_start = starts[index]
            old_end = ends[index]
            new_start = new_starts[index]
            new_end = new_ends[index]
            
            # Linear interpolation formula
            if old_end != old_start:  # Avoid division by zero
                interpolated_value = new_start + ((value - old_start) / (old_end - old_start)) * (new_end - new_start)
            else:
                interpolated_value = new_start  # If old_start == old_end, use new_start
            
            interpolated_values.append(interpolated_value)
    
    return interpolated_values

# calculate Gaussian peak with mean, sigma, and height
def gaussian(x, mean, sigma, height):
    return height * np.exp(-((x - mean) ** 2) / (2 * sigma ** 2))

# plot Gaussian peak
def plot_gaussian(mean, sigma, height):
    x = np.linspace(mean - 4 * sigma, mean + 4 * sigma, 1000)
    y = gaussian(x, mean, sigma, height)
    plt.plot(x, y, label='fitted Gaussian peak', color='orange', alpha=0.75, linestyle='-', linewidth=1)

# plot mass spectra considering that only non-eros values are included in the spectra
def plot_msSpectra(mz, intensity, norm_mz_diff = 0.0035, diffFactor = 1.3, scanIdx = 0):
    """
    Plotting mass spectra, expecting spectra where 0 intensity values were omitted.
    
    Args:
        mz: A list of mz values of the mass spectra.
        intensity: A list of intensity (non-zero) values of the mass spectra.
        norm_mz_diff: Difference between two adjacent mz at measurement point corresponding to the original sampling frequency of the mass spectra.
        diffFactor: Tolerance factor allowing to vary sampling frequency of mass spectra. It is typically set to 30% (factor 1.3).
    
    Returns:
        Figure object as fig.
    """
    spectra = {'mz': mz, 'intensity': intensity}
    newSpectra = {'mz': [], 'intensity': []}
    idxSpectra = 0
    for i in range(1, len(spectra['mz'])):
        diff = spectra['mz'][i]-spectra['mz'][i-1]
        if diff > norm_mz_diff*diffFactor:
            if diff < norm_mz_diff*2:
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i-1] + diff/2)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                print(norm_mz_diff*diffFactor, diff, norm_mz_diff*diffFactor + diff, norm_mz_diff)
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
                newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
                idxSpectra += 1
            else:
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i-1] + norm_mz_diff)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i] - norm_mz_diff)
                newSpectra['intensity'].insert(idxSpectra, 0)
                idxSpectra += 1
                newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
                newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
                idxSpectra += 1
        else:
            newSpectra['mz'].insert(idxSpectra, spectra['mz'][i])
            newSpectra['intensity'].insert(idxSpectra, spectra['intensity'][i])
            if (diff/diffFactor) > norm_mz_diff:
                norm_mz_diff = diff
            idxSpectra += 1

    fig = plt.figure(figsize=(25, 6), facecolor='white')  # Set the figure size and white background
    plt.plot(newSpectra['mz'], newSpectra['intensity'], color = 'red', marker='', linestyle='-')  # Plot mz vs. intensity
    plt.xlabel('m/z')  # Set the x-axis label
    plt.ylabel('Intensity')  # Set the y-axis label
    plt.title('Mass Spectrum of scan {}' .format(scanIdx))  # Set the title
    plt.grid(False)  # Show grid
    return fig

def plot_meshPeaks(mesh, peaks, localMax=[], size=(15,10), withFailedPeaks=False, showIndeces=False):
    """
    Plotting mesh grid with peak positions.
    
    Args:
        mesh: smoothed grid object.
        peaks: A list of peaks values.
        localMax: A list of local maxima values.
        size: A tuple defining the size of the figure (width, height).
        withFailedPeaks: If True, failed peaks are also plotted.
        showIndeces: If True, peak indeces are shown in the plot.
        
    Returns:
        Figure object as fig.
    """
    mzVec, rtVec, mzVecMax, rtVecMax, mzLocMax, rtLocMax = [], [], [], [], [], []
    # Create a mesh plot
    plot = pq.plot_mesh(mesh, transform='sqrt', figure=None)
    # print(size)
    plot['img_plot'].get_figure().set_size_inches(size[0], size[1], forward=True)
    
    # Get mesh bin ranges for filtering
    mz_min, mz_max = min(mesh.bins_mz), max(mesh.bins_mz)
    rt_min, rt_max = min(mesh.bins_rt), max(mesh.bins_rt)
    
    if not withFailedPeaks:
        for k in range(len(peaks)):
            # only include peaks with fit_failure_code == 0
            if peaks[k].fit_failure_code == 0:
                # Only include peaks within mesh bin ranges
                if (mz_min <= peaks[k].local_max_mz <= mz_max and rt_min <= peaks[k].local_max_rt <= rt_max):
                    mzVecMax.append(peaks[k].local_max_mz)
                    rtVecMax.append(peaks[k].local_max_rt)
                if (mz_min <= peaks[k].fitted_mz <= mz_max and rt_min <= peaks[k].fitted_rt <= rt_max):
                    mzVec.append(peaks[k].fitted_mz)
                    rtVec.append(peaks[k].fitted_rt)
        plot['img_plot'].scatter(mzVecMax, rtVecMax, s=100, facecolors='none', edgecolors='green', marker='o', alpha=0.7) # shows local maxima in quantified peaks
        plot['img_plot'].scatter(mzVec, rtVec, s=100, c='red', marker='.') # shows centre of quantified peaks
        if localMax != []:
            for k in range(len(localMax)):
                # Only include local maxima within mesh bin ranges
                if (mz_min <= localMax[k].mz <= mz_max and rt_min <= localMax[k].rt <= rt_max):
                    mzLocMax.append(localMax[k].mz)
                    rtLocMax.append(localMax[k].rt)
            plot['img_plot'].scatter(mzLocMax, rtLocMax, s=50, c='blue', marker='x', alpha=0.7) # shows local maxima of non quantified peaks
    else:
        mzVecnF, rtVecnF, mzVecMaxnF, rtVecMaxnF = [], [], [], []
        for k in range(len(peaks)):
            if peaks[k].fit_failure_code == 0:
                # Only include quantified peaks within mesh bin ranges
                if (mz_min <= peaks[k].fitted_mz <= mz_max and rt_min <= peaks[k].fitted_rt <= rt_max):
                    mzVec.append(peaks[k].fitted_mz)
                    rtVec.append(peaks[k].fitted_rt)
                if (mz_min <= peaks[k].local_max_mz <= mz_max and rt_min <= peaks[k].local_max_rt <= rt_max):
                    mzVecMax.append(peaks[k].local_max_mz)
                    rtVecMax.append(peaks[k].local_max_rt)
            else:
                # Only include non-quantified peaks within mesh bin ranges
                if (mz_min <= peaks[k].fitted_mz <= mz_max and rt_min <= peaks[k].fitted_rt <= rt_max):
                    mzVecnF.append(peaks[k].fitted_mz)
                    rtVecnF.append(peaks[k].fitted_rt)
                if (mz_min <= peaks[k].local_max_mz <= mz_max and rt_min <= peaks[k].local_max_rt <= rt_max):
                    mzVecMaxnF.append(peaks[k].local_max_mz)
                    rtVecMaxnF.append(peaks[k].local_max_rt)
        plot['img_plot'].scatter(mzVecMax, rtVecMax, s=100, facecolors='none', edgecolors='green', marker='o', alpha=0.7) # shows local maxima in quantified peaks
        plot['img_plot'].scatter(mzVec, rtVec, s=100, c='red', marker='.') # shows centre of quantified peaks
        plot['img_plot'].scatter(mzVecnF, rtVecnF, s=100, c='purple', marker='.') # shows centre of non quantified peaks
        plot['img_plot'].scatter(mzVecMaxnF, rtVecMaxnF, s=100,  facecolors='none', edgecolors='yellow', marker='o', alpha=0.7) # shows local maxima in non quantified peaks
    
    if showIndeces:
        for k in range(len(peaks)):
            # Only show indices for peaks within mesh bin ranges
            if not withFailedPeaks and peaks[k].fit_failure_code != 0:
                continue
            if (mz_min <= peaks[k].fitted_mz <= mz_max and rt_min <= peaks[k].fitted_rt <= rt_max):
                plot['img_plot'].text(
                    peaks[k].fitted_mz + 0.01,  # adjust offset as needed
                    peaks[k].fitted_rt - 0.01,  # adjust offset as needed
                    str(k),
                    fontsize=10,
                    color='white',
                    ha='right',
                    va='bottom',
                    bbox=dict(facecolor='none', edgecolor='none', alpha=0.5, pad=0.5)
                )
            if (mz_min <= peaks[k].local_max_mz <= mz_max and rt_min <= peaks[k].local_max_rt <= rt_max):
                plot['img_plot'].text(
                    peaks[k].local_max_mz + 0.01,  # adjust offset as needed
                    peaks[k].local_max_rt - 0.01,  # adjust offset as needed
                    str(k),
                    fontsize=10,
                    color="#1ACCEF",
                    ha='right',
                    va='bottom',
                    bbox=dict(facecolor='none', edgecolor='none', alpha=0.5, pad=0.5)
                )
    
    # Set axis ranges based on mesh bins
    plot['img_plot'].set_xlim([min(mesh.bins_mz), max(mesh.bins_mz)])
    plot['img_plot'].set_ylim([min(mesh.bins_rt), max(mesh.bins_rt)])
    
    return plot

def plot_meshPeaks_interactive(mesh, peaks, localMax=[], size=(1200, 600), withFailedPeaks=False, showIndeces=False):
    """
    Interactive Plotly version of plot_meshPeaks.
    
    Parameters:
        rt (list or array): Retention times
        mz (list or array): m/z values
        img (2D array): Smoothed intensity grid
        peaks (list of tuples): List of (rt, mz) for detected peaks
        fitted (list of tuples): List of (rt, mz) for fitted Gaussian peaks
    """

    # Prepare mesh grid data
    img = np.array(mesh.data).reshape(mesh.m, mesh.n)
    mz = np.array(mesh.bins_mz)
    rt = np.array(mesh.bins_rt)
    fig = go.Figure()

    # Add heatmap for mesh
    fig.add_trace(go.Heatmap(
        z=img,
        x=mz,
        y=rt,
        colorscale='Viridis',
        zsmooth='best',
        colorbar=dict(
            title='Intensity grid',
            x=1.02,      # Right of plot
            y=0.80,      # Push lower under legend
            xanchor='left',
            yanchor='top',
            len=0.80,    # Shorter bar
            thickness=30
        ),
        showscale=True
    ))

    if not withFailedPeaks:
        # Prepare peak data
        quant_indices = [i for i, p in enumerate(peaks) if p.fit_failure_code == 0]
        mzVec = [peaks[i].fitted_mz for i in quant_indices]
        rtVec = [peaks[i].fitted_rt for i in quant_indices]
        mzVecMax = [peaks[i].local_max_mz for i in quant_indices]
        rtVecMax = [peaks[i].local_max_rt for i in quant_indices]
        # Plot peaks
        fig.add_trace(go.Scatter(
            x=mzVec, y=rtVec, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='red', size=10, symbol='circle'),
            name='Peak centers',
            text=[str(i) for i in range(len(peaks))] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))
        fig.add_trace(go.Scatter(
            x=mzVecMax, y=rtVecMax, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='green', size=10, symbol='circle-open'),
            name='Local maxima',
            text=[str(i) for i in range(len(peaks))] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))
    else:
        # Quantified peaks (fit_failure_code == 0)
        quant_indices = [i for i, p in enumerate(peaks) if p.fit_failure_code == 0]
        mzVec = [peaks[i].fitted_mz for i in quant_indices]
        rtVec = [peaks[i].fitted_rt for i in quant_indices]
        mzVecMax = [peaks[i].local_max_mz for i in quant_indices]
        rtVecMax = [peaks[i].local_max_rt for i in quant_indices]

        fig.add_trace(go.Scatter(
            x=mzVec, y=rtVec, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='red', size=10, symbol='circle'),
            name='Quantified Peak centers',
            text=[str(i) for i in quant_indices] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))
        fig.add_trace(go.Scatter(
            x=mzVecMax, y=rtVecMax, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='green', size=10, symbol='circle-open'),
            name='Quantified Local maxima',
            text=[str(i) for i in quant_indices] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))

        # Non-quantified peaks (fit_failure_code != 0)
        nonquant_indices = [i for i, p in enumerate(peaks) if p.fit_failure_code != 0]
        mzVecnF = [peaks[i].fitted_mz for i in nonquant_indices]
        rtVecnF = [peaks[i].fitted_rt for i in nonquant_indices]
        mzVecMaxnF = [peaks[i].local_max_mz for i in nonquant_indices]
        rtVecMaxnF = [peaks[i].local_max_rt for i in nonquant_indices]

        fig.add_trace(go.Scatter(
            x=mzVecnF, y=rtVecnF, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='purple', size=10, symbol='circle'),
            name='Non-quantified Peak centers',
            text=[str(i) for i in nonquant_indices] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))
        fig.add_trace(go.Scatter(
            x=mzVecMaxnF, y=rtVecMaxnF, mode='markers+text' if showIndeces else 'markers',
            marker=dict(color='yellow', size=10, symbol='circle-open'),
            name='Non-quantified Local maxima',
            text=[str(i) for i in nonquant_indices] if showIndeces else None,
            textposition='top right',
            textfont=dict(color='white', size=10, family='Arial')
        ))

    # Update layout
    fig.update_layout(
        width=size[0], height=size[1],
        xaxis_title='m/z',
        yaxis_title='Retention time (s)',
        title=dict(text='Interactive Mesh Peaks',
            pad=dict(t=10)),
        template='plotly_dark',
        xaxis=dict(range=[min(mesh.bins_mz), max(mesh.bins_mz)]),
        yaxis=dict(range=[min(mesh.bins_rt), max(mesh.bins_rt)]),
        title_x=0.5,  # Center the title
        margin=dict(l=60, r=120, t=40, b=50),  # reduce t to bring plot up
        legend=dict(
            x=1.02,  # Place legend outside the plot (right side)
            y=1,     # Top
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(0,0,0,0.5)',
            bordercolor='white',
            borderwidth=1,
            orientation='v'  # Vertical legend
        )
    )

    fig.show()
    return fig

def plot_meshRawPeaks(mesh, lcmsData, peaks, localMax=[], size=(15,10)):
    """
    Plots a mesh and raw data with identified peaks.
    Args:
        mesh (Mesh): The mesh object containing the mesh grid data.
        rawData (raw LC-MS data object): The raw LC-MS data from pastaq.
        peaks (list): A list of Peak objects representing the identified peaks and local maxima in the mesh.
        localMax (list): A list of local maxima to be plotted.
    Returns:
        Figure: The scatter plot figure.
    Raises:
        None
    """
    
    rawData =[[], [], []] # RT, mz, Intensity
    for i in range(0,len(lcmsData.scans),1):
        rawData[0] = rawData[0] + [lcmsData.scans[i].retention_time] * (lcmsData.scans[i].num_points)
        rawData[1] = rawData[1] + lcmsData.scans[i].mz
        rawData[2] = rawData[2] + [math.sqrt(x) for x in lcmsData.scans[i].intensity]
    # print(f"Length of vectors for the raw data scatter plot: {len(rawData[0])} for retention time, {len(rawData[1])} for m/z, {len(rawData[2])} for intensity.") # comment out to check the number of elemnts in mz, rt and intensity for the scatter plot
    
    scatterRaw = plt.figure(figsize=size, facecolor='black')  # Set the figure size
    img = mesh.data
    img = np.reshape(img, (mesh.m, mesh.n))
    offset_rt = (np.array(mesh.bins_rt).max() -
                    np.array(mesh.bins_rt).min())/mesh.m / 2
    offset_mz = (np.array(mesh.bins_mz).max() -
                    np.array(mesh.bins_mz).min())/mesh.n / 2
    
    plt.pcolormesh(
                np.array(mesh.bins_mz) + offset_mz,
                np.array(mesh.bins_rt) + offset_rt,
                img,
                snap=True,
                cmap = 'viridis',
                norm=colors.PowerNorm(gamma=1./2.),
                alpha=0.6)
    cbarMesh = plt.colorbar()
    cbarMesh.set_label('Intensity grid')
    plt.xlim([np.array(mesh.bins_mz).min() - offset_mz,
                        np.array(mesh.bins_mz).max() + offset_mz])
    plt.ylim([np.array(mesh.bins_rt).min() -+ offset_rt,
                        np.array(mesh.bins_rt).max() + offset_rt])

    # axRaw = scatterRaw.add_subplot(111)
    plt.scatter(rawData[1], rawData[0], c = rawData[2], s = 5, alpha = 1, cmap = 'twilight_shifted', marker='.')

    # Add color bar to show elevation scale
    cbar = plt.colorbar()
    cbar.set_label('Intensity raw data')

    # shows local maxima
    mzVecLoc, rtVecLoc, mzVec, rtVec = [], [], [], []
    for k in range(len(peaks)):
        if peaks[k].fit_failure_code == 0:
            mzVecLoc.insert(k, peaks[k].local_max_mz)
            rtVecLoc.insert(k, peaks[k].local_max_rt)
    plt.scatter(mzVecLoc, rtVecLoc, s=100, c='green', marker='o', alpha=0.7)

    # shows identified peaks
    for k in range(len(peaks)):
        if peaks[k].fit_failure_code == 0:
            mzVec.insert(k, peaks[k].fitted_mz)
            rtVec.insert(k, peaks[k].fitted_rt)
    plt.scatter(mzVec, rtVec, s=100, c='red', marker='.')
    
    # shows local maxima
    if localMax != []:
        mzVecLocMax, rtVecLocMax = [], []
        for k in range(len(localMax)):
            mzVecLocMax.insert(k, localMax[k].mz)
            rtVecLocMax.insert(k, localMax[k].rt)
        plt.scatter(mzVecLocMax, rtVecLocMax, s=50, c='blue', marker='x', alpha=0.7)

    # Label the axes
    plt.xlabel('m/z')  # Set the x-axis label
    plt.ylabel('Retention time (s)')  # Set the y-axis label
    plt.title('Raw data of a zommed in area of the isotope cluster')  # Set the title

    # Show the plot
    plt.show()
    return scatterRaw

#print 64 bit uint as binary with bytes representation
def print_64bit_as_bytes(val):
    # Format as 64-bit binary, pad with zeros
    bin_str = f"{val:064b}"
    # Group into bytes (8 bits)
    grouped = ' '.join(bin_str[i:i+8] for i in range(0, 64, 8))
    # print(grouped)
    return grouped

def downloadTestFilesFromOSF():
    """
    Download .mzML test files used in this tutorial from OSF and save them in the 'data' directory.
    If the files already exist with correct size, skip the download.
    Only download files that are missing or have incorrect size.
    """

    fileList = [
        {"name": "3_1_extract_3400_4100_590_615.mzML", "size": 29328010, "centroid": False},
        {"name": "3_2_extract_3400_4100_590_615.mzML", "size": 26397653, "centroid": False},
        {"name": "3_1_extract_3400_4100_590_615_c.mzML", "size": 10342153, "centroid": True},
        {"name": "3_2_extract_3400_4100_590_615_c.mzML", "size": 9465509, "centroid": True},
    ]
    # Target directory for data
    target_dir = os.path.join(os.getcwd(), "data")
    os.makedirs(target_dir, exist_ok=True)

    # OSF API endpoint for files in the storage
    project_id = "b2he9"
    url = f"https://api.osf.io/v2/nodes/{project_id}/files/osfstorage/"

    print("Checking mzML files in 'data/'...")

    # Get the file list from OSF
    resp = requests.get(url)
    resp.raise_for_status()
    files = resp.json()['data']

    # Build a dict for quick lookup of OSF files
    osf_files = {file['attributes']['name']: file
                 for file in files
                 if file['attributes']['name'].endswith('.mzML')}

    # Track files to download
    files_to_download = []

    for file_info in fileList:
        fname = file_info["name"]
        expected_size = file_info["size"]
        local_path = os.path.join(target_dir, fname)
        needs_download = False

        if not os.path.exists(local_path):
            print(f"{fname} is missing, will download.")
            needs_download = True
        else:
            local_size = os.path.getsize(local_path)
            if local_size != expected_size:
                print(f"{fname} exists but has incorrect size ({local_size} != {expected_size}), will re-download.")
                needs_download = True

        if needs_download:
            if fname in osf_files:
                files_to_download.append(fname)
            else:
                print(f"{fname} not found on OSF, skipping.")

    # Download missing/incorrect files
    for fname in files_to_download:
        file = osf_files[fname]
        download_url = file['links']['download']
        out_path = os.path.join(target_dir, fname)
        print(f"Downloading {fname}...")
        r = requests.get(download_url)
        if r.status_code == 200:
            with open(out_path, 'wb') as f:
                f.write(r.content)
            print(f"Downloaded {fname}.")
        else:
            print(f"Failed to download {fname}")

    if not files_to_download:
        print("All mzML files are present and have correct size. No download needed.")

def XICVisualisation(lcmsData, grid, peak, sigmaFactor, rtWidthEIC):
    """
    Visualize the Extracted Ion Chromatogram (XIC) for a given peak.

    Parameters:
    - lcmsData: The LCMS data object.
    - grid: The grid object for the XIC visualization.
    - peak: The detected peak to be shown.
    - sigmaFactor: The factor to determine the XIC range.
    - rtWidthEIC: The width of the XIC in seconds.
    """
    # Calculate the mz and rt range for the XIC visualization
    lowBoundmz = peak.fitted_mz - sigmaFactor * peak.fitted_sigma_mz
    highBoundmz = peak.fitted_mz + sigmaFactor * peak.fitted_sigma_mz
    lowBoundrt = peak.fitted_rt - rtWidthEIC
    highBoundrt = peak.fitted_rt + rtWidthEIC
    
    # plot XIC for the raw data
    xic = pq.xic(lcmsData, lowBoundmz, highBoundmz, lowBoundrt, highBoundrt)
    
    # plot XIC for the smoothed data
    np.array(grid.bins_mz, dtype=np.float64)
    indices_mz = [i for i, val in enumerate(np.array(grid.bins_mz, dtype=np.float64)) if lowBoundmz < val < highBoundmz]
    indices_rt = [i for i, val in enumerate(np.array(grid.bins_rt, dtype=np.float64)) if lowBoundrt < val < highBoundrt]
    xicGrid = np.array(grid.data).reshape(grid.m, grid.n)
    xicGridi = xicGrid[indices_rt][:, indices_mz]
    xicGridis = np.sum(xicGridi, axis=1)
    xicGridrt = np.array(grid.bins_rt, dtype=np.float64)[indices_rt]

    plt.figure(figsize=(12, 6))  # Set the figure size
    plt.plot(xic.retention_time, xic.intensity, label='Raw data', alpha=1, color='blue')
    plt.plot(xicGridrt, xicGridis, label='Smoothed data', alpha=1, color='red')
    plt.axvline(x=peak.fitted_rt, label='Fitted peak rt', color='green', alpha=0.75, linestyle='--', linewidth=2)
    plt.plot(peak.fitted_rt, peak.fitted_height, label='Fitted Gaussian height', marker='o', markersize=7, color='red', alpha=0.75)
    
    # Plot the Gaussian peak
    plot_gaussian(peak.fitted_rt, peak.fitted_sigma_rt, peak.fitted_height)

    # Plot the XIC
    plt.title(f'Extracted Ion Chromatogram for Peak {peak.id} at {peak.fitted_rt:.0f} seconds and {peak.fitted_mz:.4f} m/z')
    plt.xlabel('Retention Time (s)')
    plt.ylabel('Intensity (cps)')
    plt.legend()
    plt.show()

def showXICRAWAligned(lcmsData1, lcmsData2, peak1, peaks2, sigmaFactor, rtWidthXIC, time_map_1r_2s):
    """
    Show the raw XIC alignment between two LCMS datasets for a specific peak in the reference chromatogram.

    Args:
        lcmsData1: The reference LCMS dataset.
        lcmsData2: The sample LCMS dataset.
        peak1: The peak in the reference dataset.
        peaks2: The peaklist of the sample dataset.
        sigmaFactor: The sigma factor for the XIC.
        rtWidthXIC: The retention time width for the XIC.
        time_map_1r_2s: The time mapping from the first to the sample dataset.
    """
    lowBoundmz = peak1.fitted_mz - sigmaFactor * peak1.fitted_sigma_mz
    highBoundmz = peak1.fitted_mz + sigmaFactor * peak1.fitted_sigma_mz
    lowBoundrt = peak1.fitted_rt - rtWidthXIC
    highBoundrt = peak1.fitted_rt + rtWidthXIC

    # get extracted ion chromatogram for file 2 using the same mz and rt range as for file 1
    xic1 = pq.xic(lcmsData1, lowBoundmz, highBoundmz, lowBoundrt, highBoundrt)
    xic2 = pq.xic(lcmsData2, lowBoundmz, highBoundmz, lowBoundrt, highBoundrt)

    # get matching peak for peaks1[peakIdx] in peaks2_1r
    matchedpeak = find_closest_peak(peak1, peaks2)

    # print target peaks and matched peaks in the other chromatograms
    print("Raw data similarity: {}" + format(peak1))
    print("Raw data similarity: {}" + format(matchedpeak))

    alignedRTs = interpolate_values(xic2.retention_time, time_map_1r_2s.sample_rt_start, time_map_1r_2s.sample_rt_end, time_map_1r_2s.rt_start, time_map_1r_2s.rt_end)

    #plot the raw EIC in the reference (blue) and sample data (red), and the aligned EIC in the sample data (green)
    plt.figure(figsize=(12, 6))  # Set the figure size
    plt.plot(xic1.retention_time, xic1.intensity, label='Raw reference', alpha=1, color='blue')
    plt.plot(xic2.retention_time, xic2.intensity, label='Raw sample', alpha=1, color='red')
    plt.plot(alignedRTs, xic2.intensity, label='Raw aligned sample', alpha=1, color='green')
    plt.axvline(x = peak1.fitted_rt, label='fitted peak rt in reference', color='blue', alpha=0.75, linestyle='--', linewidth=2)
    if matchedpeak!=None:
        plt.axvline(x = matchedpeak.fitted_rt, label='fitted peak rt in sample data', color='red', alpha=0.75, linestyle='--', linewidth=2)
        plt.axvline(x = (matchedpeak.fitted_rt + matchedpeak.rt_delta), label='fitted peak rt in sample after alignment', color='green', alpha=0.75, linestyle='-.', linewidth=2)
    plt.plot(peak1.fitted_rt, peak1.fitted_height, label='fitted Gaussian height', marker='o', markersize=7, color='blue', alpha=0.75)
    if matchedpeak!=None:
        plt.plot(matchedpeak.fitted_rt, matchedpeak.fitted_height, label='fitted Gaussian height', marker='o', markersize=7, color='red', alpha=0.75)
        plt.plot((matchedpeak.fitted_rt + matchedpeak.rt_delta), matchedpeak.fitted_height, label='fitted Gaussian height', marker='o', markersize=7, color='green', alpha=0.75)
    plt.xlabel('rt (sec)')
    plt.ylabel('intensity (cps)')
    plt.title(f'Extracted ion chromatogram of reference XIC for peak {peak1.id} at {peak1.fitted_rt} sec and sample XIC for '
              f'peak {matchedpeak.id if matchedpeak else 'N/A'} at {matchedpeak.fitted_rt if matchedpeak else 'N/A'} sec.')
    plt.legend()

    # show the plot
    plt.show()
    
    
def plot_feature_zoom(grid, feature, peaks, mzm, mzp, rtm, rtp):
    """
    Visualize a zoomed region around a feature (isotope cluster) in the smoothed 2D grid.

    Parameters:
        grid: PASTAQ grid object (with bins_mz and bins_rt)
        feature: one feature object (with .peak_ids)
        peaks: list of peak objects (with .fitted_mz and .fitted_rt)
        mzm: float, Da added to monoisotopic peak m/z (positive extension)
        mzp: float, Da subtracted from monoisotopic peak m/z (negative extension)
        rtm: float, seconds subtracted from monoisotopic peak rt (negative extension)
        rtp: float, seconds added to monoisotopic peak rt (positive extension)
    """

    # Use the first peak in the feature as the monoisotopic peak
    mono_idx = feature.peak_ids[0]
    mono_peak = peaks[mono_idx]

    plot1 = pq.plot_mesh(grid, transform='sqrt', figure=None)
    plot1['img_plot'].get_figure().set_size_inches(15, 10)
    
    # Set axis limits for zoom
    plot1['img_plot'].set_xlim(mono_peak.fitted_mz - mzm, mono_peak.fitted_mz + mzp)
    plot1['img_plot'].set_ylim(mono_peak.fitted_rt - rtm, mono_peak.fitted_rt + rtp)
    plot1['mz_plot'].set_xlim(mono_peak.fitted_mz - mzm, mono_peak.fitted_mz + mzp)
    plot1['rt_plot'].set_ylim(mono_peak.fitted_rt - rtm, mono_peak.fitted_rt + rtp)

    # Prepare vectors for isotopologue peaks
    mzVec = []
    rtVec = []
    for i in feature.peak_ids:
        mzVec.append(peaks[i].fitted_mz)
        rtVec.append(peaks[i].fitted_rt)
    plot1['img_plot'].scatter(mzVec, rtVec, s=100, c='red', marker='.')

    # Adjust color limits for zoomed region
    for artist in plot1['img_plot'].get_children():
        if hasattr(artist, 'get_array'):
            img_data = artist.get_array()
            mz_grid = grid.bins_mz
            rt_grid = grid.bins_rt
            xlim = plot1['img_plot'].get_xlim()
            ylim = plot1['img_plot'].get_ylim()
            mz_mask = (mz_grid >= xlim[0]) & (mz_grid <= xlim[1])
            rt_mask = (rt_grid >= ylim[0]) & (rt_grid <= ylim[1])
            zoomed_data = img_data[np.ix_(rt_mask, mz_mask)]
            min_zoom = np.min(zoomed_data)
            max_zoom = np.max(zoomed_data)
            artist.set_clim(vmin=min_zoom, vmax=max_zoom)
            break

    return