# helper functions for PASTAQ demonstration
import math as math
import numpy as np
import matplotlib.pyplot as plt

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