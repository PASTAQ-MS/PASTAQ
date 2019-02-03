import math

from .tapp import *

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

def calculate_mz_bins_1(raw_data, num_sampling_mz):
    current_mz = raw_data.min_mz
    mz = [current_mz]

    fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    A = fwhm_ref / ((raw_data.reference_mz) ** (3/2))
    while current_mz < raw_data.max_mz:
        # NOTE: This generates approximately 10 sampling points between the
        # peak's FWHM, not on the base.
        current_mz = current_mz + A * (current_mz) ** (3/2) / num_sampling_mz 
        mz = mz + [current_mz]
    return mz

def calculate_mz_bins_2(raw_data, num_sampling_mz):
    current_mz = raw_data.min_mz
    mz = [current_mz]

    fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    A = fwhm_ref / ((raw_data.reference_mz) ** (3/2))
    n = 0
    while current_mz < raw_data.max_mz:
        # NOTE: This generates approximately 10 sampling points between the
        # peak's FWHM, not on the base.
        current_mz = mz_at(raw_data, num_sampling_mz, n)
        n = n + 1
        mz = mz + [current_mz]
    return mz

RawData.tic = tic
RawData.mz_bins_1 = calculate_mz_bins_1
RawData.mz_bins_2 = calculate_mz_bins_2
