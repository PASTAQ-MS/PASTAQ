import tapp

# Parameters
file_names = ["file_01.mzXML", "file_02.mzXML", "file_03.mzXML"]
mz_min = 4695
mz_max = 6195
rt_min = 665
rt_max = 675.71
instrument = "ORBITRAP"
sigma_mz = (200, 0.001)
sigma_rt = 15

print("01 -- Reading raw data")
raw_data = []
for file_name in file_names:
    raw_data = raw_data + \
        [tapp.raw_data(file_name, mz_min, mz_max, rt_min, rt_max)]

print("02 -- Creating grids")
grids = []
for raw_file in raw_data:
    grids = grids + [tapp.grid(raw_file, instrument,
                               sigma_mz, sigma_rt, warped=False, parallel=True)]

print("03 -- Finding peaks")
peak_lists = []
for grid_file in grids:
    peak_lists = peak_lists + [tapp.centroid(grid_file, threshold=0,
                                             max_peaks=1000000, parallel=True)]

print("04 -- Aligning peaks")
aligned_peak_lists = tapp.warp2D(
    peak_lists, ref_index=0,
    slack=30, window_size=50, peaks_per_window=50, num_points=2000)

print("05 -- Matching relevant peaks")
metamatch = tapp.metamatch(
    aligned_peak_lists, radius_mz=0.001, radius_rt=15, fraction=0.7)

quantitative_table = metamatch.table(metric="height")
