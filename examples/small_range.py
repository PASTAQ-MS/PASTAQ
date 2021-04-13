import pastaq
import os
import matplotlib.pyplot as plt

# Only the first two chromatogram will be considered for retention time
# alignment.
input_files = [
        {'reference': True,  'raw_path': '/path/to/mzxml/1_1.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_1.mzid'},
        {'reference': True,  'raw_path': '/path/to/mzxml/1_2.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_2.mzid'},
        {'reference': False, 'raw_path': '/path/to/mzxml/1_3.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_3.mzid'},
]

# Launch the DDA pipeline for the selected files for data acquired with an
# Orbitrap instrument at 70000 resolution at 200 m/z, and with an estimated
# chromatographic peak width of 9 seconds (FWHM).
params = pastaq.default_parameters('orbitrap', 9)
params['resolution_ms1'] = 70000
params['reference_mz'] = 200

# Reading a small range of the raw data for exploration.
params['min_mz'] = 708.0
params['max_mz'] = 717.0
params['min_rt'] = 2300.0
params['max_rt'] = 3100.0
params['max_peaks'] = 1000
params['polarity'] = 'pos'

# In a small retention time range, we can configure the retention time alignment
# as follows to speed up the processing.
params['warp2d_rt_expand_factor'] = 0.5
params['warp2d_num_points'] = 500
params['warp2d_slack'] = 100
params['warp2d_window_size'] = 100
params['warp2d_peaks_per_window'] = 20

# Run the DDA pipeline. Since the mz-rt range is small, We also enable saving
# the smoothed 2D map for further visualization.
pastaq.dda_pipeline(params, input_files, 'small_range', save_grid=True)

# Visualize the smoothed map for each file in a 2D image. The images will be
# saved in the quality directory with the rest of quality control plots.
for file in input_files:
    base_name = os.path.basename(file['raw_path'])
    base_name = os.path.splitext(base_name)
    stem = base_name[0]
    grid = pastaq.read_grid('small_range/grid/{}.grid'.format(stem))
    fig = plt.figure()
    pastaq.plot_mesh(grid, figure=fig)
    plt.savefig('small_range/quality/{}_grid.png'.format(stem), dpi=300)
    plt.close(fig)

