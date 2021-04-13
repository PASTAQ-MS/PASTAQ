import pastaq

input_files = [
        {'reference': False,  'raw_path': '/path/to/mzxml/1_1.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_1.mzid'},
        {'reference': False,  'raw_path': '/path/to/mzxml/1_2.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_2.mzid'},
        {'reference': False, 'raw_path':  '/path/to/mzxml/1_3.mzXML', 'group': 'a', 'ident_path': '/path/to/mzidentml/1_3.mzid'},
]

# Launch the DDA pipeline for the selected files for data acquired with an
# Orbitrap instrument at 70000 resolution at 200 m/z, and with an estimated
# chromatographic peak width of 9 seconds (FWHM).
params = pastaq.default_parameters('orbitrap', 9)
params['resolution_ms1'] = 70000
params['reference_mz'] = 200
pastaq.dda_pipeline(params, input_files, 'full_analysis')
