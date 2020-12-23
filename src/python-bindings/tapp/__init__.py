from .tapp import *
import tapp
import os
import json
import logging
import time
import datetime

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
plt.rcParams.update({
    'font.family': 'sans-serif',
    'figure.figsize': (7.08661, 7.08661/1.618034),
    'legend.fontsize': 6,
    'font.size': 7,
})


def find_sequence_consensus(annotations, sequence_column, min_consensus_count):
    consensus = annotations[["cluster_id", "file_id", sequence_column]]
    consensus = consensus.drop_duplicates()
    consensus.columns = ["cluster_id", "file_id", "consensus_sequence"]
    consensus = consensus.groupby(["cluster_id", "consensus_sequence"]).agg(
        consensus_count=('consensus_sequence', 'count'))
    max_count = consensus.groupby("cluster_id").max().reset_index()
    max_count.columns = ["cluster_id", "consensus_count_max"]
    consensus = pd.merge(consensus.reset_index(), max_count, on="cluster_id")
    consensus = consensus[consensus["consensus_count"]
                          == consensus["consensus_count_max"]]
    consensus = consensus[consensus["consensus_count"] >= min_consensus_count]
    consensus = consensus.drop(["consensus_count_max"], axis=1)
    return consensus

def find_protein_groups(annotations, sequence_column,  protein_name_column, protein_description_column):
    seq_prot = annotations.copy()
    seq_prot = seq_prot[[sequence_column, protein_name_column, protein_description_column]]
    seq_prot = seq_prot.drop_duplicates()
    seq_prot = seq_prot.dropna()
    seq_prot = seq_prot.astype(str)
    seq_prot['combined_protein_name'] = seq_prot[protein_name_column] + "__" + seq_prot[protein_description_column]
    seq_prot_copy = seq_prot.copy()
    protein_groups = pd.DataFrame(columns=['protein_group', 'combined_protein_name'])
    protein_group_index = 0
    while not seq_prot.empty:
        cur_target_proteins = seq_prot.iloc[0]['combined_protein_name']
        cur_target_peptides = seq_prot['combined_protein_name'] == cur_target_proteins
        cur_target_peptides = seq_prot[sequence_column][cur_target_peptides]
        cur_target_peptides = cur_target_peptides.drop_duplicates()
        cur_peptides_len = cur_target_peptides.shape[0]
        prev_peptides_len = 0
        while cur_peptides_len != prev_peptides_len:
            prev_peptides_len = cur_peptides_len
            cur_target_proteins = seq_prot[seq_prot[sequence_column].isin(cur_target_peptides)]
            cur_target_proteins = cur_target_proteins['combined_protein_name'].drop_duplicates()
            cur_target_peptides = seq_prot['combined_protein_name'].isin(cur_target_proteins)
            cur_target_peptides = seq_prot[sequence_column][cur_target_peptides]
            cur_target_peptides = cur_target_peptides.drop_duplicates()
            cur_peptides_len = cur_target_peptides.shape[0]
        protein_group = pd.DataFrame({
            'protein_group': np.repeat(protein_group_index, len(cur_target_proteins)),
            'combined_protein_name': cur_target_proteins.values,
            })
        protein_groups = pd.concat([protein_groups, protein_group])
        protein_group_index += 1
        seq_prot = seq_prot[~seq_prot['combined_protein_name'].isin(cur_target_proteins)]
    seq_prot_copy = pd.merge(seq_prot_copy, protein_groups, on='combined_protein_name')
    seq_prot_copy = seq_prot_copy[[sequence_column, 'protein_group']].drop_duplicates()
    annotations = pd.merge(annotations, seq_prot_copy, on=sequence_column, how='left')
    return annotations

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


Peak.plot_xic = plot_xic


def plot_raw_points(
    peak,
    raw_data,
    img_plot=None,
    rt_plot=None,
    mz_plot=None,
        xic_method="max"):
    data_points = raw_data.raw_points(
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


# TODO: Probably we don't want avg_fwhm_rt to be a parameter being passed on
# this function. Just set it to a resonable level for the default parameters and
# modify it later as needed.
def default_parameters(instrument, avg_fwhm_rt):
    if instrument == 'orbitrap':
        tapp_parameters = {
            #
            # Instrument configuration.
            #
            'instrument_type': 'orbitrap',
            'resolution_ms1': 70000,
            'resolution_msn': 30000,
            'reference_mz': 200,
            'avg_fwhm_rt': avg_fwhm_rt,
            #
            # Meshing.
            #
            'num_samples_mz': 5,
            'num_samples_rt': 5,
            'smoothing_coefficient_mz': 0.4,
            'smoothing_coefficient_rt': 0.4,
            #
            # Warp2D.
            #
            'warp2d_slack': 30,
            'warp2d_window_size': 50,
            'warp2d_num_points': 2000,
            'warp2d_rt_expand_factor': 0.2,
            'warp2d_peaks_per_window': 100,
            #
            # MetaMatch.
            #
            'metamatch_fraction': 0.7,
            'metamatch_n_sig_mz': 1.5,
            'metamatch_n_sig_rt': 1.5,
            #
            # Feature detection.
            #
            'feature_detection_charge_states': [5, 4, 3, 2, 1],
            #
            # Other.
            #
            'max_peaks': 1000000,
            'polarity': 'both',
            'min_mz': 0,
            'max_mz': 100000,
            'min_rt': 0,
            'max_rt': 100000,
            #
            # Identification.
            #
            # Keep only the max rank PSM.
            'ident_max_rank_only': True,
            # Whether to filter the PSM that don't pass the FDR threshold.
            'ident_require_threshold': True,
            # Ignore PSMs that have been marked as decoy.
            'ident_ignore_decoy': True,
            #
            # Quality.
            #
            'similarity_num_peaks': 2000,
            #
            # Quantitative table generation.
            #
            # Options: 'height', 'volume'
            'quant_isotopes': 'height',
            # Options: 'monoisotopic_height', 'monoisotopic_volume',
            #          'total_height', 'total_volume',
            #          'max_height', 'max_volume',
            'quant_features': 'max_height',
            # Whether to remove feature annotations where the charge state of
            # the detected feature doesn't match the one given by the
            # identification engine.
            'quant_features_charge_state_filter': True,
            # Options: 'theoretical_mz', 'msms_event'
            'quant_ident_linkage': 'theoretical_mz',
            # Whether to obtain a consensus sequence and proteins on identifications.
            'quant_consensus': True,
            # Demand a minimum number of files with identification per cluster.
            'quant_consensus_min_ident': 2,
            # Whether to store all the annotations prior to cluster aggregation.
            'quant_save_all_annotations': True,
            # TODO: Protein inference method:
            #     - 'none': Don't perform protein inference.
            #     - 'general_razor': All files are considered for Occam's razor inference.
            #     - 'group_razor': Occam's razor inference is performed per group.
            #     - 'file_razor': Occam's razor inference is performed per file.
            'quant_inference_type': 'group_razor',
        }
        return tapp_parameters


def dda_pipeline_summary(tapp_parameters, input_stems, output_dir):
    summary_log = logging.getLogger('summary')
    summary_log.setLevel(logging.INFO)
    summary_log_fh = logging.FileHandler(os.path.join(output_dir, 'summary.log'))
    summary_log_fh.setLevel(logging.INFO)
    summary_log_formatter = logging.Formatter('%(message)s')
    summary_log_fh.setFormatter(summary_log_formatter)
    summary_log.addHandler(summary_log_fh)

    # Raw data
    summary_log.info('Raw data')
    for stem in input_stems:
        summary_log.info('    {}'.format(stem))

        # MS1
        in_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        if os.path.exists(in_path):
            raw_data = tapp.read_raw_data(in_path)
            summary_log.info('        MS1')
            summary_log.info('            number of scans: {}'.format(len(raw_data.scans)))
            summary_log.info('            min_mz: {}'.format(raw_data.min_mz))
            summary_log.info('            max_mz: {}'.format(raw_data.max_mz))
            summary_log.info('            min_rt: {}'.format(raw_data.min_rt))
            summary_log.info('            max_rt: {}'.format(raw_data.max_rt))

        # MS2
        in_path = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        if os.path.exists(in_path):
            raw_data = tapp.read_raw_data(in_path)
            summary_log.info('        MS2')
            summary_log.info('            number of scans: {}'.format(len(raw_data.scans)))
            summary_log.info('            min_mz: {}'.format(raw_data.min_mz))
            summary_log.info('            max_mz: {}'.format(raw_data.max_mz))
            summary_log.info('            min_rt: {}'.format(raw_data.min_rt))
            summary_log.info('            max_rt: {}'.format(raw_data.max_rt))

    # Peaks
    summary_log.info('Peaks detected')
    avg_peak_heights = []
    median_peak_heights = []
    std_peak_heights = []
    n_peaks = []
    for stem in input_stems:
        summary_log.info('    {}'.format(stem))

        in_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
        if os.path.exists(in_path):
            peaks = tapp.read_peaks(in_path)
            peak_heights = np.array([peak.fitted_height for peak in peaks])
            n_peaks += [len(peaks)]
            mean_height = peak_heights.mean()
            median_height = np.median(peak_heights)
            std_height = np.std(peak_heights)
            avg_peak_heights += [mean_height]
            median_peak_heights += [median_height]
            std_peak_heights += [std_height]
            summary_log.info('        Number of peaks: {}'.format(len(peaks)))
            summary_log.info('        Fitted height')
            summary_log.info('            mean: {}'.format(mean_height))
            summary_log.info('            median: {}'.format(median_height))
            summary_log.info('            std: {}'.format(std_height))

    if len(n_peaks) != 0:
        summary_log.info('    Overall average')
        summary_log.info('        Number of peaks: {}'.format(np.mean(n_peaks)))
        summary_log.info('        Fitted height')
        summary_log.info('            mean: {}'.format(np.mean(avg_peak_heights)))
        summary_log.info('            median: {}'.format(np.mean(median_peak_heights)))
        summary_log.info('            std: {}'.format(np.mean(std_peak_heights)))

    # Feature detection
    summary_log.info('Feature detection')
    avg_feature_monoisotopic_heights = []
    median_feature_monoisotopic_heights = []
    std_feature_monoisotopic_heights = []
    avg_feature_max_heights = []
    median_feature_max_heights = []
    std_feature_max_heights = []
    avg_feature_total_heights = []
    median_feature_total_heights = []
    std_feature_total_heights = []
    n_features = []
    for stem in input_stems:
        summary_log.info('    {}'.format(stem))

        in_path = os.path.join(output_dir, 'features', "{}.features".format(stem))
        if os.path.exists(in_path):
            features = tapp.read_features(in_path)

            feature_max_heights = np.array([feature.max_height for feature in features])
            feature_monoisotopic_heights = np.array([feature.monoisotopic_height for feature in features])
            feature_total_heights = np.array([feature.total_height for feature in features])

            feature_mean_max_height = feature_max_heights.mean()
            feature_median_max_height = np.median(feature_max_heights)
            feature_std_max_height = np.std(feature_max_heights)

            feature_mean_monoisotopic_height = feature_monoisotopic_heights.mean()
            feature_median_monoisotopic_height = np.median(feature_monoisotopic_heights)
            feature_std_monoisotopic_height = np.std(feature_monoisotopic_heights)

            feature_mean_total_height = feature_total_heights.mean()
            feature_median_total_height = np.median(feature_total_heights)
            feature_std_total_height = np.std(feature_total_heights)

            n_features += [len(features)]
            avg_feature_max_heights += [feature_mean_max_height]
            median_feature_max_heights += [feature_median_max_height]
            std_feature_max_heights += [feature_std_max_height]

            avg_feature_monoisotopic_heights += [feature_mean_monoisotopic_height]
            median_feature_monoisotopic_heights += [feature_median_monoisotopic_height]
            std_feature_monoisotopic_heights += [feature_std_monoisotopic_height]

            avg_feature_total_heights += [feature_mean_total_height]
            median_feature_total_heights += [feature_median_total_height]
            std_feature_total_heights += [feature_std_total_height]

            summary_log.info('        Number of features: {}'.format(len(features)))
            summary_log.info('        Max height')
            summary_log.info('            mean: {}'.format(feature_mean_max_height))
            summary_log.info('            median: {}'.format(feature_median_max_height))
            summary_log.info('            std: {}'.format(feature_std_max_height))
            summary_log.info('        Monoisotopic height')
            summary_log.info('            mean: {}'.format(feature_mean_monoisotopic_height))
            summary_log.info('            median: {}'.format(feature_median_monoisotopic_height))
            summary_log.info('            std: {}'.format(feature_std_monoisotopic_height))
            summary_log.info('        Total height')
            summary_log.info('            mean: {}'.format(feature_mean_total_height))
            summary_log.info('            median: {}'.format(feature_median_total_height))
            summary_log.info('            std: {}'.format(feature_std_total_height))

    if len(n_features) != 0:
        summary_log.info('    Overall average')
        summary_log.info('        Number of features: {}'.format(np.mean(n_features)))
        summary_log.info('        Max height')
        summary_log.info('            mean: {}'.format(np.mean(avg_feature_max_heights)))
        summary_log.info('            median: {}'.format(np.mean(median_feature_max_heights)))
        summary_log.info('            std: {}'.format(np.mean(std_feature_max_heights)))
        summary_log.info('        Monoisotopic height')
        summary_log.info('            mean height: {}'.format(np.mean(avg_feature_monoisotopic_heights)))
        summary_log.info('            median height: {}'.format(np.mean(median_feature_monoisotopic_heights)))
        summary_log.info('            std height: {}'.format(np.mean(std_feature_monoisotopic_heights)))
        summary_log.info('        Total height')
        summary_log.info('            mean height: {}'.format(np.mean(avg_feature_total_heights)))
        summary_log.info('            median height: {}'.format(np.mean(median_feature_total_heights)))
        summary_log.info('            std height: {}'.format(np.mean(std_feature_total_heights)))

    # Identifications and linkage
    summary_log.info('Annotations and linkage')
    for stem in input_stems:
        summary_log.info('    {}'.format(stem))

        in_path_raw_data = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_linked_msms = os.path.join(output_dir, 'linking', "{}.ms2_peak.link".format(stem))
        if os.path.exists(in_path_raw_data) and os.path.exists(in_path_linked_msms):
            raw_data = tapp.read_raw_data(in_path_raw_data)
            linked_msms = tapp.read_linked_msms(in_path_linked_msms)
            summary_log.info('        MS/MS-Peaks linkage')
            summary_log.info('            Number of ms/ms events: {}'.format(len(raw_data.scans)))
            summary_log.info('            Number of ms/ms events linked to peaks: {}'.format(len(linked_msms)))
            summary_log.info('            Linking efficiency (%): {}'.format(len(linked_msms)/len(raw_data.scans) * 100.0))

        in_path_ident_data = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        if os.path.exists(in_path_ident_data):
            ident_data = tapp.read_ident_data(in_path_ident_data)

            in_path_ident_ms2 = os.path.join(output_dir, 'linking', "{}.ident_ms2.link".format(stem))
            if os.path.exists(in_path_ident_ms2):
                ident_ms2 = tapp.read_linked_msms(in_path_ident_ms2)
                summary_log.info('        MS/MS-Identification linkage')
                summary_log.info('            Number of PSMs: {}'.format(len(ident_data.spectrum_matches)))
                summary_log.info('            Number of PSMs linked to MS/MS events: {}'.format(len(ident_ms2)))
                summary_log.info('            PSM-peaks linking efficiency (%): {}'.format(len(ident_ms2)/len(ident_data.spectrum_matches) * 100.0))

            in_path_peak_idents = os.path.join(output_dir, 'linking', "{}.ident_peak.link".format(stem))
            if os.path.exists(in_path_peak_idents):
                ident_peak = tapp.read_linked_psm(in_path_peak_idents)
                summary_log.info('        MS/MS-Identification linkage')
                summary_log.info('            Number of PSMs: {}'.format(len(ident_data.spectrum_matches)))
                summary_log.info('            Number of PSMs linked to peaks: {}'.format(len(ident_peak)))
                summary_log.info('            PSM-peaks linking efficiency (%): {}'.format(len(ident_peak)/len(ident_data.spectrum_matches) * 100.0))

    # TODO: Metamatch stats
    # TODO: Peptide stats
    # TODO: Protein group stats

def dda_pipeline(
    tapp_parameters,
    input_files,
    output_dir="TAPP",
    override_existing=False,
    save_mesh=False,
):
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
    if not os.path.exists(os.path.join(output_dir, 'time_map')):
        os.makedirs(os.path.join(output_dir, 'time_map'))
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

    # Initialize log and parameters files.
    parameters_file_name = os.path.join(output_dir, 'parameters.json')
    with open(parameters_file_name, 'w') as json_file:
        json.dump(tapp_parameters, json_file)

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
        out_path = os.path.join(output_dir, 'peaks', "{}.peaks".format(stem))
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
    out_path = os.path.join(output_dir, 'quality', 'similarity_unwarped')
    logger.info("Starting unwarped similarity matrix calculation")
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or override_existing:
        similarity_matrix = np.zeros(
            len(input_stems) ** 2).reshape(len(input_stems), len(input_stems))
        for i in range(0, len(input_stems)):
            stem_a = input_stems[i]
            logger.info("Reading peaks_a from disk: {}".format(stem_a))
            peaks_a = tapp.read_peaks(os.path.join(
                output_dir, 'peaks', '{}.peaks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'peaks', '{}.peaks'.format(stem_b)))
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
        fig.set_size_inches(7.08661, 7.08661)
        plt.savefig("{}.pdf".format(out_path), dpi=300)
        plt.savefig("{}.png".format(out_path), dpi=300)
        plt.close(fig)
        logger.info(
            "Saving similarity matrix plot: {}.pdf/png".format(out_path))
    logger.info('Finished unwarped similarity matrix calculation in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Correct retention time. If a reference sample is selected it will be
    # used, otherwise, exhaustive warping will be performed.
    # TODO: Allow usage of reference sample.
    out_path = os.path.join(output_dir, 'quality',
                            'similarity_exhaustive_warping')
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
                output_dir, 'peaks', '{}.peaks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'peaks', '{}.peaks'.format(stem_b)))
                logger.info("Warping {} peaks to {}".format(stem_b, stem_a))
                time_map = tapp.calculate_time_map(
                    peaks_a, peaks_b,
                    tapp_parameters['warp2d_slack'],
                    tapp_parameters['warp2d_window_size'],
                    tapp_parameters['warp2d_num_points'],
                    tapp_parameters['warp2d_rt_expand_factor'],
                    tapp_parameters['warp2d_peaks_per_window'])
                peaks_b = tapp.warp_peaks(peaks_b, time_map)
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
        fig.set_size_inches(7.08661, 7.08661)
        plt.savefig("{}.pdf".format(out_path), dpi=300)
        plt.savefig("{}.png".format(out_path), dpi=300)
        plt.close(fig)
        logger.info(
            "Saving similarity matrix plot: {}.pdf/png".format(out_path))
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
        output_dir, 'peaks', '{}.peaks'.format(reference_stem)))
    for stem in input_stems:
        # Check if file has already been processed.
        in_path = os.path.join(output_dir, 'peaks', "{}.peaks".format(stem))
        out_path = os.path.join(output_dir, 'warped_peaks',
                                "{}.peaks".format(stem))
        out_path_tmap = os.path.join(output_dir, 'time_map',
                                     "{}.tmap".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Warping peaks: {}".format(stem))
        if stem == reference_stem:
            tapp.write_peaks(reference_peaks, out_path)
            time_map = tapp.calculate_time_map(
                reference_peaks, reference_peaks,
                tapp_parameters['warp2d_slack'],
                tapp_parameters['warp2d_window_size'],
                tapp_parameters['warp2d_num_points'],
                tapp_parameters['warp2d_rt_expand_factor'],
                tapp_parameters['warp2d_peaks_per_window'])
            tapp.write_time_map(time_map, out_path_tmap)
        else:
            logger.info("Reading peaks from disk: {}".format(stem))
            peaks = tapp.read_peaks(in_path)
            time_map = tapp.calculate_time_map(
                reference_peaks, peaks,
                tapp_parameters['warp2d_slack'],
                tapp_parameters['warp2d_window_size'],
                tapp_parameters['warp2d_num_points'],
                tapp_parameters['warp2d_rt_expand_factor'],
                tapp_parameters['warp2d_peaks_per_window'])
            peaks = tapp.warp_peaks(peaks, time_map)
            tapp.write_peaks(peaks, out_path)
            tapp.write_time_map(time_map, out_path_tmap)
    logger.info('Finished peak warping to reference ({}) in {}'.format(
        reference_stem, datetime.timedelta(seconds=time.time()-time_start)))

    # Calculate similarity matrix after alignment, generate heatmap and save to
    # disk.
    out_path = os.path.join(output_dir, 'quality', 'similarity_warped')
    logger.info("Starting warped similarity matrix calculation")
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or override_existing:
        similarity_matrix = np.zeros(
            len(input_stems) ** 2).reshape(len(input_stems), len(input_stems))
        for i in range(0, len(input_stems)):
            stem_a = input_stems[i]
            logger.info("Reading peaks_a from disk: {}".format(stem_a))
            peaks_a = tapp.read_peaks(os.path.join(
                output_dir, 'warped_peaks', '{}.peaks'.format(stem_a)))
            for j in range(i, len(input_stems)):
                stem_b = input_stems[j]
                logger.info("Reading peaks_b from disk: {}".format(stem_b))
                peaks_b = tapp.read_peaks(os.path.join(
                    output_dir, 'warped_peaks', '{}.peaks'.format(stem_b)))
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
        fig.set_size_inches(7.08661, 7.08661)
        plt.savefig("{}.pdf".format(out_path), dpi=300)
        plt.savefig("{}.png".format(out_path), dpi=300)
        plt.close(fig)
        logger.info(
            "Saving similarity matrix plot: {}.pdf/png".format(out_path))
    logger.info('Finished warped similarity matrix calculation in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    logger.info("Starting quality control plotting")
    time_start = time.time()
    out_path_tic_bpc = os.path.join(output_dir, 'quality', 'tic_base_peak')
    out_path_rt_vs_delta = os.path.join(
        output_dir, 'quality', 'rt_vs_rt_delta')
    out_path_sigmas_density = os.path.join(
        output_dir, 'quality', 'density_sigma')
    if (not os.path.exists("{}.png".format(out_path_tic_bpc))
            or not os.path.exists("{}.png".format(out_path_rt_vs_delta))
            or not os.path.exists("{}.png".format(out_path_sigmas_density))
            or not os.path.exists("{}.pdf".format(out_path_tic_bpc))
            or not os.path.exists("{}.pdf".format(out_path_rt_vs_delta))
            or not os.path.exists("{}.pdf".format(out_path_sigmas_density))
            or override_existing):
        plt.ioff()

        fig_tic_bpc, axes = plt.subplots(2, 2, sharex=True)
        ax1, ax2 = axes[0]
        ax3, ax4 = axes[1]

        fig_rt_vs_delta, ax5 = plt.subplots(1, 1)

        fig_sigmas_density, (ax6, ax7) = plt.subplots(1, 2)

        alpha = 0.5
        for i, stem in enumerate(input_stems):
            in_path_raw_data = os.path.join(
                output_dir, 'raw', "{}.ms1".format(stem))
            in_path_tmap = os.path.join(
                output_dir, 'time_map', "{}.tmap".format(stem))
            in_path_peaks = os.path.join(
                output_dir, 'warped_peaks', "{}.peaks".format(stem))
            logger.info("Reading raw_data from disk: {}".format(stem))
            raw_data = tapp.read_raw_data(in_path_raw_data)
            logger.info("Reading tmap from disk: {}".format(stem))
            tmap = tapp.read_time_map(in_path_tmap)
            logger.info("Reading peaks from disk: {}".format(stem))
            peaks = tapp.read_peaks(in_path_peaks)

            # Plot the unwarped TIC/Base peak.
            logger.info("Plotting unwarped TIC/Base peak: {}".format(stem))
            xic = tapp.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "sum"
            )
            x = xic.retention_time
            y = xic.intensity
            x_warped = [tmap.warp(rt) for rt in x]
            ax1.plot(x, y, label=stem, alpha=alpha)
            ax3.plot(x_warped, y, label=stem, alpha=alpha)

            # Plot the warped TIC/Base peak.
            logger.info("Plotting warped TIC/Base peak: {}".format(stem))
            xic = tapp.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "max"
            )
            x = xic.retention_time
            y = xic.intensity
            x_warped = [tmap.warp(rt) for rt in x]
            ax2.plot(x, y, label=stem, alpha=alpha)
            ax4.plot(x_warped, y, label=stem, alpha=alpha)

            # Plot the warping markers for each window.
            logger.info("Plotting warping markers: {}".format(stem))
            markers_not_warped = tmap.rt_start[1:]
            markers_warped = tmap.sample_rt_start[1:]
            ax1.scatter(markers_not_warped, np.repeat(
                0, len(markers_not_warped)), alpha=alpha, marker='^', edgecolor='none')
            ax2.scatter(markers_not_warped, np.repeat(
                0, len(markers_not_warped)), alpha=alpha, marker='^', edgecolor='none')
            ax3.scatter(markers_warped, np.repeat(
                0, len(markers_warped)), alpha=alpha, marker='^', edgecolor='none')
            ax4.scatter(markers_warped, np.repeat(
                0, len(markers_warped)), alpha=alpha, marker='^', edgecolor='none')

            # Plot rt vs delta.
            logger.info("Plotting rt vs rt_delta: {}".format(stem))
            rts = np.array([peak.fitted_rt for peak in peaks])
            rt_deltas = np.array([peak.rt_delta for peak in peaks])
            idx = np.argsort(rts)
            rts = rts[idx]
            rt_deltas = rt_deltas[idx]
            ax5.plot(rts, rt_deltas, label=stem, alpha=alpha)

            # Plot sigma mz/rt ditributions.
            logger.info(
                "Plotting density of sigma_mz/sigma_rt: {}".format(stem))
            sigma_mzs = np.array([peak.fitted_sigma_mz for peak in peaks])
            sigma_rts = np.array([peak.fitted_sigma_rt for peak in peaks])
            sns.kdeplot(sigma_rts, ax=ax6)
            sns.kdeplot(sigma_mzs, ax=ax7)

        logger.info("Saving figures to disk")

        # Save TIC/Base peak figure.
        ax1.set_title('Total Ion Chromatogram (TIC)')
        ax2.set_title('Base Peak Chromatogram')
        ax2.yaxis.set_label_position("right")
        ax4.yaxis.set_label_position("right")
        ax2.set_ylabel('Unwarped', rotation=-90, labelpad=20)
        ax4.set_ylabel('Warped', rotation=-90, labelpad=20)
        fig_tic_bpc.text(0.5, 0.04, 'Retention time (s)', ha='center')
        fig_tic_bpc.text(0.04, 0.5, 'Intensity',
                         va='center', rotation='vertical')
        handles, labels = ax4.get_legend_handles_labels()
        fig_tic_bpc.legend(handles, labels, loc='upper right')
        fig_tic_bpc.set_size_inches(7.08661, 7.08661/1.618034)
        plt.figure(fig_tic_bpc.number)
        plt.savefig("{}.pdf".format(out_path_tic_bpc), dpi=300)
        plt.savefig("{}.png".format(out_path_tic_bpc), dpi=300)
        plt.close(fig_tic_bpc)

        # Save rt vs rt_delta figure.
        ax5.set_xlabel('Retention time (s)')
        ax5.set_ylabel('Retention time delta (s)')
        plt.figure(fig_rt_vs_delta.number)
        fig_rt_vs_delta.legend(handles, labels, loc='upper right')
        fig_rt_vs_delta.set_size_inches(7.08661, 7.08661/1.618034)
        plt.savefig("{}.pdf".format(out_path_rt_vs_delta), dpi=300)
        plt.savefig("{}.png".format(out_path_rt_vs_delta), dpi=300)
        plt.close(fig_rt_vs_delta)

        # Save sigma density figure.
        ax6.set_xlabel('$\\sigma_{rt}$')
        ax7.set_xlabel('$\\sigma_{mz}$')
        ax6.set_ylabel('Density')
        plt.figure(fig_sigmas_density.number)
        fig_sigmas_density.set_size_inches(7.08661, 7.08661/1.618034)
        plt.savefig("{}.pdf".format(out_path_sigmas_density), dpi=300)
        plt.savefig("{}.png".format(out_path_sigmas_density), dpi=300)
        plt.close(fig_sigmas_density)

    logger.info('Finished quality control plotting in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform feature detection using averagine
    logger.info('Starting feature detection')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.peaks".format(stem))
        out_path = os.path.join(output_dir, 'features',
                                "{}.features".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Performing feature_detection: {}".format(stem))
        features = tapp.detect_features(
            peaks, tapp_parameters['feature_detection_charge_states'])
        logger.info('Writing features: {}'.format(out_path))
        tapp.write_features(features, out_path)

    logger.info('Finished feature detection in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Read mzidentdata and save binary data to disk.
    logger.info('Starting mzIdentML parsing')
    time_start = time.time()
    for i, stem in enumerate(input_stems):
        in_path = input_ident_files[i]
        out_path = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        if in_path == 'none' or (os.path.exists(out_path) and not override_existing):
            continue
        logger.info('Reading mzIdentML: {}'.format(in_path))
        ident_data = tapp.read_mzidentml(
            in_path,
            ignore_decoy=tapp_parameters['ident_ignore_decoy'],
            require_threshold=tapp_parameters['ident_require_threshold'],
            max_rank_only=tapp_parameters['ident_max_rank_only'],
            min_mz=tapp_parameters['min_mz'],
            max_mz=tapp_parameters['max_mz'],
            min_rt=tapp_parameters['min_rt'],
            max_rt=tapp_parameters['max_rt'],
        )
        logger.info('Writing ident data: {}'.format(out_path))
        tapp.write_ident_data(ident_data, out_path)
    logger.info('Finished mzIdentML parsing in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Match ms2 events with corresponding detected peaks.
    logger.info('Starting peaks/msms linkage')
    time_start = time.time()
    for stem in input_stems:
        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.peaks".format(stem))
        out_path = os.path.join(output_dir, 'linking',
                                "{}.ms2_peak.link".format(stem))
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

    # Link ms2 events with ident information.
    logger.info('Starting ident/msms linkage')
    time_start = time.time()
    for i, stem in enumerate(input_stems):
        # Check that we had identification info.
        if input_ident_files[i] == 'none':
            continue
        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.peaks".format(stem))
        in_path_idents = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        out_path = os.path.join(output_dir, 'linking',
                                "{}.ident_ms2.link".format(stem))
        out_path_psm = os.path.join(output_dir, 'linking',
                                    "{}.ident_peak.link".format(stem))
        if os.path.exists(out_path) and not override_existing:
            continue

        logger.info("Reading raw_data from disk (MS2): {}".format(stem))
        raw_data = tapp.read_raw_data(in_path_raw)

        logger.info("Reading ident from disk: {}".format(stem))
        ident_data = tapp.read_ident_data(in_path_idents)

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Performing linkage: {}".format(stem))
        linked_idents = tapp.link_idents(ident_data, raw_data)
        logger.info('Writing linked_msms: {}'.format(out_path))
        tapp.write_linked_msms(linked_idents, out_path)
        logger.info("Performing psm linkage: {}".format(stem))
        linked_psm = tapp.link_psm(ident_data, peaks, raw_data)
        logger.info('Writing linked_psm: {}'.format(out_path))
        tapp.write_linked_psm(linked_psm, out_path_psm)

    logger.info('Finished ident/msms linkage in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform metamatch on detected peaks.
    logger.info('Starting metamatch on peaks')
    time_start = time.time()
    in_path_peaks = os.path.join(output_dir, 'warped_peaks')
    out_path = os.path.join(output_dir, 'metamatch', "peaks.clusters")
    if (not os.path.exists(out_path) or override_existing):
        logger.info("Reading peaks from disk")
        peaks = [
            tapp.read_peaks(
                os.path.join(in_path_peaks, "{}.peaks".format(input_stem)))
            for input_stem in input_stems]

        logger.info("Finding peak clusters")
        peak_clusters = tapp.find_peak_clusters(
            groups,
            peaks,
            tapp_parameters["metamatch_fraction"],
            tapp_parameters["metamatch_n_sig_mz"],
            tapp_parameters["metamatch_n_sig_rt"])

        logger.info("Writing peak clusters to disk")
        tapp.write_peak_clusters(peak_clusters, out_path)

    logger.info('Finished metamatch on peaks in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Perform metamatch on detected features.
    logger.info('Starting metamatch on features')
    time_start = time.time()
    in_path_features = os.path.join(output_dir, 'features')
    out_path = os.path.join(output_dir, 'metamatch', "features.clusters")
    if (not os.path.exists(out_path) or override_existing):
        logger.info("Reading features from disk")
        features = [
            tapp.read_features(
                os.path.join(in_path_features, "{}.features".format(input_stem)))
            for input_stem in input_stems]

        logger.info("Finding feature clusters")
        feature_clusters = tapp.find_feature_clusters(
            groups,
            features,
            tapp_parameters["metamatch_fraction"],
            tapp_parameters["metamatch_n_sig_mz"],
            tapp_parameters["metamatch_n_sig_rt"])

        logger.info("Writing feature clusters to disk")
        tapp.write_feature_clusters(feature_clusters, out_path)

    logger.info('Finished metamatch on features in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    # Create final quantitative tables.
    logger.info('Starting creation of quantitative tables')
    time_start = time.time()
    for stem in input_stems:
        # Peak quantification.
        # ====================
        in_path_peaks = os.path.join(
            output_dir, 'warped_peaks', "{}.peaks".format(stem))
        in_path_peaks_link = os.path.join(
            output_dir, 'linking', "{}.ms2_peak.link".format(stem))
        in_path_ident_link_msms = os.path.join(
            output_dir, 'linking', "{}.ident_ms2.link".format(stem))
        in_path_ident_link_theomz = os.path.join(
            output_dir, 'linking', "{}.ident_peak.link".format(stem))
        in_path_ident_data = os.path.join(
            output_dir, 'ident', "{}.ident".format(stem))
        out_path_peaks = os.path.join(output_dir, 'quant',
                                      "{}_peaks.csv".format(stem))
        out_path_peak_annotations = os.path.join(output_dir, 'quant',
                                                 "{}_peak_annotations.csv".format(stem))

        # TODO: This is probably not necessary or needs to be changed if we are
        # doing all per-peak quantification in a single loop.
        if os.path.exists(out_path_peaks) and not override_existing:
            continue

        logger.info("Reading peaks from disk: {}".format(stem))
        peaks = tapp.read_peaks(in_path_peaks)

        logger.info("Generating peaks quantitative table")
        peaks_df = pd.DataFrame({
            'peak_id': [peak.id for peak in peaks],
            'mz': [peak.fitted_mz for peak in peaks],
            'rt': [peak.fitted_rt for peak in peaks],
            'rt_delta': [peak.rt_delta for peak in peaks],
            'height': [peak.fitted_height for peak in peaks],
            'sigma_mz': [peak.fitted_sigma_mz for peak in peaks],
            'sigma_rt': [peak.fitted_sigma_rt for peak in peaks],
            'volume': [peak.fitted_volume for peak in peaks],
            'smooth_height': [peak.local_max_height for peak in peaks],
            'smooth_mz': [peak.local_max_mz for peak in peaks],
            'smooth_rt': [peak.local_max_rt for peak in peaks],
            'roi_min_mz': [peak.roi_min_mz for peak in peaks],
            'roi_max_mz': [peak.roi_max_mz for peak in peaks],
            'roi_min_rt': [peak.roi_min_rt for peak in peaks],
            'roi_max_rt': [peak.roi_max_rt for peak in peaks],
            'raw_mean_mz': [peak.raw_roi_mean_mz for peak in peaks],
            'raw_mean_rt': [peak.raw_roi_mean_rt for peak in peaks],
            'raw_std_mz': [peak.raw_roi_sigma_mz for peak in peaks],
            'raw_std_rt': [peak.raw_roi_sigma_rt for peak in peaks],
            'raw_skewness_mz': [peak.raw_roi_skewness_mz for peak in peaks],
            'raw_skewness_rt': [peak.raw_roi_skewness_rt for peak in peaks],
            'raw_kurtosis_mz': [peak.raw_roi_kurtosis_mz for peak in peaks],
            'raw_kurtosis_rt': [peak.raw_roi_kurtosis_rt for peak in peaks],
            'raw_total_intensity': [peak.raw_roi_total_intensity for peak in peaks],
            'num_points': [peak.raw_roi_num_points for peak in peaks],
            'num_scans': [peak.raw_roi_num_scans for peak in peaks],
        })

        # Peak Annotations.
        # =================
        logger.info("Reading linked peaks from disk: {}".format(stem))
        peak_annotations = peaks_df[["peak_id"]]
        linked_peaks = tapp.read_linked_msms(in_path_peaks_link)
        linked_peaks = pd.DataFrame({
            'peak_id': [linked_peak.entity_id for linked_peak in linked_peaks],
            'msms_id': [linked_peak.msms_id for linked_peak in linked_peaks],
        })
        peak_annotations = pd.merge(
            peak_annotations, linked_peaks, on="peak_id", how="left")

        if os.path.isfile(in_path_ident_data):
            logger.info("Reading ident_data from disk: {}".format(stem))
            ident_data = tapp.read_ident_data(in_path_ident_data)
            psms = pd.DataFrame({
                'psm_index': [i for i in range(0, len(ident_data.spectrum_matches))],
                'psm_id': [psm.id for psm in ident_data.spectrum_matches],
                'psm_pass_threshold': [psm.pass_threshold for psm in ident_data.spectrum_matches],
                'psm_charge_state': [psm.charge_state for psm in ident_data.spectrum_matches],
                'psm_theoretical_mz': [psm.theoretical_mz for psm in ident_data.spectrum_matches],
                'psm_experimental_mz': [psm.experimental_mz for psm in ident_data.spectrum_matches],
                'psm_retention_time': [psm.retention_time for psm in ident_data.spectrum_matches],
                'psm_rank': [psm.rank for psm in ident_data.spectrum_matches],
                'psm_peptide_id': [psm.match_id for psm in ident_data.spectrum_matches],
            })
            if not psms.empty:
                if tapp_parameters["quant_ident_linkage"] == 'theoretical_mz':
                    logger.info(
                        "Reading linked ident_peak from disk: {}".format(stem))
                    linked_idents = tapp.read_linked_psm(
                        in_path_ident_link_theomz)
                    linked_idents = pd.DataFrame({
                        'peak_id': [linked_ident.peak_id for linked_ident in linked_idents],
                        'psm_index': [linked_ident.psm_index for linked_ident in linked_idents],
                        'psm_link_distance': [linked_ident.distance for linked_ident in linked_idents],
                    })
                    linked_idents = pd.merge(
                        linked_idents, psms, on="psm_index")
                    peak_annotations = pd.merge(
                        peak_annotations, linked_idents, on="peak_id", how="left")
                elif tapp_parameters["quant_ident_linkage"] == 'msms_event':
                    logger.info(
                        "Reading linked ident_peak from disk: {}".format(stem))
                    linked_idents = tapp.read_linked_msms(
                        in_path_ident_link_msms)
                    linked_idents = pd.DataFrame({
                        'msms_id': [linked_ident.msms_id for linked_ident in linked_idents],
                        'psm_index': [linked_ident.entity_id for linked_ident in linked_idents],
                        'psm_link_distance': [linked_ident.distance for linked_ident in linked_idents],
                    })
                    linked_idents = pd.merge(
                        linked_idents, psms, on="psm_index")
                    peak_annotations = pd.merge(
                        peak_annotations, linked_idents, on="msms_id", how="left")
                else:
                    raise ValueError("unknown quant_ident_linkage parameter")

                # Get the peptide information per psm.
                def format_modification(mod):
                    ret = "monoisotopic_mass_delta: {}, ".format(
                        mod.monoisotopic_mass_delta)
                    ret += "average_mass_delta: {}, ".format(
                        mod.average_mass_delta)
                    ret += "residues: {}, ".format(mod.residues)
                    ret += "location: {}, ".format(mod.location)
                    ret += "id: {}".format("; ".join(mod.id))
                    return ret
                peptides = pd.DataFrame({
                    'psm_peptide_id': [pep.id for pep in ident_data.peptides],
                    'psm_sequence': [pep.sequence for pep in ident_data.peptides],
                    'psm_modifications_num': [len(pep.modifications) for pep in ident_data.peptides],
                    'psm_modifications_info': [" / ".join(map(format_modification, pep.modifications)) for pep in ident_data.peptides],
                })
                peak_annotations = pd.merge(
                    peak_annotations, peptides, on="psm_peptide_id", how="left")

                # Get the protein information per peptide.
                db_sequences = pd.DataFrame({
                    'db_seq_id': [db_seq.id for db_seq in ident_data.db_sequences],
                    'protein_name': [db_seq.accession for db_seq in ident_data.db_sequences],
                    'protein_description': [db_seq.description for db_seq in ident_data.db_sequences],
                })
                peptide_evidence = pd.DataFrame({
                    'db_seq_id': [pe.db_sequence_id for pe in ident_data.peptide_evidence],
                    'psm_peptide_id': [pe.peptide_id for pe in ident_data.peptide_evidence],
                    'psm_decoy': [pe.decoy for pe in ident_data.peptide_evidence],
                })
                peptide_evidence = pd.merge(
                    peptide_evidence, db_sequences, on="db_seq_id").drop(["db_seq_id"], axis=1)

                # Get the protein information per psm.
                peak_annotations = pd.merge(
                    peak_annotations, peptide_evidence, on="psm_peptide_id")

        logger.info("Saving peaks quantitative table to disk: {}".format(stem))
        peaks_df.to_csv(out_path_peaks, index=False)

        logger.info("Saving peaks annotations table to disk: {}".format(stem))
        if "msms_id" in peak_annotations:
            peak_annotations["msms_id"] = peak_annotations["msms_id"].astype(
                'Int64')
        if "psm_charge_state" in peak_annotations:
            peak_annotations["psm_charge_state"] = peak_annotations["psm_charge_state"].astype(
                'Int64')
        if "psm_rank" in peak_annotations:
            peak_annotations["psm_rank"] = peak_annotations["psm_rank"].astype(
                'Int64')
        if "psm_modifications_num" in peak_annotations:
            peak_annotations["psm_modifications_num"] = peak_annotations["psm_modifications_num"].astype(
                'Int64')
        peak_annotations = peak_annotations.sort_values("peak_id")
        peak_annotations.to_csv(out_path_peak_annotations, index=False)

        in_path_features = os.path.join(
            output_dir, 'features', "{}.features".format(stem))
        if os.path.isfile(in_path_features):
            out_path_features = os.path.join(output_dir, 'quant',
                                             "{}_features.csv".format(stem))
            out_path_feature_annotations = os.path.join(output_dir, 'quant',
                                                        "{}_feature_annotations.csv".format(stem))

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
                'charge_state': [feature.charge_state for feature in features],
                'peak_id': [feature.peak_ids for feature in features],
            })
            # Find the peak annotations that belong to each feature.
            feature_annotations = features_df[[
                "feature_id", "peak_id"]].explode("peak_id")

            # TODO: It's possible that we want to regenerate the feature table
            # without doing the same with the peak tables.
            feature_annotations = pd.merge(
                feature_annotations, peak_annotations,
                on="peak_id", how="left")

            features_df.to_csv(out_path_features, index=False)
            feature_annotations.to_csv(
                out_path_feature_annotations, index=False)

    # Matched Peaks
    # =============
    logger.info("Reading peak clusters from disk")
    in_path_peak_clusters = os.path.join(
        output_dir, 'metamatch', 'peaks.clusters')
    out_path_peak_clusters_metadata = os.path.join(output_dir, 'quant',
                                                   "peak_clusters_metadata.csv")
    out_path_peak_clusters_peaks = os.path.join(output_dir, 'quant',
                                                "peak_clusters_peaks.csv")
    out_path_peak_clusters_annotations = os.path.join(output_dir, 'quant',
                                                      "peak_clusters_annotations.csv")

    def aggregate_cluster_annotations(x):
        ret = {}
        if "psm_sequence" in x:
            ret["psm_sequence"] = ".|.".join(
                np.unique(x['psm_sequence'].dropna())).strip(".|.")
        if "psm_charge_state" in x:
            ret["psm_charge_state"] = ".|.".join(
                map(str, np.unique(x['psm_charge_state'].dropna()))).strip(".|.")
        if "psm_modifications_num" in x:
            ret["psm_modifications_num"] = ".|.".join(
                map(str, np.unique(x['psm_modifications_num'].dropna()))).strip(".|.")
        if "protein_name" in x:
            ret["protein_name"] = ".|.".join(
                np.unique(x['protein_name'].dropna())).strip(".|.")
        if "protein_description" in x:
            ret["protein_description"] = ".|.".join(
                np.unique(x['protein_description'].dropna())).strip(".|.")
        if "consensus_sequence" in x:
            ret["consensus_sequence"] = ".|.".join(
                np.unique(x['consensus_sequence'].dropna())).strip(".|.")
        if "consensus_count" in x:
            ret["consensus_count"] = ".|.".join(
                map(str, np.unique(x['consensus_count'].dropna()))).strip(".|.")
        if "consensus_protein_name" in x:
            ret["consensus_protein_name"] = ".|.".join(
                np.unique(x['consensus_protein_name'].dropna())).strip(".|.")
        if "consensus_protein_description" in x:
            ret["consensus_protein_description"] = ".|.".join(
                np.unique(x['consensus_protein_description'].dropna())).strip(".|.")
        if "protein_group" in x:
            ret["protein_group"] = ".|.".join(
                map(str, np.unique(x['protein_group'].dropna()))).strip(".|.")
        return pd.Series(ret)

    if (not os.path.exists(out_path_peak_clusters_metadata) or override_existing):
        peak_clusters = tapp.read_peak_clusters(in_path_peak_clusters)
        logger.info("Generating peak clusters quantitative table")
        peak_clusters_metadata_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in peak_clusters],
            'mz': [cluster.mz for cluster in peak_clusters],
            'rt': [cluster.rt for cluster in peak_clusters],
            'avg_height': [cluster.avg_height for cluster in peak_clusters],
        })

        logger.info("Generating peak clusters quantitative table")
        peak_clusters_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in peak_clusters],
        })
        if tapp_parameters['quant_isotopes'] == 'volume':
            out_path_peak_clusters = os.path.join(output_dir, 'quant',
                                                  "peak_clusters_volume.csv")
            for i, stem in enumerate(input_stems):
                peak_clusters_df[stem] = [cluster.file_volumes[i]
                                          for cluster in peak_clusters]
        elif tapp_parameters['quant_isotopes'] == 'height':
            out_path_peak_clusters = os.path.join(output_dir, 'quant',
                                                  "peak_clusters_height.csv")
            for i, stem in enumerate(input_stems):
                peak_clusters_df[stem] = [cluster.heights[i]
                                          for cluster in peak_clusters]
        else:
            raise ValueError("unknown quant_isotopes parameter")
        logger.info("Writing peaks quantitative table to disk")
        peak_clusters_df.to_csv(out_path_peak_clusters, index=False)

        # Peak associations.
        logger.info("Generating peak clusters peak associations table")
        cluster_peaks = [(cluster.id, cluster.peak_ids)
                         for cluster in peak_clusters]
        cluster_peaks = pd.DataFrame(cluster_peaks, columns=[
                                     "cluster_id", "peak_ids"]).explode("peak_ids")
        cluster_peaks["file_id"] = cluster_peaks["peak_ids"].map(
            lambda x: input_stems[x.file_id])
        cluster_peaks["peak_id"] = cluster_peaks["peak_ids"].map(
            lambda x: x.peak_id)
        cluster_peaks = cluster_peaks.drop(["peak_ids"], axis=1)
        logger.info("Writing cluster to peak table to disk")
        cluster_peaks.to_csv(out_path_peak_clusters_peaks, index=False)

        # Cluster annotations.
        logger.info("Generating peak clusters annotations table")
        annotations = pd.DataFrame()
        for stem in input_stems:
            logger.info("Reading peak annotations for: {}".format(stem))
            in_path_peak_annotations = os.path.join(output_dir, 'quant',
                                                    "{}_peak_annotations.csv".format(stem))
            peak_annotations = pd.read_csv(
                in_path_peak_annotations, low_memory=False)
            cluster_annotations = cluster_peaks[cluster_peaks["file_id"] == stem][[
                "cluster_id", "peak_id"]]
            cluster_annotations["file_id"] = stem
            logger.info(
                "Merging peak/clusters annotations for: {}".format(stem))
            cluster_annotations = pd.merge(
                cluster_annotations, peak_annotations, on="peak_id", how="left")
            annotations = pd.concat(
                [annotations, cluster_annotations]).reset_index(drop=True)
        # Ensure these columns have the proper type.
        if "msms_id" in annotations:
            annotations["msms_id"] = annotations["msms_id"].astype(
                'Int64')
        if "psm_charge_state" in annotations:
            annotations["psm_charge_state"] = annotations["psm_charge_state"].astype(
                'Int64')
        if "psm_rank" in annotations:
            annotations["psm_rank"] = annotations["psm_rank"].astype(
                'Int64')
        if "psm_modifications_num" in annotations:
            annotations["psm_modifications_num"] = annotations["psm_modifications_num"].astype(
                'Int64')

        if tapp_parameters['quant_consensus'] and 'psm_sequence' in annotations:
            # Find a sequence consensus
            consensus_sequence = find_sequence_consensus(
                annotations, 'psm_sequence', tapp_parameters['quant_consensus_min_ident'])
            annotations = pd.merge(
                annotations,
                consensus_sequence[[
                    "cluster_id",
                    "consensus_sequence",
                    "consensus_count",
                ]], on="cluster_id", how="left")
            # Find a consensus proteins
            proteins = annotations[annotations['psm_sequence']
                                   == annotations['consensus_sequence']]
            proteins = proteins[['cluster_id', 'protein_name',
                                 'protein_description']].drop_duplicates()
            proteins.columns = [
                'cluster_id', 'consensus_protein_name', 'consensus_protein_description']
            annotations = pd.merge(annotations, proteins,
                                   on="cluster_id", how="left")

        # Saving annotations before aggregation.
        if tapp_parameters['quant_save_all_annotations']:
            logger.info("Writing annotations to disk")
            annotations = annotations.sort_values(by=["cluster_id"])
            annotations.to_csv(out_path_peak_clusters_annotations, index=False)

        logger.info("Aggregating annotations")
        annotations = annotations.groupby(
            'cluster_id').apply(aggregate_cluster_annotations)

        # Metadata
        logger.info("Merging metadata with annotations")
        peak_clusters_metadata_df = pd.merge(
            peak_clusters_metadata_df, annotations, how="left", on="cluster_id")
        logger.info("Writing metadata to disk")
        peak_clusters_metadata_df.to_csv(
            out_path_peak_clusters_metadata, index=False)

    # Matched Features
    # ================
    logger.info("Reading feature clusters from disk")
    in_path_feature_clusters = os.path.join(
        output_dir, 'metamatch', 'features.clusters')
    out_path_feature_clusters_metadata = os.path.join(output_dir, 'quant',
                                                      "feature_clusters_metadata.csv")
    out_path_feature_clusters_features = os.path.join(output_dir, 'quant',
                                                      "feature_clusters_features.csv")
    out_path_feature_clusters_annotations = os.path.join(output_dir, 'quant',
                                                         "feature_clusters_annotations.csv")
    if (not os.path.exists(out_path_feature_clusters_metadata) or override_existing):
        feature_clusters = tapp.read_feature_clusters(
            in_path_feature_clusters)

        logger.info("Generating feature clusters quantitative table")
        metadata = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in feature_clusters],
            'mz': [cluster.mz for cluster in feature_clusters],
            'rt': [cluster.rt for cluster in feature_clusters],
            'avg_height': [cluster.avg_total_height for cluster in feature_clusters],
            'charge_state': [cluster.charge_state for cluster in feature_clusters],
        })
        data = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in feature_clusters],
        })
        if tapp_parameters['quant_features'] == 'monoisotopic_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_monoisotopic_height.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.monoisotopic_heights[i]
                                             for cluster in feature_clusters]
        elif tapp_parameters['quant_features'] == 'monoisotopic_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_monoisotopic_volume.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.monoisotopic_volumes[i]
                                             for cluster in feature_clusters]
        elif tapp_parameters['quant_features'] == 'total_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_total_height.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.total_heights[i]
                                             for cluster in feature_clusters]
        elif tapp_parameters['quant_features'] == 'total_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_total_volume.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.total_volumes[i]
                                             for cluster in feature_clusters]
        elif tapp_parameters['quant_features'] == 'max_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_max_height.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.max_heights[i]
                                             for cluster in feature_clusters]
        elif tapp_parameters['quant_features'] == 'max_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_max_volume.csv")
            for i, stem in enumerate(input_stems):
                data[stem] = [cluster.max_volumes[i]
                                             for cluster in feature_clusters]
        else:
            raise ValueError("unknown quant_features parameter")
        logger.info("Writing feature clusters quantitative table to disk")
        data.to_csv(out_path_feature_clusters, index=False)

        # Feature associations.
        logger.info("Generating feature clusters feature associations table")
        cluster_features = [(cluster.id, cluster.feature_ids)
                            for cluster in feature_clusters]
        cluster_features = pd.DataFrame(cluster_features, columns=[
                                        "cluster_id", "feature_ids"]).explode("feature_ids")
        cluster_features["file_id"] = cluster_features["feature_ids"].map(
            lambda x: input_stems[x.file_id])
        cluster_features["feature_id"] = cluster_features["feature_ids"].map(
            lambda x: x.feature_id)
        cluster_features = cluster_features.drop(["feature_ids"], axis=1)
        logger.info("Writing cluster to feature table to disk")
        cluster_features.to_csv(
            out_path_feature_clusters_features, index=False)

        # Cluster annotations.
        logger.info("Generating peak clusters annotations table")
        annotations = pd.DataFrame()
        for stem in input_stems:
            logger.info("Reading features for: {}".format(stem))
            in_path_peak_features = os.path.join(output_dir, 'features',
                                                 "{}.features".format(stem))
            in_path_peak_annotations = os.path.join(output_dir, 'quant',
                                                    "{}_peak_annotations.csv".format(stem))
            features = tapp.read_features(in_path_peak_features)
            features = [(feature.id, feature.peak_ids, feature.charge_state)
                        for feature in features]
            features = pd.DataFrame(
                features, columns=["feature_id", "peak_id", "charge_state"]).explode("peak_id")
            peak_annotations = pd.read_csv(
                in_path_peak_annotations, low_memory=False)
            cluster_annotations = cluster_features[cluster_features["file_id"] == stem][[
                "cluster_id", "feature_id"]]
            cluster_annotations["file_id"] = stem
            logger.info(
                "Merging peak/clusters annotations for: {}".format(stem))
            cluster_annotations = pd.merge(
                cluster_annotations, features, on="feature_id", how="left")
            cluster_annotations = pd.merge(
                cluster_annotations, peak_annotations, on="peak_id", how="left")
            annotations = pd.concat([annotations, cluster_annotations])

        # Ensure these columns have the proper type.
        if "msms_id" in annotations:
            annotations["msms_id"] = annotations["msms_id"].astype('Int64')
        if "charge_state" in annotations:
            annotations["charge_state"] = annotations["charge_state"].astype(
                'Int64')
        if "psm_charge_state" in annotations:
            annotations["psm_charge_state"] = annotations["psm_charge_state"].astype(
                'Int64')
        if "psm_rank" in annotations:
            annotations["psm_rank"] = annotations["psm_rank"].astype('Int64')
        if "psm_modifications_num" in annotations:
            annotations["psm_modifications_num"] = annotations["psm_modifications_num"].astype(
                'Int64')

        if tapp_parameters['quant_consensus'] and 'psm_sequence' in annotations:
            # Find a sequence consensus
            consensus_sequence = find_sequence_consensus(
                annotations, 'psm_sequence', tapp_parameters['quant_consensus_min_ident'])
            annotations = pd.merge(
                annotations,
                consensus_sequence[[
                    "cluster_id",
                    "consensus_sequence",
                    "consensus_count",
                ]], on="cluster_id", how="left")
            # Find a consensus proteins
            proteins = annotations[annotations['psm_sequence']
                                   == annotations['consensus_sequence']]
            proteins = proteins[['cluster_id', 'protein_name',
                                 'protein_description']].drop_duplicates()
            proteins.columns = [
                'cluster_id', 'consensus_protein_name', 'consensus_protein_description']
            annotations = pd.merge(annotations, proteins,
                                   on="cluster_id", how="left")

        # Calculate protein groups.
        logger.info("Calculating protein groups")
        sequence_column = 'psm_sequence'
        protein_name_column = 'protein_name'
        protein_description_column = 'protein_description'
        if tapp_parameters['quant_consensus'] and 'psm_sequence' in annotations:
            sequence_column = 'consensus_sequence'
            protein_name_column = 'consensus_protein_name'
            protein_description_column = 'consensus_protein_description'

        # Combine protein name/description in case
        if (sequence_column in annotations and
                protein_name_column in annotations and
                protein_description_column in annotations):
            annotations = find_protein_groups(
                    annotations,
                    sequence_column,
                    protein_name_column,
                    protein_description_column)

        # Saving annotations before aggregation.
        if tapp_parameters['quant_save_all_annotations']:
            logger.info("Writing annotations to disk")
            annotations = annotations.sort_values(by=["cluster_id"])
            annotations.to_csv(
                out_path_feature_clusters_annotations, index=False)

        logger.info("Aggregating annotations")
        if ("psm_charge_state" in annotations and
                tapp_parameters['quant_features_charge_state_filter']):
            annotations = annotations[annotations["psm_charge_state"]
                                      == annotations["charge_state"]]
        annotations_agg = annotations.groupby(
            'cluster_id').apply(aggregate_cluster_annotations)

        # Metadata.
        logger.info("Merging metadata with annotations")
        metadata = pd.merge(
            metadata, annotations_agg, how="left", on="cluster_id")
        logger.info("Writing metadata to disk")
        metadata.to_csv(out_path_feature_clusters_metadata, index=False)

        if 'protein_group' in metadata:
            logger.info("Aggregating protein groups")
            def aggregate_protein_group_annotations(x):
                ret = {}
                if "psm_sequence" in x:
                    ret["psm_sequence"] = ".|.".join(
                        np.unique(x['psm_sequence'].dropna())).strip(".|.")
                if "protein_name" in x:
                    ret["protein_name"] = ".|.".join(
                        np.unique(x['protein_name'].dropna())).strip(".|.")
                if "protein_description" in x:
                    ret["protein_description"] = ".|.".join(
                        np.unique(x['protein_description'].dropna())).strip(".|.")
                if "consensus_sequence" in x:
                    ret["consensus_sequence"] = ".|.".join(
                        np.unique(x['consensus_sequence'].dropna())).strip(".|.")
                if "consensus_protein_name" in x:
                    ret["consensus_protein_name"] = ".|.".join(
                        np.unique(x['consensus_protein_name'].dropna())).strip(".|.")
                if "consensus_protein_description" in x:
                    ret["consensus_protein_description"] = ".|.".join(
                        np.unique(x['consensus_protein_description'].dropna())).strip(".|.")
                return pd.Series(ret)
            prot_data = data.copy()
            prot_data = prot_data.drop(["cluster_id"], axis=1)
            prot_data['protein_group'] = pd.to_numeric(metadata['protein_group']).astype('Int64')
            prot_data = prot_data[~(prot_data['protein_group'].isna())]
            prot_data = prot_data.groupby('protein_group').agg(sum)
            prot_data = prot_data.reset_index()
            prot_metadata = annotations.copy()
            prot_metadata = prot_metadata[~prot_metadata['protein_group'].isna()]
            prot_metadata = prot_metadata.groupby('protein_group')
            prot_metadata = prot_metadata.apply(aggregate_protein_group_annotations)
            prot_metadata = prot_metadata.reset_index()
            out_path_protein_data = os.path.join(output_dir, 'quant',
                    "protein_groups.csv")
            out_path_protein_metadata = os.path.join(output_dir, 'quant',
                    "protein_groups_metadata.csv")

            logger.info("Writing protein group data/metadata to disk")
            prot_data.to_csv(out_path_protein_data, index=False)
            prot_metadata.to_csv(out_path_protein_metadata, index=False)

    logger.info('Finished creation of quantitative tables in {}'.format(
        datetime.timedelta(seconds=time.time()-time_start)))

    logger.info("Performing summary")
    dda_pipeline_summary(tapp_parameters, input_stems, output_dir)

    logger.info('Total time elapsed: {}'.format(
        datetime.timedelta(seconds=time.time()-time_pipeline_start)))

    # Stop logger.
    logger.removeHandler(logger_fh)
    logger_fh.close()
