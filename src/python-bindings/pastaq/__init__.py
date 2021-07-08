from .pastaq import *
import pastaq
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
    # NOTE: Should we ignore annotations when the consensus is ambiguous? For
    # example, if we have sequence A for 3 samples and sequence B for another 3,
    # the consensus is "unclear". Furthermore, we might choose to ignore it or
    # to choose a suitable one. For now we are keeping all annotations even if
    # that means having ambiguity about the sequence being assigned to a given
    # isotope/feature cluster.
    consensus = consensus[consensus["consensus_count"] >= min_consensus_count]
    consensus = consensus.drop(["consensus_count_max"], axis=1)
    return consensus

def find_protein_groups(
        feature_data,
        feature_annotations,
        sequence_col,
        protein_col,
        protein_description_col,
        min_peptides_number,
        remove_subset_proteins,
        ignore_ambiguous_peptides,
        protein_quant_type,
        ):


    # Create a list of annotations only with the info we need for protein groups.
    protein_annotations = feature_annotations[['cluster_id', sequence_col, protein_col, protein_description_col]].copy()
    protein_annotations = protein_annotations.drop_duplicates().dropna()

    # Assign a unique numeric id for peptide and proteins for faster comparisons.
    peptide_map = dict(zip(range(0, len(protein_annotations[sequence_col].unique())), protein_annotations[sequence_col].unique()))
    protein_map = dict(zip(range(0, len(protein_annotations[protein_col].unique())), protein_annotations[protein_col].unique()))
    peptide_map_df = pd.DataFrame({'peptide_id': peptide_map.keys(), sequence_col: peptide_map.values()})
    protein_map_df = pd.DataFrame({'protein_id': protein_map.keys(), protein_col: protein_map.values()})
    protein_annotations = pd.merge(protein_annotations, peptide_map_df)
    protein_annotations = pd.merge(protein_annotations, protein_map_df)
    protein_annotations = pd.merge(protein_annotations, protein_annotations)

    # Create a copy of the feature data for aggregation.
    protein_group_data = feature_data.copy()
    protein_group_data['protein_group_id'] = -1
    protein_group_metadata = protein_annotations.copy()
    protein_group_metadata['protein_group_id'] = -1

    # Initialize graph of protein-peptide nodes.
    protein_nodes = {} # protein_id -> [peptide_1, peptide_2...]
    peptide_nodes = {} # peptide_id -> [protein_1, protein_2...]
    for index, row in protein_annotations.iterrows():
        if row.protein_id not in protein_nodes:
            protein_nodes[row.protein_id] = set()
        protein_nodes[row.protein_id].add(row.peptide_id)
        if row.peptide_id not in peptide_nodes:
            peptide_nodes[row.peptide_id] = set()
        peptide_nodes[row.peptide_id].add(row.protein_id)

    def remove_proteins_from_graph(to_be_removed, protein_nodes, peptide_nodes):
        for protein_id in to_be_removed:
            for peptide_id in protein_nodes[protein_id]:
                peptide_nodes[peptide_id] = set([
                        this_id
                        for this_id in peptide_nodes[peptide_id]
                        if this_id != protein_id
                    ])
            if len(peptide_nodes[peptide_id]) == 0:
                del peptide_nodes[peptide_id]
            del protein_nodes[protein_id]

    # Filter out proteins that don't meet the minimum number of peptides requirement.
    if min_peptides_number > 1:
        to_be_removed = []
        for protein_id, peptide_ids in protein_nodes.items():
            if len(peptide_ids) < min_peptides_number:
                to_be_removed += [protein_id]
        remove_proteins_from_graph(to_be_removed, protein_nodes, peptide_nodes)

    # Remove fully subset proteins.
    if remove_subset_proteins:
        subset_proteins = []
        for protein_id, peptide_ids in protein_nodes.items():

            # Find proteins that share some peptide with this one.
            target_proteins = set()
            for peptide_id in peptide_ids:
                for target_protein in peptide_nodes[peptide_id]:
                    if target_protein != protein_id:
                        target_proteins.add(target_protein)

            # Check target proteins to check if peptides are fully contained within
            # another group.
            for target_protein in target_proteins:
                target_peptides = protein_nodes[target_protein]
                if set(peptide_ids).issubset(target_peptides) and len(peptide_ids) < len(target_peptides):
                    subset_proteins += [protein_id]
                    break
        remove_proteins_from_graph(subset_proteins, protein_nodes, peptide_nodes)

    #
    # Identify the type of peptide (Unique vs shared) and if razor principle is
    # applied, to which protein(s) will they be assigned.
    # - Unique: Only appears on a single protein.
    # - Shared: The peptide is shared by one or more proteins.
    # - Razor: Appears on multiple proteins, assigned to the protein with the
    #          largest number of peptides.
    # - Ambiguous: When multiple proteins have the exact same number of peptides.
    #              Has more than one assigned protein.
    #

    # Initialize protein_peptides
    protein_peptides = {}
    for protein_id in protein_nodes.keys():
        protein_peptides[protein_id] = {
                'unique': set(),
                'shared': set(),
                'razor': set(),
                'ambiguous': set(),
            }

    for peptide_id, protein_ids in peptide_nodes.items():
        if len(protein_ids) == 1:
            protein_id = list(protein_ids)[0]
            if protein_id not in protein_peptides:
                protein_peptides[protein_id] = {
                        'unique': [peptide_id],
                        'razor': [],
                        'ambiguous': [],
                }
            else:
                protein_peptides[protein_id]['unique'].add(peptide_id)
            continue

        cur_assigned_proteins = []
        cur_assigned_peptide_count = 0
        for protein_id in protein_ids:
            protein_peptides[protein_id]['shared'].add(peptide_id)
            if len(protein_nodes[protein_id]) == cur_assigned_peptide_count:
                cur_assigned_proteins += [protein_id]
            if len(protein_nodes[protein_id]) > cur_assigned_peptide_count:
                cur_assigned_proteins = [protein_id]
                cur_assigned_peptide_count = len(protein_nodes[protein_id])

        for protein_id in cur_assigned_proteins:
            if len(cur_assigned_proteins) > 1:
                protein_peptides[protein_id]['ambiguous'].add(peptide_id)
            else:
                protein_peptides[protein_id]['razor'].add(peptide_id)

    # Find protein groups that contain unique peptides
    protein_groups = {}
    protein_group_counter = 0
    unique_protein_ids = []
    non_unique_protein_ids = []
    for protein_id, peptide_ids in protein_nodes.items():
        unique_found = False
        for peptide_id in peptide_ids:
            if len(peptide_nodes[peptide_id]) == 1:
                unique_protein_ids += [protein_id]
                unique_found = True
                protein_groups[protein_group_counter] = set([protein_id])
                protein_group_counter += 1
                break
        if not unique_found:
                non_unique_protein_ids += [protein_id]

    # Remove unique protein groups from the graph.
    remove_proteins_from_graph(unique_protein_ids, protein_nodes, peptide_nodes)

    # Group proteins with shared peptides only.
    explored_proteins = set()
    for protein_id in non_unique_protein_ids:
        if protein_id in explored_proteins:
            continue

        previous_protein_group_size = 0
        protein_group = set()
        protein_group.add(protein_id)
        while previous_protein_group_size != len(protein_group):
            previous_protein_group_size = len(protein_group)
            proteins_to_explore = [id for id in protein_group if id not in explored_proteins]

            # Find all proteins associated with all peptides for unexplored proteins.
            for id in proteins_to_explore:
                for peptide_id in protein_nodes[id]:
                    for new_protein in peptide_nodes[peptide_id]:
                        protein_group.add(new_protein)
                explored_proteins.add(id)

        # Update protein group list and remove them from the available list.
        protein_groups[protein_group_counter] = protein_group
        protein_group_counter += 1
        remove_proteins_from_graph(list(protein_group), protein_nodes, peptide_nodes)

    # Depending on the selected quantification type, find which peptides should be
    # used for quantification for each protein group.
    for protein_group_id, protein_ids in protein_groups.items():
        selected_peptides = set()
        for protein_id in protein_ids:
            for peptide_id in protein_peptides[protein_id]['unique']:
                selected_peptides.add(peptide_id)
            if protein_quant_type == 'all':
                for peptide_id in protein_peptides[protein_id]['shared']:
                    selected_peptides.add(peptide_id)
            elif protein_quant_type == 'razor':
                for peptide_id in protein_peptides[protein_id]['razor']:
                    selected_peptides.add(peptide_id)
                if ignore_ambiguous_peptides:
                    for peptide_id in protein_peptides[protein_id]['ambiguous']:
                        selected_peptides.add(peptide_id)
        selected_cluster_ids = protein_annotations.cluster_id[protein_annotations.peptide_id.isin(selected_peptides)].unique()
        protein_group_data.loc[protein_group_data.cluster_id.isin(selected_cluster_ids), 'protein_group_id'] = protein_group_id
        protein_group_metadata.loc[protein_group_metadata.peptide_id.isin(selected_peptides) & protein_group_metadata.protein_id.isin(protein_ids) , 'protein_group_id'] = protein_group_id

    def aggregate_protein_group_annotations(x):
        ret = {}
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

    protein_group_data = protein_group_data[protein_group_data.protein_group_id != -1]
    del protein_group_data['cluster_id']
    protein_group_data = protein_group_data.groupby('protein_group_id').sum().reset_index()

    protein_group_metadata = protein_group_metadata[protein_group_metadata.protein_group_id != -1]
    del protein_group_metadata['cluster_id']
    del protein_group_metadata['peptide_id']
    del protein_group_metadata['protein_id']
    protein_group_metdata = protein_group_metadata.drop_duplicates()
    protein_group_metadata = protein_group_metadata.groupby('protein_group_id')
    protein_group_metadata = protein_group_metadata.apply(aggregate_protein_group_annotations)
    protein_group_metadata = protein_group_metadata.reset_index()

    return protein_group_data, protein_group_metadata

def plot_mesh(mesh, transform='sqrt', figure=None):
    plt.style.use('dark_background')

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

    plt.plot(x, y, label='peak_id = {}'.format(peak.id))
    plt.xlabel('Retention time (s)')
    plt.ylabel('Intensity')
    plt.legend()

    return figure


def plot_peak_raw_points(
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
        pastaq_parameters = {
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
            'warp2d_window_size': 100,
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
            # Annotation linking.
            #
            'link_n_sig_mz': 3,
            'link_n_sig_rt': 3,
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
            # Options: Any 'seaborn' supported palette style, like:
            #          'husl', 'crest', 'Spectral', 'flare', 'mako', etc.
            'qc_plot_palette': 'husl',
            # Options: 'png', 'pdf', 'eps'...
            'qc_plot_extension': 'png',
            # Options: 'dynamic', [0.0-1.0]
            'qc_plot_fill_alpha': 'dynamic',
            'qc_plot_line_alpha': 0.5,
            'qc_plot_scatter_alpha': 0.3,
            'qc_plot_scatter_size': 2,
            'qc_plot_min_dynamic_alpha': 0.1,
            'qc_plot_per_file': False,
            # Options: 'fill', 'line'
            'qc_plot_line_style': 'fill',
            # Plot style config.
            'qc_plot_dpi': 300,
            'qc_plot_font_family': 'sans-serif',
            'qc_plot_font_size': 7,
            'qc_plot_fig_size_x': 7.08661,
            'qc_plot_fig_size_y': 7.08661/1.618034,
            'qc_plot_fig_legend': False,
            'qc_plot_mz_vs_sigma_mz_max_peaks': 200000,
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
            'quant_ident_linkage': 'msms_event',
            # Whether to obtain a consensus sequence and proteins on identifications.
            'quant_consensus': True,
            # Demand a minimum number of files with identification per cluster.
            'quant_consensus_min_ident': 2,
            # Whether to store all the annotations prior to cluster aggregation.
            'quant_save_all_annotations': True,
            # Minimum number of peptides necessary to consider a protein for
            # quantification.
            'quant_proteins_min_peptides': 1,
            # Wether to remove proteins whose peptides are entirely contined
            # within another with a longest number of evidence peptides.
            'quant_proteins_remove_subset_proteins': True,
            # In case a peptide can't be assigned to a unique protein as 'razor'
            # we can choose to use them regardless in all instances where they
            # would if they were to be considered razor or to ignore them.
            'quant_proteins_ignore_ambiguous_peptides': True,
            # Protein inference method:
            #     - 'unique': Only unique peptides will be used for
            #                 quantification.
            #     - 'razor': Unique and peptides assigned as most likely through
            #                the Occam's razor constrain.
            #     - 'all': All peptides will be used for quantification for all
            #              protein groups. Thus shared peptides will be used more
            #              than once.
            'quant_proteins_quant_type': 'razor',
        }
        return pastaq_parameters


def _custom_log(msg, logger):
    if logger:
        logger.info(msg)
    print(msg)

def parse_raw_files(params, output_dir, logger=None, force_override=False):
    _custom_log('Starting raw data conversion', logger)
    time_start = time.time()

    for file in params['input_files']:
        raw_path = file['raw_path']
        stem = file['stem']

        # Check if file has already been processed.
        out_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        if os.path.exists(out_path) and not force_override:
            continue

        # File extension.
        file_extension = os.path.splitext(raw_path)[1]

        # Read raw files (MS1).
        _custom_log('Reading MS1: {}'.format(raw_path), logger)
        if file_extension.lower() == '.mzxml':
            raw_data = pastaq.read_mzxml(
                raw_path,
                min_mz=params['min_mz'],
                max_mz=params['max_mz'],
                min_rt=params['min_rt'],
                max_rt=params['max_rt'],
                instrument_type=params['instrument_type'],
                resolution_ms1=params['resolution_ms1'],
                resolution_msn=params['resolution_msn'],
                reference_mz=params['reference_mz'],
                fwhm_rt=params['avg_fwhm_rt'],
                polarity=params['polarity'],
                ms_level=1,
            )
        elif file_extension.lower() == '.mzml':
            raw_data = pastaq.read_mzml(
                raw_path,
                min_mz=params['min_mz'],
                max_mz=params['max_mz'],
                min_rt=params['min_rt'],
                max_rt=params['max_rt'],
                instrument_type=params['instrument_type'],
                resolution_ms1=params['resolution_ms1'],
                resolution_msn=params['resolution_msn'],
                reference_mz=params['reference_mz'],
                fwhm_rt=params['avg_fwhm_rt'],
                polarity=params['polarity'],
                ms_level=1,
            )

        # Write raw_data to disk (MS1).
        _custom_log('Writing MS1: {}'.format(out_path), logger)
        raw_data.dump(out_path)

    for file in params['input_files']:
        raw_path = file['raw_path']
        stem = file['stem']

        # Check if file has already been processed.
        out_path = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        if os.path.exists(out_path) and not force_override:
            continue

        # File extension.
        file_extension = os.path.splitext(raw_path)[1]

        # Read raw files (MS2).
        _custom_log('Reading MS2: {}'.format(raw_path), logger)
        if file_extension.lower() == '.mzxml':
            raw_data = pastaq.read_mzxml(
                raw_path,
                min_mz=params['min_mz'],
                max_mz=params['max_mz'],
                min_rt=params['min_rt'],
                max_rt=params['max_rt'],
                instrument_type=params['instrument_type'],
                resolution_ms1=params['resolution_ms1'],
                resolution_msn=params['resolution_msn'],
                reference_mz=params['reference_mz'],
                fwhm_rt=params['avg_fwhm_rt'],
                polarity=params['polarity'],
                ms_level=2,
            )
        elif file_extension.lower() == '.mzml':
            raw_data = pastaq.read_mzml(
                raw_path,
                min_mz=params['min_mz'],
                max_mz=params['max_mz'],
                min_rt=params['min_rt'],
                max_rt=params['max_rt'],
                instrument_type=params['instrument_type'],
                resolution_ms1=params['resolution_ms1'],
                resolution_msn=params['resolution_msn'],
                reference_mz=params['reference_mz'],
                fwhm_rt=params['avg_fwhm_rt'],
                polarity=params['polarity'],
                ms_level=2,
            )

        # Write raw_data to disk (MS2).
        _custom_log('Writing MS2: {}'.format(out_path), logger)
        raw_data.dump(out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished raw data parsing in {}'.format(elapsed_time), logger)

def detect_peaks(params, output_dir, save_grid=False, logger=None, force_override=False):
    # Perform resampling/smoothing and peak detection and save results to disk.
    _custom_log('Starting peak detection', logger)
    time_start = time.time()

    for file in params['input_files']:
        stem = file['stem']

        # Check if file has already been processed.
        in_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        out_path = os.path.join(output_dir, 'peaks', "{}.peaks".format(stem))
        if os.path.exists(out_path) and not force_override:
            continue

        _custom_log("Reading raw_data from disk: {}".format(stem), logger)
        raw_data = pastaq.read_raw_data(in_path)

        _custom_log("Resampling: {}".format(stem), logger)
        grid = pastaq.resample(
            raw_data,
            params['num_samples_mz'],
            params['num_samples_rt'],
            params['smoothing_coefficient_mz'],
            params['smoothing_coefficient_rt'],
        )

        if save_grid:
            mesh_path = os.path.join(output_dir, 'grid', "{}.grid".format(stem))
            _custom_log('Writing grid: {}'.format(mesh_path), logger)
            grid.dump(mesh_path)

        _custom_log("Finding peaks: {}".format(stem), logger)
        peaks = pastaq.find_peaks(raw_data, grid, params['max_peaks'])
        _custom_log('Writing peaks:'.format(out_path), logger)
        pastaq.write_peaks(peaks, out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished peak detection in {}'.format(elapsed_time), logger)

def calculate_similarity_matrix(params, output_dir, peak_dir, logger=None, force_override=False):
    out_path = os.path.join(output_dir, 'quality', 'similarity_{}.csv'.format(peak_dir))
    if os.path.exists(out_path) and not force_override:
        return

    _custom_log("Starting similarity matrix calculation for {}".format(peak_dir), logger)
    time_start = time.time()
    if not os.path.exists("{}.csv".format(out_path)) or force_override:
        input_files = params['input_files']
        n_files = len(input_files)
        similarity_matrix = np.zeros(n_files ** 2).reshape(n_files, n_files)
        for i in range(0, n_files):
            stem_a = input_files[i]['stem']
            peaks_a = pastaq.read_peaks(os.path.join(
                output_dir, peak_dir, '{}.peaks'.format(stem_a)))
            for j in range(i, n_files):
                stem_b = input_files[j]['stem']
                peaks_b = pastaq.read_peaks(os.path.join(output_dir, peak_dir, '{}.peaks'.format(stem_b)))
                _custom_log("Calculating similarity of {} vs {}".format(stem_a, stem_b), logger)
                similarity_matrix[j, i] = pastaq.find_similarity(
                    peaks_a, peaks_b,
                    params['similarity_num_peaks']).geometric_ratio
                similarity_matrix[i, j] = similarity_matrix[j, i]
        similarity_matrix = pd.DataFrame(similarity_matrix)
        similarity_matrix_names = [input_file['stem'] for input_file in input_files]
        similarity_matrix.columns = similarity_matrix_names
        similarity_matrix.rename(index=dict(zip(range(0, len(similarity_matrix_names), 1), similarity_matrix_names)), inplace=True)

        # Save similarity matrix to disk.
        _custom_log("Saving similarity matrix for {}: {}".format(peak_dir, out_path), logger)
        similarity_matrix.to_csv(out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished similarity matrix calculation from {} in {}'.format(peak_dir, elapsed_time), logger)

def perform_rt_alignment(params, output_dir, logger=None, force_override=False):
    warped_sim_path = os.path.join(output_dir, 'quality', 'similarity_warped_peaks.csv')
    if os.path.exists(warped_sim_path) and not force_override:
        return
    input_files = params['input_files']
    # Find selected reference samples.
    ref_candidates = []
    for input_file in input_files:
        if 'reference' in input_file and input_file['reference']:
            ref_candidates += [input_file]

    if len(ref_candidates) == 1:
        # If only a reference sample is marked, it will be used.
        ref = ref_candidates[0]
        _custom_log("Using selected reference: {}".format(ref['stem']), logger)
    else:
        # If no samples are selected, exhaustive search will be performed.
        if len(ref_candidates) == 0:
            _custom_log("No reference selected, performing exhaustive search", logger)
            ref_candidates = input_files

        # Find optimal reference sample from the list of candidates.
        _custom_log("Starting optimal reference search", logger)
        time_start = time.time()
        n_ref = len(ref_candidates)
        n_files = len(input_files)
        similarity_matrix = np.zeros(n_ref * n_files).reshape(n_ref, n_files)
        for i in range(0, n_ref):
            stem_a = ref_candidates[i]['stem']
            peaks_a = pastaq.read_peaks(os.path.join(
                output_dir, 'peaks', '{}.peaks'.format(stem_a)))
            for j in range(0, n_files):
                if i == j:
                    similarity_matrix[i, j] = 1
                    continue
                stem_b = input_files[j]['stem']
                peaks_b = pastaq.read_peaks(os.path.join(output_dir, 'peaks', '{}.peaks'.format(stem_b)))
                _custom_log("Warping {} peaks to {}".format(stem_b, stem_a), logger)
                time_map = pastaq.calculate_time_map(
                    peaks_a, peaks_b,
                    params['warp2d_slack'],
                    params['warp2d_window_size'],
                    params['warp2d_num_points'],
                    params['warp2d_rt_expand_factor'],
                    params['warp2d_peaks_per_window'])
                peaks_b = pastaq.warp_peaks(peaks_b, time_map)
                _custom_log("Calculating similarity of {} vs {} (warped)".format(stem_a, stem_b), logger)
                similarity_matrix[i, j] = pastaq.find_similarity(
                    peaks_a, peaks_b,
                    params['similarity_num_peaks']).geometric_ratio

        elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
        _custom_log('Finished optimal reference search in {}'.format(elapsed_time), logger)

        # Find the reference with maximum similarity overall.
        ref_index = similarity_matrix.sum(axis=1).argmax()
        ref = ref_candidates[ref_index]
        _custom_log("Selected reference: {}".format(ref['stem']), logger)

    _custom_log("Starting peak warping to reference", logger)
    time_start = time.time()
    ref_stem = ref['stem']
    ref_peaks = pastaq.read_peaks(os.path.join(output_dir, 'peaks', '{}.peaks'.format(ref_stem)))
    for input_file in input_files:
        stem = input_file['stem']
        # Check if file has already been processed.
        in_path = os.path.join(output_dir, 'peaks', "{}.peaks".format(stem))
        out_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
        out_path_tmap = os.path.join(output_dir, 'time_map', "{}.tmap".format(stem))

        peaks = pastaq.read_peaks(in_path)

        if not os.path.exists(out_path_tmap) or force_override:
            _custom_log("Calculating time_map for {}".format(stem), logger)
            time_map = pastaq.calculate_time_map(
                ref_peaks, peaks,
                params['warp2d_slack'],
                params['warp2d_window_size'],
                params['warp2d_num_points'],
                params['warp2d_rt_expand_factor'],
                params['warp2d_peaks_per_window'])
            pastaq.write_time_map(time_map, out_path_tmap)

        if os.path.exists(out_path) and not force_override:
            continue
        if stem != ref_stem:
            _custom_log("Warping {} peaks to reference {}".format(stem, ref_stem), logger)
            peaks = pastaq.warp_peaks(peaks, time_map)
        pastaq.write_peaks(peaks, out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished peak warping to reference in {}'.format(elapsed_time), logger)

def perform_feature_detection(params, output_dir, logger=None, force_override=False):
    _custom_log('Starting feature detection', logger)
    time_start = time.time()
    for input_file in params['input_files']:
        stem = input_file['stem']
        # Check if file has already been processed.
        in_path_peaks = os.path.join(output_dir, 'warped_peaks', '{}.peaks'.format(stem))
        out_path = os.path.join(output_dir, 'features', '{}.features'.format(stem))
        if os.path.exists(out_path) and not force_override:
            continue

        _custom_log("Reading peaks from disk: {}".format(stem), logger)
        peaks = pastaq.read_peaks(in_path_peaks)

        _custom_log("Performing feature_detection: {}".format(stem), logger)
        features = pastaq.detect_features(
            peaks, params['feature_detection_charge_states'])
        _custom_log('Writing features: {}'.format(out_path), logger)
        pastaq.write_features(features, out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished feature detection in {}'.format(elapsed_time), logger)

def parse_mzidentml_files(params, output_dir, logger=None, force_override=False):
    _custom_log('Starting mzIdentML parsing', logger)
    time_start = time.time()
    for input_file in params['input_files']:
        stem = input_file['stem']
        in_path = input_file['ident_path']
        out_path = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        if in_path == 'none' or (os.path.exists(out_path) and not force_override):
            continue
        _custom_log('Reading mzIdentML: {}'.format(in_path), logger)
        # TODO: We may want to add an option to pass a prefix for ignoring
        # decoys when they are not properly annotated, for example in msfragger
        # + idconvert
        ident_data = pastaq.read_mzidentml(
            in_path,
            ignore_decoy=params['ident_ignore_decoy'],
            require_threshold=params['ident_require_threshold'],
            max_rank_only=params['ident_max_rank_only'],
            min_mz=params['min_mz'],
            max_mz=params['max_mz'],
            min_rt=params['min_rt'],
            max_rt=params['max_rt'],
        )
        _custom_log('Writing ident data: {}'.format(out_path), logger)
        pastaq.write_ident_data(ident_data, out_path)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished mzIdentML parsing in {}'.format(elapsed_time), logger)

def link_peaks_msms_idents(params, output_dir, logger=None, force_override=False):
    _custom_log('Starting ident/msms linkage', logger)
    time_start = time.time()
    for input_file in params['input_files']:
        stem = input_file['stem']

        # Check if file has already been processed.
        in_path_raw = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_peaks = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
        in_path_idents = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        out_path_peak_ms2 = os.path.join(output_dir, 'linking', "{}.peak_ms2.link".format(stem))
        out_path_ident_ms2 = os.path.join(output_dir, 'linking', "{}.ident_ms2.link".format(stem))
        out_path_psm = os.path.join(output_dir, 'linking', "{}.ident_peak.link".format(stem))

        raw_data = None
        peaks = None
        ident_data = None

        if not os.path.exists(out_path_peak_ms2) or force_override:
            _custom_log("Performing peaks-msms linkage: {}".format(stem), logger)
            if raw_data is None:
                raw_data = pastaq.read_raw_data(in_path_raw)
            if peaks is None:
                peaks = pastaq.read_peaks(in_path_peaks)
            linked_msms = pastaq.link_peaks(
                    peaks,
                    raw_data,
                    params['link_n_sig_mz'],
                    params['link_n_sig_rt'],
            )
            _custom_log('Writing linked_msms: {}'.format(out_path_peak_ms2), logger)
            pastaq.write_linked_msms(linked_msms, out_path_peak_ms2)

        # Check that we had identification info.
        if input_file['ident_path'] == 'none':
            continue

        if not os.path.exists(out_path_ident_ms2) or force_override:
            _custom_log("Performing ident-msms linkage: {}".format(stem), logger)
            if ident_data is None:
                ident_data = pastaq.read_ident_data(in_path_idents)
            if raw_data is None:
                raw_data = pastaq.read_raw_data(in_path_raw)

            linked_idents = pastaq.link_idents(
                    ident_data,
                    raw_data,
                    params['link_n_sig_mz'],
                    params['link_n_sig_rt'],
            )
            _custom_log('Writing linked_msms: {}'.format(out_path_ident_ms2), logger)
            pastaq.write_linked_msms(linked_idents, out_path_ident_ms2)

        if not os.path.exists(out_path_psm) or force_override:
            _custom_log("Performing ident-peaks linkage: {}".format(stem), logger)
            if ident_data is None:
                ident_data = pastaq.read_ident_data(in_path_idents)
            if raw_data is None:
                raw_data = pastaq.read_raw_data(in_path_raw)
            if peaks is None:
                peaks = pastaq.read_peaks(in_path_peaks)
            linked_psm = pastaq.link_psm(
                    ident_data,
                    peaks,
                    raw_data,
                    params['link_n_sig_mz'],
                    params['link_n_sig_rt'],
            )
            _custom_log('Writing linked_psm: {}'.format(out_path_psm), logger)
            pastaq.write_linked_psm(linked_psm, out_path_psm)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished ident/msms linkage in {}'.format(elapsed_time), logger)

def match_peaks_and_features(params, output_dir, logger=None, force_override=False):
    input_files = params['input_files']

    # Transform groups from string into an integer list.
    groups = []
    group_map = {}
    group_counter = 0
    for input_file in input_files:
        group = input_file['group']
        if group in group_map:
            groups += [group_map[group]]
            continue
        group_map[group] = group_counter
        groups += [group_counter]
        group_counter += 1

    # Peak matching.
    _custom_log('Starting peak matching', logger)
    time_start = time.time()
    in_path_peaks = os.path.join(output_dir, 'warped_peaks')
    out_path = os.path.join(output_dir, 'metamatch', "peaks.clusters")
    if (not os.path.exists(out_path) or force_override):
        _custom_log("Reading peaks from disk", logger)
        peaks = [
            pastaq.read_peaks(os.path.join(in_path_peaks, "{}.peaks".format(input_file['stem'])))
            for input_file in input_files
        ]
        _custom_log("Finding peak clusters", logger)
        peak_clusters = pastaq.find_peak_clusters(
            groups,
            peaks,
            params["metamatch_fraction"],
            params["metamatch_n_sig_mz"],
            params["metamatch_n_sig_rt"])
        _custom_log("Writing peak clusters to disk", logger)
        pastaq.write_peak_clusters(peak_clusters, out_path)
    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished peak matching in {}'.format(elapsed_time), logger)

    # Feature matching.
    _custom_log('Starting feature matching', logger)
    time_start = time.time()
    in_path_features = os.path.join(output_dir, 'features')
    out_path = os.path.join(output_dir, 'metamatch', "features.clusters")
    if (not os.path.exists(out_path) or force_override):
        _custom_log("Reading features from disk", logger)
        features = [
            pastaq.read_features(os.path.join(in_path_features, "{}.features".format(input_file['stem'])))
            for input_file in input_files
        ]
        _custom_log("Finding feature clusters", logger)
        feature_clusters = pastaq.find_feature_clusters(
            groups,
            features,
            params["metamatch_fraction"],
            params["metamatch_n_sig_mz"],
            params["metamatch_n_sig_rt"])
        _custom_log("Writing feature clusters to disk", logger)
        pastaq.write_feature_clusters(feature_clusters, out_path)
    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished feature matching in {}'.format(elapsed_time), logger)

# NOTE: This is a giant ball of spaghetti and could use some love.
def create_quantitative_tables(params, output_dir, logger=None, force_override=False):
    input_files = params['input_files']

    _custom_log('Starting creation of quantitative tables', logger)
    time_start = time.time()
    for input_file in input_files:
        stem = input_file['stem']
        # Peak quantification.
        # ====================
        in_path_peaks = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
        in_path_peaks_link = os.path.join(output_dir, 'linking', "{}.peak_ms2.link".format(stem))
        in_path_ident_link_msms = os.path.join(output_dir, 'linking', "{}.ident_ms2.link".format(stem))
        in_path_ident_link_theomz = os.path.join(output_dir, 'linking', "{}.ident_peak.link".format(stem))
        in_path_ident_data = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        out_path_peaks = os.path.join(output_dir, 'quant',"{}_peaks.csv".format(stem))
        out_path_peak_annotations = os.path.join(output_dir, 'quant',"{}_peak_annotations.csv".format(stem))

        # TODO: This is probably not necessary or needs to be changed if we are
        # doing all per-peak quantification in a single loop.
        if os.path.exists(out_path_peaks) and not force_override:
            continue

        _custom_log("Reading peaks from disk: {}".format(stem), logger)
        peaks = pastaq.read_peaks(in_path_peaks)

        _custom_log("Generating peaks quantitative table", logger)
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
        _custom_log("Reading linked peaks from disk: {}".format(stem), logger)
        peak_annotations = peaks_df[["peak_id"]]
        linked_peaks = pastaq.read_linked_msms(in_path_peaks_link)
        linked_peaks = pd.DataFrame({
            'peak_id': [linked_peak.entity_id for linked_peak in linked_peaks],
            'msms_id': [linked_peak.msms_id for linked_peak in linked_peaks],
        })
        peak_annotations = pd.merge(
            peak_annotations, linked_peaks, on="peak_id", how="left")

        if os.path.isfile(in_path_ident_data):
            _custom_log("Reading ident_data from disk: {}".format(stem), logger)
            ident_data = pastaq.read_ident_data(in_path_ident_data)
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
                if params["quant_ident_linkage"] == 'theoretical_mz':
                    _custom_log(
                        "Reading linked ident_peak from disk: {}".format(stem), logger)
                    linked_idents = pastaq.read_linked_psm(
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
                elif params["quant_ident_linkage"] == 'msms_event':
                    _custom_log(
                        "Reading linked ident_peak from disk: {}".format(stem), logger)
                    linked_idents = pastaq.read_linked_msms(
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

        _custom_log("Saving peaks quantitative table to disk: {}".format(stem), logger)
        peaks_df.to_csv(out_path_peaks, index=False)

        _custom_log("Saving peaks annotations table to disk: {}".format(stem), logger)
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

            _custom_log("Reading features from disk: {}".format(stem), logger)
            features = pastaq.read_features(in_path_features)

            _custom_log("Generating features quantitative table", logger)
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
    _custom_log("Reading peak clusters from disk", logger)
    in_path_peak_clusters = os.path.join(
        output_dir, 'metamatch', 'peaks.clusters')
    out_path_peak_clusters_metadata = os.path.join(output_dir, 'quant',
                                                   "peak_clusters_metadata.csv")
    out_path_peak_clusters_peaks = os.path.join(output_dir, 'quant',
                                                "peak_clusters_peaks.csv")
    out_path_peak_clusters_annotations = os.path.join(output_dir, 'quant',
                                                      "peak_clusters_annotations.csv")

    # TODO: This is clearly suboptimal and will take a long time if the
    # number of clusters is very high. Needs a rewrite for optimal
    # performance.
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

    if (not os.path.exists(out_path_peak_clusters_metadata) or force_override):
        peak_clusters = pastaq.read_peak_clusters(in_path_peak_clusters)
        _custom_log("Generating peak clusters quantitative table", logger)
        peak_clusters_metadata_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in peak_clusters],
            'mz': [cluster.mz for cluster in peak_clusters],
            'rt': [cluster.rt for cluster in peak_clusters],
            'avg_height': [cluster.avg_height for cluster in peak_clusters],
        })

        _custom_log("Generating peak clusters quantitative table", logger)
        peak_clusters_df = pd.DataFrame({
            'cluster_id': [cluster.id for cluster in peak_clusters],
        })
        if params['quant_isotopes'] == 'volume':
            out_path_peak_clusters = os.path.join(output_dir, 'quant',
                                                  "peak_clusters_volume.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                peak_clusters_df[stem] = [cluster.file_volumes[i]
                                          for cluster in peak_clusters]
        elif params['quant_isotopes'] == 'height':
            out_path_peak_clusters = os.path.join(output_dir, 'quant',
                                                  "peak_clusters_height.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                peak_clusters_df[stem] = [cluster.heights[i]
                                          for cluster in peak_clusters]
        else:
            raise ValueError("unknown quant_isotopes parameter")
        _custom_log("Writing peaks quantitative table to disk", logger)
        peak_clusters_df.to_csv(out_path_peak_clusters, index=False)

        # Peak associations.
        _custom_log("Generating peak clusters peak associations table", logger)
        cluster_peaks = [(cluster.id, cluster.peak_ids)
                         for cluster in peak_clusters]
        cluster_peaks = pd.DataFrame(cluster_peaks, columns=[
                                     "cluster_id", "peak_ids"]).explode("peak_ids")
        cluster_peaks["file_id"] = cluster_peaks["peak_ids"].map(
            lambda x: input_files[x.file_id]['stem'])
        cluster_peaks["peak_id"] = cluster_peaks["peak_ids"].map(
            lambda x: x.peak_id)
        cluster_peaks = cluster_peaks.drop(["peak_ids"], axis=1)
        _custom_log("Writing cluster to peak table to disk", logger)
        cluster_peaks.to_csv(out_path_peak_clusters_peaks, index=False)

        # Cluster annotations.
        _custom_log("Generating peak clusters annotations table", logger)
        annotations = pd.DataFrame()
        for input_file in input_files:
            stem = input_file['stem']
            _custom_log("Reading peak annotations for: {}".format(stem), logger)
            in_path_peak_annotations = os.path.join(output_dir, 'quant',
                                                    "{}_peak_annotations.csv".format(stem))
            peak_annotations = pd.read_csv(
                in_path_peak_annotations, low_memory=False)
            cluster_annotations = cluster_peaks[cluster_peaks["file_id"] == stem][[
                "cluster_id", "peak_id"]]
            cluster_annotations["file_id"] = stem
            _custom_log(
                "Merging peak/clusters annotations for: {}".format(stem), logger)
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

        if params['quant_consensus'] and 'psm_sequence' in annotations:
            # Find a sequence consensus
            consensus_sequence = find_sequence_consensus(
                annotations, 'psm_sequence', params['quant_consensus_min_ident'])
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
        if params['quant_save_all_annotations']:
            _custom_log("Writing annotations to disk", logger)
            annotations = annotations.sort_values(by=["cluster_id"])
            annotations.to_csv(out_path_peak_clusters_annotations, index=False)

        _custom_log("Aggregating annotations", logger)
        # TODO: This is clearly suboptimal and will take a long time if the
        # number of clusters is very high. Needs a rewrite for optimal
        # performance.
        annotations = annotations.groupby(
            'cluster_id').apply(aggregate_cluster_annotations)

        # Metadata
        _custom_log("Merging metadata with annotations", logger)
        peak_clusters_metadata_df = pd.merge(
            peak_clusters_metadata_df, annotations, how="left", on="cluster_id")
        _custom_log("Writing metadata to disk", logger)
        peak_clusters_metadata_df.to_csv(
            out_path_peak_clusters_metadata, index=False)

    # Matched Features
    # ================
    _custom_log("Reading feature clusters from disk", logger)
    in_path_feature_clusters = os.path.join(
        output_dir, 'metamatch', 'features.clusters')
    out_path_feature_clusters_metadata = os.path.join(output_dir, 'quant',
                                                      "feature_clusters_metadata.csv")
    out_path_feature_clusters_features = os.path.join(output_dir, 'quant',
                                                      "feature_clusters_features.csv")
    out_path_feature_clusters_annotations = os.path.join(output_dir, 'quant',
                                                         "feature_clusters_annotations.csv")
    if (not os.path.exists(out_path_feature_clusters_metadata) or force_override):
        feature_clusters = pastaq.read_feature_clusters(
            in_path_feature_clusters)

        _custom_log("Generating feature clusters quantitative table", logger)
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
        if params['quant_features'] == 'monoisotopic_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_monoisotopic_height.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.monoisotopic_heights[i]
                                             for cluster in feature_clusters]
        elif params['quant_features'] == 'monoisotopic_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_monoisotopic_volume.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.monoisotopic_volumes[i]
                                             for cluster in feature_clusters]
        elif params['quant_features'] == 'total_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_total_height.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.total_heights[i]
                                             for cluster in feature_clusters]
        elif params['quant_features'] == 'total_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_total_volume.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.total_volumes[i]
                                             for cluster in feature_clusters]
        elif params['quant_features'] == 'max_height':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_max_height.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.max_heights[i]
                                             for cluster in feature_clusters]
        elif params['quant_features'] == 'max_volume':
            out_path_feature_clusters = os.path.join(output_dir, 'quant',
                                                     "feature_clusters_max_volume.csv")
            for i, input_file in enumerate(input_files):
                stem = input_file['stem']
                data[stem] = [cluster.max_volumes[i]
                                             for cluster in feature_clusters]
        else:
            raise ValueError("unknown quant_features parameter")
        _custom_log("Writing feature clusters quantitative table to disk", logger)
        data.to_csv(out_path_feature_clusters, index=False)

        # Feature associations.
        _custom_log("Generating feature clusters feature associations table", logger)
        cluster_features = [(cluster.id, cluster.feature_ids)
                            for cluster in feature_clusters]
        cluster_features = pd.DataFrame(cluster_features, columns=[
                                        "cluster_id", "feature_ids"]).explode("feature_ids")
        cluster_features["file_id"] = cluster_features["feature_ids"].map(
            lambda x: input_files[x.file_id]['stem'])
        cluster_features["feature_id"] = cluster_features["feature_ids"].map(
            lambda x: x.feature_id)
        cluster_features = cluster_features.drop(["feature_ids"], axis=1)
        _custom_log("Writing cluster to feature table to disk", logger)
        cluster_features.to_csv(
            out_path_feature_clusters_features, index=False)

        # Cluster annotations.
        _custom_log("Generating peak clusters annotations table", logger)
        annotations = pd.DataFrame()
        for input_file in input_files:
            stem = input_file['stem']
            _custom_log("Reading features for: {}".format(stem), logger)
            in_path_peak_features = os.path.join(output_dir, 'features',
                                                 "{}.features".format(stem))
            in_path_peak_annotations = os.path.join(output_dir, 'quant',
                                                    "{}_peak_annotations.csv".format(stem))
            features = pastaq.read_features(in_path_peak_features)
            features = [(feature.id, feature.peak_ids, feature.charge_state)
                        for feature in features]
            features = pd.DataFrame(
                features, columns=["feature_id", "peak_id", "charge_state"]).explode("peak_id")
            peak_annotations = pd.read_csv(
                in_path_peak_annotations, low_memory=False)
            cluster_annotations = cluster_features[cluster_features["file_id"] == stem][[
                "cluster_id", "feature_id"]]
            cluster_annotations["file_id"] = stem
            _custom_log(
                "Merging peak/clusters annotations for: {}".format(stem), logger)
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

        if params['quant_consensus'] and 'psm_sequence' in annotations:
            # Find a sequence consensus
            consensus_sequence = find_sequence_consensus(
                annotations, 'psm_sequence', params['quant_consensus_min_ident'])
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
        _custom_log("Calculating protein groups", logger)
        sequence_column = 'psm_sequence'
        protein_name_column = 'protein_name'
        protein_description_column = 'protein_description'
        if params['quant_consensus'] and 'psm_sequence' in annotations:
            sequence_column = 'consensus_sequence'
            protein_name_column = 'consensus_protein_name'
            protein_description_column = 'consensus_protein_description'

        # Combine protein name/description in case
        if (sequence_column in annotations and
                protein_name_column in annotations and
                protein_description_column in annotations):
            prot_data, prot_metadata = find_protein_groups(
                    data,
                    annotations,
                    sequence_column,
                    protein_name_column,
                    protein_description_column,
                    params['quant_proteins_min_peptides'],
                    params['quant_proteins_remove_subset_proteins'],
                    params['quant_proteins_ignore_ambiguous_peptides'],
                    params['quant_proteins_quant_type'],
                    )
            out_path_protein_data = os.path.join(output_dir, 'quant',
                    "protein_groups.csv")
            out_path_protein_metadata = os.path.join(output_dir, 'quant',
                    "protein_groups_metadata.csv")

            _custom_log("Writing protein group data/metadata to disk", logger)
            prot_data.to_csv(out_path_protein_data, index=False)
            prot_metadata.to_csv(out_path_protein_metadata, index=False)


        # Saving annotations before aggregation.
        if params['quant_save_all_annotations']:
            _custom_log("Writing annotations to disk", logger)
            annotations = annotations.sort_values(by=["cluster_id"])
            annotations.to_csv(
                out_path_feature_clusters_annotations, index=False)

        _custom_log("Aggregating annotations", logger)
        if ("psm_charge_state" in annotations and
                params['quant_features_charge_state_filter']):
            annotations = annotations[annotations["psm_charge_state"]
                                      == annotations["charge_state"]]
        # TODO: This is clearly suboptimal and will take a long time if the
        # number of clusters is very high. Needs a rewrite for optimal
        # performance.
        annotations_agg = annotations.groupby(
            'cluster_id').apply(aggregate_cluster_annotations)

        # Metadata.
        _custom_log("Merging metadata with annotations", logger)
        metadata = pd.merge(
            metadata, annotations_agg, how="left", on="cluster_id")
        _custom_log("Writing metadata to disk", logger)
        metadata.to_csv(out_path_feature_clusters_metadata, index=False)

        # Aggregate peptides.
        sequence_column = 'psm_sequence'
        if params['quant_consensus'] and 'psm_sequence' in annotations:
            sequence_column = 'consensus_sequence'

        if sequence_column in annotations:
            _custom_log("Aggregating peptide charge states", logger)
            def aggregate_peptide_annotations(x):
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
            peptide_data = data.copy()
            peptide_data = peptide_data.drop(["cluster_id"], axis=1)
            peptide_data[sequence_column] = metadata[sequence_column]
            peptide_data = peptide_data[~peptide_data[sequence_column].isna()]
            peptide_data = peptide_data[peptide_data[sequence_column] != '']
            peptide_data = peptide_data[~peptide_data[sequence_column].str.contains('\.\|\.')]
            peptide_data = peptide_data.groupby(sequence_column).agg(sum)
            peptide_data = peptide_data.reset_index()
            peptide_metadata = annotations.copy()
            peptide_metadata = peptide_metadata[peptide_metadata[sequence_column].isin(peptide_data[sequence_column])]
            peptide_metadata = peptide_metadata.groupby(sequence_column)
            peptide_metadata = peptide_metadata.apply(aggregate_peptide_annotations)
            out_path_peptide_data = os.path.join(output_dir, 'quant',
                    "peptides_data.csv")
            out_path_peptide_metadata = os.path.join(output_dir, 'quant',
                    "peptides_metadata.csv")

            _custom_log("Writing peptide data/metadata to disk", logger)
            peptide_data.to_csv(out_path_peptide_data, index=False)
            peptide_metadata.to_csv(out_path_peptide_metadata, index=False)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished creation of quantitative tables in {}'.format(elapsed_time), logger)

def dda_pipeline_summary(params, output_dir, logger):
    _custom_log("Starting summary stats", logger)
    time_start = time.time()

    input_files = params['input_files']

    summary_log = logging.getLogger('summary')
    summary_log.setLevel(logging.INFO)
    summary_log_fh = logging.FileHandler(os.path.join(output_dir, 'summary.log'))
    summary_log_fh.setLevel(logging.INFO)
    summary_log_formatter = logging.Formatter('%(message)s')
    summary_log_fh.setFormatter(summary_log_formatter)
    summary_log.addHandler(summary_log_fh)

    # Raw data
    summary_log.info('Raw data')
    for input_file in input_files:
        stem = input_file['stem']
        summary_log.info('    {}'.format(stem))

        # MS1
        in_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
        if os.path.exists(in_path):
            raw_data = pastaq.read_raw_data(in_path)
            summary_log.info('        MS1')
            summary_log.info('            number of scans: {}'.format(len(raw_data.scans)))
            summary_log.info('            min_mz: {}'.format(raw_data.min_mz))
            summary_log.info('            max_mz: {}'.format(raw_data.max_mz))
            summary_log.info('            min_rt: {}'.format(raw_data.min_rt))
            summary_log.info('            max_rt: {}'.format(raw_data.max_rt))

        # MS2
        in_path = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        if os.path.exists(in_path):
            raw_data = pastaq.read_raw_data(in_path)
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
    for input_file in input_files:
        stem = input_file['stem']
        summary_log.info('    {}'.format(stem))

        in_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
        if os.path.exists(in_path):
            peaks = pastaq.read_peaks(in_path)
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
    for input_file in input_files:
        stem = input_file['stem']
        summary_log.info('    {}'.format(stem))

        in_path = os.path.join(output_dir, 'features', "{}.features".format(stem))
        if os.path.exists(in_path):
            features = pastaq.read_features(in_path)

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
    for input_file in input_files:
        stem = input_file['stem']
        summary_log.info('    {}'.format(stem))

        in_path_raw_data = os.path.join(output_dir, 'raw', "{}.ms2".format(stem))
        in_path_linked_msms = os.path.join(output_dir, 'linking', "{}.peak_ms2.link".format(stem))
        if os.path.exists(in_path_raw_data) and os.path.exists(in_path_linked_msms):
            raw_data = pastaq.read_raw_data(in_path_raw_data)
            linked_msms = pastaq.read_linked_msms(in_path_linked_msms)
            summary_log.info('        MS/MS-Peaks linkage')
            summary_log.info('            Number of ms/ms events: {}'.format(len(raw_data.scans)))
            summary_log.info('            Number of ms/ms events linked to peaks: {}'.format(len(linked_msms)))
            if len(raw_data.scans) > 0:
                summary_log.info('            Linking efficiency (%): {}'.format(len(linked_msms)/len(raw_data.scans) * 100.0))

        in_path_ident_data = os.path.join(output_dir, 'ident', "{}.ident".format(stem))
        if os.path.exists(in_path_ident_data):
            ident_data = pastaq.read_ident_data(in_path_ident_data)

            in_path_ident_ms2 = os.path.join(output_dir, 'linking', "{}.ident_ms2.link".format(stem))
            if os.path.exists(in_path_ident_ms2):
                ident_ms2 = pastaq.read_linked_msms(in_path_ident_ms2)
                summary_log.info('        MS/MS-Identification linkage')
                summary_log.info('            Number of PSMs: {}'.format(len(ident_data.spectrum_matches)))
                summary_log.info('            Number of PSMs linked to MS/MS events: {}'.format(len(ident_ms2)))
                if len(ident_data.spectrum_matches) > 0:
                    summary_log.info('            PSM-peaks linking efficiency (%): {}'.format(len(ident_ms2)/len(ident_data.spectrum_matches) * 100.0))

            in_path_peak_idents = os.path.join(output_dir, 'linking', "{}.ident_peak.link".format(stem))
            if os.path.exists(in_path_peak_idents):
                ident_peak = pastaq.read_linked_psm(in_path_peak_idents)
                summary_log.info('        Peaks-Identification linkage')
                summary_log.info('            Number of PSMs: {}'.format(len(ident_data.spectrum_matches)))
                summary_log.info('            Number of PSMs linked to peaks: {}'.format(len(ident_peak)))
                if len(ident_data.spectrum_matches) > 0:
                    summary_log.info('            PSM-peaks linking efficiency (%): {}'.format(len(ident_peak)/len(ident_data.spectrum_matches) * 100.0))

    # TODO: Average identification linkage stats.
    # TODO: Metamatch stats
    # TODO: Peptide stats
    # TODO: Protein group stats
    summary_log.removeHandler(summary_log_fh)
    summary_log_fh.close()

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished summary stats in {}'.format(elapsed_time), logger)

def generate_qc_plots(params, output_dir, logger=None, force_override=False):
    input_files = params['input_files']

    #
    # General plot config.
    #

    # Font and figure size
    plt.rcParams.update({
        'font.family': params['qc_plot_font_family'],
        'font.size': params['qc_plot_font_size'],
        'figure.figsize': (params['qc_plot_fig_size_x'], params['qc_plot_fig_size_y']),
    })

    # Alpha parameters.
    fill_alpha = params['qc_plot_fill_alpha']
    line_alpha = params['qc_plot_line_alpha']
    scatter_alpha = params['qc_plot_scatter_alpha']
    if line_alpha == 'dynamic':
        line_alpha = max(params['qc_plot_min_dynamic_alpha'], 1.0 / len(input_files))
    if fill_alpha == 'dynamic':
        fill_alpha = max(params['qc_plot_min_dynamic_alpha'], 1.0 / len(input_files))

    # Colorscheme.
    palette = sns.color_palette(params['qc_plot_palette'], len(input_files))

    _custom_log("Starting quality control plotting", logger)
    time_start = time.time()

    #
    # Peak sigma density
    #
    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'peak_sigma_mz_rt_density.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, (ax_left, ax_right) = plt.subplots(1, 2)
            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting density of sigma_mz/sigma_rt: {}".format(stem), logger)
                peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
                peaks = pastaq.read_peaks(peaks_path)

                sigma_mzs = np.array([peak.fitted_sigma_mz for peak in peaks])
                sigma_rts = np.array([peak.fitted_sigma_rt for peak in peaks])

                if params['qc_plot_line_style'] == 'fill':
                    sns.kdeplot(sigma_rts, label=stem, ax=ax_left, fill=True, linewidth=0, alpha=fill_alpha, color=color)
                    sns.kdeplot(sigma_mzs, label=stem, ax=ax_right, fill=True, linewidth=0, alpha=fill_alpha, color=color)
                else:
                    sns.kdeplot(sigma_rts, label=stem, ax=ax_left, alpha=line_alpha, color=color)
                    sns.kdeplot(sigma_mzs, label=stem, ax=ax_right, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax_left.set_xlabel('$\\sigma_{rt}$')
            ax_right.set_xlabel('$\\sigma_{mz}$')
            ax_left.set_ylabel('Density')
            ax_right.set_ylabel('')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_peak_sigma_mz_rt_density.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, (ax_left, ax_right) = plt.subplots(1, 2)
            _custom_log("Plotting density of sigma_mz/sigma_rt: {}".format(stem), logger)
            peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
            peaks = pastaq.read_peaks(peaks_path)

            sigma_mzs = np.array([peak.fitted_sigma_mz for peak in peaks])
            sigma_rts = np.array([peak.fitted_sigma_rt for peak in peaks])

            if params['qc_plot_line_style'] == 'fill':
                sns.kdeplot(sigma_rts, label=stem, ax=ax_left, fill=True, linewidth=0, alpha=fill_alpha, color=color)
                sns.kdeplot(sigma_mzs, label=stem, ax=ax_right, fill=True, linewidth=0, alpha=fill_alpha, color=color)
            else:
                sns.kdeplot(sigma_rts, label=stem, ax=ax_left, alpha=line_alpha, color=color)
                sns.kdeplot(sigma_mzs, label=stem, ax=ax_right, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax_left.set_xlabel('$\\sigma_{rt}$')
            ax_right.set_xlabel('$\\sigma_{mz}$')
            ax_left.set_ylabel('Density')
            ax_right.set_ylabel('')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    #
    # Peak rt vs rt_delta
    #
    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'peak_rt_vs_rt_delta.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting rt vs rt_delta: {}".format(stem), logger)
                peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
                peaks = pastaq.read_peaks(peaks_path)

                rts = np.array([peak.fitted_rt for peak in peaks])
                rt_deltas = np.array([peak.rt_delta for peak in peaks])

                idx = np.argsort(rts)
                rts = rts[idx]
                rt_deltas = rt_deltas[idx]
                ax.plot(rts, rt_deltas, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Retention time delta (s)')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_peak_rt_vs_rt_delta.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting rt vs rt_delta: {}".format(stem), logger)
            peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
            peaks = pastaq.read_peaks(peaks_path)

            rts = np.array([peak.fitted_rt for peak in peaks])
            rt_deltas = np.array([peak.rt_delta for peak in peaks])

            idx = np.argsort(rts)
            rts = rts[idx]
            rt_deltas = rt_deltas[idx]
            ax.plot(rts, rt_deltas, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Retention time delta (s)')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    #
    # Peak sigma_mz vs m/z scatterplot.
    #
    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'peak_mz_vs_sigma_mz.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting mz vs sigma_mz: {}".format(stem), logger)
                peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
                peaks = pastaq.read_peaks(peaks_path)

                mz = np.array([peak.fitted_mz for peak in peaks])[0:params['qc_plot_mz_vs_sigma_mz_max_peaks']]
                sigma_mz = np.array([peak.fitted_sigma_mz for peak in peaks])[0:params['qc_plot_mz_vs_sigma_mz_max_peaks']]

                ax.scatter(mz, sigma_mz, s=params['qc_plot_scatter_size'], label=stem, edgecolors='none', alpha=scatter_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('m/z')
            ax.set_ylabel('$\\sigma_{mz}$')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_peak_mz_vs_sigma_mz.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting mz vs sigma_mz: {}".format(stem), logger)
            peaks_path = os.path.join(output_dir, 'warped_peaks', "{}.peaks".format(stem))
            peaks = pastaq.read_peaks(peaks_path)

            mz = np.array([peak.fitted_mz for peak in peaks])[0:params['qc_plot_mz_vs_sigma_mz_max_peaks']]
            sigma_mz = np.array([peak.fitted_sigma_mz for peak in peaks])[0:params['qc_plot_mz_vs_sigma_mz_max_peaks']]

            ax.scatter(mz, sigma_mz, s=params['qc_plot_scatter_size'], label=stem, edgecolors='none', alpha=scatter_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('m/z')
            ax.set_ylabel('$\\sigma_{mz}$')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    #
    # Extracted Ion Chromatogram (XIC) before and after alignment.
    #
    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'xic_unaligned.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting XIC (unaligned): {}".format(stem), logger)

                raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
                raw_data = pastaq.read_raw_data(raw_data_path)
                xic = pastaq.xic(
                    raw_data,
                    raw_data.min_mz,
                    raw_data.max_mz,
                    raw_data.min_rt,
                    raw_data.max_rt,
                    "sum"
                )
                x = xic.retention_time
                y = xic.intensity

                if params['qc_plot_line_style'] == 'fill':
                    ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
                else:
                    ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_xic_unaligned.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting XIC (unaligned): {}".format(stem), logger)

            raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
            raw_data = pastaq.read_raw_data(raw_data_path)
            xic = pastaq.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "sum"
            )
            x = xic.retention_time
            y = xic.intensity

            if params['qc_plot_line_style'] == 'fill':
                ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
            else:
                ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'xic_aligned.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting XIC (aligned): {}".format(stem), logger)

                raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
                tmap_path = os.path.join(output_dir, 'time_map', "{}.tmap".format(stem))

                raw_data = pastaq.read_raw_data(raw_data_path)
                tmap = pastaq.read_time_map(tmap_path)

                xic = pastaq.xic(
                    raw_data,
                    raw_data.min_mz,
                    raw_data.max_mz,
                    raw_data.min_rt,
                    raw_data.max_rt,
                    "sum"
                )
                x = [tmap.warp(rt) for rt in xic.retention_time]
                y = xic.intensity

                if params['qc_plot_line_style'] == 'fill':
                    ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
                else:
                    ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_xic_aligned.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting XIC (aligned): {}".format(stem), logger)

            raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
            tmap_path = os.path.join(output_dir, 'time_map', "{}.tmap".format(stem))

            raw_data = pastaq.read_raw_data(raw_data_path)
            tmap = pastaq.read_time_map(tmap_path)

            xic = pastaq.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "sum"
            )
            x = [tmap.warp(rt) for rt in xic.retention_time]
            y = xic.intensity

            if params['qc_plot_line_style'] == 'fill':
                ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
            else:
                ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    #
    # Base Peak chromatogram before and after alignment.
    #
    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'bpc_unaligned.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting Base Peak Chromatogram (unaligned): {}".format(stem), logger)

                raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
                raw_data = pastaq.read_raw_data(raw_data_path)
                xic = pastaq.xic(
                    raw_data,
                    raw_data.min_mz,
                    raw_data.max_mz,
                    raw_data.min_rt,
                    raw_data.max_rt,
                    "max"
                )
                x = xic.retention_time
                y = xic.intensity

                if params['qc_plot_line_style'] == 'fill':
                    ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
                else:
                    ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_bpc_unaligned.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting Base Peak Chromatogram (unaligned): {}".format(stem), logger)

            raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
            raw_data = pastaq.read_raw_data(raw_data_path)
            xic = pastaq.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "max"
            )
            x = xic.retention_time
            y = xic.intensity

            if params['qc_plot_line_style'] == 'fill':
                ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
            else:
                ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    if not params['qc_plot_per_file']:
        out_path = os.path.join(output_dir, 'quality', 'bpc_aligned.{}'.format(params['qc_plot_extension']))
        if not os.path.exists(out_path) or force_override:
            fig, ax = plt.subplots(1, 1)

            for input_file, color in zip(input_files, palette):
                stem = input_file['stem']
                _custom_log("Plotting Base Peak Chromatogram (aligned): {}".format(stem), logger)

                raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
                tmap_path = os.path.join(output_dir, 'time_map', "{}.tmap".format(stem))

                raw_data = pastaq.read_raw_data(raw_data_path)
                tmap = pastaq.read_time_map(tmap_path)

                xic = pastaq.xic(
                    raw_data,
                    raw_data.min_mz,
                    raw_data.max_mz,
                    raw_data.min_rt,
                    raw_data.max_rt,
                    "max"
                )
                x = [tmap.warp(rt) for rt in xic.retention_time]
                y = xic.intensity

                if params['qc_plot_line_style'] == 'fill':
                    ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
                else:
                    ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)
    else:
        for input_file, color in zip(input_files, palette):
            stem = input_file['stem']
            out_path = os.path.join(output_dir, 'quality', '{}_bpc_aligned.{}'.format(stem, params['qc_plot_extension']))
            if os.path.exists(out_path) and not force_override:
                continue
            fig, ax = plt.subplots(1, 1)
            _custom_log("Plotting Base Peak Chromatogram (aligned): {}".format(stem), logger)

            raw_data_path = os.path.join(output_dir, 'raw', "{}.ms1".format(stem))
            tmap_path = os.path.join(output_dir, 'time_map', "{}.tmap".format(stem))

            raw_data = pastaq.read_raw_data(raw_data_path)
            tmap = pastaq.read_time_map(tmap_path)

            xic = pastaq.xic(
                raw_data,
                raw_data.min_mz,
                raw_data.max_mz,
                raw_data.min_rt,
                raw_data.max_rt,
                "max"
            )
            x = [tmap.warp(rt) for rt in xic.retention_time]
            y = xic.intensity

            if params['qc_plot_line_style'] == 'fill':
                ax.fill_between(x, 0, y, lw=0, color=color, alpha=fill_alpha, label=stem)
            else:
                ax.plot(x, y, label=stem, alpha=line_alpha, color=color)

            if params['qc_plot_fig_legend']:
                plt.legend()
            ax.set_xlabel('Retention time (s)')
            ax.set_ylabel('Intensity')
            _custom_log("Saving figure: {}".format(out_path), logger)
            plt.savefig(out_path, dpi=params['qc_plot_dpi'])
            plt.close(fig)

    #
    # Similarity matrix before/after alignment.
    #
    out_path = os.path.join(output_dir, 'quality', 'similarity_unaligned.{}'.format(params['qc_plot_extension']))
    if not os.path.exists(out_path) or force_override:
        fig, ax = plt.subplots(1, 1)
        _custom_log("Plotting similarity matrix before alignment", logger)
        matrix_path = os.path.join(output_dir, 'quality', 'similarity_{}.csv'.format('peaks'))
        similarity_matrix = pd.read_csv(matrix_path, index_col=0)
        sns.heatmap(similarity_matrix, xticklabels=True, yticklabels=True, square=True, vmin=0, vmax=1)
        _custom_log("Saving figure: {}".format(out_path), logger)
        plt.savefig(out_path, dpi=params['qc_plot_dpi'])
        plt.close(fig)

    out_path = os.path.join(output_dir, 'quality', 'similarity_aligned.{}'.format(params['qc_plot_extension']))
    if not os.path.exists(out_path) or force_override:
        fig, ax = plt.subplots(1, 1)
        _custom_log("Plotting similarity matrix after alignment", logger)
        matrix_path = os.path.join(output_dir, 'quality', 'similarity_{}.csv'.format('warped_peaks'))
        similarity_matrix = pd.read_csv(matrix_path, index_col=0)
        sns.heatmap(similarity_matrix, xticklabels=True, yticklabels=True, square=True, vmin=0, vmax=1)
        _custom_log("Saving figure: {}".format(out_path), logger)
        plt.savefig(out_path, dpi=params['qc_plot_dpi'])
        plt.close(fig)

    elapsed_time = datetime.timedelta(seconds=time.time()-time_start)
    _custom_log('Finished quality control plotting in {}'.format(elapsed_time), logger)


def dda_pipeline(
    pastaq_parameters,
    input_files,
    output_dir="pastaq",
    force_override=False,
    save_grid=False,
):
    # TODO: Logger should have different levels and user can configure the
    # verbosity of output.
    # TODO: Sanitize parameters.
    # TODO: Sanitize input/outputs.
    # TODO:     - Check if there are name conflicts.

    # Create output directory and subdirectoreis if necessary.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir, 'raw')):
        os.makedirs(os.path.join(output_dir, 'raw'))
    if not os.path.exists(os.path.join(output_dir, 'quality')):
        os.makedirs(os.path.join(output_dir, 'quality'))
    if save_grid:
        if not os.path.exists(os.path.join(output_dir, 'grid')):
            os.makedirs(os.path.join(output_dir, 'grid'))
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

    # Initialize logger.
    # TODO: Log to file and cout simultaneously if the user asks for it.
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

    # Prepare input files.
    for file in input_files:
        # TODO: Check if the input path for raw files exist.
        # TODO: Check if the input path for identifications exist.

        # If there is no identification, make sure it is set as none.
        if 'ident_path' not in file:
            file['ident_path'] = 'none'

        # Obtain the stem for this file if not manually specified.
        if 'stem' not in file:
            base_name = os.path.basename(file['raw_path'])
            base_name = os.path.splitext(base_name)
            file['stem'] = base_name[0]

        # Check that all files contain a group, if not, assign the default
        # 'none' group.
        if 'group' not in file:
            file['group'] = 'none'

    # Make sure the input files are in the parameters list before saving them to
    # disk.
    pastaq_parameters['input_files'] = input_files

    # Save parameters file.
    parameters_file_name = os.path.join(output_dir, 'parameters.json')
    with open(parameters_file_name, 'w') as json_file:
        json.dump(pastaq_parameters, json_file)

    # Store current time for logging the total elapsed time for the entire run.
    time_pipeline_start = time.time()

    parse_raw_files(pastaq_parameters, output_dir, logger, force_override)
    detect_peaks(pastaq_parameters, output_dir, save_grid, logger, force_override)
    calculate_similarity_matrix(pastaq_parameters, output_dir, 'peaks', logger, force_override)
    perform_rt_alignment(pastaq_parameters, output_dir, logger, force_override)
    calculate_similarity_matrix(pastaq_parameters, output_dir, 'warped_peaks', logger, force_override)
    perform_feature_detection(pastaq_parameters, output_dir, logger, force_override)
    parse_mzidentml_files(pastaq_parameters, output_dir, logger, force_override)
    link_peaks_msms_idents(pastaq_parameters, output_dir, logger, force_override)
    match_peaks_and_features(pastaq_parameters, output_dir, logger, force_override)
    create_quantitative_tables(pastaq_parameters, output_dir, logger, force_override)
    generate_qc_plots(pastaq_parameters, output_dir, logger, force_override)
    dda_pipeline_summary(pastaq_parameters, output_dir, logger)

    logger.info('Total time elapsed: {}'.format(
        datetime.timedelta(seconds=time.time()-time_pipeline_start)))

    # Stop logger.
    logger.removeHandler(logger_fh)
    logger_fh.close()
