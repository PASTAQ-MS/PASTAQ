"""
Regression test for PeakCluster attribute names.

Validates that pastaq.PeakCluster exposes 'volumes' and 'heights' attributes
(not 'file_volumes'/'file_heights'), which are used by create_quantitative_tables
when quant_isotopes is set to 'volume' or 'height'.

Closes: https://github.com/PASTAQ-MS/PASTAQ/issues/5
"""

import pastaq


def test_peak_cluster_has_volumes_attribute():
    """Regression test: PeakCluster must expose 'volumes', not 'file_volumes'."""
    # find_peak_clusters returns PeakCluster objects; verify the attribute exists
    # on the class without needing to run a full pipeline by inspecting bindings.
    cluster_attrs = [attr for attr in dir(pastaq.PeakCluster) if not attr.startswith("_")]
    assert "volumes" in cluster_attrs, (
        "PeakCluster must have a 'volumes' attribute for quant_isotopes='volume' to work"
    )
    assert "file_volumes" not in cluster_attrs, (
        "PeakCluster must not have 'file_volumes' (use 'volumes' instead)"
    )


def test_peak_cluster_has_heights_attribute():
    """Regression test: PeakCluster must expose 'heights'."""
    cluster_attrs = [attr for attr in dir(pastaq.PeakCluster) if not attr.startswith("_")]
    assert "heights" in cluster_attrs, (
        "PeakCluster must have a 'heights' attribute for quant_isotopes='height' to work"
    )


def test_peak_cluster_expected_attributes():
    """Verify the complete set of expected PeakCluster attributes."""
    expected_attrs = {"id", "mz", "rt", "avg_height", "avg_volume", "heights", "volumes", "peak_ids"}
    cluster_attrs = set(attr for attr in dir(pastaq.PeakCluster) if not attr.startswith("_"))
    missing = expected_attrs - cluster_attrs
    assert not missing, f"PeakCluster is missing expected attributes: {missing}"


if __name__ == "__main__":
    test_peak_cluster_has_volumes_attribute()
    test_peak_cluster_has_heights_attribute()
    test_peak_cluster_expected_attributes()
    print("All PeakCluster attribute tests passed!")
