"""
Regression test for safe Peak default construction.

Validates that pastaq.Peak() creates instances with all scalar members
value-initialized (zero for numeric types, False for bool), not uninitialized garbage.
This test ensures the pybind11 lambda constructor is working correctly.
"""

import pastaq


def test_peak_default_initialization():
    """Test that Peak() initializes all numeric fields to zero."""
    peak = pastaq.Peak()
    
    # Integer fields
    assert peak.id == 0, "id should be initialized to 0"
    assert peak.raw_roi_num_points == 0, "raw_roi_num_points should be initialized to 0"
    assert peak.raw_roi_num_scans == 0, "raw_roi_num_scans should be initialized to 0"
    assert peak.fit_failure_code == 0, "fit_failure_code should be initialized to 0"
    
    # Double fields - all should be 0.0
    double_fields = [
        'local_max_mz', 'local_max_rt', 'local_max_height',
        'rt_delta',
        'roi_min_mz', 'roi_max_mz', 'roi_min_rt', 'roi_max_rt',
        'raw_roi_mean_mz', 'raw_roi_mean_rt',
        'raw_roi_sigma_mz', 'raw_roi_sigma_rt',
        'raw_roi_skewness_mz', 'raw_roi_skewness_rt',
        'raw_roi_kurtosis_mz', 'raw_roi_kurtosis_rt',
        'raw_roi_max_height', 'raw_roi_total_intensity',
        'fitted_height', 'fitted_mz', 'fitted_rt',
        'fitted_sigma_mz', 'fitted_sigma_rt', 'fitted_volume'
    ]
    
    for field in double_fields:
        value = getattr(peak, field)
        assert value == 0.0, f"{field} should be initialized to 0.0, got {value}"
    
    # Boolean field
    assert peak.peak_fit_failure == False, "peak_fit_failure should be initialized to False"


def test_peak_repr_with_default_initialization():
    """Test that Peak().__repr__() does not crash with indeterminate values."""
    peak = pastaq.Peak()
    
    # This should not raise an exception or produce NaN/inf values
    repr_str = repr(peak)
    assert "Peak" in repr_str, "Peak repr should contain 'Peak'"
    assert "id: 0" in repr_str, "Peak repr should show id: 0"
    
    # Verify no NaN or inf appear (which would indicate uninitialized garbage)
    assert "nan" not in repr_str.lower(), "repr should not contain NaN"
    assert "inf" not in repr_str.lower(), "repr should not contain inf"


if __name__ == "__main__":
    test_peak_default_initialization()
    test_peak_repr_with_default_initialization()
    print("All Peak initialization tests passed!")
