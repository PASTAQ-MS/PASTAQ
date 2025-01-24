#ifndef DE_BDAL_CPP_IO_TSFDATA_H
#define DE_BDAL_CPP_IO_TSFDATA_H

/** \file
 *
 * Definition of the public C "Mini API" for reading Bruker's timstof spectra raw-data format.
 *
 */

#include "configuration/visibility_decl.h"

#ifdef DE_BDAL_CPP_IO_TIMSDATA_BUILDING_DLL
  #define BdalTimsdataDllSpec BDAL_VIS_EXPORT_DECL
#else
  #define BdalTimsdataDllSpec BDAL_VIS_IMPORT_DECL
#endif

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

    /// Open data set.
    ///
    /// On success, returns a non-zero instance handle that needs to be passed to
    /// subsequent API calls, in particular to the required call to tims_close().
    ///
    /// On failure, returns 0, and you can use tims_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// \param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// \param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
    ///
    BdalTimsdataDllSpec uint64_t tsf_open (
        const char *analysis_directory_name,
        uint32_t use_recalibrated_state
        );

    /// Close data set.
    ///
    /// \param handle obtained by tsf_open(); passing 0 is ok and has no effect.
    ///
    BdalTimsdataDllSpec void tsf_close (uint64_t handle);

    /// Return the last error as a string (thread-local).
    ///
    /// \param buf pointer to a buffer into which the error string will be written.
    ///
    /// \param len length of the buffer
    ///
    /// \returns the actual length of the error message (including the final zero
    /// byte). If this is longer than the input parameter 'len', you know that the
    /// returned error string was truncated to fit in the provided buffer.
    ///
    BdalTimsdataDllSpec uint32_t tsf_get_last_error_string (char *buf, uint32_t len);

    /// Returns 1 if the raw data have been recalibrated after acquisition, e.g. in the
    /// DataAnalysis software. Note that masses and 1/K0 values in the raw-data SQLite
    /// file are always in the raw calibration state, not the recalibrated state.
    ///
    BdalTimsdataDllSpec uint32_t tsf_has_recalibrated_state (uint64_t handle);

    /// Read a line spectrum. Fails if no line spectrum is contained (GlobalMetatdata.HasLineSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// \param[in] handle the handle used for reading
    /// \param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// \param[out] index_array index values as double array 
    /// \param[out] intensity_array intensity values as float array 
    /// \param[in] length the length of the provided arrays
    ///
    /// \returns -1 on error, otherwise the number of entries necessary for the output arrays
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
    ///
    BdalTimsdataDllSpec int32_t tsf_read_line_spectrum_v2(
        uint64_t handle,
        int64_t spectrum_id,
        double* index_array,
        float* intensity_array,
        int32_t length
    );

    /// Read a line spectrum. Fails if no line spectrum or no peak width is contained (GlobalMetatdata.HasLineSpectra == 0 or GlobalMetatdata.HasLineSpectraPeakWidth == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// \param[in] handle the handle used for reading
    /// \param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// \param[out] index_array index values as double array 
    /// \param[out] intensity_array intensity values as float array 
    /// \param[out] width_array width values as float array 
    /// \param[in] length the length of the provided arrays
    ///
    /// \returns -1 on error, otherwise the number of entries necessary for the output arrays
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
    ///
    BdalTimsdataDllSpec int32_t tsf_read_line_spectrum_with_width_v2(
        uint64_t handle,
        int64_t spectrum_id,
        double* index_array,
        float* intensity_array,
        float* width_array,
        int32_t length
    );

    /// Read a profile spectrum. Fails if no profile is contained (GlobalMetatdata.HasProfileSpectra == 0)
    ///
    /// Note: different threads must not read spectra from the same storage handle
    /// concurrently.
    ///
    /// \param[in] handle the handle used for reading
    /// \param[in] spectrum_id from .tsf SQLite: Frames.Id
    /// \param[out] profile_array intensity values as uint32_t array, position in the array is the index
    /// \param[in] length the length of the provided array
    ///
    /// \returns -1 on error, otherwise the number of entries necessary for the output array
    /// of this call (if this is larger than the provided output array length, the result is not
    /// complete).
    ///
    BdalTimsdataDllSpec int32_t tsf_read_profile_spectrum_v2(
        uint64_t handle,
        int64_t spectrum_id,   
        uint32_t *profile_array,
        int32_t length
    );
    
    /// -----------------------------------------------------------------------------------
    ///
    /// Conversion functions coming up. All these functions share the same signature (see
    /// typedef 'BdalTimsConversionFunction'). They all return 1 on success, 0 on failure.
    ///
    /// -----------------------------------------------------------------------------------
    
    /// A function that transforms every value of the input array 'in' to a corresponding value in
    /// the output array 'out'. How many values it transforms is specified by the last argument. The
    /// individual transformations are independent of each other.
    typedef uint32_t BdalTimsConversionFunction (
        uint64_t handle,
        int64_t frame_id,      //< from .tdf SQLite: Frames.Id
        const double *in,      //< input array of values
        double *out,           //< output array of values
        uint32_t cnt           //< number of values to convert (arrays must have corresponding size)
        );
    
    /// m/z transformation: convert back and forth between (possibly non-integer) index
    /// values and m/z values.
    BdalTimsdataDllSpec BdalTimsConversionFunction tsf_index_to_mz;
    BdalTimsdataDllSpec BdalTimsConversionFunction tsf_mz_to_index;
    
    /// Set the number of threads that this DLL is allowed to use internally. [The
    /// index<->m/z transformation is internally parallelized using OpenMP; this call is
    /// simply forwarded to omp_set_num_threads()].
    ///
    /// \param n number of threads to use (n must be >= 1).
    ///
    BdalTimsdataDllSpec void tsf_set_num_threads (uint32_t n);

#ifdef __cplusplus
}
#endif

#endif //  DE_BDAL_CPP_IO_TIMSDATA_H

/* Local Variables:  */
/* mode: c           */
/* End:              */
