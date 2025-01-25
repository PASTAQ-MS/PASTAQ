#ifndef DE_BDAL_CPP_IO_TIMSDATA_H
#define DE_BDAL_CPP_IO_TIMSDATA_H

/** \file
 *
 * Definition of the public C "Mini API" for reading Bruker's TIMS raw-data format.
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
    /// Uses NoPressureCompensation.
    ///
    /// \param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// \param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
    ///
    BdalTimsdataDllSpec uint64_t tims_open (
        const char *analysis_directory_name,
        uint32_t use_recalibrated_state
        );

    typedef enum  {
        /// Apply no pressure compensation
        NoPressureCompensation = 0,

        /// If a reference pressure and pressure per frame is available,
        /// apply a pressure compensation using a single reference point from the data.
        /// All identical trafos stay identical this way.
        AnalyisGlobalPressureCompensation = 1,

        /// If a reference pressure and pressure per frame is available,
        /// apply a pressure compensation using the frame specific pressure. Each frame gets a different trafo.
        /// If a reference pressure is missing, but frame pressures are available,
        /// apply a pressure compensation using the frame specific pressure and a single reference point. Each frame gets a different trafo.
        PerFramePressureCompensation = 2,

        /// Deprecated, same behavior as PerFramePressureCompensation
        PerFramePressureCompensationWithMissingReference = 3
    } pressure_compensation_strategy;

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
    /// \param pressure_compensation_strategy the pressure compensation strategy
    ///
    BdalTimsdataDllSpec uint64_t tims_open_v2(
        const char* analysis_directory_name,
        uint32_t use_recalibrated_state,
        pressure_compensation_strategy pressure_compensation_strategy
    );

    /// Close data set.
    ///
    /// \param handle obtained by tims_open(); passing 0 is ok and has no effect.
    ///
    BdalTimsdataDllSpec void tims_close (uint64_t handle);

    /// Return the last error as a string (thread-local).
    ///
    /// \param buf pointer to a buffer into which the error string will be written.
    ///
    /// \param length length of the buffer
    ///
    /// \returns the actual length of the error message (including the final zero
    /// byte). If this is longer than the input parameter 'length', you know that the
    /// returned error string was truncated to fit in the provided buffer.
    ///
    BdalTimsdataDllSpec uint32_t tims_get_last_error_string (char *buf, uint32_t length);

    /// Returns 1 if the raw data have been recalibrated after acquisition, e.g. in the
    /// DataAnalysis software. Note that masses and 1/K0 values in the raw-data SQLite
    /// file are always in the raw calibration state, not the recalibrated state.
    ///
    BdalTimsdataDllSpec uint32_t tims_has_recalibrated_state (uint64_t handle);

    /// Read a range of scans from a single frame.
    ///
    /// Output layout: (N = scan_end - scan_begin = number of requested scans)
    ///   N x uint32_t: number of peaks in each of the N requested scans
    ///   N x (two uint32_t arrays: first indices, then intensities)
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error, otherwise the number of buffer bytes necessary for the output
    /// of this call (if this is larger than the provided buffer length, the result is not
    /// complete).
    ///
    BdalTimsdataDllSpec uint32_t tims_read_scans_v2 (
        uint64_t handle,
        int64_t frame_id,      //< from .tdf SQLite: Frames.Id
        uint32_t scan_begin,   //< first scan number to read (inclusive)
        uint32_t scan_end,     //< last scan number (exclusive)
        void *buf,             //< destination buffer allocated by user
        uint32_t length        //< length of buffer (in bytes)
        );

    /// Function type that takes a centroided peak list.
    typedef void(msms_spectrum_function)(
        int64_t id,                //< the id of the precursor or frame
        uint32_t num_peaks,        //< the number of peaks in the spectrum
        const double *mz_values,   //< all peak m/z values
        const float *area_values,  //< all peak areas
        void *user_data
    );

    /// Function type that takes a (non-centroided) profile spectrum.
    typedef void(msms_profile_spectrum_function)(
        int64_t id,                      //< the id of the precursor or frame
        uint32_t num_points,             //< number of entries in the profile spectrum
        const int32_t *intensity_values, //< the "quasi profile"
        void *user_data
    );

    /// Read peak-picked MS/MS spectra for a list of PASEF precursors.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor, performs centroiding using an algorithm and parameters
    /// suggested by Bruker, and returns the resulting MS/MS spectra (one for each
    /// precursor ID).
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_msms_v2(
        uint64_t handle,
        const int64_t *precursors,        //< list of PASEF precursor IDs; the returned spectra may be in different order
        uint32_t num_precursors,          //< number of requested spectra, must be >= 1
        msms_spectrum_function *callback, //< callback accepting the MS/MS spectra
        void *user_data                   //< will be passed to callback
    );

    /// Read peak-picked MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_msms_for_frame_v2(
        uint64_t handle,
        int64_t frame_id,                 //< frame id
        msms_spectrum_function *callback, //< callback accepting the MS/MS spectra
        void *user_data                   //< will be passed to callback
    );

    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor. These "quasi" profile spectra are passed back - one
    /// for each precursor ID.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_profile_msms_v2(
        uint64_t handle,
        const int64_t *precursors,                //< list of PASEF precursor IDs; the returned spectra may be in different order
        uint32_t num_precursors,                  //< number of requested spectra, must be >= 1
        msms_profile_spectrum_function *callback, //< callback accepting profile MS/MS spectra
        void *user_data                           //< will be passed to callback
    );

    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads for all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_profile_msms.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list. The parameter id in the callback is the
    /// precursor ID.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_profile_msms_for_frame_v2(
        uint64_t handle,
        int64_t frame_id,                         //< frame id
        msms_profile_spectrum_function *callback, //< callback accepting profile MS/MS spectra
        void *user_data                           //< will be passed to callback
    );

    /// Read peak-picked spectra for a tims frame.
    ///
    /// Given a frame ID, this function reads the frame, 
    /// sums up the corresponding scan-number ranges into a synthetic profile
    /// spectrum, performs centroiding using an algorithm and parameters
    /// suggested by Bruker, and returns the resulting spectrum (exactly one for 
    /// the frame ID).
    ///
    /// Note: Result callback identical to the tims_read_pasef_msms_v2 methods, but
    /// only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_extract_centroided_spectrum_for_frame_v2(
        uint64_t handle,
        int64_t frame_id,                 //< frame id
        uint32_t scan_begin,              //< first scan number to read (inclusive)
        uint32_t scan_end,                //< last scan number (exclusive)
        msms_spectrum_function *callback, //< callback accepting the spectra
        void *user_data                   //< will be passed to callback
    );

    /// Read peak-picked spectra for a tims frame with a custom peak picker resolution.
    ///
    /// Same as tims_extract_centroided_spectrum_for_frame_v2(), 
    /// but a user supplied resolution for the peak picker is applied.
    /// Can be used to prevent invalid split peaks in case of low ion statistics.
    /// The default suggested value in tims_extract_centroided_spectrum_for_frame_v2()
    /// is determined by the GlobalMetadata entry "PeakWidthEstimateValue" as
    /// 1 / PeakWidthEstimateValue for "PeakWidthEstimateType" = 1.
    ///
    /// Note: Result callback identical to the tims_read_pasef_msms_v2 methods, but
    /// only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_extract_centroided_spectrum_for_frame_ext(
        uint64_t handle,
        int64_t frame_id,                 /// frame id
        uint32_t scan_begin,              /// first scan number to read (inclusive)
        uint32_t scan_end,                /// last scan number (exclusive)
        double peakFinderResolution,      /// resolution for the peak picker, valid range is [1000.0 - 120,000.0]
        msms_spectrum_function* callback, /// callback accepting the spectra
        void* user_data                   /// will be passed to callback
    );

    /// Read "quasi profile" spectra for a tims frame.
    ///
    /// Given a frame ID, this function reads the frame, 
    /// and sums up the corresponding scan-number ranges into a synthetic profile
    /// spectrum. These "quasi" profile spectrum is passed back.
    ///
    /// Note: Result callback identical to the tims_read_pasef_profile_msms_v2 methods,
    /// but only returns a single result and the parameter id is the frame_id
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_extract_profile_for_frame(
        uint64_t handle,
        int64_t frame_id,                         //< frame id
        uint32_t scan_begin,                      //< first scan number to read (inclusive)
        uint32_t scan_end,                        //< last scan number (exclusive)
        msms_profile_spectrum_function *callback, //< callback accepting profile MS/MS spectra
        void *user_data                           //< will be passed to callback
    );

    /// -----------------------------------------------------------------------------------
    ///
    /// Extraction of chromatograms.
    ///
    /// -----------------------------------------------------------------------------------

    /// A chromatogram extraction job, i.e., the definition of a chromatogram trace. The value of a
    /// chromatogram trace point at a given time (i.e., for a given frame in the TDF) is determined
    /// by summing up the intensities of all peaks in that frame which fall into the specified m/z
    /// and 1/K0 window.
    typedef struct {
        int64_t id;          //< user-definable, will be passed to chromatogram_trace_sink
        double time_begin;
        double time_end;     //< (exclusive)
        double mz_min;
        double mz_max;
        double ook0_min;
        double ook0_max;
    } TimsChromatogramJob;
    
    /// A user-provided function called by the DLL whenever it is ready to process a new
    /// job. \returns 0 on error (chromatogram generation will then stop with an error), 1 when a
    /// new job has been produced, 2 when no more jobs are available
    typedef uint32_t(chromatogram_job_generator)(
        TimsChromatogramJob *,     //< members of pointed-to struct shall be filled by callback
        void *user_data            //< = whatever was passed to the main DLL call
    );

    /// Callback used to send a finished chromatogram traces from the DLL to the user. \returns 0 on
    /// error (chromatogram generation will then stop with an error), 1 when no error.
    typedef uint32_t(chromatogram_trace_sink)(
        int64_t id,                //< job ID
        uint32_t num_points,       //< number of points this chromatogram has
        const int64_t *frame_ids,  //< 'x values' = TDF frame IDs
        const uint64_t *values,    //< 'y values' = chromatogram trace
        void *user_data            //< = whatever was passed to the main DLL call
    );

    /// Extract several (MS1-only) chromatograms from an analysis.
    ///
    /// The DLL retrieves the jobs (i.e., the chromatogram definitions) from the specified generator
    /// function while iterating through the analysis. The jobs must be delivered in the order of
    /// ascending 'time_begin'.
    ///
    /// The DLL delivers chromatogram traces to the specified sink callback as soon as they are
    /// finished. When an error occurs, some of the jobs "pulled" so far might not be answered.
    ///
    /// \returns 0 on error
    BdalTimsdataDllSpec uint32_t tims_extract_chromatograms(
        uint64_t handle,
        chromatogram_job_generator *get_job,
        chromatogram_trace_sink *deliver_trace,
        void *user_data                           //< will be passed to callbacks
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
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_index_to_mz;
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_mz_to_index;

    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and 1/K0 values.
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_scannum_to_oneoverk0;
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_oneoverk0_to_scannum;

    /// mobility transformation: convert back and forth between (possibly non-integer)
    /// scan numbers and TIMS voltages.
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_scannum_to_voltage;
    BdalTimsdataDllSpec BdalTimsConversionFunction tims_voltage_to_scannum;
    
    /// Set the number of threads that this DLL is allowed to use internally. [The
    /// index<->m/z transformation is internally parallelized using OpenMP; this call is
    /// simply forwarded to omp_set_num_threads(). Has no effect on Linux].
    ///
    /// \param n number of threads to use (n must be >= 1).
    ///
    BdalTimsdataDllSpec void tims_set_num_threads (uint32_t n);

    /// Converts the 1/K0 value to CCS (in Angstrom^2) using the Mason-Shamp equation
    /// \param ook0 the 1/K0 value in Vs/cm2
    /// \param charge the charge
    /// \param mz the mz of the ion
    /// \returns the CCS value in Angstrom^2
    BdalTimsdataDllSpec double tims_oneoverk0_to_ccs_for_mz(
        const double ook0, 
        const int charge,
        const double mz
    );

    /// Converts the CCS (in Angstrom^2) to 1/K0 using the Mason-Shamp equation
    /// \param ccs the ccs value in Angstrom^2
    /// \param charge the charge
    /// \param mz the mz of the ion
    /// \returns the 1/K0 value in Vs/cm2
    BdalTimsdataDllSpec double tims_ccs_to_oneoverk0_for_mz(
        const double ccs, 
        const int charge,
        const double mz
    );


#ifdef __cplusplus
}
#endif

#endif //  DE_BDAL_CPP_IO_TIMSDATA_H

/* Local Variables:  */
/* mode: c           */
/* End:              */
