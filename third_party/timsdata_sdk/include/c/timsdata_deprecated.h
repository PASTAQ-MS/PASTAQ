#ifndef DE_BDAL_CPP_IO_TIMSDATA_DEPRECATED_H
#define DE_BDAL_CPP_IO_TIMSDATA_DEPRECATED_H

/** \file
 *
 * Definition of the public C "Mini API" for reading Bruker's TIMS raw-data format.
 *
 * DEPRECATED FUNCTIONS.
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

    /// The type of the callback functor required by tims_read_pasef_msms.
    ///
    /// \param[in] precursor_id the id of the precursor
    /// \param[in] num_peaks the number of peaks in the MS/MS spectrum
    /// \param[in] mz_values all peak m/z values as double
    /// \param[in] area_values all peak areas as float
    ///
    typedef void(*msms_spectrum_functor)(
        int64_t precursor_id,
        uint32_t num_peaks,
        double *mz_values,
        float *area_values
    );

    /// The type of the callback functor required by tims_read_pasef_profile_msms.
    ///
    /// \param[in] precursor_id the id of the precursor
    /// \param[in] num_peaks the number of peaks in the MS/MS spectrum
    /// \param[in] intensity_values as integer of "quasi profile"
    ///
    typedef void(*msms_profile_spectrum_functor)(
        int64_t precursor_id,
        uint32_t num_points,
        int32_t *intensity_values
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
    /// order in the specified precursor ID list.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_msms(
        uint64_t handle,
        int64_t *precursors,     //< list of PASEF precursor IDs; the returned spectra may be in different order
        uint32_t num_precursors, //< number of requested spectra, must be >= 1
        msms_spectrum_functor callback //< callback accepting the MS/MS spectra
    );

    /// Read peak-picked MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_msms.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_msms_for_frame(
        uint64_t handle,
        int64_t frame_id,     //< frame id
        msms_spectrum_functor callback //< callback accepting the MS/MS spectra
    );

    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a list of PASEF precursor IDs, this function reads all necessary PASEF
    /// frames, sums up the corresponding scan-number ranges into synthetic profile
    /// spectra for each precursor. These "quasi" profile spectra are passed back - one
    /// for each precursor ID.
    ///
    /// Note: the order of the returned MS/MS spectra does not necessarily match the
    /// order in the specified precursor ID list.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_profile_msms(
        uint64_t handle,
        int64_t *precursors,     //< list of PASEF precursor IDs; the returned spectra may be in different order
        uint32_t num_precursors, //< number of requested spectra, must be >= 1
        msms_profile_spectrum_functor callback //< callback accepting profile MS/MS spectra
    );

    /// Read "quasi profile" MS/MS spectra for all PASEF precursors from a given frame.
    ///
    /// Given a frame id, this function reads all contained PASEF precursors the necessary PASEF
    /// frames in the same way as tims_read_pasef_profile_msms.
    ///
    /// Note: different threads must not read scans from the same storage handle
    /// concurrently.
    ///
    /// \returns 0 on error
    ///
    BdalTimsdataDllSpec uint32_t tims_read_pasef_profile_msms_for_frame(
        uint64_t handle,
        int64_t frame_id,     //< frame id
        msms_profile_spectrum_functor callback //< callback accepting profile MS/MS spectra
    );
    
#ifdef __cplusplus
}
#endif

#endif //  DE_BDAL_CPP_IO_TIMSDATA_DEPRECATED_H

/* Local Variables:  */
/* mode: c           */
/* End:              */
