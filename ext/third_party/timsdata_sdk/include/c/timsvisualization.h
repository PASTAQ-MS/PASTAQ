#ifndef DE_BDAL_CPP_IO_TIMSVISUALIZATION_H
#define DE_BDAL_CPP_IO_TIMSVISUALIZATION_H

/** \file
 *
 * Definition of the public C "Mini API" for visualization of Bruker's TIMS raw-data.
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
    /// On failure, returns 0, and you can use tims_vis_get_last_error_string() to obtain a
    /// string describing the problem.
    ///
    /// \param analysis_directory_name the name of the directory in the file system that
    /// contains the analysis data, in UTF-8 encoding.
    ///
    /// \param use_recalibrated_state if non-zero, use the most recent recalibrated state
    /// of the analysis, if there is one; if zero, use the original "raw" calibration
    /// written during acquisition time.
    ///
    BdalTimsdataDllSpec uint64_t tims_vis_open(
        const char *analysis_directory_name,
        uint32_t use_recalibrated_state
        );

    /// Close data set.
    ///
    /// \param handle obtained by tims_vis_open(); passing 0 is ok and has no effect.
    ///
    BdalTimsdataDllSpec void tims_vis_close(uint64_t handle);

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
    BdalTimsdataDllSpec uint32_t tims_vis_get_last_error_string(char *buf, uint32_t len);

    /// Defines a filter region in the data, only the data in the given ranges will be part of the visualizations
    struct tims_vis_extraction_filter {
        double rtSecondsLower;
        double rtSecondUpper;

        double ook0Lower;
        double ook0Upper;

        double mzLower;
        double mzUpper;
    };

    /// Definition of heatmap dimensions in pixels
    struct tims_vis_heatmap_sizes {
        // X: retention time, Y: mobility
        int32_t widthRtMob;
        int32_t heightRtMob;

        // X: retention time, Y: M/Z
        int32_t widthRtMz;
        int32_t heightRtMz;

        // X: mobility, Y: M/Z
        int32_t widthMobMz;
        int32_t heightMobMz;
    };

    /// Start asynchronous calculations for the specified data filter and heatmap sizes.
    /// Cancels any currently running calculations.
    /// You can get intermediate line plots and heatmaps at any time during computation.
    /// Use tims_vis_get_state, tims_vis_wait or tims_vis_wait_complete to check/wait for
    /// the computation to be done.
    ///
    /// \param handle obtained by tims_open()
    ///
    /// \param filter the filter use for visualization
    ///
    /// \param heatmap_sizes the size definition of the 3 heatmaps
    ///
    /// \returns 1 on success, 0 on error
    ///
    BdalTimsdataDllSpec uint64_t tims_vis_calculate_async(
        uint64_t handle,
        tims_vis_extraction_filter filter,
        tims_vis_heatmap_sizes heatmap_sizes);

    /// Cancel a running calculation
    ///
    /// \param handle obtained by tims_open()
    BdalTimsdataDllSpec uint64_t tims_vis_cancel(
        uint64_t handle);

    /// Get the state of the last started computation.
    ///
    /// \param handle obtained by tims_open()
    /// \param[out] job_id the job id of the last computation (just an increasing counter)
    /// \param[out] progress the progress of the last job
    /// \param[out] complete if the last job has been completed
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_state(uint64_t handle, uint64_t& job_id, float& progress, bool& complete);

    /// Wait until the current job is finished or the timeout in milliseconds is elapsed
    ///
    /// \param handle obtained by tims_open()
    /// \param tims_in_ms the maximum wait time in milliseconds
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_wait(uint64_t handle, uint32_t time_in_ms);

    /// Wait until the current job is finished
    ///
    /// \param handle obtained by tims_open()
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_wait_complete(uint64_t handle);

    /// Possible transformations
    enum tims_vis_transformation {
        /// no transformation
        NONE,
        /// square root
        SQRT,
        /// logarithm for values > 1 and linear extrapolation for values < 1
        LIMITED_LOG
    };

    /// A line in the plot in pixel coordinates.
    /// spectrum, mobilogram and chromatogram plots consist of an array of lines
    struct tims_vis_line {
        int32_t x0;
        int32_t y0;
        int32_t x1;
        int32_t y1;
    };

    /// Get the last computed spectrum line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param pixelX the requested X dimension of the plot in pixels
    /// \param pixelY the requested Y dimension of the plot in pixels
    /// \param transformation a transformation that might be applied to the plot
    /// \param line_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the line_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_spectrum_line_plot(
        uint64_t handle,
        int32_t pixelX, int32_t pixelY,
        tims_vis_transformation transformation,
        tims_vis_line* line_array, uint32_t size);

    /// Get the last computed mobilogram line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param pixelX the requested X dimension of the plot in pixels
    /// \param pixelY the requested Y dimension of the plot in pixels
    /// \param transformation a transformation that might be applied to the plot
    /// \param line_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the line_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_mobilogram_line_plot(
        uint64_t handle,
        int32_t pixelX, int32_t pixelY,
        tims_vis_transformation transformation,
        tims_vis_line* line_array, uint32_t length);

    /// Get the last computed chromatogram line plot (may be an intermediate result).
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param pixelX the requested X dimension of the plot in pixels
    /// \param pixelY the requested Y dimension of the plot in pixels
    /// \param transformation a transformation that might be applied to the plot
    /// \param line_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the line_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_chromatogram_line_plot(
        uint64_t handle,
        int32_t pixelX, int32_t pixelY,
        tims_vis_transformation transformation,
        tims_vis_line* line_array, uint32_t length);

    /// Get the last computed (maybe intermediate) heatmap
    /// with 1/K0 on the X axis and mz on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param transformation a transformation that might be applied to the plot
    /// \param image_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the image_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_image_mob_mz (
        uint64_t handle,
        tims_vis_transformation transformation,
        float* image_array,
        uint32_t length
    );

    /// Get the last computed (maybe intermediate) heatmap
    /// with retention time in seconds on the X axis and mz on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param transformation a transformation that might be applied to the plot
    /// \param image_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the image_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_image_rt_mz (
        uint64_t handle,
        tims_vis_transformation transformation,
        float* image_array,
        uint32_t length
    );

    /// Get the last computed (maybe intermediate) heatmap
    /// with retention time in seconds on the X axis and 1/Ko on the Y axis.
    /// The dimension of the heatmap is defined in the tims_vis_calculate_async call.
    /// If the provided array is too small for the result data, no data is written.
    ///
    /// \param handle obtained by tims_open()
    /// \param transformation a transformation that might be applied to the plot
    /// \param image_array the resulting array of tims_vis_line, allocated by user
    /// \param size size of the provided array
    ///
    /// \returns 0 on error, otherwise the required size of the image_array
    ///
    BdalTimsdataDllSpec uint32_t tims_vis_get_image_rt_mob (
        uint64_t handle,
        tims_vis_transformation transformation,
        float* image_array,
        uint32_t length
    );

#ifdef __cplusplus
}
#endif

#endif

/* Local Variables:  */
/* mode: c           */
/* End:              */
