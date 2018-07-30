#include <cmath>
#include <iostream>
#include <limits>

#include "warp2d/warp2d.hpp"

typedef struct FU {
    /* Used for dynamic programming: f - cummulative function value, u -
     * position of the optimum predecessor */
    double f;
    int u;
} FU;

double similarity_func() { return 1.0; }

int Warp2D::cow_2D(std::vector<Centroid::Peak>& target_peaks,
                   std::vector<Centroid::Peak>& source_peaks, int sample_length,
                   int segment_length, int slack) {
    // The number of segments.
    int N = sample_length / segment_length;

    // Minimum time step.
    // DEBUG: Hardcoding values here.
    double rt_min = 0.0;
    double rt_max = 40.0;
    double delta_rt = (rt_max - rt_min) / (double)(sample_length - 1);

    // Allocate data structures.
    int *xstart, *xlength, *xend;
    FU **nodes, *node, *next_node;
    xstart = (int*)malloc((N + 1) * sizeof(int));
    xlength = (int*)malloc((N + 1) * sizeof(int));
    xend = (int*)malloc((N + 1) * sizeof(int));
    nodes = (FU**)malloc((N + 1) * sizeof(FU*));

    int m = segment_length;
    int t = slack;
    int Lt = sample_length;
    for (int i = 0; i <= N; i++) {
        // Compute range of possible movement for each of the (N+1) nodes.
        // Nodes 0 and N do not move, no warping possible at start or end.
        xstart[i] = std::max((i * (m - t)), (Lt - (N - i) * (m + t)));
        xend[i] = std::min((i * (m + t)), (Lt - (N - i) * (m - t)));
        xlength[i] = xend[i] - xstart[i] + 1;
        // i-th row of F,U matrix
        nodes[i] = (FU*)malloc(xlength[i] * sizeof(FU));
        // assert(nodes[i]);
    }

    // Initialize F,U matrix.
    for (int i = 0; i <= N; i++) {
        node = nodes[i];
        for (int j = 0; j < xlength[i]; j++, node++) {
            node->f = -std::numeric_limits<double>::max();
            node->u = -1;
        }
    }

    //[>
    // Find the best warping using dynamic programming
    //*/

    // initialize the only node at level N
    int level = N;
    node = nodes[level];
    node->f = 0.0;
    node->u = 0;

    // back-propagate F values
    // go through levels backwards
    // at each level, find range of allowed segments and find correlation
    for (level = N; level > 0; level--) {
        // double refTimeSegmt = delta_rt * m;
        // double refTimeStart = rt_min + (level - 1) * refTimeSegmt;
        // std::vector<Peak> target_peaksd = refWarpDB->getPeaksAtRTBand(
        // refTimeStart + refTimeSegmt / 2, refTimeSegmt / 2);
        // std::vector<Peak>& target_peaks = target_peaksd;
        std::cout << "level: " << level << std::endl;
        for (int i = 0; i < xlength[level]; i++) {
            int xi, xjmin, xjmax, jmin, jmax;
            FU* nodei;
            // Position of the node at level.
            xi = xstart[level] + i;    // x position of the i-th node at level
            nodei = nodes[level] + i;  // ith node at level
            // The position of a connected node at level-1 is constrained by the
            // position of the node at level and the window size and slack
            // constraints. find all connected nodes at level-1
            xjmin = std::max(xi - m - t, xstart[level - 1]);
            xjmax = std::min(xi - m + t, xend[level - 1]);
            jmin = xjmin - xstart[level - 1];
            jmax = xjmax - xstart[level - 1];
            std::cout << "xi: " << xi << std::endl;
            std::cout << "xstart: " << xjmin << std::endl;
            std::cout << "xend: " << xjmax << std::endl;
            for (int j = jmin; j <= jmax; j++) {
                int xj;
                FU* nodej;
                double d, sum;
                nodej =
                    nodes[level - 1] + j;  // jth node at level-1 to be filled
                                           // in and point to i if best
                // compute distance between connected nodes
                // x position of the j-th node at level-1
                xj = xstart[level - 1] + j;
                double smpTimeStart = rt_min + xj * delta_rt;
                double smpTimeSegmt = (xi - xj) * delta_rt;

                // this makes a copy of the peaks so the time can be changed in
                // the samples
                // std::vector<Peak> smpPeakd = smpWarpDB->getPeaksAtRTBand(
                // smpTimeStart + smpTimeSegmt / 2, smpTimeSegmt / 2);
                // std::vector<Peak>& source_peaks = smpPeakd;

                // d = similarity2D(refTimeStart, refTimeSegmt, target_peaks,
                // smpTimeStart, smpTimeSegmt, source_peaks);
                d = similarity_func();

                sum = d + nodei->f;
                if (sum > nodej->f) {
                    nodej->f = sum;
                    nodej->u = i;
                }
            }

            //// now need to check for ties and choose the one closest to the
            //// middle
            // int midj = (jmin + jmax) / 2;
            // int bestj = midj;
            // int bestjdist = 10000000;
            // for (int j = jmin; j <= jmax; j++) {
            // FU* nodej;
            // nodej =
            // nodes[level - 1] + j;  // jth node at level-1 to be filled
            //// in and point to i if best
            //}
        }
    }

    // printf("\n");
    //// forward-propagate U values, compute correlation/segment
    // int* warping = (int*)malloc((N + 1) * sizeof(int));
    // assert(warping);

    //// here is where the nodes are walked back
    // warping[0] = 0;
    // int i = 0;
    // totalCorr = 0;
    // for (int level = 0; level < N; level++) {
    // node = nodes[level] + i;
    // i = node->u;
    // warping[level + 1] = xstart[level + 1] + i;
    // next_node = nodes[level + 1] + i;
    // double sgCorr = node->f - next_node->f;
    // totalCorr += sgCorr;
    // segmtCorr.push_back(sgCorr);
    //}
    // printf("cow_2D()> Total Correlation Warped: %e\n", totalCorr);

    // double sampleTime = rt_min;
    // double refrncTime = rt_min;
    // origTime.push_back(sampleTime);
    // warpedTime.push_back(refrncTime);
    // for (i = 0; i < N; i++) {
    // sampleTime += (warping[i + 1] - warping[i]) * delta_rt;
    // refrncTime += m * delta_rt;
    // origTime.push_back(sampleTime);
    // warpedTime.push_back(refrncTime);
    //}

    //[>
    // Free data structures
    //*/
    for (int i = 0; i <= N; i++) free(nodes[i]);
    free(nodes);
    free(xend);
    free(xlength);
    free(xstart);
    // free(warping);

    return (1);
}

int Warp2D::warp_peaks(std::vector<Centroid::Peak>& target_peaks,
                       std::vector<Centroid::Peak>& source_peaks,
                       int sample_length, int segment_length, int slack) {
    // The number of segments.
    int N = sample_length / segment_length;

    // Minimum time step.
    // DEBUG: Hardcoding values here.
    double rt_min = 0.0;
    double rt_max = 40.0;
    double delta_rt = (rt_max - rt_min) / (double)(sample_length - 1);

    // Allocate data structures.
    int *xstart, *xlength, *xend;
    FU **nodes, *node, *next_node;
    xstart = (int*)malloc((N + 1) * sizeof(int));
    xlength = (int*)malloc((N + 1) * sizeof(int));
    xend = (int*)malloc((N + 1) * sizeof(int));
    nodes = (FU**)malloc((N + 1) * sizeof(FU*));

    int m = segment_length;
    int t = slack;
    int Lt = sample_length;
    for (int i = 0; i <= N; i++) {
        // Compute range of possible movement for each of the (N+1) nodes.
        // Nodes 0 and N do not move, no warping possible at start or end.
        xstart[i] = std::max((i * (m - t)), (Lt - (N - i) * (m + t)));
        xend[i] = std::min((i * (m + t)), (Lt - (N - i) * (m - t)));
        xlength[i] = xend[i] - xstart[i] + 1;
        // i-th row of F,U matrix
        nodes[i] = (FU*)malloc(xlength[i] * sizeof(FU));
        // assert(nodes[i]);
    }

    // Initialize F,U matrix.
    for (int i = 0; i <= N; i++) {
        node = nodes[i];
        for (int j = 0; j < xlength[i]; j++, node++) {
            node->f = -std::numeric_limits<double>::max();
            node->u = -1;
        }
    }

    //[>
    // Find the best warping using dynamic programming
    //*/

    // initialize the only node at level N
    int level = N;
    node = nodes[level];
    node->f = 0.0;
    node->u = 0;

    // back-propagate F values
    // go through levels backwards
    // at each level, find range of allowed segments and find correlation
    for (level = N; level > 0; level--) {
        // double refTimeSegmt = delta_rt * m;
        // double refTimeStart = rt_min + (level - 1) * refTimeSegmt;
        // std::vector<Peak> target_peaksd = refWarpDB->getPeaksAtRTBand(
        // refTimeStart + refTimeSegmt / 2, refTimeSegmt / 2);
        // std::vector<Peak>& target_peaks = target_peaksd;
        std::cout << "level: " << level << std::endl;
        for (int i = 0; i < xlength[level]; i++) {
            int xi, xjmin, xjmax, jmin, jmax;
            FU* nodei;
            // Position of the node at level.
            xi = xstart[level] + i;    // x position of the i-th node at level
            nodei = nodes[level] + i;  // ith node at level
            // The position of a connected node at level-1 is constrained by the
            // position of the node at level and the window size and slack
            // constraints. find all connected nodes at level-1
            xjmin = std::max(xi - m - t, xstart[level - 1]);
            xjmax = std::min(xi - m + t, xend[level - 1]);
            jmin = xjmin - xstart[level - 1];
            jmax = xjmax - xstart[level - 1];
            std::cout << "xi: " << xi << std::endl;
            std::cout << "xstart: " << xjmin << std::endl;
            std::cout << "xend: " << xjmax << std::endl;
            for (int j = jmin; j <= jmax; j++) {
                int xj;
                FU* nodej;
                double d, sum;
                nodej =
                    nodes[level - 1] + j;  // jth node at level-1 to be filled
                                           // in and point to i if best
                // compute distance between connected nodes
                // x position of the j-th node at level-1
                xj = xstart[level - 1] + j;
                double smpTimeStart = rt_min + xj * delta_rt;
                double smpTimeSegmt = (xi - xj) * delta_rt;

                // this makes a copy of the peaks so the time can be changed in
                // the samples
                // std::vector<Peak> smpPeakd = smpWarpDB->getPeaksAtRTBand(
                // smpTimeStart + smpTimeSegmt / 2, smpTimeSegmt / 2);
                // std::vector<Peak>& source_peaks = smpPeakd;

                // d = similarity2D(refTimeStart, refTimeSegmt, target_peaks,
                // smpTimeStart, smpTimeSegmt, source_peaks);
                d = similarity_func();

                sum = d + nodei->f;
                if (sum > nodej->f) {
                    nodej->f = sum;
                    nodej->u = i;
                }
            }

            //// now need to check for ties and choose the one closest to the
            //// middle
            // int midj = (jmin + jmax) / 2;
            // int bestj = midj;
            // int bestjdist = 10000000;
            // for (int j = jmin; j <= jmax; j++) {
            // FU* nodej;
            // nodej =
            // nodes[level - 1] + j;  // jth node at level-1 to be filled
            //// in and point to i if best
            //}
        }
    }

    // printf("\n");
    //// forward-propagate U values, compute correlation/segment
    // int* warping = (int*)malloc((N + 1) * sizeof(int));
    // assert(warping);

    //// here is where the nodes are walked back
    // warping[0] = 0;
    // int i = 0;
    // totalCorr = 0;
    // for (int level = 0; level < N; level++) {
    // node = nodes[level] + i;
    // i = node->u;
    // warping[level + 1] = xstart[level + 1] + i;
    // next_node = nodes[level + 1] + i;
    // double sgCorr = node->f - next_node->f;
    // totalCorr += sgCorr;
    // segmtCorr.push_back(sgCorr);
    //}
    // printf("cow_2D()> Total Correlation Warped: %e\n", totalCorr);

    // double sampleTime = rt_min;
    // double refrncTime = rt_min;
    // origTime.push_back(sampleTime);
    // warpedTime.push_back(refrncTime);
    // for (i = 0; i < N; i++) {
    // sampleTime += (warping[i + 1] - warping[i]) * delta_rt;
    // refrncTime += m * delta_rt;
    // origTime.push_back(sampleTime);
    // warpedTime.push_back(refrncTime);
    //}

    //[>
    // Free data structures
    //*/
    for (int i = 0; i <= N; i++) free(nodes[i]);
    free(nodes);
    free(xend);
    free(xlength);
    free(xstart);
    // free(warping);

    return (1);
}
