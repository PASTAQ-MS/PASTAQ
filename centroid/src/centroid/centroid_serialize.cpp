#include "centroid/centroid_serialize.hpp"
#include "utils/serialization.hpp"

bool Centroid::Serialize::read_point(std::istream &stream,
                                     Centroid::Point *point) {
    Serialization::read_uint64(stream, &point->i);
    Serialization::read_uint64(stream, &point->j);
    Serialization::read_double(stream, &point->value);
    return stream.good();
}

bool Centroid::Serialize::write_point(std::ostream &stream,
                                      const Centroid::Point &point) {
    Serialization::write_uint64(stream, point.i);
    Serialization::write_uint64(stream, point.j);
    Serialization::write_double(stream, point.value);
    return stream.good();
}

bool Centroid::Serialize::read_peak(std::istream &stream,
                                    Centroid::Peak *peak) {
    Serialization::read_uint64(stream, &peak->i);
    Serialization::read_uint64(stream, &peak->j);
    Serialization::read_double(stream, &peak->mz);
    Serialization::read_double(stream, &peak->rt);
    Serialization::read_double(stream, &peak->height);
    Serialization::read_double(stream, &peak->total_intensity);
    Serialization::read_double(stream, &peak->sigma_mz);
    Serialization::read_double(stream, &peak->sigma_rt);
    Serialization::read_double(stream, &peak->mz_centroid);
    Serialization::read_double(stream, &peak->rt_centroid);
    Serialization::read_double(stream, &peak->height_centroid);
    Serialization::read_double(stream, &peak->total_intensity_centroid);
    Serialization::read_double(stream, &peak->border_background);

    // Centroid::Peak::points
    uint64_t points_size = 0;
    Serialization::read_uint64(stream, &points_size);
    peak->points.resize(points_size);
    for (auto &point : peak->points) {
        Centroid::Serialize::read_point(stream, &point);
    }

    // Centroid::Peak::boundary
    uint64_t boundary_size = 0;
    Serialization::read_uint64(stream, &boundary_size);
    peak->boundary.resize(boundary_size);
    for (auto &point : peak->boundary) {
        Centroid::Serialize::read_point(stream, &point);
    }

    return stream.good();
}

bool Centroid::Serialize::write_peak(std::ostream &stream,
                                     const Centroid::Peak &peak) {
    Serialization::write_uint64(stream, peak.i);
    Serialization::write_uint64(stream, peak.j);
    Serialization::write_double(stream, peak.mz);
    Serialization::write_double(stream, peak.rt);
    Serialization::write_double(stream, peak.height);
    Serialization::write_double(stream, peak.total_intensity);
    Serialization::write_double(stream, peak.sigma_mz);
    Serialization::write_double(stream, peak.sigma_rt);
    Serialization::write_double(stream, peak.mz_centroid);
    Serialization::write_double(stream, peak.rt_centroid);
    Serialization::write_double(stream, peak.height_centroid);
    Serialization::write_double(stream, peak.total_intensity_centroid);
    Serialization::write_double(stream, peak.border_background);

    // Centroid::Peak::points
    Serialization::write_uint64(stream, peak.points.size());
    for (const auto &point : peak.points) {
        Centroid::Serialize::write_point(stream, point);
    }

    // Centroid::Peak::boundary
    Serialization::write_uint64(stream, peak.boundary.size());
    for (const auto &point : peak.boundary) {
        Centroid::Serialize::write_point(stream, point);
    }

    return stream.good();
}
