#include "warp2d/warp2d_serialize.hpp"
#include "utils/serialization.hpp"

bool Warp2D::Serialize::read_time_map(std::istream &stream,
                                      Warp2D::TimeMap *time_map) {
    Serialization::read_uint64(stream, &time_map->num_segments);
    Serialization::read_double(stream, &time_map->rt_min);
    Serialization::read_double(stream, &time_map->rt_max);
    Serialization::read_vector<double>(stream, &time_map->rt_start,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &time_map->rt_end,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &time_map->sample_rt_start,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &time_map->sample_rt_end,
                                       Serialization::read_double);
    return stream.good();
}

bool Warp2D::Serialize::write_time_map(std::ostream &stream,
                                       const Warp2D::TimeMap &time_map) {
    Serialization::write_uint64(stream, time_map.num_segments);
    Serialization::write_double(stream, time_map.rt_min);
    Serialization::write_double(stream, time_map.rt_max);
    Serialization::write_vector<double>(stream, time_map.rt_start,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, time_map.rt_end,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, time_map.sample_rt_start,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, time_map.sample_rt_end,
                                        Serialization::write_double);
    return stream.good();
}
