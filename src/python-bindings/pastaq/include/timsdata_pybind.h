#ifndef TIMSDATACPP_PYBIND_H
#define TIMSDATACPP_PYBIND_H

#include <pybind11/pybind11.h>
#include "timsdatacpp.h"  // Include core C++ class

namespace py = pybind11;

void register_tims_data(py::module_ &m);

#endif // TIMSDATACPP_PYBIND_H