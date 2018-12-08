#include <iostream>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

int add(int i, int j) { return i + j; }

int vector_add(std::vector<int> v) {
    int acc = 0;
    for (const auto &e : v) {
        acc += e;
    }
    return acc;
}

void say_hello() { std::cout << "hello world!" << std::endl; }

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin";

    m.def("add", &add, "A function which adds two numbers");
    m.def("vector_add", &vector_add,
          "A function which adds all elements in the given vector");
    m.def("hello", &say_hello, "Typical hello world example");
}
