#ifndef UTILS_SEARCH_HPP
#define UTILS_SEARCH_HPP

#include <vector>

// This namespace contain functions to perform search on data structures.
namespace Search {

size_t lower_bound(const std::vector<double> &haystack, double needle);

// Generalize lower_bound search that uses a custom comparison fuction.
template <class T>
struct KeySort {
    size_t index;
    T sorting_key;
};
template <typename T>
size_t lower_bound(const std::vector<KeySort<T>> &haystack, T needle) {
    size_t index = 0;
    size_t l = 0;
    size_t r = haystack.size() - 1;
    while (l <= r) {
        index = (l + r) / 2;
        if (haystack[index].sorting_key < needle) {
            l = index + 1;
        } else if (haystack[index].sorting_key > needle) {
            r = index - 1;
        } else {
            break;
        }
        if (index == 0) {
            break;
        }
    }
    return index;
}

}  // namespace Search

#endif /* UTILS_SEARCH_HPP */
