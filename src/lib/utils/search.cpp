#include "utils/search.hpp"

size_t Search::lower_bound(const std::vector<double> &haystack, double needle) {
    size_t index = 0;
    size_t l = 0;
    size_t r = haystack.size() - 1;
    while (l <= r) {
        index = (l + r) / 2;
        if (haystack[index] < needle) {
            l = index + 1;
        } else if (haystack[index] > needle) {
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
