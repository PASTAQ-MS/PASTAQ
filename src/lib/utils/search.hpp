#ifndef UTILS_SEARCH_HPP
#define UTILS_SEARCH_HPP

// This namespace contain functions to perform search on data structures.
namespace Search {

std::size_t lower_bound(const std::vector<double> &haystack, double needle);

// Generalize lower_bound search that uses a custom comparison fuction.
template <class T>
struct KeySort {
    std::size_t index;
    T sorting_key;
};
template <typename T>
std::size_t lower_bound(const std::vector<KeySort<T>> &haystack, T needle) {
    std::size_t index = 0;
    std::size_t l = 0;
    std::size_t r = haystack.size() - 1;
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
// TODO(alex): We probably want a generic lower_bound function that takes a
// predicate function and a generic type.
// size_t lower_bound(
// const std::vector<T1> &haystack, T2 needle, less_than, greater than)...

}  // namespace Search

#endif /* UTILS_SEARCH_HPP */
