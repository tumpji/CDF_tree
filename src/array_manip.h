#if not defined INCLUDED_ARRAY_MANIPULATION
#define INCLUDED_ARRAY_MANIPULATION

#include <memory>
#include <utility>
#include <array>
#include <iterator>
#include <boost/assert.hpp>

///////////////////////// UTILS /////////////////////////////
namespace utils {


template<class T, class M>
unsigned lower_or_equal_bound (const T* array, unsigned size, const M& element) {
    // returns index of element which is greater or equal than provided reference
    // output is 0-size (including borders)
    unsigned step, count = size, it = 0;
    while (count > 0) {
        step = count / 2;

        if (array[it + step] <= element) {
            it = it + step + 1;
            count -= step + 1;
        } else {
            count = step;
        }
    }
    return it;
}

template<class T, class M>
unsigned lower_bound(const T* array, unsigned size, const M& element) {
    // returns index of element which is greater or equal than provided reference
    // output is 0-size (including borders)
    unsigned step, count = size, it = 0;
    while (count > 0) {
        step = count / 2;

        if (array[it + step] < element) {
            it = it + step + 1;
            count -= step + 1;
        } else if (array[it + step] == element) {
            return it + step;
        } else {
            count = step;
        }
    }
    return it;
}

template<class T, class M>
int binary_search(const T* array, const unsigned size, const M& element) {
    // returns index of element which is greater or equal than provided reference
    // output is 0-size (including borders)
    int low = 0, up = size-1;

    while (low <= up) {
        int pivot = (up + low)/2;
        BOOST_ASSERT(pivot >= 0);
        BOOST_ASSERT(pivot < static_cast<int>(size));

        if(array[pivot] < element)
            low = pivot + 1;
        else if(array[pivot] > element) 
            up = pivot - 1;
        else
            return static_cast<int>(pivot);
    }
    return -1;
}


template<class T, long unsigned Size>
inline void insert_array_safe(std::array<T,Size>& array, int size, int index, T element) {
    std::move_backward(
        std::move_iterator<typename std::array<T,Size>::iterator>(array.begin() + index),
        std::move_iterator<typename std::array<T,Size>::iterator>(array.begin() + size),
        std::move_iterator<typename std::array<T,Size>::iterator>(array.begin() + size + 1)
        );
    array[index] = std::move(element);
}

template<class T, long unsigned Size>
inline void  two_way_array_move(
        std::array<T,Size>& input, 
        std::array<T,Size>& smaller, 
        std::array<T,Size>& bigger, 
        unsigned pivot_index) 
{
    std::move(
        std::move_iterator<typename std::array<T,Size>::iterator>(input.begin()),
        std::move_iterator<typename std::array<T,Size>::iterator>(input.begin() + pivot_index + 1),
        std::move_iterator<typename std::array<T,Size>::iterator>(smaller.begin())
        );
    std::move(
        std::move_iterator<typename std::array<T,Size>::iterator>(input.begin() + pivot_index + 1),
        std::move_iterator<typename std::array<T,Size>::iterator>(input.end()),
        std::move_iterator<typename std::array<T,Size>::iterator>(bigger.begin())
        );
}

template<class T, long unsigned Size>
inline void  one_way_array_move(
        std::array<T,Size>& input, 
        std::array<T,Size>& bigger, 
        unsigned pivot_index) 
{
    std::move(
        std::move_iterator<typename std::array<T,Size>::iterator>(input.begin() + pivot_index + 1),
        std::move_iterator<typename std::array<T,Size>::iterator>(input.end()),
        std::move_iterator<typename std::array<T,Size>::iterator>(bigger.begin())
        );
}

template<class T>
void insert_into_array(T* array, unsigned size, T element, unsigned index) {
    if (size > index) {
        std::memmove(
            reinterpret_cast<void*>(index + 1 + array), // destination 
            reinterpret_cast<void*>(index + array),     // src
            sizeof(T)*(size - index)
            );
    }
    BOOST_ASSERT(index <= size);
    array[index] = element;
}

} // end namespace utils

#endif // INCLUDED_ARRAY_MANIPULATION
