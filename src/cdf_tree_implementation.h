#if not defined INCLUDED_CDF_TREE_IMPLEMENTATION
#define INCLUDED_CDF_TREE_IMPLEMENTATION

#include <memory>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <array>
#include <algorithm>
#include <cstring>

#include <boost/assert.hpp>

#include "array_manip.h"
#include "cdf_tree_main.h"

///////////////////////////////////////////////////
///////////////// CONSTRUCTORS //////////////////// 


template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::RootNodeCluster() {
    size = 0;
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::InternalNodeCluster() : RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>() {
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::ExternalNodeCluster() {
    size = 0;
}


///////////////////////////////////////////////////
//////////////// search_PDF ///////////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_PDF(Type e) const {
    unsigned index = utils::lower_or_equal_bound(data, size, e);

    if (size == 0) {
        if (children[0] == nullptr) 
            return 0;
        return children[0]->search_PDF(e);
    }
    else 
        return children[index]->search_PDF(e);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_PDF(Type e) const {
    BOOST_ASSERT(size > 0);
    unsigned index = utils::lower_or_equal_bound(data, size, e);
    return children[index]->search_PDF(e);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_PDF(Type e) const {
    BOOST_ASSERT(size > 0);
    int index = utils::binary_search(data, size, e);
    if (index == -1)
        return 0;
    else
        return frequencies[index];
}


///////////////////////////////////////////////////
///////////////// search_CDF //////////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
CumFreqType RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_CDF(Type e) const {
    unsigned index = utils::lower_or_equal_bound(data, size, e);
    if (size == 0) {
        if (children[0] == nullptr)
            return 0;
        return children[0]->search_CDF(e);
    }
    else {
        CumFreqType sum = 0;
        for (unsigned i = 0; i < index; ++i)
            sum += cached_sums[i];
        return sum + children[index]->search_CDF(e);
    }
}
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
CumFreqType InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_CDF(Type e) const {
    BOOST_ASSERT(size > 0);
    unsigned index = utils::lower_or_equal_bound(data, size, e);

    CumFreqType sum = 0;
    for (unsigned i = 0; i < index; ++i)
        sum += cached_sums[i];
    return sum + children[index]->search_CDF(e);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
CumFreqType ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::search_CDF(Type e) const {
    BOOST_ASSERT(size > 0);
    unsigned index = utils::lower_or_equal_bound(data, size, e);
    CumFreqType sum = 0;
    for (unsigned i = 0; i < index; ++i)
        sum += frequencies[i];
    return sum;
}

///////////////////////////////////////////////////
///////////////// draw_CDF ////////////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::inverse_search_CDF(CumFreqType sum) const {
    if (children[0] == nullptr)
        throw std::runtime_error("Inverse_search_CDF on empty tree");

    unsigned index = 0;
    for (; index < size + 1; ++index) 
        if (sum > cached_sums[index]) 
            sum -= cached_sums[index];
        else
            return children[index]->inverse_search_CDF(sum);
    throw std::runtime_error("Inverse search failed");
    return data[size-1];
}
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::inverse_search_CDF(CumFreqType sum) const {
    BOOST_ASSERT(size > 0);
    unsigned index = 0;
    for (; index < size + 1; ++index)
        if (sum > cached_sums[index])
            sum -= cached_sums[index];
        else 
            return children[index]->inverse_search_CDF(sum);
    throw std::runtime_error("Inverse search failed");
    return data[size-1];
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::inverse_search_CDF(CumFreqType sum) const {
    BOOST_ASSERT(size > 0);
    unsigned index = 0;
    for (; index < size; ++index)
        if (sum <= frequencies[index])
            return data[index];
        else
            sum -= frequencies[index];
    throw std::runtime_error("Inverse search failed");
    return data[size-1];

        /*
        if (sum > frequencies[index])
            sum -= frequencies[index];
        else
            break;
        */
    //BOOST_ASSERT(sum <= frequencies[index]);
    //return data[index];
}



///////////////////////////////////////////////////
////////////// insert_sample //////////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::insert_sample(Type e, FreqType number) {
    if (number == 0)
        return number;

    if (size == 0) {
        if (children[0] == nullptr) {
            children[0] = std::make_shared<ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>>();
            cached_sums[0] = 0; 
            children[0]->parent = thisptr;
        }

        cached_sums[0] += number;
        return children[0]->insert_sample(e, number);
    } else {
        unsigned index = utils::lower_or_equal_bound(data, size, e);
        cached_sums[index] += number;
        return children[index]->insert_sample(e, number);
    }
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::insert_sample(Type e, FreqType number) {
    BOOST_ASSERT(size > 0);
    unsigned index = utils::lower_or_equal_bound(data, size, e);
    cached_sums[index] += number;
    return children[index]->insert_sample(e, number);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
FreqType ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::insert_sample(Type e, FreqType number) {
    unsigned index = utils::lower_bound(data, size, e);

    // leaf has this element already saved
    if (index < size and data[index] == e) { 
        frequencies[index] += number;
        return frequencies[index];
    } 

    utils::insert_into_array(data, size, e, index);
    utils::insert_into_array(frequencies, size, number, index);
    size += 1;

    if (size == MaxSize) // does not need to split
        split();
    return number;
}


///////////////////////////////////////////////////
/////////////////// min  + max elements ///////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::minimal_element() const {
    if (children[0]) 
        return children[0]->minimal_element();
    else
        throw std::runtime_error("Mininal Element on empty tree");
}
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::maximal_element() const {
    if (children[0]) 
        return children[size]->maximal_element();
    else
        throw std::runtime_error("Mininal Element on empty tree");
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::minimal_element() const {
    return children[0]->minimal_element();
}
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::maximal_element() const {
    return children[size]->maximal_element();
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::minimal_element() const {
    return data[0];
}
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
Type ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::maximal_element() const {
    return data[size-1];
}
///////////////////////////////////////////////////
///////////////// FACTORY       ///////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::shared_ptr<RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory() {
    auto x = std::make_shared<RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>>();
    x->thisptr = x;
    return x;
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::shared_ptr<InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory() {
    auto x = std::make_shared<InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>>();
    x->thisptr = x;
    return x;
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::shared_ptr<ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory() {
    auto x = std::make_shared<ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>>();
    x->thisptr = x;
    return x;
}

///////////////////////////////////////////////////
///////////////// Register split ///////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::register_split(
        std::shared_ptr<NodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> other, 
        Type pivot,
        CumFreqType sum) 
{
    unsigned index = utils::lower_or_equal_bound(data, size, pivot);
    BOOST_ASSERT(size < MaxSize);
    BOOST_ASSERT(index <= size);
    BOOST_ASSERT(sum > 0);

    utils::insert_into_array(data, size, pivot, index); 
    utils::insert_into_array(cached_sums, size+1, sum, index+1); 
    utils::insert_array_safe(children, size+1, index+1, other);
    size += 1;
    other->parent = thisptr;
    cached_sums[index] -= sum;

    if (size == MaxSize)
        split();
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::split() {
    constexpr unsigned size_small = (MaxSize-1)/2;
    constexpr unsigned pivot_index = (MaxSize-1)/2;
    constexpr unsigned size_big = (MaxSize)/2;

    // create greater element
    auto small_ptr = InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory();
    auto big_ptr   = InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory();

    // data
    std::memmove(small_ptr->data, &data[0],             sizeof(data[0])*size_small);
    std::memmove(big_ptr->data,   &data[pivot_index+1], sizeof(data[0])*size_big);

    // cached_sums
    CumFreqType s1 = 0, s2 = 0;
    for (unsigned i = 0; i < size_small + 1; ++i)
        s1 += cached_sums[i];
    for (unsigned i = size_small + 1; i < MaxSize; ++i)
        s2 += cached_sums[i];
    std::memmove(small_ptr->cached_sums, &cached_sums[0],            sizeof(cached_sums[0])*(size_small+1));
    std::memmove(big_ptr->cached_sums,   &cached_sums[size_small+1], sizeof(cached_sums[0])*(size_big+1));

    // ptrs
    utils::two_way_array_move(children, small_ptr->children, big_ptr->children, pivot_index);

    // size
    small_ptr->size = size_small;
    big_ptr->size = size_big;

    // parent
    small_ptr->parent = thisptr;
    big_ptr->parent = thisptr;

    // change children
    for (unsigned i = 0; i < big_ptr->size + 1; ++i) {
        BOOST_ASSERT(big_ptr->children[i] != nullptr);
        big_ptr->children[i]->parent = big_ptr;
    }
    for (unsigned i = 0; i < small_ptr->size + 1; ++i) {
        BOOST_ASSERT(big_ptr->children[i] != nullptr);
        small_ptr->children[i]->parent = small_ptr;
    }

    // roots pivot
    data[0] = data[pivot_index];
    size = 1;
    // roots ptrs
    children[0] = small_ptr;
    children[1] = big_ptr;

    cached_sums[0] = s1;
    cached_sums[1] = s2;
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::split() {
    constexpr unsigned size_less = (MaxSize)/2; 
    constexpr unsigned pivot_index = (MaxSize)/2;
    constexpr unsigned size_big = (MaxSize-1)/2;

    Type new_pivot = data[pivot_index];

    // create greater element
    std::shared_ptr<InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> big_ptr = InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory();

    // data
    std::memmove(big_ptr->data,     &data[pivot_index+1],     sizeof(data[0])*size_big);

    // cached sums
    CumFreqType sum = 0;
    for (unsigned i = pivot_index+1; i < MaxSize+1; ++i)
        sum += cached_sums[i];
    std::memmove(big_ptr->cached_sums, &cached_sums[pivot_index+1], sizeof(cached_sums[0])*(size_big+1));
    
    // ptrs
    utils::one_way_array_move(children, big_ptr->children, pivot_index);
    // size
    big_ptr->size = size_big;
    // parent
    big_ptr->parent = parent;
    // change children
    for (unsigned i = 0; i < big_ptr->size + 1; ++i)
        big_ptr->children[i]->parent = big_ptr;
    // node size
    size = size_less;
    // his pivot
    parent.lock()->register_split(big_ptr, new_pivot, sum);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::split() {
    constexpr unsigned half = (MaxSize) / 2;
    std::shared_ptr<ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>> greater_ptr = ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::factory();

    std::memmove(greater_ptr->data, &data[half], sizeof(data[0])*(MaxSize - half));
    //std::memset(&data[half], 0, sizeof(data[0])*(MaxSize - half));
    CumFreqType sum = 0;
    //for (unsigned i = half; i < MaxSize - half; ++i)
    for (unsigned i = half; i < MaxSize; ++i)
        sum += frequencies[i];
    std::memmove(greater_ptr->frequencies, &frequencies[half], sizeof(frequencies[0])*(MaxSize - half));
    //std::memset(&frequencies[half], 0, sizeof(frequencies[0])*(MaxSize - half));

    greater_ptr->size = MaxSize - half;
    size = half;

    parent.lock()->register_split(greater_ptr, greater_ptr->data[0], sum);
}


///////////////////////////////////////////////////
////////////////// SanityChecks ///////////////////

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::sanity_check() const {
    // size
    BOOST_ASSERT(size >= 0);
    BOOST_ASSERT(size < MaxSize);

    // order of elements
    for (unsigned i = 1; i < size; ++i)
        BOOST_ASSERT(data[i-1] < data[i]);

    // filled children
    if (size != 0 or children[0] != nullptr)
        for (unsigned i = 0; i < size + 1; ++i) {
            BOOST_ASSERT(children[i] != nullptr);
            BOOST_ASSERT(children[i]->parent.lock() == thisptr.lock());
            children[i]->sanity_check();
        }

    for (unsigned i = 0; i < size; ++i) {
        BOOST_ASSERT(children[i]->maximal_element() < data[i]);
        BOOST_ASSERT(children[i+1]->minimal_element() >= data[i]);
    }
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::sanity_check() const {
    BOOST_ASSERT(size >= 1);
    BOOST_ASSERT(size < MaxSize);

    // order of elements
    for (unsigned i = 1; i < size; ++i)
        BOOST_ASSERT(data[i-1] < data[i]);

    // filled children & backptrs
    for (unsigned i = 0; i < size + 1; ++i) {
        BOOST_ASSERT(children[i] != nullptr);
        BOOST_ASSERT(children[i]->parent.lock() == thisptr.lock());
        children[i]->sanity_check();
    }

    // sum of node cache are equal to parent cache
    RootNodeClusterType* parent = dynamic_cast<RootNodeClusterType*>(this->parent.lock().get());
    unsigned index_of_parent = utils::lower_or_equal_bound(parent->data, parent->size, data[0]);
    unsigned long long sum = 0;
    for(unsigned i = 0; i < size + 1; ++i) 
        sum += cached_sums[i];
    BOOST_ASSERT(parent->cached_sums[index_of_parent] == sum);
}


template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::sanity_check() const {
    // size
    BOOST_ASSERT(size >= 1);
    BOOST_ASSERT(size < MaxSize);

    // order
    for (unsigned i = 1; i < size; ++i)
        BOOST_ASSERT(data[i-1] < data[i]);

    RootNodeClusterType* parent = dynamic_cast<RootNodeClusterType*>(this->parent.lock().get());
    unsigned index_of_parent = utils::lower_or_equal_bound(parent->data, parent->size, data[0]);
    unsigned long long sum = 0;
    for(unsigned i = 0; i < size; ++i) 
        sum += frequencies[i];
    BOOST_ASSERT(parent->cached_sums[index_of_parent] == sum);

}

#endif // INCLUDED_CDF_TREE_IMPLEMENTATION
