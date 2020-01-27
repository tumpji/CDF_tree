#include <iostream>
#include <algorithm>
#include <array>

#include <map>
#include <random>

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE MyTest

#include <boost/test/unit_test.hpp>

#include "cdf_tree_main.h"

BOOST_AUTO_TEST_CASE( CDFTree_constructor ) {
    CDFTree<float> d;
}

BOOST_AUTO_TEST_CASE( CDFTree_test_all) {
    CDFTree<float> d;
    for(unsigned i = 0; i < 16; ++i) {
        d.insert_sample(i);
    }
    for(unsigned i = 0; i < 16; ++i) {
        BOOST_CHECK(d.search_count(i) == 1);
        BOOST_CHECK_CLOSE(d.search_PDF(i), 1./16, 0.001);
        BOOST_CHECK_CLOSE(d.search_CDF(i), (i+1)/16., 0.001);
    }
    float v = 0;
    for (float g = 0 + 1./32; g < 1.; g += 1./16) {
        BOOST_CHECK_CLOSE(d.inverse_search_CDF(g), v, 0.001);
        v += 1;
    }

    BOOST_CHECK_CLOSE(d.maximal_element(), 15, 0.001);
    BOOST_CHECK_CLOSE(d.minimal_element(), 0, 0.001);
    BOOST_CHECK_CLOSE(d.inverse_search_CDF(1.0), d.maximal_element(), 0.001);
    d.clear();
    d.insert_sample(-10.);
    BOOST_CHECK_CLOSE(d.search_PDF(-10.), 1., 0.001);
    BOOST_CHECK_CLOSE(d.inverse_search_CDF(1.0), d.maximal_element(), 0.001);

}

BOOST_AUTO_TEST_CASE( tree_constructor ) {
    auto root = RootNodeCluster<int>::factory();
}

BOOST_AUTO_TEST_CASE( root_insert_search_PDF ) {
    constexpr bool verbose = false;
    auto root = RootNodeCluster<int>::factory();
    unsigned suma_full = 0;

    for(int i = 0; i < 10; ++i) {
        unsigned long long a = (i+1)*(i+2)-1;
        unsigned long long t = root->insert_sample(i, a);
        suma_full += a;
        BOOST_CHECK(a == t);

        if (verbose) {
            std::cout << "-----------------------------" << std::endl;
            root->print();
        }

        //std::cout << "T: " << root->inverse_search_CDF(i+1) << "\ti: " << i << std::endl;
        BOOST_CHECK(root->inverse_search_CDF(suma_full) == i);
        BOOST_CHECK(root->inverse_search_CDF(0) == 0);
        root->sanity_check();
    }

    unsigned suma_partial = 0;
    for (int key = 0; key < 10; ++key) {
        unsigned long long true_sample_size = (key+1)*(key+2)-1;
        unsigned long long retrived_sample_size = root->search_PDF(key);
        BOOST_CHECK(true_sample_size == retrived_sample_size);

        for (unsigned i = suma_partial + 1; i <= suma_partial + true_sample_size; ++i) {
            int retrived_key_by_cdf = root->inverse_search_CDF(i);
            BOOST_CHECK(retrived_key_by_cdf == key);
        }
        suma_partial += true_sample_size;

        BOOST_CHECK(true_sample_size == retrived_sample_size);

        unsigned retrived_cdf = root->search_CDF(key);
        BOOST_CHECK(retrived_cdf == suma_partial);

        if (verbose)
            std::cout 
                << "key: " << key 
                << "\tsamples: " << true_sample_size 
                << "\tretrived_sample_size: " << retrived_sample_size
                << "\ttrue_cdf: " << suma_partial
                << "\tretrived_cdf: " << retrived_cdf
                << std::endl;
    }
    root->sanity_check();
}

BOOST_AUTO_TEST_CASE( insert_sequential ) {
    constexpr bool verbose = false;
    int size = 10000;
    auto root = RootNodeCluster<int>::factory();

    for(int i = -size; i < size; ++i) {
        unsigned long long a = (i+1)*(i+2)+1;
        unsigned long long t = root->insert_sample(i, a);
        BOOST_CHECK(a == t);
        if (verbose) {
            std::cout << "-----------------------------" << std::endl;
            root->print();
        }
        //root->sanity_check();
    }
    for (int i = -size; i < size; ++i) {
        unsigned long long a = (i+1)*(i+2)+1;
        unsigned long long g = root->search_PDF(i);
        if (verbose)
            std::cout << "key: " << i << "\trecieved: " << g << "\texpected: " << a << std::endl;
        BOOST_CHECK(g == a);
    }
    root->sanity_check();
}

BOOST_AUTO_TEST_CASE( insert_mixed ) {
    constexpr unsigned size = 200000;
    constexpr unsigned ns = 50;

    for(unsigned i = 0; i < ns; ++i) 
    {
        auto r1 = RootNodeCluster<int>::factory();
        std::map<int, unsigned long long> reference;

        std::random_device dev;
        std::mt19937 rng(dev());
        //std::mt19937 rng(42);
        std::uniform_int_distribution<std::mt19937::result_type> dist_key(-1000, 1000);
        std::uniform_int_distribution<std::mt19937::result_type> dist_value(1, 20);

        for(unsigned index = 0; index < size; ++index) {
            int key = dist_key(rng);
            unsigned long long samples = dist_value(rng);

            auto it = reference.find(key);
            if (it == reference.end()) {
                reference[key] = samples;
            } else {
                it->second += samples;
            }

            //std::cout << "Inserting sample to network: " << key << " " << samples << std::endl;
            //r2->print();
            r1->insert_sample(key, samples);
        }

        for (auto p : reference) {
            BOOST_CHECK(r1->search_PDF(p.first) == p.second);
        }
        r1->sanity_check(); 
    }
}

//BOOST_AUTO_TEST_CASE(CDF
BOOST_AUTO_TEST_CASE( lower_or_equal_bound ) {
    int array[6] = {-5,4,20,30,70,130};

    int query[13] = {-10, -5, -1, 4, 8, 20, 21, 30, 69, 70, 100, 130, 150};
    unsigned resul[13] = {0,    1,  1, 2, 2,  3,  3,  4,  4,  5,   5,   6,   6};

    for (unsigned i = 0; i < 13; ++i) 
        BOOST_CHECK( ::utils::lower_or_equal_bound(array, 6, query[i]) == resul[i] );
    // end with equal
    for (unsigned i = 0; i < 6; ++i) {
        BOOST_CHECK(::utils::lower_or_equal_bound(array, i, array[i]) == i);
        BOOST_CHECK(::utils::lower_or_equal_bound(array, i, array[i]+1) == i);
        BOOST_CHECK(::utils::lower_or_equal_bound(array, i, array[i]-1) == i);
    }


    int array2[7] = {-5,4,20,30,70,130, 300};

    int query2[15] = {-10, -5, -1, 4, 8, 20, 21, 30, 69, 70, 100, 130, 150, 300, 305};
    unsigned resul2[15] = {0,    1,  1, 2, 2,  3,  3,  4,  4,  5,   5,   6,   6,   7,   7};

    for (unsigned i = 0; i < 15; ++i) 
        BOOST_CHECK( ::utils::lower_or_equal_bound(array2, 7, query2[i]) == resul2[i] );
}

BOOST_AUTO_TEST_CASE( binary_search ) {
    int array[6] = {-5,4,20,30,70,130};

    int query[13] = {-10, -5, -1, 4, 8, 20, 21, 30, 69, 70, 100, 130, 150};
    int resul[13] = { -1,  0, -1, 1,-1,  2, -1,  3, -1,  4,  -1,   5,  -1};

    for (unsigned i = 0; i < 13; ++i) {
        BOOST_CHECK( ::utils::binary_search(array, 6, query[i]) == resul[i]);
    }


    int array2[7] = {-5,4,20,30,70,130, 300};

    int query2[15] = {-10, -5, -1, 4, 8, 20, 21, 30, 69, 70, 100, 130, 150, 300, 333};
    int resul2[15] = { -1,  0, -1, 1,-1,  2, -1,  3, -1,  4,  -1,   5,  -1,   6,   -1};

    for (unsigned i = 0; i < 15; ++i) 
        BOOST_CHECK( ::utils::binary_search(array2, 7, query2[i]) == resul2[i] );
}

/*
BOOST_AUTO_TEST_CASE( search_in_array ) {
    std::array<int, 13> odd_array ({1,2,3,4,5,6,7,8,9,10,11,12,13});

    for (int a : odd_array) {
        unsigned res = std::lower_bound(odd_array.begin(), odd_array.end(), a) - odd_array.begin();
        unsigned tes = ::utils::search_in_array(odd_array.data(), 13, a);
        BOOST_CHECK(res == tes);
        BOOST_CHECK(odd_array[tes] >= a);
    }

    std::array<int, 12> even_array ({1,2,3,4,5,6,7,8,9,10,11,120});
    for (int a : even_array) {
        unsigned res = std::lower_bound(even_array.begin(), even_array.end(), a) - even_array.begin();
        unsigned tes = ::utils::search_in_array(even_array.data(), 12, a);
        BOOST_CHECK(res == tes);
        BOOST_CHECK(even_array[tes] >= a);
    }
    unsigned a = ::utils::search_in_array(even_array.data(), 12, 0);
    BOOST_CHECK(a == 0);
    a = ::utils::search_in_array(even_array.data(), 12, 1231);
    BOOST_CHECK(a == 12);
}
*/

BOOST_AUTO_TEST_CASE( insert_into_array ) {
    int array[20] = {1,2,3,4,5,6,7};

    ::utils::insert_into_array(array, 7u, 0, 0u);
    int array_test1[20] = {0,1,2,3,4,5,6,7};

    for (int i = 0; i < 20; ++i)
        BOOST_CHECK(array[i] == array_test1[i]);

    ::utils::insert_into_array(array, 8u, 10, 1u);
    int array_test2[20] = {0,10,1,2,3,4,5,6,7};

    for (int i = 0; i < 20; ++i)
        BOOST_CHECK(array[i] == array_test2[i]);

    ::utils::insert_into_array(array, 19u, 11, 19u);
    int array_test3[20] = {0,10,1,2,3,4,5,6,7,0,0,0,0,0,0,0,0,0,0,11};

    for (int i = 0; i < 20; ++i)
        BOOST_CHECK(array[i] == array_test3[i]);
}


using REFERENCE_TYPE = std::map<int, double>;
using TESTED_TYPE = std::map<int, double>;

std::pair<std::map<int, double>, TESTED_TYPE> generate_random_reference_pair() {
    std::map<int, double> reference;
    TESTED_TYPE tested_type; 

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(-1000, 1000);

    for(int i = 0; i < 2000*2; ++i) {
        int new_sample = dist(rng);

        reference[new_sample] += 1;
        tested_type[new_sample] += 1;
    }

    return std::make_pair(reference, tested_type);
}


BOOST_AUTO_TEST_CASE( random_insert ) 
{
    std::pair<REFERENCE_TYPE, TESTED_TYPE> data = generate_random_reference_pair();
    REFERENCE_TYPE& reference = data.first;
    TESTED_TYPE& tested_type = data.second;

    for(std::pair<int, double> x : reference) {
        BOOST_CHECK(x.second == tested_type[x.first]);
    }
} 

