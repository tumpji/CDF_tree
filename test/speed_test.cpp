#include <memory>
#include <random>
#include <limits>
#include <iostream>

#include "cdf_tree_main.h"

int main () {

    ExternalNodeCluster<int> ena;
    InternalNodeCluster<int> ina;
    std::shared_ptr<RootNodeCluster<int>> rna = RootNodeCluster<int>::factory();

    std::cout << "RootCluster: " 
        << sizeof(RootNodeCluster<int,4096,unsigned,unsigned long long, false>) 
        << "\t" << "Elements: " << rna->MaxSize  << std::endl;
    std::cout << "InteCluster: " 
        << sizeof(InternalNodeCluster<int,4096,unsigned,unsigned long long, false>) 
        << "\t" << "Elements: " << rna->MaxSize  << std::endl;
    std::cout << "ExteCluster: " 
        << sizeof(ExternalNodeCluster<int,4096,unsigned,unsigned long long, false>)
        << "\t" << "Elements: " << ena.MaxSize  << std::endl;

    std::random_device rd;
    std::mt19937 key_generator(rd());
    std::uniform_int_distribution<int> distribution(
            std::numeric_limits<int>::min(),
            std::numeric_limits<int>::max()
            );

    constexpr unsigned size = 100000000; // 10**8
    for (unsigned i = 0; i < size; ++i) {
        //int key = distribution(key_generator);
        int key = key_generator();
        rna->insert_sample(key);
    }
}
