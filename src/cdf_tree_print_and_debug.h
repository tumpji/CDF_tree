#if not defined INCLUDED_CDF_TREE_PRINT_AND_DEBUG
#define INCLUDED_CDF_TREE_PRINT_AND_DEBUG

#include <iostream>
#include <ostream>
#include <sstream>
////////////////////// OSTREAM /////////////////////

// ROOT NODE CLUSTER
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (std::ostream& s, const RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>& item) {
    s << "Root Node [" << &item << "] size: " << item.size << "\n";
    for (unsigned i = 0; i < item.size; ++i) {
        s << "["<< item.children[i].get() << "]" << item.data[i];
    }
    s << "[" << item.children[item.size].get() << "]\n";
    return s;
}

// INTERNAL NODE CLUSTER
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (std::ostream& s, const InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>& item) {
    s << "Internal Node [" << item.thisptr.lock() 
        << "] parent: [" << item.parent.lock() << "] size: " << item.size << "\n";
    for (unsigned i = 0; i < item.size; ++i) {
        s << "["<< item.children[i].get() << "]" << item.data[i];
    }
    s << "[" << item.children[item.size].get() << "]\n";
    return s;
}

// EXTERNAL NODE CLUSTER
template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (std::ostream& s, const ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>& item) {
    s << "External Node [" << &item << "] parent: [" << item.parent.lock()
       <<  "] size: " << item.size << "\n";
    for (unsigned i = 0; i < item.size; ++i) {
        s << item.data[i] << " ";
    }
    s << '\n';
    return s;
}


////////////////////// PRINT /////////////////////

std::string add_tab_after_new_newline(std::string string, unsigned spaces) {
    constexpr unsigned multiplier = 2;
    std::string out;
    for (unsigned i = 0; i < spaces*multiplier; ++i) {
        out.push_back(' ');
    }

    for (auto x : string) {
        if (x == '\0')
            out.push_back('\n');
        else
            out.push_back(x);

        if (x == '\n') {
            for (unsigned i = 0; i < spaces*multiplier; ++i) {
                out.push_back(' ');
            }
        }
    }
    for (unsigned i = 0; i < spaces*multiplier; ++i) 
        out.pop_back();
    return out;
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::print(unsigned g) const {
    std::cout << *this;

    if (children[0])
    for (unsigned i = 0; i < size + 1; ++i)
        children[i]->print(g+1);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::print(unsigned g) const {
    std::stringstream s;
    s << *this;
    std::cout << add_tab_after_new_newline(s.str(), g);

    for (unsigned i = 0; i < size + 1; ++i)
        children[i]->print(g+1);
}

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
void ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>::print(unsigned i) const {
    std::stringstream s;
    s << *this;
    std::cout << add_tab_after_new_newline(s.str(), i);
}

#endif // INCLUDED_CDF_TREE_PRINT_AND_DEBUG
