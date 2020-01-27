#if not defined INCLUDED_CDF_TREE_MAIN
#define INCLUDED_CDF_TREE_MAIN

#include <ostream>
#include <memory>
#include <array>
#include <cmath>
#include <boost/assert.hpp>

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
class RootNodeCluster;

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
class InternalNodeCluster;

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
class ExternalNodeCluster;


template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (
        std::ostream&, 
        const RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (
        std::ostream&, 
        const InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);

template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
std::ostream& operator<< (
        std::ostream&, 
        const ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);


template<class Type, unsigned PageSize, class FreqType, class CumFreqType, bool overflow_check>
class NodeCluster {
    using NodeClusterType = NodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using NodeClusterPtrType = std::shared_ptr<NodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>>;

    virtual ~NodeCluster() {}

    // element -> probability
    virtual FreqType search_PDF(Type) const = 0;
    // update
    virtual FreqType insert_sample(Type, FreqType number=1) = 0;
    // element -> cummulative probability
    virtual CumFreqType search_CDF(Type e) const = 0;
    // cummulative proabibility -> element
    virtual Type inverse_search_CDF(CumFreqType cdf) const = 0;

    virtual Type minimal_element() const = 0;
    virtual Type maximal_element() const = 0;
    
    virtual void print(unsigned) const = 0;
    virtual void sanity_check () const = 0;

    virtual void register_split(NodeClusterPtrType, Type pivot, CumFreqType sum) = 0;

    std::weak_ptr<NodeClusterType> parent;
    std::weak_ptr<NodeClusterType> thisptr;
    unsigned     size;

    friend class RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>;
    friend class InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>;
    friend class ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>;

    friend std::ostream& operator<<<>(std::ostream&, const RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);
    friend std::ostream& operator<<<>(std::ostream&, const InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);
    friend std::ostream& operator<<<>(std::ostream&, const ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);

    //constexpr static unsigned MaxSize = PageSize - 2*sizeof(parent);
};
 

template<
    class Type, 
    unsigned PageSize = 4096, 
    class FreqType = unsigned, 
    class CumFreqType = unsigned long long,
    bool overflow_check = true
    >
class RootNodeCluster: public NodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check> 
{
public:
    virtual ~RootNodeCluster() {};
    RootNodeCluster();

    using NodeClusterType = NodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using NodeClusterPtrType = std::shared_ptr<NodeClusterType>;
    using RootNodeClusterType = RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using RootNodeClusterPtrType = std::shared_ptr<RootNodeClusterType>;
    using InternalNodeClusterType = InternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using InternalNodeClusterPtrType = std::shared_ptr<InternalNodeClusterType>;
    using ExternalNodeClusterType = ExternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using ExternalNodeClusterPtrType = std::shared_ptr<ExternalNodeClusterType>;

    using NodeClusterType::size;
    using NodeClusterType::thisptr;

    static constexpr unsigned MaxSize = 
        (PageSize - sizeof(NodeClusterType) - sizeof(NodeClusterPtrType) - sizeof(CumFreqType)) / 
        (sizeof(CumFreqType) + sizeof(NodeClusterPtrType) + sizeof(Type));

    //static_assert(sizeof(RootNodeClusterType) < PageSize, "Page size overflow");

    virtual FreqType        insert_sample(Type, FreqType number=1) override;
    virtual FreqType        search_PDF(Type) const override;
    virtual CumFreqType     search_CDF(Type) const override;
    virtual Type            inverse_search_CDF(CumFreqType) const override;

    static RootNodeClusterPtrType factory();

    virtual Type minimal_element() const override;
    virtual Type maximal_element() const override;

    virtual void print(unsigned x = 0) const override;
    virtual void sanity_check () const override;
protected:
    virtual void register_split(NodeClusterPtrType, Type, CumFreqType) override;
    void split();

    Type         data [MaxSize];
    CumFreqType  cached_sums[MaxSize+1];
    std::array<NodeClusterPtrType, MaxSize+1> children;

    friend std::ostream& operator<<<>(std::ostream&, const RootNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);
    friend InternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    friend ExternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
};




template<
    class Type, 
    unsigned PageSize = 4096, 
    class FreqType = unsigned, 
    class CumFreqType = unsigned long long,
    bool overflow_check = true
    >
class InternalNodeCluster: public RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check> 
{
public:
    InternalNodeCluster();
    virtual ~InternalNodeCluster() {}

protected:
    using NodeClusterType = NodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using NodeClusterPtrType = std::shared_ptr<NodeClusterType>;
    using RootNodeClusterType = RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using RootNodeClusterPtrType = std::shared_ptr<RootNodeClusterType>;
    using InternalNodeClusterType = InternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using InternalNodeClusterPtrType = std::shared_ptr<InternalNodeClusterType>;

    using RootNodeClusterType::MaxSize;
    using RootNodeClusterType::children;
    using RootNodeClusterType::data;
    using RootNodeClusterType::cached_sums;
    using NodeClusterType::parent;
    using NodeClusterType::thisptr;
    using NodeClusterType::size;

    //static_assert(sizeof(InternalNodeClusterType) < PageSize, "Page size overflow");

    virtual FreqType insert_sample(Type e, FreqType number=1) override;
    virtual FreqType search_PDF(Type e) const override;
    virtual CumFreqType search_CDF(Type e) const override;
    virtual Type inverse_search_CDF(CumFreqType cdf) const override;

    virtual Type minimal_element() const override; 
    virtual Type maximal_element() const override;

    virtual void print(unsigned) const override;
    virtual void sanity_check () const override;

    static InternalNodeClusterPtrType factory();
    void split();

    friend class RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    friend std::ostream& operator<<<>(std::ostream&, const InternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);
};





template<
    class Type, 
    unsigned PageSize = 4096, 
    class FreqType = unsigned, 
    class CumFreqType = unsigned long long,
    bool overflow_check = true
    >
class ExternalNodeCluster: public NodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check> 
{
public:
    ExternalNodeCluster();
    virtual ~ExternalNodeCluster() {};

    using NodeClusterType = NodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    static constexpr unsigned MaxSize = 
        (PageSize - sizeof(NodeClusterType)) / (sizeof(FreqType) + sizeof(Type));

protected:
    using NodeClusterPtrType = std::shared_ptr<NodeClusterType>;
    using RootNodeClusterType = RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using RootNodeClusterPtrType = std::shared_ptr<RootNodeClusterType>;
    using InternalNodeClusterType = InternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using InternalNodeClusterPtrType = std::shared_ptr<InternalNodeClusterType>;
    using ExternalNodeClusterType = ExternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    using ExternalNodeClusterPtrType = std::shared_ptr<ExternalNodeClusterType>;


    using NodeClusterType::parent;
    using NodeClusterType::size;


    //static_assert(sizeof(ExternalNodeClusterType) < PageSize, "Page size overflow");

    virtual FreqType        insert_sample(Type, FreqType) override;
    virtual FreqType        search_PDF(Type) const override;
    virtual CumFreqType     search_CDF(Type) const override;
    virtual Type            inverse_search_CDF(CumFreqType) const override;

    virtual Type minimal_element() const override;
    virtual Type maximal_element() const override;

    static ExternalNodeClusterPtrType factory();

    virtual void print(unsigned) const override;
    virtual void sanity_check() const override;

    void split();
    virtual void register_split(NodeClusterPtrType, Type, CumFreqType) override {}

    Type        data [MaxSize];
    FreqType    frequencies[MaxSize];

    friend class InternalNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    friend class RootNodeCluster<Type, PageSize, FreqType, CumFreqType, overflow_check>;
    friend std::ostream& operator<<<>(std::ostream&, const ExternalNodeCluster<Type,PageSize,FreqType,CumFreqType,overflow_check>&);
};


template<class Type>
class CDFTree {
public:
    CDFTree();

    // element -> probability
    double search_PDF(Type) const;
    unsigned search_count(Type) const;
    // update
    double insert_sample(Type);
    double insert_sample(Type, unsigned i);
    // element -> cummulative probability
    double search_CDF(Type e) const;

    // cummulative proabibility -> element
    Type inverse_search_CDF(double) const;

    Type minimal_element() const;
    Type maximal_element() const;

    void clear();

protected:
    std::shared_ptr<RootNodeCluster<Type>> root;
    unsigned long long counter;
};

template<class Type>
CDFTree<Type>::CDFTree(){
    clear();
}

template<class Type>
void CDFTree<Type>::clear() {
    counter = 0;
    root = std::make_shared<RootNodeCluster<Type>>();
    root->thisptr = root;
}

template<class Type>
inline unsigned CDFTree<Type>::search_count(Type e) const {
    return root->search_PDF(e);
}
template<class Type>
inline double CDFTree<Type>::search_PDF(Type e) const {
    unsigned s = search_count(e);
    return static_cast<double>(s) / counter;
}
template<class Type>
inline double CDFTree<Type>::insert_sample(Type e) {
    unsigned long long s = root->insert_sample(e);
    counter += 1;
    return static_cast<double>(s) / counter;
}
template<class Type>
inline double CDFTree<Type>::insert_sample(Type e, unsigned i) {
    unsigned long long s = root->insert_sample(e, i);
    counter += i;
    return static_cast<double>(s) / counter;
}
template<class Type>
inline double CDFTree<Type>::search_CDF(Type e) const {
    unsigned long long s = root->search_CDF(e);
    return static_cast<double>(s) / counter;
}

template<class Type>
inline Type CDFTree<Type>::inverse_search_CDF(double e) const {
    unsigned long long b = static_cast<unsigned long long>(std::ceil(e*counter));
    if (b <= 0) {
        throw std::runtime_error("Inversion of CDF=0 is impossible to obtain");
    }
    assert(b <= counter);
    return root->inverse_search_CDF(b);
}
template<class Type>
inline Type CDFTree<Type>::minimal_element() const {
    return root->minimal_element();
}
template<class Type>
inline Type CDFTree<Type>::maximal_element() const {
    return root->maximal_element();
}

#include "cdf_tree_implementation.h"
#include "cdf_tree_print_and_debug.h"


#endif // INCLUDED_CDF_TREE_MAIN
