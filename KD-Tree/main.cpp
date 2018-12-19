
// MIT License
//
// Copyright (c) 2018 degski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <filesystem>
#include <forward_list>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <random>
#include <stack>
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

#include <SFML/System.hpp>
#include <splitmix.hpp>
#include <plf/plf_nanotimer.h>

using Point = sf::Vector2f;

template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const Point & p_ ) noexcept {
    out_ << '<' << p_.x << ' ' << p_.y << '>';
    return out_;
}

#define nl '\n'

using Int = std::int32_t;

#if 0

// Bentley & McIlroy, "fat" partitioning scheme.
template<typename ForwardIt>
void quicksort_bentley_mcilroy ( ForwardIt first, ForwardIt last ) noexcept {

    if ( first == last ) {
        return;
    }

    auto pivot = * std::next ( first, std::distance ( first, last ) / 2 );

    ForwardIt left  = std::partition ( first, last, [ pivot ] ( const auto & em ) { return em < pivot ; } );
    ForwardIt right = std::partition (  left, last, [ pivot ] ( const auto & em ) { return pivot >= em; } );

    quicksort ( first, left );
    quicksort ( right, last );
}

template<typename ForwardIt>
void quicksort ( ForwardIt first, ForwardIt last ) noexcept {

    if ( first == last ) {
        return;
    }

    auto pivot = * std::next ( first, std::distance ( first, last ) / 2 );

    ForwardIt median = std::partition ( first, last, [ pivot ] ( const auto & em ) { return em < pivot; } );

    quicksort ( first, median );
    quicksort ( std::next ( median ),  last );
}

namespace fbt {

template<typename Data>
class BinTree;

namespace detail {

#define NODEID_INVALID_VALUE ( 0 )

struct NodeID {

    Int value;

    static const NodeID invalid;

    constexpr explicit NodeID ( ) noexcept :
        value { NODEID_INVALID_VALUE } { }
    explicit NodeID ( Int && v_ ) noexcept :
        value { std::move ( v_ ) } { }
    explicit NodeID ( const Int & v_ ) noexcept :
        value { v_ } { }
    explicit NodeID ( const std::size_t v_ ) noexcept :
        value { static_cast< Int > ( v_ ) } { }

    [[ nodiscard ]] constexpr const Int operator ( ) ( ) const noexcept {
        return value;
    }

    [[ nodiscard ]] bool operator == ( const NodeID rhs_ ) const noexcept {
        return value == rhs_.value;
    }
    [[ nodiscard ]] bool operator != ( const NodeID rhs_ ) const noexcept {
        return value != rhs_.value;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const NodeID id_ ) noexcept {
        if ( NODEID_INVALID_VALUE == id_.value ) {
            out_ << L'*';
        }
        else {
            out_ << static_cast< std::uint64_t > ( id_.value );
        }
        return out_;
    }
};

const NodeID NodeID::invalid { };


template<typename Data>
struct Node {

    NodeID left, right;

    using type = NodeID;
    using data_type = Data;

    constexpr Node ( ) noexcept {
    }
    template<typename ... Args>
    Node ( Args && ... args_ ) noexcept :
        data { std::forward<Args> ( args_ ) ... } {
    }
    template<typename ... Args>
    Node ( NodeID && l_, NodeID && r_, Args && ... args_ ) noexcept :
        left { std::move ( l_ ) },
        right { std::move ( r_ ) },
        data { std::forward<Args> ( args_ ) ... } {
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Node node_ ) noexcept {
        out_ << L'<' << node_.left << L' ' << node_.right << L'>';
        return out_;
    }

    Data data;
};

} // namespace detail.


template<typename Data>
class BinTree {

    public:

    using NodeID = detail::NodeID;
    using Node = detail::Node<Data>;
    using Nodes = std::vector<Node>;

    BinTree ( ) :
        root_node { 1 },
        m_nodes { { } } {
    }

    template<typename ... Args>
    [[ maybe_unused ]] NodeID addNode ( Args && ... args_ ) noexcept {
        const NodeID id { m_nodes.size ( ) };
        m_nodes.emplace_back ( std::forward<Args> ( args_ ) ... );
        return id;
    }

    class iterator {

        template<typename Data>
        friend class BinTree;

        typename Nodes::pointer m_ptr, m_end;

        public:

        using difference_type = typename Nodes::difference_type;
        using value_type = typename Nodes::value_type;
        using reference = typename Nodes::reference;
        using pointer = typename Nodes::pointer;
        using const_reference = typename Nodes::const_reference;
        using const_pointer = typename Nodes::const_pointer;
        using iterator_category = std::forward_iterator_tag;

        iterator ( BinTree & tree_ ) noexcept :
            m_ptr { tree_.m_nodes.data ( ) + 1 },
            m_end { tree_.m_nodes.data ( ) + tree_.m_nodes.size ( ) } {
        }

        [[ nodiscard ]] const bool is_valid ( ) const noexcept {
            return m_end != m_ptr;
        }

        [[ maybe_unused ]] iterator & operator ++ ( ) noexcept {
            ++m_ptr;
            return *this;
        }

        [[ nodiscard ]] reference operator * ( ) const noexcept {
            return *m_ptr;
        }

        [[ nodiscard ]] pointer operator -> ( ) const noexcept {
            return m_ptr;
        }
    };

    class const_iterator {

        template<typename Data>
        friend class BinTree;

        typename Nodes::pointer m_ptr, m_end;

        public:

        using difference_type = typename Nodes::difference_type;
        using value_type = typename Nodes::value_type;
        using reference = typename Nodes::reference;
        using pointer = typename Nodes::pointer;
        using const_reference = typename Nodes::const_reference;
        using const_pointer = typename Nodes::const_pointer;
        using iterator_category = std::forward_iterator_tag;

        const_iterator ( const BinTree & tree_ ) noexcept :
            m_ptr { tree_.m_nodes.data ( ) + 1 },
            m_end { tree_.m_nodes.data ( ) + tree_.m_nodes.size ( ) } {
        }

        [[ nodiscard ]] const bool is_valid ( ) const noexcept {
            return m_end != m_ptr;
        }

        [[ maybe_unused ]] const_iterator & operator ++ ( ) noexcept {
            ++m_ptr;
            return *this;
        }

        [[ nodiscard ]] const_reference operator * ( ) const noexcept {
            return *m_ptr;
        }

        [[ nodiscard ]] const_pointer operator -> ( ) const noexcept {
            return m_ptr;
        }
    };


    [[ nodiscard ]] const bool isLeaf ( const NodeID node_ ) const noexcept {
        return NodeID::invalid == m_nodes [ node_.value ].left and NodeID::invalid == m_nodes [ node_.value ].right;
    }
    [[ nodiscard ]] const bool isInternal ( const NodeID node_ ) const noexcept {
        return not ( isLeaf ( node_ ) );
    }


    [[ nodiscard ]] Node & operator [ ] ( const NodeID node_ ) noexcept {
        return m_nodes [ node_.value ];
    }
    [[ nodiscard ]] const Node & operator [ ] ( const NodeID node_ ) const noexcept {
        return m_nodes [ node_.value ];
    }

    // The number of valid nodes. This is not the same as the size of
    // the arcs-vector, which allows for some additional admin elements,
    // use nodesSize ( ) instead.
    [[ nodiscard ]] const Int nodeNum ( ) const noexcept {
        return static_cast<Int> ( m_nodes.size ( ) ) - 1;
    }

    // The size of the nodes-vector (allows for some admin elements).
    [[ nodiscard ]] const std::size_t nodesSize ( ) const noexcept {
        return m_nodes.size ( );
    }

    // Data members.

    NodeID root_node;

    private:

    Nodes m_nodes;
};

}

#endif

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

// Implicit full binary tree of size N.
template<typename T, std::size_t N, typename Int = std::int32_t>
struct Imp2DTree {

    // https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030


    using base_type = decltype ( T { }.x );
    using value_type = T;
    using pointer = T * ;
    using reference = T & ;
    using const_pointer = T const *;
    using const_reference = T const &;

    using container = std::vector<T>;
    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:

    template<typename RandomIt>
    void construct_recursive ( const pointer node_, RandomIt first_, RandomIt last_, const bool x_dim_ ) noexcept {
        if ( first_ == last_ ) {
            return;
        }
        const RandomIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        if ( x_dim_ ) {
            std::nth_element ( first_, median, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        }
        else {
            std::nth_element ( first_, median, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        }
        *node_ = *median;
        construct_recursive ( left  ( node_ ),               first_, median, not ( x_dim_ ) );
        construct_recursive ( right ( node_ ), std::next ( median ),  last_, not ( x_dim_ ) );
    }

    template<typename RandomIt>
    void construct_recursive ( const pointer node_, RandomIt first_, RandomIt last_ ) noexcept {
        return construct_recursive ( node_, first_, last_, m_dim_start );
    }

    container m_data;
    bool m_dim_start;

    public:

    Imp2DTree ( ) noexcept {
    }

    Imp2DTree ( const Imp2DTree & ) noexcept = delete;
    Imp2DTree ( Imp2DTree && ) noexcept = delete;

    Imp2DTree ( std::initializer_list<T> il_ ) noexcept :
        m_data { N, Point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_dim_start { pick_dimension ( std::begin ( il_ ), std::end ( il_ ) ) } {
        assert ( il_.size ( ) <= N );
        container points;
        std::copy ( std::begin ( il_ ), std::end ( il_ ), std::begin ( points ) );
        construct_recursive ( m_data.data ( ), std::begin ( points ), std::begin ( points ) + il_.size ( ) );
    }

    template<typename ForwardIt>
    Imp2DTree ( ForwardIt first_, ForwardIt last_ ) noexcept :
        m_data { N, Point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_dim_start { pick_dimension ( first_, last_ ) } {
        assert ( ( last_ - first_ ) <= N );
        construct_recursive ( m_data.data ( ), first_, last_ );
    }

    template<typename ForwardIt>
    [[ nodiscard ]] bool pick_dimension ( ForwardIt first_, ForwardIt last_ ) const noexcept {
        const auto x = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        const auto y = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        return ( x.second->x - x.first->x ) > ( y.second->y - y.first->y );
    }

    [[ nodiscard ]] constexpr pointer left ( const pointer p_ ) const noexcept {
        return p_ + ( p_ - m_data.data ( ) ) + 1;
    }
    [[ nodiscard ]] constexpr pointer right ( const pointer p_ ) const noexcept {
        return p_ + ( p_ - m_data.data ( ) ) + 2;
    }
    [[ nodiscard ]] constexpr pointer parent ( const pointer p_ ) const noexcept {
        return m_data.data ( ) + ( p_ - m_data.data ( ) - 1 ) / 2;
    }
    [[ nodiscard ]] constexpr const_pointer left ( const const_pointer p_ ) const noexcept {
        return p_ + ( p_ - m_data.data ( ) ) + 1;
    }
    [[ nodiscard ]] constexpr const_pointer right ( const const_pointer p_ ) const noexcept {
        return p_ + ( p_ - m_data.data ( ) ) + 2;
    }
    [[ nodiscard ]] constexpr const_pointer parent ( const const_pointer p_ ) const noexcept {
        return m_data.data ( ) + ( p_ - m_data.data ( ) - 1 ) / 2;
    }
    [[ nodiscard ]] constexpr Int left ( const Int i_ ) const noexcept {
        return 2 * i_ + 1;
    }
    [[ nodiscard ]] constexpr Int right ( const Int i_ ) const noexcept {
        return 2 * i_ + 2;
    }
    [[ nodiscard ]] constexpr Int parent ( const Int i_ ) const noexcept {
        return ( i_ - 1 ) / 2;
    }

    [[ nodiscard ]] bool is_leaf ( pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) >= std::ptrdiff_t { N / 2 };
    }
    [[ nodiscard ]] bool is_leaf ( const_pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) >= std::ptrdiff_t { N / 2 };
    }
    [[ nodiscard ]] bool is_internal ( pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) < std::ptrdiff_t { N / 2 };
    }
    [[ nodiscard ]] bool is_internal ( const_pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) < std::ptrdiff_t { N / 2 };
    }
    [[ nodiscard ]] static constexpr bool is_leaf ( const Int i_ ) noexcept {
        assert ( N >= i_ );
        return i_ >= Int { N / 2 };
    }
    [[ nodiscard ]] static constexpr bool is_internal ( const Int i_ ) noexcept {
        assert ( N >= i_ );
        return i_ < Int { N / 2 };
    }
    [[ nodiscard ]] static constexpr bool is_leaf ( const std::size_t i_ ) noexcept {
        assert ( N >= i_ );
        return i_ >= std::size_t { N / 2 };
    }
    [[ nodiscard ]] static constexpr bool is_internal ( const std::size_t i_ ) noexcept {
        assert ( N >= i_ );
        return i_ < std::size_t { N / 2 };
    }

    [[ nodiscard ]] static base_type distance_squared ( const Point & p1_, const Point & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    private:

    struct Nearest {

        Point point;
        const_pointer found;
        base_type min_distance;

        Nearest ( ) { }

        Nearest ( const Point & point_ ) noexcept :
            point { point_ },
            found { nullptr },
            min_distance { std::numeric_limits<base_type>::max ( ) } {
        }

        void init ( ) noexcept {
            found = nullptr;
            min_distance = std::numeric_limits<base_type>::max ( );
        }
    };

    void nearest_recursive0 ( const const_pointer p_, Nearest & n_, bool dim_ ) const noexcept {
        const float d = Imp2DTree::distance_squared ( n_.point, *p_ );
        if ( d < n_.min_distance ) {
            n_.min_distance = d;
            n_.found = p_;
        }
        if ( is_leaf ( p_ ) ) {
            // std::cout << "leaf " << *p_ << nl;
            return;
        }
        // std::cout << "itnl " << *p_ << nl;
        const float dx = dim_ ? p_->x - n_.point.x : p_->y - n_.point.y;
        dim_ = not ( dim_ );
        nearest_recursive0 ( dx > base_type { 0 } ? left ( p_ ) : right ( p_ ), n_, dim_ );
        if ( ( dx * dx ) >= n_.min_distance ) {
            return;
        }
        nearest_recursive0 ( dx > base_type { 0 } ? right ( p_ ) : left ( p_ ), n_, dim_ );
    }

    void nearest_recursive0 ( const const_pointer p_, bool dim_ ) const noexcept {
        float d = Imp2DTree::distance_squared ( m_nearest.point, *p_ );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d;
            m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) ) {
            return;
        }
        d = dim_ ? p_->x - m_nearest.point.x : p_->y - m_nearest.point.y;
        dim_ = not ( dim_ );
        nearest_recursive0 ( d > base_type { 0 } ? left ( p_ ) : right ( p_ ), dim_ );
        if ( ( d * d ) < m_nearest.min_distance ) {
            nearest_recursive0 ( d > base_type { 0 } ? right ( p_ ) : left ( p_ ), dim_ );
        }
    }

    const_pointer nearest_recursive0 ( const Point & p_ ) const noexcept {
        m_nearest = { p_ };
        nearest_recursive0 ( m_data.data ( ), m_dim_start );
        return m_nearest.found;
    }

    mutable Nearest m_nearest;

    struct TaggedPtr {

        std::uintptr_t ptr : 63;
        std::uintptr_t dim : 1;

        TaggedPtr ( const_pointer p_, const bool d_ ) noexcept :
            ptr ( reinterpret_cast<std::uintptr_t> ( p_ ) ),
            dim ( d_ ) {
        }
    };

    void nearest_recursive ( const TaggedPtr & p_ ) const noexcept {
        const float d = Imp2DTree::distance_squared ( m_nearest.point, * reinterpret_cast<const_pointer> ( p_.ptr ) );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d;
            m_nearest.found = reinterpret_cast<const_pointer> ( p_.ptr );
        }
        if ( is_leaf ( reinterpret_cast<const_pointer> ( p_.ptr ) ) ) {
            // std::cout << "leaf " << *p_ << nl;
            return;
        }
        // std::cout << "itnl " << *p_ << nl;
        const float dx = p_.dim ? reinterpret_cast<const_pointer> ( p_.ptr )->x - m_nearest.point.x : reinterpret_cast<const_pointer> ( p_.ptr )->y - m_nearest.point.y;
        nearest_recursive ( TaggedPtr { dx > base_type { 0 } ? left  ( reinterpret_cast<const_pointer> ( p_.ptr ) ) : right ( reinterpret_cast<const_pointer> ( p_.ptr ) ), not ( p_.dim ) } );
        if ( ( dx * dx ) >= m_nearest.min_distance ) {
            return;
        }
        nearest_recursive ( TaggedPtr { dx > base_type { 0 } ? right ( reinterpret_cast<const_pointer> ( p_.ptr ) ) : left  ( reinterpret_cast<const_pointer> ( p_.ptr ) ), not ( p_.dim ) } );
    }

    const_pointer nearest_recursive ( const Point & p_ ) const noexcept {
        m_nearest = { p_ };
        nearest_recursive ( TaggedPtr ( m_data.data ( ), m_dim_start ) );
        return m_nearest.found;
    }

    public:

    [[ nodiscard ]] Point find_nearest ( const Point & point_ ) const noexcept {
        // Fastest and correct.
        return * nearest_recursive ( point_ );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Imp2DTree & tree_ ) noexcept {
        for ( const auto p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    [[ nodiscard ]] static constexpr Int bin_tree_size ( const Int i_ ) noexcept {
        Int s = 0, p = 2;
        while ( s < i_ ) {
            s = p - 1;
            p *= 2;
        }
        return s;
    }
};

[[ nodiscard ]] constexpr Int bin_tree_size ( const Int i_ ) noexcept {
    Int s = 0, p = 2;
    while ( s < i_ ) {
        s = p - 1;
        p *= 2;
    }
    return s;
}

/*

Int wmain67878 ( ) {

    splitmix64 rng;
    std::uniform_int_distribution<std::size_t> dis { 1u, 1000u };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    plf::nanotimer timer;
    double st;

    constexpr int n = 1'000;

    using Tree = Imp2DTree<Point, bin_tree_size ( n )>;

    std::vector<Point> points;

    for ( int i = 0; i < n; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    timer.start ( );

    Tree tree ( std::begin ( points ), std::end ( points ) );

    std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::size_t rv = 0u;

    timer.start ( );

    std::cout << "dfi   " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::cout << rv << nl;

    return EXIT_SUCCESS;

    Point point { 60.0f, 20.5f };

    const auto found = tree.find_nearest_recursive ( point );

    std::cout << nl << "nearest " << found << nl;

    return EXIT_SUCCESS;
}


*/

Int wmain ( ) {


    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    plf::nanotimer timer;
    double st;

    constexpr int n = 10'000;

    std::cout << bin_tree_size ( n ) << nl;

    std::vector<Point> points;

    for ( int i = 0; i < n; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    timer.start ( );

    Imp2DTree<Point, bin_tree_size ( n )> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    // std::cout << nl << tree << nl << nl;

    Point point { disx ( rng ), disy ( rng ) };

    auto found_impl = tree.find_nearest ( point );

    constexpr int cnt = 1'000'000;

    timer.start ( );
    for ( int i = 0; i < cnt; ++i ) {
        Point point { disx ( rng ), disy ( rng ) };
        found_impl += tree.find_nearest ( { disx ( rng ), disy ( rng ) } );
    }
    std::cout << "elapsed im " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;


    std::cout << "nearest im " << found_impl << nl;

    return EXIT_SUCCESS;
}



Int wmain89879789 ( ) {

    // std::vector<Point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<Point> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    Imp2DTree<Point, bin_tree_size ( 13 )> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << nl << tree << nl << nl;

    Point point { 3.1f, 2.9f };

    std::cout << nl << nl << "nearest " << nl << tree.find_nearest ( point ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << Imp2DTree<Point, bin_tree_size ( 13 )>::distance_squared ( p, point ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}

