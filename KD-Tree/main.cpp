
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

using point = sf::Vector2<float>;

template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const point & p_ ) noexcept {
    out_ << '<' << p_.x << ' ' << p_.y << '>';
    return out_;
}

#define nl '\n'


#if 0

// Bentley & McIlroy, "fat" partitioning scheme.
template<typename forward_it>
void quicksort_bentley_mcilroy ( forward_it first, forward_it last ) noexcept {

    if ( first == last ) {
        return;
    }

    auto pivot = * std::next ( first, std::distance ( first, last ) / 2 );

    forward_it left  = std::partition ( first, last, [ pivot ] ( const auto & em ) { return em < pivot ; } );
    forward_it right = std::partition (  left, last, [ pivot ] ( const auto & em ) { return pivot >= em; } );

    quicksort ( first, left );
    quicksort ( right, last );
}

template<typename forward_it>
void quicksort ( forward_it first, forward_it last ) noexcept {

    if ( first == last ) {
        return;
    }

    auto pivot = * std::next ( first, std::distance ( first, last ) / 2 );

    forward_it median = std::partition ( first, last, [ pivot ] ( const auto & em ) { return em < pivot; } );

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

// Implicit full binary tree.
template<typename T>
struct i2dtree {

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

    struct nearest_data {
        ::point point;
        const_pointer found;
        base_type min_distance;
    };

    template<typename forward_it>
    [[ nodiscard ]] bool pick_dimension ( forward_it first_, forward_it last_ ) const noexcept {
        const std::pair x = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        const std::pair y = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        return ( x.second->x - x.first->x ) > ( y.second->y - y.first->y );
    }

    [[ nodiscard ]] pointer left ( const pointer p_ ) const noexcept {
        return ( p_ + 1 ) + ( p_ - m_data.data ( ) );
    }
    [[ nodiscard ]] pointer right ( const pointer p_ ) const noexcept {
        return ( p_ + 2 ) + ( p_ - m_data.data ( ) );
    }
    [[ nodiscard ]] const_pointer left ( const const_pointer p_ ) const noexcept {
        return ( p_ + 1 ) + ( p_ - m_data.data ( ) );
    }
    [[ nodiscard ]] const_pointer right ( const const_pointer p_ ) const noexcept {
        return ( p_ + 2 ) + ( p_ - m_data.data ( ) );
    }

    [[ nodiscard ]] bool is_leaf ( const const_pointer p_ ) const noexcept {
        return m_num_nodes_div2 <= ( p_ - m_data.data ( ) );
    }

    template<typename random_it>
    void construct_recursive ( const pointer node_, random_it first_, random_it last_, bool dim_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        if ( dim_ ) {
            std::nth_element ( first_, median, last_, [ ] ( const point & a, const point & b ) { return a.x < b.x; } );
        }
        else {
            std::nth_element ( first_, median, last_, [ ] ( const point & a, const point & b ) { return a.y < b.y; } );
        }
        *node_ = *median;
        dim_ = not ( dim_ );
        if ( first_   != median ) {
            construct_recursive ( left  ( node_ ), first_, median, dim_ );
        }
        if ( ++median !=  last_ ) {
            construct_recursive ( right ( node_ ), median,  last_, dim_ );
        }
    }

    void nearest_recursive ( const const_pointer p_, bool dim_ ) const noexcept {
        assert ( p_ >= m_data.data ( ) );
        assert ( p_ < ( m_data.data ( ) + m_data.size ( ) ) );
        base_type d = i2dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) ) {
            // std::cout << *p_ << nl;
            return;
        }
        // std::cout << *p_ << nl;
        d = dim_ ? p_->x - m_nearest.point.x : p_->y - m_nearest.point.y;
        dim_ = not ( dim_ );
        nearest_recursive ( d > base_type { 0 } ? left ( p_ ) : right ( p_ ), dim_ );
        if ( ( d * d ) < m_nearest.min_distance ) {
            nearest_recursive ( d > base_type { 0 } ? right ( p_ ) : left ( p_ ), dim_ );
        }
    }

    container m_data;
    mutable nearest_data m_nearest;
    std::ptrdiff_t m_num_nodes_div2;
    bool m_dim_start;

    public:

    i2dtree ( const i2dtree & ) = delete;
    i2dtree ( i2dtree && ) noexcept = delete;

    i2dtree ( std::initializer_list<T> il_ ) noexcept :
        m_data { bin_tree_size<std::size_t> ( il_.size ( ) ), point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_num_nodes_div2 { static_cast<std::ptrdiff_t> ( m_data.size ( ) / 2 ) },
        m_dim_start { pick_dimension ( std::begin ( il_ ), std::end ( il_ ) ) } {
        if ( il_.size ( ) ) {
            container points;
            std::copy ( std::begin ( il_ ), std::end ( il_ ), std::begin ( points ) );
            construct_recursive ( m_data.data ( ), std::begin ( points ), std::end ( points ), m_dim_start );
        }
    }

    template<typename forward_it>
    i2dtree ( forward_it first_, forward_it last_ ) noexcept :
        m_data { bin_tree_size<std::size_t> ( static_cast<std::size_t> ( std::distance ( first_, last_ ) ) ), point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_num_nodes_div2 { static_cast<std::ptrdiff_t> ( m_data.size ( ) / 2 ) },
        m_dim_start { pick_dimension ( first_, last_ ) } {
        if ( first_ != last_  ) {
            construct_recursive ( m_data.data ( ), first_, last_, m_dim_start );
        }
    }

    i2dtree & operator = ( const i2dtree & ) = delete;
    i2dtree & operator = ( i2dtree && ) noexcept = delete;

    [[ nodiscard ]] const_pointer find_nearest_ptr ( const point & point_ ) const noexcept {
        m_nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        nearest_recursive ( m_data.data ( ), m_dim_start );
        return m_nearest.found;
    }

    [[ nodiscard ]] std::ptrdiff_t find_nearest_idx ( const point & point_ ) const noexcept {
        return find_nearest_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] point find_nearest_pnt ( const point & point_ ) const noexcept {
        return * find_nearest_ptr ( point_ );
    }

    [[ nodiscard ]] static constexpr base_type distance_squared ( const point & p1_, const point & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const i2dtree & tree_ ) noexcept {
        for ( const auto p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    void add_point ( const point & p_ ) {
        auto last = std::partition ( std::begin ( m_data ), std::end ( m_data ), [ ] ( const point & p ) { return p != point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) }; } );
        std::cout << "here\n";
        if ( likely ( std::end ( m_data ) != last ) ) {
            *last = p_;
            ++last;
        }
        else {
            m_data.push_back ( p_ );
            last = std::end ( m_data );
        }
        m_num_nodes_div2 = bin_tree_size ( last - std::begin ( m_data ) );
        container c { static_cast<std::size_t> ( m_num_nodes_div2 ), point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } };
        m_num_nodes_div2 /= 2;
        construct_recursive ( c.data ( ), std::begin ( m_data ), last, m_dim_start );
        std::swap ( c, m_data );
    }

    [[ nodiscard ]] static point find_nearest_linear_pnt ( const point & point_, const std::vector<point> & points_ ) noexcept {
        // Fastest up till 50 points.
        const_pointer found = nullptr;
        float min_distance = std::numeric_limits<float>::max ( );
        for ( const auto & v : points_ ) {
            const float d = distance_squared ( point_, v );
            if ( d < min_distance ) {
                found = &v;
                min_distance = d;
            }
        }
        return * found;
    }

    // private:

    template<typename T>
    [[ nodiscard ]] static constexpr T bin_tree_size ( const T i_ ) noexcept {
        T s = 0, p = 2;
        while ( true ) {
            s = p - 1;
            if ( s >= i_ ) {
                break;
            }
            p *= 2;
        }
        return s;
    }
};

/*

int wmain67878 ( ) {

    splitmix64 rng;
    std::uniform_int_distribution<std::size_t> dis { 1u, 1000u };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    plf::nanotimer timer;
    double st;

    constexpr int n = 1'000;

    using Tree = i2dtree<point>;

    std::vector<point> points;

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

    point point { 60.0f, 20.5f };

    const auto found = tree.find_nearest_recursive ( point );

    std::cout << nl << "nearest " << found << nl;

    return EXIT_SUCCESS;
}


*/

int wmain ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    for ( int i = 0; i < 10; ++i ) {

        plf::nanotimer timer;
        double st;

        constexpr int n = 100'000;

        std::vector<point> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        // std::cout << nl << tree << nl << nl;

        //point p { disx ( rng ), disy ( rng ) };

        //std::cout << "ptf " << p << nl;

        bool result = true;

        constexpr int cnt = 100'000;

        timer.start ( );
        for ( int i = 0; i < cnt; ++i ) {
            const point p { disx ( rng ), disy ( rng ) };
            bool r = tree.find_nearest_pnt ( p ) == i2dtree<point>::find_nearest_linear_pnt ( p, points );
            if ( not ( r ) ) {
                std::cout << p << tree.find_nearest_pnt ( p ) << i2dtree<point>::find_nearest_linear_pnt ( p, points ) << nl;
                exit ( 0 );
            }
            result &= r;
        }
        // std::cout << "elapsed im " << ( std::uint64_t ) ( timer.get_elapsed_ns ( ) / cnt ) << " ns" << nl;

        std::cout << std::boolalpha << result << nl;

        // std::cout << "nearest im " << found_impl << " " << i2dtree<point>::find_nearest_linear_pnt ( p, points ) << nl;
    }

    return EXIT_SUCCESS;
}



int wmain4515 ( ) {

    // std::vector<point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<point> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << nl << tree << nl << nl;

    point ptf { 7.6f, 7.9f };

    std::cout << nl << nl << "nearest " << nl << tree.find_nearest_pnt ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << i2dtree<point>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    // tree.add_point ( ptf );

    // std::cout << nl << tree << nl << nl;

    return EXIT_SUCCESS;
}

