
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

#include "kdtree.h"

template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const point & p_ ) noexcept {
    if ( point { std::numeric_limits<decltype ( p_.x )>::max ( ), std::numeric_limits<decltype ( p_.y )>::max ( ) } != p_ ) {
        out_ << '<' << p_.x << ' ' << p_.y << '>';
    }
    else {
        out_ << "<x>";
    }
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
    [[ nodiscard ]] pointer parent ( const pointer p_ ) const noexcept {
        return const_cast<pointer> ( m_data.data ( ) + ( p_ - m_data.data ( ) - 1 ) / 2 );
    }
    [[ nodiscard ]] const_pointer left ( const const_pointer p_ ) const noexcept {
        return ( p_ + 1 ) + ( p_ - m_data.data ( ) );
    }
    [[ nodiscard ]] const_pointer right ( const const_pointer p_ ) const noexcept {
        return ( p_ + 2 ) + ( p_ - m_data.data ( ) );
    }
    [[ nodiscard ]] const_pointer parent ( const const_pointer p_ ) const noexcept {
        return m_data.data ( ) + ( p_ - m_data.data ( ) - 1 ) / 2;
    }

    [[ nodiscard ]] bool is_leaf ( const const_pointer p_ ) const noexcept {
        return m_leaf_start < p_;
    }

    template<typename random_it>
    void kd_construct_x ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const point & a, const point & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_   != median )
            kd_construct_y ( left ( p_ ), first_, median );
        if ( ++median !=  last_ )
            kd_construct_y ( right ( p_ ), median, last_ );
    }
    template<typename random_it>
    void kd_construct_y ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const point & a, const point & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_   != median )
            kd_construct_x ( left  ( p_ ), first_, median );
        if ( ++median !=  last_ )
            kd_construct_x ( right ( p_ ), median,  last_ );
    }

    void nn_search_x ( const const_pointer p_ ) const noexcept {
        base_type d = i2dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_nearest.point.x ) > base_type { 0 } ) {
            nn_search_y ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_y ( right ( p_ ) );
        }
        else {
            nn_search_y ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_y ( left ( p_ ) );
        }
    }
    void nn_search_y ( const const_pointer p_ ) const noexcept {
        base_type d = i2dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_nearest.point.y ) > base_type { 0 } ) {
            nn_search_x ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_x ( right ( p_ ) );
        }
        else {
            nn_search_x ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_x ( left ( p_ ) );
        }
    }

    container m_data;
    mutable nearest_data m_nearest;
    const_pointer m_leaf_start;
    bool m_dim_start;

    public:

    i2dtree ( const i2dtree & ) = delete;
    i2dtree ( i2dtree && ) noexcept = delete;

    i2dtree ( std::initializer_list<T> il_ ) noexcept :
        m_data { bin_tree_size<std::size_t> ( il_.size ( ) ), point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_leaf_start { m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1 },
        m_dim_start { pick_dimension ( std::begin ( il_ ), std::end ( il_ ) ) } {
        if ( il_.size ( ) ) {
            container points;
            std::copy ( std::begin ( il_ ), std::end ( il_ ), std::begin ( points ) );
            if ( m_dim_start )
                kd_construct_x ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
            else
                kd_construct_y ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
        }
    }

    template<typename forward_it>
    i2dtree ( forward_it first_, forward_it last_ ) noexcept :
        m_data { bin_tree_size<std::size_t> ( static_cast<std::size_t> ( std::distance ( first_, last_ ) ) ), point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } },
        m_leaf_start { m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1 },
        m_dim_start { pick_dimension ( first_, last_ ) } {
        if ( first_ != last_  ) {
            if ( m_dim_start )
                kd_construct_x ( m_data.data ( ), first_, last_ );
            else
                kd_construct_y ( m_data.data ( ), first_, last_ );
        }
    }

    i2dtree & operator = ( const i2dtree & ) = delete;
    i2dtree & operator = ( i2dtree && ) noexcept = delete;

    [[ nodiscard ]] const_pointer nearest_ptr ( const point & point_ ) const noexcept {
        m_nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        if ( m_dim_start )
            nn_search_x ( m_data.data ( ) );
        else
            nn_search_y ( m_data.data ( ) );
        return m_nearest.found;
    }

    [[ nodiscard ]] std::ptrdiff_t nearest_idx ( const point & point_ ) const noexcept {
        return nearest_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] point nearest_pnt ( const point & point_ ) const noexcept {
        return * nearest_ptr ( point_ );
    }

    [[ nodiscard ]] static constexpr base_type distance_squared ( const point & p1_, const point & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const i2dtree & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    void add_point ( const point & p_ ) {
        if ( m_data.size ( ) == m_data.capacity ( ) ) {
            m_data.resize ( ( m_data.size ( ) + 1 ) * 2 - 1, point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
            m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
        }
        auto last = std::partition ( std::begin ( m_data ), std::end ( m_data ), [ ] ( const point & p ) { return p != point { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) }; } );
        *last++ = p_;
        i2dtree tmp { std::begin ( m_data ), last };
        std::swap ( m_data, tmp.m_data );
        m_dim_start = tmp.m_dim_start;
    }

    [[ nodiscard ]] static point nearest_linear_pnt ( const point & point_, const std::vector<point> & points_ ) noexcept {
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

    const auto found = tree.nearest_recursive ( point );

    std::cout << nl << "nearest " << found << nl;

    return EXIT_SUCCESS;
}


*/

int wmain8797 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    for ( int i = 0; i < 100; ++i ) {

        plf::nanotimer timer;
        double st;

        constexpr int n = 1'000;

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
            bool r = tree.nearest_pnt ( p ) == i2dtree<point>::nearest_linear_pnt ( p, points );
            if ( not ( r ) ) {
                const point p1 = tree.nearest_pnt ( p ), p2 = i2dtree<point>::nearest_linear_pnt ( p, points );
                const float f1 = i2dtree<point>::distance_squared ( p, p1 ), f2 = i2dtree<point>::distance_squared ( p, p2 );
                if ( f1 == f2 ) {
                    continue;
                }
                std::cout << p1 << f1 << p2 << f2 << nl;
            }
            result &= r;
        }
        // std::cout << "elapsed im " << ( std::uint64_t ) ( timer.get_elapsed_ns ( ) / cnt ) << " ns" << nl;

        std::cout << std::boolalpha << result << nl;

        // std::cout << "nearest im " << found_impl << " " << i2dtree<point>::nearest_linear_pnt ( p, points ) << nl;
    }

    return EXIT_SUCCESS;
}



int wmain897 ( ) {

    // std::vector<point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<point> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << nl << tree << nl << nl;

    point ptf { 7.6f, 7.9f };

    std::cout << nl << nl << "nearest " << nl << tree.nearest_pnt ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << i2dtree<point>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    tree.add_point ( ptf );

    std::cout << nl << tree << nl << nl;

    std::cout << nl << nl << "nearest " << nl << tree.nearest_pnt ( point { 7.7f, 8.0f } ) << nl;

    point ptf2 { 1.6f, 9.9f };

    tree.add_point ( ptf2 );

    std::cout << nl << tree << nl << nl;

    point ptf3 { 4.1f, 6.0f };

    tree.add_point ( ptf3 );

    std::cout << nl << tree << nl << nl;

    return EXIT_SUCCESS;
}


using namespace std;

const int k = 2;

// A structure to represent node of kd tree
struct Node {
    int point [ k ]; // To store k dimensional point
    Node *left, *right;
};

// A method to create a node of K D tree
struct Node* newNode ( int arr [ ] ) {
    struct Node* temp = new Node;

    for ( int i = 0; i < k; i++ )
        temp->point [ i ] = arr [ i ];

    temp->left = temp->right = NULL;
    return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec ( Node *root, int point [ ], unsigned depth ) {
    // Tree is empty?
    if ( root == NULL )
        return newNode ( point );

    // Calculate current dimension (cd) of comparison
    unsigned cd = depth % k;

    // Compare the new point with root on current dimension 'cd'
    // and decide the left or right subtree
    if ( point [ cd ] < ( root->point [ cd ] ) )
        root->left = insertRec ( root->left, point, depth + 1 );
    else
        root->right = insertRec ( root->right, point, depth + 1 );

    return root;
}

// Function to insert a new point with given point in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node* insert ( Node *root, int point [ ] ) {
    return insertRec ( root, point, 0 );
}

// A utility function to find minimum of three integers
Node *minNode ( Node *x, Node *y, Node *z, int d ) {
    Node *res = x;
    if ( y != NULL && y->point [ d ] < res->point [ d ] )
        res = y;
    if ( z != NULL && z->point [ d ] < res->point [ d ] )
        res = z;
    return res;
}

// Recursively finds minimum of d'th dimension in KD tree
// The parameter depth is used to determine current axis.
Node *findMinRec ( Node* root, int d, unsigned depth ) {
    // Base cases
    if ( root == NULL )
        return NULL;

    // Current dimension is computed using current depth and total
    // dimensions (k)
    unsigned cd = depth % k;

    // Compare point with root with respect to cd (Current dimension)
    if ( cd == d ) {
        if ( root->left == NULL )
            return root;
        return findMinRec ( root->left, d, depth + 1 );
    }

    // If current dimension is different then minimum can be anywhere
    // in this subtree
    return minNode ( root,
        findMinRec ( root->left, d, depth + 1 ),
        findMinRec ( root->right, d, depth + 1 ), d );
}

// A wrapper over findMinRec(). Returns minimum of d'th dimension
Node *findMin ( Node* root, int d ) {
    // Pass current level or depth as 0
    return findMinRec ( root, d, 0 );
}

// A utility method to determine if two Points are same
// in K Dimensional space
bool arePointsSame ( int point1 [ ], int point2 [ ] ) {
    // Compare individual pointinate values
    for ( int i = 0; i < k; ++i )
        if ( point1 [ i ] != point2 [ i ] )
            return false;

    return true;
}

// Copies point p2 to p1
void copyPoint ( int p1 [ ], int p2 [ ] ) {
    for ( int i = 0; i < k; i++ )
        p1 [ i ] = p2 [ i ];
}

// Function to delete a given point 'point[]' from tree with root
// as 'root'.  depth is current depth and passed as 0 initially.
// Returns root of the modified tree.
Node *deleteNodeRec ( Node *root, int point [ ], int depth ) {
    // Given point is not present
    if ( root == NULL )
        return NULL;

    // Find dimension of current node
    int cd = depth % k;

    // If the point to be deleted is present at root
    if ( arePointsSame ( root->point, point ) ) {
        // 2.b) If right child is not NULL
        if ( root->right != NULL ) {
            // Find minimum of root's dimension in right subtree
            Node *min = findMin ( root->right, cd );

            // Copy the minimum to root
            copyPoint ( root->point, min->point );

            // Recursively delete the minimum
            root->right = deleteNodeRec ( root->right, min->point, depth + 1 );
        }
        else if ( root->left != NULL ) // same as above
        {
            Node *min = findMin ( root->left, cd );
            copyPoint ( root->point, min->point );
            root->right = deleteNodeRec ( root->left, min->point, depth + 1 );
        }
        else // If node to be deleted is leaf node
        {
            delete root;
            return NULL;
        }
        return root;
    }

    // 2) If current node doesn't contain point, search downward
    if ( point [ cd ] < root->point [ cd ] )
        root->left = deleteNodeRec ( root->left, point, depth + 1 );
    else
        root->right = deleteNodeRec ( root->right, point, depth + 1 );
    return root;
}

// Function to delete a given point from K D Tree with 'root'
Node* deleteNode ( Node *root, int point [ ] ) {
  // Pass depth as 0
    return deleteNodeRec ( root, point, 0 );
}

// Driver program to test above functions
int main8798797 ( ) {
    struct Node *root = NULL;
    int points [ ] [ k ] = { {30, 40}, {5, 25}, {70, 70},
                      {10, 12}, {50, 30}, {35, 45} };

    int n = sizeof ( points ) / sizeof ( points [ 0 ] );

    for ( int i = 0; i < n; i++ )
        root = insert ( root, points [ i ] );

    // Delet (30, 40);
    root = deleteNode ( root, points [ 0 ] );

    cout << "Root after deletion of (30, 40)\n";
    cout << root->point [ 0 ] << ", " << root->point [ 1 ] << endl;

    return 0;
}


struct KDTree {

    kdtree *ptree;

    template<typename forward_it>
    KDTree ( forward_it first_, forward_it last_ ) noexcept :
        ptree { kd_create ( 2 ) } {
        std::for_each ( first_, last_, [ this ] ( point & pos ) { kd_insertf ( ptree, ( const float * ) & pos, NULL ); } );
    }

    ~KDTree ( ) noexcept {
        kd_free ( ptree );
    }

    [[ nodiscard ]] point nearest_pnt ( const point & pos_ ) const noexcept {
        point pos;
        struct kdres * res = kd_nearestf ( ptree, ( const float * ) & pos_ );
        kd_res_itemf ( res, ( float * ) & pos );
        kd_res_free ( res );
        return pos;
    }
};

int main ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    {
        plf::nanotimer timer;
        double st;

        constexpr int n = 100'000;

        std::vector<point> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        point ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    {
        plf::nanotimer timer;
        double st;

        constexpr int n = 100'000;

        std::vector<point> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        point ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    {
        plf::nanotimer timer;
        double st;

        constexpr int n = 100'000;

        std::vector<point> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        point ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    {
        plf::nanotimer timer;
        double st;

        constexpr int n = 100'000;

        std::vector<point> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        i2dtree<point> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        point ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }
    return EXIT_SUCCESS;
}