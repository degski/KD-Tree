
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
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

#include <SFML/System.hpp>

using Point = sf::Vector2f;

template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const Point & p_ ) noexcept {
    out_ << '<' << p_.x << ' ' << p_.y << '>';
    return out_;
}

#define nl '\n'

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


using Int = std::int32_t;

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

template<typename Tree, typename ForwardIt>
typename Tree::NodeID construct_kdtree_impl ( ForwardIt first_, ForwardIt last_, const bool x_dim_, Tree & tree_ ) noexcept {
    if ( first_ == last_ ) {
        return Tree::NodeID::invalid;
    }
    if ( x_dim_ ) { // <7 2><5 4><2 3><4 7><9 6><8 1>
        std::sort ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
    }
    else {
        std::sort ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
    }
    const ForwardIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
    const typename Tree::NodeID node { tree_.addNode ( * median ) };
    std::cout << ( int ) ( node.value - 1 ) << ' ' << ( *median ) << nl;
    tree_ [ node ].left   = construct_kdtree_impl (               first_, median, not ( x_dim_ ), tree_ );
    tree_ [ node ].right  = construct_kdtree_impl ( std::next ( median ),  last_, not ( x_dim_ ), tree_ );
    return node;
}

template<typename Tree, typename Points>
void construct_kdtree ( Tree & tree_, Points points_ ) noexcept {
    construct_kdtree_impl ( std::begin ( points_ ), std::end ( points_ ), true, tree_ );
}


template<typename T>
[[ nodiscard ]] T distance_squared ( const sf::Vector2<T> & p1_, const sf::Vector2<T> & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
}


template<typename Tree, typename Points, typename T = typename Points::Point::value_type>
void find_nearest (
    const Tree & tree_,
    typename Tree::NodeID node_,
    const Point & point_,                      // Looking for closest node to this point.
    Point & closest_,   // Closest node (so far).
    T & min_distance_,
    const bool x_dim_ ) noexcept {

    if ( tree_.isLeaf ( node_ ) ) {
        const T dist = distance_squared ( point_, tree_ [ node_ ].data );
        if ( dist < min_distance_ ) {
            closest_ = tree_ [ node_ ].data;
            min_distance_ = dist;
        }
    }
    else {
        const T value { x_dim_ ? point_.x : point_.y }, pivot { x_dim_ ? tree_ [ node_ ].data.x : tree_ [ node_ ].data.y };
        if ( value < pivot ) {
            // Search left first.
            find_nearest ( tree_, tree_ [ node_ ].left, point_, closest_, min_distance_, not ( x_dim_ ) );
            if ( value + min_distance_ >= pivot ) {
                find_nearest ( tree_, tree_ [ node_ ].right, point_, closest_, min_distance_, not ( x_dim_ ) );
            }
        }
        else {
            // Search right first.
            find_nearest ( tree_, tree_ [ node_ ].right, point_, closest_, min_distance_, not ( x_dim_ ) );
            if ( value - min_distance_ <= pivot ) {
                find_nearest ( tree_, tree_ [ node_ ].left, point_, closest_, min_distance_, not ( x_dim_ ) );
            }
        }
    }
}



// Implicit full binary tree of size N.
template<typename T, std::size_t N, typename Int = std::int32_t>
struct Imp2DTree {

    // https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030


    using base_type = decltype ( T { }.x );
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using const_pointer = T const *;
    using const_reference = T const &;

    using container = std::array<T, N>;
    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:

    template<typename RandomIt>
    void construct ( const pointer node_, RandomIt first_, RandomIt last_, const bool x_dim_ ) noexcept {
        if ( first_ == last_ ) {
            return;
        }
        if ( x_dim_ ) {
            std::sort ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        }
        else {
            std::sort ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        }
        const RandomIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        *node_ = *median;
        construct ( left  ( node_ ),               first_, median, not ( x_dim_ ) );
        construct ( right ( node_ ), std::next ( median ),  last_, not ( x_dim_ ) );
    }

    container m_data;

    public:

    Imp2DTree ( ) noexcept {
    }
    Imp2DTree ( Imp2DTree && ) noexcept = delete;
    Imp2DTree ( const Imp2DTree & t_ ) noexcept {
        m_data = t_.m_data;
    }
    Imp2DTree ( std::initializer_list<T> il_ ) noexcept {
        assert ( il_.size ( ) <= N );
        container points;
        std::copy ( std::begin ( il_ ), std::end ( il_ ), std::begin ( points ) );
        construct ( m_data.data ( ), std::begin ( points ), std::begin ( points ) + il_.size ( ), true );
    }
    template<typename ForwardIt>
    Imp2DTree ( ForwardIt first_, ForwardIt last_ ) noexcept {
        assert ( ( last_ - first_ ) <= N );
        construct ( m_data.data ( ), first_, last_, true );
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
        return ( p_ - m_data.data ( ) ) / ( N / 2 );
    }
    [[ nodiscard ]] bool is_leaf ( const_pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) / ( N / 2 );
    }
    [[ nodiscard ]] static constexpr bool is_leaf ( const Int i_ ) noexcept {
        assert ( N >= i_ );
        return i_ / ( N / 2 );
    }

    [[ nodiscard ]] static base_type distance_squared ( const T & p1_, const T & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }
    [[ nodiscard ]] static base_type distance_squared ( const T & p1_, const_pointer p2_ ) noexcept {
        return distance_squared ( p1_, *p2_ );
    }

    private:

    struct Nearest {

        Point point, found;
        base_type min_distance;

        Nearest ( const Point & point_ ) noexcept :
            point { point_ },
            found { },
            min_distance { std::numeric_limits<base_type>::max ( ) } {
        }
    };

    void find_nearest_impl ( const const_pointer node_, const bool x_dim_, Nearest & nearest_ ) const noexcept {
        if ( is_leaf ( node_ ) ) {
            const base_type distance = Imp2DTree::distance_squared ( nearest_.point, node_ );
            if ( distance < nearest_.min_distance ) {
                nearest_.found = * node_;
                nearest_.min_distance = distance;
            }
        }
        else {
            const base_type value { x_dim_ ? nearest_.point.x : nearest_.point.y }, pivot { x_dim_ ? node_->x : node_->y };
            if ( value < pivot ) { // Search left first.
                find_nearest_impl ( left ( node_ ), not ( x_dim_ ), nearest_ );
                if ( value + nearest_.min_distance >= pivot ) {
                    find_nearest_impl ( right ( node_ ), not ( x_dim_ ), nearest_ );
                }
            }
            else { // Search right first.
                find_nearest_impl ( right ( node_ ), not ( x_dim_ ), nearest_ );
                if ( value - nearest_.min_distance <= pivot ) {
                    find_nearest_impl ( left ( node_ ), not ( x_dim_ ), nearest_ );
                }
            }
        }
    }

    public:

    [[ nodiscard ]] Point find_nearest ( const Point & point_ ) const noexcept {
        Nearest nearest { point_ };
        find_nearest_impl ( m_data.data ( ), true, nearest );
        std::cout << nearest.min_distance << ' ';
        return nearest.found;
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



Int wmain ( ) {

    // std::vector<Point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<Point> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    Imp2DTree<Point, bin_tree_size ( 13 )> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << tree << nl;

    Point ptf { 6, 5 };

    std::cout << tree.find_nearest ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << Imp2DTree<Point, bin_tree_size ( 13 )>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}



Int wmain7867867 ( ) {

    std::vector<Point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };

    for ( auto p : points ) {
        std::cout << '<' << p.x << ' ' << p.y << '>';
    }
    std::cout << '\n';

    fbt::BinTree<Point> tree;

    using It = typename fbt::BinTree<Point>::iterator;

    construct_kdtree ( tree, points );

    std::cout << tree.nodeNum ( ) << '\n';

    for ( It it { tree }; it.is_valid ( ); ++it ) {
        std::cout << '<' << it->data.x << ' ' << it->data.y << '>';
    }

    return EXIT_SUCCESS;
}
