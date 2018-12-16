
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
#include <splitmix.hpp>
#include <plf/plf_nanotimer.h>

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


#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

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
        const RandomIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        if ( x_dim_ ) {
            std::nth_element ( first_, median, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        }
        else {
            std::nth_element ( first_, median, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        }
        *node_ = *median;
        construct ( left  ( node_ ),               first_, median, not ( x_dim_ ) );
        if ( 1 == std::distance ( first_, median ) ) {
            *right ( node_ ) = Point { std::numeric_limits<base_type>::min ( ), std::numeric_limits<base_type>::min ( ) };
        }
        construct ( right ( node_ ), std::next ( median ),  last_, not ( x_dim_ ) );
    }

    container m_data;
    const bool m_dim_start;
    const std::size_t m_dim_mask;

    public:

    Imp2DTree ( ) noexcept {
    }
    Imp2DTree ( Imp2DTree && ) noexcept = delete;
    Imp2DTree ( const Imp2DTree & t_ ) noexcept {
        m_data = t_.m_data;
        m_dim_mask = t_.m_dim_mask;
    }
    Imp2DTree ( std::initializer_list<T> il_ ) noexcept :
        m_dim_start { pick_dimension ( std::begin ( il_ ), std::end ( il_ ) ) },
        m_dim_mask { m_dim_start ? 0x5555'5555'5555'5555 : 0xAAAA'AAAA'AAAA'AAAA } {
        assert ( il_.size ( ) <= N );
        container points;
        std::copy ( std::begin ( il_ ), std::end ( il_ ), std::begin ( points ) );
        construct ( m_data.data ( ), std::begin ( points ), std::begin ( points ) + il_.size ( ), dim ( 1u ) );
    }
    template<typename ForwardIt>
    Imp2DTree ( ForwardIt first_, ForwardIt last_ ) noexcept :
        m_dim_start { pick_dimension ( first_, last_ ) },
        m_dim_mask { m_dim_start ? 0x5555'5555'5555'5555 : 0xAAAA'AAAA'AAAA'AAAA } {
        assert ( ( last_ - first_ ) <= N );
        construct ( m_data.data ( ), first_, last_, dim ( 1u ) );
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
        return ( p_ - m_data.data ( ) ) / ( N / 2 );
    }
    [[ nodiscard ]] bool is_leaf ( const_pointer p_ ) const noexcept {
        assert ( N >= ( p_ - m_data.data ( ) ) );
        return ( p_ - m_data.data ( ) ) / ( N / 2 );
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

        const Point point;
        const_pointer found;
        base_type min_distance;

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

    void nearest_impl ( const const_pointer p_, Nearest & n_, bool dim_ ) const noexcept {
        const float d = Imp2DTree::distance_squared ( n_.point, *p_ );
        if ( d < n_.min_distance ) {
            n_.min_distance = d;
            n_.found = p_;
        }
        if ( is_leaf ( p_ ) ) {
            std::cout << "leaf " << *p_ << nl;
            return;
        }
        std::cout << "itnl " << *p_ << nl;
        const float dx = dim_ ? p_->x - n_.point.x : p_->y - n_.point.y;
        dim_ = not ( dim_ );
        nearest_impl ( dx > base_type { 0 } ? left ( p_ ) : right ( p_ ), n_, dim_ );
        if ( ( dx * dx ) >= n_.min_distance ) {
            return;
        }
        nearest_impl ( dx > base_type { 0 } ? right ( p_ ) : left ( p_ ), n_, dim_ );
    }

    public:

    [[ nodiscard ]] static bool dim_from_index ( std::size_t x ) noexcept {
        if ( likely ( x ) )
            return not ( __builtin_clzll ( x ) & int { 1 } );
        return false;
    }

    void nns_impl ( Nearest & n_ ) const noexcept {
        std::size_t index { 1u }, level_index { 0u };

        do {
            const std::size_t node = index - 1;

            const float dx = Imp2DTree::dim_from_index ( index ) ? m_data [ node ].x - n_.point.x : m_data [ node ].y - n_.point.y;
            const bool left_first { dx > base_type { 0 } };

            const float d = Imp2DTree::distance_squared ( n_.point, m_data [ node ] );
            if ( d < n_.min_distance ) {
                n_.min_distance = d;
                n_.found = m_data.data ( ) + node;
            }

            if ( is_leaf ( node ) ) {
                std::cout << "leaf " << m_data [ node ] << " idx " << index << " dim " << Imp2DTree::dim_from_index ( index ) << " li " << level_index << '\n';
            }
            else {

                std::cout << "itnl " << m_data [ node ] << " idx " << index << " dim " << Imp2DTree::dim_from_index ( index ) << " li " << level_index << '\n';

                // test children of node
                // if any accepted then
                {
                    index <<= 1;
                    level_index <<= 1;

                    // if right child ﬁrst then
                    if ( not ( left_first ) ) {
                        ++index;
                    }
                    // if rejected one child then
                    //if ( ( dx * dx ) >= n_.min_distance ) {
                    //    ++level_index;
                     //   continue;
                    //}
                    continue;
                }
            }

            const int up = __builtin_ctzll ( ++level_index );
            index >>= up;
            index += ( 1 - 2 * ( index & std::size_t { 1u } ) );
            level_index >>= up;

        } while ( index > 1u );
    }

    void nns ( const const_pointer p_, Nearest & n_, bool dim_ ) const noexcept {
        std::size_t level_start { 1u }, level_index { 0u };
        do {
            const std::size_t node = level_start + level_index - 1;
            if ( is_leaf ( node ) ) {
                std::cout << "leaf " << m_data [ node ] << '\n';
            }
            else {
                std::cout << "itnl " << m_data [ node ] << '\n';
                level_start <<= 1;
                level_index <<= 1;
                continue;
            }
            ++level_index;
            const int up = __builtin_ctzll ( level_index );
            level_start >>= up;
            level_index >>= up;

        } while ( level_start > 1u );
    }

    [[ nodiscard ]] bool dim ( const std::size_t level_start_ ) const noexcept {
        return level_start_ & m_dim_mask;
    }

    void nnsa ( Nearest & n_ ) const noexcept {
        std::size_t level_start { 1u }, level_index { 0u }, level_index_dynamic { 0u };
        static splitmix64 at_random;
        static std::bernoulli_distribution select_right;
        do {
            const std::size_t node = level_start + level_index_dynamic - 1;
            if ( is_leaf ( node ) ) {
                std::cout << "leaf " << m_data [ node ] << " ls " << level_start << " dim " << dim ( level_start ) << " li " << level_index << " lid " << level_index_dynamic << '\n';
            }
            else {
                std::cout << "itnl " << m_data [ node ] << " ls " << level_start << " dim " << dim ( level_start ) << " li " << level_index << " lid " << level_index_dynamic << '\n';
                // test children of node
                // if any accepted then
                {
                    level_start <<= 1;
                    level_index <<= 1;
                    level_index_dynamic <<= 1;

                    // if right child ﬁrst then
                    if ( select_right ( at_random ) ) {
                        ++level_index_dynamic;
                    }
                    // if rejected one child then
                    {
                    //    ++level_index;
                    }
                    continue;
                }
            }
            ++level_index;
            const int up = __builtin_ctzll ( level_index );
            level_start >>= up;
            level_index >>= up;
            level_index_dynamic >>= up;
            level_index_dynamic += ( 1 - 2 * ( level_index_dynamic & std::size_t { 1 } ) );

        } while ( level_start > 1u );
    }

    public:

    [[ nodiscard ]] Point find_nearest ( const Point & point_ ) const noexcept {
        Nearest nearest { point_ };
        std::cout << nl << nl;
        nearest_impl ( m_data.data (), nearest, m_dim_start );
        nearest.init ( );
        std::cout << nl << nl;
        nns_impl ( nearest );
        std::cout << nl << nl;
        return {}; // *nearest.found;
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

    for ( std::size_t i = 0u; i < std::size_t { 1'000'000u }; ++i ) {
        rv += Tree::dim_from_index ( i );
    }

    std::cout << "dfi   " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::cout << rv << nl;

    return EXIT_SUCCESS;

    Point ptf { 60.0f, 20.5f };

    timer.start ( );

    const auto found = tree.find_nearest ( ptf );

    std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::cout << nl << "nearest " << found << nl;

    return EXIT_SUCCESS;
}


Int wmain ( ) {

    splitmix64 rng;
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    plf::nanotimer timer;
    double st;

    constexpr int n = 15;

    std::vector<Point> points;

    for ( int i = 0; i < n; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    timer.start ( );

    Imp2DTree<Point, bin_tree_size ( n )> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    Point ptf { 60.0f, 20.5f };

    timer.start ( );

    const auto found = tree.find_nearest ( ptf );

    std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::cout << nl << "nearest " << found << nl;

    return EXIT_SUCCESS;
}


Int wmain766878 ( ) {

    std::vector<Point> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    // std::vector<Point> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    Imp2DTree<Point, bin_tree_size ( 6 )> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << nl << tree << nl << nl;

    Point ptf { 2.1f, 2.9f };

    std::cout << nl << nl << "nearest " << nl << tree.find_nearest ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << Imp2DTree<Point, bin_tree_size ( 6 )>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}
