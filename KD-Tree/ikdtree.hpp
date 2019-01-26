
// MIT License
//
// Copyright (c) 2018, 2019 degski
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

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <limits>
#include <type_traits>


namespace kd {

template<typename T>
struct Point2 {

    using value_type = T;

    value_type x, y;

    Point2 ( ) noexcept = default;
    Point2 ( const Point2 & ) noexcept = default;
    Point2 ( Point2 && ) noexcept = default;
    Point2 ( value_type && x_, value_type && y_ ) noexcept :
        x { std::move ( x_ ) }, y { std::move ( y_ ) } {
    }

    [[ maybe_unused ]] Point2 & operator = ( const Point2 & ) noexcept = default;
    [[ maybe_unused ]] Point2 & operator = ( Point2 && ) noexcept = default;

    [[ nodiscard ]] bool operator == ( const Point2 & p_ ) const noexcept {
        return x == p_.x and y == p_.y;
    }
    [[ nodiscard ]] bool operator != ( const Point2 & p_ ) const noexcept {
        return x != p_.x or y != p_.y;
    }

    [[ maybe_unused ]] Point2 & operator += ( const Point2 & p_ ) noexcept {
        x += p_.x; y += p_.y;
        return *this;
    }
    [[ maybe_unused ]] Point2 & operator -= ( const Point2 & p_ ) noexcept {
        x -= p_.x; y -= p_.y;
        return *this;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Point2 & p_ ) noexcept {
        if ( Point2 { std::numeric_limits<value_type>::max ( ), std::numeric_limits<value_type>::max ( ) } != p_ )
            out_ << '<' << p_.x << ' ' << p_.y << '>';
        else
            out_ << "<* *>";
        return out_;
    }
};

template<typename T>
struct Point3 {

    using value_type = T;

    value_type x, y, z;

    Point3 ( ) noexcept = default;
    Point3 ( const Point3 & ) noexcept = default;
    Point3 ( Point3 && ) noexcept = default;
    Point3 ( value_type && x_, value_type && y_, value_type && z_ ) noexcept :
        x { std::move ( x_ ) }, y { std::move ( y_ ) }, z { std::move ( z_ ) } {
    }

    [[ maybe_unused ]] Point3 & operator = ( const Point3 & ) noexcept = default;
    [[ maybe_unused ]] Point3 & operator = ( Point3 && ) noexcept = default;

    [[ nodiscard ]] bool operator == ( const Point3 & p_ ) const noexcept {
        return x == p_.x and y == p_.y and z == p_.z;
    }
    [[ nodiscard ]] bool operator != ( const Point3 & p_ ) const noexcept {
        return x != p_.x or y != p_.y or z != p_.z;
    }

    [[ maybe_unused ]] Point3 & operator += ( const Point3 & p_ ) noexcept {
        x += p_.x; y += p_.y; z += p_.z;
        return *this;
    }
    [[ maybe_unused ]] Point3 & operator -= ( const Point3 & p_ ) noexcept {
        x -= p_.x; y -= p_.y; z -= p_.z;
        return *this;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Point3 & p_ ) noexcept {
        if ( Point3 { std::numeric_limits<value_type>::max ( ), std::numeric_limits<value_type>::max ( ), std::numeric_limits<value_type>::max ( ) } != p_ )
            out_ << '<' << p_.x << ' ' << p_.y << ' ' << p_.z << '>';
        else
            out_ << "<* * *>";
        return out_;
    }
};


using Point2f = Point2<float>;
using Point2d = Point2<double>;

using Point3f = Point3<float>;
using Point3d = Point3<double>;


// Implicit KD full binary tree of dimension 2, P can be substituted by sf::Vector2<>.
template<typename T, typename P = Point2<T>>
struct Tree2D {

    using value_type = P;
    using base_type = T;
    using pointer = value_type *;
    using reference = value_type &;
    using const_pointer = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<value_type>;

    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    struct nearest_data {
        value_type to;
        const_pointer point;
        base_type distance;
    };

    private:

    template<typename forward_it>
    [[ nodiscard ]] std::size_t get_dimensions_order ( forward_it first_, forward_it last_ ) const noexcept {
        const std::pair x = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        const std::pair y = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        return ( x.second->x - x.first->x ) < ( y.second->y - y.first->y );
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
        return ( m_leaf_start < p_ ) or ( std::numeric_limits<base_type>::max ( ) == left ( p_ )->x );
    }

    template<typename random_it>
    void kd_construct_xy ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const const_pointer p_ ) const noexcept {
        base_type d = Tree2D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - nearest.to.x ) > base_type { 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const const_pointer p_ ) const noexcept {
        base_type d = Tree2D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - nearest.to.y ) > base_type { 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_linear ( ) const noexcept {
        for ( auto && v : m_data ) {
            const base_type d = distance_squared ( nearest.to, v );
            if ( d < nearest.distance ) {
                nearest.point = &v;
                nearest.distance = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start;
    std::size_t m_dim;

    static constexpr std::size_t m_linear_bound = 44u;

    public:

    mutable nearest_data nearest;

    Tree2D ( const Tree2D & ) = delete;
    Tree2D ( Tree2D && rhs_ ) noexcept :
        m_data { std::move ( rhs_.m_data ) },
        m_leaf_start { rhs_.m_leaf_start },
        m_dim { rhs_.mdim } {
    }

    Tree2D ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                m_dim = get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) );
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( m_dim ) {
                case 0: kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
                case 1: kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
                }
            }
            else {
                m_data.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    template<typename forward_it>
    Tree2D ( forward_it first_, forward_it last_ ) noexcept {
        initialize ( first_, last_ );
    }

    [[ nodiscard ]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[ nodiscard ]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[ nodiscard ]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[ nodiscard ]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[ nodiscard ]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[ nodiscard ]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[ nodiscard ]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    Tree2D & operator = ( const Tree2D & ) = delete;
    Tree2D & operator = ( Tree2D && rhs_ ) noexcept {
        m_data = std::move ( rhs_.m_data );
        m_leaf_start = rhs_.m_leaf_start;
        m_dim = rhs_.m_dim;
        return * this;
    }

    template<typename size_type>
    [[ nodiscard ]] reference operator [ ] ( const size_type i_ ) noexcept { return m_data [ i_ ]; }
    template<typename size_type>
    [[ nodiscard ]] const_reference operator [ ] ( const size_type i_ ) const noexcept { return m_data [ i_ ]; }

    template<typename forward_it>
    void initialize ( forward_it first_, forward_it last_ ) noexcept {
        if ( first_ < last_ ) {
            const std::size_t n = std::distance ( first_, last_ );
            if ( n > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( static_cast<std::size_t> ( n ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                m_dim = get_dimensions_order ( first_, last_ );
                switch ( m_dim ) {
                case 0: kd_construct_xy ( m_data.data ( ), first_, last_ ); break;
                case 1: kd_construct_yx ( m_data.data ( ), first_, last_ ); break;
                }
            }
            else {
                m_data.reserve ( n );
                std::copy ( first_, last_, std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    [[ nodiscard ]] const_pointer nearest_ptr ( const value_type & point_ ) const noexcept {
        nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        switch ( m_dim ) {
        case 0: nn_search_xy ( m_data.data ( ) ); break;
        case 1: nn_search_yx ( m_data.data ( ) ); break;
        case 2: nn_search_linear ( ); break;
        }
        return nearest.point;
    }

    [[ nodiscard ]] std::ptrdiff_t nearest_idx ( const value_type & point_ ) const noexcept {
        return nearest_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] value_type nearest_pnt ( const value_type & point_ ) const noexcept {
        return *nearest_ptr ( point_ );
    }

    [[ nodiscard ]] static constexpr base_type distance_squared ( const value_type & p1_, const value_type & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Tree2D & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    private:

    template<typename U>
    [[ nodiscard ]] static constexpr U bin_tree_size ( const U i_ ) noexcept {
        assert ( i_ > 0 );
        if ( i_ > m_linear_bound ) {
            U p = 1;
            while ( p < i_ ) {
                p += p + 1;
            }
            return p;
        }
        else {
            return i_;
        }
    }
};


template<typename T, typename K, typename M>
struct TreeMap2D {

    using base_type = T;
    using key_type = K;
    using mapped_type = M;

    using value_type = std::pair<K, M>;
    using pointer = value_type *;
    using reference = value_type &;
    using const_pointer = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<value_type>;

    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    struct nearest_data {
        key_type to;
        const_pointer point;
        base_type distance;
    };

    private:

    template<typename forward_it>
    [[ nodiscard ]] std::size_t get_dimensions_order ( forward_it first_, forward_it last_ ) const noexcept {
        const std::pair x = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.first.x < b.first.x; } );
        const std::pair y = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.first.y < b.first.y; } );
        return ( x.second->first.x - x.first->first.x ) < ( y.second->first.y - y.first->first.y );
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
        return ( m_leaf_start < p_ ) or ( std::numeric_limits<base_type>::max ( ) == left ( p_ )->first.x );
    }

    template<typename random_it>
    void kd_construct_xy ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.first.x < b.first.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.first.y < b.first.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const const_pointer p_ ) const noexcept {
        base_type d = TreeMap2D::distance_squared ( ( *p_ ).first, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->first.x - nearest.to.x ) > base_type { 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const const_pointer p_ ) const noexcept {
        base_type d = TreeMap2D::distance_squared ( ( *p_ ).first, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->first.y - nearest.to.y ) > base_type { 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_linear ( ) const noexcept {
        for ( auto && v : m_data ) {
            const base_type d = distance_squared ( nearest.to, v.first );
            if ( d < nearest.distance ) {
                nearest.point = &v;
                nearest.distance = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start;
    std::size_t m_dim;

    static constexpr std::size_t m_linear_bound = 1u;

    public:

    mutable nearest_data nearest;

    TreeMap2D ( ) { }
    TreeMap2D ( const TreeMap2D & ) = delete;
    TreeMap2D ( TreeMap2D && rhs_ ) noexcept :
        m_data { std::move ( rhs_.m_data ) },
        m_leaf_start { rhs_.m_leaf_start },
        m_dim { rhs_.m_dim } {
    }

    TreeMap2D ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), value_type { key_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) }, mapped_type { } } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                m_dim = get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) );
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( m_dim ) {
                case 0: kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
                case 1: kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
                }
            }
            else {
                m_data.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    template<typename forward_it>
    TreeMap2D ( forward_it first_, forward_it last_ ) noexcept {
        initialize ( first_, last_ );
    }

    [[ nodiscard ]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[ nodiscard ]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[ nodiscard ]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[ nodiscard ]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[ nodiscard ]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[ nodiscard ]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[ nodiscard ]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    TreeMap2D & operator = ( const TreeMap2D & ) = delete;
    TreeMap2D & operator = ( TreeMap2D && rhs_ ) noexcept {
        m_data = std::move ( rhs_.m_data );
        m_leaf_start = rhs_.m_leaf_start;
        m_dim = rhs_.m_dim;
        return *this;
    }

    template<typename size_type>
    [[ nodiscard ]] reference operator [ ] ( const size_type i_ ) noexcept { return m_data [ i_ ]; }
    template<typename size_type>
    [[ nodiscard ]] const_reference operator [ ] ( const size_type i_ ) const noexcept { return m_data [ i_ ]; }

    template<typename forward_it>
    void initialize ( forward_it first_, forward_it last_ ) noexcept {
        if ( m_data.empty ( ) ) {
            if ( first_ < last_ ) {
                const std::size_t n = std::distance ( first_, last_ );
                if ( n > m_linear_bound ) {
                    m_data.resize ( bin_tree_size<std::size_t> ( static_cast<std::size_t> ( n ) ), value_type { key_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) }, mapped_type { } } );
                    m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                    m_dim = get_dimensions_order ( first_, last_ );
                    switch ( m_dim ) {
                    case 0: kd_construct_xy ( m_data.data ( ), first_, last_ ); break;
                    case 1: kd_construct_yx ( m_data.data ( ), first_, last_ ); break;
                    }
                }
                else {
                    m_data.reserve ( n );
                    std::copy ( first_, last_, std::back_inserter ( m_data ) );
                    m_dim = 2u;
                }
            }
        }
    }

    [[ nodiscard ]] const_pointer nearest_ptr ( const key_type & point_ ) const noexcept {
        nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        switch ( m_dim ) {
        case 0: nn_search_xy ( m_data.data ( ) ); break;
        case 1: nn_search_yx ( m_data.data ( ) ); break;
        case 2: nn_search_linear ( ); break;
        }
        return nearest.point;
    }

    [[ nodiscard ]] std::ptrdiff_t nearest_idx ( const key_type & point_ ) const noexcept {
        return nearest_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] value_type nearest_pnt ( const key_type & point_ ) const noexcept {
        return *nearest_ptr ( point_ );
    }

    [[ nodiscard ]] mapped_type nearest_val ( const key_type & point_ ) const noexcept {
        return nearest_ptr ( point_ )->second;
    }

    [[ nodiscard ]] static constexpr base_type distance_squared ( const key_type & p1_, const key_type & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const TreeMap2D & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            if constexpr ( 1u == sizeof ( mapped_type ) ) {
                out_ << '<' << p.first << ' ' << static_cast<int> ( p.second ) << '>';
            }
            else {
                out_ << '<' << p.first << ' ' << p.second << '>';
            }
        }
        return out_;
    }

    private:

    template<typename U>
    [[ nodiscard ]] static constexpr U bin_tree_size ( const U i_ ) noexcept {
        assert ( i_ > 0 );
        if ( i_ > m_linear_bound ) {
            U p = 1;
            while ( p < i_ ) {
                p += p + 1;
            }
            return p;
        }
        else {
            return i_;
        }
    }
};


// Implicit KD full binary tree of dimension 3, P can be substituted by sf::Vector3<>.
template<typename T, typename P = Point3<T>>
struct Tree3D {

    using value_type = P;
    using base_type = T;
    using pointer = value_type *;
    using reference = value_type &;
    using const_pointer = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<value_type>;

    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    struct nearest_data {
        value_type to;
        const_pointer point;
        base_type distance;
    };

    private:

    template<typename forward_it>
    [[ nodiscard ]] std::size_t get_dimensions_order ( forward_it first_, forward_it last_ ) const noexcept {
        const std::pair x = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.x < b.x; } );
        const std::pair y = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.y < b.y; } );
        const std::pair z = std::minmax_element ( first_, last_, [ ] ( const auto & a, const auto & b ) { return a.z < b.z; } );
        std::pair<base_type, std::int32_t> dx { x.second->x - x.first->x, 0 }, dy { y.second->y - y.first->y, 1 }, dz { z.second->z - z.first->z, 2 };
        // sort list of 3.
        if ( dx.first < dy.first )
            std::swap ( dx, dy );
        if ( dx.first < dz.first )
            std::swap ( dx, dz );
        if ( dy.first < dz.first )
            std::swap ( dy, dz );
        // decide xyz- or xzy-order.
        return ( ( dx.second == 0 and dy.second == 1 ) or ( dx.second == 1 and dy.second == 2 ) or ( dx.second == 2 and dy.second == 0 ) ) ? dx.second : 3 + dx.second;
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
        return ( m_leaf_start < p_ ) or ( std::numeric_limits<base_type>::max ( ) == left ( p_ )->x );
    }

    template<typename random_it>
    void kd_construct_xy ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yz ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yz ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yz ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_zx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_zx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_zx ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.z < b.z; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    template<typename random_it>
    void kd_construct_xz ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_zy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_zy ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xz ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xz ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_zy ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.z < b.z; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - nearest.to.x ) > base_type { 0 } ) {
            nn_search_yz ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yz ( right ( p_ ) );
        }
        else {
            nn_search_yz ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yz ( left ( p_ ) );
        }
    }
    void nn_search_yz ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - nearest.to.y ) > base_type { 0 } ) {
            nn_search_zx ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_zx ( right ( p_ ) );
        }
        else {
            nn_search_zx ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_zx ( left ( p_ ) );
        }
    }
    void nn_search_zx ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - nearest.to.z ) > base_type { 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_xz ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - nearest.to.x ) > base_type { 0 } ) {
            nn_search_zy ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_zy ( right ( p_ ) );
        }
        else {
            nn_search_zy ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_zy ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - nearest.to.y ) > base_type { 0 } ) {
            nn_search_xz ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xz ( right ( p_ ) );
        }
        else {
            nn_search_xz ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_xz ( left ( p_ ) );
        }
    }
    void nn_search_zy ( const const_pointer p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, nearest.to );
        if ( d < nearest.distance ) {
            nearest.distance = d; nearest.point = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - nearest.to.z ) > base_type { 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < nearest.distance ) )
                nn_search_yx ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const const_pointer ) const noexcept {
        for ( const auto & v : m_data ) {
            const base_type d = distance_squared ( nearest.to, v );
            if ( d < nearest.distance ) {
                nearest.point = &v;
                nearest.distance = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start;
    void ( Tree3D::*nn_search ) ( const const_pointer ) const noexcept;

    static constexpr std::size_t m_linear_bound = 44u;

    public:

    mutable nearest_data nearest;

    Tree3D ( const Tree3D & ) = delete;
    Tree3D ( Tree3D && ) noexcept = delete;

    Tree3D ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) ) ) {
                case 0: kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_xy; break;
                case 1: kd_construct_yz ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_yz; break;
                case 2: kd_construct_zx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_zx; break;
                case 3: kd_construct_xz ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_xz; break;
                case 4: kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_yx; break;
                case 5: kd_construct_zy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & Tree3D::nn_search_zy; break;
                }
            }
            else {
                m_data.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                nn_search = & Tree3D::nn_search_linear;
            }
        }
    }

    template<typename forward_it>
    Tree3D ( forward_it first_, forward_it last_ ) noexcept {
        initialize ( first_, last_ );
    }

    [[ nodiscard ]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[ nodiscard ]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[ nodiscard ]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[ nodiscard ]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[ nodiscard ]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[ nodiscard ]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[ nodiscard ]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    Tree3D & operator = ( const Tree3D & ) = delete;
    Tree3D & operator = ( Tree3D && ) noexcept = delete;

    template<typename size_type>
    [[ nodiscard ]] reference operator [ ] ( const size_type i_ ) noexcept { return m_data [ i_ ]; }
    template<typename size_type>
    [[ nodiscard ]] const_reference operator [ ] ( const size_type i_ ) const noexcept { return m_data [ i_ ]; }

    template<typename forward_it>
    void initialize ( forward_it first_, forward_it last_ ) noexcept {
        if ( first_ < last_ ) {
            const std::size_t n = std::distance ( first_, last_ );
            if ( n > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( static_cast<std::size_t> ( n ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                case 0: kd_construct_xy ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_xy; break;
                case 1: kd_construct_yz ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_yz; break;
                case 2: kd_construct_zx ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_zx; break;
                case 3: kd_construct_xz ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_xz; break;
                case 4: kd_construct_yx ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_yx; break;
                case 5: kd_construct_zy ( m_data.data ( ), first_, last_ ); nn_search = &Tree3D::nn_search_zy; break;
                }
            }
            else {
                m_data.reserve ( n );
                std::copy ( first_, last_, std::back_inserter ( m_data ) );
                nn_search = &Tree3D::nn_search_linear;
            }
        }
    }

    [[ nodiscard ]] const_pointer nearest_ptr ( const value_type & point_ ) const noexcept {
        nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        ( this->*nn_search ) ( m_data.data ( ) );
        return nearest.point;
    }

    [[ nodiscard ]] std::ptrdiff_t nearest_idx ( const value_type & point_ ) const noexcept {
        return nearest_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] value_type nearest_pnt ( const value_type & point_ ) const noexcept {
        return *nearest_ptr ( point_ );
    }

    [[ nodiscard ]] static constexpr base_type distance_squared ( const value_type & p1_, const value_type & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) + ( ( p1_.z - p2_.z ) * ( p1_.z - p2_.z ) ) );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const Tree3D & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    private:

    template<typename U>
    [[ nodiscard ]] static constexpr U bin_tree_size ( const U i_ ) noexcept {
        assert ( i_ > 0 );
        if ( i_ > m_linear_bound ) {
            U p = 1;
            while ( p < i_ ) {
                p += p + 1;
            }
            return p;
        }
        else {
            return i_;
        }
    }
};

namespace detail {

template<std::size_t S>
struct message { // needs fixing.

    template<typename ... Args>
    message ( Args ... ) {
        static_assert ( not ( 2 == S or 3 == S ), "2 or 3 dimensions only" );
    }
};
}

template<typename base_type, std::size_t S>
using ikdtree = typename std::conditional<2 == S, Tree2D<base_type>, typename std::conditional<3 == S, Tree3D<base_type>, detail::message<S>>::type>::type;

}
