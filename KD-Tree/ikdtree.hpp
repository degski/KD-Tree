
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

#include <iostream>
#include <limits>
#include <type_traits>

#include <SFML/System.hpp>


// https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030

namespace kdt {

using point2f = sf::Vector2<float>;
using point3f = sf::Vector3<float>;
}

template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const kdt::point2f & p_ ) noexcept {
    if ( kdt::point2f { std::numeric_limits<decltype ( p_.x )>::max ( ), std::numeric_limits<decltype ( p_.y )>::max ( ) } != p_ ) {
        out_ << '<' << p_.x << ' ' << p_.y << '>';
    }
    else {
        out_ << "<* *>";
    }
    return out_;
}
template<typename Stream>
[[ maybe_unused ]] Stream & operator << ( Stream & out_, const kdt::point3f & p_ ) noexcept {
    if ( kdt::point3f { std::numeric_limits<decltype ( p_.x )>::max ( ), std::numeric_limits<decltype ( p_.y )>::max ( ), std::numeric_limits<decltype ( p_.z )>::max ( ) } != p_ ) {
        out_ << '<' << p_.x << ' ' << p_.y << ' ' << p_.z << '>';
    }
    else {
        out_ << "<* * *>";
    }
    return out_;
}

namespace kdt {

// Implicit full binary tree of dimension 2.
template<typename T>
struct i2dtree {

    // https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030

    using base_type = T;
    using value_type = sf::Vector2<T>;
    using pointer = value_type * ;
    using reference = value_type & ;
    using const_pointer = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<value_type>;

    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:

    static_assert ( std::is_same<value_type, point2f>::value, "point is not consistently defined" );

    struct nearest_data {
        value_type point;
        const_pointer found;
        base_type min_distance;
    };

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
        return m_leaf_start < p_ or std::numeric_limits<base_type>::max ( ) == left ( p_ )->x;
    }

    template<typename random_it>
    void kd_construct_x ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_y ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_y ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_y ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( const value_type & a, const value_type & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_x ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_x ( right ( p_ ), median, last_ );
        }
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

    void nn_search_linear ( ) const noexcept {
        // Fastest up till 44 points.
        for ( auto && v : m_data ) {
            const base_type d = distance_squared ( m_nearest.point, v );
            if ( d < m_nearest.min_distance ) {
                m_nearest.found = &v;
                m_nearest.min_distance = d;
            }
        }
    }

    container m_data;
    mutable nearest_data m_nearest;
    const_pointer m_leaf_start;
    std::size_t m_dim;

    static constexpr std::size_t linear = 44u;

    public:

    i2dtree ( const i2dtree & ) = delete;
    i2dtree ( i2dtree && ) noexcept = delete;

    i2dtree ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > linear ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                m_dim = get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) );
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( m_dim ) {
                case 0: kd_construct_x ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
                case 1: kd_construct_y ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); break;
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
    i2dtree ( forward_it first_, forward_it last_ ) noexcept {
        if ( first_ < last_ ) {
            const std::size_t n = std::distance ( first_, last_ );
            if ( n > linear ) {
                m_data.resize ( bin_tree_size<std::size_t> ( static_cast<std::size_t> ( n ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                m_dim = get_dimensions_order ( first_, last_ );
                switch ( m_dim ) {
                case 0: kd_construct_x ( m_data.data ( ), first_, last_ ); break;
                case 1: kd_construct_y ( m_data.data ( ), first_, last_ ); break;
                }
            }
            else {
                m_data.reserve ( n );
                std::copy ( first_, last_, std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    i2dtree & operator = ( const i2dtree & ) = delete;
    i2dtree & operator = ( i2dtree && ) noexcept = delete;

    [[ nodiscard ]] const_pointer nearest_ptr ( const value_type & point_ ) const noexcept {
        m_nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        switch ( m_dim ) {
        case 0: nn_search_x ( m_data.data ( ) ); break;
        case 1: nn_search_y ( m_data.data ( ) ); break;
        case 2: nn_search_linear ( ); break;
        }
        return m_nearest.found;
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
    [ [ maybe_unused ] ] friend Stream & operator << ( Stream & out_, const i2dtree & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    private:

    template<typename U>
    [[ nodiscard ]] static constexpr U bin_tree_size ( const U i_ ) noexcept {
        assert ( i_ > 0 );
        if ( i_ > linear ) {
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


// Implicit full binary tree of dimension 3.
template<typename T>
struct i3dtree {

    // https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030

    using base_type = T;
    using value_type = sf::Vector3<T>;
    using pointer = value_type *;
    using reference = value_type &;
    using const_pointer = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<value_type>;

    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:

    static_assert ( std::is_same<value_type, point3f>::value, "point is not consistently defined" );

    struct nearest_data {
        value_type point;
        const_pointer found;
        base_type min_distance;
    };

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
        return m_leaf_start < p_ or std::numeric_limits<base_type>::max ( ) == left ( p_ )->x;
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
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_nearest.point.x ) > base_type { 0 } ) {
            nn_search_yz ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_yz ( right ( p_ ) );
        }
        else {
            nn_search_yz ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_yz ( left ( p_ ) );
        }
    }
    void nn_search_yz ( const const_pointer p_ ) const noexcept {
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_nearest.point.y ) > base_type { 0 } ) {
            nn_search_zx ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_zx ( right ( p_ ) );
        }
        else {
            nn_search_zx ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_zx ( left ( p_ ) );
        }
    }
    void nn_search_zx ( const const_pointer p_ ) const noexcept {
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - m_nearest.point.z ) > base_type { 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_xz ( const const_pointer p_ ) const noexcept {
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_nearest.point.x ) > base_type { 0 } ) {
            nn_search_zy ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_zy ( right ( p_ ) );
        }
        else {
            nn_search_zy ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_zy ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const const_pointer p_ ) const noexcept {
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_nearest.point.y ) > base_type { 0 } ) {
            nn_search_xz ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_xz ( right ( p_ ) );
        }
        else {
            nn_search_xz ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_xz ( left ( p_ ) );
        }
    }
    void nn_search_zy ( const const_pointer p_ ) const noexcept {
        base_type d = i3dtree::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - m_nearest.point.z ) > base_type { 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < m_nearest.min_distance ) )
                nn_search_yx ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const const_pointer ) const noexcept {
        // Fastest up till 44 points.
        for ( const auto & v : m_data ) {
            const base_type d = distance_squared ( m_nearest.point, v );
            if ( d < m_nearest.min_distance ) {
                m_nearest.found = &v;
                m_nearest.min_distance = d;
            }
        }
    }

    container m_data;
    mutable nearest_data m_nearest;
    const_pointer m_leaf_start;
    void ( i3dtree::*nn_search ) ( const const_pointer ) const noexcept;

    static constexpr std::size_t linear = 44;

    public:

    i3dtree ( const i3dtree & ) = delete;
    i3dtree ( i3dtree && ) noexcept = delete;

    i3dtree ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > linear ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) ) ) {
                case 0: kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_xy; break;
                case 1: kd_construct_yz ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_yz; break;
                case 2: kd_construct_zx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_zx; break;
                case 3: kd_construct_xz ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_xz; break;
                case 4: kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_yx; break;
                case 5: kd_construct_zy ( m_data.data ( ), std::begin ( points ), std::end ( points ) ); nn_search = & i3dtree::nn_search_zy; break;
                }
            }
            else {
                m_data.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                nn_search = & i3dtree::nn_search_linear;
            }
        }
    }

    template<typename forward_it>
    i3dtree ( forward_it first_, forward_it last_ ) noexcept {
        if ( first_ < last_ ) {
            const std::size_t n = std::distance ( first_, last_ );
            if ( n > linear ) {
                m_data.resize ( bin_tree_size<std::size_t> ( static_cast< std::size_t > ( n ) ), value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_leaf_start = m_data.data ( ) + ( m_data.size ( ) / 2 ) - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                case 0: kd_construct_xy ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_xy; break;
                case 1: kd_construct_yz ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_yz; break;
                case 2: kd_construct_zx ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_zx; break;
                case 3: kd_construct_xz ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_xz; break;
                case 4: kd_construct_yx ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_yx; break;
                case 5: kd_construct_zy ( m_data.data ( ), first_, last_ ); nn_search = & i3dtree::nn_search_zy; break;
                }
            }
            else {
                m_data.reserve ( n );
                std::copy ( first_, last_, std::back_inserter ( m_data ) );
                nn_search = & i3dtree::nn_search_linear;
            }
        }
    }

    i3dtree & operator = ( const i3dtree & ) = delete;
    i3dtree & operator = ( i3dtree && ) noexcept = delete;

    [[ nodiscard ]] const_pointer nearest_ptr ( const value_type & point_ ) const noexcept {
        m_nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        ( this->*nn_search ) ( m_data.data ( ) );
        return m_nearest.found;
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
    [ [ maybe_unused ] ] friend Stream & operator << ( Stream & out_, const i3dtree & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    private:

    template<typename U>
    [[ nodiscard ]] static constexpr U bin_tree_size ( const U i_ ) noexcept {
        assert ( i_ > 0 );
        if ( i_ > linear ) {
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

template<typename base_type, std::size_t S>
using ikdtree = typename std::conditional<2 == S, i2dtree<base_type>, typename std::conditional<3 == S, i3dtree<base_type>, std::false_type>::type>::type;

}
