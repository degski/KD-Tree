
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
#include <sax/iostream.hpp>
#include <type_traits>

#include <sax/stl.hpp>

namespace sax {
namespace detail {

constexpr std::size_t linear_bound = 1u;

template<typename>
struct signed_double_width_integer {};
template<>
struct signed_double_width_integer<std::uint8_t> {
    using type = std::int16_t;
};
template<>
struct signed_double_width_integer<std::uint16_t> {
    using type = std::int32_t;
};
template<>
struct signed_double_width_integer<std::uint32_t> {
    using type = std::int64_t;
};
template<>
struct signed_double_width_integer<std::uint64_t> {
    using type = std::int64_t;
}; // Inshallah.
template<>
struct signed_double_width_integer<std::int8_t> {
    using type = std::int16_t;
};
template<>
struct signed_double_width_integer<std::int16_t> {
    using type = std::int32_t;
};
template<>
struct signed_double_width_integer<std::int32_t> {
    using type = std::int64_t;
};
template<>
struct signed_double_width_integer<std::int64_t> {
    using type = std::int64_t;
}; // Inshallah.

template<typename T>
struct same_sized_int {
    using type = std::make_signed_t<T>;
};
template<>
struct same_sized_int<float> {
    using type = std::int32_t;
};
template<>
struct same_sized_int<double> {
    using type = std::int64_t;
};

// Integer LogN.
template<int Base, typename T, typename sfinae = std::enable_if_t<std::conjunction_v<std::is_integral<T>, std::is_unsigned<T>>>>
constexpr T iLog ( T const n_, T const p_ = T ( 0 ) ) noexcept {
    return n_ < Base ? p_ : iLog<Base, T, sfinae> ( n_ / Base, p_ + 1 );
}

// Integer Log2.
template<typename T, typename = std::enable_if_t<std::conjunction_v<std::is_integral<T>, std::is_unsigned<T>>>>
constexpr T ilog2 ( T const n_ ) noexcept {
    return iLog<2, T> ( n_ );
}

template<typename T, typename = std::enable_if_t<std::conjunction_v<std::is_integral<T>, std::is_unsigned<T>>>>
constexpr T next_power_2 ( T const n_ ) noexcept {
    return n_ > 2 ? T ( 1 ) << ( ilog2<T> ( n_ - 1 ) + 1 ) : n_;
}

template<std::size_t N>
constexpr std::size_t array_size ( ) noexcept {
    return N > detail::linear_bound ? detail::next_power_2 ( N + 1 ) - 1 : N;
}

template<std::size_t S>
struct message { // needs fixing.

    template<typename... Args>
    message ( Args... ) {
        static_assert ( not( 2 == S or 3 == S ), "2 or 3 dimensions only" );
    }
};
} // namespace detail

template<typename T>
struct Point2 {

    using value_type = T;

    value_type x, y;

    Point2 ( ) noexcept : x{ std::numeric_limits<value_type>::quiet_NaN ( ) } {};
    Point2 ( Point2 const & ) noexcept = default;
    Point2 ( Point2 && ) noexcept      = default;
    Point2 ( value_type && x_, value_type && y_ ) noexcept : x{ std::move ( x_ ) }, y{ std::move ( y_ ) } {}

    // template<typename SfmlVec>
    // Point2 ( SfmlVec && v_ ) noexcept : x{ std::move ( v_.x ) }, y{ std::move ( v_.y ) } {}

    [[maybe_unused]] Point2 & operator= ( Point2 const & ) noexcept = default;
    [[maybe_unused]] Point2 & operator= ( Point2 && ) noexcept = default;

    [[nodiscard]] bool operator== ( Point2 const & p_ ) const noexcept { return x == p_.x and y == p_.y; }
    [[nodiscard]] bool operator!= ( Point2 const & p_ ) const noexcept { return x != p_.x or y != p_.y; }

    [[maybe_unused]] Point2 & operator+= ( Point2 const & p_ ) noexcept {
        x += p_.x;
        y += p_.y;
        return *this;
    }
    [[maybe_unused]] Point2 & operator-= ( Point2 const & p_ ) noexcept {
        x -= p_.x;
        y -= p_.y;
        return *this;
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, Point2 const & p_ ) noexcept {
        if ( not std::isnan ( p_.x ) )
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

    Point3 ( ) noexcept : x{ std::numeric_limits<value_type>::quiet_NaN ( ) } {};
    Point3 ( Point3 const & ) noexcept = default;
    Point3 ( Point3 && ) noexcept      = default;
    Point3 ( value_type && x_, value_type && y_, value_type && z_ ) noexcept :
        x{ std::move ( x_ ) }, y{ std::move ( y_ ) }, z{ std::move ( z_ ) } {}

    //  template<typename SfmlVec>
    //  Point3 ( SfmlVec && v_ ) noexcept : x{ std::move ( v_.x ) }, y{ std::move ( v_.y ) }, z{ std::move ( v_.z ) } {}

    [[maybe_unused]] Point3 & operator= ( Point3 const & ) noexcept = default;
    [[maybe_unused]] Point3 & operator= ( Point3 && ) noexcept = default;

    [[nodiscard]] bool operator== ( Point3 const & p_ ) const noexcept { return x == p_.x and y == p_.y and z == p_.z; }
    [[nodiscard]] bool operator!= ( Point3 const & p_ ) const noexcept { return x != p_.x or y != p_.y or z != p_.z; }

    [[maybe_unused]] Point3 & operator+= ( Point3 const & p_ ) noexcept {
        x += p_.x;
        y += p_.y;
        z += p_.z;
        return *this;
    }
    [[maybe_unused]] Point3 & operator-= ( Point3 const & p_ ) noexcept {
        x -= p_.x;
        y -= p_.y;
        z -= p_.z;
        return *this;
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, Point3 const & p_ ) noexcept {
        if ( not std::isnan ( p_.x ) )
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

struct vector_tag_t {};
struct array_tag_t {};
/*
// Implicit KD full binary tree of dimension 2.
template<typename T, typename P = Point2<T>>
struct Tree2D {

    using value_type = P;
    using base_type  = T;
    using dist_type =
        std::conditional_t<std::is_floating_point_v<base_type>, base_type, detail::signed_double_width_integer<base_type>>;
    using pointer         = value_type *;
    using reference       = value_type &;
    using const_pointer   = value_type const *;
    using const_reference = value_type const &;

    using container = std::vector<P>;

    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:
    template<typename forward_it>
    [[nodiscard]] std::size_t get_dimensions_order ( forward_it const first_, forward_it const last_ ) const noexcept {
        auto const [ min_x, max_x ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.x < b.x; } );
        auto const [ min_y, max_y ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.y < b.y; } );
        return ( max_x->x - min_x->x ) < ( max_y->y - min_y->y );
    }

    [[nodiscard]] pointer left ( pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] pointer right ( pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer left ( const_pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer right ( const_pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }

    [[nodiscard]] bool is_leaf ( const_pointer const p_ ) const noexcept {
        return m_leaf_start < p_ or std::isnan ( left ( p_ )->x );
    }

    template<typename random_it>
    void kd_construct_xy ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        std::swap ( *p_, *median );
        std::cout << "x " << *this << nl;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        std::swap ( *p_, *median );
        std::cout << "y " << *this << nl;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const_pointer const p_ ) const noexcept {
        dist_type d = Tree2D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > dist_type{ 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const_pointer const p_ ) const noexcept {
        dist_type d = Tree2D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_to.y ) > dist_type{ 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const_pointer const ) const noexcept {
        for ( value_type const & v : m_data ) {
            dist_type const d = distance_squared ( v, m_to );
            if ( d < m_min_distance_squared ) {
                m_point                = &v;
                m_min_distance_squared = d;
            }
        }
    }

    container m_data;
    std::size_t m_size;
    const_pointer m_leaf_start = nullptr;
    void ( Tree2D::*nn_search ) ( const_pointer const ) const noexcept;
    container m_recently_added;

    // These mutable types are class global result types.
    mutable const_pointer m_point = nullptr;
    mutable value_type m_to;
    mutable dist_type m_min_distance_squared = std::numeric_limits<dist_type>::max ( );

    public:
    Tree2D ( ) noexcept {}
    Tree2D ( Tree2D const & rhs_ ) :
        m_data{ rhs_.m_data }, m_size{ rhs_.m_size }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {}
    Tree2D ( Tree2D && rhs_ ) noexcept :
        m_data{ std::move ( rhs_.m_data ) }, m_size{ rhs_.m_size }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {
    }

    Tree2D ( container const & c_ ) : m_data{ c_ } { initialize ( std::begin ( m_data ), std::end ( m_data ) ); }
    Tree2D ( container && c_ ) : m_data{ std::move ( c_ ) } { initialize ( std::begin ( m_data ), std::end ( m_data ) ); }

    [[nodiscard]] std::size_t size ( ) noexcept { return m_size; }
    [[nodiscard]] std::size_t capacity ( ) noexcept { return m_data.size ( ); }

    [[nodiscard]] bool contains ( value_type const & point_ ) const noexcept {
        m_to                   = point_;
        m_min_distance_squared = std::numeric_limits<dist_type>::max ( );
        ( this->*nn_search ) ( m_data.data ( ) );
        return *m_point == point_;
    }

    template<typename... Args>
    void emplace ( Args &&... args_ ) noexcept {
        value_type point = { std::forward<Args &&> ( args_ )... };
        if ( std::end ( m_recently_added ) ==
             std::find ( std::begin ( m_recently_added ), std::end ( m_recently_added ), point ) ) // If not recently added.
            if ( not contains ( point ) )                                                          // If not already in tree.
                m_recently_added.emplace_back ( std::move ( point ) );
    }

    void rebalance ( ) noexcept {
        if ( m_size > detail::linear_bound ) {
            if ( m_recently_added.size ( ) ) {
                {
                    // Grow container, iff required.
                    auto const skip_higher = m_data.size ( ) / 2 - 1;
                    std::cout << "size " << m_size << nl;
                    m_size += m_recently_added.size ( );
                    std::cout << "size " << m_size << nl;
                    if ( m_size > m_data.size ( ) ) {
                        std::cout << "resizing" << nl;
                        m_data.resize ( capacity ( m_size ) );
                        m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                    }
                    std::cout << *this << nl;
                    // Fill invalid points with recent emplacements.
                    auto it_nan = std::find_if ( std::begin ( m_data ) + skip_higher, std::end ( m_data ),
                                                 [] ( auto const & p ) noexcept { return std::isnan ( p.x ); } );
                    while ( m_recently_added.size ( ) ) {
                        *it_nan = m_recently_added.back ( );
                        m_recently_added.pop_back ( );
                        it_nan = std::find_if ( it_nan, std::end ( m_data ),
                                                [] ( auto const & p ) noexcept { return std::isnan ( p.x ); } );
                    }
                    std::cout << *this << nl;
                    if ( it_nan < std::prev ( std::end (
                                      m_data ) ) ) { // Check if invalid points (that is not already at the back) are present.
                        // Move invalid points to the back.
                        auto it_non_nan = std::find_if ( std::rbegin ( m_data ), std::rend ( m_data ),
                                                         [] ( auto const & p ) noexcept { return not std::isnan ( p.x ); } );
                        while ( &*it_nan < &*it_non_nan ) {
                            std::swap ( *it_nan, *it_non_nan );
                            it_nan     = std::find_if ( it_nan, std::end ( m_data ),
                                                    [] ( auto const & p ) noexcept { return std::isnan ( p.x ); } );
                            it_non_nan = std::find_if ( std::rbegin ( m_data ), it_non_nan + 1,
                                                        [] ( auto const & p ) noexcept { return not std::isnan ( p.x ); } );
                        }
                    }
                }
                {
                    // Re-balance.
                    switch ( get_dimensions_order ( std::begin ( m_data ), std::begin ( m_data ) + m_size ) ) {
                        case 0:
                            kd_construct_xy ( m_data.data ( ), std::begin ( m_data ), std::begin ( m_data ) + m_size );
                            nn_search = &Tree2D::nn_search_xy;
                            break;
                        case 1:
                            kd_construct_yx ( m_data.data ( ), std::begin ( m_data ), std::begin ( m_data ) + m_size );
                            nn_search = &Tree2D::nn_search_yx;
                            break;
                    }
                }
            }
        }
    }

    [[nodiscard]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[nodiscard]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[nodiscard]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[nodiscard]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[nodiscard]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[nodiscard]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[nodiscard]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    [[nodiscard]] bool is_valid ( iterator it_ ) noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_valid ( const_iterator it_ ) const noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( iterator it_ ) noexcept { return std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( const_iterator it_ ) const noexcept { return std::isnan ( it_->x ); }

    [[nodiscard]] static bool is_valid ( const_reference value_type_ ) noexcept { return not std::isnan ( value_type_.x ); }
    [[nodiscard]] static bool is_not_valid ( const_reference value_type_ ) noexcept { return std::isnan ( value_type_.x ); }

    Tree2D & operator= ( Tree2D const & rhs_ ) {
        m_data       = rhs_.m_data;
        m_size       = rhs_.m_size;
        m_leaf_start = rhs_.m_leaf_start;
        nn_search    = rhs_.nn_search;
        return *this;
    }
    Tree2D & operator= ( Tree2D && rhs_ ) noexcept {
        m_data       = std::move ( rhs_.m_data );
        m_size       = rhs_.m_size;
        m_leaf_start = rhs_.m_leaf_start;
        nn_search    = rhs_.nn_search;
        return *this;
    }

    template<typename size_type>
    [[nodiscard]] reference operator[] ( size_type const i_ ) noexcept {
        return m_data[ i_ ];
    }
    template<typename size_type>
    [[nodiscard]] const_reference operator[] ( size_type const i_ ) const noexcept {
        return m_data[ i_ ];
    }

    template<typename forward_it>
    void initialize ( forward_it first_, forward_it last_ ) {
        if ( first_ < last_ ) {
            m_size = static_cast<std::size_t> ( std::distance ( first_, last_ ) );
            if ( m_size > detail::linear_bound ) {
                m_data.resize ( capacity ( m_size ) );
                first_       = std::begin ( m_data );
                last_        = first_ + m_size;
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree2D::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yx ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree2D::nn_search_yx;
                        break;
                }
            }
            else {
                nn_search = &Tree2D::nn_search_linear;
            }
        }
    }

    [[nodiscard]] const_pointer nn_pointer ( value_type const & point_ ) const noexcept {
        m_to                   = point_;
        m_min_distance_squared = std::numeric_limits<dist_type>::max ( );
        ( this->*nn_search ) ( m_data.data ( ) );
        return m_point;
    }
    [[nodiscard]] sax::pair<const_pointer, base_type> nn_pointer_distance ( value_type const & point_ ) const noexcept {
        return { nn_pointer ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] std::ptrdiff_t nn_index ( value_type const & point_ ) const noexcept {
        return nn_pointer ( point_ ) - m_data.data ( );
    }
    [[nodiscard]] sax::pair<std::ptrdiff_t, base_type> nn_index_distance ( value_type const & point_ ) const noexcept {
        return { nn_index ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] base_type nn_distance ( value_type const & point_ ) const noexcept {
        assert ( m_min_distance_squared != std::numeric_limits<dist_type>::max ( ) );
        return static_cast<base_type> ( m_min_distance_squared );
    }

    [[nodiscard]] static constexpr dist_type distance_squared ( value_type const & p1_, value_type const & p2_ ) noexcept {
        return ( ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) *
                 ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) ) +
               ( ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) *
                 ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) );
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, Tree2D const & tree_ ) noexcept {
        for ( auto const & p : tree_.m_data )
            out_ << p;
        return out_;
    }

    private:
    template<typename U>
    [[nodiscard]] static constexpr U capacity ( U const i_ ) noexcept {
        assert ( i_ > 0 );
        return i_ > detail::linear_bound ? detail::next_power_2 ( i_ + 1 ) - 1 : i_;
    }
}; // namespace sax

*/

// Implicit KD full binary tree of dimension 2.
template<typename T, typename P = Point2<T>, typename Type = vector_tag_t, std::size_t N = 0>
struct Tree2D {

    using value_type = P;
    using base_type  = T;
    using dist_type =
        std::conditional_t<std::is_floating_point_v<base_type>, base_type, detail::signed_double_width_integer<base_type>>;
    using pointer         = value_type *;
    using reference       = value_type &;
    using const_pointer   = value_type const *;
    using const_reference = value_type const &;

    using container_type = Type;
    using container =
        std::conditional_t<std::is_same_v<container_type, array_tag_t>, std::array<P, detail::array_size<N> ( )>, std::vector<P>>;

    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:
    template<typename forward_it>
    [[nodiscard]] std::size_t get_dimensions_order ( forward_it const first_, forward_it const last_ ) const noexcept {
        auto const [ min_x, max_x ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.x < b.x; } );
        auto const [ min_y, max_y ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.y < b.y; } );
        return ( max_x->x - min_x->x ) < ( max_y->y - min_y->y );
    }

    [[nodiscard]] pointer left ( pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] pointer right ( pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer left ( const_pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer right ( const_pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }

    [[nodiscard]] bool is_leaf ( const_pointer const p_ ) const noexcept {
        return m_leaf_start < p_ or std::isnan ( left ( p_ )->x );
    }

    template<typename random_it>
    void kd_construct_xy ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const_pointer const p_ ) const noexcept {
        dist_type d = Tree2D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > dist_type{ 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const_pointer const p_ ) const noexcept {
        dist_type d = Tree2D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_to.y ) > dist_type{ 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const_pointer const ) const noexcept {
        for ( value_type const & v : m_data ) {
            dist_type const d = distance_squared ( m_to, v );
            if ( d < m_min_distance_squared ) {
                m_point                = &v;
                m_min_distance_squared = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start = nullptr;
    void ( Tree2D::*nn_search ) ( const_pointer const ) const noexcept;

    // These mutable types are class global result types.
    mutable const_pointer m_point = nullptr;
    mutable value_type m_to;
    mutable dist_type m_min_distance_squared = std::numeric_limits<dist_type>::max ( );

    public:
    Tree2D ( ) noexcept {}
    Tree2D ( Tree2D const & ) = delete;
    Tree2D ( Tree2D && rhs_ ) noexcept :
        m_data{ std::move ( rhs_.m_data ) }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {}

    Tree2D ( std::initializer_list<value_type> il_ ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
            assert ( il_.size ( ) == N );
        }
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{} );
                }
                else {
                    m_data.resize ( capacity<std::size_t> ( il_.size ( ) ) );
                }
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree2D::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree2D::nn_search_yx;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::copy_n ( std::begin ( il_ ), il_.size ( ), std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( il_.size ( ) );
                    std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                }
                nn_search = &Tree2D::nn_search_linear;
            }
        }
    }

    template<typename forward_it>
    Tree2D ( forward_it first_, forward_it last_ ) noexcept {
        initialize ( first_, last_ );
    }

    // Returns (constexpr) the size of the std::array, or the class template parameter N ( = 0).
    [[nodiscard]] static constexpr std::size_t size ( ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
            return detail::array_size<N> ( );
        }
        else {
            return N;
        }
    }

    [[nodiscard]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[nodiscard]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[nodiscard]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[nodiscard]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[nodiscard]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[nodiscard]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[nodiscard]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    [[nodiscard]] bool is_valid ( iterator it_ ) noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_valid ( const_iterator it_ ) const noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( iterator it_ ) noexcept { return std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( const_iterator it_ ) const noexcept { return std::isnan ( it_->x ); }

    [[nodiscard]] static bool is_valid ( const_reference value_type_ ) noexcept { return not std::isnan ( value_type_.x ); }
    [[nodiscard]] static bool is_not_valid ( const_reference value_type_ ) noexcept { return std::isnan ( value_type_.x ); }

    Tree2D & operator= ( Tree2D const & ) = delete;
    Tree2D & operator                     = ( Tree2D && rhs_ ) noexcept {
        m_data       = std::move ( rhs_.m_data );
        m_leaf_start = rhs_.m_leaf_start;
        nn_search    = rhs_.nn_search;
        return *this;
    }

    template<typename size_type>
    [[nodiscard]] reference operator[] ( size_type const i_ ) noexcept {
        return m_data[ i_ ];
    }
    template<typename size_type>
    [[nodiscard]] const_reference operator[] ( size_type const i_ ) const noexcept {
        return m_data[ i_ ];
    }

    template<typename forward_it>
    void initialize ( forward_it const first_, forward_it const last_ ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
            assert ( std::distance ( first_, last_ ) == N );
        }
        if ( first_ < last_ ) {
            auto const n = std::distance ( first_, last_ );
            if ( n > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{} );
                }
                else {
                    m_data.resize ( capacity<std::size_t> ( static_cast<std::size_t> ( n ) ) );
                }
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree2D::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yx ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree2D::nn_search_yx;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::copy_n ( first_, n, std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( n );
                    std::copy ( first_, last_, std::back_inserter ( m_data ) );
                }
                nn_search = &Tree2D::nn_search_linear;
            }
        }
    }

    [[nodiscard]] const_pointer nn_pointer ( value_type const & point_ ) const noexcept {
        m_to                   = point_;
        m_min_distance_squared = std::numeric_limits<dist_type>::max ( );
        ( this->*nn_search ) ( m_data.data ( ) );
        return m_point;
    }
    [[nodiscard]] sax::pair<const_pointer, base_type> nn_pointer_distance ( value_type const & point_ ) const noexcept {
        return { nn_pointer ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] std::ptrdiff_t nn_index ( value_type const & point_ ) const noexcept {
        return nn_pointer ( point_ ) - m_data.data ( );
    }
    [[nodiscard]] sax::pair<std::ptrdiff_t, base_type> nn_index_distance ( value_type const & point_ ) const noexcept {
        return { nn_index ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] base_type nn_distance ( value_type const & point_ ) const noexcept {
        assert ( m_min_distance_squared != std::numeric_limits<dist_type>::max ( ) );
        return static_cast<base_type> ( m_min_distance_squared );
    }

    [[nodiscard]] static constexpr dist_type distance_squared ( value_type const & p1_, value_type const & p2_ ) noexcept {
        return ( ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) *
                 ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) ) +
               ( ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) *
                 ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) );
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, Tree2D const & tree_ ) noexcept {
        for ( auto const & p : tree_.m_data )
            out_ << p;
        return out_;
    }

    private:
    template<typename U>
    [[nodiscard]] static constexpr U capacity ( U const i_ ) noexcept {
        assert ( i_ > 0 );
        return i_ > detail::linear_bound ? detail::next_power_2 ( i_ + 1 ) - 1 : i_;
    }
};

// Implicit KD full binary tree of dimension 3.
template<typename T, typename P = Point3<T>, typename Type = vector_tag_t, std::size_t N = 0>
struct Tree3D {

    using value_type = P;
    using base_type  = T;
    using dist_type =
        std::conditional_t<std::is_floating_point_v<base_type>, base_type, detail::signed_double_width_integer<base_type>>;
    using pointer         = value_type *;
    using reference       = value_type &;
    using const_pointer   = value_type const *;
    using const_reference = value_type const &;

    using container_type = Type;
    using container =
        std::conditional_t<std::is_same_v<container_type, array_tag_t>, std::array<P, detail::array_size<N> ( )>, std::vector<P>>;

    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:
    template<typename forward_it>
    [[nodiscard]] std::size_t get_dimensions_order ( forward_it const first_, forward_it const last_ ) const noexcept {
        auto const [ min_x, max_x ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.x < b.x; } );
        auto const [ min_y, max_y ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.y < b.y; } );
        auto const [ min_z, max_z ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.z < b.z; } );
        sax::pair<base_type, detail::same_sized_int<base_type>> dx{ max_x->x - min_x->x, 0 }, dy{ max_y->y - min_y->y, 1 },
            dz{ max_z->z - min_z->z, 2 };
        // sort list of 3.
        if ( dx.first < dy.first )
            std::swap ( dx, dy );
        if ( dx.first < dz.first )
            std::swap ( dx, dz );
        if ( dy.first < dz.first )
            std::swap ( dy, dz );
        // decide xyz- or xzy-order.
        return ( ( dx.second == 0 and dy.second == 1 ) or ( dx.second == 1 and dy.second == 2 ) or
                 ( dx.second == 2 and dy.second == 0 ) )
                   ? dx.second
                   : 3 + dx.second;
    }

    [[nodiscard]] pointer left ( pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] pointer right ( pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer left ( const_pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer right ( const_pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }

    [[nodiscard]] bool is_leaf ( const_pointer const p_ ) const noexcept {
        return m_leaf_start < p_ or std::isnan ( left ( p_ )->x );
    }

    template<typename random_it>
    void kd_construct_xy ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yz ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yz ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yz ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_zx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_zx ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_zx ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.z < b.z; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    template<typename random_it>
    void kd_construct_xz ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_zy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_zy ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xz ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xz ( right ( p_ ), median, last_ );
        }
    }
    template<typename random_it>
    void kd_construct_zy ( pointer const p_, random_it const first_, random_it const last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.z < b.z; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > base_type{ 0 } ) {
            nn_search_yz ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yz ( right ( p_ ) );
        }
        else {
            nn_search_yz ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yz ( left ( p_ ) );
        }
    }
    void nn_search_yz ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_to.y ) > base_type{ 0 } ) {
            nn_search_zx ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_zx ( right ( p_ ) );
        }
        else {
            nn_search_zx ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_zx ( left ( p_ ) );
        }
    }
    void nn_search_zx ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - m_to.z ) > base_type{ 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_xz ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > base_type{ 0 } ) {
            nn_search_zy ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_zy ( right ( p_ ) );
        }
        else {
            nn_search_zy ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_zy ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_to.y ) > base_type{ 0 } ) {
            nn_search_xz ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xz ( right ( p_ ) );
        }
        else {
            nn_search_xz ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xz ( left ( p_ ) );
        }
    }
    void nn_search_zy ( const_pointer const p_ ) const noexcept {
        base_type d = Tree3D::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->z - m_to.z ) > base_type{ 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const_pointer const ) const noexcept {
        for ( value_type const & v : m_data ) {
            dist_type const d = distance_squared ( m_to, v );
            if ( d < m_min_distance_squared ) {
                m_point                = &v;
                m_min_distance_squared = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start = nullptr;
    void ( Tree3D::*nn_search ) ( const_pointer const ) const noexcept;

    // These mutable types are class global result types.
    mutable const_pointer m_point = nullptr;
    mutable value_type m_to;
    mutable dist_type m_min_distance_squared = std::numeric_limits<dist_type>::max ( );

    public:
    Tree3D ( ) noexcept {}
    Tree3D ( Tree3D const & ) = delete;
    Tree3D ( Tree3D && rhs_ ) noexcept :
        m_data{ std::move ( rhs_.m_data ) }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {}

    Tree3D ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{} );
                }
                else {
                    m_data.resize ( capacity<std::size_t> ( il_.size ( ) ) );
                }
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yz ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_yz;
                        break;
                    case 2:
                        kd_construct_zx ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_zx;
                        break;
                    case 3:
                        kd_construct_xz ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_xz;
                        break;
                    case 4:
                        kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_yx;
                        break;
                    case 5:
                        kd_construct_zy ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &Tree3D::nn_search_zy;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::copy_n ( std::begin ( il_ ), il_.size ( ), std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( il_.size ( ) );
                    std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                }
                nn_search = &Tree3D::nn_search_linear;
            }
        }
    }

    template<typename forward_it>
    Tree3D ( forward_it first_, forward_it last_ ) noexcept {
        initialize ( first_, last_ );
    }

    [[nodiscard]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[nodiscard]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[nodiscard]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[nodiscard]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[nodiscard]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[nodiscard]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[nodiscard]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    [[nodiscard]] bool is_valid ( iterator it_ ) noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_valid ( const_iterator it_ ) const noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( iterator it_ ) noexcept { return std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( const_iterator it_ ) const noexcept { return std::isnan ( it_->x ); }

    [[nodiscard]] static bool is_valid ( const_reference value_type_ ) noexcept { return not std::isnan ( value_type_.x ); }
    [[nodiscard]] static bool is_not_valid ( const_reference value_type_ ) noexcept { return std::isnan ( value_type_.x ); }

    Tree3D & operator= ( Tree3D const & ) = delete;
    Tree3D & operator                     = ( Tree3D && rhs_ ) noexcept {
        m_data       = std::move ( rhs_.m_data );
        m_leaf_start = rhs_.m_leaf_start;
        nn_search    = rhs_.nn_search;
        return *this;
    }

    template<typename size_type>
    [[nodiscard]] reference operator[] ( size_type const i_ ) noexcept {
        return m_data[ i_ ];
    }
    template<typename size_type>
    [[nodiscard]] const_reference operator[] ( size_type const i_ ) const noexcept {
        return m_data[ i_ ];
    }

    template<typename forward_it>
    void initialize ( forward_it const first_, forward_it const last_ ) noexcept {
        if ( first_ < last_ ) {
            auto const n = std::distance ( first_, last_ );
            if ( n > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{} );
                }
                else {
                    m_data.resize ( capacity<std::size_t> ( static_cast<std::size_t> ( n ) ) );
                }
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yz ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_yz;
                        break;
                    case 2:
                        kd_construct_zx ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_zx;
                        break;
                    case 3:
                        kd_construct_xz ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_xz;
                        break;
                    case 4:
                        kd_construct_yx ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_yx;
                        break;
                    case 5:
                        kd_construct_zy ( m_data.data ( ), first_, last_ );
                        nn_search = &Tree3D::nn_search_zy;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag_t> ) {
                    std::copy_n ( first_, n, std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( n );
                    std::copy ( first_, last_, std::back_inserter ( m_data ) );
                }
                nn_search = &Tree3D::nn_search_linear;
            }
        }
    }

    [[nodiscard]] const_pointer nn_pointer ( value_type const & point_ ) const noexcept {
        m_to                   = point_;
        m_min_distance_squared = std::numeric_limits<dist_type>::max ( );
        ( this->*nn_search ) ( m_data.data ( ) );
        return m_point;
    }
    [[nodiscard]] sax::pair<const_pointer, base_type> nn_pointer_distance ( value_type const & point_ ) const noexcept {
        return { nn_pointer ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] std::ptrdiff_t nn_index ( value_type const & point_ ) const noexcept {
        return nn_pointer ( point_ ) - m_data.data ( );
    }
    [[nodiscard]] sax::pair<std::ptrdiff_t, base_type> nn_index_distance ( value_type const & point_ ) const noexcept {
        return { nn_index ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] base_type nn_distance ( value_type const & point_ ) const noexcept {
        assert ( m_min_distance_squared != std::numeric_limits<dist_type>::max ( ) );
        return static_cast<base_type> ( m_min_distance_squared );
    }

    [[nodiscard]] static constexpr base_type distance_squared ( value_type const & p1_, value_type const & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) +
               ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) + ( ( p1_.z - p2_.z ) * ( p1_.z - p2_.z ) ) );
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, Tree3D const & tree_ ) noexcept {
        for ( auto const & p : tree_.m_data )
            out_ << p;
        return out_;
    }

    private:
    template<typename U>
    [[nodiscard]] static constexpr U capacity ( U const i_ ) noexcept {
        assert ( i_ > 0 );
        return i_ > detail::linear_bound ? detail::next_power_2 ( i_ + 1 ) - 1 : i_;
    }
};

template<typename base_type, std::size_t S>
using ikdtree = typename std::conditional<2 == S, Tree2D<base_type>,
                                          typename std::conditional<3 == S, Tree3D<base_type>, detail::message<S>>::type>::type;

} // namespace sax
