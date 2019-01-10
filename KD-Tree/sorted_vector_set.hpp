
// MIT License
//
// Copyright (c) 2019 degski
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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ciso646>
#include <cassert>

#include <algorithm>
#include <functional>
#include <iostream>
#include <utility>


template<typename T, typename Container, typename Compare = std::less<T>>
class sorted_vector_set {

public:

    using reference = T&;
    using rv_reference = T&&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using size_type = std::size_t;
    using value_type = T;
    using container = Container;
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;

    sorted_vector_set ( ) noexcept { }
    sorted_vector_set ( std::initializer_list<T> init ) : m_data { std::forward<std::initializer_list<T>> ( init ) } { }
    sorted_vector_set ( const sorted_vector_set & svs ) : m_data ( svs.m_data ) { }
    sorted_vector_set ( sorted_vector_set && svs ) noexcept : m_data ( std::forward<container> ( svs.m_data ) ) { std::cout << "svs move constructed\n"; }

    ~sorted_vector_set ( ) {

        if ( not ( m_data.empty ( ) ) ) { // Has not been moved from.
            std::cout << "svs destroyed" << std::endl;
        }
        else {
            std::cout << "svs was moved from" << std::endl;
        }
    }

    [[ maybe_unused ]] sorted_vector_set & operator = ( const sorted_vector_set & svs ) { m_data = svs.m_data; return *this;  }
    [[ maybe_unused ]] sorted_vector_set & operator = ( sorted_vector_set && svs ) noexcept { m_data = std::move ( svs.m_data ); return *this; }

    // Return iterator to an element with key equivalent to key. If no such
    // element is found, past-the-end iterator is returned.
    [[ nodiscard ]] iterator find ( const_reference t_ ) noexcept {
        const iterator it = lower_bound ( t_ );
        if ( std::end ( m_data ) == it or Compare ( ) ( *it, t_ ) ) {
            return std::end ( m_data );
        }
        return it;
    }

    [[ nodiscard ]] const_iterator find ( const_reference t_ ) const noexcept {
        const const_iterator it = lower_bound ( t_ );
        if ( std::cend ( m_data ) == it or Compare ( ) ( *it, t_ ) ) {
            return std::cend ( m_data );
        }
        return it;
    }

    [[ maybe_unused ]] iterator insert ( const_reference t_ ) noexcept {
        const iterator it = lower_bound ( t_ );
        if ( it == m_data.end ( ) or Compare ( ) ( t_, *it ) ) {
            return m_data.insert ( it, t_ );
        }
        return it;
    }

    [[ maybe_unused ]] iterator insert ( const_iterator it_, const_reference t_ ) noexcept {
        return m_data.insert ( it_, t_ );
    }

    [[ maybe_unused ]] iterator insert ( const_iterator it_, rv_reference t_ ) noexcept {
        return m_data.insert ( it_, std::move ( t_ ) );
    }

    // IMPORTANT: Changes the value, i.e. the key, if this does not respect the sort-order,
    // the set will no longer be sorted (the part of the key used for sorting should be
    // equal to the value of that part of the key of the update).
    [[ maybe_unused ]] iterator insert_or_update_unsafe ( const_reference t_ ) noexcept {
        iterator it = lower_bound ( t_ );
        if ( it == m_data.end ( ) or Compare ( ) ( t_, *it ) ) {
            it = m_data.insert ( it, t_ );
        }
        else {
            *it = t_;
        }
        return it;
    }

    // IMPORTANT: Changes the value, i.e. the key AND ASSUMES THE VALUE EXISTS.
    void update_unsafe ( const_reference t_ ) noexcept {
        *lower_bound ( t_ ) = t_;
    }

    // Gives direct acces to the underlying container, irrespective of ordering.
    [[ nodiscard ]] container & container_ref ( ) noexcept { return m_data; }
    [[ nodiscard ]] const container & container_ref ( ) const noexcept { return m_data; }

    [[ nodiscard ]] bool member ( const_reference t_ ) const noexcept { return not ( not_member ( t_ ) ); }
    [[ nodiscard ]] bool not_member ( const_reference t_ ) const noexcept { return Compare ( ) ( t_, *lower_bound ( t_ ) ); }

    [[ nodiscard ]] iterator erase ( iterator p_ ) noexcept { return m_data.erase ( p_ ); }
    [[ nodiscard ]] iterator erase ( const_iterator p_ ) noexcept { return m_data.erase ( p_ ); }
    [[ nodiscard ]] iterator erase ( reference t_ ) noexcept { return _Erase ( t_ ); }
    [[ nodiscard ]] iterator erase ( const_reference t_ ) noexcept { return _Erase ( t_ ); }

    [[ nodiscard ]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[ nodiscard ]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[ nodiscard ]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[ nodiscard ]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[ nodiscard ]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[ nodiscard ]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[ nodiscard ]] const_reference front ( ) const noexcept { return m_data.front ( ); }
    [[ nodiscard ]] const_reference back ( ) const noexcept { return m_data.back ( ); }

    void reserve ( const size_type r_ ) { m_data.reserve ( r_ ); }
    void clear ( ) noexcept { m_data.clear ( ); }

    [[ nodiscard ]] size_type size ( ) const noexcept { return m_data.size ( ); }
    [[ nodiscard ]] size_type capacity ( ) const noexcept { return m_data.capacity ( ); }
    [[ nodiscard ]] bool empty ( ) const noexcept { return m_data.empty ( ); }

    // This non-const operator allows assignment to the key (i.e. the key itself is modified),
    // better make sure that the assignment does not change the sort order!!!
    [[ nodiscard ]] reference operator [ ] ( const size_type i_ ) noexcept { return m_data [ i_ ]; }
    [[ nodiscard ]] const_reference operator [ ] ( const size_type i_ ) const noexcept { return m_data [ i_ ]; }
    [[ nodiscard ]] reference nth ( const size_type i_ ) noexcept { return m_data [ i_ ]; }
    [[ nodiscard ]] const_reference nth ( const size_type i_ ) const noexcept { return m_data [ i_ ]; }
    [[ nodiscard ]] reference at ( const size_type i_ ) noexcept { return m_data.at ( i_ ); }
    [[ nodiscard ]] const_reference at ( const size_type i_ ) const noexcept { return m_data.at ( i_ ); }

    [[ nodiscard ]] iterator lower_bound ( const_reference t_ ) noexcept {
        return std::lower_bound ( std::begin ( m_data ), std::end ( m_data ), t_, Compare ( ) );
    }

    [[ nodiscard ]] const_iterator lower_bound ( const_reference t_ ) const noexcept {
        return std::lower_bound ( std::cbegin ( m_data ), std::cend ( m_data ), t_, Compare ( ) );
    }

    [[ nodiscard ]] bool operator == ( const_reference rhs ) const noexcept {
        if ( size ( ) != rhs.size ( ) ) {
            return false;
        }
        if constexpr ( std::is_pod<T>::value ) {
            std::memcmp ( m_data.data ( ), rhs.m_data.data ( ), sizeof ( T ) * size ( ) ) == 0;
        }
        else {
            std::equal ( begin ( ), end ( ), rhs.begin ( ) );
        }
    }

    [[ nodiscard ]] bool operator != ( const_reference rhs ) const noexcept {
        return not ( *this == rhs );
    }

    private:

    [[ nodiscard ]] iterator _Erase ( const_reference t_ ) noexcept {
        const iterator it = lower_bound ( t_ );
        return Compare ( ) ( t_, *it ) ? std::end ( m_data ) : m_data.erase ( it );
    }

    container m_data;
};
