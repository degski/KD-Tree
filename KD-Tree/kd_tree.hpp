
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

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <array>
#include <limits>

#include <SFML/System.hpp>


// https://stackoverflow.com/questions/1627305/nearest-neighbor-k-d-tree-wikipedia-proof/37107030#37107030


namespace sf {

template<typename T = float>
struct KDNode2 {
    KDNode2 ( ) = delete;
    KDNode2 ( const KDNode2 & ) = delete;
    // Branch node.
    KDNode2 ( const T split_, std::unique_ptr<const KDNode2> & lhs_, std::unique_ptr<const KDNode2> & rhs_ ) :
        m_left ( std::move ( lhs_ ) ),
        m_right ( std::move ( rhs_ ) ),
        m_pivot ( split_ ) {
    };
    // Leaf node.
    KDNode2 ( std::shared_ptr<const sf::Vector2<T>> n_ ) :
        m_leaf { n_ },
        m_pivot { T { 0 } } {
    };

    ~KDNode2 ( ) = default;

    KDNode2 & operator = ( const KDNode2 & ) = delete;

    [[ nodiscard ]] bool is_leaf ( ) const noexcept {
        return static_cast<bool> ( m_leaf );
    }

    const std::unique_ptr<const KDNode2> m_left, m_right;
    const std::shared_ptr<const sf::Vector2<T>> m_leaf;
    const T m_pivot;
};


template<typename T>
struct KDTree2 {
    KDTree2 ( ) = delete;
    KDTree2 ( const KDTree2 & ) = delete;
    KDTree2 ( KDTree2 && t_ ) :
        m_root ( std::move ( const_cast<std::unique_ptr<const KDNode2<T>>&> ( t_.m_root ) ) ) {
    };

    ~KDTree2 ( ) = default;

    KDTree2 & operator = ( const KDTree2 & ) = delete;
    [[ maybe_unused ]] KDTree2 & operator = ( KDTree2 && t_ ) noexcept {
        m_root = std::move ( const_cast<std::unique_ptr<const KDNode2<T>>&> ( t_.m_root ) );
    }

    [[ nodiscard ]] const sf::Vector2<T> find_nearest ( const sf::Vector2<T> & point_ ) const noexcept;

    // Nearest neighbour search - runs in O(log n).
    void find_nearest (
        const std::unique_ptr<const KDNode2<T>> & node_,
        const sf::Vector2<T> & point_,
        const sf::Vector2<T> & closest_,
        T & min_dist_,
        const std::size_t depth_ = 0u ) const noexcept;

    const std::unique_ptr<const KDNode2<T>> m_root;
};



template<typename T>
[[ nodiscard ]] T distance_squared ( const sf::Vector2<T> & p1, const sf::Vector2<T> & p2 ) noexcept {
    return ( ( p1.x - p2.x ) * ( p1.x - p2.x ) ) + ( ( p1.y - p2.y ) * ( p1.y - p2.y ) );
}


template<typename T>
void KDTree2<T>::find_nearest (
    const std::unique_ptr<const KDNode2<T>> & node_,
    const sf::Vector2<T> & point_,                      // Looking for closest node to this point.
    const sf::Vector2<T> & closest_,   // Closest node (so far).
    T & min_dist_,
    const std::size_t depth_ ) const noexcept {
    if ( node_->is_leaf ( ) ) {
        const T dist = distance_squared ( point, node_->m_leaf->m_point );
        if ( dist < min_dist_ ) {
            closest_ = *( node_->m_leaf );
            min_dist_ = dist;
        }
    }
    else {
        const T val { depth_ % 2 ? point_.y : point_.x };
        if ( val < node_->m_pivot ) {
            // Search left first.
            find_nearest ( node_->m_left, point_, closest_, min_dist_, depth_ + 1 );
            if ( val + min_dist_ >= node_->m_pivot ) {
                find_nearest ( node_->m_right, point_, closest_, min_dist_, depth_ + 1 );
            }
        }
        else {
            // Search right first.
            find_nearest ( node_->m_right, point_, closest_, min_dist_, depth_ + 1 );
            if ( val - min_dist_ <= node_->m_pivot ) {
                find_nearest ( node_->m_left, point_, closest_, min_dist_, depth_ + 1 );
            }
        }
    }
}

// Nearest neighbour.
template<typename T>
[[ nodiscard ]] const sf::Vector2<T> KDTree2<T>::find_nearest ( const sf::Vector2<T> & point_ ) const noexcept {
    const sf::Vector2<T> closest;
    T min_dist = std::numeric_limits<T>::max ( );
    find_nearest ( m_root, point_, closest, min_dist );
    return closest;
}

}


/*

http://lith.me/code/2015/06/08/Nearest-Neighbor-Search-with-KDTree

private static Node findNearest(final Node current, final Node target, final int depth) {
        final int axis = depth % K;
        final int direction = getComparator(axis).compare(target, current);
        final Node next = (direction < 0) ? current.left : current.right;
        final Node other = (direction < 0) ? current.right : current.left;
        Node best = (next == null) ? current : findNearest(next, target, depth + 1);
        if (current.euclideanDistance(target) < best.euclideanDistance(target)) {
            best = current;
        }
        if (other != null) {
            if (current.verticalDistance(target, axis) < best.euclideanDistance(target)) {
                final Node possibleBest = findNearest(other, target, depth + 1);
                if (possibleBest.euclideanDistance(target) < best.euclideanDistance(target)) {
                    best = possibleBest;
                }
            }
        }
        return best;
    }


*/