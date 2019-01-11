
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

#include <spatial/idle_point_multiset.hpp>
#include <spatial/neighbor_iterator.hpp>

#include <splitmix.hpp>
#include <plf/plf_nanotimer.h>

#include "ikdtree.hpp"
#include "kdtree.h"
#include "kdtree2.hpp"

#define nl '\n'

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)


[[ nodiscard ]] constexpr auto nn_distance_squared ( const kd::Point2f & p1_, const kd::Point2f & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
}
[[ nodiscard ]] constexpr auto nn_distance_squared ( const kd::Point3f & p1_, const kd::Point3f & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) ) + ( ( p1_.z - p2_.z ) * ( p1_.z - p2_.z ) );
}

template<typename forward_it, typename value_type>
[[ nodiscard ]] value_type nn_search_linear ( forward_it first_, forward_it last_, const value_type & p_ ) noexcept {
    using base_type = decltype ( p_.x );
    base_type min_distance = std::numeric_limits<base_type>::max ( );
    forward_it found = last_;
    while ( first_ != last_ ) {
        const base_type d = nn_distance_squared ( p_, * first_ );
        if ( d < min_distance ) {
            min_distance = d;
            found = first_;
        }
        ++first_;
    }
    return * found;
}


bool test ( const int n_ ) noexcept {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };

    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 20.0f, 40.0f };

    std::vector<kd::Point2f> points;
    points.reserve ( n_ );

    for ( int i = 0; i < n_; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

    bool rv = true;

    for ( int i = 0; i < 100'000; ++i ) {
        const kd::Point2f ptf { disx ( rng ), disy ( rng ) };
        rv = rv and ( tree.nearest_pnt ( ptf ) == nn_search_linear ( std::begin ( points ), std::end ( points ), ptf ) );
        if ( not ( rv ) ) {
            std::cout << "fail\n";
            exit ( 0 );
        }
    }

    return rv;
}

int wmain68461 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };
    std::uniform_int_distribution<int> dis { 16, 10'000 };

    for ( int i = 0; i < 1'000; ++i ) {
        std::cout << std::boolalpha << test ( dis ( rng ) ) << nl;
    }

    return EXIT_SUCCESS;
}


int wmain879879 ( ) {

    // std::vector<kd::Point2f> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<kd::Point2f> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    kd::Tree2D<float> tree { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    std::cout << nl << tree << nl << nl;

    kd::Point2f ptf { 7.6f, 7.9f };

    std::cout << nl << nl << "nearest " << nl << tree.nearest_pnt ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << kd::Tree2D<float>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}


using namespace std;

constexpr int k = 2;

// A structure to represent node of kd tree
struct Node {
    int Point2f [ k ]; // To store k dimensional Point2f
    Node *left, *right;
};

// A method to create a node of K D tree
struct Node* newNode ( int arr [ ] ) {
    struct Node* temp = new Node;

    for ( int i = 0; i < k; i++ )
        temp->Point2f [ i ] = arr [ i ];

    temp->left = temp->right = NULL;
    return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec ( Node *root, int Point2f [ ], unsigned depth ) {
    // Tree is empty?
    if ( root == NULL )
        return newNode ( Point2f );

    // Calculate current dimension (cd) of comparison
    unsigned cd = depth % k;

    // Compare the new Point2f with root on current dimension 'cd'
    // and decide the left or right subtree
    if ( Point2f [ cd ] < ( root->Point2f [ cd ] ) )
        root->left = insertRec ( root->left, Point2f, depth + 1 );
    else
        root->right = insertRec ( root->right, Point2f, depth + 1 );

    return root;
}

// Function to insert a new Point2f with given Point2f in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node* insert ( Node *root, int Point2f [ ] ) {
    return insertRec ( root, Point2f, 0 );
}

// A utility function to find minimum of three integers
Node *minNode ( Node *x, Node *y, Node *z, int d ) {
    Node *res = x;
    if ( y != NULL && y->Point2f [ d ] < res->Point2f [ d ] )
        res = y;
    if ( z != NULL && z->Point2f [ d ] < res->Point2f [ d ] )
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

    // Compare Point2f with root with respect to cd (Current dimension)
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
bool arePointsSame ( int point1 [ ], int Point2 [ ] ) {
    // Compare individual pointinate values
    for ( int i = 0; i < k; ++i )
        if ( point1 [ i ] != Point2 [ i ] )
            return false;

    return true;
}

// Copies Point2f p2 to p1
void copyPoint ( int p1 [ ], int p2 [ ] ) {
    for ( int i = 0; i < k; i++ )
        p1 [ i ] = p2 [ i ];
}

// Function to delete a given Point2f 'Point2f[]' from tree with root
// as 'root'.  depth is current depth and passed as 0 initially.
// Returns root of the modified tree.
Node *deleteNodeRec ( Node *root, int Point2f [ ], int depth ) {
    // Given Point2f is not present
    if ( root == NULL )
        return NULL;

    // Find dimension of current node
    int cd = depth % k;

    // If the Point2f to be deleted is present at root
    if ( arePointsSame ( root->Point2f, Point2f ) ) {
        // 2.b) If right child is not NULL
        if ( root->right != NULL ) {
            // Find minimum of root's dimension in right subtree
            Node *min = findMin ( root->right, cd );

            // Copy the minimum to root
            copyPoint ( root->Point2f, min->Point2f );

            // Recursively delete the minimum
            root->right = deleteNodeRec ( root->right, min->Point2f, depth + 1 );
        }
        else if ( root->left != NULL ) // same as above
        {
            Node *min = findMin ( root->left, cd );
            copyPoint ( root->Point2f, min->Point2f );
            root->right = deleteNodeRec ( root->left, min->Point2f, depth + 1 );
        }
        else // If node to be deleted is leaf node
        {
            delete root;
            return NULL;
        }
        return root;
    }

    // 2) If current node doesn't contain Point2f, search downward
    if ( Point2f [ cd ] < root->Point2f [ cd ] )
        root->left = deleteNodeRec ( root->left, Point2f, depth + 1 );
    else
        root->right = deleteNodeRec ( root->right, Point2f, depth + 1 );
    return root;
}

// Function to delete a given Point2f from K D Tree with 'root'
Node* deleteNode ( Node *root, int Point2f [ ] ) {
  // Pass depth as 0
    return deleteNodeRec ( root, Point2f, 0 );
}

// Driver program to test above functions
int main8798797 ( ) {
    struct Node *root = nullptr;
    int points [ ] [ k ] = { {30, 40}, {5, 25}, {70, 70}, {10, 12}, {50, 30}, {35, 45} };

    constexpr int n = sizeof ( points ) / sizeof ( points [ 0 ] );

    for ( int i = 0; i < n; i++ ) {
        root = insert ( root, points [ i ] );
    }

    // Delete (30, 40);
    root = deleteNode ( root, points [ 0 ] );

    cout << "Root after deletion of (30, 40)\n";
    cout << root->Point2f [ 0 ] << ", " << root->Point2f [ 1 ] << endl;

    return 0;
}



using PointArray = std::array<float, 2>;
using PointToID = spatial::idle_point_multiset<2, PointArray>;

[[ nodiscard ]] PointArray toArray ( const kd::Point2f & v_ ) noexcept {
    return *reinterpret_cast<const PointArray*> ( &v_ );
}

[[ nodiscard ]] kd::Point2f fromArray ( const PointArray & p_ ) noexcept {
    return *reinterpret_cast<const kd::Point2f*> ( &p_ );
}


struct KDTree {

    kdtree *ptree;

    template<typename forward_it>
    KDTree ( forward_it first_, forward_it last_ ) noexcept :
        ptree { kd_create ( 2 ) } {
        std::for_each ( first_, last_, [ this ] ( kd::Point2f & pos ) { kd_insertf ( ptree, ( const float * ) & pos, NULL ); } );
    }

    ~KDTree ( ) noexcept {
        kd_free ( ptree );
    }

    [[ nodiscard ]] kd::Point2f nearest_pnt ( const kd::Point2f & pos_ ) const noexcept {
        kd::Point2f pos;
        struct kdres * res = kd_nearestf ( ptree, ( const float * ) & pos_ );
        kd_res_itemf ( res, ( float * ) & pos );
        kd_res_free ( res );
        return pos;
    }
};

int wmain77897897 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f,  40.0f };

    constexpr int n = 100'000;

    {
        plf::nanotimer timer;
        double st;

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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




int main877989 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };

    std::uniform_real_distribution<float> disy { 0.0f, 550.0f };
    std::uniform_real_distribution<float> disz { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f,  40.0f };

    constexpr int n = 100'000;

    {
        plf::nanotimer timer;
        double st;

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point3f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ), disz ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree3D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point3f ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ), disz ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    {
        plf::nanotimer timer;
        double st;

        std::vector<kd::Point3f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ), disz ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kd::Tree3D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point3f ptf;

        timer.start ( );

        for ( int i = 0; i < 1'000'000; ++i ) {
            ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ), disz ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    return EXIT_SUCCESS;
}


int wmain89879 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100'000.0f };
    std::uniform_real_distribution<float> disx { 0.0f,  40'000.0f };

    constexpr int n = 100'000;

    {
        plf::nanotimer timer;
        double st;

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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

        std::vector<kd::Point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kd::Point2f ptf;

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


template<typename T, typename = std::enable_if_t<std::conjunction_v<std::is_integral<T>, std::is_unsigned<T>>>>
void print_bits ( const T n ) noexcept {
    int c = 0;
    T i = T ( 1 ) << ( sizeof ( T ) * 8 - 1 );
    while ( i ) {
        putchar ( int ( ( n & i ) > 0 ) + int ( 48 ) );
        if ( 0 == c or 8 == c ) {
            putchar ( 32 );
        }
        i >>= 1;
        ++c;
    }
}


union tagged_float {

    explicit tagged_float ( const bool tag_ = false ) noexcept {
        if ( tag_ )
            i = 0b0000'0000'0000'0000'0000'0000'0000'0001;
        else
            value = 0.0f;
    }
    explicit tagged_float ( const float & f_ ) noexcept :
        value { f_ } {
    }
    explicit tagged_float ( float && f_ ) noexcept :
        value { std::move ( f_ ) } {
    }
    explicit tagged_float ( const float & f_, const bool tag_ ) noexcept :
        value { f_ } {
        if ( tag_ )
            i |= 0b0000'0000'0000'0000'0000'0000'0000'0001;
        else
            i &= 0b1111'1111'1111'1111'1111'1111'1111'1110;
    }
    explicit tagged_float ( float && f_, const bool tag_ ) noexcept :
        value { std::move ( f_  ) } {
        if ( tag_ )
            i |= 0b0000'0000'0000'0000'0000'0000'0000'0001;
        else
            i &= 0b1111'1111'1111'1111'1111'1111'1111'1110;
    }

    [[ nodiscard ]] bool is_tagged ( ) const noexcept {
        return i & 0b0000'0000'0000'0000'0000'0000'0000'0001;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const tagged_float & tf_ ) noexcept {
        out_ << tf_.value << ( tf_.is_tagged ( ) ? '#' : ' ' );
        return out_;
    }

    float value;

    private:

    std::uint32_t i;
};

union tagged_double {

    explicit tagged_double ( const bool tag_ = false ) noexcept {
        if ( tag_ )
            i = 0b0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0001;
        else
            value = 0.0;
    }
    explicit tagged_double ( const double & d_ ) noexcept :
        value { d_ } {
    }
    explicit tagged_double ( double && d_ ) noexcept :
        value { std::move ( d_ ) } {
    }
    explicit tagged_double ( const double & d_, const bool tag_ ) noexcept :
        value { d_ } {
        if ( tag_ )
            i |= 0b0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0001;
        else
            i &= 0b1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1110;
    }
    explicit tagged_double ( double && d_, const bool tag_ ) noexcept :
        value { std::move ( d_ ) } {
        if ( tag_ )
            i |= 0b0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0001;
        else
            i &= 0b1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1111'1110;
    }

    [[ nodiscard ]] bool is_tagged ( ) const noexcept {
        return i & 0b0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0001;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const tagged_double & td_ ) noexcept {
        out_ << td_.value << ( td_.is_tagged ( ) ? '#' : ' ' );
        return out_;
    }

    double value;

    private:

    std::uint64_t i;
};


template<typename R>
using tagged_real = typename std::conditional<std::is_same<float, R>::value, tagged_float, tagged_double>::type;


int wmain ( ) {

    tagged_real<float> f1 { 7.1f, true };

    std::cout << f1 << nl;

    tagged_float f2 { 1.9f, false };

    std::cout << f2 << nl;

    return EXIT_SUCCESS;
}




#include "kdtree2.hpp"

#include <boost/multi_array.hpp>

splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };
std::uniform_real_distribution<float> disy { 0.0f, 100'000.0f };
std::uniform_real_distribution<float> disx { 0.0f,  40'000.0f };

//
// define, for convenience a 2d array of floats.
//
typedef boost::multi_array<float, 2> array2dfloat;


double time_a_search ( kdtree2::KDTree* tree, int nsearch ) {

    const int dim = tree->dim;
    std::vector<float> query ( dim );
    kdtree2::KDTreeResultVector result;

    plf::nanotimer timer;

    timer.start ( );

    for ( int i = 0; i < nsearch; i++ ) {
        query [ 0 ] = disx ( rng ); query [ 1 ] = disy ( rng );
        tree->n_nearest ( query, 1, result );
    }

    return timer.get_elapsed_us ( ) / 1'000'000.0;
}

double time_a_search ( kd::Tree2D<float> & tree, int nsearch ) {

    plf::nanotimer timer;

    kd::Point2f ptf;

    timer.start ( );

    for ( int i = 0; i < nsearch; i++ ) {
        ptf += tree.nearest_pnt ( { disx ( rng ), disy ( rng ) } );
    }

    return timer.get_elapsed_us ( ) / 1'000'000.0;
}

void time_random_searches ( kdtree2::KDTree* tree ) {

  // emit the number of searches per second.

    int nsearch;

    nsearch = 50;

    while ( true ) {
        double t = time_a_search ( tree, nsearch );
        if ( t < 1.0 ) {
            nsearch *= 5;
            continue;
        }
        else {
            std::uint64_t sps = double ( nsearch ) / t;
            std::cout << "C++ impl, for nn=" << '1' << " searches/sec = " << sps << "\n";
            return;
        }
    }
}

void time_random_searches ( kd::Tree2D<float> & tree ) {

  // emit the number of searches per second.

    int nsearch;

    nsearch = 50;

    while ( true ) {
        double t = time_a_search ( tree, nsearch );
        if ( t < 1.0 ) {
            nsearch *= 5;
            continue;
        }
        else {
            std::uint64_t sps = double ( nsearch ) / t;
            std::cout << "C++ impl, for nn=" << '1' << " searches/sec = " << sps << "\n";
            return;
        }
    }
}

int main684984 ( ) {

    int N = 100'000'000, dim = 2;

    {
        plf::nanotimer timer;

        array2dfloat realdata;
        realdata.resize ( boost::extents [ N ] [ dim ] );
        for ( int i = 0; i < N; i++ ) {
            realdata [ i ] [ 0 ] = disx ( rng );
            realdata [ i ] [ 1 ] = disy ( rng );
        }

        kdtree2::KDTree* tree;
        kdtree2::KDTreeResultVector res;

        timer.start ( );
        tree = new kdtree2::KDTree ( realdata, true );
        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        time_random_searches ( tree );

        std::cout << nl;

        free ( tree );
    }

    {
        plf::nanotimer timer;

        std::vector<kd::Point2f> points;
        for ( int i = 0; i < N; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );
        kd::Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );
        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        time_random_searches ( tree );

        std::cout << nl;
    }

    return EXIT_SUCCESS;
}


#include <queue>

template<typename Point>
struct PQType {

    using base_type = decltype ( Point { }.x );

    base_type value = std::numeric_limits<base_type>::max ( );
    Point point;

    PQType ( ) noexcept = default;
    PQType ( const PQType & ) noexcept = default;
    PQType ( PQType && ) noexcept = default;
    PQType ( const base_type & value_, const Point & p_ ) noexcept :
        value { value_ }, point { p_ } {
    }
    PQType ( base_type && value_, Point && p_ ) noexcept :
        value { std::move ( value_ ) }, point { std::move ( p_ ) } {
    }
    PQType ( base_type && value_, base_type && x_, base_type && y_ ) noexcept :
        value { std::move ( value_ ) }, point { std::move ( x_ ), std::move ( y_ ) } {
    }

    [[ maybe_unused ]] PQType & operator = ( const PQType & ) noexcept = default;
    [[ maybe_unused ]] PQType & operator = ( PQType && ) noexcept = default;

    [[ nodiscard ]] bool operator < ( const PQType & p_ ) const noexcept {
        return value < p_.value;
    }
    [[ nodiscard ]] bool operator > ( const PQType & p_ ) const noexcept {
        return value > p_.value;
    }
};

#include "sorted_vector_set.hpp"


template<typename Point>
using PQueue = sorted_vector_set<PQType<Point>>;


template<typename Point>
struct KNearest : public sorted_vector_set<PQType<Point>> {

    using base = sorted_vector_set<PQType<Point>>;

    using typename base::reference;
    using typename base::rv_reference;
    using typename base::const_reference;
    using typename base::pointer;
    using typename base::const_pointer;
    using typename base::size_type;
    using typename base::value_type;
    using typename base::container;
    using iterator = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    using base_type = typename value_type::base_type;

    KNearest ( const std::size_t s_ ) : base::sorted_vector_set ( s_ ) { }

    template<typename ... Args>
    void emplace ( base_type && value_, Args && ... args_ ) noexcept {
        if ( base::top ( ).value > value_ ) {
            base::pop ( );
            base::emplace ( std::move ( value_ ), std::forward<Args> ( args_ ) ... );
        }
    }

    [[ nodiscard ]] const_reference bottom ( ) const noexcept { return base::bottom ( ); }
    [[ nodiscard ]] const_reference top ( ) const noexcept { return base::top ( ); }

    [[ nodiscard ]] base_type value ( ) const noexcept { return base::top ( ).value; }
};


int wmain879808 ( ) {

    std::uniform_real_distribution<float> pdisv { 0.0f,  100.0f };
    std::uniform_real_distribution<float> pdisy { 0.0f, 100'000.0f };
    std::uniform_real_distribution<float> pdisx { 0.0f,  40'000.0f };

    KNearest<kd::Point2f> knn ( 50u );

    for ( int i = 0; i < 500; ++i ) {
        knn.emplace ( pdisv ( rng ), pdisx ( rng ), pdisy ( rng ) );
    }

    std::cout << knn.top ( ).point << ' ' << knn.top ( ).value << nl;
    std::cout << knn.bottom ( ).point << ' ' << knn.bottom ( ).value << nl;
    std::cout << knn.value ( ) << nl;

    return EXIT_SUCCESS;
}


int wmain64631 ( ) {

    sorted_vector_multiset<int> ms { 1, 2, 3, 4, 5 };

    for ( auto & v : ms ) {
        std::cout << v << ' ';
    }
    std::cout << nl;

    std::cout << *ms.lower_bound ( 3 ) << nl;

    ms.insert ( 3 );

    for ( auto & v : ms ) {
        std::cout << v << ' ';
    }
    std::cout << nl;

    ms.insert ( 0 );

    for ( auto & v : ms ) {
        std::cout << v << ' ';
    }
    std::cout << nl;

    return EXIT_SUCCESS;
}
