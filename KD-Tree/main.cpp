
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


#include <spatial/idle_point_multimap.hpp>
#include <spatial/idle_point_multiset.hpp>
#include <spatial/neighbor_iterator.hpp>


#include <SFML/System.hpp>
#include <splitmix.hpp>
#include <plf/plf_nanotimer.h>

#include "ikdtree.hpp"
#include "kdtree.h"

#define nl '\n'

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)


[[ nodiscard ]] constexpr auto nn_distance_squared ( const kdt::point2f & p1_, const kdt::point2f & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
}
[[ nodiscard ]] constexpr auto nn_distance_squared ( const kdt::point3f & p1_, const kdt::point3f & p2_ ) noexcept {
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

    std::vector<kdt::point2f> points;

    for ( int i = 0; i < n_; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

    bool rv = true;

    for ( int i = 0; i < 1'000'000; ++i ) {
        const kdt::point2f ptf { disx ( rng ), disy ( rng ) };
        rv = rv and ( tree.nearest_pnt ( ptf ) == nn_search_linear ( std::begin ( points ), std::end ( points ), ptf ) );
        if ( not ( rv ) ) {
            std::cout << "fail\n";
            exit ( 0 );
        }
    }

    return rv;
}

int main8978978 ( ) {

    std::cout << std::boolalpha << test ( 1'000'000 ) << nl;

    return EXIT_SUCCESS;
}



int wmain8797 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f, 40.0f };

    for ( int i = 0; i < 100; ++i ) {

        plf::nanotimer timer;
        double st;

        constexpr int n = 1'000;

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        // std::cout << nl << tree << nl << nl;

        //kdt::point2f p { disx ( rng ), disy ( rng ) };

        //std::cout << "ptf " << p << nl;

        bool result = true;

        constexpr int cnt = 100'000;

        timer.start ( );
        for ( int i = 0; i < cnt; ++i ) {
            /*
            const kdt::point2f p { disx ( rng ), disy ( rng ) };
            bool r = tree.nearest_pnt ( p ) == kdt::i2dtree<float>::nearest_linear_pnt ( p, points );
            if ( not ( r ) ) {
                const kdt::point2f p1 = tree.nearest_pnt ( p ), p2 = kdt::i2dtree<float>::nearest_linear_pnt ( p, points );
                const float f1 = kdt::i2dtree<float>::distance_squared ( p, p1 ), f2 = kdt::i2dtree<float>::distance_squared ( p, p2 );
                if ( f1 == f2 ) {
                    continue;
                }
                std::cout << p1 << f1 << p2 << f2 << nl;
            }
            result &= r;
            */
        }
        // std::cout << "elapsed im " << ( std::uint64_t ) ( timer.get_elapsed_ns ( ) / cnt ) << " ns" << nl;

        std::cout << std::boolalpha << result << nl;

        // std::cout << "nearest im " << found_impl << " " << kdt::i2dtree<kdt::point2f>::nearest_linear_pnt ( p, points ) << nl;
    }

    return EXIT_SUCCESS;
}



int wmain ( ) {

    // std::vector<kdt::point2f> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
    std::vector<kdt::point2f> points { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    kdt::i2dtree<float> tree { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

    std::cout << nl << tree << nl << nl;

    kdt::point2f ptf { 7.6f, 7.9f };

    std::cout << nl << nl << "nearest " << nl << tree.nearest_pnt ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << kdt::i2dtree<float>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}


using namespace std;

const int k = 2;

// A structure to represent node of kd tree
struct Node {
    int point2f [ k ]; // To store k dimensional point2f
    Node *left, *right;
};

// A method to create a node of K D tree
struct Node* newNode ( int arr [ ] ) {
    struct Node* temp = new Node;

    for ( int i = 0; i < k; i++ )
        temp->point2f [ i ] = arr [ i ];

    temp->left = temp->right = NULL;
    return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec ( Node *root, int point2f [ ], unsigned depth ) {
    // Tree is empty?
    if ( root == NULL )
        return newNode ( point2f );

    // Calculate current dimension (cd) of comparison
    unsigned cd = depth % k;

    // Compare the new kdt::point2f with root on current dimension 'cd'
    // and decide the left or right subtree
    if ( point2f [ cd ] < ( root->point2f [ cd ] ) )
        root->left = insertRec ( root->left, point2f, depth + 1 );
    else
        root->right = insertRec ( root->right, point2f, depth + 1 );

    return root;
}

// Function to insert a new kdt::point2f with given point2f in
// KD Tree and return new root. It mainly uses above recursive
// function "insertRec()"
Node* insert ( Node *root, int point2f [ ] ) {
    return insertRec ( root, point2f, 0 );
}

// A utility function to find minimum of three integers
Node *minNode ( Node *x, Node *y, Node *z, int d ) {
    Node *res = x;
    if ( y != NULL && y->point2f [ d ] < res->point2f [ d ] )
        res = y;
    if ( z != NULL && z->point2f [ d ] < res->point2f [ d ] )
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

    // Compare point2f with root with respect to cd (Current dimension)
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

// Copies point2f p2 to p1
void copyPoint ( int p1 [ ], int p2 [ ] ) {
    for ( int i = 0; i < k; i++ )
        p1 [ i ] = p2 [ i ];
}

// Function to delete a given point2f 'point2f[]' from tree with root
// as 'root'.  depth is current depth and passed as 0 initially.
// Returns root of the modified tree.
Node *deleteNodeRec ( Node *root, int point2f [ ], int depth ) {
    // Given point2f is not present
    if ( root == NULL )
        return NULL;

    // Find dimension of current node
    int cd = depth % k;

    // If the point2f to be deleted is present at root
    if ( arePointsSame ( root->point2f, point2f ) ) {
        // 2.b) If right child is not NULL
        if ( root->right != NULL ) {
            // Find minimum of root's dimension in right subtree
            Node *min = findMin ( root->right, cd );

            // Copy the minimum to root
            copyPoint ( root->point2f, min->point2f );

            // Recursively delete the minimum
            root->right = deleteNodeRec ( root->right, min->point2f, depth + 1 );
        }
        else if ( root->left != NULL ) // same as above
        {
            Node *min = findMin ( root->left, cd );
            copyPoint ( root->point2f, min->point2f );
            root->right = deleteNodeRec ( root->left, min->point2f, depth + 1 );
        }
        else // If node to be deleted is leaf node
        {
            delete root;
            return NULL;
        }
        return root;
    }

    // 2) If current node doesn't contain point2f, search downward
    if ( point2f [ cd ] < root->point2f [ cd ] )
        root->left = deleteNodeRec ( root->left, point2f, depth + 1 );
    else
        root->right = deleteNodeRec ( root->right, point2f, depth + 1 );
    return root;
}

// Function to delete a given point2f from K D Tree with 'root'
Node* deleteNode ( Node *root, int point2f [ ] ) {
  // Pass depth as 0
    return deleteNodeRec ( root, point2f, 0 );
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
    cout << root->point2f [ 0 ] << ", " << root->point2f [ 1 ] << endl;

    return 0;
}



using PointArray = std::array<float, 2>;
using PointToID = spatial::idle_point_multiset<2, PointArray>;

[[ nodiscard ]] PointArray toArray ( const kdt::point2f & v_ ) noexcept {
    return *reinterpret_cast<const PointArray*> ( &v_ );
}

[[ nodiscard ]] kdt::point2f fromArray ( const PointArray & p_ ) noexcept {
    return *reinterpret_cast<const kdt::point2f*> ( &p_ );
}


struct KDTree {

    kdtree *ptree;

    template<typename forward_it>
    KDTree ( forward_it first_, forward_it last_ ) noexcept :
        ptree { kd_create ( 2 ) } {
        std::for_each ( first_, last_, [ this ] ( kdt::point2f & pos ) { kd_insertf ( ptree, ( const float * ) & pos, NULL ); } );
    }

    ~KDTree ( ) noexcept {
        kd_free ( ptree );
    }

    [[ nodiscard ]] kdt::point2f nearest_pnt ( const kdt::point2f & pos_ ) const noexcept {
        kdt::point2f pos;
        struct kdres * res = kd_nearestf ( ptree, ( const float * ) & pos_ );
        kd_res_itemf ( res, ( float * ) & pos );
        kd_res_free ( res );
        return pos;
    }
};

int main676786 ( ) {

    splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 0.0f,  40.0f };

    constexpr int n = 43;

    {
        plf::nanotimer timer;
        double st;

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        KDTree tree ( std::begin ( points ), std::end ( points ) );
        // kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point3f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ), disz ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i3dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point3f ptf;

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

        std::vector<kdt::point3f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ), disz ( rng ) );
        }

        timer.start ( );

        // KDTree tree ( std::begin ( points ), std::end ( points ) );
        kdt::i3dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point3f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

        std::vector<kdt::point2f> points;

        for ( int i = 0; i < n; ++i ) {
            points.emplace_back ( disx ( rng ), disy ( rng ) );
        }

        timer.start ( );

        kdt::i2dtree<float> tree ( std::begin ( points ), std::end ( points ) );

        std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        kdt::point2f ptf;

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

constexpr std::uint32_t bit1 = 0b0000'0000'0000'0000'0000'0000'0000'0001;
constexpr std::uint32_t bits31 = ~bit1;

union float_type {

    float_type ( const float & num_ ) noexcept :
        f ( num_ ) {
    }
    float_type ( float && num_ = 0.0f ) noexcept :
        f ( std::move ( num_ ) ) {
    }
    float_type ( std::uint32_t && num_ = 0u ) noexcept :
        i ( std::move ( num_ ) ) {
    }

    // Portable extraction of components.

    bool is_negative ( ) const noexcept { return i & ( 1u << 31 ); }
    std::uint32_t exponent ( ) const noexcept { return ( i >> 23 ) & 0xFF; }
    std::uint32_t mantissa ( ) const noexcept { return i & ( ( 1 << 23 ) - 1 ); }

    void tag ( ) noexcept {
        i |= bit1;
    }

    void clear_low_bit ( ) noexcept {
        i &= bits31;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const float_type & f_ ) noexcept {
        print_bits ( f_.i );
        std::cout << ' ' << f_.f;
        return out_;
    }

    std::uint32_t i;
    float f;
};



int main789789 ( ) {

    float_type f { -22661.6516651175f };

    std::cout << f.exponent ( ) << ' ' << f.mantissa ( ) << nl;

    std::cout << f << nl;

    f.tag ( );

    std::cout << f.exponent ( ) << ' ' << f.mantissa ( ) << nl;

    std::cout << f << nl;

    float_type f2 { 0b1100'0001'1011'0110'0000'0000'0000'0001 };

    std::cout << f2.exponent ( ) << ' ' << f2.mantissa ( ) << nl;

    std::cout << f2 << nl;

    return EXIT_SUCCESS;
}
