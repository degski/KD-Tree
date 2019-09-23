
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
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
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
#include <iterator>
#include <list>
#include <map>
#include <random>
#include <sax/iostream.hpp>
#include <stack>
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

#include <spatial/idle_point_multiset.hpp>
#include <spatial/neighbor_iterator.hpp>

#include <sax/prng_jsf.hpp>
#include <sax/splitmix.hpp>
#include <sax/uniform_int_distribution.hpp>

#include <plf/plf_nanotimer.h>

#include "ikdtree.hpp"
#include "kdtree.h"
#include "kdtree2.hpp"

#define likely( x ) __builtin_expect ( !!( x ), 1 )
#define unlikely( x ) __builtin_expect ( !!( x ), 0 )

int main6576756 ( ) {

    // std::vector<kd::Point2f> points { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, {
    // 8, 1 }, { 7, 2 } };
    std::vector<kd::Point2f> points{ { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 },
                                     { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 },  { 8, 5 }, { 9, 9 } };

    for ( auto p : points ) {
        std::cout << p;
    }
    std::cout << nl;

    kd::Tree2D<float> tree{ { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 },
                            { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 },  { 8, 5 }, { 9, 9 } };

    std::cout << nl << tree << nl << nl;

    kd::Point2f ptf{ 7.6f, 7.9f };

    std::cout << nl << nl << "nearest " << nl << *tree.nn_pointer ( ptf ) << nl;

    std::cout << nl;

    for ( auto p : points ) {
        std::cout << kd::Tree2D<float>::distance_squared ( p, ptf ) << ' ' << p << nl;
    }

    return EXIT_SUCCESS;
}

void handleEptr ( std::exception_ptr eptr ) { // Passing by value is ok.
    try {
        if ( eptr )
            std::rethrow_exception ( eptr );
    }
    catch ( const std::exception & e ) {
        std::cout << "Caught exception \"" << e.what ( ) << "\"\n";
    }
}

template<typename T>
struct Point2 {

    using value_type = T;

    value_type x, y;

    Point2 ( ) noexcept                = default;
    Point2 ( Point2 const & ) noexcept = default;
    Point2 ( Point2 && ) noexcept      = default;
    Point2 ( value_type && x_, value_type && y_ ) noexcept : x{ std::move ( x_ ) }, y{ std::move ( y_ ) } {}

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
        if ( Point2{ std::numeric_limits<value_type>::max ( ), std::numeric_limits<value_type>::max ( ) } != p_ )
            out_ << '<' << p_.x << ' ' << p_.y << '>';
        else
            out_ << "<* *>";
        return out_;
    }
};

[[nodiscard]] constexpr auto nn_distance_squared ( Point2<float> const & p1_, Point2<float> const & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
}
[[nodiscard]] constexpr auto nn_distance_squared ( kd::Point2f const & p1_, kd::Point2f const & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
}
[[nodiscard]] constexpr auto nn_distance_squared ( kd::Point3f const & p1_, kd::Point3f const & p2_ ) noexcept {
    return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) ) +
           ( ( p1_.z - p2_.z ) * ( p1_.z - p2_.z ) );
}

template<typename forward_it, typename value_type>
[[nodiscard]] value_type * nn_search_linear ( forward_it first_, forward_it last_, value_type const & p_ ) noexcept {
    using base_type        = decltype ( p_.x );
    base_type min_distance = std::numeric_limits<base_type>::max ( );
    forward_it found       = last_;
    while ( first_ != last_ ) {
        auto const d = nn_distance_squared ( p_, *first_ );
        if ( d < min_distance ) {
            min_distance = d;
            found        = first_;
        }
        ++first_;
    }
    return &*found;
}

using fran = std::uniform_real_distribution<float>;

bool test ( int const n_ ) {
    sax::Rng rng{ std::uint64_t ( n_ + 1 ) };
    fran disy{ fran ( 0.0f, 4'999.0f ) ( rng ), fran ( 5'000.0f, 9'999.0f ) ( rng ) };
    fran disx{ fran ( 0.0f, 4'999.0f ) ( rng ), fran ( 5'000.0f, 9'999.0f ) ( rng ) };
    std::vector<kd::Point2f> points;
    points.reserve ( n_ );
    for ( int i = 0; i < n_; ++i )
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    kd::ikdtree<float, 2> tree ( std::begin ( points ), std::end ( points ) );
    for ( int i = 0; i < 1'000; ++i ) {
        kd::Point2f const ptf{ disx ( rng ), disy ( rng ) }, kdt{ *tree.nn_pointer ( ptf ) },
            ls{ *nn_search_linear ( std::begin ( points ), std::end ( points ), ptf ) };
        if ( nn_distance_squared ( ptf, kdt ) != nn_distance_squared ( ptf, ls ) ) {
            std::cout << "fail\n";
            std::cout << "n = " << n_ << " point " << ptf << nl;
            std::cout << "values " << kdt << " vs " << ls << nl;
            std::cout << "dist " << nn_distance_squared ( ptf, kdt ) << " vs " << nn_distance_squared ( ptf, ls ) << nl;
            std::cout << std::boolalpha << ( nn_distance_squared ( ptf, kdt ) == nn_distance_squared ( ptf, ls ) ) << nl;
            return false;
        }
    }
    return true;
}

int main ( ) {
    std::exception_ptr eptr;
    sax::Rng rng{ sax::fixed_seed ( ) };
    sax::uniform_int_distribution<int> dis{ 16, 10'000 };
    try {
        for ( int i = 0; i < 10'000; ++i )
            std::cout << i << ' ' << std::boolalpha << test ( dis ( rng ) ) << nl;
    }
    catch ( ... ) {
        eptr = std::current_exception ( ); // Capture.
    }
    handleEptr ( eptr );
    return EXIT_SUCCESS;
}

#if 0


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

[[ nodiscard ]] PointArray toArray ( kd::Point2f const & v_ ) noexcept {
    return *reinterpret_cast<const PointArray*> ( &v_ );
}

[[ nodiscard ]] kd::Point2f fromArray ( const PointArray & p_ ) noexcept {
    return *reinterpret_cast<kd::Point2f const*> ( &p_ );
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

    [[ nodiscard ]] kd::Point2f nn_pointer ( kd::Point2f const & pos_ ) const noexcept {
        kd::Point2f pos;
        struct kdres * res = kd_nearestf ( ptree, ( const float * ) & pos_ );
        kd_res_itemf ( res, ( float * ) & pos );
        kd_res_free ( res );
        return pos;
    }
};

int main77897 ( ) {

    sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    return EXIT_SUCCESS;
}




int main877989 ( ) {

    sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };

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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ), disz ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ), disz ( rng ) } );
        }

        std::cout << "elapsed search " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

        std::cout << nl << nl << "nearest " << ptf << nl;

        std::cout << nl;
    }

    return EXIT_SUCCESS;
}


int wmain89879 ( ) {

    sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) ); } ( ) };
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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
            ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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


template<typename T = float>
union tagged {

    using value_type = T;
    using uintv_type = typename std::conditional<std::is_same<T, float>::value, std::uint32_t, std::uint64_t>::type;

    static constexpr uintv_type one { 1u };
    static constexpr uintv_type not_one { ~one };

    explicit tagged ( const bool tag_ = false ) noexcept {
        if ( tag_ )
            i = one;
        else
            value = value_type { 0 };
    }
    explicit tagged ( const tagged & tr_ ) noexcept :
        value { tr_.value } { }
    explicit tagged ( tagged && tr_ ) noexcept :
        value { std::move ( tr_.value ) } { }
    explicit tagged ( const tagged & tr_, const bool tag_ ) noexcept :
        value { tr_.value } {
        if ( tag_ )
            i |= one;
        else
            i &= not_one;
    }
    explicit tagged ( tagged && tr_, const bool tag_ ) noexcept :
        value { std::move ( tr_.value ) } {
        if ( tag_ )
            i |= one;
        else
            i &= not_one;
    }
    explicit tagged ( value_type const & d_ ) noexcept :
        value { d_ } {
        i &= not_one;
    }
    explicit tagged ( value_type && d_ ) noexcept :
        value { std::move ( d_ ) } {
        i &= not_one;
    }
    explicit tagged ( value_type const & d_, const bool tag_ ) noexcept :
        value { d_ } {
        if ( tag_ )
            i |= one;
        else
            i &= not_one;
    }
    explicit tagged ( value_type && d_, const bool tag_ ) noexcept :
        value { std::move ( d_ ) } {
        if ( tag_ )
            i |= one;
        else
            i &= not_one;
    }

    [[ nodiscard ]] bool is_tagged ( ) const noexcept {
        return i & one;
    }

    void tag ( ) noexcept {
        i |= one;
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const tagged & td_ ) noexcept {
        out_ << td_.value << ( td_.is_tagged ( ) ? '#' : ' ' );
        return out_;
    }

    [[ nodiscard ]] bool operator == ( const tagged & t_ ) const noexcept {
        return ( i >> 1 ) == ( t_.i >> 1 );
    }
    [[ nodiscard ]] bool operator == ( value_type const & v_ ) const noexcept {
        return *this == static_cast<tagged> ( v_ );
    }
    [[ nodiscard ]] bool operator != ( const tagged & t_ ) const noexcept {
        return ( i >> 1 ) != ( t_.i >> 1 );
    }
    [[ nodiscard ]] bool operator != ( value_type const & v_ ) const noexcept {
        return *this != static_cast< tagged > ( v_ );
    }

    [[ nodiscard ]] bool operator < ( const tagged & t_ ) const noexcept {
        return ( i & not_one ) < ( t_.i & not_one );
    }
    [[ nodiscard ]] bool operator < ( value_type const & v_ ) const noexcept {
        return *this < static_cast<tagged> ( v_ );
    }

    [[ nodiscard ]] value_type operator - ( const tagged & t_ ) const noexcept {
        return static_cast<value_type> ( i & not_one ) - static_cast<value_type> ( t_.i & not_one );
    }
    [[ nodiscard ]] value_type operator + ( const tagged & t_ ) const noexcept {
        return static_cast<value_type> ( i & not_one ) + static_cast<value_type> ( t_.i & not_one );
    }
    /*
    [[ nodiscard ]] value_type operator - ( value_type const & v_ ) const noexcept {
        return ( static_cast< tagged > ( static_cast<value_type> ( i & not_one ) - static_cast<value_type> ( static_cast<tagged> ( v_ ).i & not_one ) & not_one ).value;
    }
    [[ nodiscard ]] value_type operator + ( value_type const & v_ ) const noexcept {
        return ( static_cast<tagged> ( static_cast<value_type> ( i & not_one ) + static_cast<value_type> ( static_cast<tagged> ( v_ ).i & not_one ) ) & not_one ).value;
    }
    */

    value_type value;

    private:

    uintv_type i;
};




template<typename T = float>
struct TaggedPoint2 {

    using value_type = T;

    tagged<value_type> x; value_type y;

    TaggedPoint2 ( ) noexcept = default;
    TaggedPoint2 ( const TaggedPoint2 & ) noexcept = default;
    TaggedPoint2 ( TaggedPoint2 && ) noexcept = default;
    TaggedPoint2 ( value_type && x_, value_type && y_ ) noexcept :
        x { std::move ( x_ ) }, y { std::move ( y_ ) } {
    }
    TaggedPoint2 ( tagged<value_type> && x_, value_type && y_ ) noexcept :
        x { std::move ( x_ ) }, y { std::move ( y_ ) } {
    }

    [[ maybe_unused ]] TaggedPoint2 & operator = ( const TaggedPoint2 & p_ ) noexcept {
        x.value = p_.x.value; y = p_.y;
        return * this;
    }
    [[ maybe_unused ]] TaggedPoint2 & operator = ( TaggedPoint2 && p_ ) noexcept {
        x.value = std::move ( p_.x.value ); y = std::move ( p_.y );
        return * this;
    }
    [[ maybe_unused ]] TaggedPoint2 & operator = ( Point2 const<T> & p_ ) noexcept {
        x.value = p_.x; y = p_.y;
        return *this;
    }
    [[ maybe_unused ]] TaggedPoint2 & operator = ( Point2<T> && p_ ) noexcept {
        x.value = std::move ( p_.x ); y = std::move ( p_.y );
        return *this;
    }

    [[ nodiscard ]] bool operator == ( const TaggedPoint2 & p_ ) const noexcept {
        return x.value == p_.x.value and y == p_.y;
    }
    [[ nodiscard ]] bool operator != ( const TaggedPoint2 & p_ ) const noexcept {
        return x.value != p_.x.value or y != p_.y;
    }

    void tag ( ) noexcept {
        x.tag ( );
    }

    template<typename Stream>
    [[ maybe_unused ]] friend Stream & operator << ( Stream & out_, const TaggedPoint2 & p_ ) noexcept {
        if ( TaggedPoint2 { std::numeric_limits<value_type>::max ( ), std::numeric_limits<value_type>::max ( ) } != p_ )
            out_ << '<' << p_.x << ' ' << p_.y << '>';
        else
            out_ << "<* *>";
        return out_;
    }
};

template<typename T = float, typename P = Point2<T>, typename TP = TaggedPoint2<T>>
struct Tree2D {

    using value_type = P;
    using tagged_value_type = TP;
    using base_type = T;
    using pointer = tagged_value_type * ;
    using reference = tagged_value_type & ;
    using const_pointer = tagged_value_type const *;
    using const_reference = tagged_value_type const &;

    using container = std::vector<value_type>;
    using tagged_container = std::vector<tagged_value_type>;

    using iterator = typename tagged_container::iterator;
    using const_iterator = typename tagged_container::const_iterator;

    private:

    struct nn_data {
        value_type point;
        const_pointer found;
        base_type min_distance;
    };

    template<typename forward_it>
    [ [ nodiscard ] ] std::size_t get_dimensions_order ( forward_it first_, forward_it last_ ) const noexcept {
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
        return p_->x.is_tagged ( );
    }

    template<typename random_it>
    void kd_construct_xy ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
        else {
            p_->tag ( );
        }
    }
    template<typename random_it>
    void kd_construct_yx ( const pointer p_, random_it first_, random_it last_ ) noexcept {
        random_it median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [ ] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
        else {
            p_->tag ( );
        }
    }

    void nn_search_xy ( const const_pointer p_ ) const noexcept {
        base_type d = Tree2D::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_nearest.point.x ) > base_type { 0 } ) {
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
    void nn_search_yx ( const const_pointer p_ ) const noexcept {
        base_type d = Tree2D::distance_squared ( *p_, m_nearest.point );
        if ( d < m_nearest.min_distance ) {
            m_nearest.min_distance = d; m_nearest.found = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_nearest.point.y ) > base_type { 0 } ) {
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

    void nn_search_linear ( ) const noexcept {
        /*
        for ( auto && v : m_data ) {
            auto const d = distance_squared ( m_nearest.point, v );
            if ( d < m_nearest.min_distance ) {
                m_nearest.found = &v;
                m_nearest.min_distance = d;
            }
        }
        */
    }

    tagged_container m_data;
    std::size_t m_dim;
    mutable nn_data m_nearest;

    static constexpr std::size_t m_linear_bound = 4u;

    public:

    Tree2D ( const Tree2D & ) = delete;
    Tree2D ( Tree2D && ) noexcept = delete;

    Tree2D ( std::initializer_list<value_type> il_ ) noexcept {
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( il_.size ( ) ), tagged_value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
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
                //std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    template<typename forward_it>
    Tree2D ( forward_it first_, forward_it last_ ) noexcept {
        if ( first_ < last_ ) {
            std::size_t const n = std::distance ( first_, last_ );
            if ( n > m_linear_bound ) {
                m_data.resize ( bin_tree_size<std::size_t> ( static_cast< std::size_t > ( n ) ), tagged_value_type { std::numeric_limits<base_type>::max ( ), std::numeric_limits<base_type>::max ( ) } );
                m_dim = get_dimensions_order ( first_, last_ );
                switch ( m_dim ) {
                case 0: kd_construct_xy ( m_data.data ( ), first_, last_ ); break;
                case 1: kd_construct_yx ( m_data.data ( ), first_, last_ ); break;
                }
            }
            else {
                m_data.reserve ( n );
                //std::copy ( first_, last_, std::back_inserter ( m_data ) );
                m_dim = 2u;
            }
        }
    }

    Tree2D & operator = ( const Tree2D & ) = delete;
    Tree2D & operator = ( Tree2D && ) noexcept = delete;

    [[ nodiscard ]] const_pointer nn_ptr ( value_type const & point_ ) const noexcept {
        m_nearest = { point_, nullptr, std::numeric_limits<base_type>::max ( ) };
        switch ( m_dim ) {
        case 0: nn_search_xy ( m_data.data ( ) ); break;
        case 1: nn_search_yx ( m_data.data ( ) ); break;
        case 2: nn_search_linear ( ); break;
        }
        return m_nearest.found;
    }

    [[ nodiscard ]] std::ptrdiff_t nn_idx ( value_type const & point_ ) const noexcept {
        return nn_ptr ( point_ ) - m_data.data ( );
    }

    [[ nodiscard ]] value_type nn_pointer ( value_type const & point_ ) const noexcept {
        auto [ x, y ] { *nn_ptr ( point_ ) };
        return { std::move ( x.value ), std::move ( y ) };
    }

    template<typename T1, typename T2>
    [[ nodiscard ]] static constexpr base_type distance_squared ( const T1 & p1_, const T2 & p2_ ) noexcept {
        return ( ( p1_.x - p2_.x ) * ( p1_.x - p2_.x ) ) + ( ( p1_.y - p2_.y ) * ( p1_.y - p2_.y ) );
    }

    template<typename Stream>
    [ [ maybe_unused ] ] friend Stream & operator << ( Stream & out_, const Tree2D & tree_ ) noexcept {
        for ( const auto & p : tree_.m_data ) {
            out_ << p;
        }
        return out_;
    }

    private:

    template<typename U>
    [ [ nodiscard ] ] static constexpr U bin_tree_size ( const U i_ ) noexcept {
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

bool test2 ( int const n_ ) noexcept {

    sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };

    std::uniform_real_distribution<float> disy { 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx { 20.0f, 40.0f };

    std::vector<Point2<float>> points;
    points.reserve ( n_ );

    for ( int i = 0; i < n_; ++i ) {
        points.emplace_back ( disx ( rng ), disy ( rng ) );
    }

    Tree2D<float> tree ( std::begin ( points ), std::end ( points ) );

    bool rv = true;

    for ( int i = 0; i < 100'000; ++i ) {
        Point2<float> const ptf { disx ( rng ), disy ( rng ) };
        auto p1 = tree.nn_pointer ( ptf );
        auto p2 = nn_search_linear ( std::begin ( points ), std::end ( points ), ptf );
        std::cout << p1 << p2 << nl;
        rv = rv and ( p1 == p2 );
        if ( not ( rv ) ) {
            std::cout << "fail\n";
            exit ( 0 );
        }
    }

    return rv;
}

int wmain65152 ( ) {

    sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };
    sax::uniform_int_distribution<int> dis { 6, 20 };

    for ( int i = 0; i < 1'000; ++i ) {
        std::cout << std::boolalpha << test2 ( dis ( rng ) ) << nl;
    }

    return EXIT_SUCCESS;
}


int wmain ( ) {

    {

        Tree2D<float> tree { { 1, 3 }, { 1, 8 }, { 2, 2 }, { 2, 10 }, { 3, 6 }, { 4, 1 }, { 5, 4 }, { 6, 8 }, { 7, 4 }, { 7, 7 }, { 8, 2 }, { 8, 5 }, { 9, 9 } };

        std::cout << nl << tree << nl << nl;

        Point2 ptf { 7.6f, 7.9f };

        std::cout << nl << nl << "nearest " << nl << tree.nn_pointer ( ptf ) << nl;
    }


    return EXIT_SUCCESS;
}

#    include "kdtree2.hpp"

#    include <boost/multi_array.hpp>

sax::splitmix64 rng { [ ] ( ) { std::random_device rdev; return ( static_cast< std::size_t > ( rdev ( ) ) << 32 ) | static_cast< std::size_t > ( rdev ( ) ); } ( ) };
std::uniform_real_distribution<float> disy { 0.0f, 100'000.0f };
std::uniform_real_distribution<float> disx { 0.0f,  40'000.0f };

//
// define, for convenience a 2d array of floats.
//
typedef boost::multi_array<float, 2> array2dfloat;


double time_a_search ( kdtree2::KDTree* tree, int nsearch ) {

    int const dim = tree->dim;
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
        ptf += tree.nn_pointer ( { disx ( rng ), disy ( rng ) } );
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

#    include <queue>

template<typename Point>
struct PQType {

    using base_type = decltype ( Point { }.x );

    base_type value = std::numeric_limits<base_type>::max ( );
    Point point;

    PQType ( ) noexcept = default;
    PQType ( const PQType & ) noexcept = default;
    PQType ( PQType && ) noexcept = default;
    PQType ( auto const & value_, const Point & p_ ) noexcept :
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

#    include "sorted_vector_set.hpp"


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

    KNearest ( std::size_t const s_ ) : base::sorted_vector_set ( s_ ) { }

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

#endif
