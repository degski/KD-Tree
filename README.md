
# KD-Tree

A static KD-Tree in 2 and 3 dimensions. Insert and delete are not implemented and the Tree should simply be reconstructed, which is fairly fast. The Tree is implicit, i.e. flat with no overhead. Uses either a std::array or a std::vector as container type. Falls back to linear search iff data is small. The size of the container is next ( N ^ 2 ) - 1, so on average there is 25% wasted space.

The implementation uses recursion for both construction and nn-search.

    namespace kd

    struct vector { };
    struct array { };
    
    template<typename T, typename P = Point2<T>, typename Type = vector, std::size_t N = 0>
    struct Tree2D;

A Point2 and Point3 class are provided, or just drop in SFML's sf::Vector2 or sf::Vector3. The parameter N has no effect in case of a vector, which is the default container type.

No map is implemented, but one can return an index of the node that's closest. Simply make lookup tables for any "mapped" data. 

TODO: The case of a power of 2 will have to be looked at, as this might be common, but actually has the maximum space wastage (50%). The root of the tree could be taken out of the data-strucutre to accomodate that, which will give no waste as a result.