
# KD-Tree

Implicit KD-Tree in 2 and 3 dimensions. Uses either a std::array or a std::vector as container type. Falls back to linear search iff data is small.

    namespace kd

    struct vector { };
    struct array { };
    
    template<typename T, typename P = Point2<T>, typename Type = vector, std::size_t N = 0>
    struct Tree2D;

A Point2 and Point3 class are provided, or just drop in SFML's sf::Vector2 or sf::Vector3. The parameter N has no effect in case of a vector, which is the default container type.
