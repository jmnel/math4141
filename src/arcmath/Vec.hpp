#pragma once

#include <cassert>
#include <iostream>

#include "CoreMath.hpp"

namespace arc {
    
    namespace details {
        template<bool>
        struct rule_is_greater_than_4;
        
        template<>
        struct rule_is_greater_than_4<true> {};
        
        template<>
        struct rule_is_greater_than_4<false> {};
        
        template<class T, size_t N, size_t M>
        class inner_storage : rule_is_greater_than_4<( M > 4 )> {
        public:
            T x, y, z, w;
        private:
            T more_data[ N - 4 ];
        };
        
        template<class T, size_t N>
        class inner_storage<T, 2, N> {
        public:
            T x, y;
        };
        
        template<class T, size_t N>
        class inner_storage<T, 3, N> {
        public:
            T x, y, z;
        };
        
        template<class T, size_t N>
        class inner_storage<T, 4, N> {
        public:
            T x, y, z, w;
        };
    };
    
    /**
     *
     *  @brief <#description#>
     *
     */
    template<class T, size_t N>
    class Vec : public details::inner_storage<T, N, N> {
    public:
        static_assert( N > 1, "Size of 'n-tuple' must be > 1." );
        Vec();
        Vec( T k );
        Vec( T x, T y );
        Vec( T x, T y, T z );
        Vec( T x, T y, T z, T w );
        
        template<class S>
        Vec( Vec<S, N> const &v );
        
//        Vec( Vec<T, 2> v, T z );
//        Vec( Vec<T, 2> v, T z, T w );
//        Vec( Vec<T, 3> v, T w );

// const size_t size();
        size_t size() const;
        T &operator[]( const size_t i );
        T const &operator[]( const size_t i ) const;
        Vec<T, N> &operator=( Vec<T, N> const &t );
//        friend std::ostream &operator<<( std::ostream &os, Vec<T, N> const &t );
       
        std::string toString() const;

        static const Vec zero() { return Vec( 0.0 ); }
        static const Vec<T, 2> left() { return Vec( -1.0, 0.0 ); }
        static const Vec<T, 2> right() { return Vec( 1.0, 0.0 ); }
        static const Vec<T, 2> down() { return Vec( 0.0, -1.0 ); }
        static const Vec<T, 2> up() { return Vec( 0.0, 1.0 ); }
        
    };
    
    template<class T, size_t N>
    Vec<T, N> operator+( Vec<T, N> const &t );
    
    template<class T, size_t N>
    Vec<T, N> operator-( Vec<T, N> const &t );
    
    template<class T, size_t N>
    Vec<T, N> operator+( Vec<T, N> const &t1, Vec<T, N> const &t2 );
    
    template<class T, size_t N>
    Vec<T, N> operator-( Vec<T, N> const &t1, Vec<T, N> const &t2 );

    template <class S, class T, size_t N>
    Vec<T, N> operator*(S const &s, Vec<T, N> const &t2);

    template <class S, class T, size_t N>
    Vec<T, N> operator*( Vec<T, N> const &t1, S const &s );
    
    template <class S, class T, size_t N>
    Vec<T, N> operator/( Vec<T, N> const &t1, S const &s );
    
    template<class T, size_t N>
    bool operator==( Vec<T, N> const &t1, Vec<T, N> const &t2 );
    
    template<class T, size_t N>
    bool operator!=( Vec<T, N> const &t1, Vec<T, N> const &t2 );

    template<class T, size_t N>
    inline std::ostream &operator<<( std::ostream &os, Vec<T, N> const &t );

    typedef Vec<float, 2> Vec2f;
    typedef Vec<double, 2> Vec2d;
    typedef Vec<float, 3> Vec3f;
    typedef Vec<double, 3> Vec3d;
    typedef Vec<float, 4> Vec4f;
    typedef Vec<double, 4> Vec4d;
    typedef Vec<Real, 2> Vec2;
    typedef Vec<Real, 3> Vec3;
    typedef Vec<Real, 4> Vec4;

    /**
     *  @brief Calculate vector dot product.
     *
     *  @details Standard inner dot product of vectors \f$\vec{a}\f$ and
     * \f$\vec{b}\f$.
     *
     *  @param v1 Vector a
     *  @param v2 Vector b
     *  @return Result of 2D dot product
     *
     */
    template <class T, size_t N>
    T dot(Vec<T, N> const &v1, Vec<T, N> const &v2) {
        T product = T(0.0);
        for (size_t i = 0; i < N; i++) {
            product = product + v1[i] * v2[i];
        }
        return product;
    }

    template <class T>
    T cross(Vec<T, 2> const &v1, Vec<T, 2> const &v2) {
        return v1.x * v2.y - v1.y * v2.x;
    }

    template <class T>
    Vec<T,3> cross(Vec<T,3> const &v1, Vec<T,3> const &v2 ) {
        
        Vec<T,3> result =
        {
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x
        };
        
        return result;
    }

    /**
     *  @brief Calculate length squared of 2D vector.
     *
     *  @details Simply returns dot product of vector with itself. Uses function
     * [dot](@ref dot).
     *
     *  @param v Vector to calculate length squared of
     *
     *  @return Returns length squared of vector
     *
     */
    template <class T, size_t N>
    T lengthSquared(Vec<T, N> const &v) {
        T sum = T(0.0);
        for (size_t i = 0; i < N; i++) {
            sum = sum + v[i] * v[i];
        }
        return sum;
    }

    /**
     *  @brief Calculate length of vector.
     *
     *  @details Gives the square root of [lengthSquared](@ref lengthSquared) to
     * calculate length of vector.
     *
     *  @param v Vector to calculate length of
     *
     *  @return Returns length of vector
     *
     */
    template <class T, size_t N>
    T length(Vec<T, N> const &v) {
        ///@todo Implement as universal sqrt function
        return std::sqrt(lengthSquared(v));
    }

    /**
     *  @brief Calculates normalized unit vector of input vector.
     *
     *  @details Divides vector by length.
     *
     *  @param v Vector to normalize
     *
     *  @return Returns normalized vector
     *
     */
    template <class T, size_t N>
    Vec<T, N> normalize(Vec<T, N> const &v) {
        T dist = length(v);
        assert(dist > 0.0);

        return v / dist;
    }
}

#include "Vec.inl"
