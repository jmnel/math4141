#pragma once

#include "Vec.hpp"

namespace arc {
    
    /**
     *
     *  @brief <#description#>
     *
     */
    template<class T, size_t M, size_t N>
    class Mat {
    private:
//        Vec<T, N> rows[ M ] = { T( 0.0 ) };
        T elements[ M * N ] = { T( 0.0 ) };
    public:
        static const Mat identity();
        static const Mat zero();
        Mat();
        Mat( T const &k );
        Mat( Vec<T, N> const &row );
        Mat( Vec<T, N> const ( &rows )[ M ] );
//        Mat( const T ( &rows )[ M ][ N ] );
        Mat( std::initializer_list<std::initializer_list<T>> values );
        
        virtual ~Mat();

        // const size_t size();
        size_t size() const;
        Vec<T, N> &operator[]( const size_t i );
        Vec<T, N> const &operator[]( const size_t i ) const;
        Mat &operator=( const Mat &m );
//        friend std::ostream &operator<<( std::ostream &os, Mat &m );
    };
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator+( Mat<T, M, N> const &m );
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator-( Mat<T, M, N> const &m );
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator+( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 );
    
    template<class T, size_t M, size_t N, size_t K>
    Mat<T, M, K> operator*( Mat<T, M, N> const &m1, Mat<T, N, K> const &m2 );
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator*( T const &s, Mat<T, M, N> const &m );
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator*( Mat<T, M, N> const &m, T const &s );
    
    template<class T, size_t M, size_t N>
    Vec<T, N> operator*( Mat<T, M, N> const &m, Vec<T, M> const &v );
    
    template<class T, size_t M, size_t N>
    bool operator==( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 );
    
    template<class T, size_t M, size_t N>
    bool operator!=( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 );
    
    template<class T, size_t M, size_t N>
    inline std::ostream &operator<<( std::ostream &os, Mat<T, M, N> const &m );

    using Mat4x4f = Mat<float, 4, 4>;
}

#include "Mat.inl"
