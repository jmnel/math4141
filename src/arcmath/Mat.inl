// include "FunctionsGeometric.hpp"

namespace arc {
    template<class T, size_t M, size_t N>
    const Mat<T, M, N> Mat<T, M, N>::identity() {
        Mat<T, M, N> identity;
        for( size_t i = 0; i < M; i++ ) {
            for( size_t j = 0; j < N; j++ ) {
                identity[ i ][ j ] = i == j ? T( 1.0 ) : T( 0.0 );
            }
        }
        return identity;
    }
    
    template<class T, size_t M, size_t N>
    const Mat<T, M, N> Mat<T, M, N>::zero() {
        return Mat<T, M, N>();
    }
    
    // -- Constructors --
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::Mat() {
        for( size_t i = 0; i < M * N; i++ ) {
            elements[ i ] = T( 0.0 );
        }
//        for( size_t i = 0; i < M; i++ ) {
//            for( size_t j = 0; j < N; j++ ) {
//                rows[ i ][ j ] = T( 0.0 );
//            }
//        }
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::Mat( T const &k ) {
        for( size_t i = 0; i < M; i++ ) {
            for( size_t j = 0; j < N; j++ ) {
                elements[ i * M + j ] = i == j ? k : T( 0.0 );
            }
        }
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::Mat( Vec<T, N> const &row ) {
        for( size_t i = 0; i < M; i++ ) {
            for( size_t j = 0; j < N; j++ ) {
//                rows[ i ][ j ] = row[ j ];
                elements[ i * M + j ] = row[ j ];
            }
        }
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::Mat( Vec<T, N> const ( &rows )[ M ] ) {
        for( size_t i = 0; i < M; i++ ) {
            for( size_t j = 0; j < N; j++ ) {
                this->elements[ i * M + j ] = rows[ i ][ j ];
            }
        }
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::Mat( std::initializer_list<std::initializer_list<T>> values ) {
        assert( values.size() == M );
        for( size_t i = 0; i < M; i++ ) {
            auto row = values.begin() + i;
            assert( row->size() == N );
            for( size_t j = 0; j < N; j++ ) {
                elements[ i * M + j ] = *( row->begin() + j );
            }
        }
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N>::~Mat() {
        
    }

    /*
    template<class T, size_t M, size_t N>
    const size_t Mat<T, M, N>::size() {
        return M;
    }
    */

    template<class T, size_t M, size_t N>
    size_t Mat<T, M, N>::size() const {
        return M;
    }
    
    template<class T, size_t M, size_t N>
    Vec<T, N> &Mat<T, M, N>::operator[]( const size_t i ) {
        assert( i < M );
//        return rows[ i ];
        return *(Vec<T, N> *)&elements[ i * M ];
    }
    
    template<class T, size_t M, size_t N>
    Vec<T, N> const &Mat<T, M, N>::operator[]( const size_t i ) const {
        assert( i < M );
//        return rows[ i ];
        return *(Vec<T, N> *)&elements[ i * M ];
//        return &rows[ i ];
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> &Mat<T, M, N>::operator=( const Mat<T, M, N> &m ) {
//        for( size_t i = 0; i < M; i++ ) {
//            this->rows[ i ] = m.rows[ i ];
//        }
        for( size_t i = 0; i < M * N; i++ ) {
            this->elements[ i ] = m.elements[ i ];
        }
        return *this;
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator+( Mat<T, M, N> const &m ) {
        return m;
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator-( Mat<T, M, N> const &m ) {
        return m * T( -1.0 );
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator+( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 ) {
        Mat<T, M, N> sum;
        for( size_t i = 0; i < M; i++ ) {
            sum[ i ] = m1[ i ] + m2[ i ];
        }
        return sum;
    }

    template<class T, size_t M, size_t N, size_t K>
    Mat<T, M, K> operator*( Mat<T, M, N> const &m1, Mat<T, N, K> const &m2 ) {
        Mat<T, M, K> product;
        for( size_t i = 0; i < M; i++ ) {
            for( size_t j = 0; j < K; j++ ) {
                Vec<T, N> col;
                for( size_t q = 0; q < N; q++ ) {
                    col[ q ] = m2[ q ][ j ];
                }
                product[ i ][ j ] = dot( m1[ i ], col );
            }
        }
        return product;
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator*( T const &s, Mat<T, M, N> const &m ) {
        Mat<T, M, N> product;
        for( size_t i = 0; i < M; i++ ) {
            product[ i ] = s * m[ i ];
        }
        return product;
    }
    
    template<class T, size_t M, size_t N>
    Mat<T, M, N> operator*( Mat<T, M, N> const &m, T const &s ) {
        return s * m;
    }
    
    template<class T, size_t M, size_t N>
    Vec<T, N> operator*( Mat<T, M, N> const &m, Vec<T, M> const &v ) {
        assert( "not implemented" );
    }
    
    template<class T, size_t M, size_t N>
    bool operator==( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 ) {
        for( size_t i = 0; i < M; i++ ) {
            if( m1[ i ] != m2[ i ] ) { return false; }
        }
        return true;
    }
    
    template<class T, size_t M, size_t N>
    bool operator!=( Mat<T, M, N> const &m1, Mat<T, M, N> const &m2 ) {
        return !( m1 == m2 );
    }
    
    template<class T, size_t M, size_t N>
    inline std::ostream &operator<<( std::ostream &os, Mat<T, M, N> const &m ) {
        for( size_t i = 0; i < M; i++ ) {
            os << m[ i ] << std::endl;
        }
        return os;
    }
    
}
