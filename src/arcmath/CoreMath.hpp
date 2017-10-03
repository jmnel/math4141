#pragma once

//#include <assert.h>

//#include "Mat2x2.hpp"

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstdint>

//#include "Fixed.hpp"

namespace arc
{
    
//#define CORE_PRECISION_DOUBLE
    
#ifdef CORE_PRECISION_DOUBLE
    
    /**
     *
     *  @brief Real-valued scalar with double precision.
     *
     */
    typedef double Real;
//    typedef int64_t Fixed;
    
#else
    
    /**
     *
     *  @brief Real-valued scalar with single precision.
     *
     */
    typedef float Real;
    
#endif
    
    typedef unsigned int uint;
    
    template<class T>
    double kronecker( T i, T j ) {
        return i == j ? 1.0f : 0.0f;
    };
    
    template <class T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
    
    //template <class T = Real>
	const double Pi = std::acos( -1.0 );
    //const T Pi = T( M_PI );
    
    template <class T>
    T abs( T x ) { return std::abs( x ); };
    
    template <class T>
    T copysign( T x, T y ) { return std::copysign( x, y ); }
    
    template <class T>
    T sqrt( T x ) { return std::sqrt( x ); };
    
    template <class T>
    T cos( T x ) { return std::cos( x ); };
    
    template <class T>
    T sin( T x ) { return std::sin( x ); };
    
    template <class T>
    T atan2( T y, T x ) { return std::atan2( y, x ); }
    
    template <class T>
    T fmod( T x, T y ) { return std::fmod( x, y ); }

    template <class T, class S>
    T pow(T x, S y) {
        return std::pow(x, y);
    }

    template <class T>
    T max(T x, T y) {
        return std::max(x, y);
    };

    template <class T>
    T min(T x, T y) {
        return std::min(x, y);
    };

    template <class T>
    T clamp( T x, T min, T max ) {
        x = x > max ? max : x;
        x = x < min ? min : x;
        return x;
    }
    
    template <class T>
    T smoothertep( T edge0, T edge1, T x ) {
        x = clamp( ( x - edge0 ) / ( edge1 - edge0 ), T( 0.0 ), T( 1.0 ) );
        x = x * x * x * ( x * ( x * T( 6.0 ) - T( 15.0 ) ) + T( 10.0 ) );
        return x;
    }
    
    uint nextPowerOf2( uint x );
    
}
