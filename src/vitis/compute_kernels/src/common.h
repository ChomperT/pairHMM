#ifndef __COMMON_H__
#define __COMMON_H__

#define F64 0
#define CORE 6
#define MAX_SEQ_LEN 351
#define GROUP 19

#ifdef __SYNTHESIS__
#include <hls_stream.h>
#include <hls_math.h>
#include <ap_int.h>
typedef ap_uint<128> u128;
typedef ap_uint<256> u256;
typedef ap_uint<512> u512;
typedef ap_uint<10> u10;
typedef ap_uint<10> btype;
typedef ap_uint<8> len_c_t;
#define MIN(a, b) hls::min(a, b)
#else
#include <stdio.h>
#include <stdint.h>
#include <cmath>
#include <stddef.h>
#include <cstdint>
#include <string>
#include <immintrin.h>

typedef __m256i u256;
typedef unsigned char len_c_t;
typedef unsigned short u10;
typedef unsigned short btype;
#define MIN(a, b) std::min(a, b)
#endif

#include <limits>

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long u64;

typedef float f32;
typedef double f64;
typedef u64 intype;

#if F64
typedef f64 precision;
#include "constantd.h"
#ifndef __SYNTHESIS__
#define LOG10(x) log10(x)
#define POW(x, y) pow(x, y)
#else
#define LOG10(x) hls::log10(x)
#define POW(x, y) hls::pow(x, y)
#endif
#define LMAX ((precision)0X1.330CF3D4EDA85P+8)
constexpr precision MAX_VALUE = std::numeric_limits<precision>::max();
#else


typedef f32 precision;
#include "constantf.h"
#ifndef __SYNTHESIS__
#define LOG10(x) log10(x)
#define POW(x, y) powf(x, y)
#else
#define LOG10(x) hls::log10(x)
#define POW(x, y) hls::powf(x, y)
#endif
#define LMAX ((precision)0x1.2a9f2b5e2e739p+5)
constexpr precision MAX_VALUE = std::numeric_limits<precision>::max()/16.0f;
#endif


#endif
