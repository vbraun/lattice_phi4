// Misc. Utilities

#ifndef UTIL__HPP
#define UTIL__HPP

#ifdef __GNUC__
#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#endif

#ifdef __GNUC__
#define noinline __attribute__ ((noinline))
#else
#define noinline
#endif

#endif /* UTIL__HPP */
