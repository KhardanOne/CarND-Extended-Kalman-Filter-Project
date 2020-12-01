#ifndef DEBUG_PRINT_
#define DEBUG_PRINT_

// use this macro to turn on/off extra couts:
#define VERBOSE_COUT_ENABLED 0


#if VERBOSE_COUT_ENABLED
#define debug_print(x) std::cout << x
#else
#define debug_print(x)
#endif

#endif  // DEBUG_PRINT_
