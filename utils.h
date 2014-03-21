#ifndef UTILS
#define UTILS

#define _min_(a, b) (a < b ? a : b)
#define _max_(a, b) (a < b ? b : a)
#define _minmod_(a, b) (a > 0 ? (b > 0 ? _min_(a, b) : 0) : (b > 0 ? 0 : - _min_(-a, -b)))
/* Int2Type */
template<int I>
struct Int2Type{ enum {value = I};};


/* Compile time error check */
template<int> struct CompileTimeError;
template<> struct CompileTimeError<true> {};
#define _static_check_(expr, msg) \
    { CompileTimeError<((expr) != 0)> ERROR_##msg; (void)ERROR_##msg; } 
#endif

/* define error messages */
#define ALARM(message)  std::cout << message << std::endl 
#define ASSERT_FILE_NOT_FOUND(AFILE) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "can not open file \"" << AFILE << "\" ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        exit(0) 
#define ASSERT_VARIABLE_OUT_OF_RANGE(AVAR) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "variable \"" << AVAR << "\" out of range ..."<< std::endl; \
        std::cerr << "*****************" << std::endl; \
        exit(0) 
#define ASSERT_DIMENSION_MISMATCH(AVAR1, AVAR2) \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "dimension mismatch for \"" << AVAR1 << "\", \"" << AVAR2 << std::endl; \
        std::cerr << "*****************" << std::endl; \
        exit(0) 
#define ASSERT_NOT_SUPPORTED \
        std::cerr << "***** ERROR *****" << std::endl; \
        std::cerr << "method not supported!" << std::endl; \
        std::cerr << "*****************" << std::endl; \
        exit(0) 
#define BREAKPOINT exit(0)

