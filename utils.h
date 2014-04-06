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

/* set up initial value from nc file */
#define _ncInitialize1_(var) var.initialize(setups::ncvar[#var].domain()); var = setups::ncvar[#var] 
#define _ncInitialize2_(var1, var2) _ncInitialize1_(var1); _ncInitialize1_(var2)
#define _ncInitialize3_(var1, var2, var3) _ncInitialize1_(var1); _ncInitialize2_(var2, var3)
#define _ncInitialize4_(var1, var2, var3, var4) \
    _ncInitialize1_(var1); _ncInitialize3_(var2, var3, var4)
#define _ncInitialize5_(var1, var2, var3, var4, var5) \
    _ncInitialize1_(var1); _ncInitialize4_(var2, var3, var4, var5)
#define _ncInitialize6_(var1, var2, var3, var4, var5, var6) \
    _ncInitialize1_(var1); _ncInitialize5_(var2, var3, var4, var5, var6)
#define _ncInitialize7_(var1, var2, var3, var4, var5, var6, var7) \
    _ncInitialize1_(var1); _ncInitialize6_(var2, var3, var4, var5, var6, var7)
#define _ncInitialize8_(var1, var2, var3, var4, var5, var6, var7, var8) \
    _ncInitialize1_(var1); _ncInitialize7_(var2, var3, var4, var5, var6, var7, var8)
#define _ncInitialize9_(var1, var2, var3, var4, var5, var6, var7, var8, var9) \
    _ncInitialize1_(var1); _ncInitialize8_(var2, var3, var4, var5, var6, var7, var8, var9)
#define _ncInitialize10_(var1, var2, var3, var4, var5, var6, var7, var8, var9, var10) \
    _ncInitialize1_(var1); _ncInitialize9_(var2, var3, var4, var5, var6, var7, var8, var9, var10)
