
#ifndef UTILSLIB_EXPORTS_H
#define UTILSLIB_EXPORTS_H

#ifdef UtilsLib_EXPORTS_BUILT_AS_STATIC
#  define UTILSLIB_EXPORTS
#  define UTILSLIB_NO_EXPORT
#else
#  ifndef UTILSLIB_EXPORTS
#    ifdef UtilsLib_EXPORTS
        /* We are building this library */
#      define UTILSLIB_EXPORTS __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define UTILSLIB_EXPORTS __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef UTILSLIB_NO_EXPORT
#    define UTILSLIB_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef UTILSLIB_DEPRECATED
#  define UTILSLIB_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef UTILSLIB_DEPRECATED_EXPORT
#  define UTILSLIB_DEPRECATED_EXPORT UTILSLIB_EXPORTS UTILSLIB_DEPRECATED
#endif

#ifndef UTILSLIB_DEPRECATED_NO_EXPORT
#  define UTILSLIB_DEPRECATED_NO_EXPORT UTILSLIB_NO_EXPORT UTILSLIB_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef UTILSLIB_NO_DEPRECATED
#    define UTILSLIB_NO_DEPRECATED
#  endif
#endif

#endif /* UTILSLIB_EXPORTS_H */
