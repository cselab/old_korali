#ifndef ENGINE_CMAES_UTILS_HPP
#define ENGINE_CMAES_UTILS_HPP

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
	#define ChangeDir _chdir
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
	#define ChangeDir chdir
#endif

#endif /* ENGINE_CMAES_UTILS_HPP */
