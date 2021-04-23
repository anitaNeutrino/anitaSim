cmake_minimum_required(VERSION 3.1.0)
project(anitaSim)
set(libname "AnitaSimIcemc")

### Site-specific options
option(ANITA3_EVENTREADER "Load ANITA3 EventReader" ON) 

## ROOT 
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathMore MathCore RIO Hist Tree Net)

message("ROOT_VERSION is set to ${ROOT_VERSION}")
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include )

set (CMAKE_CXX_STANDARD 11)

FILE (GLOB SOURCE_FILES "src/*.cc")
FILE (GLOB INCLUDE_FILES "include/*.h") 

set(LIBS ${ROOT_LIBRARIES} MathMore ) 

if (ANITA3_EVENTREADER) 
  add_definitions(-DANITA3_EVENTREADER)
  set(LIBS ${LIBS} AnitaEvent ) 
endif()


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
    message("ANITA_UTIL_INSTALL_DIR is found and set to $ENV{ANITA_UTIL_INSTALL_DIR}")
    add_definitions(-DANITA_UTIL_EXISTS)
    set (LIBS ${LIBS} RootFftwWrapper ) 
    set(ANITA_UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(ANITA_UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(LD_ANITA_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
    set(INC_ANITA_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)
    set(ANITA_UTIL_ETC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/etc)
    link_directories(${ANITA_UTIL_LIB_DIR})
    include_directories(${ANITA_UTIL_LIB_DIR})
    include_directories(${ANITA_UTIL_INC_DIR})	
 endif()

set(CMAKE_CXX_FLAGS_DEFAULT 
  "-O2 -g -pipe -m64 -pthread -W -Wall -Wextra -Woverloaded-virtual" 
  CACHE STRING "c++ Flags used during default icemc builds" 
  FORCE ) 

#We need to compile stuff into a library unless we want it to be compiled a million times over
set(DICTNAME G__${libname})
ROOT_GENERATE_DICTIONARY(${DICTNAME} ${INCLUDE_FILES} LINKDEF LinkDef.h)
if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()



#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ZLIB_LIBRARIES} RootFftwWrapper AnitaEvent AnitaIceMC ${ROOT_LIBRARIES} ${FFTW_LIBRARIES})


if(DEFINED ENV{ANITA_UTIL_INSTALL_DIR})
  set(UTIL_LIB_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
  set(UTIL_INC_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/include)
  set(UTIL_BIN_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/bin)
  set(UTIL_SHARE_DIR $ENV{ANITA_UTIL_INSTALL_DIR}/share)
  set(LD_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/lib)
  set(INC_UTIL $ENV{ANITA_UTIL_INSTALL_DIR}/include)

  install (FILES ${INCLUDE_FILES} DESTINATION ${UTIL_INC_DIR})
  install (TARGETS ${libname}
    ARCHIVE DESTINATION ${UTIL_LIB_DIR}
    LIBRARY DESTINATION ${UTIL_LIB_DIR}
    RUNTIME DESTINATION ${UTIL_BIN_DIR})
  #Only needed for ROOT6
  if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    install (FILES ${PROJECT_BINARY_DIR}/${DICTNAME}_rdict.pcm DESTINATION ${UTIL_LIB_DIR})
    #  install (FILES ${PROJECT_BINARY_DIR}/lib${libname}.rootmap DESTINATION ${UTIL_LIB_DIR})
  endif()
endif()

macro(add_binary bin) 
  add_executable(${bin} test/${bin}.cc)
  #target_link_libraries(${bin} icemc_objs ${LIBS})
  target_link_libraries(${bin} ${LIBS} ${libname})
endmacro() 

add_binary(simulateAnita)
