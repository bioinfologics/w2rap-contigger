find_path(JEMALLOC_INCLUDE_DIR
          NAMES jemalloc/jemalloc.h
          HINTS ${JEMALLOC_ROOT_DIR}/include)

find_library(JEMALLOC_LIBRARIES
             NAMES jemalloc
             HINTS ${JEMALLOC_ROOT_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(jemalloc DEFAULT_MSG JEMALLOC_LIBRARIES JEMALLOC_INCLUDE_DIR)

mark_as_advanced(
        JEMALLOC_LIBRARIES
        JEMALLOC_INCLUDE_DIR)
