# FindGenFit.cmake
# Locates the GenFit tracking library (https://github.com/GenFit/GenFit).
#
# GenFit does not ship a CMakeConfig file, so this module finds it via
# find_path / find_library, honouring CMAKE_PREFIX_PATH (set by the
# key4hep nightlies) and the optional GENFIT_INSTALL_DIR hint.
#
# Imported target
#   GenFit::genfit2
#
# Result variables (also set for projects that don't use the imported target)
#   GENFIT_FOUND
#   GENFIT_INCLUDE_DIRS
#   GENFIT_LIBRARIES

find_path(GENFIT_INCLUDE_DIR
    NAMES AbsBField.h
    PATH_SUFFIXES include
    HINTS ENV GENFIT_INSTALL_DIR
)

find_library(GENFIT_LIBRARY
    NAMES genfit2
    PATH_SUFFIXES lib lib64
    HINTS ENV GENFIT_INSTALL_DIR
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GenFit
    REQUIRED_VARS GENFIT_LIBRARY GENFIT_INCLUDE_DIR
    VERSION_VAR   GENFIT_VERSION
)

if(GENFIT_FOUND AND NOT TARGET GenFit::genfit2)
    add_library(GenFit::genfit2 SHARED IMPORTED)
    set_target_properties(GenFit::genfit2 PROPERTIES
        IMPORTED_LOCATION             "${GENFIT_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GENFIT_INCLUDE_DIR}"
    )
endif()

set(GENFIT_INCLUDE_DIRS "${GENFIT_INCLUDE_DIR}")
set(GENFIT_LIBRARIES    "${GENFIT_LIBRARY}")

mark_as_advanced(GENFIT_INCLUDE_DIR GENFIT_LIBRARY)
