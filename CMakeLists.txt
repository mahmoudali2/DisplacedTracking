cmake_minimum_required(VERSION 3.12)
project(CombinatorialTrackFinder LANGUAGES CXX)

# Standalone build of the CombinatorialTrackFinder (formerly DisplacedTracking /
# StandaloneMuonTracking). The sources mirror Tracking/ in k4RecTracker; only the
# files this component needs are kept here.

# Find required packages
find_package(DD4hep REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Gaudi REQUIRED)
find_package(ROOT REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(podio REQUIRED)

# GenFit (no upstream CMakeConfig — use our module in cmake/)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(GenFit REQUIRED)

# Define sources
set(sources
    components/CombinatorialTrackFinder.cpp
    components/MergeTrackerHitPlanes.cpp
    src/TrackingPrimitives.cpp
    src/LayerGeometry.cpp
    src/HelixExtrapolator.cpp
    src/SeedFinder.cpp
    src/TrackQuality.cpp
    src/genfit_interfaces/GenfitField.cpp
    src/genfit_interfaces/GenfitMaterialInterface.cpp
)

# Create the module
gaudi_add_module(CombinatorialTrackFinder
    SOURCES ${sources}
    LINK
        podio::podioIO
        Gaudi::GaudiKernel
        k4FWCore::k4FWCore
        k4FWCore::k4Interface
        EDM4HEP::edm4hep
        DD4hep::DDRec
        DD4hep::DDCore
        ROOT::Core
        ROOT::MathCore
        ROOT::Geom
        Eigen3::Eigen
        GenFit::genfit2
)

# Include directories
target_include_directories(CombinatorialTrackFinder PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include/genfit_interfaces>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install the target
install(TARGETS CombinatorialTrackFinder
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"
    COMPONENT dev
)

# Set up configuration generation directory
set(GAUDI_GENCONF_DIR "genConfDir")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${GAUDI_GENCONF_DIR})

# Unit tests (no Gaudi, no input files): `ctest` in the build directory.
enable_testing()

# Analytic helix crossings vs an independent RK4 Lorentz-force integration.
add_executable(test_HelixExtrapolator
    test/testHelixExtrapolator/test_HelixExtrapolator.cpp
    src/HelixExtrapolator.cpp)
target_include_directories(test_HelixExtrapolator PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(test_HelixExtrapolator PRIVATE Eigen3::Eigen)
add_test(NAME test_HelixExtrapolator COMMAND test_HelixExtrapolator)

# Toy calibration of the multiple-scattering-aware kink chi2.
add_executable(test_TrackQuality
    test/testTrackQuality/test_TrackQuality.cpp
    src/TrackQuality.cpp
    src/HelixExtrapolator.cpp)
target_include_directories(test_TrackQuality PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
target_link_libraries(test_TrackQuality PRIVATE Eigen3::Eigen)
add_test(NAME test_TrackQuality COMMAND test_TrackQuality)
