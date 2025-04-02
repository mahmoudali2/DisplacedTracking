cmake_minimum_required(VERSION 3.12)
project(KalmanTracking LANGUAGES CXX)

# Find required packages
find_package(DD4hep REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Gaudi REQUIRED)
find_package(ROOT REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(podio REQUIRED)

#include(cmake/Key4hepConfig.cmake)

# Define sources
set(sources
    components/KalmanTracking.cpp
)

# Create the module
gaudi_add_module(KalmanTracking
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
)

# Include directories
target_include_directories(KalmanTracking PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install the target
install(TARGETS KalmanTracking
    EXPORT ${CMAKE_PROJECT_NAME}Targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"
    COMPONENT dev
)

# Set up configuration generation directory
set(GAUDI_GENCONF_DIR "genConfDir")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${GAUDI_GENCONF_DIR})
