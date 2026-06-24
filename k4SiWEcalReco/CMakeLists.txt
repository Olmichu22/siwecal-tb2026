cmake_minimum_required(VERSION 3.16)
project(k4SiWEcalReco LANGUAGES CXX)

# -- Gaudi/key4hep plugin package: turns the per-event `ecal` tree into EDM4hep
#    (CalorimeterHit + a per-event Cluster carrying the particle-ID shower
#    variables in its shapeParameters). Build under key4hep (source setup.sh).

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

find_package(Gaudi REQUIRED)
find_package(k4FWCore REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(podio REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
find_package(DD4hep REQUIRED COMPONENTS DDCore)

gaudi_add_module(k4SiWEcalRecoPlugins
  SOURCES
    src/components/EcalToEDM4hep.cpp
    src/components/EcalPidTransformer.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podio
    podio::podioRootIO
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::TreePlayer
    DD4hep::DDCore
)

target_include_directories(k4SiWEcalRecoPlugins
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

install(TARGETS k4SiWEcalRecoPlugins
  EXPORT k4SiWEcalRecoTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY options/ DESTINATION "${CMAKE_INSTALL_DATADIR}/k4SiWEcalReco"
  FILES_MATCHING PATTERN "*.py")
