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
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Spectrum)
find_package(DD4hep REQUIRED COMPONENTS DDCore)
# LCIO: only EcalLcioDecoder needs it (reads the EUDAQ ROC_run_*_tp.slcio files).
find_package(LCIO REQUIRED)
# ACTS tracking stage (ACTSGeoSvc, ShowerTagger, SiPadMeasConverter,
# ACTSProtoTracker) -- same dependency pin as siwecal_k4sim's CMakeLists.txt,
# both repos share the key4hep release (.key4hep-release).
find_package(Acts 44.3.0 REQUIRED COMPONENTS Core PluginDD4hep)
find_package(k4ActsTracking REQUIRED)

gaudi_add_module(k4SiWEcalRecoPlugins
  SOURCES
    src/components/EcalToEDM4hep.cpp
    src/components/EcalPidTransformer.cpp
    src/components/EcalEventCut.cpp
    src/components/EcalRawDecoder.cpp
    src/components/EcalLcioDecoder.cpp
    src/components/PedestalMipCalibrator.cpp
    src/components/EcalEventBuilder.cpp
    src/components/ACTSGeoSvc.cpp
    src/components/ShowerTagger.cpp
    src/components/SiPadMeasConverter.cpp
    src/components/ACTSProtoTracker.cpp
  LINK
    Gaudi::GaudiKernel
    k4FWCore::k4FWCore
    EDM4HEP::edm4hep
    podio::podio
    podio::podioRootIO
    ROOT::Core ROOT::RIO ROOT::Tree ROOT::TreePlayer ROOT::Spectrum
    DD4hep::DDCore
    LCIO::lcio
    Acts::Core
    Acts::PluginDD4hep
    k4ActsTracking::k4ActsTracking
)

target_include_directories(k4SiWEcalRecoPlugins
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ACTS headers: derived from the Acts::Core target find_package(Acts) sets up,
# not hardcoded to one key4hep release's install path -- the release (and so
# the hashed acts/<version>-<hash>/include directory) is pinned in
# .key4hep-release and changes between stacks. Same workaround as
# siwecal_k4sim's top-level CMakeLists.txt.
get_target_property(_acts_include_dirs Acts::Core INTERFACE_INCLUDE_DIRECTORIES)
if(_acts_include_dirs)
    target_include_directories(k4SiWEcalRecoPlugins PRIVATE ${_acts_include_dirs})
endif()

install(TARGETS k4SiWEcalRecoPlugins
  EXPORT k4SiWEcalRecoTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY options/ DESTINATION "${CMAKE_INSTALL_DATADIR}/k4SiWEcalReco"
  FILES_MATCHING PATTERN "*.py")
