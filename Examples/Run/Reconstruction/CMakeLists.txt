# shared code 
add_subdirectory(Common)

# Truth fitting for Generic detector 
add_executable(
  ActsExampleTruthTracksGeneric
  GenericRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

# Truth fitting for Telescope detector
add_executable(
  ActsExampleTruthTracksTelescope
  TelescopeRecTruthTracks.cpp)
target_link_libraries(
  ActsExampleTruthTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# CKF reco for Generic detector 
add_executable(
  ActsExampleCKFTracksGeneric
  GenericRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

# CKF reco for Telescope detector 
add_executable(
  ActsExampleCKFTracksTelescope
  TelescopeRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# Chi2 for Generic Detector
add_executable(
  ActsExampleChi2TracksGeneric
  GenericRecChi2Tracks.cpp)
target_link_libraries(
  ActsExampleChi2TracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)
# chi2 reco for Telescope detector 
add_executable(
  ActsExampleChi2TracksTelescope
  TelescopeRecChi2Tracks.cpp)
target_link_libraries(
  ActsExampleChi2TracksTelescope
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorTelescope)

# Seeding for Generic detector
add_executable(
  ActsExampleSeedingGeneric
  GenericSeedingExample.cpp)
target_link_libraries(
  ActsExampleSeedingGeneric
  PRIVATE
    ActsExamplesRecTracksCommon
    ActsExamplesDetectorGeneric
)

install(
  TARGETS
    ActsExampleTruthTracksGeneric
    ActsExampleTruthTracksTelescope
    ActsExampleCKFTracksGeneric
    ActsExampleCKFTracksTelescope
    ActsExampleChi2TracksGeneric
    ActsExampleChi2TracksTelescope
    ActsExampleSeedingGeneric 
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
