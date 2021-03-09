add_library(
  ActsPluginMinuit2 SHARED
  src/Minuit2Minimizer.cpp
  # src/Minuit2LayerBuilder.cpp
  # src/Minuit2Parser.cpp
  # src/Minuit2PrimitivesHelper.cpp
  # src/Minuit2SurfaceConverter.cpp
  )
target_include_directories(
  ActsPluginMinuit2
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginMinuit2
  PUBLIC ActsCore ROOT::Minuit2)

install(
  TARGETS ActsPluginMinuit2
  EXPORT ActsPluginMinuit2Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
