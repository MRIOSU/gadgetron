if (WIN32)
  add_definitions(-D__BUILD_GADGETRON_CARTESIAN__)
endif ()

include_directories(
  ${CMAKE_SOURCE_DIR}/gadgets/mri_core
  )

add_library(gadgetron_cartesian SHARED 
  gadgetron_cartesian_export.h
  CartesianToGenericGadget.h
  CartesianToGenericGadget.cpp
)

set_target_properties(gadgetron_cartesian PROPERTIES VERSION ${GADGETRON_VERSION_STRING} SOVERSION ${GADGETRON_SOVERSION})

target_link_libraries(gadgetron_cartesian gadgetron_toolbox_cpucore 
  gadgetron_gadgetbase
  gadgetron_toolbox_log
  ${ISMRMRD_LIBRARIES}
  optimized ${ACE_LIBRARIES} debug ${ACE_DEBUG_LIBRARY}
  )

if(ARMADILLO_FOUND)
    target_link_libraries(gadgetron_cartesian gadgetron_toolbox_cpucore_math )
endif()

install (FILES  gadgetron_cartesian_export.h
                CartesianToGenericGadget.h
                DESTINATION ${GADGETRON_INSTALL_INCLUDE_PATH} COMPONENT main)

install (TARGETS gadgetron_cartesian DESTINATION lib COMPONENT main)
