if(PKG_MESSAGE)
  if(LAMMPS_SIZES STREQUAL BIGBIG)
    message(FATAL_ERROR "The MESSAGE Package is not compatible with -DLAMMPS_BIGBIG")
  endif()
  option(MESSAGE_ZMQ "Use ZeroMQ in MESSAGE package" OFF)
  file(GLOB_RECURSE cslib_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.F
      ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.c
      ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/[^.]*.cpp)

  add_library(cslib STATIC ${cslib_SOURCES})
  install(TARGETS cslib EXPORT LAMMPS_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  if(BUILD_MPI)
    target_compile_definitions(cslib PRIVATE -DMPI_YES)
    set_target_properties(cslib PROPERTIES OUTPUT_NAME "csmpi")
  else()
    target_compile_definitions(cslib PRIVATE -DMPI_NO)
    target_include_directories(cslib PRIVATE ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src/STUBS_MPI)
    set_target_properties(cslib PROPERTIES OUTPUT_NAME "csnompi")
  endif()

  if(MESSAGE_ZMQ)
    target_compile_definitions(cslib PRIVATE -DZMQ_YES)
    find_package(ZMQ REQUIRED)
    target_include_directories(cslib PRIVATE ${ZMQ_INCLUDE_DIRS})
    target_link_libraries(cslib PUBLIC ${ZMQ_LIBRARIES})
  else()
    target_compile_definitions(cslib PRIVATE -DZMQ_NO)
    target_include_directories(cslib PRIVATE ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src/STUBS_ZMQ)
  endif()

  target_link_libraries(lammps PRIVATE cslib)
  target_include_directories(lammps PRIVATE ${LAMMPS_LIB_SOURCE_DIR}/message/cslib/src)
endif()
