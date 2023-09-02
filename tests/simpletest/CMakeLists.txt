add_library(simpletest
  simpletest.h
  simpletest.cpp
  )

target_include_directories(simpletest
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )