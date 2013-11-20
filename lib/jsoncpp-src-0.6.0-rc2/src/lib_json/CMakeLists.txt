set(JSON_FILES json_reader.cpp json_value.cpp json_writer.cpp)

set(CMAKE_CXX_FLAGS "-g -O2 ${CMAKE_CXX_FLAGS}")

add_library(json STATIC ${JSON_FILES})
set_property(TARGET json PROPERTY COMPILE_DEFINITIONS)
