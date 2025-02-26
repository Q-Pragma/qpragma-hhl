# PYBIND LINALG QPU
message("-- Linking pybind Linalg QPU for demo")
set(qpu qpragma-pybindlinalg)

# DISPLAY QPU
#message("-- Linking Display QPU for demo")
#set(qpu displayqpu)


set(cpp-source src/main.cpp)

# Compile using Q-Pragma plugin
set_property(SOURCE ${cpp-source} PROPERTY COMPILE_FLAGS "-fplugin=$<TARGET_FILE:clang-qpragma> -I /utils")

# Get target name
#get_filename_component(cpp-target ${cpp-source} NAME_WE)

# Create executable
add_executable(hybrid_hhl EXCLUDE_FROM_ALL ${cpp-source})
add_dependencies(hybrid_hhl clang-qpragma)

# Link QPU
target_link_libraries(hybrid_hhl PUBLIC qpragma ${qpu})
