# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/cmake/372/bin/cmake

# The command to remove a file.
RM = /snap/cmake/372/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/thomas/Documents/github/Raytracer/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/thomas/Documents/github/Raytracer/src

# Include any dependencies generated for this target.
include CMakeFiles/lray.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/lray.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lray.dir/flags.make

CMakeFiles/lray.dir/Vector.cpp.o: CMakeFiles/lray.dir/flags.make
CMakeFiles/lray.dir/Vector.cpp.o: Vector.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thomas/Documents/github/Raytracer/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lray.dir/Vector.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lray.dir/Vector.cpp.o -c /home/thomas/Documents/github/Raytracer/src/Vector.cpp

CMakeFiles/lray.dir/Vector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lray.dir/Vector.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thomas/Documents/github/Raytracer/src/Vector.cpp > CMakeFiles/lray.dir/Vector.cpp.i

CMakeFiles/lray.dir/Vector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lray.dir/Vector.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thomas/Documents/github/Raytracer/src/Vector.cpp -o CMakeFiles/lray.dir/Vector.cpp.s

CMakeFiles/lray.dir/Raytracer.cpp.o: CMakeFiles/lray.dir/flags.make
CMakeFiles/lray.dir/Raytracer.cpp.o: Raytracer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/thomas/Documents/github/Raytracer/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/lray.dir/Raytracer.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lray.dir/Raytracer.cpp.o -c /home/thomas/Documents/github/Raytracer/src/Raytracer.cpp

CMakeFiles/lray.dir/Raytracer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lray.dir/Raytracer.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/thomas/Documents/github/Raytracer/src/Raytracer.cpp > CMakeFiles/lray.dir/Raytracer.cpp.i

CMakeFiles/lray.dir/Raytracer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lray.dir/Raytracer.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/thomas/Documents/github/Raytracer/src/Raytracer.cpp -o CMakeFiles/lray.dir/Raytracer.cpp.s

# Object files for target lray
lray_OBJECTS = \
"CMakeFiles/lray.dir/Vector.cpp.o" \
"CMakeFiles/lray.dir/Raytracer.cpp.o"

# External object files for target lray
lray_EXTERNAL_OBJECTS =

lray: CMakeFiles/lray.dir/Vector.cpp.o
lray: CMakeFiles/lray.dir/Raytracer.cpp.o
lray: CMakeFiles/lray.dir/build.make
lray: CMakeFiles/lray.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/thomas/Documents/github/Raytracer/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable lray"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lray.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lray.dir/build: lray

.PHONY : CMakeFiles/lray.dir/build

CMakeFiles/lray.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lray.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lray.dir/clean

CMakeFiles/lray.dir/depend:
	cd /home/thomas/Documents/github/Raytracer/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/thomas/Documents/github/Raytracer/src /home/thomas/Documents/github/Raytracer/src /home/thomas/Documents/github/Raytracer/src /home/thomas/Documents/github/Raytracer/src /home/thomas/Documents/github/Raytracer/src/CMakeFiles/lray.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lray.dir/depend

