# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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
CMAKE_COMMAND = /n/lanl_linux/packages/cmake28/bin/cmake

# The command to remove a file.
RM = /n/lanl_linux/packages/cmake28/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /n/lanl_linux/packages/cmake28/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/quanb/arctic_polyGen/mkdual-1.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quanb/arctic_polyGen/mkdual-1.1/build

# Include any dependencies generated for this target.
include CMakeFiles/mkdual.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mkdual.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mkdual.dir/flags.make

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o: CMakeFiles/mkdual.dir/flags.make
CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o: ../src/MESH_MakeDual.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/quanb/arctic_polyGen/mkdual-1.1/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o   -c /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual.c

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.i"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual.c > CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.i

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.s"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual.c -o CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.s

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.requires:
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.requires

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.provides: CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.requires
	$(MAKE) -f CMakeFiles/mkdual.dir/build.make CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.provides.build
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.provides

CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.provides.build: CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.provides.build

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o: CMakeFiles/mkdual.dir/flags.make
CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o: ../src/MESH_MakeDual2.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/quanb/arctic_polyGen/mkdual-1.1/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o   -c /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual2.c

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.i"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual2.c > CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.i

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.s"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual2.c -o CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.s

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.requires:
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.requires

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.provides: CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.requires
	$(MAKE) -f CMakeFiles/mkdual.dir/build.make CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.provides.build
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.provides

CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.provides.build: CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.provides.build

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o: CMakeFiles/mkdual.dir/flags.make
CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o: ../src/MESH_MakeDual3.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/quanb/arctic_polyGen/mkdual-1.1/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o   -c /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual3.c

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.i"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual3.c > CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.i

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.s"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_MakeDual3.c -o CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.s

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.requires:
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.requires

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.provides: CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.requires
	$(MAKE) -f CMakeFiles/mkdual.dir/build.make CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.provides.build
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.provides

CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.provides.build: CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o
.PHONY : CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.provides.build

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o: CMakeFiles/mkdual.dir/flags.make
CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o: ../src/MESH_HasPlanarFaces.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/quanb/arctic_polyGen/mkdual-1.1/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o   -c /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_HasPlanarFaces.c

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.i"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_HasPlanarFaces.c > CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.i

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.s"
	/home/quanb/local/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/quanb/arctic_polyGen/mkdual-1.1/src/MESH_HasPlanarFaces.c -o CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.s

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.requires:
.PHONY : CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.requires

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.provides: CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.requires
	$(MAKE) -f CMakeFiles/mkdual.dir/build.make CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.provides.build
.PHONY : CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.provides

CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.provides.build: CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o
.PHONY : CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.provides.build

# Object files for target mkdual
mkdual_OBJECTS = \
"CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o" \
"CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o" \
"CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o" \
"CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o"

# External object files for target mkdual
mkdual_EXTERNAL_OBJECTS =

libmkdual.a: CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o
libmkdual.a: CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o
libmkdual.a: CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o
libmkdual.a: CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o
libmkdual.a: CMakeFiles/mkdual.dir/build.make
libmkdual.a: CMakeFiles/mkdual.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C static library libmkdual.a"
	$(CMAKE_COMMAND) -P CMakeFiles/mkdual.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mkdual.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mkdual.dir/build: libmkdual.a
.PHONY : CMakeFiles/mkdual.dir/build

CMakeFiles/mkdual.dir/requires: CMakeFiles/mkdual.dir/src/MESH_MakeDual.c.o.requires
CMakeFiles/mkdual.dir/requires: CMakeFiles/mkdual.dir/src/MESH_MakeDual2.c.o.requires
CMakeFiles/mkdual.dir/requires: CMakeFiles/mkdual.dir/src/MESH_MakeDual3.c.o.requires
CMakeFiles/mkdual.dir/requires: CMakeFiles/mkdual.dir/src/MESH_HasPlanarFaces.c.o.requires
.PHONY : CMakeFiles/mkdual.dir/requires

CMakeFiles/mkdual.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mkdual.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mkdual.dir/clean

CMakeFiles/mkdual.dir/depend:
	cd /home/quanb/arctic_polyGen/mkdual-1.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quanb/arctic_polyGen/mkdual-1.1 /home/quanb/arctic_polyGen/mkdual-1.1 /home/quanb/arctic_polyGen/mkdual-1.1/build /home/quanb/arctic_polyGen/mkdual-1.1/build /home/quanb/arctic_polyGen/mkdual-1.1/build/CMakeFiles/mkdual.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mkdual.dir/depend

