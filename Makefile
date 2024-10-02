# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/deborah/Desktop/cppCodes/PhaseField

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/deborah/Desktop/cppCodes/PhaseField/03PhaseFieldCpp

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/deborah/Desktop/cppCodes/PhaseField/03PhaseFieldCpp/CMakeFiles /home/deborah/Desktop/cppCodes/PhaseField/03PhaseFieldCpp//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/deborah/Desktop/cppCodes/PhaseField/03PhaseFieldCpp/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named solid_library

# Build rule for target.
solid_library: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 solid_library
.PHONY : solid_library

# fast build rule for target.
solid_library/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/build
.PHONY : solid_library/fast

#=============================================================================
# Target rules for targets named s

# Build rule for target.
s : cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 s
.PHONY : s

# fast build rule for target.
s/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/s.dir/build.make CMakeFiles/s.dir/build
.PHONY : s/fast

main.o: main.cpp.o
.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/s.dir/build.make CMakeFiles/s.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i
.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/s.dir/build.make CMakeFiles/s.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s
.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/s.dir/build.make CMakeFiles/s.dir/main.cpp.s
.PHONY : main.cpp.s

mesh_interface/sources/Geometry.o: mesh_interface/sources/Geometry.cpp.o
.PHONY : mesh_interface/sources/Geometry.o

# target to build an object file
mesh_interface/sources/Geometry.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Geometry.cpp.o
.PHONY : mesh_interface/sources/Geometry.cpp.o

mesh_interface/sources/Geometry.i: mesh_interface/sources/Geometry.cpp.i
.PHONY : mesh_interface/sources/Geometry.i

# target to preprocess a source file
mesh_interface/sources/Geometry.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Geometry.cpp.i
.PHONY : mesh_interface/sources/Geometry.cpp.i

mesh_interface/sources/Geometry.s: mesh_interface/sources/Geometry.cpp.s
.PHONY : mesh_interface/sources/Geometry.s

# target to generate assembly for a file
mesh_interface/sources/Geometry.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Geometry.cpp.s
.PHONY : mesh_interface/sources/Geometry.cpp.s

mesh_interface/sources/Inclusion.o: mesh_interface/sources/Inclusion.cpp.o
.PHONY : mesh_interface/sources/Inclusion.o

# target to build an object file
mesh_interface/sources/Inclusion.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Inclusion.cpp.o
.PHONY : mesh_interface/sources/Inclusion.cpp.o

mesh_interface/sources/Inclusion.i: mesh_interface/sources/Inclusion.cpp.i
.PHONY : mesh_interface/sources/Inclusion.i

# target to preprocess a source file
mesh_interface/sources/Inclusion.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Inclusion.cpp.i
.PHONY : mesh_interface/sources/Inclusion.cpp.i

mesh_interface/sources/Inclusion.s: mesh_interface/sources/Inclusion.cpp.s
.PHONY : mesh_interface/sources/Inclusion.s

# target to generate assembly for a file
mesh_interface/sources/Inclusion.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Inclusion.cpp.s
.PHONY : mesh_interface/sources/Inclusion.cpp.s

mesh_interface/sources/Line.o: mesh_interface/sources/Line.cpp.o
.PHONY : mesh_interface/sources/Line.o

# target to build an object file
mesh_interface/sources/Line.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Line.cpp.o
.PHONY : mesh_interface/sources/Line.cpp.o

mesh_interface/sources/Line.i: mesh_interface/sources/Line.cpp.i
.PHONY : mesh_interface/sources/Line.i

# target to preprocess a source file
mesh_interface/sources/Line.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Line.cpp.i
.PHONY : mesh_interface/sources/Line.cpp.i

mesh_interface/sources/Line.s: mesh_interface/sources/Line.cpp.s
.PHONY : mesh_interface/sources/Line.s

# target to generate assembly for a file
mesh_interface/sources/Line.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Line.cpp.s
.PHONY : mesh_interface/sources/Line.cpp.s

mesh_interface/sources/LineLoop.o: mesh_interface/sources/LineLoop.cpp.o
.PHONY : mesh_interface/sources/LineLoop.o

# target to build an object file
mesh_interface/sources/LineLoop.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/LineLoop.cpp.o
.PHONY : mesh_interface/sources/LineLoop.cpp.o

mesh_interface/sources/LineLoop.i: mesh_interface/sources/LineLoop.cpp.i
.PHONY : mesh_interface/sources/LineLoop.i

# target to preprocess a source file
mesh_interface/sources/LineLoop.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/LineLoop.cpp.i
.PHONY : mesh_interface/sources/LineLoop.cpp.i

mesh_interface/sources/LineLoop.s: mesh_interface/sources/LineLoop.cpp.s
.PHONY : mesh_interface/sources/LineLoop.s

# target to generate assembly for a file
mesh_interface/sources/LineLoop.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/LineLoop.cpp.s
.PHONY : mesh_interface/sources/LineLoop.cpp.s

mesh_interface/sources/MeshFactor.o: mesh_interface/sources/MeshFactor.cpp.o
.PHONY : mesh_interface/sources/MeshFactor.o

# target to build an object file
mesh_interface/sources/MeshFactor.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/MeshFactor.cpp.o
.PHONY : mesh_interface/sources/MeshFactor.cpp.o

mesh_interface/sources/MeshFactor.i: mesh_interface/sources/MeshFactor.cpp.i
.PHONY : mesh_interface/sources/MeshFactor.i

# target to preprocess a source file
mesh_interface/sources/MeshFactor.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/MeshFactor.cpp.i
.PHONY : mesh_interface/sources/MeshFactor.cpp.i

mesh_interface/sources/MeshFactor.s: mesh_interface/sources/MeshFactor.cpp.s
.PHONY : mesh_interface/sources/MeshFactor.s

# target to generate assembly for a file
mesh_interface/sources/MeshFactor.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/MeshFactor.cpp.s
.PHONY : mesh_interface/sources/MeshFactor.cpp.s

mesh_interface/sources/PlaneSurface.o: mesh_interface/sources/PlaneSurface.cpp.o
.PHONY : mesh_interface/sources/PlaneSurface.o

# target to build an object file
mesh_interface/sources/PlaneSurface.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/PlaneSurface.cpp.o
.PHONY : mesh_interface/sources/PlaneSurface.cpp.o

mesh_interface/sources/PlaneSurface.i: mesh_interface/sources/PlaneSurface.cpp.i
.PHONY : mesh_interface/sources/PlaneSurface.i

# target to preprocess a source file
mesh_interface/sources/PlaneSurface.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/PlaneSurface.cpp.i
.PHONY : mesh_interface/sources/PlaneSurface.cpp.i

mesh_interface/sources/PlaneSurface.s: mesh_interface/sources/PlaneSurface.cpp.s
.PHONY : mesh_interface/sources/PlaneSurface.s

# target to generate assembly for a file
mesh_interface/sources/PlaneSurface.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/PlaneSurface.cpp.s
.PHONY : mesh_interface/sources/PlaneSurface.cpp.s

mesh_interface/sources/Point.o: mesh_interface/sources/Point.cpp.o
.PHONY : mesh_interface/sources/Point.o

# target to build an object file
mesh_interface/sources/Point.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Point.cpp.o
.PHONY : mesh_interface/sources/Point.cpp.o

mesh_interface/sources/Point.i: mesh_interface/sources/Point.cpp.i
.PHONY : mesh_interface/sources/Point.i

# target to preprocess a source file
mesh_interface/sources/Point.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Point.cpp.i
.PHONY : mesh_interface/sources/Point.cpp.i

mesh_interface/sources/Point.s: mesh_interface/sources/Point.cpp.s
.PHONY : mesh_interface/sources/Point.s

# target to generate assembly for a file
mesh_interface/sources/Point.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Point.cpp.s
.PHONY : mesh_interface/sources/Point.cpp.s

mesh_interface/sources/Surface.o: mesh_interface/sources/Surface.cpp.o
.PHONY : mesh_interface/sources/Surface.o

# target to build an object file
mesh_interface/sources/Surface.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Surface.cpp.o
.PHONY : mesh_interface/sources/Surface.cpp.o

mesh_interface/sources/Surface.i: mesh_interface/sources/Surface.cpp.i
.PHONY : mesh_interface/sources/Surface.i

# target to preprocess a source file
mesh_interface/sources/Surface.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Surface.cpp.i
.PHONY : mesh_interface/sources/Surface.cpp.i

mesh_interface/sources/Surface.s: mesh_interface/sources/Surface.cpp.s
.PHONY : mesh_interface/sources/Surface.s

# target to generate assembly for a file
mesh_interface/sources/Surface.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/mesh_interface/sources/Surface.cpp.s
.PHONY : mesh_interface/sources/Surface.cpp.s

solid/sources/Element.o: solid/sources/Element.cpp.o
.PHONY : solid/sources/Element.o

# target to build an object file
solid/sources/Element.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Element.cpp.o
.PHONY : solid/sources/Element.cpp.o

solid/sources/Element.i: solid/sources/Element.cpp.i
.PHONY : solid/sources/Element.i

# target to preprocess a source file
solid/sources/Element.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Element.cpp.i
.PHONY : solid/sources/Element.cpp.i

solid/sources/Element.s: solid/sources/Element.cpp.s
.PHONY : solid/sources/Element.s

# target to generate assembly for a file
solid/sources/Element.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Element.cpp.s
.PHONY : solid/sources/Element.cpp.s

solid/sources/Node.o: solid/sources/Node.cpp.o
.PHONY : solid/sources/Node.o

# target to build an object file
solid/sources/Node.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Node.cpp.o
.PHONY : solid/sources/Node.cpp.o

solid/sources/Node.i: solid/sources/Node.cpp.i
.PHONY : solid/sources/Node.i

# target to preprocess a source file
solid/sources/Node.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Node.cpp.i
.PHONY : solid/sources/Node.cpp.i

solid/sources/Node.s: solid/sources/Node.cpp.s
.PHONY : solid/sources/Node.s

# target to generate assembly for a file
solid/sources/Node.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Node.cpp.s
.PHONY : solid/sources/Node.cpp.s

solid/sources/Solid.o: solid/sources/Solid.cpp.o
.PHONY : solid/sources/Solid.o

# target to build an object file
solid/sources/Solid.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Solid.cpp.o
.PHONY : solid/sources/Solid.cpp.o

solid/sources/Solid.i: solid/sources/Solid.cpp.i
.PHONY : solid/sources/Solid.i

# target to preprocess a source file
solid/sources/Solid.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Solid.cpp.i
.PHONY : solid/sources/Solid.cpp.i

solid/sources/Solid.s: solid/sources/Solid.cpp.s
.PHONY : solid/sources/Solid.s

# target to generate assembly for a file
solid/sources/Solid.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/solid_library.dir/build.make CMakeFiles/solid_library.dir/solid/sources/Solid.cpp.s
.PHONY : solid/sources/Solid.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... s"
	@echo "... solid_library"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... mesh_interface/sources/Geometry.o"
	@echo "... mesh_interface/sources/Geometry.i"
	@echo "... mesh_interface/sources/Geometry.s"
	@echo "... mesh_interface/sources/Inclusion.o"
	@echo "... mesh_interface/sources/Inclusion.i"
	@echo "... mesh_interface/sources/Inclusion.s"
	@echo "... mesh_interface/sources/Line.o"
	@echo "... mesh_interface/sources/Line.i"
	@echo "... mesh_interface/sources/Line.s"
	@echo "... mesh_interface/sources/LineLoop.o"
	@echo "... mesh_interface/sources/LineLoop.i"
	@echo "... mesh_interface/sources/LineLoop.s"
	@echo "... mesh_interface/sources/MeshFactor.o"
	@echo "... mesh_interface/sources/MeshFactor.i"
	@echo "... mesh_interface/sources/MeshFactor.s"
	@echo "... mesh_interface/sources/PlaneSurface.o"
	@echo "... mesh_interface/sources/PlaneSurface.i"
	@echo "... mesh_interface/sources/PlaneSurface.s"
	@echo "... mesh_interface/sources/Point.o"
	@echo "... mesh_interface/sources/Point.i"
	@echo "... mesh_interface/sources/Point.s"
	@echo "... mesh_interface/sources/Surface.o"
	@echo "... mesh_interface/sources/Surface.i"
	@echo "... mesh_interface/sources/Surface.s"
	@echo "... solid/sources/Element.o"
	@echo "... solid/sources/Element.i"
	@echo "... solid/sources/Element.s"
	@echo "... solid/sources/Node.o"
	@echo "... solid/sources/Node.i"
	@echo "... solid/sources/Node.s"
	@echo "... solid/sources/Solid.o"
	@echo "... solid/sources/Solid.i"
	@echo "... solid/sources/Solid.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

