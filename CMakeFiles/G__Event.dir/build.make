# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/peter/Alpha-g/PG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/peter/Alpha-g/PG

# Utility rule file for G__Event.

# Include the progress variables for this target.
include CMakeFiles/G__Event.dir/progress.make

CMakeFiles/G__Event: G__Event.cxx
CMakeFiles/G__Event: libEvent_rdict.pcm
CMakeFiles/G__Event: libEvent.rootmap


G__Event.cxx: /home/peter/rootBuild/include/TVector.h
G__Event.cxx: /home/peter/rootBuild/include/TVector.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/peter/Alpha-g/PG/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating G__Event.cxx, libEvent_rdict.pcm, libEvent.rootmap"
	/usr/bin/cmake -E env LD_LIBRARY_PATH=/home/peter/rootBuild/lib:/home/peter/rootBuild/lib /home/peter/rootBuild/bin/rootcling -v2 -f G__Event.cxx -s /home/peter/Alpha-g/PG/libEvent.so -rml libEvent.so -rmf /home/peter/Alpha-g/PG/libEvent.rootmap -I/home/peter/rootBuild/include -I/home/peter/Alpha-g/PG TVector.h

libEvent_rdict.pcm: G__Event.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libEvent_rdict.pcm

libEvent.rootmap: G__Event.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libEvent.rootmap

G__Event: CMakeFiles/G__Event
G__Event: G__Event.cxx
G__Event: libEvent_rdict.pcm
G__Event: libEvent.rootmap
G__Event: CMakeFiles/G__Event.dir/build.make

.PHONY : G__Event

# Rule to build all files generated by this target.
CMakeFiles/G__Event.dir/build: G__Event

.PHONY : CMakeFiles/G__Event.dir/build

CMakeFiles/G__Event.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/G__Event.dir/cmake_clean.cmake
.PHONY : CMakeFiles/G__Event.dir/clean

CMakeFiles/G__Event.dir/depend:
	cd /home/peter/Alpha-g/PG && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/peter/Alpha-g/PG /home/peter/Alpha-g/PG /home/peter/Alpha-g/PG /home/peter/Alpha-g/PG /home/peter/Alpha-g/PG/CMakeFiles/G__Event.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/G__Event.dir/depend

