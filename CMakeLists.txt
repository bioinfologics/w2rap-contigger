cmake_minimum_required (VERSION 2.8)

if (POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif (POLICY CMP0074)

IF (MY_COMPILER_OPTION STREQUAL "Intel")
  FIND_PROGRAM(ICPC_PROGRAM icpc ...)
  IF(ICPC_PROGRAM)
    MESSAGE(STATUS "** Compiling with Intel settings **")
    IF(ENV{CXX})
      MESSAGE(WARNING "Overwriting CXX envvar")
    ENDIF()
    SET(ENV{CXX} "${ICPC_PROGRAM}")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ELSE()
    MESSAGE(FATAL_ERROR "Intel compiler not found")
  ENDIF()
ENDIF ()

project (w2rap-contigger)

# Set a default build type if none was specified
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Get the remote url
execute_process(
        COMMAND git remote get-url origin
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_ORIGIN_URL
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Fill the version variables in case the compilation isn't done from a git repo
if("${GIT_ORIGIN_URL}" STREQUAL "")
    add_definitions("-DGIT_ORIGIN_URL=\"not set\"")
else()
    add_definitions("-DGIT_ORIGIN_URL=\"${GIT_ORIGIN_URL}\"")
endif()
if("${GIT_COMMIT_HASH}" STREQUAL "")
add_definitions("-DGIT_COMMIT_HASH=\"experimental\"")
else()
add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
endif()
if("${GIT_BRANCH}" STREQUAL "")
add_definitions("-DGIT_COMMIT_HASH=\"not from git\"")
else()
add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")
endif()

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(DEFINED PROFILE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
endif()

if(DEFINED COVERAGE)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

## Compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fpermissive -w -Wno-deprecated -std=c++11 -gdwarf")

## To move the bins to the final location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

set(MALLOC_LIBRARY "" CACHE PATH "Malloc library path")

if (MALLOC_LIBRARY)
    ## MALLOC_LIBRARY, path to libtbbmalloc.so
    get_filename_component(LIBSO_PATH ${MALLOC_LIBRARY} DIRECTORY)
    get_filename_component(LIBSO_NAME ${MALLOC_LIBRARY} NAME)
    ## Link directories
    link_directories(${LIBSO_PATH})
endif()

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

find_package(BZip2)
if (BZIP2_FOUND)
    add_definitions(-DHAVE_BZIP2)
    include_directories(${BZIP2_INCLUDE_DIRS})
endif()

## Include directories
include_directories("src/")

## Mac compilation issue
include_directories("/usr/include/malloc")

## Set sources variables
add_library(base_libs OBJECT
        src/Charvector.cc
        src/Equiv.cc
        src/FastIfstream.cc
        src/Fastavector.cc
        src/IntPairVec.cc
        src/Intvector.cc
        src/MemberOf.cc
        src/PackAlign.cc
        src/PairsManager.cc
        src/ParseRange.cc
        src/ParseSet.cc
        src/Qualvector.cc
        src/TokenizeString.cc
        src/Vec.cc
        src/dna/Bases.cc
        src/efasta/EfastaTools.cc
        src/feudal/BaseVec.cc
        src/feudal/BinaryStream.cc
        src/feudal/BitVec.cc
        src/feudal/CharString.cc
        src/feudal/FeudalControlBlock.cc
        src/feudal/FeudalFileReader.cc
        src/feudal/FeudalFileWriter.cc
        src/feudal/FeudalTools.cc
        src/feudal/FieldVec.cc
        src/feudal/Generic.cc
        src/feudal/Mempool.cc
        src/feudal/Oob.cc
        src/feudal/PQVec.cc
        src/graph/Digraph.cc
        src/graph/FindCells.cc
        src/kmers/BigKPather.cc
        src/kmers/KMerContext.cc
        src/kmers/KmerRecord.cc
        src/kmers/ReadPather.cc
        src/lookup/SAM.cc
        src/math/Functions.cc
        src/math/HoInterval.cc
        src/math/Matrix.cc
        src/math/Permutation.cc
        src/math/PowerOf2.cc
        src/paths/FindErrorsCore.cc
        src/paths/HyperBasevector.cc
        src/paths/HyperEfasta.cc
        src/paths/KmerBaseBroker.cc
        src/paths/KmerPath.cc
        src/paths/KmerPathInterval.cc
        src/paths/Sepdev.cc
        src/paths/UnipathScaffold.cc
        src/paths/long/Correct1Pre.cc
        src/paths/long/CorrectByStack.cc
        src/paths/long/CorrectPairs1.cc
        src/paths/long/DiscovarTools.cc
        src/paths/long/EMEC3.cc
        src/paths/long/FillPairs.cc
        src/paths/long/FriendAligns.cc
        src/paths/long/HBVFromEdges.cc
        src/paths/long/LargeKDispatcher.cc
        src/paths/long/LoadCorrectCore.cc
        src/paths/long/LongProtoTools.cc
        src/paths/long/LongReadsToPaths.cc
        src/paths/long/MakeAlignments.cc
        src/paths/long/PreCorrectAlt1.cc
        src/paths/long/PreCorrectOldNew.cc
        src/paths/long/ReadPathTools.cc
        src/paths/long/ReadStack.cc
        src/random/Bernoulli.cc
        src/random/RNGen.cc
        src/reporting/PerfStat.cc
        src/system/Assert.cc
        src/system/ErrNo.cc
        src/system/Exit.cc
        src/system/HostName.cc
        src/system/ProcBuf.cc
        src/system/SysConf.cc
        src/system/System.cc
        src/system/Thread.cc
        src/system/ThreadsafeIO.cc
        src/system/TraceVal.cc
        src/system/WorklistUtils.cc
        src/system/file/File.cc
        src/system/file/FileReader.cc
        src/system/file/FileWriter.cc
        src/system/file/TempFile.cc
        src/util/Logger.cc
        src/kmers/kmer_parcels/KmerParcelsBuilder.cc
        src/math/IntDistribution.cc
        src/pairwise_aligners/MaxMutmerFromMer.cc
        src/paths/MakeAlignsPathsParallelX.cc
        src/paths/ReadsToPathsCoreX.cc
        src/paths/RemodelGapTools.cc
        src/paths/UnibaseUtils.cc
        src/paths/Unipath.cc
        src/paths/long/CreateGenome.cc
        src/paths/long/KmerCount.cc
        src/paths/long/LongHyper.cc
        src/paths/long/SupportedHyperBasevector.cc
        src/paths/long/SupportedHyperBasevector3.cc
        src/paths/long/SupportedHyperBasevector4.cc
        src/paths/long/SupportedHyperBasevector5.cc
        src/paths/long/SupportedHyperBasevector6.cc
        src/paths/long/SupportedHyperBasevector7.cc
        src/paths/long/large/CN1PeakFinder.cc
        src/paths/long/large/GapToyTools.cc
        src/paths/long/large/GapToyTools3.cc
        src/paths/long/large/Lines.cc
        src/random/NormalDistribution.cc
        src/GFADump.cc
        src/util/TextTable.cc
        src/SMR.h
        src/KMerFreqFactory.h
        src/paths/long/ReadPath.cc
        src/util/OutputLog.h
        src/util/OutputLog.cc
        src/paths/long/KMerNodeFreq_s.cc
        )

add_library(specific_w2rap-contigger OBJECT
        src/paths/long/large/GapToyTools4.cc
        src/paths/long/large/GapToyTools5.cc
        src/kmers/MakeLookup.cc
        src/kmers/kmatch/KMatch.cc
        src/paths/long/OverlapReads.cc
        src/paths/long/KmerAlign.cc
        src/paths/long/Friends.cc
        src/math/Combinatorics.cc
        src/kmers/LongReadPather.cc
        src/feudal/QualNibbleVec.cc
        src/Basevector.cc
        src/VecString.cc
        src/BasevectorTools.cc
        src/Misc.cc
        src/paths/long/large/AssembleGaps.cc
        src/paths/long/large/GapToyTools2.cc
        src/paths/long/large/Unsat.cc
        src/paths/long/BuildReadQGraph.cc
        src/paths/long/ExtendReadPath.cc
        src/paths/long/KmerAlign.h
        src/paths/long/ShortKmerReadPather.cc
        src/paths/long/large/Clean200.cc
        src/paths/long/large/ExtractReads.cc
        src/paths/long/large/ReadNameLookup.cc
        src/paths/long/large/Repath.cc
        src/system/Crash.cc
        src/util/PeakFinder.h
        src/paths/long/large/FinalFiles.cc
        src/paths/long/large/MakeGaps.cc
        src/paths/long/large/Simplify.cc
        src/paths/long/large/ImprovePath.cc
        src/util/w2rap_timers.h
        src/paths/PathFinder.cc
        src/paths/PathFinder.h
        src/paths/long/large/OverlapValidator.cc
        src/paths/long/large/OverlapValidator.h
        src/paths/long/large/GraphImprover.cc
        src/paths/long/large/GraphImprover.h
        src/paths/long/large/ConsensusChecker.cc
        src/paths/long/large/ConsensusChecker.h
        )

## add executables and respective sources
add_executable(w2rap-contigger src/modules/w2rap-contigger.cc
        $<TARGET_OBJECTS:specific_w2rap-contigger>
        $<TARGET_OBJECTS:base_libs>
        )

add_executable(hbv2gfa src/modules/hbv2gfa.cc
        $<TARGET_OBJECTS:base_libs>
        )

##Zlib link
if (ZLIB_FOUND)
  set(ZLIB libz.so)
  target_link_libraries(w2rap-contigger ${ZLIB_LIBRARIES})
  target_link_libraries(hbv2gfa ${ZLIB_LIBRARIES})
endif()

#Have the malloc library linked at the end, for compatibility issues with gperftools/tcmalloc
if (DEFINED MALLOC_LIBRARY)
    ## Link libraries
    if (LIBSO_NAME MATCHES ".*tcmalloc*")
        target_link_libraries(w2rap-contigger ${LIBSO_NAME} profiler)
    else ()
        target_link_libraries(w2rap-contigger ${LIBSO_NAME})
    endif()
endif()
