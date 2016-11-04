cmake_minimum_required (VERSION 2.8)
project (w2rap-contigger)

## Compiler flags
SET(CMAKE_CXX_FLAGS "-g -fopenmp -fpermissive -w -Wno-deprecated -std=c++11 -Ofast -DNDEBUG")

#SET(CMAKE_CXX_LINK_EXECUTABLE "/usr/bin/ld")
add_definitions(${CMAKE_CXX_FLAGS})

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

#if (DEFINED ZLIB_PATH)
#  include_directories("${ZLIB_PATH}/include")
#  link_directories("${ZLIB_PATH}/lib")
#else()
#    set(ZLIB_PATH "/usr/lib" CACHE FILEPATH "ZLIB Library path")
#endif()

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
endif()

## Include directories
include_directories("src/")

## Mac compilation issue
include_directories("/usr/include/malloc")

## Set sources variables

add_library(base_libs OBJECT
        src/Alignment.cc
        src/Basevector.cc
        src/Charvector.cc
        src/Equiv.cc
        src/FastIfstream.cc
        src/Fastavector.cc
        src/FetchReads.cc
        src/Floatvector.cc
        src/IndexedAlignmentPlusVector.cc
        src/IntPairVec.cc
        src/Intvector.cc
        src/LinkTimestamp.cc
        src/MemberOf.cc
        src/PackAlign.cc
        src/PairsManager.cc
        src/ParseRange.cc
        src/ParseSet.cc
        src/PrintAlignment.cc
        src/Qualvector.cc
        src/ScoreAlignment.cc
        src/TokenizeString.cc
        src/Vec.cc
        src/VecAlignmentPlus.cc
        src/VecString.cc
        src/dna/Bases.cc
        src/efasta/AmbiguityScore.cc
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
        src/feudal/QualNibbleVec.cc
        src/graph/Digraph.cc
        src/graph/FindCells.cc
        src/kmers/BigKPather.cc
        src/kmers/KMerContext.cc
        src/kmers/KmerRecord.cc
        src/kmers/LongReadPather.cc
        src/kmers/ReadPather.cc
        src/lookup/SAM.cc
        src/math/Combinatorics.cc
        src/math/Functions.cc
        src/math/HoInterval.cc
        src/math/Matrix.cc
        src/math/Permutation.cc
        src/math/PowerOf2.cc
        src/pairwise_aligners/SmithWatBandedA.cc
        src/pairwise_aligners/SmithWatFree.cc
        src/paths/BigMapTools.cc
        src/paths/FindErrorsCore.cc
        src/paths/HyperBasevector.cc
        src/paths/HyperEfasta.cc
        src/paths/KmerBaseBroker.cc
        src/paths/KmerPath.cc
        src/paths/KmerPathInterval.cc
        src/paths/LongReadTools.cc
        src/paths/ProcessGap.cc
        src/paths/Sepdev.cc
        src/paths/UnipathScaffold.cc
        src/paths/Uniseq.cc
        src/paths/long/Correct1Pre.cc
        src/paths/long/CorrectByStack.cc
        src/paths/long/CorrectPairs1.cc
        src/paths/long/DiscovarTools.cc
        src/paths/long/EMEC3.cc
        src/paths/long/FillPairs.cc
        src/paths/long/FriendAligns.cc
        src/paths/long/Friends.cc
        src/paths/long/HBVFromEdges.cc
        src/paths/long/KmerAlign.cc
        src/paths/long/LargeKDispatcher.cc
        src/paths/long/LoadCorrectCore.cc
        src/paths/long/LongProtoTools.cc
        src/paths/long/LongReadsToPaths.cc
        src/paths/long/MakeAlignments.cc
        src/paths/long/OverlapReads.cc
        src/paths/long/PreCorrectAlt1.cc
        src/paths/long/PreCorrectOldNew.cc
        src/paths/long/ReadPathTools.cc
        src/paths/long/ReadStack.cc
        src/paths/long/ultra/ConsensusScoreModel.cc
        src/paths/long/ultra/FounderAlignment.cc
        src/paths/long/ultra/GetFriendsAndAlignsInitial.cc
        src/paths/long/ultra/MakeBlocks.cc
        src/paths/long/ultra/MultipleAligner.cc
        src/paths/long/ultra/Prefab.cc
        src/paths/long/ultra/ThreadedBlocks.cc
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
        src/system/gzstream/gzstream.C
        src/util/Logger.cc
        src/fastg/FastgGraph.cc
        src/graphics/BasicGraphics.cc
        src/kmers/kmer_parcels/KmerParcelsBuilder.cc
        src/kmers/MakeLookup.cc
        src/math/IntDistribution.cc
        src/pairwise_aligners/MaxMutmerFromMer.cc
        src/pairwise_aligners/SmithWatAffine.cc
        src/paths/MakeAlignsPathsParallelX.cc
        src/paths/ReadsToPathsCoreX.cc
        src/paths/RemodelGapTools.cc
        src/paths/UnibaseUtils.cc
        src/paths/Unipath.cc
        src/paths/long/AssessBestAlignCore.cc
        src/paths/long/CleanEfasta.cc
        src/paths/long/CreateGenome.cc
        src/paths/long/DigraphFromWords.cc
        #src/paths/long/EvalAssembly.cc
        src/paths/long/EvalByReads.cc
        src/paths/long/KmerCount.cc
        src/paths/long/LongHyper.cc
        src/paths/long/ReadOriginTracker.cc
        src/paths/long/RefTrace.cc
        src/paths/long/RefTraceControl.cc
        src/paths/long/RefTraceTools.cc
        src/paths/long/SupportedHyperBasevector.cc
        src/paths/long/SupportedHyperBasevector3.cc
        src/paths/long/SupportedHyperBasevector4.cc
        src/paths/long/SupportedHyperBasevector5.cc
        src/paths/long/SupportedHyperBasevector6.cc
        src/paths/long/SupportedHyperBasevector7.cc
        src/paths/long/VariantCallTools.cc
        src/paths/long/VariantFilters.cc
        src/paths/long/VariantPostProcess.cc
        src/paths/long/VariantReadSupport.cc
        src/paths/long/Variants.cc
        src/paths/long/large/CN1PeakFinder.cc
        src/paths/long/large/GapToyTools.cc
        src/paths/long/large/GapToyTools3.cc
        src/paths/long/large/GapToyTools4.cc
        src/paths/long/large/GapToyTools5.cc
        src/paths/long/large/Lines.cc
        src/paths/simulation/VCF.cc
        src/random/NormalDistribution.cc
        src/util/TextTable.cc
        src/util/w2rap_timers.h
        src/paths/long/ReadPath.cc)

add_library(hb_base_libs OBJECT
        src/Alignment.cc
        src/Charvector.cc
        src/Equiv.cc
        src/FastIfstream.cc
        src/Fastavector.cc
        src/Floatvector.cc
        src/IndexedAlignmentPlusVector.cc
        src/IntPairVec.cc
        src/Intvector.cc
        src/LinkTimestamp.cc
        src/MemberOf.cc
        src/PackAlign.cc
        src/PairsManager.cc
        src/ParseRange.cc
        src/ParseSet.cc
        src/PrintAlignment.cc
        src/Qualvector.cc
        src/ScoreAlignment.cc
        src/TokenizeString.cc
        src/Vec.cc
        src/VecAlignmentPlus.cc
        src/dna/Bases.cc
        src/efasta/AmbiguityScore.cc
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
        src/pairwise_aligners/SmithWatBandedA.cc
        src/pairwise_aligners/SmithWatFree.cc
        src/paths/BigMapTools.cc
        src/paths/FindErrorsCore.cc
        src/paths/HyperBasevector.cc
        src/paths/HyperEfasta.cc
        src/paths/KmerBaseBroker.cc
        src/paths/KmerPath.cc
        src/paths/KmerPathInterval.cc
        src/paths/LongReadTools.cc
        src/paths/ProcessGap.cc
        src/paths/Sepdev.cc
        src/paths/UnipathScaffold.cc
        src/paths/Uniseq.cc
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
        src/fastg/FastgGraph.cc
        src/graphics/BasicGraphics.cc
        src/kmers/kmer_parcels/KmerParcelsBuilder.cc
        src/math/IntDistribution.cc
        src/pairwise_aligners/MaxMutmerFromMer.cc
        src/pairwise_aligners/SmithWatAffine.cc
        src/paths/MakeAlignsPathsParallelX.cc
        src/paths/ReadsToPathsCoreX.cc
        src/paths/RemodelGapTools.cc
        src/paths/UnibaseUtils.cc
        src/paths/Unipath.cc
        src/paths/long/AssessBestAlignCore.cc
        src/paths/long/CleanEfasta.cc
        src/paths/long/CreateGenome.cc
        src/paths/long/DigraphFromWords.cc
        src/paths/long/EvalByReads.cc
        src/paths/long/KmerCount.cc
        src/paths/long/LongHyper.cc
        src/paths/long/ReadOriginTracker.cc
        src/paths/long/RefTrace.cc
        src/paths/long/RefTraceControl.cc
        src/paths/long/RefTraceTools.cc
        src/paths/long/SupportedHyperBasevector.cc
        src/paths/long/SupportedHyperBasevector3.cc
        src/paths/long/SupportedHyperBasevector4.cc
        src/paths/long/SupportedHyperBasevector5.cc
        src/paths/long/SupportedHyperBasevector6.cc
        src/paths/long/SupportedHyperBasevector7.cc
        src/paths/long/VariantCallTools.cc
        src/paths/long/VariantFilters.cc
        src/paths/long/VariantPostProcess.cc
        src/paths/long/VariantReadSupport.cc
        src/paths/long/Variants.cc
        src/paths/long/large/CN1PeakFinder.cc
        src/paths/long/large/GapToyTools.cc
        src/paths/long/large/GapToyTools3.cc
        src/paths/long/large/Lines.cc
        src/paths/simulation/VCF.cc
        src/random/NormalDistribution.cc
        src/util/TextTable.cc
        src/GFADump.cc
        src/paths/long/ReadPath.cc)

add_library(specific_w2rap-contigger OBJECT
        src/BasevectorTools.cc
        src/CompressedSequence.cc
        src/FastaConverter.cc
        src/FastaFileset.cc
        src/FastaFilestream.cc
        src/FastaFilestreamPreview.cc
        src/FastaNameParser.cc
        src/FastaVerifier.cc
        src/Misc.cc
        src/bam/ReadBAM.cc
        src/pairwise_aligners/SmithWatBanded.cc
        src/paths/MuxGraph.cc
        src/paths/OffsetTracker.cc
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
        src/paths/long/DisplayTools.cc
        src/paths/long/large/FinalFiles.cc
        src/paths/long/large/MakeGaps.cc
        src/paths/long/large/Simplify.cc
        src/paths/long/large/ImprovePath.cc
        src/GFADump.cc
        src/paths/PathFinder.cc
        src/paths/PathFinder.h)

## add executables and respective sources
add_executable(w2rap-contigger src/modules/w2rap-contigger.cc
        $<TARGET_OBJECTS:specific_w2rap-contigger>
        $<TARGET_OBJECTS:base_libs>
        )

add_executable(hbv2gfa src/modules/hbv2gfa.cc
        $<TARGET_OBJECTS:hb_base_libs>
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

