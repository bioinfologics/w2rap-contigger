cmake_minimum_required (VERSION 2.8)
project (A-Discovar)

## Compiler flags
SET(CMAKE_CXX_FLAGS "-g -fopenmp -fpermissive -w -Wno-deprecated -std=c++11 -Ofast -lz")
#SET(CMAKE_CXX_LINK_EXECUTABLE "/usr/bin/ld")
add_definitions(${CMAKE_CXX_FLAGS})

## To move the bins to the final location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

## Include directories
include_directories("src/")

## Link directories
link_directories("/tgac/software/production/intel_studio/2015.2/x86_64/composerxe/tbb/lib/intel64/gcc4.4")


## Set sources variables
SET(SOURCE_FILES_UNIPATH_01 src/modules/unipaths_01.cc src/paths/long/large/GapToyTools6.cc src/paths/long/DiscovarTools.cc src/system/System.cc
        src/feudal/CharString.cc src/paths/long/large/Samples.cc src/paths/long/large/ExtractReads.cc src/system/Assert.cc
        src/system/file/FileReader.cc src/system/file/FileWriter.cc src/Charvector.cc src/Basevector.cc src/feudal/PQVec.cc
        src/system/file/File.cc src/system/Exit.cc src/system/SysConf.cc src/FastIfstream.cc src/paths/long/fosmid/Fosmids.cc
        src/ParseSet.cc src/lookup/SAM.cc src/Fastavector.cc src/util/Logger.cc src/system/ErrNo.cc src/system/file/TempFile.cc
        src/feudal/Mempool.cc src/TokenizeString.cc src/feudal/BinaryStream.cc src/feudal/FeudalTools.cc src/bam/ReadBAM.cc
        src/feudal/FieldVec.cc src/paths/long/large/ReadNameLookup.cc src/Qualvector.cc src/feudal/BaseVec.cc src/feudal/BitVec.cc
        src/feudal/QualNibbleVec.cc src/FetchReads.cc src/paths/long/LoadCorrectCore.cc src/dna/Bases.cc src/random/Shuffle.cc
        src/PairsManager.cc src/feudal/FeudalFileReader.cc src/VecString.cc src/feudal/FeudalFileWriter.cc src/system/WorklistUtils.cc
        src/feudal/Generic.cc src/feudal/Oob.cc src/system/Thread.cc src/lookup/SAM2CRD.cc src/system/ProcBuf.cc src/math/HoInterval.cc
        src/random/RNGen.cc src/graph/Digraph.cc src/Equiv.cc src/feudal/FeudalControlBlock.cc src/Floatvector.cc src/paths/MergeReadSetsCore.cc
        src/lookup/LookAlign.cc src/paths/long/LongProtoTools.cc src/paths/long/PreCorrectAlt1.cc src/paths/long/PreCorrectOldNew.cc
        src/paths/long/FillPairs.cc src/paths/KmerPath.cc src/paths/long/Correct1Pre.cc src/paths/FindErrorsCore.cc src/efasta/EfastaTools.cc
        src/paths/long/LongReadsToPaths.cc src/paths/long/CorrectPairs1.cc src/kmers/KmerRecord.cc src/paths/long/Correct1.cc
        src/paths/HyperBasevector.cc src/random/MersenneTwister.cc src/Vec.cc src/system/ThreadsafeIO.cc src/math/Matrix.cc
        src/PackAlign.cc src/ScoreAlignment.cc src/PrintAlignment.cc src/BlockAlign.cc src/system/ParsedArgs.cc src/paths/long/ReadStack.cc
        src/IntPairVec.cc src/Intvector.cc src/paths/OffsetTracker.cc src/kmers/LongReadPather.cc src/reporting/PerfStat.cc
        src/ParseRange.cc src/math/PowerOf2.cc src/system/TraceVal.cc src/pairwise_aligners/ClusterAligner.cc src/paths/KmerPathInterval.cc
        src/paths/long/ultra/ThreadedBlocks.cc src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/paths/long/ultra/ConsensusScoreModel.cc
        src/math/Permutation.cc src/paths/long/ultra/Prefab.cc src/paths/long/ultra/MakeBlocks.cc src/paths/long/Friends.cc src/math/Combinatorics.cc
        src/paths/long/OverlapReads.cc src/random/Bernoulli.cc src/paths/long/ultra/FounderAlignment.cc src/paths/long/KmerAlign.h
        src/paths/KmerBaseBroker.cc src/paths/long/ShortKmerReadPather.cc src/MemberOf.cc src/paths/long/ultra/MultipleAligner.cc
        src/paths/long/FriendAligns.cc src/paths/long/ReadPathTools.cc src/paths/long/ExtendReadPath.cc src/paths/long/CorrectByStack.cc
        src/kmers/BigKPather.cc src/paths/long/HBVFromEdges.cc src/paths/long/ReadPath.cc src/paths/long/KmerAlign.cc
        src/math/Functions.cc src/paths/BigMapTools.cc src/paths/UnipathScaffold.cc src/paths/long/EMEC3.cc src/paths/HyperEfasta.cc
        src/kmers/KMerContext.cc src/paths/long/LargeKDispatcher.cc src/system/HostName.cc src/LinkTimestamp.cc src/system/UseGDB.cc
        src/paths/MuxGraph.cc src/paths/LongReadTools.cc src/pairwise_aligners/SmithWatBandedA.cc src/pairwise_aligners/SmithWatFree.cc
        src/paths/ProcessGap.cc src/paths/Uniseq.cc src/graph/FindCells.cc src/paths/Sepdev.cc src/Superb.cc src/kmers/ReadPather.cc
        src/efasta/AmbiguityScore.cc src/paths/long/MakeAlignments.cc src/AnnotatedContig.cc src/Misc.cc src/system/Crash.cc
        src/Alignment.cc src/VecAlignmentPlus.cc src/IndexedAlignmentPlusVector.cc)

SET(SOURCE_FILES_QGRAPH_02 src/modules/qgraph_02.cc src/paths/long/large/GapToyTools6.cc src/paths/long/DiscovarTools.cc
        src/system/System.cc src/feudal/CharString.cc src/paths/HyperBasevector.cc src/paths/long/large/GapToyTools.cc
        src/paths/long/large/Repath.cc src/system/Assert.cc src/feudal/BinaryStream.cc src/system/file/FileReader.cc
        src/system/file/FileWriter.cc src/feudal/Mempool.cc src/feudal/FieldVec.cc src/Charvector.cc
        src/paths/long/ReadPath.cc src/Basevector.cc src/feudal/PQVec.cc src/IntPairVec.cc src/BasevectorTools.cc
        src/Intvector.cc src/feudal/FeudalFileReader.cc src/system/file/File.cc src/feudal/FeudalFileWriter.cc
        src/system/SysConf.cc src/FastIfstream.cc src/paths/long/fosmid/Fosmids.cc src/ParseSet.cc src/lookup/SAM.cc
        src/Fastavector.cc src/util/Logger.cc src/system/Exit.cc src/system/ErrNo.cc src/system/file/TempFile.cc
        src/TokenizeString.cc src/feudal/BaseVec.cc src/paths/KmerPath.cc src/Equiv.cc src/graph/Digraph.cc src/dna/Bases.cc
        src/paths/KmerBaseBroker.cc src/VecString.cc src/paths/long/CreateGenome.cc src/paths/long/EvalAssembly.cc
        src/paths/ReadsToPathsCoreX.cc src/paths/Unipath.cc src/efasta/EfastaTools.cc src/paths/long/LongHyper.cc
        src/paths/long/SupportedHyperBasevector6.cc src/PairsManager.cc src/paths/long/LoadCorrectCore.cc
        src/paths/long/SupportedHyperBasevector3.cc src/paths/long/large/LocalLayout.cc src/lookup/LookAlign.cc
        src/paths/long/large/GapToyTools3.cc src/math/PowerOf2.cc src/system/WorklistUtils.cc src/kmers/KmerRecord.cc
        src/paths/KmerPathInterval.cc src/FastaFileset.cc src/kmers/KMerContext.cc src/feudal/Generic.cc
        src/kmers/ReadPather.cc src/paths/long/HBVFromEdges.cc src/system/ProcBuf.cc src/feudal/FieldVec.cc
        src/paths/long/BuildReadQGraph.cc src/feudal/FeudalControlBlock.cc src/feudal/Oob.cc src/system/Thread.cc
        src/paths/long/ShortKmerReadPather.cc src/paths/long/ExtendReadPath.cc src/system/ThreadsafeIO.cc
        src/FastaFilestream.cc src/FastaFilestreamPreview.cc src/FastaVerifier.cc src/FastaConverter.cc
        src/FastaNameParser.cc src/BlockAlign.cc src/PrintAlignment.cc src/paths/long/SupportedHyperBasevector.cc
        src/paths/long/Correct1.cc src/paths/long/FriendAligns.cc src/math/Matrix.cc src/graphics/BasicGraphics.cc
        src/paths/long/SupportedHyperBasevector4.cc src/paths/long/ReadStack.cc src/paths/HyperEfasta.cc
        src/system/ParsedArgs.cc src/feudal/BitVec.cc src/paths/MakeAlignsPathsParallelX.cc src/ParseRange.cc src/MemberOf.cc
        src/system/TraceVal.cc src/FetchReads.cc src/pairwise_aligners/MaxMutmerFromMer.cc
        src/kmers/kmer_parcels/KmerParcelsBuilder.cc src/pairwise_aligners/SmithWatAffine.cc src/Floatvector.cc
        src/paths/long/EvalByReads.cc src/PackAlign.cc src/reporting/PerfStat.cc src/paths/long/LongProtoTools.cc
        src/paths/MergeReadSetsCore.cc src/paths/long/KmerCount.cc src/paths/long/Friends.cc src/random/RNGen.cc
        src/VecString.cc src/paths/long/ultra/ThreadedBlocks.cc src/paths/long/ultra/ConsensusScoreModel.cc
        src/paths/BigMapTools.cc src/paths/long/ultra/MultipleAligner.cc src/paths/LongReadTools.cc
        src/paths/UnipathScaffold.cc src/pairwise_aligners/SmithWatFree.cc src/paths/long/ultra/Prefab.cc
        src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/Vec.cc src/paths/long/KmerAlign.cc src/paths/long/OverlapReads.cc
        src/paths/long/SupportedHyperBasevector5.cc src/math/HoInterval.cc src/Alignment.cc src/random/Bernoulli.cc
        src/math/Combinatorics.cc src/paths/long/ultra/FounderAlignment.cc src/pairwise_aligners/SmithWatBandedA.cc
        src/math/Permutation.cc src/paths/long/RefTraceTools.cc src/IndexedAlignmentPlusVector.cc src/paths/ProcessGap.cc
        src/paths/long/large/Lines.cc src/paths/Uniseq.cc src/kmers/BigKPather.cc src/paths/long/SupportedHyperBasevector7.cc
        src/paths/Sepdev.cc src/math/IntDistribution.cc src/math/Functions.cc src/paths/UnibaseUtils.cc src/CompressedSequence.cc
        src/paths/long/fosmid/FosmidPool.cc src/feudal/FeudalTools.cc src/paths/long/LongReadsToPaths.cc src/random/Shuffle.cc
        src/paths/long/RefTraceControl.cc src/paths/long/PreCorrectAlt1.cc src/lookup/SAM2CRD.cc src/paths/long/FillPairs.cc
        src/kmers/LongReadPather.cc src/random/MersenneTwister.cc src/fastg/FastgGraph.cc src/paths/long/Correct1Pre.cc
        src/paths/long/DigraphFromWords.cc src/paths/long/PreCorrectOldNew.cc src/paths/long/CleanEfasta.cc src/Superb.cc
        src/util/PeakFinder.h src/graph/FindCells.cc src/efasta/AmbiguityScore.cc src/VecAlignmentPlus.cc src/system/UseGDB.cc
        src/paths/long/PreCorrectOldNew.cc src/paths/long/Variants.cc src/paths/long/large/CN1PeakFinder.cc
        src/paths/long/LargeKDispatcher.cc src/LinkTimestamp.cc src/paths/long/EMEC3.cc src/pairwise_aligners/ClusterAligner.cc
        src/paths/long/RefTrace.cc src/ScoreAlignment.cc src/paths/long/ultra/MakeBlocks.cc src/paths/long/VariantPostProcess.cc
        src/system/HostName.cc src/paths/long/CorrectPairs1.cc src/paths/long/CorrectByStack.cc src/paths/long/AssessBestAlignCore.cc
        src/paths/long/MakeAlignments.cc src/paths/long/ReadPath.cc src/paths/long/VariantReadSupport.cc src/paths/RemodelGapTools.cc
        src/paths/long/ReadPathTools.cc src/paths/long/VariantCallTools.cc src/Qualvector.cc src/paths/long/VariantFilters.cc
        src/paths/FindErrorsCore.cc src/paths/RemodelGapTools.cc src/paths/long/ReadOriginTracker.cc src/Qualvector.cc
        src/random/NormalDistribution.cc src/paths/simulation/VCF.cc src/feudal/QualNibbleVec.cc src/util/TextTable.cc
        src/AnnotatedContig.cc src/system/Crash.cc src/Misc.cc)

SET(SOURCE_FILES_CLEAN_03 src/modules/clean_03.cc src/paths/long/DiscovarTools.cc src/system/System.cc src/paths/long/large/GapToyTools6.cc
        src/ParseSet.cc src/feudal/CharString.cc src/system/SysConf.cc src/paths/long/fosmid/Fosmids.cc src/TokenizeString.cc src/Fastavector.cc
        src/system/Thread.cc src/FastaFilestream.cc src/Charvector.cc src/system/file/FileReader.cc src/system/file/File.cc
        src/system/ErrNo.cc src/graph/Digraph.cc src/Basevector.cc src/FastaFilestream.cc src/paths/HyperBasevector.cc src/paths/KmerPath.cc
        src/math/HoInterval.cc src/paths/KmerPathInterval.cc src/paths/long/large/Clean200.cc src/Intvector.cc src/feudal/FieldVec.cc
        src/Floatvector.cc src/system/file/FileWriter.cc src/math/Matrix.cc src/dna/Bases.cc src/feudal/Mempool.cc src/Equiv.cc src/system/file/TempFile.cc
        src/feudal/BaseVec.cc src/FastaFilestreamPreview.cc src/Vec.cc src/feudal/Generic.cc src/paths/long/large/GapToyTools3.cc
        src/paths/long/large/GapToyTools.cc src/paths/long/CreateGenome.cc src/feudal/FeudalFileWriter.cc src/feudal/FeudalFileReader.cc
        src/FetchReads.cc src/graphics/BasicGraphics.cc src/paths/long/large/LocalLayout.cc src/system/ThreadsafeIO.cc src/paths/KmerBaseBroker.cc
        src/paths/long/ReadStack.cc src/IntPairVec.cc src/PairsManager.cc src/feudal/FeudalControlBlock.cc src/paths/long/SupportedHyperBasevector3.cc
        src/paths/long/ReadPath.cc src/efasta/EfastaTools.cc src/paths/long/FriendAligns.cc src/feudal/BitVec.cc src/random/Bernoulli.cc
        src/paths/long/EvalAssembly.cc src/lookup/LookAlign.cc src/feudal/PQVec.cc src/PackAlign.cc src/math/PowerOf2.cc src/random/RNGen.cc
        src/system/WorklistUtils.cc src/kmers/KMerContext.cc src/feudal/BinaryStream.cc src/paths/long/RefTraceTools.cc src/paths/long/MakeAlignments.cc
        src/ScoreAlignment.cc src/paths/long/LoadCorrectCore.cc src/paths/RemodelGapTools.cc src/VecString.cc src/paths/long/Correct1Pre.cc
        src/paths/MergeReadSetsCore.cc src/paths/long/LongProtoTools.cc src/kmers/KmerRecord.cc src/system/Assert.cc
        src/kmers/ReadPather.cc src/util/Logger.cc src/system/ProcBuf.cc src/reporting/PerfStat.cc src/system/ParsedArgs.cc src/math/IntDistribution.cc
        src/paths/long/ultra/ThreadedBlocks.cc src/Qualvector.cc src/FastaConverter.cc src/feudal/Oob.cc src/pairwise_aligners/SmithWatBanded.cc
        src/paths/long/SupportedHyperBasevector4.cc src/lookup/SAM.cc src/paths/long/EvalByReads.cc src/paths/HyperEfasta.cc
        src/MemberOf.cc src/pairwise_aligners/ClusterAligner.cc src/paths/long/PreCorrectAlt1.cc src/feudal/FeudalTools.cc
        src/paths/LongReadTools.cc src/paths/long/LongReadsToPaths.cc src/pairwise_aligners/SmithWatAffine.cc src/ParseRange.cc
        src/paths/ReadsToPathsCoreX.cc src/PrintAlignment.cc src/paths/long/ReadPathTools.cc
        src/paths/FindErrorsCore.cc src/kmers/BigKPather.cc src/paths/BigMapTools.cc src/paths/long/fosmid/FosmidPool.cc src/paths/long/RefTrace.cc
        src/pairwise_aligners/SmithWatFree.cc src/paths/long/HBVFromEdges.cc src/paths/long/PreCorrectOldNew.cc src/paths/long/Variants.cc
        src/pairwise_aligners/SmithWatBandedA.cc src/paths/simulation/VCF.cc src/util/TextTable.cc src/paths/UnipathScaffold.cc
        src/feudal/QualNibbleVec.cc src/paths/long/ultra/MultipleAligner.cc src/kmers/LongReadPather.cc src/paths/long/LargeKDispatcher.cc
        src/paths/long/SupportedHyperBasevector6.cc src/FastIfstream.cc src/paths/long/ultra/MakeBlocks.cc src/paths/long/ultra/ConsensusScoreModel.cc
        src/paths/long/ReadOriginTracker.cc src/paths/long/ultra/ConsensusScoreModel.cc src/paths/long/VariantCallTools.cc
        src/CompressedSequence.cc src/paths/long/ReadOriginTracker.cc src/paths/long/RefTraceControl.cc src/system/TraceVal.cc
        src/FastaVerifier.cc src/math/Functions.cc src/lookup/SAM2CRD.cc src/paths/long/SupportedHyperBasevector7.cc src/paths/long/Friends.cc
        src/system/HostName.cc src/paths/long/Correct1.cc src/paths/long/AssessBestAlignCore.cc src/paths/long/CleanEfasta.cc
        src/graph/Digraph.cc src/paths/long/ultra/Prefab.cc src/paths/long/LongProtoTools.cc src/paths/long/VariantFilters.cc
        src/paths/MakeAlignsPathsParallelX.cc src/random/NormalDistribution.cc src/paths/long/VariantPostProcess.cc src/paths/long/CorrectPairs1.cc
        src/paths/long/CorrectByStack.cc src/paths/long/KmerAlign.cc src/paths/long/VariantReadSupport.cc src/pairwise_aligners/MaxMutmerFromMer.cc
        src/paths/long/KmerCount.cc src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/paths/long/large/Lines.cc src/paths/long/FillPairs.cc
        src/paths/long/LongHyper.cc src/random/Shuffle.cc src/graph/FindCells.cc src/random/MersenneTwister.cc src/math/Combinatorics.cc
        src/paths/UnibaseUtils.cc src/paths/long/ultra/FounderAlignment.cc src/fastg/FastgGraph.cc src/paths/long/OverlapReads.cc
        src/paths/long/SupportedHyperBasevector5.cc src/kmers/kmer_parcels/KmerParcelsBuilder.cc src/system/Exit.cc src/BlockAlign.cc
        src/math/Permutation.cc src/paths/long/EMEC3.cc src/paths/ProcessGap.cc src/paths/Unipath.cc src/efasta/AmbiguityScore.cc
        src/paths/Uniseq.cc src/efasta/AmbiguityScore.cc src/paths/long/large/CN1PeakFinder.cc src/paths/Sepdev.cc src/Alignment.cc
        src/Superb.cc src/IndexedAlignmentPlusVector.cc src/VecAlignmentPlus.cc src/paths/long/SupportedHyperBasevector.cc src/LinkTimestamp.cc
        src/system/UseGDB.cc src/paths/long/DigraphFromWords.cc)

SET(SOURCE_FILES_PATCHING_04 src/modules/patching_04.cc src/paths/HyperBasevector.cc src/Intvector.cc src/paths/long/ReadPath.cc src/VecString.cc
        src/paths/KmerPath.cc src/graph/Digraph.cc src/Charvector.cc src/Basevector.cc src/feudal/PQVec.cc src/feudal/FieldVec.cc
        src/paths/KmerBaseBroker.cc src/Equiv.cc src/feudal/FeudalFileReader.cc src/feudal/Mempool.cc src/system/file/FileWriter.cc
        src/system/file/FileReader.cc src/system/Assert.cc src/feudal/BaseVec.cc src/math/Matrix.cc src/system/ErrNo.cc src/paths/KmerPathInterval.cc
        src/feudal/CharString.cc src/feudal/FeudalControlBlock.cc src/system/file/File.cc src/dna/Bases.cc src/feudal/Generic.cc src/system/Exit.cc
        src/feudal/Oob.cc src/feudal/FeudalFileWriter.cc src/feudal/BinaryStream.cc src/math/Permutation.cc src/paths/long/large/GapToyTools3.cc
        src/system/System.cc src/paths/long/ReadStack.cc src/PairsManager.cc src/paths/long/DiscovarTools.cc src/paths/long/large/GapToyTools4.cc
        src/random/Bernoulli.cc src/paths/long/large/GapToyTools.cc src/paths/long/large/GapToyTools6.cc src/Vec.cc src/math/HoInterval.cc
        src/paths/long/large/AssembleGaps.cc src/ParseSet.cc src/paths/long/large/GapToyTools5.cc src/graphics/BasicGraphics.cc
        src/paths/long/FriendAligns.cc src/system/TraceVal.cc src/efasta/EfastaTools.cc src/feudal/BitVec.cc src/paths/long/SupportedHyperBasevector3.cc
        src/lookup/LookAlign.cc src/system/file/TempFile.cc src/math/PowerOf2.cc src/paths/long/CreateGenome.cc src/paths/long/large/GapToyTools2.cc
        src/paths/long/large/LocalLayout.cc src/paths/long/EvalAssembly.cc src/paths/long/fosmid/Fosmids.cc src/paths/long/LongProtoTools.cc
        src/paths/long/MakeAlignments.cc src/ScoreAlignment.cc src/system/SysConf.cc src/paths/long/LoadCorrectCore.cc src/paths/RemodelGapTools.cc
        src/paths/long/fosmid/Fosmids.cc src/paths/long/LongReadsToPaths.cc src/TokenizeString.cc src/PackAlign.cc src/paths/long/KmerCount.cc
        src/paths/HyperEfasta.cc src/paths/long/Friends.cc src/FetchReads.cc src/paths/long/ultra/MakeBlocks.cc src/paths/long/Correct1Pre.cc
        src/paths/MergeReadSetsCore.cc src/paths/long/RefTraceTools.cc src/paths/long/CorrectPairs1.cc src/paths/long/PreCorrectAlt1.cc
        src/paths/long/PreCorrectOldNew.cc src/FastIfstream.cc src/Fastavector.cc src/paths/ReadsToPathsCoreX.cc src/system/Thread.cc
        src/kmers/MakeLookup.cc src/paths/long/ultra/Prefab.cc src/pairwise_aligners/SmithWatAffine.cc src/PrintAlignment.cc
        src/paths/long/ReadPathTools.cc src/paths/long/OverlapReads.cc src/paths/long/ultra/ThreadedBlocks.cc src/math/IntDistribution.cc
        src/system/ParsedArgs.cc src/paths/long/ultra/ConsensusScoreModel.cc src/Floatvector.cc src/IntPairVec.cc src/system/ThreadsafeIO.cc
        src/reporting/PerfStat.cc src/system/ProcBuf.cc src/paths/long/ultra/MultipleAligner.cc src/util/Logger.cc src/kmers/ReadPather.cc
        src/random/RNGen.cc src/kmers/KmerRecord.cc src/kmers/LongReadPather.cc src/paths/long/LargeKDispatcher.cc src/paths/long/SupportedHyperBasevector4.cc
        src/paths/long/SupportedHyperBasevector6.cc src/pairwise_aligners/SmithWatBandedA.cc src/system/WorklistUtils.cc src/paths/long/ReadOriginTracker.cc
        src/BlockAlign.cc src/kmers/KMerContext.cc src/paths/long/EMEC3.cc src/random/Shuffle.cc src/paths/long/LongHyper.cc src/paths/long/FillPairs.cc
        src/paths/long/RefTrace.cc src/paths/long/RefTraceControl.cc src/paths/long/SupportedHyperBasevector5.cc
        src/paths/long/Correct1.cc src/paths/Unipath.cc src/lookup/SAM.cc src/lookup/SAM2CRD.cc src/paths/long/SupportedHyperBasevector.cc
        src/paths/LongReadTools.cc src/paths/long/large/Unsat.cc src/paths/long/large/GapToyTools3.cc src/math/Functions.cc src/paths/RemodelGapTools.cc
        src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/feudal/FeudalTools.cc src/paths/MakeAlignsPathsParallelX.cc
        src/paths/long/fosmid/FosmidPool.cc src/paths/BigMapTools.cc src/paths/UnipathScaffold.cc src/Alignment.cc src/kmers/kmer_parcels/KmerParcelsBuilder.cc
        src/paths/long/KmerAlign.cc src/paths/ProcessGap.cc src/random/NormalDistribution.cc src/paths/long/EvalByReads.cc src/paths/Uniseq.cc
        src/IndexedAlignmentPlusVector.cc src/graph/FindCells.cc src/paths/long/large/Lines.cc src/kmers/BigKPather.cc src/paths/long/ultra/FounderAlignment.cc
        src/paths/FindErrorsCore.cc src/pairwise_aligners/MaxMutmerFromMer.cc src/paths/long/CleanEfasta.cc src/efasta/AmbiguityScore.cc
        src/paths/UnibaseUtils.cc src/ParseRange.cc src/feudal/CharString.cc src/paths/long/Variants.cc src/paths/long/DigraphFromWords.cc
        src/random/MersenneTwister.cc src/math/Combinatorics.cc src/paths/long/AssessBestAlignCore.cc src/pairwise_aligners/ClusterAligner.cc
        src/efasta/AmbiguityScore.cc src/MemberOf.cc src/pairwise_aligners/SmithWatFree.cc src/Qualvector.cc src/system/HostName.cc src/paths/Sepdev.cc
        src/feudal/QualNibbleVec.cc src/Superb.cc src/VecAlignmentPlus.cc src/paths/long/VariantCallTools.cc src/paths/simulation/VCF.cc
        src/util/TextTable.cc src/paths/long/SupportedHyperBasevector7.cc src/paths/long/large/CN1PeakFinder.cc src/paths/long/VariantFilters.cc
        src/paths/long/HBVFromEdges.cc src/paths/long/VariantReadSupport.cc src/paths/long/CorrectByStack.cc src/paths/long/VariantPostProcess.cc
        src/LinkTimestamp.cc src/LinkTimestamp.cc src/fastg/FastgGraph.cc src/system/UseGDB.cc)

SET(SOURCE_FILES_SIMPLIFY_05 src/modules/simplify_05.cc src/paths/HyperBasevector.cc src/Intvector.cc src/Charvector.cc src/paths/long/ReadPath.cc
        src/VecString.cc src/graph/Digraph.cc src/paths/KmerPath.cc src/Basevector.cc src/feudal/PQVec.cc src/feudal/FieldVec.cc
        src/paths/KmerBaseBroker.cc src/Equiv.cc src/feudal/FeudalFileReader.cc src/feudal/Mempool.cc src/feudal/CharString.cc
        src/system/file/FileWriter.cc src/system/file/FileReader.cc src/feudal/FieldVec.cc src/feudal/Mempool.cc src/feudal/BaseVec.cc
        src/system/Assert.cc src/system/file/File.cc src/dna/Bases.cc src/system/Exit.cc src/feudal/Generic.cc src/feudal/Oob.cc
        src/Equiv.cc src/feudal/FeudalFileWriter.cc src/feudal/BinaryStream.cc src/paths/long/large/Lines.cc src/paths/long/large/Simplify.cc
        src/paths/long/large/GapToyTools3.cc src/system/System.cc src/paths/long/large/GapToyTools5.cc src/paths/long/large/GapToyTools3.cc
        src/system/ThreadsafeIO.cc src/paths/long/ReadStack.cc src/paths/long/FriendAligns.cc src/math/Matrix.cc src/efasta/EfastaTools.cc
        src/system/ErrNo.cc src/feudal/FeudalControlBlock.cc src/paths/KmerPathInterval.cc src/PairsManager.cc src/math/PowerOf2.cc
        src/kmers/KmerRecord.cc src/kmers/ReadPather.cc src/system/WorklistUtils.cc src/paths/long/large/GapToyTools.cc src/paths/long/large/CN1PeakFinder.cc
        src/PairsManager.cc src/kmers/KMerContext.cc src/paths/long/EMEC3.cc src/TokenizeString.cc src/paths/long/large/GapToyTools2.cc
        src/graphics/BasicGraphics.cc src/math/Functions.cc src/system/SysConf.cc src/paths/long/ReadPathTools.cc src/PrintAlignment.cc
        src/math/HoInterval.cc src/paths/long/large/GapToyTools4.cc src/kmers/MakeLookup.cc src/system/Thread.cc src/paths/long/LongReadsToPaths.cc
        src/paths/RemodelGapTools.cc src/ParseSet.cc src/paths/long/LongProtoTools.cc src/paths/long/large/tools/NhoodInfoCore.cc
        src/graph/FindCells.cc src/paths/long/large/GapToyTools6.cc src/paths/long/large/ImprovePath.cc src/random/Bernoulli.cc
        src/paths/long/DiscovarTools.cc src/paths/HyperBasevector.cc src/PackAlign.cc src/paths/HyperEfasta.cc src/paths/long/KmerCount.cc
        src/paths/long/Friends.cc src/kmers/BigKPather.cc src/feudal/BitVec.cc src/paths/long/SupportedHyperBasevector3.cc
        src/paths/long/ultra/ThreadedBlocks.cc src/system/ParsedArgs.cc src/lookup/LookAlign.cc src/system/file/TempFile.cc
        src/reporting/PerfStat.cc src/util/Logger.cc src/random/RNGen.cc src/kmers/LongReadPather.cc src/paths/long/LargeKDispatcher.cc
        src/paths/long/SupportedHyperBasevector6.cc src/paths/long/large/tools/NhoodInfoState.cc src/FastIfstream.cc src/system/ParsedArgs.cc
        src/Fastavector.cc src/math/Permutation.cc src/paths/long/LongHyper.cc src/paths/Unipath.cc src/random/NormalDistribution.cc
        src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/paths/ReadsToPathsCoreX.cc src/paths/long/LoadCorrectCore.cc
        src/paths/long/fosmid/Fosmids.cc src/Vec.cc src/Qualvector.cc src/paths/long/EvalAssembly.cc src/system/UseGDB.cc src/system/ProcBuf.cc
        src/system/TraceVal.cc src/paths/long/SupportedHyperBasevector7.cc src/paths/long/ultra/ConsensusScoreModel.cc src/math/IntDistribution.cc
        src/paths/long/RefTraceTools.cc src/paths/long/ultra/ConsensusScoreModel.cc src/feudal/CharString.cc src/lookup/SAM.cc
        src/paths/long/ultra/MultipleAligner.cc src/fastg/FastgGraph.cc src/paths/long/SupportedHyperBasevector.cc src/pairwise_aligners/SmithWatBandedA.cc
        src/BlockAlign.cc src/random/Shuffle.cc src/paths/long/large/tools/NhoodInfoStuff.cc src/paths/long/KmerAlign.cc src/paths/long/FillPairs.cc
        src/paths/long/Correct1.cc src/lookup/SAM2CRD.cc src/paths/LongReadTools.cc src/feudal/FeudalTools.cc src/paths/MakeAlignsPathsParallelX.cc
        src/pairwise_aligners/SmithWatAffine.cc src/paths/BigMapTools.cc src/paths/long/ultra/FounderAlignment.cc src/paths/FindErrorsCore.cc
        src/paths/long/LongProtoTools.cc src/paths/long/CleanEfasta.cc src/paths/UnipathScaffold.cc src/feudal/QualNibbleVec.cc
        src/paths/long/SupportedHyperBasevector4.cc src/IntPairVec.cc src/paths/long/SupportedHyperBasevector5.cc src/Alignment.cc
        src/kmers/kmer_parcels/KmerParcelsBuilder.cc src/paths/long/ReadOriginTracker.cc src/paths/long/RefTraceControl.cc
        src/paths/long/RefTrace.cc src/paths/ProcessGap.cc src/paths/Uniseq.cc src/paths/Sepdev.cc src/IndexedAlignmentPlusVector.cc
        src/IndexedAlignmentPlusVector.cc src/Superb.cc src/VecAlignmentPlus.cc src/paths/long/ultra/Prefab.cc src/pairwise_aligners/MaxMutmerFromMer.cc
        src/paths/long/OverlapReads.cc src/paths/MergeReadSetsCore.cc src/efasta/AmbiguityScore.cc src/paths/long/HBVFromEdges.cc
        src/paths/UnibaseUtils.cc src/ParseRange.cc src/Fastavector.cc src/paths/long/Variants.cc src/paths/long/PreCorrectOldNew.cc
        src/paths/long/DigraphFromWords.cc src/random/MersenneTwister.cc src/paths/long/EvalByReads.cc src/math/Combinatorics.cc
        src/paths/long/PreCorrectAlt1.cc src/paths/long/VariantCallTools.cc src/paths/simulation/VCF.cc src/util/TextTable.cc
        src/paths/long/VariantFilters.cc src/paths/long/VariantReadSupport.cc src/ScoreAlignment.cc src/pairwise_aligners/SmithWatFree.cc
        src/paths/long/ultra/MakeBlocks.cc src/paths/long/DisplayTools.cc src/paths/long/MakeAlignments.cc src/paths/long/AssessBestAlignCore.cc
        src/paths/long/CorrectByStack.cc src/pairwise_aligners/ClusterAligner.cc src/MemberOf.cc src/paths/long/CorrectPairs1.cc
        src/system/HostName.cc src/paths/long/large/LocalLayout.cc src/paths/long/Correct1Pre.cc src/paths/long/CreateGenome.cc
        src/paths/long/VariantPostProcess.cc src/LinkTimestamp.cc src/paths/long/fosmid/FosmidPool.cc src/FetchReads.cc src/Floatvector.cc)

SET(SOURCE_FILES_SCAFFOLDING_06 src/modules/scaffolding_06.cc src/paths/HyperBasevector.cc src/graph/Digraph.cc src/paths/KmerPath.cc src/Intvector.cc
        src/Charvector.cc src/paths/long/ReadPath.cc src/IntPairVec.cc src/Basevector.cc src/feudal/FieldVec.cc src/math/Matrix.cc
        src/paths/KmerBaseBroker.cc src/Equiv.cc src/paths/KmerPathInterval.cc src/feudal/FeudalFileReader.cc src/feudal/Mempool.cc
        src/feudal/BaseVec.cc src/feudal/CharString.cc src/system/file/FileWriter.cc src/system/file/FileReader.cc src/feudal/PQVec.cc
        src/system/Assert.cc src/system/file/File.cc src/dna/Bases.cc src/system/Exit.cc src/feudal/Generic.cc src/feudal/Oob.cc
        src/system/ThreadsafeIO.cc src/system/ErrNo.cc src/feudal/FeudalControlBlock.cc src/random/RNGen.cc src/feudal/FeudalFileWriter.cc
        src/feudal/BinaryStream.cc src/math/Permutation.cc src/paths/long/large/MakeGaps.cc src/paths/long/large/Samples.cc
        src/FastIfstream.cc src/system/ProcBuf.cc src/TokenizeString.cc src/system/System.cc src/paths/long/large/FinalFiles.cc
        src/paths/long/DisplayTools.cc src/paths/long/large/GapToyTools6.cc src/feudal/BitVec.cc src/system/file/TempFile.cc src/kmers/KmerRecord.cc src/paths/long/large/GapToyTools.cc
        src/paths/long/large/Lines.cc src/efasta/EfastaTools.cc src/feudal/CharString.cc src/paths/long/SupportedHyperBasevector3.cc
        src/math/PowerOf2.cc src/paths/long/SupportedHyperBasevector6.cc src/system/WorklistUtils.cc src/paths/long/KmerCount.cc src/Charvector.cc
        src/paths/long/DiscovarTools.cc src/paths/long/CreateGenome.cc src/paths/long/large/LocalLayout.cc src/ParseSet.cc
        src/paths/long/EvalAssembly.cc src/paths/long/fosmid/Fosmids.cc src/Vec.cc src/paths/long/large/GapToyTools5.cc src/PackAlign.cc
        src/system/TraceVal.cc src/lookup/SAM.cc src/paths/long/SupportedHyperBasevector4.cc src/pairwise_aligners/SmithWatBandedA.cc
        src/paths/long/RefTraceTools.cc src/paths/long/large/CN1PeakFinder.cc src/PairsManager.cc src/paths/long/LongHyper.cc
        src/paths/Unipath.cc src/paths/long/large/GapToyTools4.cc src/FetchReads.cc src/paths/long/EvalByReads.cc src/graph/FindCells.cc
        src/paths/long/LongProtoTools.cc src/pairwise_aligners/ClusterAligner.cc src/paths/HyperEfasta.cc src/paths/long/Friends.cc
        src/VecString.cc src/paths/long/ultra/ThreadedBlocks.cc src/math/IntDistribution.cc src/system/ParsedArgs.cc src/Fastavector.cc
        src/lookup/LookAlign.cc src/Floatvector.cc src/system/UseGDB.cc src/paths/long/SupportedHyperBasevector7.cc src/paths/long/ultra/ConsensusScoreModel.cc
        src/reporting/PerfStat.cc src/paths/long/ultra/MultipleAligner.cc src/util/Logger.cc src/kmers/ReadPather.cc src/fastg/FastgGraph.cc
        src/kmers/BigKPather.cc src/paths/long/SupportedHyperBasevector5.cc src/paths/long/ReadOriginTracker.cc src/paths/long/RefTraceControl.cc
        src/BlockAlign.cc src/kmers/KMerContext.cc src/paths/long/large/GapToyTools3.cc src/paths/LongReadTools.cc src/math/Functions.cc
        src/paths/long/ultra/GetFriendsAndAlignsInitial.cc src/paths/long/ReadStack.cc src/paths/long/FriendAligns.cc src/Alignment.cc
        src/paths/long/LargeKDispatcher.cc src/paths/long/ultra/MakeBlocks.cc src/efasta/AmbiguityScore.cc src/Superb.cc
        src/random/Bernoulli.cc src/graphics/BasicGraphics.cc src/system/HostName.cc src/pairwise_aligners/SmithWatFree.cc src/math/HoInterval.cc
        src/MemberOf.cc src/paths/long/MakeAlignments.cc src/ScoreAlignment.cc src/system/SysConf.cc src/math/Combinatorics.cc
        src/paths/long/LoadCorrectCore.cc src/paths/RemodelGapTools.cc src/paths/long/HBVFromEdges.cc src/paths/long/LongReadsToPaths.cc
        src/kmers/LongReadPather.cc src/IndexedAlignmentPlusVector.cc src/paths/long/EMEC3.cc src/random/Shuffle.cc src/paths/long/KmerAlign.cc
        src/paths/long/FillPairs.cc src/paths/long/RefTrace.cc src/paths/long/Correct1.cc src/lookup/SAM2CRD.cc src/paths/long/SupportedHyperBasevector.cc
        src/random/NormalDistribution.cc src/feudal/FeudalTools.cc src/pairwise_aligners/SmithWatAffine.cc src/paths/long/fosmid/FosmidPool.cc
        src/PrintAlignment.cc src/paths/UnibaseUtils.cc src/paths/ReadsToPathsCoreX.cc src/ParseRange.cc src/LinkTimestamp.cc
        src/paths/long/Correct1Pre.cc src/paths/MergeReadSetsCore.cc src/paths/long/CorrectPairs1.cc src/paths/long/CorrectByStack.cc
        src/paths/long/AssessBestAlignCore.cc src/random/MersenneTwister.cc src/paths/long/PreCorrectAlt1.cc src/paths/long/DigraphFromWords.cc
        src/paths/long/PreCorrectOldNew.cc src/paths/long/Variants.cc src/feudal/CharString.cc src/ParseRange.cc src/paths/long/CleanEfasta.cc
        src/kmers/MakeLookup.cc src/paths/long/ultra/Prefab.cc src/paths/long/LongProtoTools.cc src/paths/long/ReadPathTools.cc
        src/VecAlignmentPlus.cc src/paths/long/OverlapReads.cc src/paths/long/VariantCallTools.cc src/paths/simulation/VCF.cc
        src/util/TextTable.cc src/paths/MakeAlignsPathsParallelX.cc src/paths/long/VariantFilters.cc src/paths/long/ultra/FounderAlignment.cc
        src/paths/FindErrorsCore.cc src/system/Thread.cc src/feudal/QualNibbleVec.cc src/kmers/kmer_parcels/KmerParcelsBuilder.cc
        src/paths/BigMapTools.cc src/pairwise_aligners/MaxMutmerFromMer.cc src/paths/long/VariantReadSupport.cc src/Qualvector.cc
        src/paths/long/VariantPostProcess.cc src/paths/UnipathScaffold.cc src/paths/ProcessGap.cc src/paths/Uniseq.cc src/paths/Sepdev.cc)

## add executables and respective sources
add_executable(01_unipaths ${SOURCE_FILES_UNIPATH_01})
add_executable(02_qgraph ${SOURCE_FILES_QGRAPH_02})
add_executable(03_clean ${SOURCE_FILES_CLEAN_03})
add_executable(04_patching ${SOURCE_FILES_PATCHING_04})
add_executable(05_simplify ${SOURCE_FILES_SIMPLIFY_05})
add_executable(06_scaffolding ${SOURCE_FILES_SCAFFOLDING_06})

## Link libraries
target_link_libraries(01_unipaths libtbbmalloc.so)
target_link_libraries(02_qgraph libtbbmalloc.so)
target_link_libraries(03_clean libtbbmalloc.so)
target_link_libraries(04_patching libtbbmalloc.so)
target_link_libraries(05_simplify libtbbmalloc.so)
target_link_libraries(06_scaffolding libtbbmalloc.so)
