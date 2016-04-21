cmake_minimum_required (VERSION 2.8)
project (A-Discovar)

## Compiler flags
SET(CMAKE_CXX_FLAGS "-g -fopenmp -fpermissive -w -Wno-deprecated -std=c++11 -Ofast -lz")
#SET(CMAKE_CXX_LINK_EXECUTABLE "/usr/bin/ld")
add_definitions(${CMAKE_CXX_FLAGS})

## To move the bins to the final location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

## Include directories
include_directories("../")

## Link directories
link_directories("/tgac/software/production/intel_studio/2015.2/x86_64/composerxe/tbb/lib/intel64/gcc4.4")


## Set sources variables
SET(SOURCE_FILES_UNIPATH_01 unipaths_01.cc ../paths/long/large/GapToyTools6.cc ../paths/long/DiscovarTools.cc ../system/System.cc
        ../feudal/CharString.cc ../paths/long/large/Samples.cc ../paths/long/large/ExtractReads.cc ../system/Assert.cc
        ../system/file/FileReader.cc ../system/file/FileWriter.cc ../Charvector.cc ../Basevector.cc ../feudal/PQVec.cc
        ../system/file/File.cc ../system/Exit.cc ../system/SysConf.cc ../FastIfstream.cc ../paths/long/fosmid/Fosmids.cc
        ../ParseSet.cc ../lookup/SAM.cc ../Fastavector.cc ../util/Logger.cc ../system/ErrNo.cc ../system/file/TempFile.cc
        ../feudal/Mempool.cc ../TokenizeString.cc ../feudal/BinaryStream.cc ../feudal/FeudalTools.cc ../bam/ReadBAM.cc
        ../feudal/FieldVec.cc ../paths/long/large/ReadNameLookup.cc ../Qualvector.cc ../feudal/BaseVec.cc ../feudal/BitVec.cc
        ../feudal/QualNibbleVec.cc ../FetchReads.cc ../paths/long/LoadCorrectCore.cc ../dna/Bases.cc ../random/Shuffle.cc
        ../PairsManager.cc ../feudal/FeudalFileReader.cc ../VecString.cc ../feudal/FeudalFileWriter.cc ../system/WorklistUtils.cc
        ../feudal/Generic.cc ../feudal/Oob.cc ../system/Thread.cc ../lookup/SAM2CRD.cc ../system/ProcBuf.cc ../math/HoInterval.cc
        ../random/RNGen.cc ../graph/Digraph.cc ../Equiv.cc ../feudal/FeudalControlBlock.cc ../Floatvector.cc ../paths/MergeReadSetsCore.cc
        ../lookup/LookAlign.cc ../paths/long/LongProtoTools.cc ../paths/long/PreCorrectAlt1.cc ../paths/long/PreCorrectOldNew.cc
        ../paths/long/FillPairs.cc ../paths/KmerPath.cc ../paths/long/Correct1Pre.cc ../paths/FindErrorsCore.cc ../efasta/EfastaTools.cc
        ../paths/long/LongReadsToPaths.cc ../paths/long/CorrectPairs1.cc ../kmers/KmerRecord.cc ../paths/long/Correct1.cc
        ../paths/HyperBasevector.cc ../random/MersenneTwister.cc ../Vec.cc ../system/ThreadsafeIO.cc ../math/Matrix.cc
        ../PackAlign.cc ../ScoreAlignment.cc ../PrintAlignment.cc ../BlockAlign.cc ../system/ParsedArgs.cc ../paths/long/ReadStack.cc
        ../IntPairVec.cc ../Intvector.cc ../paths/OffsetTracker.cc ../kmers/LongReadPather.cc ../reporting/PerfStat.cc
        ../ParseRange.cc ../math/PowerOf2.cc ../system/TraceVal.cc ../pairwise_aligners/ClusterAligner.cc ../paths/KmerPathInterval.cc
        ../paths/long/ultra/ThreadedBlocks.cc ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../paths/long/ultra/ConsensusScoreModel.cc
        ../math/Permutation.cc ../paths/long/ultra/Prefab.cc ../paths/long/ultra/MakeBlocks.cc ../paths/long/Friends.cc ../math/Combinatorics.cc
        ../paths/long/OverlapReads.cc ../random/Bernoulli.cc ../paths/long/ultra/FounderAlignment.cc ../paths/long/KmerAlign.h
        ../paths/KmerBaseBroker.cc ../paths/long/ShortKmerReadPather.cc ../MemberOf.cc ../paths/long/ultra/MultipleAligner.cc
        ../paths/long/FriendAligns.cc ../paths/long/ReadPathTools.cc ../paths/long/ExtendReadPath.cc ../paths/long/CorrectByStack.cc
        ../kmers/BigKPather.cc ../paths/long/HBVFromEdges.cc ../paths/long/ReadPath.cc ../paths/long/KmerAlign.cc
        ../math/Functions.cc ../paths/BigMapTools.cc ../paths/UnipathScaffold.cc ../paths/long/EMEC3.cc ../paths/HyperEfasta.cc
        ../kmers/KMerContext.cc ../paths/long/LargeKDispatcher.cc ../system/HostName.cc ../LinkTimestamp.cc ../system/UseGDB.cc
        ../paths/MuxGraph.cc ../paths/LongReadTools.cc ../pairwise_aligners/SmithWatBandedA.cc ../pairwise_aligners/SmithWatFree.cc
        ../paths/ProcessGap.cc ../paths/Uniseq.cc ../graph/FindCells.cc ../paths/Sepdev.cc ../Superb.cc ../kmers/ReadPather.cc
        ../efasta/AmbiguityScore.cc ../paths/long/MakeAlignments.cc ../AnnotatedContig.cc ../Misc.cc ../system/Crash.cc
        ../Alignment.cc ../VecAlignmentPlus.cc ../IndexedAlignmentPlusVector.cc)

SET(SOURCE_FILES_QGRAPH_02 qgraph_02.cc ../paths/long/large/GapToyTools6.cc ../paths/long/DiscovarTools.cc
        ../system/System.cc ../feudal/CharString.cc ../paths/HyperBasevector.cc ../paths/long/large/GapToyTools.cc
        ../paths/long/large/Repath.cc ../system/Assert.cc ../feudal/BinaryStream.cc ../system/file/FileReader.cc
        ../system/file/FileWriter.cc ../feudal/Mempool.cc ../feudal/FieldVec.cc ../Charvector.cc
        ../paths/long/ReadPath.cc ../Basevector.cc ../feudal/PQVec.cc ../IntPairVec.cc ../BasevectorTools.cc
        ../Intvector.cc ../feudal/FeudalFileReader.cc ../system/file/File.cc ../feudal/FeudalFileWriter.cc
        ../system/SysConf.cc ../FastIfstream.cc ../paths/long/fosmid/Fosmids.cc ../ParseSet.cc ../lookup/SAM.cc
        ../Fastavector.cc ../util/Logger.cc ../system/Exit.cc ../system/ErrNo.cc ../system/file/TempFile.cc
        ../TokenizeString.cc ../feudal/BaseVec.cc ../paths/KmerPath.cc ../Equiv.cc ../graph/Digraph.cc ../dna/Bases.cc
        ../paths/KmerBaseBroker.cc ../VecString.cc ../paths/long/CreateGenome.cc ../paths/long/EvalAssembly.cc
        ../paths/ReadsToPathsCoreX.cc ../paths/Unipath.cc ../efasta/EfastaTools.cc ../paths/long/LongHyper.cc
        ../paths/long/SupportedHyperBasevector6.cc ../PairsManager.cc ../paths/long/LoadCorrectCore.cc
        ../paths/long/SupportedHyperBasevector3.cc ../paths/long/large/LocalLayout.cc ../lookup/LookAlign.cc
        ../paths/long/large/GapToyTools3.cc ../math/PowerOf2.cc ../system/WorklistUtils.cc ../kmers/KmerRecord.cc
        ../paths/KmerPathInterval.cc ../FastaFileset.cc ../kmers/KMerContext.cc ../feudal/Generic.cc
        ../kmers/ReadPather.cc ../paths/long/HBVFromEdges.cc ../system/ProcBuf.cc ../feudal/FieldVec.cc
        ../paths/long/BuildReadQGraph.cc ../feudal/FeudalControlBlock.cc ../feudal/Oob.cc ../system/Thread.cc
        ../paths/long/ShortKmerReadPather.cc ../paths/long/ExtendReadPath.cc ../system/ThreadsafeIO.cc
        ../FastaFilestream.cc ../FastaFilestreamPreview.cc ../FastaVerifier.cc ../FastaConverter.cc
        ../FastaNameParser.cc ../BlockAlign.cc ../PrintAlignment.cc ../paths/long/SupportedHyperBasevector.cc
        ../paths/long/Correct1.cc ../paths/long/FriendAligns.cc ../math/Matrix.cc ../graphics/BasicGraphics.cc
        ../paths/long/SupportedHyperBasevector4.cc ../paths/long/ReadStack.cc ../paths/HyperEfasta.cc
        ../system/ParsedArgs.cc ../feudal/BitVec.cc ../paths/MakeAlignsPathsParallelX.cc ../ParseRange.cc ../MemberOf.cc
        ../system/TraceVal.cc ../FetchReads.cc ../pairwise_aligners/MaxMutmerFromMer.cc
        ../kmers/kmer_parcels/KmerParcelsBuilder.cc ../pairwise_aligners/SmithWatAffine.cc ../Floatvector.cc
        ../paths/long/EvalByReads.cc ../PackAlign.cc ../reporting/PerfStat.cc ../paths/long/LongProtoTools.cc
        ../paths/MergeReadSetsCore.cc ../paths/long/KmerCount.cc ../paths/long/Friends.cc ../random/RNGen.cc
        ../VecString.cc ../paths/long/ultra/ThreadedBlocks.cc ../paths/long/ultra/ConsensusScoreModel.cc
        ../paths/BigMapTools.cc ../paths/long/ultra/MultipleAligner.cc ../paths/LongReadTools.cc
        ../paths/UnipathScaffold.cc ../pairwise_aligners/SmithWatFree.cc ../paths/long/ultra/Prefab.cc
        ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../Vec.cc ../paths/long/KmerAlign.cc ../paths/long/OverlapReads.cc
        ../paths/long/SupportedHyperBasevector5.cc ../math/HoInterval.cc ../Alignment.cc ../random/Bernoulli.cc
        ../math/Combinatorics.cc ../paths/long/ultra/FounderAlignment.cc ../pairwise_aligners/SmithWatBandedA.cc
        ../math/Permutation.cc ../paths/long/RefTraceTools.cc ../IndexedAlignmentPlusVector.cc ../paths/ProcessGap.cc
        ../paths/long/large/Lines.cc ../paths/Uniseq.cc ../kmers/BigKPather.cc ../paths/long/SupportedHyperBasevector7.cc
        ../paths/Sepdev.cc ../math/IntDistribution.cc ../math/Functions.cc ../paths/UnibaseUtils.cc ../CompressedSequence.cc
        ../paths/long/fosmid/FosmidPool.cc ../feudal/FeudalTools.cc ../paths/long/LongReadsToPaths.cc ../random/Shuffle.cc
        ../paths/long/RefTraceControl.cc ../paths/long/PreCorrectAlt1.cc ../lookup/SAM2CRD.cc ../paths/long/FillPairs.cc
        ../kmers/LongReadPather.cc ../random/MersenneTwister.cc ../fastg/FastgGraph.cc ../paths/long/Correct1Pre.cc
        ../paths/long/DigraphFromWords.cc ../paths/long/PreCorrectOldNew.cc ../paths/long/CleanEfasta.cc ../Superb.cc
        ../util/PeakFinder.h ../graph/FindCells.cc ../efasta/AmbiguityScore.cc ../VecAlignmentPlus.cc ../system/UseGDB.cc
        ../paths/long/PreCorrectOldNew.cc ../paths/long/Variants.cc ../paths/long/large/CN1PeakFinder.cc
        ../paths/long/LargeKDispatcher.cc ../LinkTimestamp.cc ../paths/long/EMEC3.cc ../pairwise_aligners/ClusterAligner.cc
        ../paths/long/RefTrace.cc ../ScoreAlignment.cc ../paths/long/ultra/MakeBlocks.cc ../paths/long/VariantPostProcess.cc
        ../system/HostName.cc ../paths/long/CorrectPairs1.cc ../paths/long/CorrectByStack.cc ../paths/long/AssessBestAlignCore.cc
        ../paths/long/MakeAlignments.cc ../paths/long/ReadPath.cc ../paths/long/VariantReadSupport.cc ../paths/RemodelGapTools.cc
        ../paths/long/ReadPathTools.cc ../paths/long/VariantCallTools.cc ../Qualvector.cc ../paths/long/VariantFilters.cc
        ../paths/FindErrorsCore.cc ../paths/RemodelGapTools.cc ../paths/long/ReadOriginTracker.cc ../Qualvector.cc
        ../random/NormalDistribution.cc ../paths/simulation/VCF.cc ../feudal/QualNibbleVec.cc ../util/TextTable.cc
        ../AnnotatedContig.cc ../system/Crash.cc ../Misc.cc)

SET(SOURCE_FILES_CLEAN_03 clean_03.cc ../paths/long/DiscovarTools.cc ../system/System.cc ../paths/long/large/GapToyTools6.cc
        ../ParseSet.cc ../feudal/CharString.cc ../system/SysConf.cc ../paths/long/fosmid/Fosmids.cc ../TokenizeString.cc ../Fastavector.cc
        ../system/Thread.cc ../FastaFilestream.cc ../Charvector.cc ../system/file/FileReader.cc ../system/file/File.cc
        ../system/ErrNo.cc ../graph/Digraph.cc ../Basevector.cc ../FastaFilestream.cc ../paths/HyperBasevector.cc ../paths/KmerPath.cc
        ../math/HoInterval.cc ../paths/KmerPathInterval.cc ../paths/long/large/Clean200.cc ../Intvector.cc ../feudal/FieldVec.cc
        ../Floatvector.cc ../system/file/FileWriter.cc ../math/Matrix.cc ../dna/Bases.cc ../feudal/Mempool.cc ../Equiv.cc ../system/file/TempFile.cc
        ../feudal/BaseVec.cc ../FastaFilestreamPreview.cc ../Vec.cc ../feudal/Generic.cc ../paths/long/large/GapToyTools3.cc
        ../paths/long/large/GapToyTools.cc ../paths/long/CreateGenome.cc ../feudal/FeudalFileWriter.cc ../feudal/FeudalFileReader.cc
        ../FetchReads.cc ../graphics/BasicGraphics.cc ../paths/long/large/LocalLayout.cc ../system/ThreadsafeIO.cc ../paths/KmerBaseBroker.cc
        ../paths/long/ReadStack.cc ../IntPairVec.cc ../PairsManager.cc ../feudal/FeudalControlBlock.cc ../paths/long/SupportedHyperBasevector3.cc
        ../paths/long/ReadPath.cc ../efasta/EfastaTools.cc ../paths/long/FriendAligns.cc ../feudal/BitVec.cc ../random/Bernoulli.cc
        ../paths/long/EvalAssembly.cc ../lookup/LookAlign.cc ../feudal/PQVec.cc ../PackAlign.cc ../math/PowerOf2.cc ../random/RNGen.cc
        ../system/WorklistUtils.cc ../kmers/KMerContext.cc ../feudal/BinaryStream.cc ../paths/long/RefTraceTools.cc ../paths/long/MakeAlignments.cc
        ../ScoreAlignment.cc ../paths/long/LoadCorrectCore.cc ../paths/RemodelGapTools.cc ../VecString.cc ../paths/long/Correct1Pre.cc
        ../paths/MergeReadSetsCore.cc ../paths/long/LongProtoTools.cc ../kmers/KmerRecord.cc ../system/Assert.cc
        ../kmers/ReadPather.cc ../util/Logger.cc ../system/ProcBuf.cc ../reporting/PerfStat.cc ../system/ParsedArgs.cc ../math/IntDistribution.cc
        ../paths/long/ultra/ThreadedBlocks.cc ../Qualvector.cc ../FastaConverter.cc ../feudal/Oob.cc ../pairwise_aligners/SmithWatBanded.cc
        ../paths/long/SupportedHyperBasevector4.cc ../lookup/SAM.cc ../paths/long/EvalByReads.cc ../paths/HyperEfasta.cc
        ../MemberOf.cc ../pairwise_aligners/ClusterAligner.cc ../paths/long/PreCorrectAlt1.cc ../feudal/FeudalTools.cc
        ../paths/LongReadTools.cc ../paths/long/LongReadsToPaths.cc ../pairwise_aligners/SmithWatAffine.cc ../ParseRange.cc
        ../paths/ReadsToPathsCoreX.cc ../PrintAlignment.cc ../paths/long/ReadPathTools.cc
        ../paths/FindErrorsCore.cc ../kmers/BigKPather.cc ../paths/BigMapTools.cc ../paths/long/fosmid/FosmidPool.cc ../paths/long/RefTrace.cc
        ../pairwise_aligners/SmithWatFree.cc ../paths/long/HBVFromEdges.cc ../paths/long/PreCorrectOldNew.cc ../paths/long/Variants.cc
        ../pairwise_aligners/SmithWatBandedA.cc ../paths/simulation/VCF.cc ../util/TextTable.cc ../paths/UnipathScaffold.cc
        ../feudal/QualNibbleVec.cc ../paths/long/ultra/MultipleAligner.cc ../kmers/LongReadPather.cc ../paths/long/LargeKDispatcher.cc
        ../paths/long/SupportedHyperBasevector6.cc ../FastIfstream.cc ../paths/long/ultra/MakeBlocks.cc ../paths/long/ultra/ConsensusScoreModel.cc
        ../paths/long/ReadOriginTracker.cc ../paths/long/ultra/ConsensusScoreModel.cc ../paths/long/VariantCallTools.cc
        ../CompressedSequence.cc ../paths/long/ReadOriginTracker.cc ../paths/long/RefTraceControl.cc ../system/TraceVal.cc
        ../FastaVerifier.cc ../math/Functions.cc ../lookup/SAM2CRD.cc ../paths/long/SupportedHyperBasevector7.cc ../paths/long/Friends.cc
        ../system/HostName.cc ../paths/long/Correct1.cc ../paths/long/AssessBestAlignCore.cc ../paths/long/CleanEfasta.cc
        ../graph/Digraph.cc ../paths/long/ultra/Prefab.cc ../paths/long/LongProtoTools.cc ../paths/long/VariantFilters.cc
        ../paths/MakeAlignsPathsParallelX.cc ../random/NormalDistribution.cc ../paths/long/VariantPostProcess.cc ../paths/long/CorrectPairs1.cc
        ../paths/long/CorrectByStack.cc ../paths/long/KmerAlign.cc ../paths/long/VariantReadSupport.cc ../pairwise_aligners/MaxMutmerFromMer.cc
        ../paths/long/KmerCount.cc ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../paths/long/large/Lines.cc ../paths/long/FillPairs.cc
        ../paths/long/LongHyper.cc ../random/Shuffle.cc ../graph/FindCells.cc ../random/MersenneTwister.cc ../math/Combinatorics.cc
        ../paths/UnibaseUtils.cc ../paths/long/ultra/FounderAlignment.cc ../fastg/FastgGraph.cc ../paths/long/OverlapReads.cc
        ../paths/long/SupportedHyperBasevector5.cc ../kmers/kmer_parcels/KmerParcelsBuilder.cc ../system/Exit.cc ../BlockAlign.cc
        ../math/Permutation.cc ../paths/long/EMEC3.cc ../paths/ProcessGap.cc ../paths/Unipath.cc ../efasta/AmbiguityScore.cc
        ../paths/Uniseq.cc ../efasta/AmbiguityScore.cc ../paths/long/large/CN1PeakFinder.cc ../paths/Sepdev.cc ../Alignment.cc
        ../Superb.cc ../IndexedAlignmentPlusVector.cc ../VecAlignmentPlus.cc ../paths/long/SupportedHyperBasevector.cc ../LinkTimestamp.cc
        ../system/UseGDB.cc ../paths/long/DigraphFromWords.cc)

SET(SOURCE_FILES_PATCHING_04 patching_04.cc ../paths/HyperBasevector.cc ../Intvector.cc ../paths/long/ReadPath.cc ../VecString.cc
        ../paths/KmerPath.cc ../graph/Digraph.cc ../Charvector.cc ../Basevector.cc ../feudal/PQVec.cc ../feudal/FieldVec.cc
        ../paths/KmerBaseBroker.cc ../Equiv.cc ../feudal/FeudalFileReader.cc ../feudal/Mempool.cc ../system/file/FileWriter.cc
        ../system/file/FileReader.cc ../system/Assert.cc ../feudal/BaseVec.cc ../math/Matrix.cc ../system/ErrNo.cc ../paths/KmerPathInterval.cc
        ../feudal/CharString.cc ../feudal/FeudalControlBlock.cc ../system/file/File.cc ../dna/Bases.cc ../feudal/Generic.cc ../system/Exit.cc
        ../feudal/Oob.cc ../feudal/FeudalFileWriter.cc ../feudal/BinaryStream.cc ../math/Permutation.cc ../paths/long/large/GapToyTools3.cc
        ../system/System.cc ../paths/long/ReadStack.cc ../PairsManager.cc ../paths/long/DiscovarTools.cc ../paths/long/large/GapToyTools4.cc
        ../random/Bernoulli.cc ../paths/long/large/GapToyTools.cc ../paths/long/large/GapToyTools6.cc ../Vec.cc ../math/HoInterval.cc
        ../paths/long/large/AssembleGaps.cc ../ParseSet.cc ../paths/long/large/GapToyTools5.cc ../graphics/BasicGraphics.cc
        ../paths/long/FriendAligns.cc ../system/TraceVal.cc ../efasta/EfastaTools.cc ../feudal/BitVec.cc ../paths/long/SupportedHyperBasevector3.cc
        ../lookup/LookAlign.cc ../system/file/TempFile.cc ../math/PowerOf2.cc ../paths/long/CreateGenome.cc ../paths/long/large/GapToyTools2.cc
        ../paths/long/large/LocalLayout.cc ../paths/long/EvalAssembly.cc ../paths/long/fosmid/Fosmids.cc ../paths/long/LongProtoTools.cc
        ../paths/long/MakeAlignments.cc ../ScoreAlignment.cc ../system/SysConf.cc ../paths/long/LoadCorrectCore.cc ../paths/RemodelGapTools.cc
        ../paths/long/fosmid/Fosmids.cc ../paths/long/LongReadsToPaths.cc ../TokenizeString.cc ../PackAlign.cc ../paths/long/KmerCount.cc
        ../paths/HyperEfasta.cc ../paths/long/Friends.cc ../FetchReads.cc ../paths/long/ultra/MakeBlocks.cc ../paths/long/Correct1Pre.cc
        ../paths/MergeReadSetsCore.cc ../paths/long/RefTraceTools.cc ../paths/long/CorrectPairs1.cc ../paths/long/PreCorrectAlt1.cc
        ../paths/long/PreCorrectOldNew.cc ../FastIfstream.cc ../Fastavector.cc ../paths/ReadsToPathsCoreX.cc ../system/Thread.cc
        ../kmers/MakeLookup.cc ../paths/long/ultra/Prefab.cc ../pairwise_aligners/SmithWatAffine.cc ../PrintAlignment.cc
        ../paths/long/ReadPathTools.cc ../paths/long/OverlapReads.cc ../paths/long/ultra/ThreadedBlocks.cc ../math/IntDistribution.cc
        ../system/ParsedArgs.cc ../paths/long/ultra/ConsensusScoreModel.cc ../Floatvector.cc ../IntPairVec.cc ../system/ThreadsafeIO.cc
        ../reporting/PerfStat.cc ../system/ProcBuf.cc ../paths/long/ultra/MultipleAligner.cc ../util/Logger.cc ../kmers/ReadPather.cc
        ../random/RNGen.cc ../kmers/KmerRecord.cc ../kmers/LongReadPather.cc ../paths/long/LargeKDispatcher.cc ../paths/long/SupportedHyperBasevector4.cc
        ../paths/long/SupportedHyperBasevector6.cc ../pairwise_aligners/SmithWatBandedA.cc ../system/WorklistUtils.cc ../paths/long/ReadOriginTracker.cc
        ../BlockAlign.cc ../kmers/KMerContext.cc ../paths/long/EMEC3.cc ../random/Shuffle.cc ../paths/long/LongHyper.cc ../paths/long/FillPairs.cc
        ../paths/long/RefTrace.cc ../paths/long/RefTraceControl.cc ../paths/long/SupportedHyperBasevector5.cc
        ../paths/long/Correct1.cc ../paths/Unipath.cc ../lookup/SAM.cc ../lookup/SAM2CRD.cc ../paths/long/SupportedHyperBasevector.cc
        ../paths/LongReadTools.cc ../paths/long/large/Unsat.cc ../paths/long/large/GapToyTools3.cc ../math/Functions.cc ../paths/RemodelGapTools.cc
        ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../feudal/FeudalTools.cc ../paths/MakeAlignsPathsParallelX.cc
        ../paths/long/fosmid/FosmidPool.cc ../paths/BigMapTools.cc ../paths/UnipathScaffold.cc ../Alignment.cc ../kmers/kmer_parcels/KmerParcelsBuilder.cc
        ../paths/long/KmerAlign.cc ../paths/ProcessGap.cc ../random/NormalDistribution.cc ../paths/long/EvalByReads.cc ../paths/Uniseq.cc
        ../IndexedAlignmentPlusVector.cc ../graph/FindCells.cc ../paths/long/large/Lines.cc ../kmers/BigKPather.cc ../paths/long/ultra/FounderAlignment.cc
        ../paths/FindErrorsCore.cc ../pairwise_aligners/MaxMutmerFromMer.cc ../paths/long/CleanEfasta.cc ../efasta/AmbiguityScore.cc
        ../paths/UnibaseUtils.cc ../ParseRange.cc ../feudal/CharString.cc ../paths/long/Variants.cc ../paths/long/DigraphFromWords.cc
        ../random/MersenneTwister.cc ../math/Combinatorics.cc ../paths/long/AssessBestAlignCore.cc ../pairwise_aligners/ClusterAligner.cc
        ../efasta/AmbiguityScore.cc ../MemberOf.cc ../pairwise_aligners/SmithWatFree.cc ../Qualvector.cc ../system/HostName.cc ../paths/Sepdev.cc
        ../feudal/QualNibbleVec.cc ../Superb.cc ../VecAlignmentPlus.cc ../paths/long/VariantCallTools.cc ../paths/simulation/VCF.cc
        ../util/TextTable.cc ../paths/long/SupportedHyperBasevector7.cc ../paths/long/large/CN1PeakFinder.cc ../paths/long/VariantFilters.cc
        ../paths/long/HBVFromEdges.cc ../paths/long/VariantReadSupport.cc ../paths/long/CorrectByStack.cc ../paths/long/VariantPostProcess.cc
        ../LinkTimestamp.cc ../LinkTimestamp.cc ../fastg/FastgGraph.cc ../system/UseGDB.cc)

SET(SOURCE_FILES_SIMPLIFY_05 simplify_05.cc ../paths/HyperBasevector.cc ../Intvector.cc ../Charvector.cc ../paths/long/ReadPath.cc
        ../VecString.cc ../graph/Digraph.cc ../paths/KmerPath.cc ../Basevector.cc ../feudal/PQVec.cc ../feudal/FieldVec.cc
        ../paths/KmerBaseBroker.cc ../Equiv.cc ../feudal/FeudalFileReader.cc ../feudal/Mempool.cc ../feudal/CharString.cc
        ../system/file/FileWriter.cc ../system/file/FileReader.cc ../feudal/FieldVec.cc ../feudal/Mempool.cc ../feudal/BaseVec.cc
        ../system/Assert.cc ../system/file/File.cc ../dna/Bases.cc ../system/Exit.cc ../feudal/Generic.cc ../feudal/Oob.cc
        ../Equiv.cc ../feudal/FeudalFileWriter.cc ../feudal/BinaryStream.cc ../paths/long/large/Lines.cc ../paths/long/large/Simplify.cc
        ../paths/long/large/GapToyTools3.cc ../system/System.cc ../paths/long/large/GapToyTools5.cc ../paths/long/large/GapToyTools3.cc
        ../system/ThreadsafeIO.cc ../paths/long/ReadStack.cc ../paths/long/FriendAligns.cc ../math/Matrix.cc ../efasta/EfastaTools.cc
        ../system/ErrNo.cc ../feudal/FeudalControlBlock.cc ../paths/KmerPathInterval.cc ../PairsManager.cc ../math/PowerOf2.cc
        ../kmers/KmerRecord.cc ../kmers/ReadPather.cc ../system/WorklistUtils.cc ../paths/long/large/GapToyTools.cc ../paths/long/large/CN1PeakFinder.cc
        ../PairsManager.cc ../kmers/KMerContext.cc ../paths/long/EMEC3.cc ../TokenizeString.cc ../paths/long/large/GapToyTools2.cc
        ../graphics/BasicGraphics.cc ../math/Functions.cc ../system/SysConf.cc ../paths/long/ReadPathTools.cc ../PrintAlignment.cc
        ../math/HoInterval.cc ../paths/long/large/GapToyTools4.cc ../kmers/MakeLookup.cc ../system/Thread.cc ../paths/long/LongReadsToPaths.cc
        ../paths/RemodelGapTools.cc ../ParseSet.cc ../paths/long/LongProtoTools.cc ../paths/long/large/tools/NhoodInfoCore.cc
        ../graph/FindCells.cc ../paths/long/large/GapToyTools6.cc ../paths/long/large/ImprovePath.cc ../random/Bernoulli.cc
        ../paths/long/DiscovarTools.cc ../paths/HyperBasevector.cc ../PackAlign.cc ../paths/HyperEfasta.cc ../paths/long/KmerCount.cc
        ../paths/long/Friends.cc ../kmers/BigKPather.cc ../feudal/BitVec.cc ../paths/long/SupportedHyperBasevector3.cc
        ../paths/long/ultra/ThreadedBlocks.cc ../system/ParsedArgs.cc ../lookup/LookAlign.cc ../system/file/TempFile.cc
        ../reporting/PerfStat.cc ../util/Logger.cc ../random/RNGen.cc ../kmers/LongReadPather.cc ../paths/long/LargeKDispatcher.cc
        ../paths/long/SupportedHyperBasevector6.cc ../paths/long/large/tools/NhoodInfoState.cc ../FastIfstream.cc ../system/ParsedArgs.cc
        ../Fastavector.cc ../math/Permutation.cc ../paths/long/LongHyper.cc ../paths/Unipath.cc ../random/NormalDistribution.cc
        ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../paths/ReadsToPathsCoreX.cc ../paths/long/LoadCorrectCore.cc
        ../paths/long/fosmid/Fosmids.cc ../Vec.cc ../Qualvector.cc ../paths/long/EvalAssembly.cc ../system/UseGDB.cc ../system/ProcBuf.cc
        ../system/TraceVal.cc ../paths/long/SupportedHyperBasevector7.cc ../paths/long/ultra/ConsensusScoreModel.cc ../math/IntDistribution.cc
        ../paths/long/RefTraceTools.cc ../paths/long/ultra/ConsensusScoreModel.cc ../feudal/CharString.cc ../lookup/SAM.cc
        ../paths/long/ultra/MultipleAligner.cc ../fastg/FastgGraph.cc ../paths/long/SupportedHyperBasevector.cc ../pairwise_aligners/SmithWatBandedA.cc
        ../BlockAlign.cc ../random/Shuffle.cc ../paths/long/large/tools/NhoodInfoStuff.cc ../paths/long/KmerAlign.cc ../paths/long/FillPairs.cc
        ../paths/long/Correct1.cc ../lookup/SAM2CRD.cc ../paths/LongReadTools.cc ../feudal/FeudalTools.cc ../paths/MakeAlignsPathsParallelX.cc
        ../pairwise_aligners/SmithWatAffine.cc ../paths/BigMapTools.cc ../paths/long/ultra/FounderAlignment.cc ../paths/FindErrorsCore.cc
        ../paths/long/LongProtoTools.cc ../paths/long/CleanEfasta.cc ../paths/UnipathScaffold.cc ../feudal/QualNibbleVec.cc
        ../paths/long/SupportedHyperBasevector4.cc ../IntPairVec.cc ../paths/long/SupportedHyperBasevector5.cc ../Alignment.cc
        ../kmers/kmer_parcels/KmerParcelsBuilder.cc ../paths/long/ReadOriginTracker.cc ../paths/long/RefTraceControl.cc
        ../paths/long/RefTrace.cc ../paths/ProcessGap.cc ../paths/Uniseq.cc ../paths/Sepdev.cc ../IndexedAlignmentPlusVector.cc
        ../IndexedAlignmentPlusVector.cc ../Superb.cc ../VecAlignmentPlus.cc ../paths/long/ultra/Prefab.cc ../pairwise_aligners/MaxMutmerFromMer.cc
        ../paths/long/OverlapReads.cc ../paths/MergeReadSetsCore.cc ../efasta/AmbiguityScore.cc ../paths/long/HBVFromEdges.cc
        ../paths/UnibaseUtils.cc ../ParseRange.cc ../Fastavector.cc ../paths/long/Variants.cc ../paths/long/PreCorrectOldNew.cc
        ../paths/long/DigraphFromWords.cc ../random/MersenneTwister.cc ../paths/long/EvalByReads.cc ../math/Combinatorics.cc
        ../paths/long/PreCorrectAlt1.cc ../paths/long/VariantCallTools.cc ../paths/simulation/VCF.cc ../util/TextTable.cc
        ../paths/long/VariantFilters.cc ../paths/long/VariantReadSupport.cc ../ScoreAlignment.cc ../pairwise_aligners/SmithWatFree.cc
        ../paths/long/ultra/MakeBlocks.cc ../paths/long/DisplayTools.cc ../paths/long/MakeAlignments.cc ../paths/long/AssessBestAlignCore.cc
        ../paths/long/CorrectByStack.cc ../pairwise_aligners/ClusterAligner.cc ../MemberOf.cc ../paths/long/CorrectPairs1.cc
        ../system/HostName.cc ../paths/long/large/LocalLayout.cc ../paths/long/Correct1Pre.cc ../paths/long/CreateGenome.cc
        ../paths/long/VariantPostProcess.cc ../LinkTimestamp.cc ../paths/long/fosmid/FosmidPool.cc ../FetchReads.cc ../Floatvector.cc)

SET(SOURCE_FILES_SCAFFOLDING_06 scaffolding_06.cc ../paths/HyperBasevector.cc ../graph/Digraph.cc ../paths/KmerPath.cc ../Intvector.cc
        ../Charvector.cc ../paths/long/ReadPath.cc ../IntPairVec.cc ../Basevector.cc ../feudal/FieldVec.cc ../math/Matrix.cc
        ../paths/KmerBaseBroker.cc ../Equiv.cc ../paths/KmerPathInterval.cc ../feudal/FeudalFileReader.cc ../feudal/Mempool.cc
        ../feudal/BaseVec.cc ../feudal/CharString.cc ../system/file/FileWriter.cc ../system/file/FileReader.cc ../feudal/PQVec.cc
        ../system/Assert.cc ../system/file/File.cc ../dna/Bases.cc ../system/Exit.cc ../feudal/Generic.cc ../feudal/Oob.cc
        ../system/ThreadsafeIO.cc ../system/ErrNo.cc ../feudal/FeudalControlBlock.cc ../random/RNGen.cc ../feudal/FeudalFileWriter.cc
        ../feudal/BinaryStream.cc ../math/Permutation.cc ../paths/long/large/MakeGaps.cc ../paths/long/large/Samples.cc
        ../FastIfstream.cc ../system/ProcBuf.cc ../TokenizeString.cc ../system/System.cc ../paths/long/large/FinalFiles.cc
        ../paths/long/DisplayTools.cc ../paths/long/large/GapToyTools6.cc ../feudal/BitVec.cc ../system/file/TempFile.cc ../kmers/KmerRecord.cc ../paths/long/large/GapToyTools.cc
        ../paths/long/large/Lines.cc ../efasta/EfastaTools.cc ../feudal/CharString.cc ../paths/long/SupportedHyperBasevector3.cc
        ../math/PowerOf2.cc ../paths/long/SupportedHyperBasevector6.cc ../system/WorklistUtils.cc ../paths/long/KmerCount.cc ../Charvector.cc
        ../paths/long/DiscovarTools.cc ../paths/long/CreateGenome.cc ../paths/long/large/LocalLayout.cc ../ParseSet.cc
        ../paths/long/EvalAssembly.cc ../paths/long/fosmid/Fosmids.cc ../Vec.cc ../paths/long/large/GapToyTools5.cc ../PackAlign.cc
        ../system/TraceVal.cc ../lookup/SAM.cc ../paths/long/SupportedHyperBasevector4.cc ../pairwise_aligners/SmithWatBandedA.cc
        ../paths/long/RefTraceTools.cc ../paths/long/large/CN1PeakFinder.cc ../PairsManager.cc ../paths/long/LongHyper.cc
        ../paths/Unipath.cc ../paths/long/large/GapToyTools4.cc ../FetchReads.cc ../paths/long/EvalByReads.cc ../graph/FindCells.cc
        ../paths/long/LongProtoTools.cc ../pairwise_aligners/ClusterAligner.cc ../paths/HyperEfasta.cc ../paths/long/Friends.cc
        ../VecString.cc ../paths/long/ultra/ThreadedBlocks.cc ../math/IntDistribution.cc ../system/ParsedArgs.cc ../Fastavector.cc
        ../lookup/LookAlign.cc ../Floatvector.cc ../system/UseGDB.cc ../paths/long/SupportedHyperBasevector7.cc ../paths/long/ultra/ConsensusScoreModel.cc
        ../reporting/PerfStat.cc ../paths/long/ultra/MultipleAligner.cc ../util/Logger.cc ../kmers/ReadPather.cc ../fastg/FastgGraph.cc
        ../kmers/BigKPather.cc ../paths/long/SupportedHyperBasevector5.cc ../paths/long/ReadOriginTracker.cc ../paths/long/RefTraceControl.cc
        ../BlockAlign.cc ../kmers/KMerContext.cc ../paths/long/large/GapToyTools3.cc ../paths/LongReadTools.cc ../math/Functions.cc
        ../paths/long/ultra/GetFriendsAndAlignsInitial.cc ../paths/long/ReadStack.cc ../paths/long/FriendAligns.cc ../Alignment.cc
        ../paths/long/LargeKDispatcher.cc ../paths/long/ultra/MakeBlocks.cc ../efasta/AmbiguityScore.cc ../Superb.cc
        ../random/Bernoulli.cc ../graphics/BasicGraphics.cc ../system/HostName.cc ../pairwise_aligners/SmithWatFree.cc ../math/HoInterval.cc
        ../MemberOf.cc ../paths/long/MakeAlignments.cc ../ScoreAlignment.cc ../system/SysConf.cc ../math/Combinatorics.cc
        ../paths/long/LoadCorrectCore.cc ../paths/RemodelGapTools.cc ../paths/long/HBVFromEdges.cc ../paths/long/LongReadsToPaths.cc
        ../kmers/LongReadPather.cc ../IndexedAlignmentPlusVector.cc ../paths/long/EMEC3.cc ../random/Shuffle.cc ../paths/long/KmerAlign.cc
        ../paths/long/FillPairs.cc ../paths/long/RefTrace.cc ../paths/long/Correct1.cc ../lookup/SAM2CRD.cc ../paths/long/SupportedHyperBasevector.cc
        ../random/NormalDistribution.cc ../feudal/FeudalTools.cc ../pairwise_aligners/SmithWatAffine.cc ../paths/long/fosmid/FosmidPool.cc
        ../PrintAlignment.cc ../paths/UnibaseUtils.cc ../paths/ReadsToPathsCoreX.cc ../ParseRange.cc ../LinkTimestamp.cc
        ../paths/long/Correct1Pre.cc ../paths/MergeReadSetsCore.cc ../paths/long/CorrectPairs1.cc ../paths/long/CorrectByStack.cc
        ../paths/long/AssessBestAlignCore.cc ../random/MersenneTwister.cc ../paths/long/PreCorrectAlt1.cc ../paths/long/DigraphFromWords.cc
        ../paths/long/PreCorrectOldNew.cc ../paths/long/Variants.cc ../feudal/CharString.cc ../ParseRange.cc ../paths/long/CleanEfasta.cc
        ../kmers/MakeLookup.cc ../paths/long/ultra/Prefab.cc ../paths/long/LongProtoTools.cc ../paths/long/ReadPathTools.cc
        ../VecAlignmentPlus.cc ../paths/long/OverlapReads.cc ../paths/long/VariantCallTools.cc ../paths/simulation/VCF.cc
        ../util/TextTable.cc ../paths/MakeAlignsPathsParallelX.cc ../paths/long/VariantFilters.cc ../paths/long/ultra/FounderAlignment.cc
        ../paths/FindErrorsCore.cc ../system/Thread.cc ../feudal/QualNibbleVec.cc ../kmers/kmer_parcels/KmerParcelsBuilder.cc
        ../paths/BigMapTools.cc ../pairwise_aligners/MaxMutmerFromMer.cc ../paths/long/VariantReadSupport.cc ../Qualvector.cc
        ../paths/long/VariantPostProcess.cc ../paths/UnipathScaffold.cc ../paths/ProcessGap.cc ../paths/Uniseq.cc ../paths/Sepdev.cc)

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
