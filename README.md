# w2rap-contigger

An Illumina PE genome contig assembler, can handle large (17Gbp) complex (hexaploid) genomes.

http://bioinfologics.github.io/the-w2rap-contigger/ 

## Installation
### Pre-requisites

* Cmake 2.8.0+  
* GCC 5.2.0+ (you can also use ICC)
* (Optional) jemalloc or another malloc library (intel's work too).

### Compilation instructions
This is a big-ish codebase to compile, so we recomment using the `-j` flag on make to use multiplr processors (the examples use 4, but more is better).

```
git clone https://github.com/gonzalogacc/w2rap-contigger.git
cd w2rap-contigger
cmake -D CMAKE_CXX_COMPILER=g++ .  
make -j 4
```

If you want to link a particular malloc library set the `MALLOC_LIBRARY` variable during cmake:

```
git clone https://github.com/gonzalogacc/w2rap-contigger.git
cd w2rap-contigger
cmake -D CMAKE_CXX_COMPILER=g++ -D MALLOC_LIBRARY=<path_to_library.so> .  
make -j 4
```

Right now Intel's tbbmalloc have the best performance in our systems. Jemalloc ha also improved performance in the past. In both cases it is a small gain and varies from system to system, it can even turn into a loss, so beware.

*Note:* for older versions, due to inneficient allocation patterns, jemalloc used to have quite a positive impact in some scenarios, so if you're using any of the "6 binaries" versions, do link with jemalloc.

## Running w2rap-contigger

You need to create a new directory for the intermediate and output files.

The current version of the w2rap contigger runs in 7 steps. By default the contigger will run each of these steps in order, not dumping unnecessary intermediate files. You can use the `--from_step` and `--to_step` options to start from and finish after particular steps. When run this way, the contigger will read the output files from the previous step and dump the necessary files for the next step to run. If you want to dump the output of every step you can use the `--dump_all 1` option.


Step # | Description | Outputs
:---|---|---
1 | Read loading | binary-formatted reads
2 | Build small k (k=60) graph from reads | small k graph, read paths
3 | Build large K graph from small k graph and reads | large K graph, read paths
4 | Clean large K graph | large K cleaned graph, read paths
5 | Local assemblies on the large K graph "gaps" | large K completed graph, read paths
6 | Graph simplification and PathFinder | large K simplified graph, read paths, raw/contig-lines GFA and fasta
7 | PE-scale scaffolding across gaps in the large K graph | large K simplified graph with jumps, read paths, raw/lines GFA and fasta


###Parallel performance considerations

The code has been optimised with local process binning. You should make sure your system will pin threads (as per openmp definition of) to make the best use of memory locality. This can be achieved setting the `OMP_PROC_BIND` or `GOMP_CPU_AFFINITY`/ `KMP_AFFINITY` variables. Particular optimal settings will depend on your system. If you run your software through a scheduler such as SLURM or PBS, the scheduler should set all variables if correctly configured.

In most systems (specially most NUMA systems), using thread-local allocation should have a positive impact on performance. Whilst many systems will use thread-local by default, or have some smart policy, you should consider setting the `MALLOC_PER_THREAD=1` variable if that improves performance on your system (i.e. linux's default malloc can have a good gain from this).


###Examples
Example run with input bam file, K=260:

```
mkdir test_k260
#optionally, remember to set MALLOC_PER_THREAD and/or OMP_PROC_BIND
export OMP_PROC_BIND=spread
export MALLOC_PER_THREAD=1
./w2rap_contigger -o test_k260 -p example -r example.bam -K 260
```


Example run with input fastq files, K=260 dumping all intermediate files:

```
mkdir test_k260
#optionally, remember to set MALLOC_PER_THREAD and/or OMP_PROC_BIND
export OMP_PROC_BIND=spread
export MALLOC_PER_THREAD=1
./w2rap_contigger -o test_k260 -p example -r example_r1.fastq,example_r2.fastq -K 260
```

Example run with input fastq files, K=260, runing step2 independently (this steps usually has usually the highest memory peak):

```
mkdir test_k260
#optionally, remember to set MALLOC_PER_THREAD and/or OMP_PROC_BIND
export OMP_PROC_BIND=spread
export MALLOC_PER_THREAD=1
./w2rap_contigger -o test_k260 -p example -r example_r1.fastq,example_r2.fastq -K 260 --to_step 1
./w2rap_contigger -o test_k260 -p example -r example_r1.fastq,example_r2.fastq -K 260 --from_step 2 --to_step 2
./w2rap_contigger -o test_k260 -p example -r example_r1.fastq,example_r2.fastq -K 260 --from_step 3
```

Your assembly will end up in the `a.lines.fasta` file.

## Running **old versions** of w2rap-contigger

Previous versions of the contigger have 6 binaries for the different steps of the assembly. You need to create a new directory for the intermediate and output files, then run each of the 6 steps sequentially.

Each step requires the path to the output directory and a prefix. Step 01 also requires a path to the input reads (either a bam file or 2 fastq files separated by a ','). Step 02 accepts an optional -K parameter to change the K of the last stage DBG. All steps accept a -t flag to set the thread count for parallel sections, although not all sections will use all processors.

In most systems (specially most NUMA systems), using thread-local allocation will have a positive impact on performance. We are not aware of cases where it produced a significant negative impact, so we recommend setting the `MALLOC_PER_THREAD=1` variable.

Example run with input bam file:

```
mkdir test_k260
export MALLOC_PER_THREAD=1
./01_unipaths -o test_k260 -p example -r example.bam
./02_qgraph -o test_k260 -p example -K 260
./03_clean -o test_k260 -p example
./04_patching -o test_k260 -p example
./05_simplify -o test_k260 -p example
./06_scaffolding -o test_k260 -p example
```


Example run with input fastq files:

```
mkdir test_k260
export MALLOC_PER_THREAD=1
./01_unipaths -o test_k260 -p example -r example_r1.fastq,example_r2.fastq
./02_qgraph -o test_k260 -p example -K 260
./03_clean -o test_k260 -p example
./04_patching -o test_k260 -p example
./05_simplify -o test_k260 -p example
./06_scaffolding -o test_k260 -p example
```

Your assembly will end up in the `a.lines.fasta` file.
