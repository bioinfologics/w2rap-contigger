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

If you want to link a particular malloc library (jemalloc gives the best performance right now) set the `MALLOC_LIBRARY` variable during cmake:

```
git clone https://github.com/gonzalogacc/w2rap-contigger.git
cd w2rap-contigger
cmake -D CMAKE_CXX_COMPILER=g++ -D MALLOC_LIBRARY=<path_to_library.so> .  
make -j 4
```

## Running w2rap-contigger

You need to create a new directory for the intermediate and output files, then run each of the 6 steps sequentially.

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
