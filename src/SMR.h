//
// Created by Luis Yanes (EI) on 24/07/2017.
//

#ifndef W2RAP_CONTIGGER_SMR_H
#define W2RAP_CONTIGGER_SMR_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>

#include <fcntl.h>
#include <cstdlib>

#ifndef O_DIRECT
#define O_DIRECT O_RDONLY
#endif


template<typename FileRecord>
class FastBReader {
public:
    FastBReader(const vecbvec &_reads, const std::vector<uint16_t> &_len, const size_t _readsBegin,
                const size_t _readsEnd) :
            reads(_reads),
            len(_len),
            readsBegin(_readsBegin),
            readsEnd(_readsEnd),
            pos(_readsBegin),
            total(_readsEnd - _readsBegin),
            progress( (_readsEnd-_readsBegin) / 100 )
    { }

    bool done() {return pos >= readsEnd;}
    bool next_record(FileRecord &rec) {
        if ( (pos-readsBegin) % progress == 0) {
//#pragma omp critical(reader_progress)
//            std::cout << "Processed " << 100 * (pos-readsBegin) / total << "% reads on thread " << omp_get_thread_num() << std::endl;
        }

        if (pos < readsEnd) {
            rec.first = reads[pos];
            rec.second = len[pos];
            pos++;
            return true;
        }
        return false;
    }
private:
    const vecbvec &reads;
    const std::vector<uint16_t> &len;
    const size_t readsBegin, readsEnd, total;
    const uint64_t progress;
    size_t pos;
};

template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ParamStruct >
class SMR {
public:
    SMR(ParamStruct p, uint64_t maxMem, uint min = 0, const std::string &Otmp = "") :
            parameters(p),
            sizeElements(10000),
            myBatches(0),
            nRecs(0),
            nKmers(0),
            tKmers(0),
            merged_kmers(0),
            tmp2(Otmp),
            minCount(min),
            maxThreads((unsigned int) omp_get_max_threads()) {
        numElementsPerBatch = maxMem/sizeof(RecordType) / maxThreads;
        mergeCount=4;
    }

    void read_from_file(const vecbvec &reads, const std::vector<uint16_t> &len, int numThreads) {
        uint64_t numReadsReduction(0);
        OutputLog(3) << "Begin parallel reduction of " << reads.size() << " reads to batches\n";
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

#pragma omp parallel reduction(+: numReadsReduction)
            {
                size_t from, to;
                from = (reads.size() / omp_get_max_threads()) * (omp_get_thread_num());
                to = (reads.size() / omp_get_max_threads()) * (omp_get_thread_num() + 1);

                if (omp_get_thread_num() + 1 == omp_get_num_threads()) to = reads.size();

                FileReader myFileReader(reads, len, from, to);
                mapElementsToBatches(myFileReader, numReadsReduction);
            }
            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            OutputLog(3) << "Done parallel reduction to batches in " << elapsed_seconds.count() << "s" << std::endl;
            this->nRecs = numReadsReduction;
    };

    void getRecords(RecordType* records) {
        OutputLog(3) << "Merging final files" << std::endl;
        std::vector<std::string> threadFiles(maxThreads);

        for (auto threadID = 0; threadID<maxThreads; threadID++) {
            std::string name(tmp + "thread_" + std::to_string(threadID) + ".tmc");
            //std::cout << "Opening file " << name << std::endl;
            threadFiles[threadID] = name;
        }


        nKmers = merge("final.kc", threadFiles);
        OutputLog(3) << "Done merging final files" << std::endl;
        std::ifstream outf("final.kc", std::ios_base::binary);
        elements = (RecordType*) malloc(nKmers*sizeof(RecordType));
        outf.read((char*) elements, sizeof(RecordType)*nKmers);

        std::swap(records,elements);

        std::cout << "Total Kmers " << tKmers << "\n";
        std::cout << "Final nKmers " << nKmers << "\n";
        std::cout << "Final reads " << nRecs << "\n";
    }

private:

    uint64_t merge(const std::string &tmpName, const std::vector<std::string> &files, std::vector<KMerNodeFreq_s> &elements) {
        auto files_size(files.size());
        auto numMemoryElements(elements.size());
        uint64_t memoryElement(0);
        std::vector<int> in_fds( files_size );
        int out = ::open(tmpName.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0777);

        std::vector<bool> active_file(files_size,true);

        std::vector<KMerNodeFreq_s*> next_element_from_file(files.size());
        std::vector<uint> count_element_from_file(files_size,0);
        std::vector<uint> size_element_from_file(files_size,0);
        std::vector<uint> seen_element_from_file(files_size,0);
        uint64_t outCount(0), count0(0), bufferSize(10000000);
        uint64_t outBufferSize(50000000);

        ::write(out, &outCount, sizeof(outCount));

        KMerNodeFreq_s *outfreqs=(KMerNodeFreq_s*) malloc(outBufferSize* sizeof(KMerNodeFreq_s));

        for (auto i = 0; i < files_size; i++) {
            in_fds[i] = ::open(files[i].c_str(), O_RDONLY);
            uint64_t size;
            ::read(in_fds[i], &size, sizeof(size));
            //std::cout << files[i] << " size " << size << " in disk " << sizeof(size) + size* sizeof(KMerNodeFreq_s) << std::endl;
            next_element_from_file[i] = (KMerNodeFreq_s*) malloc(bufferSize* sizeof(KMerNodeFreq_s));
            count_element_from_file[i] = 0;
            size_element_from_file[i]=
                    ::read(in_fds[i], next_element_from_file[i], bufferSize* sizeof(KMerNodeFreq_s)) /
                    sizeof(KMerNodeFreq_s);
            active_file[i] = true;
        }

        uint min = 0;
        for (uint i = 1; i < files_size; ++i) {
            if (next_element_from_file[i][0] < next_element_from_file[i-1][0]) min = i;
        }

        KMerNodeFreq_s current_element;

        current_element = std::min(elements[0],next_element_from_file[min][0]);

        current_element.count=0;


        KMerNodeFreq_s min_element;
        bool active;
        do {
            active=false;
            min_element.kdata[0] = std::numeric_limits<uint64_t>::max();
            min_element.kdata[1] = std::numeric_limits<uint64_t>::max();

            for (int i=0; i<files_size; ++i) {
                if (size_element_from_file[i] <= count_element_from_file[i]) continue;

                if (current_element == next_element_from_file[i][count_element_from_file[i]]) {
                    current_element.merge(next_element_from_file[i][count_element_from_file[i]]);
                    count_element_from_file[i]++;
                    if (count_element_from_file[i] == size_element_from_file[i]) {
                        size_element_from_file[i] =
                                ::read(in_fds[i], (char *) next_element_from_file[i], bufferSize*sizeof(KMerNodeFreq_s))
                                / sizeof(KMerNodeFreq_s);
                        seen_element_from_file[i]+=count_element_from_file[i];
                        count_element_from_file[i] = 0;
                    }
                }
                if (count_element_from_file[i] < size_element_from_file[i]) {
                    if (!(next_element_from_file[i][count_element_from_file[i]] > min_element)){
                        active=true;
                        min_element = next_element_from_file[i][count_element_from_file[i]];
                    }
                }
            }
            if (memoryElement < numMemoryElements) {
                if (current_element == elements[memoryElement]) {
                    current_element.merge(elements[memoryElement]);
                    memoryElement++;
                }
            }
            if (memoryElement<numMemoryElements and !(elements[memoryElement] > min_element)) {
                active=true;
                min_element = elements[memoryElement];
            }

            outfreqs[count0] = current_element;
            ++count0;
            if (count0 == outBufferSize) {
                ::write(out, (char *) outfreqs, outBufferSize * sizeof(KMerNodeFreq_s));
                outCount += count0;
                count0 = 0;
            }
            current_element = min_element;
            current_element.count=0;
        } while ( active );

        if (count0 > 0) {
            ::write(out, (char *) outfreqs, count0 * sizeof(KMerNodeFreq_s));
            outCount += count0;
        }
        ::lseek(out, 0, SEEK_SET);
        ::write(out, (char *) &outCount, sizeof(outCount));
        free(outfreqs);
        ::close(out);

        for (auto i = 0; i < files_size; i++) {
            free(next_element_from_file[i]);
            std::remove(files[i].c_str());
        }
        return outCount;
    }

    uint64_t merge(const std::string &tmpName, const std::vector<std::string> &files) {
        auto files_size(files.size());
        std::vector<int> in_fds( files_size );
        int out = ::open(tmpName.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0777);

        std::vector<bool> active_file(files_size,true);

        std::vector<KMerNodeFreq_s*> next_element_from_file(files.size());
        std::vector<uint> count_element_from_file(files_size,0);
        std::vector<uint> size_element_from_file(files_size,0);
        std::vector<uint> seen_element_from_file(files_size,0);
        uint64_t outCount(0), count0(0), bufferSize(10000000);
        uint64_t outBufferSize(50000000);

        ::write(out, &outCount, sizeof(outCount));

        KMerNodeFreq_s *outfreqs=(KMerNodeFreq_s*) malloc(outBufferSize* sizeof(KMerNodeFreq_s));

        // Open files and initialise batch arrays
        for (auto i = 0; i < files_size; i++) {
            in_fds[i] = ::open(files[i].c_str(), O_RDONLY);
            uint64_t size;
            ::read(in_fds[i], &size, sizeof(size));
            //std::cout << files[i] << " size " << size << " in disk " << sizeof(size) + size* sizeof(KMerNodeFreq_s) << std::endl;
            next_element_from_file[i] = (KMerNodeFreq_s*) malloc(bufferSize* sizeof(KMerNodeFreq_s));
            count_element_from_file[i] = 0;
            size_element_from_file[i]=
                    ::read(in_fds[i], next_element_from_file[i], bufferSize* sizeof(KMerNodeFreq_s)) /
                    sizeof(KMerNodeFreq_s);
            active_file[i] = true;
        }

        // Find the file holding the minimum element
        uint min = 0;
        for (uint i = 1; i < files_size; ++i) {
            if (next_element_from_file[i][0] < next_element_from_file[i-1][0]) min = i;
        }

        // Set the current element to the minimum from the files
        KMerNodeFreq_s current_element(next_element_from_file[min][0]);
        current_element.count=0;
        current_element.kc=0;

        KMerNodeFreq_s min_element;
        bool active;
        do {
            active=false;

            // Set the minimum to a very large value (ensure the min is found)
            min_element.kdata[0] = std::numeric_limits<uint64_t>::max();
            min_element.kdata[1] = std::numeric_limits<uint64_t>::max();

            for (int i=0; i<files_size; ++i) {
                // If this file is depleted continue
                if (size_element_from_file[i] <= count_element_from_file[i]) continue;

                // If this file's element == current element, merge onto the current element
                // and advance the batch-file record
                if (current_element == next_element_from_file[i][count_element_from_file[i]]) {
                    current_element.merge(next_element_from_file[i][count_element_from_file[i]]);
                    count_element_from_file[i]++;

                    // If this batch ran out of elements, get a new batch
                    if (count_element_from_file[i] == size_element_from_file[i]) {
                        // If there are no more elements, size_element_from_file == 0 so it is deactivated
                        size_element_from_file[i] =
                                ::read(in_fds[i], (char *) next_element_from_file[i], bufferSize*sizeof(KMerNodeFreq_s))
                                / sizeof(KMerNodeFreq_s);
                        seen_element_from_file[i]+=count_element_from_file[i];
                        count_element_from_file[i] = 0;
                    }
                }
                // If there are any elements left on this file, check if it is the minimum
                if (count_element_from_file[i] < size_element_from_file[i]) {
                    // If this is the minimum, we need to merge it with the rest
                    if (!(next_element_from_file[i][count_element_from_file[i]] > min_element)){
                        active=true;
                        min_element = next_element_from_file[i][count_element_from_file[i]];
                    }
                }
            }

            // Check if we need to output current element
            if (current_element.count >= minCount) {
                outfreqs[count0] = current_element;
                ++count0;
            }
            merged_kmers++;
            // Check if there's space left on the 'output' batch
            if (count0 == outBufferSize) {
                ::write(out, (char *) outfreqs, outBufferSize * sizeof(KMerNodeFreq_s));
                outCount += count0;
                count0 = 0;
            }
            // Set the new current_element to the min_element and it's count to 0 for the next iteration
            current_element = min_element;
            current_element.count=0;
            current_element.kc=0;
        } while ( active );

        // If there were any elements left on the 'output' batch dump them
        if (count0 > 0) {
            ::write(out, (char *) outfreqs, count0 * sizeof(KMerNodeFreq_s));
            outCount += count0;
        }

        // Write the number of elements on the file at the header
        ::lseek(out, 0, SEEK_SET);
        ::write(out, (char *) &outCount, sizeof(outCount));
        free(outfreqs);
        ::close(out);

        // Cleanup
        for (auto i = 0; i < files_size; i++) {
            free(next_element_from_file[i]);
            //std::remove(files[i].c_str());
        }
        OutputLog(3) << "Total kmers " << outCount << "/" << merged_kmers << " with freq >= " << minCount << std::endl;
        return outCount;
    }

    void mapElementsToBatches(FileReader &myFileReader, uint64_t &numReadsReduction) {
        std::vector<RecordType> _elements;
        uint totalBatches(0);
        uint currentBatch(0);
        _elements.reserve(numElementsPerBatch);
        RecordFactory myRecordFactory(parameters);
        FileRecord frecord;

        while (not myFileReader.done()) {
            RecordType record;
            if (myFileReader.next_record(frecord)) numReadsReduction++;
            myRecordFactory.setFileRecord(frecord);
            while (myRecordFactory.next_element(record)) {
                _elements.push_back(record);
                if (_elements.size() >= numElementsPerBatch) {
                    //std::cout << "Dumped " << totalBatches << " on thread "<< omp_get_thread_num() << std::endl;
                    totalBatches++;
                    tKmers+=_elements.size();
                    if (currentBatch==mergeCount) {
                        mergeBatches(currentBatch, _elements);
                        _elements.clear();
                        currentBatch=1;
                    } else {
                        dumpBatch(currentBatch, _elements);
                        _elements.clear();
                        currentBatch++;
                    }
                }
            }
        }
        tKmers+=_elements.size();

        if (totalBatches) mergeBatches(currentBatch, _elements);
        else {
            dumpBatch(currentBatch, _elements);
        }
        _elements.clear();
        rename(std::string(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_batch_0.tmc").c_str(), std::string(tmp+"thread_"+ std::to_string(omp_get_thread_num()) +".tmc").c_str());
    }

    void dumpBatch(const uint currentBatch, std::vector<RecordType> &_elements) {
        collapse(_elements);
        // Simply dump the file
        std::string outBatchName(std::string(tmp
                                             + "thread_"
                                             + std::to_string(omp_get_thread_num())
                                             + "_batch_" + std::to_string(currentBatch)
                                             + ".tmc"));
        // Simply dump the file
        std::ofstream outBatch(outBatchName.data()
                , std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        auto size = _elements.size();
        outBatch.write((char *)&size, sizeof(size));
        outBatch.write((char *)_elements.data(), size*sizeof(RecordType));
        outBatch.close();
        _elements.clear();
    }

    void mergeBatches(uint currentCount, std::vector<RecordType> &_elements) {
        collapse(_elements);
        uint numFilesToMerge = currentCount;

        std::vector<std::string> inBatches(numFilesToMerge);
        for (auto batch=0; batch < numFilesToMerge; batch++) {
            std::string filename(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_batch_" + std::to_string(batch)+ ".tmc");
            inBatches[batch] = filename;
        }
        std::string threadFile(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_tmpBatch1.tmc");
        merge(threadFile, inBatches, _elements);
        rename(std::string(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_tmpBatch1.tmc").c_str(), std::string(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_batch_0.tmc").c_str());
    }

    void collapse(std::vector<RecordType> &_elements){
        if (_elements.size()==0) return;
        std::sort(_elements.begin(), _elements.end());
        typename std::vector<RecordType>::iterator writePtr;
        writePtr = _elements.begin();
        typename std::vector<RecordType>::const_iterator endPtr;
        endPtr = _elements.end();

        // Keep the total in the first equal KMer
        typename std::vector<RecordType>::const_iterator readPtr;
        for (readPtr = _elements.begin(); readPtr != endPtr; ++writePtr) {
            *writePtr = *readPtr;
            ++readPtr;
            while (readPtr != endPtr and *writePtr == *readPtr) {
                writePtr->merge(*readPtr);
                ++readPtr;
            }
        }
        // After accumulating the values on the first KMer, remove all duplicates using unique
        // this operation is safe because the first copy is kept for all unique values
        // resize the object by erasing the duplicates
        _elements.resize(std::distance(_elements.begin(),writePtr));
    }


    uint minCount;
    ParamStruct parameters;
    RecordType* elements;
    uint64_t numElementsPerBatch;
    std::atomic<uint64_t> myBatches;
    uint64_t nRecs;
    uint64_t nKmers;
    std::atomic<uint64_t> tKmers;
    uint64_t merged_kmers;
    const std::string tmp2;
    const int unsigned maxThreads;
    int mergeCount; // How many batches to keep rolling before merging
    std::string tmp;

    uint64_t sizeElements;
};

#endif //W2RAP_CONTIGGER_SMR_H
