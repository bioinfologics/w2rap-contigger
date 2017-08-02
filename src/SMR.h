//
// Created by Luis Yanes (EI) on 24/07/2017.
//

#ifndef W2RAP_CONTIGGER_SMR_H
#define W2RAP_CONTIGGER_SMR_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <parallel/algorithm>

struct FastqRecord{
    std::string name,seq,qual;
};

template<typename FileRecord>
class FastQReader {
public:
    // Single-end reads constructor
    explicit FastQReader(std::string filepath){
        std::cout << "Opening: " << filepath << "\n";
        buffer = (char *) malloc(bufsize*sizeof(char));
        read.rdbuf()->pubsetbuf(buffer,bufsize);
        read.open(filepath.data());
    }
    bool next_record(FileRecord& rec){
        std::string dummy;
        std::getline(read, record.name);
        std::getline(read, record.seq);
        std::getline(read, dummy);
        std::getline(read, record.qual);
        std::swap(record,rec);
        return !rec.name.empty();
    }

    bool done() {
        return read.eof();
    }

    ~FastQReader(){
        free(buffer);
    }
private:
    std::ifstream read;
    FileRecord record;
    static const size_t bufsize=4*1024*1024;
    char *buffer;
};

template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ParamStruct >
class SMR2 {
public:
    SMR2(ParamStruct p, uint64_t maxMem, const std::string &Otmp = "") :
            parameters(p),
            numMersPerThread(maxMem/sizeof(RecordType)),
            myBatches(0),
            nRecs(0),
            nKmers(0),
            tKmers(0),
            tmp2(Otmp),
            maxThreads(omp_get_max_threads())

    {};

    RecordType *getRecords(FileReader &myFileReader){

        std::cout<<"MAP Threads = "<<omp_get_max_threads()<<std::endl;
        mapElementsToBatches(myFileReader);

        omp_set_num_threads(1);
        reduceElementsFromBatches();
        omp_set_num_threads(maxThreads);

        std::cout << "Total Kmers " << tKmers << "\n";
        std::cout << "Final nKmers " << nKmers << "\n";
        std::cout << "Final reads " << nRecs << "\n";
        return elements;
    }

    void read_from_file(RecordType *records, std::string file1, std::string file2) {
        FileReader myFileReader(file1, file2);
        std::replace(file1.begin(), file1.end(), '/', '.');
        std::replace(file2.begin(), file2.end(), '/', '.');
        tmp = tmp2+file1+"_"+file2;
        std::swap(elements,records);
        getRecords(myFileReader);
        std::swap(records,elements);
    };

    void read_from_file(RecordType *records, std::string file) {
        FileReader myFileReader(file);
        std::replace(file.begin(), file.end(), '/', '.');
        tmp = tmp2+file;
        std::swap(elements,records);
        getRecords(myFileReader);
        std::swap(records,elements);
    };


private:

    void reduceElementsFromBatches() {
        std::vector<std::ifstream> finalMerge(myBatches);
        for (auto batchID = 0; batchID < myBatches; ++batchID) {
            finalMerge[batchID].open(tmp + "_batch_" + std::to_string(batchID + 1)
                                     + ".tmc", std::ios_base::in | std::ios_base::binary);
        }
        std::vector<uint64_t > fileSizes(myBatches);
        for (auto i=0; i < myBatches; ++i) {
            finalMerge[i].read((char *) &fileSizes[i], sizeof(uint64_t));
        }

        if (myBatches > 1) {

            // TODO: Multi-threaded n-way merge
            // Open all the files, find the boundary elements (file1[size/threadID]
            std::vector<RecordType> batchStarts(omp_get_max_threads());
            std::cout<<"MERGE Threads = "<<omp_get_max_threads()<<std::endl;
            RecordType k;
            batchStarts[0] = k;
            for (auto threadID = 1; threadID < omp_get_max_threads(); ++threadID) {
                RecordType k;
                size_t numRecs = (fileSizes[0]/(omp_get_max_threads()-threadID))* sizeof(RecordType);
                finalMerge[0].seekg( sizeof(uint64_t) + numRecs,
                                     std::ios_base::beg);
                finalMerge[0].read((char *) &k, sizeof(RecordType));
                batchStarts[threadID]=k;
            }

#pragma omp parallel
            {
                // Merge the subsets using different threads
                // Join all the merged kmers into elements ordered by threadID
                std::vector<std::ifstream> threadMerge(myBatches);
                for (auto batchID = 0; batchID < myBatches; ++batchID) {
                    threadMerge[batchID].open(tmp + "_batch_" + std::to_string(batchID + 1)
                                              + ".tmc", std::ios_base::in | std::ios_base::binary);
                    threadMerge[batchID].seekg(sizeof(uint64_t));
                }
                // Place the file pointer at the first record this thread needs to merge
                for (auto i = 0; i < myBatches; ++i) {
                    fileBinarySearch(threadMerge[i], batchStarts[omp_get_thread_num()],
                                     0, fileSizes[i]);
                }
                std::vector<RecordType> thread_elements;
                priorityQueueMerge(threadMerge, thread_elements);

#pragma omp critical (SMR_FinalMerge)
                {
                    nKmers += thread_elements.size();
                    elements = (RecordType*) realloc(elements, nKmers * sizeof(RecordType));
                    memcpy(elements, thread_elements.data(), thread_elements.size()*sizeof(RecordType));
                }
            }

            for (auto batchID = 0; batchID < myBatches; batchID++) {
                std::remove(std::string(tmp + "_batch_" + std::to_string(batchID + 1) + ".tmc").data());
            }

            std::ofstream threadMerged(tmp + ".kc",
                                       std::ios_base::trunc | std::ios_base::out | std::ios_base::binary);
            threadMerged.write((char *) &nKmers, sizeof(nKmers));
            threadMerged.write((char *) elements, sizeof(RecordType) * nKmers);
        } else {
            // Just one batch so simply rename the _batch_ file
            std::rename(std::string(tmp + "_batch_" + std::to_string(1) + ".tmc").data(), std::string(tmp + ".kc").data());
        }

    }

    bool fileBinarySearch(std::ifstream &file, RecordType &mVal, size_t lo, size_t hi) {
        file.seekg(sizeof(uint64_t), std::ios_base::beg);
        RecordType k;
        if (k == mVal) {
            return true;
        }

        while (lo <= hi and hi > 0) {
            size_t mid = (size_t) floor((lo + (hi - lo) / 2.));
            file.seekg(mid * sizeof(RecordType)+ sizeof(uint64_t), std::ios_base::beg);
            RecordType k;
            file.read((char *) &k, sizeof(RecordType ));
            if (k == mVal) {
                return true;
            }
            k<mVal ? (lo = mid+1) : (hi = mid - 1);
        }
        return false;
    }

    void priorityQueueMerge(std::vector<std::ifstream> &finalMerge, std::vector<RecordType> &_elements) {
        typedef std::pair<int, RecordType> fileRecPair;
        auto pCmp = [](const fileRecPair &left, const fileRecPair &right) { return (left.second > right.second); };

        std::priority_queue<fileRecPair, std::vector<fileRecPair>, decltype(pCmp)> fQueue(pCmp);
        for (auto i = 0; i < myBatches; ++i) {
            RecordType k;
            finalMerge[i].read((char *) &k, sizeof(RecordType));
            fQueue.emplace(i, k);
        }

        RecordType helper;
        fileRecPair current, next;
        current = fQueue.top();
        fQueue.pop();
        do {
            // Add an element from the same file to the PrtyQueue
            finalMerge[current.first].read((char *) &helper, sizeof(RecordType));
            if (finalMerge[current.first].gcount() == sizeof(RecordType))
                fQueue.emplace(current.first, helper);

            next = fQueue.top();
            fQueue.pop();
            // While every next element is equal to the current element, merge them
            while (current.second == next.second and !fQueue.empty()) {
                current.second.merge(next.second);
                finalMerge[next.first].read((char *) &helper, sizeof(RecordType));
                if (finalMerge[next.first].gcount() == sizeof(RecordType))
                    fQueue.emplace(next.first, helper);
                next = fQueue.top();
                fQueue.pop();
            }
            // Add the element to the final list
            _elements.emplace_back(current.second);
            swap(current, next);

        } while (!fQueue.empty());
    }


    void mapElementsToBatches(FileReader &myFileReader) {
#pragma omp parallel reduction(+: nRecs, tKmers, nKmers)
        {
            std::vector<RecordType> _elements;
            _elements.reserve(numMersPerThread);
            RecordFactory myRecordFactory(parameters);
            FileRecord frecord;
            // Make the next_record function get a chunk of elements
            // The chunk of elements could be implemented as a double buffer
            // Divide the chunk onto threads
            // Generate mers from each chunk and dump batches from each thread

            // This gets the filerecord metadata and keeps track of progress over the file
            while (not myFileReader.done()) {
                RecordType record;
#pragma omp critical (SMR_NextFileRecord)
                {
                    myFileReader.next_record(frecord);
                }
                nRecs++;
                // A bunch of records should be sent to each thread
                // MyRecordFactory should be threadable
                // would make setFileRecord a "process" batch
                myRecordFactory.setFileRecord(frecord);

                while (myRecordFactory.next_element(record)) {
                    tKmers++;
                    _elements.emplace_back(record);
                    if (_elements.size() >= numMersPerThread) {
                        dumpBatch(_elements);
                    }
                }
            }

            if (_elements.size() > 0) {
                if (myBatches > 1) {
                    dumpBatch(_elements);
                } else {
                    collapse(_elements);
#pragma omp critical (mergeBatchesInMem)
                    {
                        nKmers += _elements.size();
                        elements = (RecordType*) realloc(elements, nKmers * sizeof(RecordType));
                        memcpy(elements, _elements.data(), _elements.size()*sizeof(RecordType));
                    }
                }
            }
        };
        if (myBatches < 1) {
            std::vector<RecordType> _t(elements, elements + sizeof elements / sizeof elements[0]);
            collapse(_t);
            elements = realloc(elements, nKmers* sizeof(RecordType));
            memcpy(elements, _t.data(), _t.size()*sizeof(RecordType));
        }

    }

    void dumpBatch(std::vector<RecordType> &_elements) {
        collapse(_elements);
        std::ofstream of;
#pragma omp critical (dumpBatchOutput)
        {
            std::cout << "Dumping batch #" << myBatches << "\n";
            myBatches++;
            of.open(std::string(tmp + "_batch_" + std::to_string(myBatches) + ".tmc"),
                    std::ios::out | std::ios::trunc | std::ios::binary);
        }
        uint64_t size(_elements.size());

        // TODO: Make this write async
        of.write((const char *) &size, sizeof(size));
        of.write((const char *) _elements.data(), sizeof(RecordType)*size);
        _elements.clear();
    }

    void collapse(std::vector<RecordType> &_elements){

        std::sort(_elements.begin(), _elements.end());
        typename std::vector<RecordType>::iterator writePtr;
        writePtr = _elements.begin();
        typename std::vector<RecordType>::const_iterator endPtr;
        endPtr = _elements.end();

        // Keep the total in the first equal KMer
        typename std::vector<RecordType>::iterator readPtr;
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
        _elements.erase(writePtr, _elements.end());
    };

    ParamStruct parameters;
    RecordType *elements;
    uint64_t numMersPerThread;
    std::atomic<uint64_t> myBatches;
    uint64_t nRecs;
    uint64_t nKmers;
    uint64_t tKmers;
    const std::string tmp2;
    const int maxThreads;
    std::string tmp;
};


template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ParamStruct >
class SMR {
public:
    SMR(ParamStruct p, uint64_t maxMem, const std::string Otmp = "./") :
            parameters(p),
            numMersPerThread((maxMem/sizeof(RecordType)) / 2),
            myBatches(0),
            nRecs(0),
            tKmers(0),
            tmp(Otmp)
    {};

    RecordType* getRecords(FileReader &myFileReader, const std::string &name){
        RecordFactory myRecordFactory(parameters);
        FileRecord frecord;
        std::vector<RecordType> bElements;
        // This gets the filerecord metadata and keeps track of progress over the file
        while (not myFileReader.done()) {
            RecordType record;
            myFileReader.next_record(frecord);
            myRecordFactory.setFileRecord(frecord);
            nRecs++;
            while (myRecordFactory.next_element(record)) {
                tKmers++;
                bElements.push_back(record);
                if (bElements.size() >= numMersPerThread) {
                    dumpBatch(bElements);
                }
            }
        }
        if (bElements.size() > 0) {
            dumpBatch(bElements);
        }

        std::vector<std::ifstream*> finalMerge(myBatches);
        std::ofstream threadMerged(tmp+name+".kc", std::ios::trunc | std::ios::out | std::ios::binary);
        for (auto batchID = 0; batchID < myBatches; ++batchID) {
            finalMerge[batchID] = new std::ifstream (tmp+"_batch_" + std::to_string(batchID+1)
                                                     + ".tmc", std::ios::in | std::ios::binary);
        }
        std::vector<uint64_t > fileSizes(myBatches);
        for (auto i=0; i < myBatches; ++i) {
            finalMerge[i]->read((char *) &fileSizes[i], sizeof(uint64_t));
        }

        std::vector<std::pair<unsigned int, RecordType>> kCursor;
        for (auto i=0; i < myBatches; ++i) {
            RecordType k;
            finalMerge[i]->read((char *) &k, sizeof(RecordType));
            kCursor.push_back(std::make_pair(i, k));
        }

        std::vector<bool> finished(finalMerge.size(), false);
        uint64_t nKmers(numMersPerThread);
        elements = (RecordType *) malloc(nKmers * sizeof(RecordType));
        uint64_t countKmers(0);
        do {
            RecordType k;
            std::pair<unsigned int, RecordType> minFile;
            minFile.second = kCursor.begin()->second;
            for (auto itr = kCursor.begin(); itr != kCursor.cend(); ++itr) {
                if (itr->second < minFile.second and !finished[itr->first]) {
                    minFile.first = itr->first;
                    minFile.second = itr->second;
                }
            }
            auto currentMer(minFile.second);
            for (auto itr=kCursor.cbegin(); itr != kCursor.cend(); ++itr){
                if (currentMer == itr->second and !finished[itr->first]) {
                    k.merge(itr->second);
                    finalMerge[itr->first]->read((char *) &kCursor[itr->first].second, sizeof(RecordType));
                    if (finalMerge[itr->first]->eof()) {
                        finished[itr->first] = true;
                    }
                }
            }

            k = currentMer;
            if (countKmers >= nKmers) {
                nKmers = (uint64_t) (countKmers * 1.2);
                elements = (RecordType *) realloc(elements, nKmers * sizeof(RecordType));
            }
            elements[countKmers] = k;
            ++countKmers;

        } while (std::count_if(finished.begin(), finished.end(), [&](const bool &file) -> bool { return !file;})>0);

        for (auto batchID = 0; batchID<myBatches; batchID++) {
            std::remove(std::string(tmp+"_batch_" + std::to_string(batchID+1)+ ".tmc").data());
        }
        threadMerged.write((char *) &countKmers, sizeof(countKmers));
        threadMerged.write((char *) elements, sizeof(RecordType)*countKmers);
        std::cout << "Total Kmers " << tKmers << "\n";
        std::cout << "Final nKmers " << countKmers << "\n";
        std::cout << "Final reads " << nRecs << "\n";
        return elements;
    }

    RecordType * read_from_file(std::string file1, std::string file2) {
        myBatches = 0;
        nRecs = 0;
        tKmers = 0;
        FileReader myFileReader(file1, file2);
        std::replace(file1.begin(), file1.end(), '/', '.');
        std::replace(file2.begin(), file2.end(), '/', '.');
        return getRecords(myFileReader, file1+"_"+file2);
    };

    RecordType * read_from_file(std::string file) {
        myBatches = 0;
        nRecs = 0;
        tKmers = 0;
        FileReader myFileReader(file);
        std::replace(file.begin(), file.end(), '/', '.');
        return getRecords(myFileReader, file);
    };


private:

    void dumpBatch(std::vector<RecordType> &bElements) {
        myBatches++;
        collapse(bElements);
        std::ofstream of(std::string(tmp+"_batch_" + std::to_string(myBatches) + ".tmc"),
                         std::ios::out | std::ios::trunc | std::ios::binary);
        uint64_t size(bElements.size());
        of.write((const char *) &size, sizeof(size));
        of.write((const char *) bElements.data(), sizeof(RecordType)*size);
        bElements.clear();
    }
    void collapse(std::vector<RecordType> &bElements){

        __gnu_parallel::sort(bElements.begin(), bElements.end());
        typename std::vector<RecordType>::iterator writePtr = bElements.begin();
        typename std::vector<RecordType>::const_iterator endPtr = bElements.end();

        // Keep the total in the first equal KMer
        for (typename std::vector<RecordType>::iterator readPtr = bElements.begin(); readPtr != endPtr; ++writePtr) {
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
        bElements.erase(writePtr, bElements.end());
    };

    ParamStruct parameters;
    RecordType* elements;
    uint64_t numMersPerThread;
    uint64_t myBatches;
    uint64_t nRecs;
    uint64_t tKmers;
    const std::string &tmp;
};

#endif //W2RAP_CONTIGGER_SMR_H
