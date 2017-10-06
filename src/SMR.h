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
    FastQReader(std::string filepath){
        std::cout << "Opening: " << filepath << "\n";
        read.open(filepath.data());
    }

    void next_record(FileRecord& rec){
        std::string dummy;
        std::getline(read, record.name);
        std::getline(read, record.seq);
        std::getline(read, dummy);
        std::getline(read, record.qual);
        std::swap(record,rec);
    }

    bool done() {
        return read.eof();
    }

private:
    std::ifstream read;
    FileRecord record;
};

#include <vector>
#include <iostream>
#include <algorithm>

template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ParamStruct >
class SMR {
public:
    SMR(ParamStruct p, uint64_t maxMem, const std::string Otmp = "") :
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
        return getRecords(myFileReader, "/"+file);
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
