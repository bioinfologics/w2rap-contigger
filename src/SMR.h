//
// Created by Luis Yanes (EI) on 24/07/2017.
//

#ifndef W2RAP_CONTIGGER_SMR_H
#define W2RAP_CONTIGGER_SMR_H

#include <vector>
#include <iostream>
#include <algorithm>

#include <fcntl.h>
#include <stdlib.h>
#include <malloc.h>

#ifndef O_DIRECT
#define O_DIRECT O_RDONLY
#endif

std::istream& operator>> (std::istream& is, KMerNodeFreq_s& kmer) {
    is.read((char*)&kmer, sizeof(kmer));
    return is;
}

std::ostream& operator<<(std::ostream& os, const KMerNodeFreq_s& kmer) {
    os << kmer.kdata[0] << '-' << kmer.kdata[1] << "\n";
    os << kmer.count << ' ' << kmer.kc << "\n";
    return os;
}

template<typename FileRecord>
class FastBReader {
public:
    FastBReader(const vecbvec &_reads, const size_t _readsBegin, const size_t _readsEnd): reads(_reads), readsBegin(_readsBegin), readsEnd(_readsEnd) { pos = readsBegin; }
    bool done() {return pos >= readsEnd;}
    bool next_record(FileRecord &rec){
        if (pos < readsEnd){
            rec = reads[pos];
            pos++;
            return true;
        }
        return false;
    }
private:
    const vecbvec &reads;
    const size_t readsBegin, readsEnd;
    size_t pos;
};

template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ParamStruct >
class SMR2 {
public:
    SMR2(ParamStruct p, uint64_t maxMem, const std::string &Otmp = "") :
            parameters(p),
            myBatches(0),
            nRecs(0),
            nKmers(0),
            tKmers(0),
            tmp2(Otmp),
            maxThreads((unsigned int) omp_get_max_threads())
    {
        numElementsPerBatch = maxMem/sizeof(RecordType) / maxThreads;
        mergeCount=2;
    }

    void read_from_file(const vecbvec &reads, int numThreads) {
        uint64_t numReadsReduction(0), numKmersReduction(0);
        OutputLog(3) << "Begin parallel reduction of " << reads.size() << " reads to batches\n";
#pragma omp parallel reduction(+: numReadsReduction, numKmersReduction)
        {
            size_t from, to;
            from = (reads.size()/omp_get_max_threads()) * (omp_get_thread_num());
            to = (reads.size()/omp_get_max_threads()) * (omp_get_thread_num()+1);

            if(omp_get_thread_num()+1 == omp_get_num_threads()) to = reads.size();

            FileReader myFileReader(reads, from, to);
#pragma omp critical (mapStart)
            {
                std::cout << "Thread = " << omp_get_thread_num() << " " << from << " - " << to << std::endl;
            }
            mapElementsToBatches2(myFileReader, numReadsReduction, numKmersReduction);
        }
        OutputLog(3) << "Done parallel reduction to batches\n";
        OutputLog(3) << "Reduction Reads = " << std::to_string(numReadsReduction) << std::endl;
        OutputLog(3) << "Reduction Kmers = " << std::to_string(numKmersReduction) << std::endl;
        this->nRecs = numReadsReduction;
    };

    void getRecords(RecordType *records){
        std::swap(elements,records);
        reduceElementsFromBatches2();
        std::swap(records,elements);

        std::cout << "Total Kmers " << tKmers << "\n";
        std::cout << "Final nKmers " << nKmers << "\n";
        std::cout << "Final reads " << nRecs << "\n";
    }

    void reduceElementsFromBatches2() {
        OutputLog(3) << "Merging thread files" << std::endl;
        std::vector<std::ifstream> threadFiles(maxThreads);
        for (auto threadID = 0; threadID<maxThreads; threadID++) {
            std::string name(tmp + "thread_" + std::to_string(threadID) + ".tmc");
            //std::cout << "Opening file " << name << std::endl;
            threadFiles[threadID].open(name.data(), std::ios_base::in|std::ios_base::binary);
        }
        std::vector<RecordType> elements;
        nKmers = merge("final.kc", threadFiles);
        OutputLog(3) << "Done merging thread files" << std::endl;
    }

private:
    typedef std::pair<int, RecordType> fileRecPair;

    struct FileRecPairMoreCmp
    {
        bool operator()(const fileRecPair& lhs, const fileRecPair& rhs)
        {
            return lhs.second > rhs.second;
        }
    };

    typedef std::priority_queue<fileRecPair, std::vector<fileRecPair>, FileRecPairMoreCmp> fileRecPQ;

    bool addToQueue(std::vector<std::ifstream> &inf, fileRecPQ &q, const fileRecPair &pair) const
    {
        RecordType helper;
        inf[pair.first].read((char *) &helper, sizeof(RecordType));
        if (inf[pair.first].gcount() == sizeof(RecordType)) {
//            if (isEqual(helper)) std::cerr << "Read from file " << pair.first << " on thread " << omp_get_thread_num() << std::endl;
            q.emplace(pair.first, helper);
            return true;
        }
        return false;
    }

    void outputElement(std::ofstream &outf, const RecordType &smallest) const {
        outf.write((char*) &smallest, sizeof(RecordType));
    }

    bool isEqual(const RecordType &rec) const {
        return false;
        uint64_t a(18446744072652575747LLU);
        uint64_t b(13835058055282163712LLU);
        KMerNodeFreq_s tk;
        tk.kdata[0] = a;
        tk.kdata[1] = b;
        if (rec == tk) {
            return true;
        }
        return false;
    }

    uint64_t merge(const std::string &tmpName, std::vector<std::ifstream> &inf) {
        //std::cout << "Begin merge " << tmpName << "\n";
        std::ofstream outf(tmpName.data(), std::ios_base::binary|std::ios_base::out|std::ios_base::trunc);
        std::vector<RecordType> outvec;
        uint64_t numElements(0);
        outf.write((char*)&numElements, sizeof(uint64_t));

        RecordType helper;

        fileRecPQ fQueue;
        for (size_t i = 0; i < inf.size(); ++i) {
            uint64_t tmp;
            inf[i].read((char *) &tmp, sizeof(uint64_t));
            inf[i].read((char *) &helper, sizeof(RecordType));
            fQueue.emplace(i, helper);
        }
        fileRecPair smallest;
        if (!fQueue.empty()) smallest = fQueue.top();fQueue.pop();
        addToQueue(inf, fQueue, smallest);
        fileRecPair next;
        if (!fQueue.empty()) next = fQueue.top();fQueue.pop();
        while (!fQueue.empty()) {
            if (smallest.second < next.second) {
                outputElement(outf, smallest.second);
                numElements++;
                addToQueue(inf, fQueue, smallest);
                if (! (fQueue.top().second < next.second) ) {
                    smallest = next;
                    addToQueue(inf, fQueue, next);
                    next = fQueue.top();fQueue.pop();
                }
            } else if (smallest.second == next.second) {
                smallest.second.merge(next.second);
                addToQueue(inf, fQueue, next);
                next = fQueue.top();fQueue.pop();
            }
        }
        if (smallest.second < next.second) {
            outputElement(outf, smallest.second);
            numElements++;
            outputElement(outf, next.second);
            numElements++;
        }
        outf.seekp(0, std::ios_base::beg);
        outf.write((char*)&numElements, sizeof(uint64_t));

        return numElements;
    }

    uint64_t merge(const std::string &tmpName, std::vector<std::ifstream> &inf, std::vector<RecordType> &elements) {
        std::cout << "Begin merge+memory " << tmpName << "\n";
        std::ofstream outf(tmpName.data(), std::ios_base::binary|std::ios_base::out|std::ios_base::trunc);
        std::vector<RecordType> outvec;
        uint64_t numElements(0), memElement(0);
        outf.write((char*)&numElements, sizeof(uint64_t));

        RecordType helper;
        fileRecPQ fQueue;
        for (size_t i = 0; i < inf.size(); ++i) {
            uint64_t tmp;
            inf[i].read((char *) &tmp, sizeof(uint64_t));
            inf[i].read((char *) &helper, sizeof(RecordType));
            fQueue.emplace(i, helper);
        }

        fileRecPair fileSmallest = fQueue.top();fQueue.pop();
        if (memElement < elements.size()) {
            do {
                if (fileSmallest.second == elements[memElement]) {
                    // Merge
                    elements[memElement].merge(fileSmallest.second);
                    addToQueue(inf, fQueue, fileSmallest);
                    fileSmallest = fQueue.top();fQueue.pop();
                } else if (elements[memElement] < fileSmallest.second) {
                    // Output memElement
                    outputElement(outf, elements[memElement]);
                    memElement++;
                    numElements++;
                } else {
                    // Output fileSmallest
                    outputElement(outf, fileSmallest.second);
                    numElements++;
                    addToQueue(inf, fQueue, fileSmallest);
                    fileSmallest = fQueue.top();fQueue.pop();
                }
            } while (memElement < elements.size() and !fQueue.empty());
        }

        addToQueue(inf, fQueue, fileSmallest);
        if (!fQueue.empty()) {
            while (!fQueue.empty()) {
                fileRecPair next = fQueue.top();
                fQueue.pop();
                if (fileSmallest.second == next.second) {
                    fileSmallest.second.merge(next.second);
                    addToQueue(inf, fQueue, next);
                } else {
                    outputElement(outf, fileSmallest.second);
                    numElements++;
                    addToQueue(inf, fQueue, fileSmallest);
                    fileSmallest = next;
                }
            }
            outputElement(outf, fileSmallest.second);
            numElements++;
        } else {
            auto leftToWrite = elements.size()-memElement;
            if (leftToWrite > 0) {
//                std::copy(elements.begin()+memElement, elements.end(), std::ostream_iterator<RecordType>(std::cout,"\n"));
                outf.write((char*) &elements[memElement], leftToWrite*sizeof(RecordType));
                numElements+=leftToWrite;
            }
        }

        outf.seekp(0, std::ios_base::beg);
        outf.write((char*)&numElements, sizeof(uint64_t));
        std::cout << "Done with " << tmpName << std::endl;

        return numElements;
    }

    void mapElementsToBatches2(FileReader &myFileReader, uint64_t &numReadsReduction, uint64_t &numKmersReduction) {
        std::vector<RecordType> _elements;
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
                numKmersReduction++;
                if (_elements.size() >= numElementsPerBatch) {
                    tKmers+=_elements.size();
                    if (currentBatch==mergeCount) {
                        dumpBatch(currentBatch, _elements);
                        mergeBatches(currentBatch);
                        currentBatch=1;
                    } else {
                        dumpBatch(currentBatch, _elements);
                        currentBatch++;
                    }
                }
            }
        }
        tKmers+=_elements.size();
        dumpBatch(currentBatch, _elements);

        mergeBatches(currentBatch);
        rename(std::string(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_batch_0.tmc").c_str(), std::string(tmp+"thread_"+ std::to_string(omp_get_thread_num()) +".tmc").c_str());
    }

    void collapsedElementsToFile(const std::string &tmpName, std::vector<RecordType> &_elements) const {
        // Simply dump the file
        //std::cout << "Dumping file: " << tmpName << std::endl;
        std::ofstream outBatch(tmpName.data()
                    , std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        auto size = _elements.size();
        outBatch.write((char *)&size, sizeof(size));
        outBatch.write((char *)_elements.data(), size*sizeof(RecordType));
        //std::cout << "Done dumping file: " << tmpName << std::endl;
        outBatch.close();
        _elements.clear();
    }

    void dumpBatch(const uint currentBatch, std::vector<RecordType> &_elements) {
        collapse(_elements);
        // Simply dump the file
        //std::cout << "Thread " << omp_get_thread_num() << " dumping batch #" << currentBatch << "\n";
        std::string outBatchName(std::string(tmp
                                             + "thread_"
                                             + std::to_string(omp_get_thread_num())
                                             + "_batch_" + std::to_string(currentBatch)
                                             + ".tmc"));
        collapsedElementsToFile(outBatchName, _elements);
        isOrdered(outBatchName);
    }

    void mergeBatches(uint currentCount) {
        uint numFilesToMerge = currentCount+1;
        std::vector<std::ifstream> inBatches(numFilesToMerge);
        for (auto batch=0; batch < numFilesToMerge; batch++) {
            std::string filename(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_batch_" + std::to_string(batch)+ ".tmc");
            //std::cout << "Opening " << filename << std::endl;
            inBatches[batch].open(filename, std::ios_base::binary | std::ios_base::in);
        }
        std::string threadFile(tmp+ "thread_"+ std::to_string(omp_get_thread_num())+ "_tmpBatch1.tmc");
        OutputLog(3) << "Merging " << numFilesToMerge << " in thread " << omp_get_thread_num() << std::endl;
        merge(threadFile, inBatches);
        OutputLog(3) << "DONE - Merging " << numFilesToMerge << " in thread " << omp_get_thread_num() << std::endl;
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

    bool isOrdered(const std::string &name) {
        return true;
        std::ifstream inf(name.data(), std::ios_base::binary|std::ios_base::in);
        if (!inf) {
            return false;
        }
        uint64_t totalKmers(0), fileKmers(0);
        inf.read((char *) &totalKmers, sizeof(totalKmers));

#pragma omp critical
        {
            std::cout << name << " - Final count " << totalKmers << "\n";
        }
        KMerNodeFreq_s mer, prev;
        inf.read((char*) &mer, sizeof(RecordType));
        fileKmers++;
//        std::cout << mer;
        prev = mer;
        while ( inf.read((char*) &mer, sizeof(RecordType)) ) {
            if(isEqual(mer)) std::cerr << "Found it on file " << name << " at " << fileKmers << std::endl;
//            std::cout << mer;
            if (! (mer > prev)) {
                std::cerr << "Input error at " << inf.tellg() << " prev >= current\n" << "Prev: "<< prev << "Curr: " << mer;
                std::cerr << "File kmer = " << fileKmers << std::endl;
                return false;
            }
            prev = mer;
            fileKmers++;
        }

        if (totalKmers!=fileKmers) {
            std::cerr << "Total kmers (" << totalKmers << ") is different from read kmers " << fileKmers << std::endl;
            return false;
        }

        std::cout << totalKmers << " = " << fileKmers << std::endl;
        return true;
    }

    ParamStruct parameters;
    RecordType *elements;
    uint64_t numElementsPerBatch;
    std::atomic<uint64_t> myBatches;
    uint64_t nRecs;
    uint64_t nKmers;
    std::atomic<uint64_t> tKmers;
    const std::string tmp2;
    const int unsigned maxThreads;
    int mergeCount; // How many batches to keep rolling before merging
    std::string tmp;
};

#endif //W2RAP_CONTIGGER_SMR_H
