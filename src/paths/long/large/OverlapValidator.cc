//
// Created by Bernardo Clavijo (TGAC) on 12/12/2016.
//

#include <assert.h>
#include "OverlapValidator.h"

InformativePair::InformativePair(uint64_t index, ReadPath &r1Path, ReadPath &r2Path, vec<int> &inv) {
    //copy paths and generate inverse paths
    //std::cout<<r1Path.size()<<" "<<r2Path.size()<<std::endl;
    pathIndex=index;
    for (auto &e:r1Path) {
        uint64_t eid = e;
        r1path.push_back(eid);
    }
    for(auto eid=r1path.rbegin();eid!=r1path.rend();++eid){
        if (inv[*eid]>0) {
            uint64_t rid=inv[*eid];
            r1rpath.push_back(rid);
        }
        else r1rpath.push_back(*eid);
    }
    int64_t overlap=-1;
    for (auto &e:r2Path){
        uint64_t eid=e;
        r2path.push_back(eid);
    }
    for(auto eid=r2path.rbegin();eid!=r2path.rend();++eid){
        if (inv[*eid]>0) {
            uint64_t rid=inv[*eid];
            r2rpath.push_back(rid);
        }
        else r2rpath.push_back(*eid);
        if (r1path.size()>0 and r2rpath.back()==r1path.back()) overlap=r2rpath.size()-1;
    }
    //evaluate overlap and combine if possible
    if (overlap!=-1){
        combined_path=r1path;
        combined_path.pop_back(); //avoid the dupplication of edge
        for (auto e=r2rpath.begin()+overlap; e!=r2rpath.end();++e) combined_path.push_back(*e);
        for(auto eid=combined_path.rbegin();eid!=combined_path.rend();++eid){
            if (inv[*eid]>0) {
                uint64_t rid=inv[*eid];
                combined_rpath.push_back(rid);
            }
            else combined_rpath.push_back(*eid);
        }
        //std::cout<<"combined into single path of size"<<combined_path.size()<<std::endl;
    }

}

std::vector<uint64_t> InformativePair::overlaps_crossed(vec<int> & toLeft, vec<int> &toRight) {
    //returns vertex indexes on the argument HBV crossed by the pair
    std::vector<uint64_t> crossed;
    std::vector<std::vector<uint64_t>> paths_to_analyze;
    //if combined path only consider that (forward and reverse)
    if (is_combined()){
        paths_to_analyze.push_back(combined_path);
        paths_to_analyze.push_back(combined_rpath);
    }
    //else consider each reads path on its own (forward and reverse)
    else{
        paths_to_analyze.push_back(r1path);
        paths_to_analyze.push_back(r1rpath);
        paths_to_analyze.push_back(r2path);
        paths_to_analyze.push_back(r2rpath);
    }
    for (auto &p:paths_to_analyze){
        for (int i=0;i<((int)p.size())-1;++i){
            if (toRight[p[i]]==toLeft[p[i+1]]){
                uint64_t vid=toRight[p[i]];
                crossed.push_back(vid);
            }
        }
    }
    std::sort(crossed.begin(),crossed.end());
    crossed.erase(std::unique(crossed.begin(),crossed.end()),crossed.end());
    return crossed;
}

std::vector<uint64_t> InformativePair::overlaps_jumped(vec<int> & toLeft, vec<int> &toRight, uint8_t max_dist) {
    //returns vertex indexes on the argument HBV jumped by the pair with up to max_dist intermediate overlaps
    assert(max_dist==0); //so far we are only implementing pairs direcly jumping through a single vertex
    std::vector<uint64_t> crossed;
    //if combined path there is no jumps
    if (!is_combined() && r1path.size()>0 and r2path.size()>0) {
        //else check if there is a single vertex between r1 and r2
        auto r1_end=r1path.back();
        auto r2_start=r2rpath.front();
        if (toRight[r1_end]==toLeft[r2_start]) {
            uint64_t vid=toRight[r1_end];
            crossed.push_back(vid);
        }
    }
    return crossed;
}



void OverlapValidator::compute_overlap_support() {
    find_informative_pairs();
    //Go through each vertex, consider the possible overlap of each
}

void OverlapValidator::find_informative_pairs() {
    //Go through all paths, when multi-part paths found, find the corresponding vertex and insert path number for later analysis
    //process each path as pair
    vec<int> toLeft,toRight;
    mHBV.ToLeft(toLeft);
    mHBV.ToRight(toRight);
    mCross.clear();
    mCross.resize(mHBV.N());
    mJump.clear();
    mJump.resize(mHBV.N());
    std::atomic_uint_fast64_t cross(0), jump(0), overlap(0);
    std::cout<<Date()<<": analysing "<<mPaths.size()<<" paths for information"<<std::endl;
//#pragma omp parallel for
    for(int64_t i=0;i<mPaths.size();i+=2){
        //if both ends map single-part to the same edge, discard
        if (mPaths[i].size()==1 && mPaths[i+1].size()==1 and mPaths[i][0]==mInv[mPaths[i+1][0]]) continue;

        InformativePair new_path(i,mPaths[i],mPaths[i+1],mInv);
        if (new_path.is_combined()) overlap++;
        uint64_t pid;
#pragma omp critical(OVLVAL_ipaths)
        {
            mInformativePairs.push_back(new_path);
            pid=mInformativePairs.size()-1;
        }
        for (auto &oj:new_path.overlaps_crossed(toLeft,toRight)){
#pragma omp critical(OVLVAL_cross)
            mCross[oj].push_back(pid);
            ++cross;
        }
        for (auto &oj:new_path.overlaps_jumped(toLeft,toRight,0)){
#pragma omp critical(OVLVAL_jump)
            mJump[oj].push_back(pid);
            ++jump;
        }

    }
    std::cout<<Date()<<": "<<mInformativePairs.size()<<" informative pairs ( "<<overlap<<" combined )"<<std::endl;
    uint64_t vcross=0,vjump=0,unsupported=0,complex=0;

    for (auto &vc:mCross) if (vc.size()>0) vcross++;
    for (auto &vj:mJump) if (vj.size()>0) vjump++;
    for (uint64_t i=0;i<mCross.size();++i){
        if (mCross[i].size()==0 and mJump[i].size()==0) ++unsupported;
        else {
            if (mHBV.From(i).size()>1 and mHBV.From(i).size()>1) complex++;
        }
    }
    std::cout<<Date()<<": "<<vcross<<" vertices crossed by a total of "<<cross<<" pairs"<<std::endl;
    std::cout<<Date()<<": "<<vjump<<" vertices jumped by a total of "<<jump<<" pairs"<<std::endl;
    std::cout<<Date()<<": "<<complex<<" complex vertices with support"<<std::endl;
    std::cout<<Date()<<": "<<unsupported<<" vertices unsupported"<<std::endl;

}

