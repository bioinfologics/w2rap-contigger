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

bool InformativePair::crosses_transition(uint64_t e1, uint64_t e2) {
    if (is_combined()){
        for (auto i=0;i<combined_path.size()-1;++i) if (combined_path[i]==e1 and combined_path[i+1]==e2) return true;
        for (auto i=0;i<combined_rpath.size()-1;++i) if (combined_rpath[i]==e1 and combined_rpath[i+1]==e2) return true;
    } else {
        if (r1path.size()>1) {
            for (auto i = 0; i < r1path.size() - 1; ++i) if (r1path[i] == e1 and r1path[i + 1] == e2) return true;
            for (auto i = 0; i < r1rpath.size() - 1; ++i) if (r1rpath[i] == e1 and r1rpath[i + 1] == e2) return true;
        }
        if (r2path.size()>1) {
            for (auto i = 0; i < r2path.size() - 1; ++i) if (r2path[i] == e1 and r2path[i + 1] == e2) return true;
            for (auto i = 0; i < r2rpath.size() - 1; ++i) if (r2rpath[i] == e1 and r2rpath[i + 1] == e2) return true;
        }
    }
    return false;
}

bool InformativePair::jumps_transition(uint64_t e1, uint64_t e2) {
    if (!is_combined() and r1path.size()>0 and r2path.size()>0){
        if (r1path.back()==e1 and r2rpath.front()==e2) return true;
        if (r2path.back()==e1 and r1rpath.front()==e2) return true;
    }
    return false;
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
    uint64_t vcross=0,vjump=0,unsupported=0,complex=0,in_eq_out=0;

    for (auto &vc:mCross) if (vc.size()>0) vcross++;
    for (auto &vj:mJump) if (vj.size()>0) vjump++;
    for (uint64_t i=0;i<mCross.size();++i){
        if (mCross[i].size()==0 and mJump[i].size()==0) ++unsupported;
        else {
            if (mHBV.From(i).size()>1 and mHBV.To(i).size()>1) {
                complex++;
                if (mHBV.From(i).size() == mHBV.To(i).size()) in_eq_out++;
            }
        }
    }
    std::cout<<Date()<<": "<<vcross<<" vertices crossed by a total of "<<cross<<" pairs"<<std::endl;
    std::cout<<Date()<<": "<<vjump<<" vertices jumped by a total of "<<jump<<" pairs"<<std::endl;
    std::cout<<Date()<<": "<<complex<<" complex vertices with support ( "<<in_eq_out<<" with #in=#out )"<<std::endl;
    std::cout<<Date()<<": "<<unsupported<<" vertices unsupported"<<std::endl;

}

void OverlapValidator::analyse_complex_overlaps() {
    struct transition_support{
        uint64_t e1;
        uint64_t e2;
        uint64_t cross;
        uint64_t jump;
    };
    uint64_t to_expand=0;
    for (uint64_t vi=0;vi<mCross.size();++vi){
        if (mHBV.From(vi).size()>1 and mHBV.To(vi).size()>1){
            //std::cout<<" analyssing complex vertex "<<vi<<std::endl;
            std::vector<struct transition_support> transitions;//generate all possible combinations of in-outs
            std::vector<struct transition_support> valid_transitions;
            //std::cout<<" generating transitions "<<vi<<std::endl;
            for (auto i1=0;i1<mHBV.To(vi).size();++i1)
                for (auto i2=0;i2<mHBV.From(vi).size();++i2){
                    struct transition_support ts;
                    ts.e1=mHBV.EdgeObjectIndexByIndexTo(vi,i1);
                    ts.e2=mHBV.EdgeObjectIndexByIndexFrom(vi,i2);
                    ts.cross=0;
                    ts.jump=0;
                    //std::cout<<ts.e1<<"->"<<ts.e2<<std::endl;
                    for (auto &p:mCross[vi]){
                        //std::cout<<"analysing crossing for informative pair "<<p<<std::endl;
                        if (mInformativePairs[p].crosses_transition(ts.e1,ts.e2)) ++ts.cross;
                        //std::cout<<"done!!!"<<std::endl;
                    }
                    //std::cout<<"cross: "<<ts.cross<<std::endl;
                    for (auto &p:mJump[vi]) if (mInformativePairs[p].jumps_transition(ts.e1,ts.e2)) ++ts.jump;
                    //std::cout<<"jump: "<<ts.jump<<std::endl;
                    transitions.push_back(ts);
                }
            //std::cout<<" testing for possible solutions "<<vi<<std::endl;
            if (mHBV.From(vi).size()==mHBV.To(vi).size()) {//simplest case, find 1-to-1s
                uint64_t total_cross = mCross[vi].size(), total_jump = mJump[vi].size(), total=total_cross+total_jump;
                //todo consider percentajes and such
                std::vector<uint64_t> used_ins,used_outs;
                bool conflict=false;
                for (auto t:transitions) {
                    if (t.cross+t.jump>1) {
                        if (std::find(used_ins.begin(),used_ins.end(),t.e1)!=used_ins.end()) conflict=true;
                        if (std::find(used_outs.begin(),used_outs.end(),t.e2)!=used_outs.end()) conflict=true;
                        used_ins.push_back(t.e1);
                        used_outs.push_back(t.e2);
                        valid_transitions.push_back(t);
                    }
                }
                std::cout<<"Vertex "<<vi<<" IN: ";
                for (auto i1=0;i1<mHBV.To(vi).size();++i1) std::cout<<mHBV.EdgeObjectIndexByIndexTo(vi,i1)<<" ";
                std::cout<<"  OUT: ";
                for (auto i1=0;i1<mHBV.From(vi).size();++i1) std::cout<<mHBV.EdgeObjectIndexByIndexFrom(vi,i1)<<" ";
                std::cout<<std::endl<<"Supported transitions: ";
                for (auto &t:valid_transitions) std::cout<<"  "<<t.e1<<"->"<<t.e2<<":"<<t.cross<<","<<t.jump;
                std::cout<<std::endl;
                if (valid_transitions.size()==mHBV.From(vi).size() and conflict==false) {
                    std::cout<<"Solved!"<<std::endl;
                    to_expand++;
                }

            }
            //std::cout<<" done "<<vi<<std::endl;


        }
    }
    std::cout<<Date()<<": "<<to_expand<<" overlaps could be trivially expanded"<<std::endl;
}



std::vector<uint64_t> OverlapValidator::find_perfect_tips(uint16_t max_size) {
    std::set<uint64_t> tip1, tip2;
    std::vector<uint64_t> tips;
    //Find vertices with only one input
    for (auto vi=0;vi<mHBV.N();++vi) if (mHBV.To(vi).size()==1 and mHBV.From(vi).size()==0) {
            //check edges and validate their lack of support on the output transition
            uint64_t tip=mHBV.EdgeObjectIndexByIndexTo(vi,0);
            auto vfork=mHBV.To(vi)[0];
            //check Y topology
            if (mHBV.From(vfork).size()==2 and mHBV.To(vfork).size()==1) {
                auto prev=mHBV.EdgeObjectIndexByIndexTo(vfork,0);
                auto other=mHBV.EdgeObjectIndexByIndexFrom(vfork,0);
                if (other==tip) other=mHBV.EdgeObjectIndexByIndexFrom(vfork,1);

                if (collect_all_support(vi,prev,tip)*10<collect_all_support(vi,prev,other))tip1.insert(tip);
                //add to tip1

            }
        }
    //Find vertices with only one output
    for (auto vi=0;vi<mHBV.N();++vi) if (mHBV.From(vi).size()==1 and mHBV.To(vi).size()==0) {
            //check edges and validate their lack of support on the input transition
            uint64_t tip=mHBV.EdgeObjectIndexByIndexFrom(vi,0);
            auto vfork=mHBV.From(vi)[0];
            //check Y topology
            if (mHBV.To(vfork).size()==2 and mHBV.From(vfork).size()==1) {
                auto next=mHBV.EdgeObjectIndexByIndexFrom(vfork,0);
                auto other=mHBV.EdgeObjectIndexByIndexTo(vfork,0);
                if (other==tip) other=mHBV.EdgeObjectIndexByIndexTo(vfork,1);

                if (collect_all_support(vi,tip,next)*10<collect_all_support(vi,other,next))tip2.insert(tip);

            }
        }

    //edge on tip1, inverse on tip2, size >max_size -> ADD BOTH.
    for (auto &tip:tip1) if (tip2.find(mInv[tip])!=tip2.end() and mHBV.EdgeObject(tip).size()<=max_size) tips.push_back(tip);
    return tips;
}

uint64_t OverlapValidator::collect_all_support(uint64_t vi, uint64_t e1, uint64_t e2) {
    uint64_t cross=0,jump=0;
    for (auto &p:mCross[vi]) if (mInformativePairs[p].crosses_transition(e1,e2)) ++cross;
    for (auto &p:mJump[vi]) if (mInformativePairs[p].jumps_transition(e1,e2)) ++jump;
    return cross+jump;
}
