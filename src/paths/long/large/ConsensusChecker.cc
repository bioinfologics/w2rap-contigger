//
// Created by Bernardo Clavijo (TGAC) on 23/12/2016.
//

#include "ConsensusChecker.h"

bool ConsensusChecker::consensus_OK(uint64_t edge) {
    collect_edge_calls(edge);
    collect_read_calls(edge);


    //Base-By-Base, check how many reads agree on consensus, if less than 80%, flag as problematic
    bool OK=true;
    for (auto i=0;i<edge_call.size();++i){
        if (read_calls[i].total()>9 and read_calls[i].support(edge_call[i])<.7) {
            //std::cout<<"consensus not strong enough on position "<<i<<std::endl;
            OK=false;
            break;
        }
    }


    return OK;
}

void ConsensusChecker::print_detail(){
    std::cout<<"Position, Call, A, C, G, T, N, correct %"<<std::endl;
    char b[5]={'A','C','G','T','N'};
    for (auto i=0;i<edge_call.size();++i) {
        if (read_calls[i].total()>9 and read_calls[i].support(edge_call[i])<.7)
        std::cout<<i<<", "<<b[edge_call[i]]<<", "<<read_calls[i].A<<", "<<read_calls[i].C<<", "<<read_calls[i].G<<", "<<read_calls[i].T<<", "<<read_calls[i].N<<","<<read_calls[i].support(edge_call[i])<<std::endl;
    }
}

void ConsensusChecker::collect_read_calls(uint64_t edge) {
    read_calls.clear();
    read_calls.resize(mHBV.EdgeObject(edge).size(),{0,0,0,0,0});
    //Condense all reads placed into the edges
    auto invedge=mInv[edge];
    uint64_t start_pos;
    std::vector<char> read_trimmed_seq;
    int64_t pid=-1;
    for (auto p:mPaths) {//TODO: path inversion (i.e. dict from path to edge would be much more useful and fast here
        pid++;
        bool used=false;
        uint edge_pos=p.getOffset();
        uint left_trim=0;
        for (auto e:p) {
            if (e==edge) { //TODO: check rc and rc sequences, etc
                used = true;
                break;

            }
            left_trim+=mHBV.EdgeObject(e).size()-edge_pos+28;//reads are hardcoded to be mapped at 60-mers.
            edge_pos=0;
        }
        if (not used) continue;

        if (left_trim) continue;//just to make sure
        //
        std::string readseq=mBases[pid].ToString();

        //std::cout<<">path_"<<pid<<"_pos_"<<edge_pos<<std::endl<<readseq<<std::endl;
        //now add the read's consensus to read_calls.
        auto cp=edge_pos;
        for (auto i=0;i<readseq.size()-left_trim and cp<read_calls.size();++i,++cp){
            switch (readseq[left_trim+i]){
                case 'A':
                    read_calls[cp].A++;
                    break;
                case 'C':
                    read_calls[cp].C++;
                    break;
                case 'G':
                    read_calls[cp].G++;
                    break;
                case 'T':
                    read_calls[cp].T++;
                    break;
                default:
                    read_calls[cp].N++;
                    break;

            }

        }



    }

}

void ConsensusChecker::collect_edge_calls(uint64_t edge) {
    edge_call.clear();
    auto s=mHBV.EdgeObject(edge).ToString();
    //std::cout<<">edge"<<edge<<std::endl<<s<<std::endl;
    edge_call.reserve(s.size());
    for (auto c:s){
        switch (c){
            case 'A':
                edge_call.push_back(0);
                break;
            case 'C':
                edge_call.push_back(1);
                break;
            case 'G':
                edge_call.push_back(2);
                break;
            case 'T':
                edge_call.push_back(3);
                break;
            default:
                edge_call.push_back(4);
                break;
        }
    }

}
