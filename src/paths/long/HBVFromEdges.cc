///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2016) by the     //
//  Earlham Institute.  All rights are reserved.  This software is supplied  //
//   without any warranty or guaranteed support whatsoever. The Earlham      //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * HBVFromEdges.cc
 *
 *  Created on: Sept 26, 2016
 *      Author: bjclavijo/kebarr
 */

#include "paths/long/HBVFromEdges.h"
#include "dna/CanonicalForm.h"
#include "paths/KmerPathInterval.h"
#include "system/ID.h"
#include <list>



class EdgeEnd{
public:
    std::string seq;
    bool distal,rc;
    uint64_t edge_index;
    EdgeEnd(const uint64_t _edge_index, const bvec & bvseq,const bool &_rc,const bool &_distal, const unsigned K) {//K is the HBV's k, overlaps are k-1
        edge_index=_edge_index;
        distal=_distal;
        rc=_rc;
        auto b=bvseq;
        if (rc) b.ReverseComplement();
        seq=b.ToString();
        if (_distal)
            seq=seq.substr(seq.length()-K+1,seq.length());
        else
            seq.resize(K-1);
    };
    bool operator<(const EdgeEnd & other){return seq<other.seq;};//WARNING, for simplicity
    bool operator==(const EdgeEnd & other){return seq==other.seq;};//WARNING, for simplicity

};

typedef struct {
    int64_t fw_v1,fw_v2,rc_v1,rc_v2;
} edge_vertex_list;


void buildHBVFromEdges( vecbvec const& edges, unsigned K, HyperBasevector* pHBV,
                            std::vector<int> &pFwdEdgeXlat, std::vector<int> &pRevEdgeXlat )
{
    pHBV->Clear();
    pFwdEdgeXlat.clear();
    pRevEdgeXlat.clear();
    if ( edges.empty() )
    {
        return;
    }

    std::vector<EdgeEnd> ends;
    ends.reserve(2*edges.size());

    for (uint64_t i=0;i<edges.size();++i){
        ends.push_back(EdgeEnd(i,edges[i],false,false,K));
        ends.push_back(EdgeEnd(i,edges[i],false,true,K));
        if ( edges[i].getCanonicalForm() != CanonicalForm::PALINDROME ) {
            ends.push_back(EdgeEnd(i,edges[i],true,false,K));
            ends.push_back(EdgeEnd(i,edges[i],true,true,K));
        }
    }
    std::sort(ends.begin(),ends.end());

    //number the vertices and fill a list of vertex numbers for the edges
    std::vector<edge_vertex_list> edge_vertices;
    edge_vertices.resize(edges.size(),{-1,-1,-1,-1});
    uint64_t vID=0;
    //Can't be parallel, because vID depends on previous iteration
    for (auto i=0; i<ends.size(); ++i){
        if (i>0 and not (ends[i-1]==ends[i])) vID++;
        if (!ends[i].rc) {
            if (!ends[i].distal) edge_vertices[ends[i].edge_index].fw_v1=vID;
            else edge_vertices[ends[i].edge_index].fw_v2=vID;
        } else {
            if (!ends[i].distal) edge_vertices[ends[i].edge_index].rc_v1=vID;
            else edge_vertices[ends[i].edge_index].rc_v2=vID;
        }
    }
    ++vID;

    pHBV->SetK(K);
    pHBV->AddVertices(vID);
    pHBV->EdgesMutable().reserve(2*edges.size());

    pFwdEdgeXlat.resize(edges.size(),-1);
    pRevEdgeXlat.resize(edges.size(),-1);
    //TODO: Stupidly enough, we need to sort the edges???

    //Now add the edges, and their rcs to the graph, this probably shouldn't be parallel neither (data corruption)
    for (uint64_t i=0;i<edges.size();++i){
        auto fwEdgeId = pHBV->EdgeObjectCount();
        if (edge_vertices[i].fw_v1==-1 or edge_vertices[i].fw_v2==-1) {
            std::cout<<"Edge "<<i<<" has fw_v1="<<edge_vertices[i].fw_v1<<" and fw_v2="<<edge_vertices[i].fw_v2<<"!!!"<<std::endl;
            FatalErr("Vertex not set!");
        }
        pHBV->AddEdge(edge_vertices[i].fw_v1,edge_vertices[i].fw_v2,edges[i]);
        pFwdEdgeXlat[i]=fwEdgeId;
        if ( edges[i].getCanonicalForm() == CanonicalForm::PALINDROME ) {
            pRevEdgeXlat[i] = fwEdgeId;
        }
        else {
            auto bwEdgeId = pHBV->EdgeObjectCount();
            auto rcedge = edges[i];
            rcedge.ReverseComplement();
            if (edge_vertices[i].rc_v1 == -1 or edge_vertices[i].rc_v2 == -1) {
                std::cout << "Edge " << i << " has rc_v1=" << edge_vertices[i].rc_v1 << " and rc_v2="
                          << edge_vertices[i].rc_v2 << "!!!" << std::endl;
                FatalErr("Vertex not set!");
            }
            pHBV->AddEdge(edge_vertices[i].rc_v1, edge_vertices[i].rc_v2, rcedge);
            pRevEdgeXlat[i] = bwEdgeId;
        }
    }
    for (auto i=0;i<pFwdEdgeXlat.size();++i) {
        if (pFwdEdgeXlat[i]==-1) std::cout<<"Forward translation not set for edge "<<i<<" with canonical form "<<edges[i].getCanonicalForm()<<std::endl;
    }
    for (auto i=0;i<pRevEdgeXlat.size();++i) {
        if (pRevEdgeXlat[i]==-1) std::cout<<"Reverse translation not set for edge "<<i<<" with canonical form "<<edges[i].getCanonicalForm()<<std::endl;
    }

}


void buildHKPFromHBV( HyperBasevector const& hbv,
                        std::vector<int> const& fwdEdgeXlat,
                        std::vector<int> const& revEdgeXlat,
                        HyperKmerPath* pHKP )
{
    unsigned K = hbv.K();
    pHKP->Clear();
    pHKP->SetK(K);
    pHKP->ToMutable() = hbv.To();
    pHKP->FromMutable() = hbv.From();
    pHKP->ToEdgeObjMutable() = hbv.ToEdgeObj();
    pHKP->FromEdgeObjMutable() = hbv.FromEdgeObj();
    pHKP->EdgesMutable().resize(hbv.EdgeObjectCount());
    long nextPalindrome = first_palindrome;
    long nextKmerId = 1;
    auto oItr = pHKP->EdgesMutable().begin();
    vec<bvec> const edges = hbv.Edges();//TODO: if the graph is any large, this is going to kill the machine
    for ( auto itr=edges.begin(),end=edges.end(); itr != end; ++itr,++oItr )
    {
        size_t len = itr->size();
        switch ( itr->getCanonicalForm() )
        {
        case CanonicalForm::FWD:
            oItr->Assign(nextKmerId,nextKmerId+len-K);
            nextKmerId += len;
            break;
        case CanonicalForm::PALINDROME:
            oItr->Assign(nextPalindrome,nextPalindrome+len-K);
            nextPalindrome += len;
            break;
        case CanonicalForm::REV:
            break;
        }
    }
    auto rItr = revEdgeXlat.begin();
    auto end = fwdEdgeXlat.end();
    for ( auto itr=fwdEdgeXlat.begin(); itr != end; ++itr,++rItr )
    {
        if ( pHKP->EdgeObject(*itr).empty() )
        {
            KmerPath& edge = pHKP->EdgeObjectMutable(*itr);
            edge = pHKP->EdgeObject(*rItr);
            edge.ReverseNoConcatenate();
        }
        else if ( pHKP->EdgeObject(*rItr).empty() )
        {
            KmerPath& edge = pHKP->EdgeObjectMutable(*rItr);
            edge = pHKP->EdgeObject(*itr);
            edge.ReverseNoConcatenate();
        }
        else if (*itr != *rItr) {
            int eid=itr-fwdEdgeXlat.begin();
            std::cout<<"Edge translation for old edge "<<eid<<" is "<<*itr<<" / "<<*rItr<<" and both kmer translations were already populated"<<std::endl;
            ForceAssertEq(*itr, *rItr);
        }
    }
}
