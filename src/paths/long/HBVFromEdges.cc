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
#include "math/Hash.h"

class EdgeEnd
{
public:
    EdgeEnd()=default;
    EdgeEnd( bvec const* pBV, bool rc, bool distal, unsigned K )
            : mItr(pBV,distal ? pBV->size()-K:0,rc), mK(K), mHash(FNV1a(mItr,mItr+K)) {}

    size_t getHash() const { return mHash; }
    bvec::SwitchHitterIter begin() const { return mItr; }
    bvec::SwitchHitterIter end() const { return mItr+mK; }

    struct Hasher
    { typedef EdgeEnd argument_type;
        size_t operator()( EdgeEnd const& ee ) const { return ee.getHash(); } };

    friend bool operator<( EdgeEnd const& ee1, EdgeEnd const& ee2 )
    { if ( ee1.mHash != ee2.mHash ) return ee1.mHash < ee2.mHash;
        for ( auto i1=ee1.begin(),e1=ee1.end(),i2=ee2.begin(); i1!=e1; ++i1,++i2 )
            if ( *i1 != *i2 ) return *i1 < *i2;
        return false; }

    friend bool operator==( EdgeEnd const& ee1, EdgeEnd const& ee2 )
    { if ( ee1.mHash != ee2.mHash ) return false;
        for ( auto i1=ee1.begin(),e1=ee1.end(),i2=ee2.begin(); i1!=e1; ++i1,++i2 )
            if ( *i1 != *i2 ) return false;
        return true; }

private:
    bvec::SwitchHitterIter mItr;
    unsigned mK;
    size_t mHash;
};

class EdgeEndWrapper{
public:
    EdgeEnd ee;
    bool distal,rc;
    uint64_t edge_index;
    EdgeEndWrapper(const uint64_t _edge_index, const bvec & bvseq,const bool &_rc,const bool &_distal, const unsigned overlap) {
        edge_index=_edge_index;
        distal=_distal;
        rc=_rc;
        ee=EdgeEnd(&bvseq,_rc,_distal,overlap);
    };
    bool operator<(const EdgeEndWrapper & other){return ee<other.ee;};//WARNING, for simplicity
    bool operator==(const EdgeEndWrapper & other){return ee==other.ee;};//WARNING, for simplicity

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

    std::vector<EdgeEndWrapper> ends;
    ends.reserve(2*edges.size());

    for (uint64_t i=0;i<edges.size();++i){
        ends.push_back(EdgeEndWrapper(i,edges[i],false,false,K-1));
        ends.push_back(EdgeEndWrapper(i,edges[i],false,true,K-1));
        if ( edges[i].getCanonicalForm() != CanonicalForm::PALINDROME ) {
            ends.push_back(EdgeEndWrapper(i,edges[i],true,false,K-1));
            ends.push_back(EdgeEndWrapper(i,edges[i],true,true,K-1));
        }
    }
    std::sort(ends.begin(),ends.end());


    //number the vertices and fill a list of vertex numbers for the edges
    std::vector<edge_vertex_list> edge_vertices;
    edge_vertices.resize(edges.size(),{-1,-1,-1,-1});
    uint64_t vID=0;
    //uint64_t ecount[5]={0,0,0,0,0};
    //uint64_t vcount=0;
    //Can't be parallel, because vID depends on previous iteration
    for (auto i=0; i<ends.size(); ++i){
        if (i>0 and not (ends[i-1]==ends[i])) {
            vID++;
            //if (vcount>5) vcount=5;
            //++ecount[vcount-1];
            //vcount=0;
        }
        //++vcount;
        if (!ends[i].rc) {
            if (!ends[i].distal) edge_vertices[ends[i].edge_index].fw_v1=vID;
            else edge_vertices[ends[i].edge_index].fw_v2=vID;
        } else {
            if (!ends[i].distal) edge_vertices[ends[i].edge_index].rc_v1=vID;
            else edge_vertices[ends[i].edge_index].rc_v2=vID;
        }
    }
    ++vID;

    //std::cout<<Date()<<": vertices by edge count:  1:"<<ecount[0]<<" 2:"<<ecount[1]
    //         <<" 3:"<<ecount[2]<<" 4:"<<ecount[3]<<" 5+:"<<ecount[4]<<std::endl;
    pHBV->SetK(K);
    pHBV->AddVertices(vID);
    pHBV->EdgesMutable().reserve(2*edges.size());

    pFwdEdgeXlat.resize(edges.size(),-1);
    pRevEdgeXlat.resize(edges.size(),-1);

    //Now add the edges, and their rcs to the graph, this probably shouldn't be parallel neither (data corruption)
    for (uint64_t i=0;i<edges.size();++i){
        auto fwEdgeId = pHBV->EdgeObjectCount();
        bvec edge = edges[i];
        pHBV->AddEdge(edge_vertices[i].fw_v1,edge_vertices[i].fw_v2,edge);
        pFwdEdgeXlat[i]=fwEdgeId;
        if ( edges[i].getCanonicalForm() == CanonicalForm::PALINDROME ) {
            pRevEdgeXlat[i] = fwEdgeId;
        }
        else {
            auto bwEdgeId = pHBV->EdgeObjectCount();
            pHBV->AddEdge(edge_vertices[i].rc_v1, edge_vertices[i].rc_v2, edge);
            pHBV->EdgeObjectMutable(bwEdgeId).ReverseComplement();
            pRevEdgeXlat[i] = bwEdgeId;
        }
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
