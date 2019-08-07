#ifndef GAP_TOY_TOOLS_H
#define GAP_TOY_TOOLS_H

#include "Bitvector.h"
#include "CoreTools.h"
#include "Intvector.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "system/SpinLockedData.h"
#include <memory>
#include <fstream>


void BasesToGraph( vecbasevector& bpathsx, const int K, HyperBasevector& hb );


void MakeLocalAssembly2( VecEFasta& corrected, const vec<int>& lefts, const vec<int>& rights,
     SupportedHyperBasevector& shb, const int K2_FLOOR,
     vecbasevector& creads, /*LongProtoTmpDirManager& tmp_mgr,*/ vec<int>& cid,
     vec<pairing_info>& cpartner );


void LayoutReads( const HyperBasevector& hb, const vec<int>& inv, 
     const vecbasevector& bases, const ReadPathVec& paths, 
     std::vector<std::vector<int>>& layout_pos, std::vector<std::vector<int64_t>>& layout_id,
     std::vector<std::vector<bool>>& layout_or );


void RemoveHangs( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_del );

void Degloop( const int mode, HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const vecbasevector& bases, const VecPQVec& quals, const double min_dist,
     const int verbosity = 0 );

// DegloopCore: H = HyperBasevector or HyperBasevectorX.

template<class H> void DegloopCore( const int mode, H& hb, vec<int>& inv, 
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals,
     const VecULongVec& paths_index, const int v, const int pass,
     const double min_dist, vec<int>& EDELS, const int verbosity,
     const vec<int>* ids = NULL );

void Patch( HyperBasevector& hb,
     vec<HyperBasevector>& mhbp,
     vecbvec& new_stuff );

void CleanupCore( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void Cleanup( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void CleanupLoops( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths );

void RemoveUnneededVertices( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

// For two-edge loops only:

void RemoveUnneededVerticesLoopsOnly( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

void RemoveUnneededVertices2( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths, Bool debug = false );

void RemoveSmallComponents3( HyperBasevector& hb, 
     const Bool remove_small_cycles = False );


void AddNewStuff( vecbvec& new_stuff, HyperBasevector& hb, vec<int>& inv2, 
     ReadPathVec& paths2, const vecbasevector& bases, const VecPQVec& quals, 
     const int MIN_GAIN, const int EXT_MODE );

void ExtendPath( ReadPath& p, const int64_t i, const HyperBasevector& hb, 
     const vec<int>& to_right, const bvec& bases,
     const QualVec& quals, const int min_gain, const Bool verbose,
     const int mode );

void ExtendPath2( ReadPath& p, const int64_t i, const HyperBasevector& hb, 
     const vec<int>& to_left, const vec<int>& to_right, const bvec& bases,
     const QualVec& quals, const int min_gain, const Bool verbose,
     const int mode );

// a class for dealing with bubbles and support
// the basic usage is
// 1) bubble_logger(graph,inv)
// 2) for each read, call log_read(read,qual,read-path)
// 3) call getData() to analyze data
class bubble_logger{
public:
    //stores the state of the bubble
    struct bubble_data_t{
        struct support_t{
            support_t(int a, int b):read_branch(a),branch_branch(b){}
            int read_branch;   // qsum between read and path
            int branch_branch; // qsum difference against the losing branch
        };
        bubble_data_t(const int branch0_edge,const int branch1_edge)
            :branch_edges({branch0_edge,branch1_edge})
            ,branch_supports(2)
            ,lock_ptr(new SpinLockedData)
            {}
        bubble_data_t(const int branch0_edge,const int branch1_edge
                   ,const int branch0_edge_rc, const int branch1_edge_rc)
            :branch_edges({branch0_edge,branch1_edge,branch0_edge_rc,branch1_edge_rc})
            ,branch_supports(4)
            ,lock_ptr(new SpinLockedData)
            {}
        bubble_data_t(bubble_data_t && other)noexcept
            :branch_edges(std::move(other.branch_edges))
            ,branch_supports(std::move(other.branch_supports))
            ,lock_ptr(std::move(other.lock_ptr)) { }

        void addSupport(size_t branch, support_t const&weight);
        vec<support_t> getSupport(size_t branch)const{
            SpinLocker locker(*lock_ptr);
            return branch_supports[branch];
        };
        vec<int>const& getEdges()const{ return branch_edges;};

    private:
        bubble_data_t();
        bubble_data_t(bubble_data_t const&other);
        const vec< int > branch_edges;
        vec< vec<support_t> > branch_supports;
        std::unique_ptr<SpinLockedData> lock_ptr;
    };

    // given the problem statement this is the only meaningful constructor
    bubble_logger(const HyperBasevector& hb, const vec<int>& inv);

    // if one or more edges in rp is part of a bubble, perform gap-free alignment on the path as the alternate path
    // and collect the result
    // returns true if at any point (qsum of alt path) < (qsum of orig path)
    bool log_read(basevector const&read, QualVec const&qual, ReadPath const&rp, bool bVerbose=false);

    //do a gap-free alignment of the read against the graph according to rp
    //returns the sum of read-quality-score at the position the sequence mismatches
    int getQ(basevector const&read, QualVec const&qual, ReadPath const&rp, const QualVec::value_type min_q=4);

    // assign weight to the bubble-branch of an edge
    void addWeight(int edge,bubble_data_t::support_t const&weight){
        const auto& b_b=edge_bubble_branch_[edge];
        int bubble_data_size = bubble_data_.size();
        ForceAssertLt(b_b.first, bubble_data_size);
        bubble_data_[ b_b.first ].addSupport(b_b.second,weight);
    }

    // if edge is part of a bubble, return the edge index corresponding to another branch, or -1 otherwise
    int alt(int edge)const{ return edge_alt_[edge];}

    // if an edge is in a bubble
    bool inBubble(int edge)const{ return edge_bubble_branch_[edge].first>=0;}

    // return the current list of bubble data
    const std::vector<bubble_data_t>& getData()const{return bubble_data_;};

private:
    bubble_logger();
    const HyperBasevector& hb_;
    vec<int> edge_alt_;                            // [edge_id] -> edge id of another branch
    vec< std::pair<int,int> > edge_bubble_branch_; // [edge_id] -> (bubble,branch)
    std::vector< bubble_data_t > bubble_data_;             // [bubble_id] -> data

};
std::ostream& operator<<(std::ostream&os, bubble_logger const& in);
std::ostream& operator<<(std::ostream& os, bubble_logger::bubble_data_t const&in);
void PopBubbles( HyperBasevector& hb , const vec<int>& inv2
               , const vecbasevector& bases, const VecPQVec& quals, const ReadPathVec& paths2);

// HIGHLY INCOMPLETE:

void Validate( const HyperBasevector& hb, const vec<int>& inv, 
     const ReadPathVec& paths );


void TestInvolution( const HyperBasevector& hb, const vec<int>& inv );

void DeleteFunkyPathPairs( const HyperBasevector& hb, const vec<int>& inv,
     const vecbasevector& bases, ReadPathVec& paths, const Bool verbose );



void ReroutePaths( const HyperBasevector& hb, const vec<int>& inv,
     ReadPathVec& paths, const vecbasevector& bases, const VecPQVec& quals );

void Tamp( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     const int max_shift );


void LogTime( const double clock, const String& what, const String& work_dir = "" );


void PartnersToEnds( const HyperBasevector& hb, ReadPathVec& paths,
                        const vecbasevector& bases, const VecPQVec& quals );



void FragDist( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const String out_file );

void UnwindThreeEdgePlasmids(HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths);

// Compute fraction of edges in the assembly whose copy number is within 10% of
// the nearest integer value.

double CNIntegerFraction(const HyperBasevector& hb, const vec<vec<covcount>>& covs,
			 const double frac = 0.1, const int min_edge_size = 2000);

void RemoveUnneededVerticesGeneralizedLoops( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths );

void BuildAll( vecbasevector& all, const HyperBasevector& hb, 
     const int64_t extra = 0 );

void TranslatePaths( ReadPathVec& paths2, const HyperBasevector& hb3,
     const vec<vec<int>>& to3, const vec<int>& left3 );


#endif
