///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"

#include "paths/long/EvalByReads.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/RefTrace.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "random/Bernoulli.h"
#include "reporting/PerfStat.h"
#include "util/NullOStream.h"
#include "paths/long/ReadOriginTracker.h"

namespace { // open anonymous namespace

void FixPath( vec<int>& p, const vec< triple<int,int,int> >& merges,
     const vec< vec<int> >& merges_index, const HyperBasevector& hb_orig, 
     int& left_add, int& right_add )
{    left_add = 0, right_add = 0;
     while(1)
     {    Bool changed = False;
          vec<int> mids;
          for ( int l = 0; l < p.isize( ); l++ )
               if ( p[l] >= 0 ) mids.append( merges_index[ p[l] ] );
          UniqueSort(mids);
          for ( int mj = 0; mj < mids.isize( ); mj++ )
          {    int j = mids[mj];
               int e1 = merges[j].first, e2 = merges[j].second;
               ForceAssert( e1 != e2 );
               int enew = merges[j].third;
               for ( int l = 0; l < p.isize( ); l++ )
               {    if ( l < p.isize( ) - 1 && p[l] == e1 && p[l+1] == e2 )
                    {    p[l] = enew;
                         for ( int m = l+2; m < p.isize( ); m++ )
                              p[m-1] = p[m];
                         p.pop_back( );
                         changed = True;    }
                    else if ( l == p.isize( ) - 1 && p[l] == e1 )
                    {    right_add += hb_orig.EdgeLengthKmers(e2);
                         p[l] = enew;
                         changed = True;    }
                    else if ( l == 0 && p[l] == e2 )
                    {    left_add += hb_orig.EdgeLengthKmers(e1);
                         p[l] = enew;
                         changed = True;    }
                    if (changed) break;    }
               if (changed) break;    }
          if ( !changed ) break;    }    }

} // close anonymous namespace

void SupportedHyperBasevector::RemoveUnneededVertices0( 
     vec< triple<int,int,int> >& merges )
{    vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    int e1 = EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = EdgeObjectIndexByIndexFrom( i, 0 );
               basevector p = Cat( e1, e2 );
               int enew = EdgeObjectCount( );
               int re1 = Inv(e1), re2 = Inv(e2);
               ForceAssert( ( re1 < 0 && re2 < 0 ) || ( re1 >= 0 && re2 >= 0 ) );
               int v = To(i)[0], w = From(i)[0];
               Bool loop = ( v == w && From(v).solo( ) && To(v).solo( ) );
               // v --e1--> i --e2--> w
               merges.push( e1, e2, enew );
               JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);
               if ( re1 < 0 ) InvMutable( ).push_back(-1);
               else if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *
                    InvMutable( ).push_back(enew);    }
               else if ( re2 == e2 && re1 != e1 )
               {    // * --e1--> * --e2=re2--> * --re1--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( enew, re1, enew2 );
                    basevector p2 = TrimCat( K( ), p, EdgeObject(re1) );
                    to_left.push_back(v), to_right.push_back( to_right[re1] );
                    JoinEdges( w, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 != e2 )
               {    // * --re2--> * --e1=re1--> * --e2--> *
                    int enew2 = EdgeObjectCount( );
                    merges.push( re2, enew, enew2 );
                    basevector p2 = TrimCat( K( ), EdgeObject(re2), p );
                    to_left.push_back( to_left[re2] ), to_right.push_back(w);
                    JoinEdges( v, p2 );
                    InvMutable( ).push_back( -1, enew2 );    }
               else if ( re1 == e1 && re2 == e2 )
               {    if (loop) InvMutable( ).push_back(-1);
                    else
                    {    // not sure if this can happen
                         ForceAssert( 0 == 1 );    }    }
               else
               {    // e1, e2, re1, re2 all different
                    int renew = EdgeObjectCount( );
                    basevector rp = Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    InvMutable( ).push_back(renew, enew);    }    }    }    }

void SupportedHyperBasevector::RemoveUnneededVertices( )
{    vec< triple<int,int,int> > merges;
     RemoveUnneededVertices0(merges);
     vec< vec<int> > merges_index( EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     #pragma omp parallel for
     for ( int i = 0; i < NPaths( ); i++ )
     {    int left_add, right_add;
          FixPath( PathMutable(i), merges, merges_index, *this,
               left_add, right_add );    }
     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    int left_add1, right_add1, left_add2, right_add2;
          FixPath( PairLeftMutable(i), merges, merges_index, *this,
               left_add1, right_add1 );
          FixPath( PairRightMutable(i), merges, merges_index, *this,
               left_add2, right_add2 );
          AddTrim( i, right_add1 + left_add2 );    }
     UniqueOrderPaths( );
     RemoveEdgelessVertices( );    }

void SupportedHyperBasevector::DeleteUnusedPaths( )
{    vec<Bool> used, to_delete( NPaths( ), False );
     Used(used);
     for ( int i = 0; i < NPaths( ); i++ )
     {    for ( int j = 0; j < Path(i).isize( ); j++ )
               if ( Path(i,j) >= 0 && !used[ Path(i,j) ] ) to_delete[i] = True;    }
     EraseIf( PathsMutable( ), to_delete );
     EraseIf( WeightsFwMutable( ), to_delete );    
     EraseIf( WeightsRcMutable( ), to_delete );    
     // The following if is temporary - until origins fully implemented.
     if ( to_delete.size( ) == WeightsFwOrigin( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );    
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
                    if ( p[j] >= 0 && !used[ p[j] ] ) to_delete[i] = True;    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }

void SupportedHyperBasevector::RemoveDeadEdgeObjects0( )
{    vec<Bool> used;
     Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( !InvDef(i) ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ Inv(i) ] );    }
     InvMutable( ) = inv2;
     HyperBasevector::RemoveDeadEdgeObjects( );    }

void SupportedHyperBasevector::RemoveDeadEdgeObjects( )
{    vec<Bool> used;
     Used(used);
     vec<int> to_new_id( used.size( ), -1 );
     {    int count = 0;
          for ( int i = 0; i < used.isize( ); i++ )
               if ( used[i] ) to_new_id[i] = count++;    }
     vec<int> inv2;
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    if ( !used[i] ) continue;
          if ( !InvDef(i) ) inv2.push_back(-1);
          else inv2.push_back( to_new_id[ Inv(i) ] );    }
     InvMutable( ) = inv2;

     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    vec<int>& p = PathMutable(i);
          for ( int j = 0; j < Path(i).isize( ); j++ )
          {    if ( Path(i,j) >= 0 )
               {    int n = to_new_id[ Path(i,j) ];
                    if ( n < 0 ) to_delete[i] = True;
                    else PathMutable(i)[j] = n;    }    }    }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( p[j] >= 0 )
                    {    int n = to_new_id[ p[j] ];
                         if ( n < 0 ) to_delete[i] = True;
                         else p[j] = n;    }    }    }    }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );

     HyperBasevector::RemoveDeadEdgeObjects( );    }

void SupportedHyperBasevector::UniqueOrderPaths( )
{    
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ),
               WeightsFwOriginMutable( ), WeightsRcOriginMutable( ) );    }
     else SortSync( PathsMutable( ), WeightsFwMutable( ), WeightsRcMutable( ) );
     vec<Bool> to_delete( NPaths( ), False );
     for ( int i = 0; i < NPaths( ); i++ )
     {    int j = Paths( ).NextDiff(i);
          fix64_6 cfw = 0.0, crc = 0.0;
          for ( int k = i; k < j; k++ )
               cfw += WeightFw(k);
          WeightFwMutable(i) = cfw;
          for ( int k = i; k < j; k++ )
               crc += WeightRc(k);
          WeightRcMutable(i) = crc;
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    for ( int k = i+1; k < j; k++ )
               {    WeightFwOriginMutable(i).append(
                         WeightFwOriginMutable(k) );
                    WeightRcOriginMutable(i).append(
                         WeightRcOriginMutable(k) );    }
               Sort( WeightFwOriginMutable(i) );
               Sort( WeightRcOriginMutable(i) );    }
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( cfw + crc == 0 ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PathsMutable( ), to_delete );
     // Temporary if.
     if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
     {    EraseIf( WeightsFwOriginMutable( ), to_delete );
          EraseIf( WeightsRcOriginMutable( ), to_delete );    }
     EraseIf( WeightsFwMutable( ), to_delete );
     EraseIf( WeightsRcMutable( ), to_delete );

     // Now do pairs.

     SortSync( PairsMutable( ), PairDataMutable( ) );
     to_delete.resize_and_set( NPairs( ), False );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          int j = Pairs( ).NextDiff(i);
          vec<pair_point> x;
          for ( int k = i; k < j; k++ )
               x.append( PairData(k) );
          Sort(x);
          PairDataMutable(i) = x;
          for ( int k = i+1; k < j; k++ )
               to_delete[k] = True;
          if ( PairData(i).empty( ) ) to_delete[i] = True;
          i = j - 1;   }
     EraseIf( PairsMutable( ), to_delete );
     EraseIf( PairDataMutable( ), to_delete );    }


namespace
{

template<int K> void OrientCore( SupportedHyperBasevector& shb, 
     const vecbasevector& genome )
{    vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup1( genome, kmers_plus );
     int fw = 0, rc = 0;
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          kmers[i] = kmers_plus[i].first;
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     vec< vec<int> > comp;
     shb.ComponentsE(comp);
     vec<Bool> reversed( shb.EdgeObjectCount( ), False );
     for ( int c = 0; c < comp.isize( ); c++ )
     {    if ( comp[c].empty( ) ) continue;
          if ( shb.Inv( comp[c][0] ) >= 0 ) continue;
          int fw = 0, rc = 0;
          for ( int i = 0; i < comp[c].isize( ); i++ )
          {    int e = comp[c][i];
               const basevector& E = shb.EdgeObject(e);
               vec<Bool> fwv( E.isize( ) - K + 1, False ); 
               vec<Bool> rcv( E.isize( ) - K + 1, False );
               #pragma omp parallel for
               for ( int j = 0; j <= E.isize( ) - K; j++ )
               {    kmer<K> x;
                    x.SetToSubOf( E, j );
                    if ( BinMember( kmers, x ) ) fwv[j] = True;
                    x.ReverseComplement( );
                    if ( BinMember( kmers, x ) ) rcv[j] = True;    }
               fw += Sum(fwv);
               rc += Sum(rcv);    }
          if ( rc > fw )
          {    int v = to_left[ comp[c][0] ];
               for ( int i = 0; i < comp[c].isize( ); i++ )
               {    int e = comp[c][i];
                    shb.EdgeObjectMutable(e).ReverseComplement( );
                    reversed[e] = True;    }
               shb.ReverseComponent(v);    }    }
     for ( int i = 0; i < shb.NPaths( ); i++ )
          if ( reversed[ shb.Path(i)[0] ] ) shb.PathMutable(i).ReverseMe( );

     // Reverse pairs.  Note that if a pair goes from one component to another,
     // and only one of the components is reversed, then the pair no longer makes
     // sense.  Not sure what to do about that.

     for ( int i = 0; i < shb.NPairs( ); i++ )
     {    Bool rev1 = reversed[ shb.PairLeft(i)[0] ];
          Bool rev2 = reversed[ shb.PairRight(i)[0] ];
          if (rev1) shb.PairLeftMutable(i).ReverseMe( );
          if (rev2) shb.PairRightMutable(i).ReverseMe( );
          if (rev1)
          {    shb.PairMutable(i) = std::make_pair(
                    shb.PairRight(i), shb.PairLeft(i) );    }    }    }

template <int K>
struct OrientCoreFunctor
{
    void operator()( SupportedHyperBasevector& shb,
                        const vecbasevector& genome )
    { OrientCore<K>(shb,genome); }
};

}

void OrientToReference( SupportedHyperBasevector& shb, const vecbasevector& genome,
     const long_logging& logc )
{    double clock = WallClockTime( );
     BigK::dispatch<OrientCoreFunctor>(shb.K(),shb,genome);
     shb.UniqueOrderPaths( );
     shb.TestValid(logc);
     REPORT_TIME( clock, "used orienting to reference" );    }

void AssessAssembly( const String& SAMPLE, const SupportedHyperBasevector& shb, 
     const HyperEfasta& he, const vec<Bool>& hide, const String& TMP, 
     const ref_data& ref, const String& HUMAN_CONTROLS,
     const long_logging& logc, const uint NUM_THREADS, RefTraceControl RTCtrl )
{
     const vecbasevector& G = ref.G;
     const vec<HyperBasevector>& GH = ref.GH; 
     const vec<bool>& is_circular = ref.is_circular;

     if (logc.STATUS_LOGGING) DATE_MSG( "assessing assembly" );
     if ( logc.COUNT_COV > 0 ) CountCov( shb, TMP, logc.COUNT_COV );
     int assembly_count = 0, reference_count = 0;
     if ( logc.READ_EVAL == "True" )
     {    double eclock = WallClockTime( );
          Bool print_a = False;
          vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
          vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
          vec<basevector> R;
          for ( int g = 0; g < (int) G.size( ); g++ )
               R.push_back( G[g] );
          HyperBasevector hb_R( shb.K( ), R );
          EvalByReads( shb, hb_R, bases, quals, assembly_count, reference_count, 
               print_a, logc.PRINT_FAVORING_REF );
          REPORT_TIME( eclock, "used evaluating by reads" );    }
     std::cout << "\n=================================================================="
          << "==================\n\n";
     std::cout << "SUMMARY STATS  --  dexter longread stats shown as -[...]-\n\n";
     if ( SAMPLE != "unknown" && logc.REFTRACE == "True" )
     {    RTCtrl.ReadySampleLookup( );
          ReadOriginTracker read_tracker(RTCtrl);
          RefTraceHeuristics rth = RefTraceHeuristics( );
          Bool fix_bug = False;
          NullOStream nullOS;
          if ( HUMAN_CONTROLS != "" )
          {    fix_bug = True;
               int np = 3;
               vec<int> penalty(np), gaps(np), meta_events(np);
               std::vector<std::ostringstream> out(np);
               for ( int pass = 0; pass < np; pass++ )
               {    if ( pass == 0 )
                    {    rth.min_group_frac = 0.1;
                         rth.min_group_save = 200;    }
                    if ( pass == 1 )
                    {    rth.max_offset_diff = 30;
                         rth.max_error_rate = 0.31;
                         rth.offset_add = 5;
                         rth.min_group_frac = 0.1;
                         rth.max_twiddle = 5;    }
                    if ( pass == 2 )
                    {    rth.max_offset_diff = 350;
                         rth.max_error_rate = 0.31;
                         rth.offset_add = 5;
                         rth.min_group_frac = 0.75;
                         rth.max_twiddle = 120;    }
                    RefTraceResults res = RefTrace( ref, shb, shb.Inv( ),
                            logc.verb[ "REFTRACE" ], logc, out[pass], rth, "",
                            fix_bug);
                    //RefTraceResults res = RefTrace( ref, shb, shb.Inv( ),
                    //        logc.verb[ "REFTRACE" ], logc, out[pass], nullOS,
                    //        rth, "", fix_bug, RTCtrl, &read_tracker);
                    penalty[pass] = res.penalty;
                    gaps[pass] = res.gaps;
                    meta_events[pass] = res.meta_events;    }
               int pix = -1;
               for ( int pi = 0; pi < np; pi++ )
               {    if ( penalty[pi] == Min(penalty) )
                    {    pix = pi;
                         std::cout << out[pi].str( );
                         if ( logc.SHOW_REFTRACE_EVENTS && logc.PERF_STATS )
                         {    PerfStat::log( ) << std::fixed << std::setprecision(0)
                                   << PerfStat( "error_meta_events", 
                                   "error meta-events", meta_events[pi] );
                              PerfStat::log( ) << std::fixed << std::setprecision(0)
                                   << PerfStat( "gaps", "gaps", gaps[pi] );    }
                         break;    }    }    
               for ( int pi = 0; pi < np; pi++ )
               {    if ( pi == pix ) continue;
                    std::istringstream in( out[pi].str( ) );
                    String line;
                    while(1)
                    {    getline( in, line );
                         if ( in.fail( ) ) break;
                         if ( line.Contains( " seconds " ) 
                              || line.Contains( " minutes " )
                              || line.Contains( " hours " ) 
                              || line.Contains( " days " ) )
                         {    std::cout << line << "\n";    }    }    }    }
          else {
            RefTraceAndCallVaraint( ref, shb, shb.Inv( ), logc.verb[ "REFTRACE"],
                    logc, std::cout, nullOS, rth, "", fix_bug, RTCtrl, &read_tracker); }
     }
     int hidden = 0;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          if ( shb.Inv(e) >= 0 && shb.Inv(e) < e ) hidden++;
     std::cout << shb.EdgeObjectCount( ) << " edges (" << shb.EdgeObjectCount( ) - hidden
          << " visible, " << hidden << " hidden)\n";
     std::cout << he.EdgeObjectCount( ) << " edges in efasta" 
          << " (" << he.EdgeObjectCount( ) - Sum(hide) << " visible, "
          << Sum(hide) << " hidden)" << std::endl;
     if ( !shb.Acyclic( ) ) std::cout << "assembly has a cycle\n";
     std::cout << "assembly has " << shb.NComponents( ) << " components\n";
     if ( G.size( ) > 0 )
     {    int64_t gsize = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               gsize += G[g].size( );
          int excess_edges = shb.EdgeObjectCount( ) - (int) G.size( );
          double excess_edges_per = 1000000.0 * double(excess_edges) / double(gsize);
          std::cout << "-[" << std::setiosflags(std::ios::fixed) << std::setprecision(0) 
               << excess_edges_per 
               << std::resetiosflags(std::ios::fixed) << "]- excess edges per Mb" 
               << " (total excess edges = " << excess_edges << ")" << std::endl;
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << std::setprecision(0) 
                    << PerfStat( "excess_edges", "excess edges per Mb",
                    excess_edges_per );    }    }
     if ( SAMPLE != "unknown" ) ReportExcessEdges( he, G, logc.PERF_STATS );
     if ( logc.READ_EVAL == "True" )
     {    double eclock = WallClockTime( );
          int64_t gsize = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               gsize += G[g].size( );
          double bad_reads_per = 1000000.0 * double(reference_count) / double(gsize);
          std::cout << "-[" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << bad_reads_per
               << std::resetiosflags(std::ios::fixed) << "]- reads favoring reference per Mb" 
               << " (favors ref = " << reference_count << ")\n";
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << std::setprecision(1) << PerfStat( 
                    "bad_reads_per", "bad reads per Mb", bad_reads_per );    }
          vec<int> sources, sinks;
          shb.Sources(sources), shb.Sinks(sinks);
          int expected = 0;
          for ( int g = 0; g < (int) G.size( ); g++ )
               if ( !is_circular[g] ) expected += 2;
          int excess_ends = sources.isize( ) + sinks.isize( ) - expected;
          double excess_ends_per = double(excess_ends)/(double(gsize)/1000000.0);
          std::cout << "-[" << std::setiosflags(std::ios::fixed) << std::setprecision(2) 
               << excess_ends_per << std::resetiosflags(std::ios::fixed) 
               << "]- excess ends per Mb (excess ends = " << excess_ends << ")\n";
          if ( logc.PERF_STATS )
          {    PerfStat::log( ) << std::fixed << std::setprecision(2) 
                    << PerfStat( "excess_ends_per", "excess ends per Mb",
                    excess_ends_per );    }
          REPORT_TIME( eclock, "used evaluating by reads tail" );    }    }

void CountCov( const SupportedHyperBasevector& shb, const String& TMP,
     const int gp1 )
{
     // Determine file count.

     int fcount = 0;
     for ( fcount = 0; ; fcount++ )
          if ( !IsRegularFile( TMP + "/" + ToString(fcount) + ".fastb" ) ) break;

     // Set up counts.

     vec< vec<int> > cov( shb.EdgeObjectCount( ), vec<int>(2,0) );

     // Build data structures.

     const int L = 12;
     HyperBasevector hb_fw(shb), hb_rc(shb);
     hb_rc.Reverse( );
     vec<int> to_right_fw, to_right_rc;
     hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);
     vecbasevector x_fw, x_rc;
     for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ )
          x_fw.push_back( hb_fw.EdgeObject(i) );
     for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ )
          x_rc.push_back( hb_rc.EdgeObject(i) );
     VecIntPairVec locs_fw, locs_rc;
     CreateGlocs( x_fw, L, locs_fw );
     CreateGlocs( x_rc, L, locs_rc );

     // Heuristics.

     const double prox = 20 * 1000;
     const int max_qual = 100 * 1000;
     const int min_pos_evidence = 8;

     // Go through the files.

     for ( int f = 0; f < fcount; f++ )
     {    vecbasevector bases( TMP + "/" +ToString(f) + ".fastb" );
          vecqualvector quals( TMP + "/" +ToString(f) + ".qualb" );
          int gid = ( f < gp1 ? 0 : 1 );

          // Align reads.

          #pragma omp parallel for
          for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
          {    int n = KmerId( bases[id], L, 0 );
               vec<read_place> places;
               const int infinity = 1000000000;
               int qual_sum = infinity;
               const int min_qual = 1;
               FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw,
                    to_right_rc, locs_fw, locs_rc, places, qual_sum, min_qual,
                    prox );

               vec<Bool> to_delete( places.size( ), False );
               for ( int i = 0; i < places.isize( ); i++ )
                    if ( places[i].Qsum( ) > max_qual ) to_delete[i] = True;
               EraseIf( places, to_delete );
               Bool OK = False;
               if ( places.solo( ) ) OK = True;
               if ( places.size( ) == 2 && shb.InvDef( places[0].E( ).front( ) ) )
               {    vec<int> e = places[0].E( );
                    for ( int j = 0; j < e.isize( ); j++ )
                         e[j] = shb.Inv( e[j] );
                    // Doing this two ways, not sure which is right.
                    if ( e == places[1].E( ) ) OK = True;
                    e.ReverseMe( );
                    if ( e == places[1].E( ) ) OK = True;    }

               if (OK)
               {    for ( int i = 0; i < places.isize( ); i++ )
                    {    const read_place& p = places[i];
                         #pragma omp critical
                         {    for ( int j = 0; j < p.N( ); j++ )
                                   cov[ p.E(j) ][gid]++;    }    }    }    }    }

     // Announce results.

     std::cout << "\ncoverage:\n";
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    if ( cov[e][0] < min_pos_evidence ) continue;
          if ( cov[e][1] > 0 ) continue;
          std::cout << e << " --> " << printSeq( cov[e] ) << "\n";    }    }


