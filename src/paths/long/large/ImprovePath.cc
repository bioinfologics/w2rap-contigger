///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Could be generalized to allow seeding on larger kmers, which would find more
// alignments.  Separate problem: find a good set of seeds.
//
// Also could treat the existing alignment as a seed.
//
// Run time for lookup table generation: 2 minutes.
//
// Run time for main loop: 1.02 seconds / million reads.
// (Extrapolates to ~13 minutes.)

/*

review after "Compare core to existing alignment".  How are backups handled?

memory for 40-mer table
10 + 8 --> 24 bytes * 10G = 240 GB

could we just make this run, and leave in place logging machinery so we
could continue the investigation later?

====================================================================================

*/

// Look for read path placement.  Take the first 20-mer in the read.  Suppose
// that it has at most 10 placements on the assembly.  Suppose that each extension
// of one of these placements goes far enough to align the read, and suppose that
// there are at most 100 extensions (including partials).  Score each extension.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "math/Hash.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/MakeGaps.h"

namespace { // open anonymous namespace

//[GONZA] TODO: why this is istringstream if the passes an ostringstream to this function every time!? changed here
//void FinalPrint( std::istringstream* pout )
void FinalPrint( std::ostringstream* pout )
{
     if ( pout != NULL )
     {
          #pragma omp critical
          {    std::cout << pout->str( );    }
          delete pout;    }    }

template<int K> void MakeKmerLookup0H( const vecbasevector& unibases,
     vec< triple<uint64_t,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               kmers_plus[r].first = FNV1a( x.Ints( ), x.Ints( ) + (K+15)/16 );
               kmers_plus[r].second = i;
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);    }

} // close anonymous namespace

template<int L> void ImprovePath( const vec<int>& rstarts,
     const vec< std::pair<int64_t, std::pair<int,int> > >& locsx,
     const ReadPathVec& paths, const int64_t xi, ReadPath& p, const int id,
     const HyperBasevector& hb, const vec<int>& inv, const vec<int>& to_left,
     const vec<int>& to_right, const basevector& b, const qualvector& q,
     const vec< triple<kmer<L>,int,int> >& kmers_plus, const path_improver& pimp,
     path_improver::path_status& status )
{
     // Logging.

     //[GONZA] TODO: check this, this is required to be a istringstream instead of ostringstream by some functions FinalPrint(pout)
     //std::ostringstream* pout = NULL;
     std::ostringstream* pout = NULL;
     //[GONZA] TODO: commented this line
     //if ( pimp.Logging( ) ) pout = new std::istringstream;

     // Heuristics.  See also ImprovePaths.

     // const int L = 20;                     // seed kmer size
     const int max_locs1 = 10;                // max seed locs
     const int max_locs2 = 100;               // max extensions
     const int window = 60;                   // window for junk assessment
     const int max_mis = 6;                   // max mismatches in window
     const int min_gain = 5;                  // min required q score improvement
     const int flank = 10;                    // for raising q2
     const int max_ok = 150;                  // max quality score sum
     const int min_back_over = 60;            // to back up

     // Look for prelocated seeds.

     vec< std::pair<int,int> > locs;
     int64_t low = LowerBound1( locsx, xi ), high = UpperBound1( locsx, xi );
     for ( int64_t j = low; j < high; j++ )
          locs.push_back( locsx[j].second );

     // Find seeds.  We try the first four nonoverlapping 20-mers on the read
     // (spanning 80 bases in total).

     kmer<L> x;
     for ( int ri = 0; ri < rstarts.isize( ); ri++ )
     {    int rstart = rstarts[ri];
          if ( rstart + L > b.isize( ) ) continue;
          x.SetToSubOf( b, rstart );
          int64_t low = LowerBound1( kmers_plus, x );
          int64_t high = UpperBound1( kmers_plus, x );
          if ( high - low <= max_locs1 )
          {    for ( int64_t li = low; li < high; li++ )
               {    int e = kmers_plus[li].second;
                    int start = kmers_plus[li].third - rstart;
                    if ( start >= 0 ) locs.push( e, start );
                    else
                    {    int w = to_left[e];
                         for ( int i = 0; i < hb.To(w).isize( ); i++ )
                         {    e = hb.ITo( w, i );
                              int pstart = start + hb.EdgeLengthKmers(e);
                              if ( pstart >= 0 )
                                   locs.push( e, pstart );    }    }    }    }    }
     UniqueSort(locs);
     if ( locs.empty( ) )
     {    if ( p.size( ) > 0 )
          {    status = path_improver::old_better;
               if (pimp.show_old_better)
               {    *pout << "\nid = " << id << ", old better" << std::endl;
                    *pout << "was placed, but can't find seed" << std::endl;    }    }
          else
          {    status = path_improver::same;
               if (pimp.show_same) *pout << "\nid = " << id << ", same" << std::endl;    }
          FinalPrint(pout);
          return;    }

     // If there is only one possible path and it agrees with the existing path,
     // don't do anything.

     if ( locs.solo( ) && p.size( ) == 1 && locs[0].first == p[0]
          && locs[0].second == p.getOffset( )
          && hb.EdgeLengthBases( p[0] ) - p.getOffset( ) >= b.isize( ) )
     {    status = path_improver::same;
          if (pimp.show_same) *pout << "\nid = " << id << ", same" << std::endl;
          FinalPrint(pout);
          return;    }

     // Extend seeds.

     vec<vec<int>> exts;
     vec<int> starts;
     vec<int> exts_len;  // bases left to place (<=0 if all used)
     for ( int i = 0; i < locs.isize( ); i++ )
     {    int e = locs[i].first, start = locs[i].second;
          int v = to_right[e];
          vec<int> e1 = {e};
          exts.push_back(e1);
          exts_len.push_back( b.isize( ) - (hb.EdgeLengthBases(e) - start) );
          starts.push_back(start);    }
     for ( int j = 0; j < exts.isize( ); j++ )
     {    if ( j > max_locs2 )  // too many extensions
                                // - abandon attempts to find new paths
	  {    if ( p.size( ) > 0 )
               {    status = path_improver::old_better;
                    if (pimp.show_old_better)
                    {    *pout << "\nid = " << id << ", old better" << std::endl;
                         *pout << "was placed, but too many extensions"
                              << std::endl;    }    }
               else
               {    status = path_improver::same;
                    if (pimp.show_same)
                    {    *pout << "\nid = " << id << ", same" << std::endl;
                         *pout << "wasn't placed, and hit too many extensions"
                              << std::endl;    }    }
               FinalPrint(pout);
               return;    }
          if ( exts_len[j] <= 0 )   // full alignment found - don't extend further
               continue;
          int y = to_right[ exts[j].back( ) ];
          if ( hb.From(y).empty( ) )  // graph ends, can't extend further
          {    if ( p.size( ) > 0 )
               {    status = path_improver::old_better;
                    if (pimp.show_old_better)
                    {    *pout << "\nid = " << id << ", old better" << std::endl;
                         *pout << "was placed, but found dead end" << std::endl;    }    }
               else
               {    status = path_improver::same;
                    if (pimp.show_same)
                    {    *pout << "\nid = " << id << ", same" << std::endl;
                         *pout << "wasn't placed, and found dead end"
                              << std::endl;    }    }
               FinalPrint(pout);
               return;    }
	  // continue extending to the right from current edge
          for ( int l = 0; l < hb.From(y).isize( ); l++ )
          {    vec<int> e = exts[j];
               int n = hb.EdgeObjectIndexByIndexFrom( y, l );
               e.push_back(n);
               exts.push_back(e);
               exts_len.push_back( exts_len[j] - hb.EdgeLengthKmers(n) );
               starts.push_back( starts[j] );    }    }
     // discard partial alignments
     vec<Bool> to_del( exts.size( ), False );
     for ( int j = 0; j < exts.isize( ); j++ )
          if ( exts_len[j] > 0 ) to_del[j] = True;
     EraseIf(exts, to_del), EraseIf(exts_len, to_del), EraseIf(starts, to_del);

     // Nothing to do if there are no paths.

     if ( exts.empty( ) )
     {    if ( p.size( ) > 0 )
          {    status = path_improver::old_better;
               if (pimp.show_old_better)
               {    *pout << "\nid = " << id << ", old better" << std::endl;
                    *pout << "was placed, but found no extensions" << std::endl;    }    }
          else
          {    status = path_improver::same;
               if (pimp.show_same)
               {    *pout << "\nid = " << id << ", same" << std::endl;
                    *pout << "wasn't placed, and found no extensions"
                         << std::endl;    }    }
          FinalPrint(pout);
          return;    }

     // Evaluate extensions by scoring using quality scores

     int nexts = exts.size( );
     vec<int> qsum( nexts, 0 ), qsum2( nexts, 0 ), qmax( nexts, 0 );
     for ( int j = 0; j < exts.isize( ); j++ )
     {    const vec<int>& e = exts[j];
          basevector E = hb.Cat(e);
          int start = starts[j];
          for ( int m = 0; m < b.isize( ); m++ )
          {    if ( b[m] != E[start+m] )
               {    qsum[j] += q[m];
                    if ( q[m] > 2 ) qsum2[j] += q[m];
                    qmax[j] = Max( qmax[j], (int) q[m] );    }    }    }
     SortSync( qsum, qsum2, qmax, exts, starts );

     // Delete weak extensions.

     for ( int j = 1; j < exts.isize( ); j++ )
     {    if ( qsum[j] - qsum[0] >= min_gain )
          {    qsum.resize(j), qsum2.resize(j), qmax.resize(j);
               exts.resize(j), starts.resize(j);
               break;    }    }

     // If there are two extensions, differing at only a single Q2 base, and the
     // read agrees with one in a large enough neighborhood, delete the other.

     if ( exts.size( ) == 2 && qsum[1] - qsum[0] == 2 )
     {    vec<int> diffs;
          basevector E1 = hb.Cat( exts[0] ), E2 = hb.Cat( exts[1] );
          for ( int m = 0; m < b.isize( ); m++ )
          {    if ( b[m] == E1[starts[0]+m] && b[m] != E2[starts[1]+m] )
                    diffs.push_back(m);    }
          if ( diffs.solo( ) )
          {    int d = diffs[0];
               if ( d >= flank && d < b.isize( ) - flank )
               {    int mis = 0;
                    for ( int m = d - flank; m <= d + flank; m++ )
                         if ( b[m] != E1[starts[0]+m] ) mis++;
                    if ( mis == 0 )
                    {    qsum.resize(1), qsum2.resize(1), qmax.resize(1);
                         exts.resize(1), starts.resize(1);    }    }    }    }

     // Check to see if the best path agrees with what we started with.

     ReadPath pnew( starts[0], exts[0] );
     if ( p == pnew )
     {    status = path_improver::same;
          if (pimp.show_same)
          {    *pout << "\nid = " << id << ", same" << std::endl;
               *pout << "best path agrees with old path" << std::endl;    }
          FinalPrint(pout);
          return;    }

     // Check alignment for goodness.  Require window of 60 with at most
     // 10% mismatches.

     Bool good = False;
     {    const vec<int>& e = exts[0];
          basevector E = hb.Cat( exts[0] );
          int start = starts[0], mis = 0, count = 0;
          for ( int m = 0; m < b.isize( ); m++ )
          {    count++;
               if ( b[m] != E[start+m] ) mis++;
               if ( count > window )
               {    if ( b[m-window] != E[start+m-window] ) mis--;
                    count--;    }
               if ( count == window && mis <= max_mis )
               {    good = True;
                    break;    }    }    }
     if ( !good && p.size( ) == 0 )
     {    status = path_improver::same;
          if (pimp.show_same)
          {    *pout << "\nid = " << id << ", same" << std::endl;
               *pout << "wasn't placed, but new placement looks bad" << std::endl;    }
          FinalPrint(pout);
          return;    }

     // If extensions are inconsistent, give up.

     for ( int m = 1; m < exts.isize( ); m++ )
     {    if ( starts[m] != starts[0] || exts[m][0] != exts[0][0] )
          {    if ( p.size( ) > 0 )
               {    status = path_improver::old_better;
                    if (pimp.show_old_better)
                    {    *pout << "\nid = " << id << ", old better" << std::endl;
                         *pout << "was placed, but extensions inconsistent"
                              << std::endl;    }    }
               else
               {    status = path_improver::same;
                    if (pimp.show_same)
                    {    *pout << "\nid = " << id << ", same" << std::endl;
                         *pout << "wasn't placed, but new extensions are "
                              << "inconsistent" << std::endl;    }    }
               FinalPrint(pout);
               return;    }    }

     // Compute core extension.

     vec<int> core;
     {    int m;
          for ( m = 1; m < exts.isize( ); m++ )
               if ( qsum[m] - qsum[0] >= min_gain ) break;
          for ( int j = 0; j < exts[0].isize( ); j++ )
          {    Bool mis = False;
               for ( int l = 1; l < m; l++ )
               {    if ( j >= exts[l].isize( ) || exts[l][j] != exts[0][j] )
                    {    mis = True;
                         break;    }    }
               if (mis) break;
               core.push_back( exts[0][j] );    }    }
     pnew = ReadPath( starts[0], core );

     // Try extending core extension backwards.

     int v = to_left[ core[0] ];
     if ( hb.To(v).solo( ) )
     {    int e = hb.ITo( v, 0 );
          int ne = hb.EdgeLengthKmers(e);
          int start = starts[0] + ne;
          if ( start <= hb.EdgeLengthBases(e) - min_back_over )
          {    core.push_front(e);
               starts[0] = start;    }    }

     // Decision criteria:
     // 1. Don't do anything unless qsum for full alignment <= 150.
     // 2. If there is an existing alignment, core must go at least as far.
     // 2. If there is an existing path, the path for the core up to the same
     //    point has to be at least as good.
     //
     // If core extension agrees with p, then nothing to do.

     vec<int> py;
     for ( int j = 0; j < (int) p.size( ); j++ )
          py.push_back( p[j] );
     if ( core == py && starts[0] == p.getOffset( ) )
     {    status = path_improver::same;
          if (pimp.show_same)
          {    *pout << "\nid = " << id << ", same" << std::endl;
               *pout << "new = old" << std::endl;    }
          FinalPrint(pout);
          return;    }

     // Other cases where done.

     if ( qsum[0] > max_ok )
     {    if ( p.size( ) == 0 )
          {    status = path_improver::same;
               if (pimp.show_same)
               {    *pout << "\nid = " << id << ", same" << std::endl;
                    *pout << "wasn't placed, but new has high qsum"
                         << std::endl;    }    }
          else // calling old better, but not necessarily true
          {    status = path_improver::old_better;
               if (pimp.show_old_better)
               {    *pout << "\nid = " << id << ", old better" << std::endl;
                    *pout << "new path has high qsum (but perhaps the "
                         << "old one does too)" << std::endl;    }    }
          FinalPrint(pout);
          return;    }

     if ( p.size( ) == 0 )
     {    status = path_improver::new_better;
          p = pnew;
          if (pimp.show_new_better)
          {    *pout << "\nid = " << id << ", new better" << std::endl;
               *pout << "old: empty" << std::endl;
               *pout << "new: " << starts[0] << ": " << printSeq(core) << std::endl;    }
          FinalPrint(pout);
          return;    }

     // Compare core to existing alignment.

     int old_stop = -1, old_qsum = 0, new_start = 0, new_stop = -1, new_qsum = 0;
     int old_start = p.getOffset( ) >= 0 ? 0 : -p.getOffset( );
     {    vec<int> e;
          for ( int j = 0; j < (int) p.size( ); j++ )
               e.push_back( p[j] );
          basevector E = hb.Cat(e);
          int start = p.getOffset( ), m = 0;
          for ( m = 0; m < b.isize( ); m++ )
          {    if ( start+m < 0 ) continue;
               if ( start+m == E.isize( ) ) break;
               if ( b[m] != E[start+m] ) old_qsum += q[m];    }
          old_stop = m;    }
     {    vec<int> e;
          for ( int j = 0; j < core.isize( ); j++ )
               e.push_back( core[j] );
          basevector E = hb.Cat(e);
          int start = starts[0], m = 0;
          for ( m = 0; m < b.isize( ); m++ )
          {    if ( start+m == E.isize( ) ) break;
               if ( b[m] != E[start+m] )
                    if ( old_start <= m && m < old_stop ) new_qsum += q[m];    }
          new_stop = m;    }

     // Make decision.  Note that this doesn't test 'extra' edges tacked on
     // at beginning.

     Bool print = False;
     if ( new_start > old_start || new_stop < old_stop )
     {    status = path_improver::old_better; // Not necessarily!
          if (pimp.show_old_better)
               *pout << "\nid = " << id << ", old better" << std::endl;    }
     else if ( new_start == old_start && new_stop == old_stop )
     {    if ( new_qsum < old_qsum )
          {    status = path_improver::new_better;
               p = pnew;
               if (pimp.show_new_better)
               {    *pout << "\nid = " << id << ", new better" << std::endl;
                    print = True;    }    }
          else if ( old_qsum < new_qsum )
          {    status = path_improver::old_better;
               if (pimp.show_old_better)
               {    *pout << "\nid = " << id << ", old better" << std::endl;
                    print = True;    }    }
          else
          {    status = path_improver::same;
               if (pimp.show_same)
               {    *pout << "\nid = " << id << ", same" << std::endl;
                    print = True;    }    }    }
     else if ( new_qsum <= old_qsum )
     {    status = path_improver::new_better;
          p = pnew;
          if (pimp.show_new_better)
          {    *pout << "\nid = " << id << ", new better" << std::endl;
               print = True;    }    }
     else
     {    status = path_improver::indet;
          if (pimp.show_indet)
          {    *pout << "\nid = " << id << ", indeterminate" << std::endl;
               print = True;    }    }
     if (print)
     {    *pout << "old: " << p.getOffset( ) << ":" << printSeq(p) << std::endl;
          vec<int> e;
          for ( int j = 0; j < (int) p.size( ); j++ )
               e.push_back( p[j] );
          basevector E = hb.Cat(e);
          int start = p.getOffset( );
          {    int qsum = 0, qmax = 0;
               for ( int m = 0; m < b.isize( ); m++ )
               {    if ( start+m < 0 ) continue;
                    if ( start+m == E.isize( ) ) break;
                    if ( b[m] != E[start+m] )
                    {    qsum += q[m];
                         qmax = Max( qmax, (int) q[m] );    }    }
               *pout << "old has qsum = " << qsum << ", qmax = "
               << qmax << std::endl;    }
          int xip = ( xi % 2 == 0 ? xi + 1 : xi - 1 );
          const ReadPath& pp = paths[xip];
          vec<int> ppx;
          for ( int j = (int) pp.size( ) - 1; j >= 0; j-- )
               ppx.push_back( inv[ pp[j] ] );
          *pout << "inv(partner(old)): " << printSeq(ppx) << std::endl;
          vec<int> core_inv;
          for ( int j = (int) core.size( ) - 1; j >= 0; j-- )
               core_inv.push_back( inv[ core[j] ] );
          *pout << "new: " << starts[0] << ": " << printSeq(core)
               << " (inv = " << printSeq(core_inv) << ")" << std::endl;
          *pout << "new qsum = " << printSeq(qsum) << " (qsum2 = "
               << printSeq(qsum2) << ", qmax = " << printSeq(qmax) << ")" << std::endl;
          *pout << "comparison qsums: old = " << old_qsum << ", new = "
               << new_qsum << std::endl;
          if ( !good ) *pout << "new alignment looks like junk" << std::endl;
          if (pimp.print_align)
          {    basevector E = hb.Cat(core);
               avector<int> gaps(1), lengths(1);
               gaps(0) = 0, lengths(0) = b.size( );
               align a( 0, starts[0], gaps, lengths );
               PrintVisualAlignment( True, *pout, b, E, a, q );    }    }
     FinalPrint(pout);    }

template<int L> void ImprovePathsCoreCore( const vec<int>& to_left,
     const vec<int>& to_right, const vec< triple<kmer<L>,int,int> >& kmers_plus,
     const vec< std::pair<int64_t, std::pair<int,int> > >& locsx,
     const vec<int>& rstarts,
     ReadPathVec& paths, const HyperBasevector& hb,
     const vec<int>& inv, const vecbasevector& bases, const VecPQVec& quals,
     const vec<int64_t>& ids, const path_improver& pimp )
{
     // Improve paths.

     Bool track_results = pimp.Logging( );
     int count_old_better = 0, count_new_better = 0;
     int count_same = 0, count_indet = 0;
     #pragma omp parallel for
     for ( int64_t id = 0; id < (int64_t) bases.size( ); id++ )
     {    path_improver::path_status status;
          int64_t true_id = ( ids.empty( ) ? id : ids[id] );

          ImprovePath( rstarts, locsx, paths, id, paths[id], true_id, hb, inv,
               to_left, to_right, bases[id], quals.begin()[id], kmers_plus, pimp,
               status );

          if (track_results) // slow
          {
               #pragma omp critical
               {    if ( status == path_improver::old_better ) count_old_better++;
                    if ( status == path_improver::new_better ) count_new_better++;
                    if ( status == path_improver::same ) count_same++;
                    if ( status == path_improver::indet ) count_indet++;
                         }    }    }

     // Report results.

     if (track_results)
     {    std::cout << "\n";
          PRINT(count_old_better);
          PRINT(count_new_better);
          PRINT(count_same);
          PRINT(count_indet);    }    }

template<int L> void BuildLookup( vec< triple<kmer<L>,int,int> >& kmers_plus,
     const HyperBasevector& hb )
{
     // Build L-mer lookup table.

     const int K = hb.K( );
     vecbasevector edges( hb.E( ) );
     for ( int e = 0; e < hb.E( ); e++ )
     {    edges[e] = hb.EdgeObject(e);
          int ne = edges[e].size( );
          if ( ne > 0 ) edges[e].resize( ne + L - K );    }
     MakeKmerLookup0( edges, kmers_plus );    }

void ImprovePaths( ReadPathVec& paths, const HyperBasevector& hb,
     const vec<int>& inv, const vecbasevector& bases, const VecPQVec& quals,
     const vec<int64_t>& ids, const path_improver& pimp,
     const Bool IMPROVE_PATHS_LARGE, const Bool BETSYBOB )
{
     // Build indices.

     //std::cout << Date( ) << ": improving paths" << std::endl;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);

     // First try large K.  Output is seeds.

     vec< std::pair<int64_t, std::pair<int,int> > > locsx;
     const int K = hb.K( );
     if (IMPROVE_PATHS_LARGE)
     {
          // Hash large kmers.

          const int L = 200;
          ForceAssertEq( L, K );
          vec< triple<uint64_t,int,int> > kmers_plus;
          std::cout << Date( ) << ": forming edges" << std::endl;
          vecbasevector edges( hb.E( ) );
          for ( int e = 0; e < hb.E( ); e++ )
          {    edges[e] = hb.EdgeObject(e);
               int ne = edges[e].size( );
               if ( ne > 0 ) edges[e].resize( ne + L - K );    }
          std::cout << Date( ) << ": hashing" << std::endl;
          MakeKmerLookup0H<L>( edges, kmers_plus );

          // Find unmapped reads.

          std::cout << Date( ) << ": finding unmapped reads" << std::endl;
          vec<int64_t> un;
          for ( int64_t i = 0; i < (int64_t) bases.size( ); i++ )
               if ( paths[i].size( ) == 0 ) un.push_back(i);

          // For unmapped reads, find seeds based on large kmers.

          std::cout << Date( ) << ": traversing reads" << std::endl;
          #pragma omp parallel for
          for ( int64_t q = 0; q < un.isize( ); q++ )
          {    int64_t id = un[q];
               for ( int j = 0; j <= bases[id].isize( ) - L; j++ )
               {    if ( j > 0 ) break;
                    kmer<L> x;
                    x.SetToSubOf( bases[id], j );
                    uint64_t hash = FNV1a( x.Ints( ), x.Ints( ) + (K+15)/16 );
                    int64_t low = LowerBound1( kmers_plus, hash );
                    int64_t high = UpperBound1( kmers_plus, hash );
                    vec<int64_t> ms;
                    for ( int64_t m = low; m < high; m++ )
                    {    int e = kmers_plus[m].second;
                         int epos = kmers_plus[m].third;
                         Bool mismatch = False;
                         for ( int l = 0; l < L; l++ )
                         {    if ( hb.EdgeObject(e)[epos+l] != bases[id][l] )
                              {    mismatch = True;
                                   break;     }    }
                         if ( !mismatch ) ms.push_back(m);    }
                    if ( ms.solo( ) )
                    {    int64_t m = ms[0];
                         int e = kmers_plus[m].second, epos = kmers_plus[m].third;
                         #pragma omp critical
                         {     locsx.push(
                                   id, std::make_pair( e, epos ) );    }    }    }    }
          ParallelSort(locsx);
          std::cout << Date( ) << ": done" << std::endl;    }

     // Two passes of improvement.

     int npasses = 1;
     if (BETSYBOB) npasses = 3;
     for ( int pass = 1; pass <= npasses; pass++ )
     {    if ( pass == 1 )
          {    const int L = 20;
               const vec<int> rstarts = {0,20,40,60};
               vec< triple<kmer<L>,int,int> > kmers_plus;
               BuildLookup( kmers_plus, hb );
               ImprovePathsCoreCore( to_left, to_right, kmers_plus, locsx,
                    rstarts, paths, hb, inv, bases, quals, ids, pimp );    }
          if ( pass == 2 )
          {    const int L = 40;
               const vec<int> rstarts = {0};
               vec< triple<kmer<L>,int,int> > kmers_plus;
               BuildLookup( kmers_plus, hb );
               ImprovePathsCoreCore( to_left, to_right, kmers_plus, locsx,
                    rstarts, paths, hb, inv, bases, quals, ids, pimp );    }
          if ( pass == 3 )
          {    const int L = 80;
               const vec<int> rstarts = {0,80};
               vec< triple<kmer<L>,int,int> > kmers_plus;
               BuildLookup( kmers_plus, hb );
               ImprovePathsCoreCore( to_left, to_right, kmers_plus, locsx,
                    rstarts, paths, hb, inv, bases, quals, ids, pimp );    }    }
     //std::cout << Date( ) << ": done" << std::endl;
}
