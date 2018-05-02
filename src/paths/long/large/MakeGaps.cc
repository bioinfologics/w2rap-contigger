///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/MakeGaps.h"
#include "paths/long/large/GapToyTools.h"

void MakeGaps( HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths,
     VecULongVec& edgeToPathIds, const int MIN_LINE, const int MIN_LINK_COUNT, 
     const String& work_dir, const String& FIN, const Bool verbose,
     const Bool GAP_CLEANUP )
{
     // Set up data structures.

     double clock = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( work_dir + "/" + FIN + ".fin.lines", &lines );
     vec<int> llens, npairs;
     GetLineLengths( hb, lines, llens );
     BinaryReader::readFile( work_dir + "/" + FIN + ".fin.lines.npairs", &npairs );
     vec<double> cov( lines.size( ) );
     for ( int i = 0; i < lines.isize( ); i++ )
          cov[i] = 100.0 * double(npairs[i]) / double(llens[i]);
     vec<int> tol;
     GetTol( hb, lines, tol );

     // Heuristics.

     const int max_hang = 800;
     const int max_depth = 2;
     const int max_int = 1500; // could try lowering
     const int passes = 3;
     const double max_cov_pc_off = 20.0;
     const int max_line_to_ignore = 500;

     // Define edge groups.  These are bunches of edges that are near sinks and
     // sources.  Via 'tom', edges near such ends are mapped back to 'primary'
     // edges that are farther from the ends, their distance to the end is noted
     // via dist_to_end, and they are flagged via sink_like or source_like.

     // std::cout << Date( ) << ": defining edge groups" << std::endl;
     int nobj = hb.EdgeObjectCount( );
     vec<int> tom( nobj, vec<int>::IDENTITY );
     vec<Bool> sink_like( nobj, False ), source_like( nobj, False );
     vec<int> dist_to_end( nobj, 0 );
     for ( int e = 0; e < nobj; e++ )
     {
         if ( hb.From( to_right[e] ).empty( ) ) {
             sink_like[e] = True;
         }
         if ( hb.To( to_left[e] ).empty( ) )  {
              source_like[e] = True;
         }
     }
     for ( int pass = 1; pass <= passes; pass++ )
     {    for ( int zpass = 1; zpass <= 2; zpass++ )
          {    hb.Reverse( );
               vec<int> to_left, to_right;
               hb.ToLeft(to_left), hb.ToRight(to_right);

               /*
               // trailing cycle???
               for ( int e = 0; e < nobj; e++ )
               {    int v = to_right[e];
                    if ( hb.To(v).size( ) != 2 ) continue;
                    if ( hb.From(v).size( ) != 1 ) continue;
                    if ( hb.From(v)[0] != v ) continue;
                    int ec = hb.EdgeObjectIndexByIndexFrom( v, 0 );
                    if ( hb.Kmers(ec) > max_hang ) continue;
               */

               for ( int e = 0; e < nobj; e++ )
               {
                    int v = to_right[e];
                    if ( hb.From(v).size( ) != 2 || hb.To(v).size( ) != 1 ) continue;

                    int e1 = hb.EdgeObjectIndexByIndexFrom( v, 0 );
                    int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1 );
                    int w1 = to_right[e1], w2 = to_right[e2];

                    if ( zpass == 2 && ( !sink_like[e1] || !sink_like[e2] ) ) continue;
                    if ( zpass == 1 && ( !source_like[e1] || !source_like[e2] ) ) continue;
                    if ( w1 == w2 && hb.To(w1).size( ) != 2 ) continue;
                    if ( w1 != w2 && ( !hb.To(w1).solo( ) || !hb.To(w2).solo( ) ) ) continue;

                    int d1 = hb.Kmers(e1) + dist_to_end[e1];
                    int d2 = hb.Kmers(e2) + dist_to_end[e2];
                    if ( d1 > max_hang || d2 > max_hang ) continue;
                    if ( zpass == 2 ) sink_like[e] = True;
                    else source_like[e] = True;
                    dist_to_end[e] = Max( d1, d2 );
                    tom[e1] = tom[e], tom[e2] = tom[e];
               }

               for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               {    int v = to_right[e];
                    if ( hb.From(v).size( ) != 2 || hb.To(v).size( ) != 1 ) continue;
                    int e1 = hb.EdgeObjectIndexByIndexFrom( v, 0 );
                    int e2 = hb.EdgeObjectIndexByIndexFrom( v, 1 );
                    int w1 = to_right[e1], w2 = to_right[e2];
                    if ( w1 != w2 ) continue;
                    if ( hb.To(w1).size( ) != 2 || !hb.From(w1).solo( ) ) continue;
                    int z = hb.From(w1)[0];
                    if ( !hb.To(z).solo( ) ) continue;
                    int e3 = hb.EdgeObjectIndexByIndexFrom( w1, 0 );
                    if ( zpass == 2 && !sink_like[e3] ) continue;
                    if ( zpass == 1 && !source_like[e3] ) continue;
                    int d1 = hb.Kmers(e1) + hb.Kmers(e3) + dist_to_end[e3];
                    int d2 = hb.Kmers(e2) + hb.Kmers(e3) + dist_to_end[e3];
                    if ( d1 > max_hang || d2 > max_hang ) continue;
                    if ( zpass == 2 ) sink_like[e] = True;
                    else source_like[e] = True;
                    dist_to_end[e] = Max( d1, d2 );
                    tom[e1] = tom[e], tom[e2] = tom[e], tom[e3] = tom[e];    
                         }    }    }

     // Define edges that are near each other.

     // std::cout << Date( ) << ": computing nears" << std::endl;
     vec< std::pair<int,int> > nears;
     vec< vec<int> > nears1(nobj), nears2(nobj);
     nears.reserve( int64_t( round( 0.7 * paths.size( ) ) ) );
     // PRINT( nears.capacity( ) );
     const int nbatches = 100;
     int64_t npids = paths.size( ) / 2;
     #pragma omp parallel for
     for ( int bi = 0; bi < nbatches; bi++ )
     {    vec< std::pair<int,int> > n;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int64_t pid1 = ( bi * npids ) / nbatches;
               int64_t pid2 = ( (bi+1) * npids ) / nbatches;
               for ( int64_t pid = pid1; pid < pid2; pid++ )
               {    int64_t id1 = 2*pid, id2 = 2*pid+1;
                    if ( paths[id1].size( ) == 0 || paths[id2].size( ) == 0 ) 
                    continue;
                    vec<int> x, y;
                    for ( int64_t j = 0; j < (int64_t) paths[id1].size( ); j++ )
                         x.push_back( paths[id1][j] );
                    for ( int64_t j = 0; j < (int64_t) paths[id2].size( ); j++ )
                         y.push_back( paths[id2][j] );
                    y.ReverseMe( );
                    for ( int j = 0; j < y.isize( ); j++ )
                         y[j] = inv[ y[j] ];
                    if ( pass == 2 )
                    {    swap( x, y );
                         x.ReverseMe( ), y.ReverseMe( );
                         for ( int j = 0; j < x.isize( ); j++ )
                              x[j] = inv[ x[j] ];
                         for ( int j = 0; j < y.isize( ); j++ )
                              y[j] = inv[ y[j] ];    }
     
                    // Apply tom.
     
                    for ( int j = 0; j < x.isize( ); j++ )
                         x[j] = tom[ x[j] ];
                    for ( int j = 0; j < y.isize( ); j++ )
                         y[j] = tom[ y[j] ];
     
                    // Compress consecutive duplicates.

                    int t = 0;
                    for ( int j = 1; j < x.isize( ); j++ )
                         if ( x[j] != x[t] ) x[++t] = x[j];
                    x.resize(t+1);
                    t = 0;
                    for ( int j = 1; j < y.isize( ); j++ )
                         if ( y[j] != y[t] ) y[++t] = y[j];
                    y.resize(t+1);

                    // Ignore little stuff.

                    t = 0;
                    for ( int j = 0; j < x.isize( ); j++ )
                         if ( llens[tol[x[j]]] > max_line_to_ignore ) x[t++] = x[j];
                    x.resize(t);
                    t = 0;
                    for ( int j = 0; j < y.isize( ); j++ )
                         if ( llens[tol[y[j]]] > max_line_to_ignore ) y[t++] = y[j];
                    y.resize(t);
               
                    // Generate nears.
     
                    vec<int> ys(y);
                    UniqueSort(ys);
                    for ( int j1 = 0; j1 < x.isize( ); j1++ )
                    {    if ( BinMember( ys, x[j1] ) ) continue;
                         for ( int j2 = 0; j2 < y.isize( ); j2++ )
                         {    int e1 = x[j1], e2 = y[j2];
                              if ( e1 == e2 ) continue;
                              n.push( e1, e2 );    }    }    }    }
          #pragma omp critical
          {    nears.append(n);
               for ( int i = 0; i < n.isize( ); i++ )
               {    int e1 = n[i].first, e2 = n[i].second;
                    nears1[e1].push_back(e2);
                    nears2[e2].push_back(e1);    }    }    }
     // PRINT2( nears.size( ), paths.size( ) );
     // std::cout << Date( ) << ": sorting nears" << std::endl;
     ParallelSort(nears);
     #pragma omp parallel for
     for ( int e = 0; e < nobj; e++ )
     {    Sort( nears1[e] ), Sort( nears2[e] );    }

     // Find good links.

     // std::cout << Date( ) << ": finding links" << std::endl;
     int events = 0;
     vec< std::pair<int,int> > links;
     vec<int> counts;
     const int batches = 1000;
     vec<int64_t> bstart(batches+1);
     for ( int64_t i = 0; i <= batches; i++ )
          bstart[i] = ( (int64_t) nears.size( ) * i ) / batches;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t i = 1; i < batches; i++ )
     {    int64_t& s = bstart[i];
          while( s > 0 && nears[s].first == nears[s-1].first )
          {    s--;    }    }
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches; bi++ )
     {    for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j = nears.NextDiff(i);
               {    int e1 = nears[i].first, e2 = nears[i].second;
     
                    // Check to see if e2 is within max_depth of e1 or its 
                    // predecessor.

                    Bool close = False;
                    vec<int> x = {e1}, d = {-1}, k = {0};
                    if ( hb.To( to_left[e1] ).solo( ) )
                    {    int e = hb.EdgeObjectIndexByIndexTo( to_left[e1], 0 );
                         x.push_back(e), d.push_back(-1), k.push_back(0);    }
                    for ( int j = 0; j < x.isize( ); j++ )
                    {    int e = x[j];
                         if ( e == e2 )
                         {    close = True;
                              break;    }
                         if ( k[j] > max_int || d[j] == max_depth ) continue;
                         int v = to_right[e], w = to_left[e];
                         for ( int l = 0; l < hb.From(v).isize( ); l++ )
                         {    int e = hb.EdgeObjectIndexByIndexFrom( v, l );
                              x.push_back(e), d.push_back( d[j] + 1 );
                              k.push_back( k[j] + hb.EdgeLengthKmers(e) );    }    
                         for ( int l = 0; l < hb.To(w).isize( ); l++ )
                         {    int e = hb.EdgeObjectIndexByIndexTo( w, l );
                              x.push_back(e), d.push_back( d[j] + 1 );
                              k.push_back( k[j] + hb.EdgeLengthKmers(e) );    }    }
                    if ( !close )
                    {    
                         #pragma omp critical
                         {    links.push( tom[e1], tom[e2] );
                              counts.push_back(j-i);    }    }    }
               i = j - 1;    }    }

     // Sort links and note objects that appear on the right of two different links.

     // std::cout << Date( ) << ": sorting links" << std::endl;
     ParallelSortSync( links, counts );

     // Compute edge groups.

     // std::cout << Date( ) << ": computing edge groups" << std::endl;

     vec<vec<int>> gp(nobj);
    if (verbose) {
        for (int e = 0; e < nobj; e++)
            gp[tom[e]].push_back(e);
    }
     // Finalize links.

     // std::cout << Date( ) << ": finalizing links" << std::endl;
     vec< std::pair<int,int> > accepted;
     for ( int i = 0; i < links.isize( ); i++ )
     {    int e1 = links[i].first, e2 = links[i].second;

          // Require enough links, and some length on both sides.

          if ( counts[i] < MIN_LINK_COUNT ) continue;
          if ( llens[tol[e1]] < MIN_LINE || llens[tol[e2]] < MIN_LINE ) continue;

          // Require about the same coverage on both sides.

          double c1 = cov[tol[e1]], c2 = cov[tol[e2]];
          if ( c1 < c2 ) std::swap( c1, c2 );
          if ( c1/c2 - 1.0 > max_cov_pc_off/100.0 ) continue;

          // Must be winner.

          int max_alt = 0;
          for ( int l = 0; l < nears1[e1].isize( ); l++ )
          {    int m = nears1[e1].NextDiff(l);
               max_alt = Max( max_alt, m - l );
               l = m - 1;    }
          for ( int l = 0; l < nears2[e2].isize( ); l++ )
          {    int m = nears2[e2].NextDiff(l);
               max_alt = Max( max_alt, m - l );
               l = m - 1;    }
          if ( max_alt > counts[i] ) continue;

          // Advance past simple bubbles, and require canonical.

          int e1x = e1, e2x = e2;
          for ( int p = 0; p < passes; p++ )
          {    int v = to_right[e1x];
               if ( !hb.To(v).solo( ) || hb.From(v).size( ) != 2 ) break;
               if ( hb.From(v)[0] != hb.From(v)[1] ) break;
               int w = hb.From(v)[0];
               if ( hb.To(w).size( ) != 2 || !hb.From(w).solo( ) ) break;
               e1x = hb.IFrom( w, 0 );    }
          for ( int p = 0; p < passes; p++ )
          {    int v = to_left[e2x];
               if ( !hb.From(v).solo( ) || hb.To(v).size( ) != 2 ) break;
               if ( hb.To(v)[0] != hb.To(v)[1] ) break;
               int w = hb.To(v)[0];
               if ( hb.From(w).size( ) != 2 || !hb.To(w).solo( ) ) break;
               e2x = hb.ITo( w, 0 );    }
          int l1 = tol[e1x], l2 = tol[e2x];
          if ( lines[l1].back( )[0][0] != e1x ) continue;
          if ( lines[l2].front( )[0][0] != e2x ) continue;

          // Mark as accepted, but note that some of these are discarded later.

          accepted.push( e1, e2 );

          // Print.

          events++;
          if (verbose)
          {    std::cout << "\n" << events << ". [" << counts[i] << "] " << e1 
                    << " --> " << e2 << std::endl;
               std::cout << "lines: " << tol[e1] << "[len=" << llens[tol[e1]] << "], " 
                    << tol[e2] << "[len=" << llens[tol[e2]] << "]" << std::endl;
               std::cout << "lefts: ";
               for ( int l = 0; l < nears1[e1].isize( ); l++ )
               {    int m = nears1[e1].NextDiff(l);
                    if ( l > 0 ) std::cout << ",";
                    std::cout << nears1[e1][l];
                    if ( m - l > 1 ) std::cout << "^" << m - l;
                    l = m - 1;    }
               std::cout << "\nrights: ";
               for ( int l = 0; l < nears2[e2].isize( ); l++ )
               {    int m = nears2[e2].NextDiff(l);
                    if ( l > 0 ) std::cout << ",";
                    std::cout << nears2[e2][l];
                    if ( m - l > 1 ) std::cout << "^" << m - l;
                    l = m - 1;    }
               std::cout << "\nleft edge group: " << printSeq( gp[e1] ) << std::endl;
               std::cout << "right edge group: " << printSeq( gp[e2] ) << std::endl;    }    }

     // Unaccept links that are not one-to-one.

     vec<int> a1, a2;
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    a1.push_back( accepted[i].first ), a2.push_back( accepted[i].second );    }
     Sort(a1), Sort(a2);
     vec<Bool> adel( accepted.size( ), False );
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    int low1 = LowerBound( a1, accepted[i].first );
          int high1 = UpperBound( a1, accepted[i].first );
          int low2 = LowerBound( a2, accepted[i].second );
          int high2 = UpperBound( a2, accepted[i].second );
          if ( high1 - low1 != 1 || high2 - low2 != 1 ) 
          {    if (verbose)
               {    std::cout << "deleting " << accepted[i].first << " --> "
                         << accepted[i].second << "\n";    }
               adel[i] = True;    }    }
     EraseIf( accepted, adel );

     // Advance past simple bubbles.

     for ( int i = 0; i < accepted.isize( ); i++ )
     {    int &e1 = accepted[i].first, &e2 = accepted[i].second;
          for ( int p = 0; p < passes; p++ )
          {    int v = to_right[e1];
               if ( !hb.To(v).solo( ) || hb.From(v).size( ) != 2 ) break;
               if ( hb.From(v)[0] != hb.From(v)[1] ) break;
               int w = hb.From(v)[0];
               if ( hb.To(w).size( ) != 2 || !hb.From(w).solo( ) ) break;
               e1 = hb.IFrom( w, 0 );    }
          for ( int p = 0; p < passes; p++ )
          {    int v = to_left[e2];
               if ( !hb.From(v).solo( ) || hb.To(v).size( ) != 2 ) break;
               if ( hb.To(v)[0] != hb.To(v)[1] ) break;
               int w = hb.To(v)[0];
               if ( hb.From(w).size( ) != 2 || !hb.To(w).solo( ) ) break;
               e2 = hb.ITo( w, 0 );    }    }

     // Force symmetry.  Maybe more complicated than it has to be.

     UniqueSort(accepted);
     int na = accepted.size( );
     vec<Bool> xa1( hb.E( ), False ), xa2( hb.E( ), False );
     for ( int i = 0; i < na; i++ )
     {    int e1 = accepted[i].first, e2 = accepted[i].second;
          xa1[e1] = True, xa2[e2] = True;    }
     vec< std::pair<int,int> > accepted2;
     vec<Bool> acdel( accepted.size( ), False );
     for ( int i = 0; i < na; i++ )
     {    int e1 = accepted[i].first, e2 = accepted[i].second;
          int re1 = inv[e1], re2 = inv[e2];
          if ( !BinMember( accepted, std::make_pair( re2, re1 ) ) )
          {    if ( !xa1[re2] && !xa2[re1] ) accepted2.push( re2, re1 );
               else acdel[i] = True;    }    }
     EraseIf( accepted, acdel );
     accepted.append(accepted2);
     UniqueSort(accepted);
     std::cout << Date( ) << ": deleting " << Sum(acdel) << " gaps and adding " 
          << accepted.isize( ) - na << " gaps to force symmetry" << std::endl;

     // Fix problem with overlinked edges.

     vec<int> cleft( hb.E( ), 0 ), cright( hb.E( ), 0 );
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    cleft[accepted[i].first]++;
          cright[accepted[i].second]++;    }
     vec<Bool> overlinked( hb.E( ) );
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( cleft[e] > 1 ) overlinked[e] = True;
          if ( cright[e] > 1 ) overlinked[e] = True;    }
     vec<Bool> del2( accepted.size( ), False );
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    int e1 = accepted[i].first, e2 = accepted[i].second;
          if ( overlinked[e1] || overlinked[e2] ) del2[i] = True;    }
     EraseIf( accepted, del2 );

     // Edit graph to add gap edges.

     for ( int i = 0; i < accepted.isize( ); i++ )
     {
         int e1 = accepted[i].first, e2 = accepted[i].second;
          int N = hb.N( );
          hb.AddVertices(2);
          hb.GiveEdgeNewToVx( e1, to_right[e1], N );
          to_right[e1] = N;
          hb.GiveEdgeNewFromVx( e2, to_left[e2], N+1 );
          to_left[e2] = N+1;
          hb.AddEdge( N, N+1, basevector( ) );

          // need to truncate paths
          std::map<size_t, vec<size_t>> to_remove;
          for ( auto const& pathid : edgeToPathIds[e1] ) {
              auto& path = paths[pathid];
              auto pos = std::find( path.begin(), path.end(), e1 );
              if ( pos != path.end() ) {
                  for ( auto itr = pos+1; itr != path.end(); ++itr )
                      to_remove[*itr].push_back(pathid);
                  path.erase(pos+1, path.end());
              }
          }

          for ( auto const& pathid : edgeToPathIds[e2] ) {
              auto& path = paths[pathid];
              auto pos = path.end();
              auto tmp = path.begin();
              // find the last occurrence of e2 in path
              while ( (tmp = std::find(tmp, path.end(), e2)) != path.end() ) {
                  pos = tmp++;
              }
              int offset = path.getOffset();
              if ( pos != path.end() && pos != path.begin() ) {
                  for ( auto itr = path.begin(); itr != pos; ++itr ) {
                      to_remove[*itr].push_back(pathid);
                      offset -= hb.EdgeLengthKmers(*itr);
                  }
                  path.erase( path.begin(), pos );
                  path.setOffset(offset);
              }
          }

          // need to update paths index
          for ( auto itr = to_remove.cbegin(); itr != to_remove.cend(); ++itr ) {
              size_t edge = itr->first;
              vec<size_t> const& skip = itr->second;
              auto& orig = edgeToPathIds[edge];
              ULongVec temp;
              std::copy_if( orig.begin(), orig.end(), std::back_inserter(temp),
                      [&skip]( unsigned long el ) { return !Member( skip, el ); } );
              orig = temp;
          }
     }

     // Fix the inversion.

     int nold = inv.size( );
     inv.resize( hb.E( ) );
     for ( int i = 0; i < accepted.isize( ); i++ )
     {    int e1 = accepted[i].first, e2 = accepted[i].second;
          int re1 = inv[e1], re2 = inv[e2];
          int ri = BinPosition( accepted, std::make_pair( re2, re1 ) );
          ForceAssertGe( ri, 0 );
          inv[ nold + i ] = nold + ri;    }

     // Clean up.

     if (GAP_CLEANUP) 
     {
         RemoveSmallComponents4( hb, True );
         Cleanup( hb, inv, paths );
         CleanupLoops( hb, inv, paths );
     }

     // Done.

     if (verbose) std::cout << "\n";
     // PRINT(events);
     if (verbose) std::cout << "\n";
     std::cout << Date( ) << ": done making gaps, time used = "
          << TimeSince(clock) << std::endl;
          }
