// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/long/CreateGenome.h"
//#include "paths/long/EvalByReads.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/DiscovarTools.h"
#include "random/Bernoulli.h"

void RemoveNegatives( vec<int>& p )
{    vec<Bool> to_delete( p.size( ), False );
     for ( int j = 0; j < p.isize( ); j++ )
          if ( p[j] == -1 ) to_delete[j] = True;
     EraseIf( p, to_delete );    }

Bool Compati( const vec<int>& X, const vec<int>& v, const int r, const int pos )
{    Bool compati = False;
     for ( int rr = 1; rr < v.isize( ) - 1; rr++ )
     {    if ( v[rr] != v[r] ) continue;
          Bool compatj = True;
          for ( int s = rr + 1; s < v.isize( ); s++ )
          {    int p = pos + s - rr;
               if ( p >= X.isize( ) ) break;
               if ( v[s] != X[p] )
               {    compatj = False;
                    break;    }    }
          for ( int s = rr - 1; s >= 0; s-- )
          {    int p = pos + s - rr;
               if ( p < 0 ) break;
               if ( v[s] != X[p] )
               {    compatj = False;
                    break;    }    }
          if (compatj)
          {    compati = True;
               break;    }    }
     return compati;    }

#if 0
void SupportedHyperBasevector::PullApart2( const double min_weight_split,
     const long_logging& logc )
{
     // Logging stuff.

     if (logc.STATUS_LOGGING) std::cout << Date( ) << ": starting PullApart2" << std::endl;
     int verbosity = logc.verb[ "PULL_APART2" ];
     double clock = WallClockTime( );

     // Iterate until no improvement.

     while(1)
     {    Bool progress = False;
          // ---->

     // Set up indices.

     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     vec< vec< std::pair<int,int> > > paths_index( EdgeObjectCount( ) );
     for (int i = 0; i < Paths( ).isize(); i++)    
     for (int j = 0; j < Path(i).isize(); j++)
          paths_index[ Path(i,j) ].push(i, j);

     // Keep track of proposed joins.

     vec< triple< vec<int>, vec<int>, vec<int> > > joins;
     vec< vec< vec<int> > > paths1, paths2;

     // Start with an edge e1, ending in a vertex v.

     int count = 0;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int e1 = 0; e1 < EdgeObjectCount( ); e1++ )
     {    
          if (logc.STATUS_LOGGING)
          {
               #pragma omp critical
               {    int dots1 = (100*count)/EdgeObjectCount( );
                    count++;
                    int dots2 = (100*count)/EdgeObjectCount( );
                    if ( dots2 > dots1 )
                    {    for ( int j = dots1; j < dots2; j++ )
                         {    std::cout << ".";
                              if ( (j+1) % 50 == 0 ) std::cout << "\n";
                              else if ( (j+1) % 10 == 0 ) std::cout << " ";    }
                         flush(std::cout);    }    }    }

          int v = to_right[e1];

          // Form an edge neighborhood of v.

          const int radius = 4;
          vec< std::pair<int,int> > ed;
          ed.push( e1, 0 );
          while(1)
          {    Bool progress = False;
               for ( int i = 1; i < ed.isize( ); i++ )
               {    if ( ed[i].second == radius ) continue;
                    int w = to_left[ ed[i].first ];
                    for ( int j = 0; j < To(w).isize( ); j++ )
                    {    int e = EdgeObjectIndexByIndexTo( w, j );
                         int d = ed[i].second + 1;
                         Bool found = False;
                         for ( int l = 0; l < ed.isize( ); l++ )
                         {    if ( ed[l].first == e )
                              {    found = True;
                                   if ( ed[l].second > d )
                                   {    ed[l].second = d;
                                        progress = True;    }    }    }
                         if ( !found )
                         {    ed.push( e, d );
                              progress = True;    }    }    }
               for ( int i = 0; i < ed.isize( ); i++ )
               {    if ( ed[i].second == radius ) continue;
                    int w = to_right[ ed[i].first ];
                    for ( int j = 0; j < From(w).isize( ); j++ )
                    {    int e = EdgeObjectIndexByIndexFrom( w, j );
                         int d = ed[i].second + 1;
                         Bool found = False;
                         for ( int l = 0; l < ed.isize( ); l++ )
                         {    if ( ed[l].first == e )
                              {    found = True;
                                   if ( ed[l].second > d )
                                   {    ed[l].second = d;
                                        progress = True;    }    }    }
                         if ( !found )
                         {    ed.push( e, d );
                              progress = True;    }    }    }
               if ( !progress ) break;    }

          // Define (e1,f1) and (e2,f2) pairs.

          vec< std::pair<int,int> > s1, s2;
          {    vec<int> x, s;
               for ( int l = 0; l < ed.isize( ); l++ )
                    x.push_back( ed[l].first );
               UniqueSort(x);
               s.push_back(e1);
               for ( int l = 0; l < radius; l++ )
               {    int n = s.size( );
                    for ( int m = 0; m < n; m++ )
                         s.append( FromEdgeObj( to_right[ s[m] ] ) );
                    UniqueSort(s);    }
               for ( int l = 0; l < s.isize( ); l++ )
                    s1.push( e1, s[l] );
               for ( int l = 0; l < x.isize( ); l++ )
               {    int e2 = x[l];
                    vec<int> t;
                    t.push_back(e2);
                    for ( int k = 0; k < radius; k++ )
                    {    int n = t.size( );
                         for ( int m = 0; m < n; m++ )
                              t.append( FromEdgeObj( to_right[ t[m] ] ) );
                         UniqueSort(t);    }
                    t = Intersection( t, x );
                    for ( auto m : t ) s2.push( e2, m );    }    }

          // Define e2, f1, f2.
          
          for ( int me = 0; me < s1.isize( ); me++ )
          for ( int mf = 0; mf < s2.isize( ); mf++ )
          {    int f1 = s1[me].second, e2 = s2[mf].first, f2 = s2[mf].second;
               vec<int> all;
               all.push_back( e1, e2, f1, f2 );
               UniqueSort(all);
               if ( all.size( ) != 4 ) continue;

               // Find edges reachable from e1 and e2, and bridge paths.

               vec<int> reach;
               Bool found_f1 = False, found_f2 = False;
               vec< vec<int> > bridge11, bridge22, bridge12, bridge21;
               Bool bad = False;
               for ( int z = 1; z <= 2; z++ )
               {    int e = ( z == 1 ? e1 : e2 );
                    for ( int i = 0; i < paths_index[e].isize( ); i++ )
                    {    int id = paths_index[e][i].first; 
                         int p = paths_index[e][i].second;
                         for ( int j = p + 1; j < Path(id).isize( ); j++ )
                         {    int g = Path(id)[j];
                              if ( e == e1 && g == e2 ) bad = True;
                              if ( e == e2 && g == e1 ) bad = True;
                              if ( e == e1 && g == f1 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge11.push_back(b);    }
                              if ( e == e2 && g == f2 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge22.push_back(b);    }
                              if ( e == e1 && g == f2 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge12.push_back(b);    }
                              if ( e == e2 && g == f1 )
                              {    vec<int> b;
                                   for ( int l = p; l <= j; l++ )
                                        b.push_back( Path(id)[l] );
                                   bridge21.push_back(b);    }
                              if ( g == f1 ) 
                              {    found_f1 = True;
                                   for ( int k = j + 1; k < Path(id).isize( ); k++ )
                                   {    if ( e == e1 && Path(id)[k] == e2 ) 
                                             bad = True;
                                        if ( e == e2 && Path(id)[k] == e1 ) 
                                             bad = True;    }
                                   break;    }
                              if ( g == f2 ) 
                              {    found_f2 = True;
                                   for ( int k = j + 1; k < Path(id).isize( ); k++ )
                                   {    if ( e == e1 && Path(id)[k] == e2 ) 
                                             bad = True;
                                        if ( e == e2 && Path(id)[k] == e1 ) 
                                             bad = True;    }
                                   break;    }
                              reach.push_back(g);    }    }    }
               if ( bad || !found_f1 || !found_f2 ) continue;
               UniqueSort(reach), UniqueSort(bridge11), UniqueSort(bridge22);
               if ( BinMember( reach, e1 ) || BinMember( reach, e2 ) ) continue;

               // Check for paths that are consistent with no bridge.

               vec< vec<int> > all_bridges;
               all_bridges.append(bridge11);
               all_bridges.append(bridge22);
               all_bridges.append(bridge12);
               all_bridges.append(bridge21);
               for ( int z = 1; z <= 2; z++ )
               {    int e = ( z == 1 ? e1 : e2 );
                    for ( int i = 0; i < paths_index[e].isize( ); i++ )
                    {    int id = paths_index[e][i].first; 
                         int p = paths_index[e][i].second;
                         Bool term = False;
                         for ( int j = p + 1; j < Path(id).isize( ); j++ )
                         {    if ( Path(id)[j] == f1 || Path(id)[j] == f2 )
                                   term = True;    }
                         if (term) continue;
                         vec<int> b;
                         for ( int j = p; j < Path(id).isize( ); j++ )
                              b.push_back( Path(id)[j] );
                         Bool normal = False;
                         for ( int i = 0; i < all_bridges.isize( ); i++ )
                              if ( all_bridges[i].Contains( b, 0 ) ) normal = True;
                         if ( !normal) bad = True;    }    }
               if (bad) continue;

               // See if {e1,e2}{reach}{f1,f2} makes sense.  Note that we're 
               // requiring all intermediate edges to be reachable by paths 
               // starting at e1 or e2, which is probably too stringent.

               for ( int i = 0; i < reach.isize( ); i++ )
               {    int x = to_left[ reach[i] ], y = to_right[ reach[i] ];
                    for ( int j = 0; j < To(x).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexTo( x, j );
                         if ( g != e1 && g != e2 && !BinMember( reach, g ) )
                              bad = True;    }
                    for ( int j = 0; j < From(y).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexFrom( y, j );
                         if ( g != f1 && g != f2 && !BinMember( reach, g ) )
                              bad = True;    }    }
               for ( int i = 0; i < 2; i++ )
               {    int x = ( i == 0 ? to_left[f1] : to_left[f2] );
                    for ( int j = 0; j < To(x).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexTo( x, j );
                         if ( g != e1 && g != e2 && !BinMember( reach, g ) )
                              bad = True;    }
                    int y = ( i == 0 ? to_right[e1] : to_right[e2] );
                    for ( int j = 0; j < From(y).isize( ); j++ )
                    {    int g = EdgeObjectIndexByIndexFrom( y, j );
                         if ( g != f1 && g != f2 && !BinMember( reach, g ) )
                              bad = True;    }    }
               if (bad) continue;

               // The reach edges must comprise a connected subgraph.

               if ( reach.nonempty( ) )
               {    vec<int> r;
                    r.push_back( reach[0] );
                    for ( int i = 0; i < r.isize( ); i++ )
                    {    int v = to_left[ r[i] ], w = to_right[ r[i] ];
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    int x = ( pass == 1 ? v : w );
                              for ( int j = 0; j < From(x).isize( ); j++ )
                              {    int e = EdgeObjectIndexByIndexFrom( x, j );
                                   if ( BinMember( reach, e ) && !Member( r, e ) )
                                        r.push_back(e);    }
                              for ( int j = 0; j < To(x).isize( ); j++ )
                              {    int e = EdgeObjectIndexByIndexTo( x, j );
                                   if ( BinMember( reach, e ) && !Member( r, e ) )
                                        r.push_back(e);    }    }    }
                    if ( r.size( ) < reach.size( ) ) continue;    }

               // Compute maximum bridge length.

               int B = 0;
               for ( int i = 0; i < bridge11.isize( ); i++ )
               {    int b = 0;
                    for ( int j = 1; j < bridge11[i].isize( ) - 1; j++ )
                         b += EdgeLengthKmers( bridge11[i][j] );
                    B = Max( B, b );    }
               for ( int i = 0; i < bridge22.isize( ); i++ )
               {    int b = 0;
                    for ( int j = 1; j < bridge22[i].isize( ) - 1; j++ )
                         b += EdgeLengthKmers( bridge22[i][j] );
                    B = Max( B, b );    }

               // See if cross-paths satisfy hypotheses.

               fix64_6 weight_11 = 0, weight_12 = 0, weight_21 = 0, weight_22 = 0;
               for ( int i = 0; i < paths_index[e1].isize( ); i++ )
               {    int id = paths_index[e1][i].first, p = paths_index[e1][i].second;
                    for ( int j = p + 1; j < Path(id).isize( ); j++ )
                    {    if ( Path(id)[j] == f1 ) weight_11 += Weight(id);
                         if ( Path(id)[j] == f2 ) weight_12 += Weight(id);    }    }
               for ( int i = 0; i < paths_index[e2].isize( ); i++ )
               {    int id = paths_index[e2][i].first, p = paths_index[e2][i].second;
                    for ( int j = p + 1; j < Path(id).isize( ); j++ )
                    {    if ( Path(id)[j] == f1 ) weight_21 += Weight(id);
                         if ( Path(id)[j] == f2 ) weight_22 += Weight(id);    }    }
               const int min_weight_split_low = 2;
               Bool OK = False;
               if ( weight_11 >= min_weight_split_low
                    && weight_22 >= min_weight_split_low
                    && weight_12 == 0 && weight_21 == 0 )
               {    OK = True;    } 
               if ( weight_11 >= min_weight_split && weight_22 >= min_weight_split
                    && weight_12 + weight_21 <= 2
                    && B <= MedianCorrectedReadLengthFudge( ) )
               {    OK = True;    }
               if ( weight_11 >= min_weight_split/2 
                    && weight_22 >= min_weight_split/2
                    && weight_12 + weight_21 < 2
                    && B <= MedianCorrectedReadLengthFudge( ) )
               {    OK = True;    }
               if ( ( OK && verbosity >= 2 ) || verbosity >= 3 )
               {    
                    #pragma omp critical
                    {    std::cout << "\n";
                         PRINT4( e1, e2, f1, f2 );
                         std::cout << "reach = " << printSeq(reach) << std::endl;
                         std::cout << "bridge11:\n";
                         for ( int i = 0; i < bridge11.isize( ); i++ )
                         {    std::cout << "[" << i+1 << "] " << printSeq( bridge11[i] ) 
                                   << "\n";    }
                         std::cout << "bridge22:\n";
                         for ( int i = 0; i < bridge22.isize( ); i++ )
                         {    std::cout << "[" << i+1 << "] " << printSeq( bridge22[i] ) 
                                   << "\n";    }
                         if ( verbosity >= 3 && !OK ) 
                         {    PRINT4( weight_11, weight_22, weight_12, weight_21 );
                              std::cout << "rejecting" << std::endl;    }    }    }
               if ( !OK ) continue;
     
               // Save join.

               vec<int> e12, f12;
               e12.push_back(e1,e2), f12.push_back(f1,f2);
               #pragma omp critical
               {    joins.push( e12, reach, f12 );
                    paths1.push_back(bridge11); 
                    paths2.push_back(bridge22);    }    }    }

     // Process joins.

     ParallelSortSync( joins, paths1, paths2 );
     if (logc.STATUS_LOGGING)
     {    std::cout << Date( ) << ": processing " << joins.size( ) 
               << " potential joins" << std::endl;    }
     if ( verbosity >= 2 ) std::cout << "\n";
     vec<Bool> touched( EdgeObjectCount( ), False );
     for ( int i = 0; i < joins.isize( ); i++ )
     {    Bool overlap = False;
          for ( int j = 0; j < joins[i].first.isize( ); j++ )
               if ( touched[ joins[i].first[j] ] ) overlap = True;
          for ( int j = 0; j < joins[i].second.isize( ); j++ )
               if ( touched[ joins[i].second[j] ] ) overlap = True;
          for ( int j = 0; j < joins[i].third.isize( ); j++ )
               if ( touched[ joins[i].third[j] ] ) overlap = True;
          if (overlap) continue;

          vec< triple< vec<int>, vec<int>, vec<int> > > proc;
          vec< vec< vec<int> > > proc1, proc2;
          proc.push_back( joins[i] );
          proc1.push_back( paths1[i] ), proc2.push_back( paths2[i] );

          if ( Inv( joins[i].first[0] ) >= 0 )
          {    vec<int> a, b, c;
               a.push_back( Inv( joins[i].third[0] ), Inv( joins[i].third[1] ) );
               for ( int j = 0; j < joins[i].second.isize( ); j++ )
                    b.push_back( Inv( joins[i].second[j] ) );
               Sort(b);
               c.push_back( Inv( joins[i].first[0] ), Inv( joins[i].first[1] ) );
               vec<int> all1, all2;
               all1.append( joins[i].first );
               all1.append( joins[i].second );
               all1.append( joins[i].third );
               all2.append(a), all2.append(b), all2.append(c);
               Sort(all1), Sort(all2);
               if ( Meet( all1, all2 ) ) continue;
               for ( int j = 0; j < all2.isize( ); j++ )
                    if ( touched[ all2[j] ] ) overlap = True;
               if (overlap) continue;
               proc.push( a, b, c );   
               vec< vec<int> > p1 = paths1[i], p2 = paths2[i];
               for ( int j = 0; j < p1.isize( ); j++ )
               {    p1[j].ReverseMe( );
                    for ( int l = 0; l < p1[j].isize( ); l++ )
                         p1[j][l] = Inv( p1[j][l] );   }
               for ( int j = 0; j < p2.isize( ); j++ )
               {    p2[j].ReverseMe( );
                    for ( int l = 0; l < p2[j].isize( ); l++ )
                         p2[j][l] = Inv( p2[j][l] );   }
               proc1.push_back(p1), proc2.push_back(p2);    }

          vec<fix64_6> weight( EdgeObjectCount( ), 0 );
          for ( int i = 0; i < NPaths( ); i++ )
          for ( int j = 0; j < Path(i).isize( ); j++ )
               weight[ Path(i)[j] ] += Weight(i);
          vec<int> count(2, 0);
          for ( int p = 0; p < proc.isize( ); p++ )
          {    if ( verbosity >= 2 )
               {    int e1 = proc[p].first[0], e2 = proc[p].first[1];
                    int f1 = proc[p].third[0], f2 = proc[p].third[1];
                    std::cout << "joining: ";
                    PRINT4(e1,e2,f1,f2);    }
               for ( int j = 0; j < proc[p].first.isize( ); j++ )
                    touched[ proc[p].first[j] ] = True;
               for ( int j = 0; j < proc[p].second.isize( ); j++ )
                    touched[ proc[p].second[j] ] = True;
               for ( int j = 0; j < proc[p].third.isize( ); j++ )
                    touched[ proc[p].third[j] ] = True;
               vec<int> dels;
               dels.append( proc[p].first );
               dels.append( proc[p].second );
               dels.append( proc[p].third );
               const int max_del_weight = 4;
               Bool bad = False;
               for ( int i = 0; i < dels.isize( ); i++ )
               {    Bool used = False;
                    for ( int j = 0; j < proc1[p].isize( ); j++ )
                         if ( Member( proc1[p][j], dels[i] ) ) used = True;
                    for ( int j = 0; j < proc2[p].isize( ); j++ )
                         if ( Member( proc2[p][j], dels[i] ) ) used = True;
                    if (used) continue;
                    if ( weight[ dels[i] ] > max_del_weight ) bad = True;    }
               if (bad)
               {    if ( verbosity >= 2 ) std::cout << "aborting join" << std::endl;
                    continue;    }
               progress = True;
               int v1 = to_left[ proc[p].first[0] ]; 
               int v2 = to_left[ proc[p].first[1] ];
               int w1 = to_right[ proc[p].third[0] ]; 
               int w2 = to_right[ proc[p].third[1] ];
               vec< vec<int> > e( proc1[p] );
               e.append( proc2[p] );
               vec<int> f;
               for ( int j = 0; j < proc1[p].isize( ); j++ )
               {    f.push_back( EdgeObjectCount( ) );
                    int rid = -1;
                    if ( p == 1 ) rid = EdgeObjectCount( ) - count[0];
                    InvMutable( ).push_back(rid);
                    if ( p == 1 )
                    {    InvMutable( EdgeObjectCount( ) - count[0] ) 
                              = EdgeObjectCount( );    }
                    AddEdge( v1, w1, Cat( proc1[p][j] ) );
                    count[p]++;    }
               for ( int j = 0; j < proc2[p].isize( ); j++ )
               {    f.push_back( EdgeObjectCount( ) );
                    int rid = -1;
                    if ( p == 1 ) rid = EdgeObjectCount( ) - count[0];
                    InvMutable( ).push_back(rid);
                    if ( p == 1 )
                    {    InvMutable( EdgeObjectCount( ) - count[0] ) 
                              = EdgeObjectCount( );    }
                    AddEdge( v2, w2, Cat( proc2[p][j] ) );
                    count[p]++;    }
               DeleteEdges(dels);
               TransformPaths( e, f );    }    }
     if ( verbosity >= 2 ) std::cout << "\n";

     // Clean up.

     UniqueOrderPaths( );
     RemoveEdgelessVertices( );
     REPORT_TIME( clock, "used in PullApart2" );
     RemoveDeadEdgeObjects( );
     TestValid(logc);
     DeleteReverseComplementComponents(logc);    

     if ( !progress ) break;    }    }
#endif

// Note that we're not yet testing pairs for symmetry, and it's not 100% clear
// what the tests would be.

void SupportedHyperBasevectorPrivate::TestValid( const long_logging& logc,
     const Bool test_paths ) const
{    double clock = WallClockTime( );

     // Test for empty.

     if ( EdgeObjectCount( ) == 0 ) DiscovarTools::ExitAssemblyEmpty( );

     // Test HyperBasevector structure.

     HyperBasevector::TestValid( );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);

     // Test involution.

     if ( inv_.isize( ) != EdgeObjectCount( ) )
     {    std::cout << "\n";
          PRINT2( EdgeObjectCount( ), inv_.size( ) );
          std::cout << "SupportedHyperBasevector is invalid.\n"
               << "Involution has wrong size.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }
     for ( int e = 0; e < EdgeObjectCount( ); e++ )
     {    if ( !InvDef(e) ) continue;
          if ( Inv(e) < 0 || Inv(e) >= EdgeObjectCount( ) )
          {    std::cout << "\n";
               PRINT3( e, Inv(e), EdgeObjectCount( ) );
               std::cout << "SupportedHyperBasevector is invalid.\n"
                    << "Illegal involution value.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }
          basevector b = EdgeObject(e);
          b.ReverseComplement( );
          if ( b != EdgeObject( Inv(e) ) )
          {    std::cout << "\n";
               int re = Inv(e);
               PRINT4( e, re, b.size( ), EdgeObject(re).size( ) );
               std::cout << "SupportedHyperBasevector is invalid.\n"
                    << "Involution value not rc.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }
          if ( Inv(Inv(e)) != e )
          {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Involution is not an involution.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i1 = 0; i1 < To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < From(v).isize( ); i2++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = EdgeObjectIndexByIndexFrom( v, i2 );
               if ( InvDef(e1) != InvDef(e2) )
               {    std::cout << "\n";
                    PRINT3( e1, EdgeLengthKmers(e1), Inv(e1) );
                    PRINT3( e2, EdgeLengthKmers(e2), Inv(e2) );
                    std::cout << "SupportedHyperBasevector is invalid.\n"
                         << "Some edges in the same component have the inversion "
                         << "defined and some do not.\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i1 = 0; i1 < To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < From(v).isize( ); i2++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = EdgeObjectIndexByIndexFrom( v, i2 );
               int re1 = Inv(e1), re2 = Inv(e2);
               if ( InvDef(e1) && InvDef(e2) && to_right[re2] != to_left[re1] )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Involution does not preserve graph structure.\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }
     // Note sure the following can happen:
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < From(v).isize( ); i++ )
          {    int w = From(v)[i];
               int e = EdgeObjectIndexByIndexFrom( v, i );
               if ( !InvDef(e) ) continue;
               int re = Inv(e);
               int rv = to_left[re], rw = to_right[re];
               if ( From(rv).size( ) != To(w).size( )
                    || To(rw).size( ) != From(v).size( ) )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Graph structure is asymmetric.\n" << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }    

     // Return if we're not testing paths.

     if ( !test_paths )
     {    REPORT_TIME( clock, "used testing valid" );
          return;    }

     for ( int j = 0; j < NPaths( ); j++ )
     {    if ( Path(j).empty( ) )
          {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Path is empty.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }
          for ( int i = 0; i < Path(j).isize( ); i++ )
          {    if ( Path(j,i) < 0 )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal negative path value.\n" << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < Path(j).isize( ) - 1; i++ )
          {    if ( Path(j,i) < 0 || Path(j,i+1) < 0 ) continue;
               if ( to_right[ Path(j,i) ] != to_left[ Path(j,i+1) ] )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in path " << j << ": " << Path(j,i) 
                         << ", " << Path(j,i+1) << "\n" << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }    

     if ( Paths( ).size( ) != WeightsFw( ).size( ) )
     {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Paths and weights have different sizes.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }
     if ( WeightsFw( ).size( ) != WeightsRc( ).size( ) )
     {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Fw and rc weights have different sizes.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }

     // Test for no support.

     if ( NPaths( ) == 0 )
     {    std::cout << "\nSupportedHyperBasevector is invalid.\n";
          std::cout << "There are no paths.\n" << "Abort." << std::endl;
          DiscovarTools::ExitPathsEmpty( ); }
//          TracebackThisProcess( );    }

     // Test paths for unique ordering.

     if ( !Paths( ).UniqueOrdered( ) )
     {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Paths are not unique-ordered.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }

     // Test paths and weights for symmetry.

     for ( int i1 = 0; i1 < NPaths( ); i1++ )
     {    const vec<int>& p1 = Path(i1);
          if ( p1[0] < 0 || !InvDef( p1[0] ) ) continue;
          vec<int> p2;
          for ( int j = 0; j < p1.isize( ); j++ )
               p2.push_back( Inv( p1[j] ) );
          p2.ReverseMe( );
          int i2 = BinPosition( Paths( ), p2 );
          if ( i2 < 0 )
          {    std::cout << "\nPath(" << i1 << ") = " << printSeq( Path(i1) ) << "\n"
                    << "SupportedHyperBasevector is invalid.\n"
                    << "The reverse complement of path " << i1 << " is missing.\n" 
                    << "Abort." << std::endl;
                    TracebackThisProcess( );    }
          if ( WeightFw(i1) != WeightRc(i2) || WeightRc(i1) != WeightFw(i2) )
          {    std::cout << "\n";
               PRINT2( WeightFw(i1), WeightRc(i2) );
               PRINT2( WeightRc(i1), WeightFw(i2) );
               std::cout << "SupportedHyperBasevector is invalid.\n"
                    << "Weights are asymmetric.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }    }

     // Test pairs.

     if ( Pairs( ).size( ) != PairData( ).size( ) )
     {    std::cout << "SupportedHyperBasevector is invalid.\n"
               << "Pairs and PairData have different sizes.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }
     if ( !Pairs( ).UniqueOrdered( ) )
     {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
               << "Pairs are not unique-ordered.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }
     for ( int i = 0; i < NPairs( ); i++ )
     {    if ( PairLeft(i).empty( ) || PairRight(i).empty( ) )
          {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "Pair side is empty.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }    }
     for ( int j = 0; j < NPairs( ); j++ )
     {    for ( int i = 0; i < PairLeft(j).isize( ); i++ )
          {    if ( PairLeft(j,i) < 0 || PairLeft(j,i) >= EdgeObjectCount( ) )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal edge id in left side of pair.\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairRight(j).isize( ); i++ )
          {    if ( PairRight(j,i) < 0 || PairRight(j,i) >= EdgeObjectCount( ) )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Illegal edge id in right side of pair.\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairLeft(j).isize( ) - 1; i++ )
          {    if ( PairLeft(j,i) < 0 || PairLeft(j,i+1) < 0 ) continue;
               if ( to_right[ PairLeft(j,i) ] != to_left[ PairLeft(j,i+1) ] )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in left side of pair " << j << ": " 
                         << PairLeft(j,i) << ", " << PairLeft(j,i+1) << "\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }
          for ( int i = 0; i < PairRight(j).isize( ) - 1; i++ )
          {    if ( PairRight(j,i) < 0 || PairRight(j,i+1) < 0 ) continue;
               if ( to_right[ PairRight(j,i) ] != to_left[ PairRight(j,i+1) ] )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Non-adjacent edges in right side of pair " << j << ": " 
                         << PairRight(j,i) << ", " << PairRight(j,i+1) << "\n" 
                         << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }
          if ( PairData(j).empty( ) )
          {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                    << "No pair data." << std::endl << "Abort." << std::endl;
               TracebackThisProcess( );    }
          for ( int i = 0; i < PairData(j).isize( ); i++ )
          {    const pair_point& p = PairData(j,i);
               if ( p.Trim( ) < -1000000000 || p.Trim( ) > 1000000000 )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Insane value " << p.Trim( ) << " for trim." 
                         << std::endl << "Abort." << std::endl;
                    TracebackThisProcess( );    }
               if ( p.Lib( ) < 0 )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Negative library id." << std::endl << "Abort." << std::endl;
                    TracebackThisProcess( );    }
               if ( p.Weight( ) < 0 )
               {    std::cout << "\nSupportedHyperBasevector is invalid.\n"
                         << "Pair has negative weight." << std::endl << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }

     REPORT_TIME( clock, "used testing valid" );    }

void SupportedHyperBasevector::DeleteReverseComplementComponents( const long_logging& logc , const int64_t iDirSortFactor)
{    double clock = WallClockTime( );
     vec< vec<int> > components;
     ComponentEdges(components);
     
     if(iDirSortFactor==0){
         Sort(components);
     }
     else{
         vec<int> to_left, to_right;
         ToLeft(to_left), ToRight(to_right);
         vec<std::tuple<int64_t,int64_t,vec<int>>> order;
         order.reserve(components.size());
         
         for(const auto& entry: components){
             int64_t diff=0,bdiff=0;
             for(const auto& edge: entry){
                 if( Source(to_left[edge])){
                     bdiff-=iDirSortFactor*EdgeObject(edge).size();
                     diff-=iDirSortFactor;
                 }
                 if( Sink(to_right[edge])){
                     bdiff+=iDirSortFactor*EdgeObject(edge).size();
                     diff+=iDirSortFactor;
                 }
             }
             order.emplace_back(diff,bdiff,entry);
         }
         Sort(order);
         
         components.clear();
         for(const auto&entry: order){ components.push_back(std::get<2>(entry)); }
     }
     
     vec<int> rc_to_delete;
     for ( size_t i = 0; i < components.size(); i++ )
     {    vec<int> rc;
          for ( size_t j = 0; j < components[i].size(); j++ )
               rc.push_back( Inv( components[i][j] ) );
          Sort(rc);
          int p = (iDirSortFactor==0)?BinPosition( components, rc ):Position( components, rc );
          if ( p > (int)i )
          {    for ( size_t j = 0; j < components[p].size(); j++ )
                    rc_to_delete.push_back( components[p][j] );
               for ( size_t j = 0; j < components[i].size(); j++ )
                    InvMutable( components[i][j] ) = -1;    }    }
     Sort(rc_to_delete);
     if ( logc.verb[ "DELETE_RC" ] >= 1 )
     {    std::cout << Date( ) << ": rc component deletion -- deleting " 
               << rc_to_delete.size( ) << " edges" << std::endl;    }
     DeleteEdges(rc_to_delete);
     DeleteUnusedPaths( );
     REPORT_TIME( clock, "used deleting reverse complement components" );
     double eclock = WallClockTime( );
     RemoveEdgelessVertices( );
     RemoveUnneededVertices( );
     REPORT_TIME( eclock, "used cleaning up after deleting rc components" );
     RemoveDeadEdgeObjects( );
     TestValid(logc);    }

void SupportedHyperBasevector::Reverse( )
{    HyperBasevector::Reverse( );
     for ( int i = 0; i < NPaths( ); i++ )
          PathMutable(i).ReverseMe( );
     for ( int i = 0; i < NPairs( ); i++ )
     {    PairLeftMutable(i).ReverseMe( );
          PairRightMutable(i).ReverseMe( );
          PairMutable(i) = std::make_pair( PairRight(i), PairLeft(i) );    }    }

