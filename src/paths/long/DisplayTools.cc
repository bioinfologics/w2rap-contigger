///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "ParseSet.h"
#include "paths/HyperBasevector.h"
#include "paths/long/DisplayTools.h"
#include "random/Random.h"

void ParseSeeds( const HyperBasevector& hb, const vec<int>& to_right,
     const String& SEEDS, const int RANDOM_SEED, const String& SEEDS_MINUS,
     vec<int>& seeds )
{
     if ( RANDOM_SEED >= 0 ) srandomx(RANDOM_SEED);
     seeds.clear( );
     vec<int> seeds_minus;
     ParseIntSet( SEEDS_MINUS, seeds_minus );
     if ( SEEDS == "all" )
     {    for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
               seeds.push_back(e);    }
     else 
     {    vec<String> seedsx;
          ParseStringSet( SEEDS, seedsx );
          for ( int i = 0; i < seedsx.isize( ); i++ )
          {    String s = seedsx[i];

               if ( s.Contains( "random:", 0 ) )
               {    int l = 0;
                    s = s.After( "random:" );
                    if ( s.Contains( ":" ) )
                    {    l = s.After( ":" ).Int( );
                         s = s.Before( ":" );    }
                    int n = s.Int( );
                    for ( int j = 0; j < n; j++ )
                    {    int e = randomx( ) % hb.EdgeObjectCount( );
                         if ( Member( seeds, e ) || hb.EdgeLengthKmers(e) < l )
                         {    j--;
                              continue;    }
                         seeds.push_back(e);    }    }

               else if ( s.Contains( "trandom:", 0 ) )
               {    uint seed;
                    timeval t;
                    gettimeofday( &t, NULL );
                    seed = t.tv_sec + t.tv_usec;
                    srandomx(seed);
                    int l = 0;
                    s = s.After( "trandom:" );
                    if ( s.Contains( ":" ) )
                    {    l = s.After( ":" ).Int( );
                         s = s.Before( ":" );    }
                    int n = s.Int( );
                    for ( int j = 0; j < n; j++ )
                    {    int e = randomx( ) % hb.EdgeObjectCount( );
                         if ( Member( seeds, e ) || hb.EdgeLengthKmers(e) < l )
                         {    j--;
                              continue;    }
                         seeds.push_back(e);    }    }

               else if ( s.IsInt( ) ) seeds.push_back( s.Int( ) );
               else
               {    if ( !s.Contains( ".." ) )
                    {    std::cout << "Illegal seed notation." << std::endl;
                         Scram(1);    }
                    int s1 = s.Before( ".." ).Int( ), s2 = s.After( ".." ).Int( );
                    vec<int> x;
                    x.push_back(s1);
                    for ( int j = 0; j < x.isize( ); j++ )
                    {    int v = to_right[ x[j] ];
                         for ( int l = 0; l < hb.From(v).isize( ); l++ )
                         {    int e = hb.EdgeObjectIndexByIndexFrom( v, l );
                              if ( e != s2 && !Member( x, e ) ) 
                                   x.push_back(e);    }    }
                    x.push_back(s2);
                    UniqueSort(x);
                    seeds.append(x);    }    }    }
     UniqueSort(seeds), UniqueSort(seeds_minus);
     vec<Bool> sdel( seeds.size( ) );
     for ( int i = 0; i < seeds.isize( ); i++ )
          if ( BinMember( seeds_minus, seeds[i] ) ) sdel[i] = True;
     EraseIf( seeds, sdel );
     for ( int i = 0; i < seeds.isize( ); i++ )
     {    if ( seeds[i] < 0 || seeds[i] >= hb.EdgeObjectCount( ) )
          {    std::cout << seeds[i] << " is an illegal seed value" << std::endl;
               Scram(1);    }    }    }

Bool ParseSeeds( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec< triple<kmer<20>,int,int> >& kmers_plus,
     const vec<vec<vec<vec<int>>>>& lines, const vec<String>& genome_names,
     const vec< std::pair<int,ho_interval> >& ambint, Bool& ambflag,
     const vec< vec< std::pair<int,int> > >& hits, const String& SEEDS,
     const int RANDOM_SEED, const String& SEEDS_MINUS, vec<int>& seeds,
     const int max_seeds, std::ostream& tout )
{
     if ( RANDOM_SEED >= 0 ) srandomx(RANDOM_SEED);
     seeds.clear( );
     vec<int> seeds_minus;
     ParseIntSet( SEEDS_MINUS, seeds_minus );
     if ( SEEDS == "all" )
     {    for ( int e = 0; e < hb.E( ); e++ )
               seeds.push_back(e);    }
     else 
     {    vec<String> seedsx;
          ParseStringSet( SEEDS, seedsx );
          for ( int i = 0; i < seedsx.isize( ); i++ )
          {    String s = seedsx[i];

               // Handle sequence case.

               vec<char> content, ACGT = {'A', 'C', 'G', 'T'};
               for ( int j = 0; j < s.isize( ); j++ )
                    content.push_back( s[j] );
               UniqueSort(content);
               if ( s.size( ) > 0 && BinSubset( content, ACGT ) )
               {    basevector b(s);
                    const int L = 20;
                    tout << "query has length " << b.size( ) << std::endl;
                    if ( b.isize( ) < L )
                    {    tout << "Sorry, need at least " << L << " bases." << std::endl;
                         return False;    }
                    basevector rb(b);
                    rb.ReverseComplement( );
                    vec<int> all;
                    if ( kmers_plus.size( ) > 0 )
                    {    kmer<L> x;
                         x.SetToSubOf( b, 0 );
                         int64_t low = LowerBound1( kmers_plus, x );
                         int64_t high = UpperBound1( kmers_plus, x );
                         for ( int64_t j = low; j < high; j++ )
                         {    int e = kmers_plus[j].second;
                              int epos = kmers_plus[j].third;
                              Bool mismatch = False;
                              for ( int l = L; l < b.isize( ); l++ )
                              {    if ( epos+l >= hb.EdgeObject(e).isize( )
                                        || b[l] != hb.EdgeObject(e)[epos+l] )
                                   {    mismatch = True;
                                        break;    }    }
                              if ( !mismatch) all.push_back(e);    }
                         x.ReverseComplement( );
                         low = LowerBound1( kmers_plus, x );
                         high = UpperBound1( kmers_plus, x );
                         for ( int64_t j = low; j < high; j++ )
                         {    int re = kmers_plus[j].second;
                              int repos = kmers_plus[j].third;
                              Bool mismatch = False;
                              int n = rb.size( );
                              for ( int l = 0; l < n-L; l++ )
                              {    if ( repos - (n-L) + l < 0
                                        || rb[l] != hb.EdgeObject(re)
                                           [ repos - (n-L) + l ] )
                                   {    mismatch = True;
                                        break;    }    }
                              if ( !mismatch) all.push_back( inv[re] );    }    }
                    else
                    {
                         #pragma omp parallel for
                         for ( int e = 0; e < hb.E( ); e++ )
                         {    String t = hb.EdgeObject(e).ToString( );
                              if ( t.Contains(s) )

                                   #pragma omp critical
                                   {    seeds.push_back(e);    
                                        all.push_back(e);    }    }    }
                    UniqueSort(all);
                    seeds.append(all);
                    tout << "using seeds";
                    if ( all.size( ) <= 40 ) tout << " " << printSeq(all) << std::endl;
                    else
                    {    for ( int i = 0; i < 20; i++ )
                              tout << " " << all[i];
                         tout << "... [" << all.size( ) << " in total]"
                              << std::endl;    }    }

               // Handle random: cases.

               else if ( s.Contains( "random:", 0 ) )
               {    int l = 0;
                    s = s.After( "random:" );
                    if ( s.Contains( ":" ) )
                    {    if ( !s.After( ":" ).IsInt( ) )
                         {    tout << "Illegal seed notation " << s << "." << std::endl;
                              return False;    }
                         l = s.After( ":" ).Int( );
                         s = s.Before( ":" );    }
                    int n = s.Int( );
                    for ( int j = 0; j < n; j++ )
                    {    int e = randomx( ) % hb.E( );
                         if ( Member( seeds, e ) || hb.Kmers(e) < l )
                         {    j--;
                              continue;    }
                         seeds.push_back(e);    }    }

               // Handle trandom: cases.

               else if ( s.Contains( "trandom:", 0 ) )
               {    uint seed;
                    timeval t;
                    gettimeofday( &t, NULL );
                    seed = t.tv_sec + t.tv_usec;
                    srandomx(seed);
                    int l = 0;
                    s = s.After( "trandom:" );
                    if ( s.Contains( ":" ) )
                    {    if ( !s.After( ":" ).IsInt( ) )
                         {    tout << "Illegal seed notation " << s << "." << std::endl;
                              return False;    }
                         l = s.After( ":" ).Int( );
                         s = s.Before( ":" );    }
                    int n = s.Int( );
                    for ( int j = 0; j < n; j++ )
                    {    int e = randomx( ) % hb.E( );
                         if ( Member( seeds, e ) || hb.Kmers(e) < l )
                         {    j--;
                              continue;    }
                         tout << "using seed " << e << std::endl;
                         seeds.push_back(e);    }    }
               
               // Handle simple edge case.

               else if ( s.IsInt( ) ) seeds.push_back( s.Int( ) );

               // Handle line case.

               else if ( s.Contains( "L", 0 ) && s.After( "L" ).IsInt( ) )
               {    int l = s.After( "L" ).Int( );
                    if ( l < 0 || l >= lines.isize( ) )
                    {    tout << "Illegal seed notation " << s << "." << std::endl;
                         return False;    }
                    int s1 = lines[l].front( )[0][0], s2 = lines[l].back( )[0][0];
                    if ( s1 == s2 ) seeds.push_back(s1);
                    else
                    {    vec<int> x;
                         x.push_back(s1);
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    int v = hb.ToRight( x[j] );
                              for ( int l = 0; l < (int) hb.From(v).size( ); l++ )
                              {    int e = hb.IFrom( v, l );
                                   if ( e != s2 && !Member( x, e ) ) 
                                        x.push_back(e);    }    }
                         x.push_back(s2);
                         UniqueSort(x);
                         seeds.append(x);    }    }
     
               // Handle edge range case.

               else if ( s.Contains( ".." ) )
               {    if ( !s.Before( ".." ).IsInt( ) || !s.After( ".." ).IsInt( ) )
                    {    tout << "Illegal seed notation " << s << "." << std::endl;
                         return False;    }
                    int s1 = s.Before( ".." ).Int( ), s2 = s.After( ".." ).Int( );
                    if ( s1 == s2 ) seeds.push_back(s1);
                    else
                    {    vec<int> x;
                         x.push_back(s1);
                         for ( int j = 0; j < x.isize( ); j++ )
                         {    if ( x.isize( ) > max_seeds ) break;
                              int v = hb.ToRight( x[j] );
                              for ( int l = 0; l < (int) hb.From(v).size( ); l++ )
                              {    int e = hb.IFrom( v, l );
                                   if ( e != s2 && !Member( x, e ) ) 
                                        x.push_back(e);    }    }
                         x.push_back(s2);
                         UniqueSort(x);
                         seeds.append(x);    }    }
     
               // Handle alignment case.

               else
               {    if ( genome_names.size( ) == 0 )
                    {    tout << "Can't parse " << s << " without having "
                              << "a genome.names file. Make sure you ran "
			      << "DiscovarDeNovo with the REFHEAD option."
			      << std::endl;
                         return False;    }
                    if ( s.Contains( ":" ) && !s.Contains( "-" )
                         && s.After( ":" ).IsInt( ) )
                    {    s = s + "-" + ToString( s.After( ":" ).Int( ) + 1 );    }
                    if ( !s.Contains( ":" ) || !s.After( ":" ).Contains( "-" ) )
                    {    tout << "Illegal seed notation " << s << "." << std::endl;
                         return False;    }
                    String gs = s.Before( ":" );
                    int g = Position( genome_names, gs );
                    if ( g < 0 ) 
                    {    tout << "Unknown genome record name." << std::endl;
                         return False;    }
                    if ( !s.Between( ":", "-" ).IsInt( ) 
                         || !s.After( "-" ).IsInt( ) )
                    {    tout << "Illegal seed notation " << s << "." << std::endl;
                         return False;    }
                    int start = s.Between( ":", "-" ).Int( );
                    int stop = s.After( "-" ).Int( );
                    if ( start >= stop )
                    {    tout << "Illegal seed notation " << s << "." << std::endl;
                         return False;    }
                    for ( int j = 0; j < ambint.isize( ); j++ )
                    {    if ( g == ambint[j].first && Subset( 
                              ho_interval(start,stop), ambint[j].second ) )
                         {    tout << "You've specified a reference locus "
                                   << "consisting "
                                   << "entirely of ambiguous bases." << std::endl;
                              ambflag = True;
                              break;    }    }
                    for ( int e = 0; e < hits.isize( ); e++ )
                    for ( int m = 0; m < hits[e].isize( ); m++ )
                    {    if ( hits[e][m].first != g ) continue;
                         {    int xstart = hits[e][m].second;
                              int xstop = xstart + hb.Bases(e);
                              if ( IntervalOverlap(start, stop, xstart, xstop) > 0 )
                                   seeds.push_back(e);    }    }    }    }    }

     UniqueSort(seeds), UniqueSort(seeds_minus);
     vec<Bool> sdel( seeds.size( ) );
     for ( int i = 0; i < seeds.isize( ); i++ )
          if ( BinMember( seeds_minus, seeds[i] ) ) sdel[i] = True;
     EraseIf( seeds, sdel );
     for ( int i = 0; i < seeds.isize( ); i++ )
     {    if ( seeds[i] < 0 || seeds[i] >= hb.E( ) )
          {    tout << seeds[i] << " is an illegal seed value" << std::endl;
               return False;    }    }
     return True;    }

void BuildNhood( const HyperBasevector& hb, const vec<int>& to_left,
     const vec<int>& to_right, const vec<int>& seeds, const int DEPTH, 
     vec<Bool>& invisible )
{
     invisible.resize( hb.EdgeObjectCount( ), True );
     for ( int i = 0; i < seeds.isize( ); i++ )
          invisible[ seeds[i] ] = False;
     // std::cout << Date( ) << ": building nhood" << std::endl;
     for ( int d = 0; d < DEPTH; d++ )
     {    vec<Bool> inv(invisible);
          for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          {    if ( inv[e] ) continue;
               int v = to_left[e], w = to_right[e];
               for ( int u = 0; u < 2; u++ )
               {    int x = ( u == 0 ? v : w );
                    for ( int j = 0; j < hb.To(x).isize( );  j++ )
                         invisible[ hb.EdgeObjectIndexByIndexTo( x, j ) ] = False;
                    for ( int j = 0; j < hb.From(x).isize( );  j++ )
                    {    invisible[ hb.EdgeObjectIndexByIndexFrom( x, j ) ]
                              = False;    }    }    }    }
     int vis = 0;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          if ( !invisible[e] ) vis++;
     std::cout << "found " << vis << " edges" << std::endl;    }

void BuildNhood( const HyperBasevectorX& hb, const vec<int>& seeds, 
     const int DEPTH, vec<Bool>& invisible, const int max_edges )
{    invisible.resize( hb.E( ), True );
     int nedges = seeds.size( );
     for ( int i = 0; i < seeds.isize( ); i++ )
          invisible[ seeds[i] ] = False;
     for ( int d = 0; d < DEPTH; d++ )
     {    if ( nedges > max_edges ) break;
          vec<Bool> inv(invisible);
          for ( int e = 0; e < hb.E( ); e++ )
          {    if ( nedges > max_edges ) break;
               if ( inv[e] ) continue;
               int v = hb.ToLeft(e), w = hb.ToRight(e);
               for ( int u = 0; u < 2; u++ )
               {    int x = ( u == 0 ? v : w );
                    for ( int j = 0; j < (int) hb.To(x).size( );  j++ )
                    {    if ( invisible[ hb.ITo( x, j ) ] )
                         {    invisible[ hb.ITo( x, j ) ] = False;
                              nedges++;    }    }
                    for ( int j = 0; j < (int) hb.From(x).size( );  j++ )
                    {    if ( invisible[ hb.IFrom( x, j ) ] )
                         {    invisible[ hb.IFrom( x, j ) ] = False;    
                              nedges++;    }    }    }    }    }    }
