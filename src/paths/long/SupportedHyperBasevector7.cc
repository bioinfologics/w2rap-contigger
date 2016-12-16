///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "fastg/FastgGraph.h"
#include "paths/HyperBasevector.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/LargeKDispatcher.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/SupportedHyperBasevector.h"

namespace { // open anonymous namespace

template<int K2> void MapAssembly( const HyperBasevector& hb, 
     const HyperBasevector& hb_orig, vec< vec<int> >& content, const int verbosity )
{    if ( verbosity >= 1 ) std::cout << Date( ) << ": mapping hb onto hb_orig" << std::endl;
     vecbasevector B;
     for ( int e = 0; e < hb_orig.EdgeObjectCount( ); e++ )
          B.push_back_reserve( hb_orig.EdgeObject(e) );
     vec< triple<kmer<K2>,int,int> > kmers_plus;
     MakeKmerLookup1( B, kmers_plus );
     vec< kmer<K2> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     content.resize( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    const basevector& b = hb.EdgeObject(e);
          kmer<K2> x;
          int pos = 0;
          while( pos < b.isize( ) - (K2-1) )
          {    x.SetToSubOf( b, pos );
               int64_t low = LowerBound( kmers, x );
               int64_t high = UpperBound( kmers, x );

               // A kmer in hb cannot occur more than once in hb_orig.  However,
               // note that it is possible (but rare, and perhaps never occurred) 
               // to have a kmer in hb that is not present in hb_orig.  The reason 
               // for this is that at some point two edges in hb_orig might have 
               // been joined along a (K-1)-base overlap, resulting in the creation 
               // of new kmers.  It is not entirely clear what the implications of 
               // this are here.

               if ( high == low )
               {    pos++;
                    continue;    }
               ForceAssertEq( high, low + 1 );

               int id = kmers_plus[low].second;
               content[e].push_back(id);
               if ( verbosity >= 1 ) PRINT2( e, id );

               // There used to be an assert here, but it doesn't work when there
               // are instances of high = low, as discussed above.

               // ForceAssertEq( kmers_plus[low].third, 0 );

               pos += B[id].isize( ) - (K2-1) - kmers_plus[low].third;    }    }    }

template <int K>
struct MapAsmFunctor
{
    void operator()( const HyperBasevector& hb,  const HyperBasevector& hb_orig,
                        vec< vec<int> >& content, const int verbosity )
    { MapAssembly<K>(hb,hb_orig,content,verbosity); }
};

void MapOne( const vec<int>& u, const vec< vec< std::pair<int,int> > >& cindex,
     const vec< vec<int> >& content, const HyperBasevector& hb,
     const vec<int>& to_right, vec< vec<int> >& fulls, vec<int>& fulls_pos )
{
     fulls.clear( );
     fulls_pos.clear( );
     vec< triple< vec<int>, int, int > > partials;
     const vec< std::pair<int,int> >& locs = cindex[ u[0] ];
     for ( int j = 0; j < locs.isize( ); j++ )
     {    int x = locs[j].first, y = locs[j].second;
          Bool mismatch = False;
          int l;
          for ( l = 0; l < u.isize( ); l++ )
          {    if ( y+l == content[x].isize( ) ) break;
               if ( content[x][y+l] != u[l] ) 
               {    mismatch = True;
                    break;    }    }
          if (mismatch) continue;
          vec<int> p(1);
          p[0] = x;
          partials.push( p, l, y );    }
     while( partials.nonempty( ) )
     {    vec<int> p = partials.back( ).first;
          int pos = partials.back( ).second;
          int y = partials.back( ).third;
          partials.pop_back( );
          if ( pos == u.isize( ) )
          {    fulls.push_back(p);
               fulls_pos.push_back(y);
               continue;    }
          int v = to_right[ p.back( ) ];
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int l, e = hb.EdgeObjectIndexByIndexFrom( v, j ); 
               Bool mismatch = False;
               for ( l = pos; l < u.isize( ); l++ )
               {    if ( l-pos == content[e].isize( ) ) break;
                    if ( content[e][l-pos] != u[l] ) 
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               vec<int> q(p);
               q.push_back(e);
               partials.push( q, l, y );    }    }    }

} // close anonymous namespace

// TransformPaths: replace each occurrence of x[i] by y[i].

void SupportedHyperBasevector::TransformPaths( const vec< vec<int> >& x, 
     const vec<int>& y, vec< vec< std::pair<int,int> > >& paths_index )
{    if ( paths_index.nonempty( ) ) paths_index.resize( EdgeObjectCount( ) );
     vec<int> all;
     for ( int j = 0; j < x.isize( ); j++ )
          all.append( x[j] );
     UniqueSort(all);

     vec<int> used_paths;
     if ( paths_index.nonempty( ) )
     {    for ( int i = 0; i < all.isize( ); i++ )
          {    int e = all[i];
               for ( int j = 0; j < paths_index[e].isize( ); j++ )
                    used_paths.push_back( paths_index[e][j].first );    }
          UniqueSort(used_paths);    }
     int to_use = ( paths_index.empty( ) ? NPaths( ) : used_paths.isize( ) );

     const int infty = 1000000000;
     #pragma omp parallel for
     for ( int pi = 0; pi < to_use; pi++ )
     {    int i = ( paths_index.empty( ) ? pi : used_paths[pi] );
          vec<int>& p = PathMutable(i);
          vec<int> p_old = p;
          for ( int j = 0; j < p.isize( ); j++ )
          {    if ( !BinMember( all, p[j] ) ) continue;

               // Find all the proper overlaps between p and some x[l], that use
               // the jth entry of p.

               vec< std::pair<int,int> > pos;
               for ( int l = 0; l < x.isize( ); l++ )
               {    for ( int m = 0; m < x[l].isize( ); m++ )
                    {    if ( x[l][m] != p[j] ) continue;
                         if ( Overlap( p, x[l], j - m ) ) pos.push(l,m);    }    }
               if ( !pos.solo( ) )
               {    p.clear( );
                    WeightFwMutable(i) = 0.0;
                    WeightRcMutable(i) = 0.0;
                    break;    }
               int l = pos[0].first, m = pos[0].second;
               int start = j;
               int stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
               for ( int r = start; r < stop - 1; r++ )
                    p[r] = -1;
               p[stop-1] = y[l];    }
          RemoveNegatives(p);

          if ( paths_index.nonempty( ) && p != p_old )
          {    for ( int l = 0; l < p_old.isize( ); l++ )
               {    int e = p_old[l];
                    #pragma omp critical
                    {    auto low = LowerBound( paths_index[e], {i,0} );
                         auto high = UpperBound( paths_index[e], {i,infty} );
                         paths_index[e].erase( 
                              paths_index[e].begin( ) + low,
                              paths_index[e].begin( ) + high );    }    }
               #pragma omp critical
               {    for ( int l = 0; l < p.isize( ); l++ )
                    {    int e = p[l];
                         paths_index[e].push( i, l );
                         Sort( paths_index[e] );    }    }    }    }

     #pragma omp parallel for
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> &p1 = PairLeftMutable(i), &p2 = PairRightMutable(i);
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<int>& p = ( pass == 1 ? p1 : p2 );
               for ( int j = 0; j < p.isize( ); j++ )
               {    if ( !BinMember( all, p[j] ) ) continue;
                    vec< std::pair<int,int> > pos;
                    for ( int l = 0; l < x.isize( ); l++ )
                    {    for ( int m = 0; m < x[l].isize( ); m++ )
                         {    if ( x[l][m] != p[j] ) continue;
                              if ( Overlap( p, x[l], j - m ) ) 
                                   pos.push(l,m);    }    }
                    if ( !pos.solo( ) )
                    {    p1.clear( ), p2.clear( );
                         PairDataMutable(i).clear( );
                         break;    }
                    int l = pos[0].first, m = pos[0].second;
                    int start = j, stop = Min( p.isize( ), x[l].isize( ) + (j-m) );
                    int add = 0;
                    if ( pass == 1 && j - m + x[l].isize( ) > p.isize( ) )
                    {    for ( int r = p.isize( ); r < x[l].isize( ) + j-m; r++ )
                              add += EdgeLengthKmers( x[l][ r - (j-m) ] );    }
                    if ( pass == 2 && j - m < 0 )
                    {    for ( int r = 0; r < -(j-m); r++ )
                              add += EdgeLengthKmers( x[l][r] );    }
                    AddTrim( i, add );
                    for ( int r = start; r < stop - 1; r++ )
                         p[r] = -1;
                    p[stop-1] = y[l];    }
               RemoveNegatives(p);    }    }    }


void SupportedHyperBasevector::writeBinary( BinaryWriter& writer ) const
{    HyperBasevector::writeBinary(writer);
     writer.write( Inv( ) );
     writer.write( Paths( ) );
     writer.write( WeightsFw( ) );
     writer.write( WeightsRc( ) );
     writer.write( WeightsFwOrigin( ) );
     writer.write( WeightsRcOrigin( ) );
     writer.write( Pairs( ) );
     writer.write( PairData( ) );
     writer.write( ReadCount( ) );
     writer.write( ReadLengthDist( ) );
     writer.write( FudgeMult( ) );    }

void SupportedHyperBasevector::readBinary( BinaryReader& reader )
{    HyperBasevector::readBinary(reader);
     reader.read( &InvMutable( ) );
     reader.read( &PathsMutable( ) );
     reader.read( &WeightsFwMutable( ) );
     reader.read( &WeightsRcMutable( ) );
     reader.read( &WeightsFwOriginMutable( ) );
     reader.read( &WeightsRcOriginMutable( ) );
     reader.read( &PairsMutable( ) );
     reader.read( &PairDataMutable( ) );
     reader.read( &ReadCountMutable( ) );
     reader.read( &ReadLengthDistMutable( ) );
     reader.read( &FudgeMultMutable( ) );    }

void SupportedHyperBasevector::DumpDot( const String& head,
     const vec<Bool>& invisible, const vec<String>& edge_color,
     const long_logging& logc, const Bool hide_inv, 
     const vec<String>& edge_names ) const
{    vec<Bool> hide;
     if (hide_inv) FlagEdgesForHiding( *this, Inv( ), hide, logc );
     else hide.resize( EdgeObjectCount( ), False );
     const Bool DOT_LABEL_CONTIGS = True;
     const Bool DOT_LABEL_VERTICES = False;
     vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
          lengths[i] = EdgeLengthKmers(i);
     vec<String> edge_id_names( EdgeObjectCount( ) );
     if ( edge_names.size( ) > 0 ) edge_id_names = edge_names;
     else
     {    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          {    if ( !InvDef(i) ) edge_id_names[i] = ToString(i);
               else 
               {    edge_id_names[i] = ToString(i) 
                         + "=" + ToString( Inv(i) ) + "'";    }    }    }
     Ofstream( dout, head + ".dot" );
     PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
          HyperBasevector::edge_label_info::DIRECT, &edge_id_names ),
          DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, &hide,
          &invisible, &edge_color, NULL, logc.LAYOUT );    }

void SupportedHyperBasevector::DumpFilesStandard( 
     const long_logging_control& log_control, const long_logging& logc,
     const int id ) const
{    if ( log_control.OUT_INT_HEAD != "" )
     {    DumpFiles( log_control.OUT_INT_HEAD + "." + ToString(id), 
               log_control, logc );    }    }

void SupportedHyperBasevector::DumpFiles( const String& head,
     const long_logging_control& log_control, const long_logging& logc ) const
{    double clock = WallClockTime( );
     
     // Output .shbv and .dot.

     BinaryWriter::writeFile( head + ".shbv", *this );
     vec<Bool> invisible( EdgeObjectCount( ), False );
     vec<String> edge_color( EdgeObjectCount( ), "" );
     DumpDot( head, invisible, edge_color, logc );

     // Output .support.

     Ofstream( sout, head + ".support" );
     for ( int i = 0; i < NPaths( ); i++ )
     {    int kmers = 0;
          for ( int j = 0; j < Path(i).isize( ); j++ )
               kmers += EdgeLengthKmers( Path(i)[j] );
          sout << "[" << i << "," << std::setiosflags(std::ios::fixed) << std::setprecision(1) 
               << Weight(i) << std::resetiosflags(std::ios::fixed) << "x] "
               << printSeq( Path(i) ) << " [" << kmers << " kmers]";
          // Temporary if.
          if ( WeightsFwOrigin( ).size( ) == WeightsFw( ).size( ) )
          {    vec<fix64_6> support;
               for ( int j = 0; j < WeightFwOrigin(i).isize( ); j++ )
                    support.push_back( WeightFwOrigin(i)[j].first );
               for ( int j = 0; j < WeightRcOrigin(i).isize( ); j++ )
                    support.push_back( WeightRcOrigin(i)[j].first );
               ReverseSort(support);
               sout << "; weights=";
               int count = 0;
               for ( int j = 0; j < support.isize( ); j++ )
               {    if ( j > 0 ) sout << ",";
                    if ( count++ > 5 ) 
                    {    sout << "...";
                         break;    }
                    int k = support.NextDiff(j);
                    sout << support[j] << "[" << k-j << "x]";
                    j = k - 1;    }    }
          sout << "\n";    }

     // Output .pair.

     Ofstream( pout, head + ".pairs" );
     for ( int i = 0; i < NPairs( ); i++ )
     {    vec<int> trims;
          for ( int j = 0; j < PairData(i).isize( ); j++ )
               trims.push_back( PairData(i,j).Trim( ) );
          Sort(trims);
          int Q1 = trims[ trims.size( ) / 4 ];
          int Q2 = trims[ trims.size( ) / 2 ];
          int Q3 = trims[ ( 3 * trims.size( ) ) / 4 ];
          pout << "\n[" << i << "] " << printSeqExp( PairLeft(i) ) << " --> "
               << printSeqExp( PairRight(i) ) << " (count=" << PairData(i).size( ) 
               << ", trim Q123 = " << Q1 << "," << Q2 << "," << Q3 << ")\n";    }
     pout << "\n";

     // Output .fasta.

     Ofstream( fout, head + ".fasta" );
     vec<int> to_left, to_right;
     ToLeft(to_left), ToRight(to_right);
     for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
     {    int v = to_left[e], w = to_right[e];
          fout << ">edge_" << e << " " << v << ":" << w 
               << " K=" << K( ) << " kmers=" << EdgeLengthKmers(e) << "\n";
          EdgeObject(e).Print(fout);    }    

     if(logc.MAKE_FASTG) fastg::WriteFastg(head+".fastg", *this);

     // Output .glocs and .gpaths.

     if ( logc.USE_GENOME_FOR_DUMP && log_control.G->size( ) > 0 )
     {    Ofstream( out, head + ".glocs" );
          Ofstream( pout, head + ".gpaths" );

          // Compute placements of edges.

          vec< vec<placementy> > places( EdgeObjectCount( ) );
          #pragma omp parallel for
          for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
               places[e] = log_control.FindGenomicPlacements( EdgeObject(e) );
          for ( int e = 0; e < EdgeObjectCount( ); e++ ) 
          {    int v = to_left[e], w = to_right[e];
               out << "edge_" << e << "[v=" << v << ",w=" << w << ",l=" 
                    << EdgeObject(e).size( ) << "]:";
               for ( int j = 0; j < places[e].isize( ); j++ )
               {    placementy& p = places[e][j];
                    int n = (*log_control.G)[p.g].size( );
                    if ( p.Pos < K( ) - 1 ) p.Pos += n;
                    out << " " << ( p.fw ? "fw" : "rc" ) << p.g << ".";
                    if (p.fw) out << p.pos << "-" << p.Pos;
                    else out << n - p.Pos << "-" << n - p.pos;    }
               out << "\n";    }    

          // Define a data structure P that assigns to each oriented edge the
          // vector of its genomic placements.

          int n = EdgeObjectCount( );
          vec< vec< triple<int,int,int> > > P(2*n);
          const vec<bool>& is_circular = *log_control.is_circular;
          for ( int e = 0; e < n; e++ )
          {    for ( int j = 0; j < places[e].isize( ); j++ )
               {    const placementy& p = places[e][j];
                    int start = p.pos;
                    int stop = start + EdgeObject(e).isize( );
                    int n = (*log_control.G)[p.g].size( );
                    if ( is_circular[p.g] && stop >= n + K( ) - 1 ) stop -= n;
                    stop -= ( K( ) - 1 );
                    P[ 2*e + ( p.fw ? 0 : 1 ) ].push( p.g, start, stop );    }    }

          // Build graph.

          // int N = 0;
          vec< triple<int,int,int> > V;
          vec<int> PS;
          PS.push_back(0);
          for ( int j = 0; j < P.isize( ); j++ )
          {    // N += P[j].size( );
               PS.push_back( PS.back( ) + P[j].isize( ) );
               for ( int l = 0; l < P[j].isize( ); l++ )
                    V.push( j, P[j][l].first, P[j][l].second );    }
          int N = PS.back( );
          vec< vec<int> > from(N), to(N);
          const int batch = 100;
          #pragma omp parallel for
          for ( int bz1 = 0; bz1 < this->N( ); bz1 += batch )
          {    vec< std::pair<int,int> > edges;
               for ( int z1 = bz1; z1 < Min( bz1 + batch, this->N( ) ); z1++ )
               {    for ( int j1 = 0; j1 < From(z1).isize( ); j1++ )
                    {    int z2 = From(z1)[j1];
                         int v1 = EdgeObjectIndexByIndexFrom( z1, j1 );
                         if ( P[2*v1].empty( ) && P[1+2*v1].empty( ) ) continue;
                         for ( int j2 = 0; j2 < From(z2).isize( ); j2++ )
                         {    int z3 = From(z2)[j2];
                              int v2 = EdgeObjectIndexByIndexFrom( z2, j2 );
                              if ( P[2*v2].empty( ) && P[1+2*v2].empty( ) ) continue;
                              int x1, x2;
                              for ( int pass = 1; pass <= 2; pass++ )
                              {    if ( pass == 1 )
                                   {    x1 = 2*v1, x2 = 2*v2;    }
                                   else
                                   {    x2 = 2*v1 + 1, x1 = 2*v2 + 1;    }
                                   vec< triple< std::pair<int,int>, int, int > > X;
                                   for ( int l1 = 0; l1 < P[x1].isize( ); l1++ )
                                   {    X.push( std::make_pair( P[x1][l1].first,
                                             P[x1][l1].third ), 0, l1 );    }
                                   for ( int l2 = 0; l2 < P[x2].isize( ); l2++ )
                                   {    X.push( std::make_pair( P[x2][l2].first,
                                             P[x2][l2].second ), 1, l2 );    }
                                   Sort(X);
                                   for ( int i = 0; i < X.isize( ); i++ )
                                   {    int j, mid;
                                        for ( j = i + 1; j < X.isize( ); j++ )
                                             if ( X[j].first != X[i].first ) break;
                                        for ( mid = i + 1; mid < j; mid++ )
                                             if ( X[mid].second == 1 ) break;
                                        for ( int u1 = i; u1 < mid; u1++ )
                                        for ( int u2 = mid; u2 < j; u2++ )
                                        {    int l1 = X[u1].third, l2 = X[u2].third;
                                             int m1 = l1, m2 = l2;
                                             m1 += PS[x1];
                                             m2 += PS[x2];
                                             edges.push( m1, m2 );    }    
                                        i = j - 1;    }    }    }    }    }    
               #pragma omp critical
               {    for ( int i = 0; i < edges.isize( ); i++ )
                    {    from[ edges[i].first ].push_back( edges[i].second );
                         to[ edges[i].second ].push_back( 
                              edges[i].first );    }    }   }
          #pragma omp parallel for
          for ( int i = 0; i < N; i++ )
          {    Sort(from[i]), Sort(to[i]);    }
          digraph H( from, to );
          vec< triple< int, std::pair<int,int>, vec<int> > > matches;
          vec< vec<int> > paths;
          H.AllPaths( -1, -1, paths );

          if ( Sum(is_circular) > 0 ) // wrong, need granular test
          {    vec< vec<int> > all_loops;
               for ( int v = 0; v < N; v++ )
               {    vec< vec<int> > loops;
                    H.AllPaths( v, v, loops, -1, True );
                    for ( int j = 0; j < loops.isize( ); j++ )
                    {    vec<int> x = loops[j];
                         if ( x.solo( ) ) continue;
                         int m = Min(x);
                         int p;
                         for ( p = 0; p < x.isize( ); p++ )
                              if ( x[p] == m ) break;
                         vec<int> y;
                         for ( int i = p; i < x.isize( ) - 1; i++ )
                              y.push_back( x[i] );
                         for ( int i = 0; i < p; i++ )
                              y.push_back( x[i] );
                         if ( !Member( all_loops, y ) ) 
                              all_loops.push_back(y);    }    }
               paths.append(all_loops);    }

          for ( int j = 0; j < paths.isize( ); j++ )
          {    const vec<int>& p = paths[j];
               int g = V[ p.front( ) ].second, pos = V[ p.front( ) ].third;
               int Pos = V[ p.back( ) ].third
                    + EdgeObject( V[ p.back( ) ].first / 2 ).isize( );
               vec<int> q;
               for ( int l = 0; l < paths[j].isize( ); l++ )
               {    int e = V[ paths[j][l] ].first;
                    if ( e % 2 == 0 ) q.push_back(e/2);
                    else q.push_back( -e/2-1 );    }
               matches.push( g, std::make_pair( pos, Pos ), q );    }
          Sort(matches);
          for ( int j = 0; j < matches.isize( ); j++ )
          {    int g = matches[j].first;
               int n = (*log_control.G)[g].isize( );
               int start = matches[j].second.first, stop = matches[j].second.second;
               if ( is_circular[g] && stop >= n + K( ) - 1 ) stop -= n;
               pout << g << "." << start << "-" << stop;
               if ( stop - start == K( ) - 1 ) pout << " (perfect circle)";
               pout << " :: ";
               const vec<int>& v = matches[j].third;
               for ( int l = 0; l < v.isize( ); l++ )
               {    if ( l > 0 ) pout << ",";
                    if ( v[l] >= 0 ) pout << v[l] << "fw";
                    else pout << -v[l]-1 << "rc";    }
               pout << "\n";    }    }

     REPORT_TIME( clock, "used dumping files" );    }


