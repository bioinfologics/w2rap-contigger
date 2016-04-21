///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EdgeInfo. Get info associated to particular edges in a GapToy assembly.

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "paths/long/ReadPath.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_String_OrDefault_Doc(E, "", 
          "boolean expression involving edge ids, in "
          "either the format e1|...|en or e1&...&en; blank space allowed for "
          "readability");
     CommandArgument_Int_OrDefault_Doc(R, -1, "look up this read id");
     CommandArgument_Int_OrDefault_Doc(P, -1, "look up this pair id");
     CommandArgument_Bool_OrDefault_Doc(PIDS, False, "show pids");
     CommandArgument_Bool_OrDefault_Doc(WRITE_READ_IDS, False, "write read ids");
     CommandArgument_Bool_OrDefault_Doc(SCRUT, False, 
          "scrutinize placements of reads on the given edges; rudimentary at best");
     CommandArgument_String_OrDefault_Doc(OUT_READ_HEAD, "", 
          "if specified, dump reads to this.{fasta,qualb}");
     CommandArgument_Bool_OrDefault_Doc(OUT_ORIENT, False, 
          "if OUT_READ_HEAD specified, swap order of reads within a pair as needed "
          "to make the first read forward, if possible");
     EndCommandArguments;
         
     // Parse E.

     vec<int> ors, ands;
     if ( E.nonempty( ) )
     {    E.GlobalReplaceBy( " ", "" );
          if ( E.IsInt( ) ) ors.push_back( E.Int( ) );
          else if ( E.Contains( "|" ) )
          {    E.GlobalReplaceBy( "|", "," );
               ParseIntSet( "{" + E + "}", ors );    }
          else if ( E.Contains( "&" ) )
          {    E.GlobalReplaceBy( "&", "," );
               ParseIntSet( "{" + E + "}", ands );    }
          else
          {    cout << "Illegal value for E." << std::endl;
               Scram(1);    }    }

     // Load inversion.

     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );

     // Look up edges.

     vec<int> e(ors);
     e.append(ands);
     for ( int i = 0; i < e.isize( ); i++ )
     {    if ( e[i] >= inv.isize( ) )
          {    cout << e[i] << " is not a valid edge id." << std::endl;
               Scram(1);    }    }
     int ne = e.size( );
     vec<int> e2(e);
     for ( int i = 0; i < ne; i++ )
               if ( inv[ e[i] ] >= 0 ) e2.push_back( inv[ e[i] ] );
     UniqueSort(e2);
     VecULongVec x0;
     x0.Read( DIR + "/a.paths.inv", e2 );
     vec<vec<int64_t>> x( x0.size( ) );
     for ( int i = 0; i < x.isize( ); i++ )
     for ( int j = 0; j < (int) x0[i].size( ); j++ )
          x[i].push_back( x0[i][j] );

     // Define pids.

     vec<int64_t> pids;
     if ( ors.nonempty( ) )
     {    for ( int j = 0; j < e2.isize( ); j++ )
          for ( int i = 0; i < (int) x[j].size( ); i++ )
               pids.push_back( x[j][i]/2 );
          UniqueSort(pids);    }
     else if ( ands.nonempty( ) )
     {    vec<vec<int64_t>> pidsx( ands.size( ) );
          for ( int i = 0; i < ands.isize( ); i++ )
          {    pidsx[i] = x[ BinPosition( e2, ands[i] ) ];
               if ( inv[ ands[i] ] >= 0 )
                    pidsx[i].append( x[ BinPosition( e2, inv[ ands[i] ] ) ] );
               for ( int j = 0; j < pidsx[i].isize( ); j++ )
                    pidsx[i][j] /= 2;
               UniqueSort(pidsx[i]);    }
          Intersection( pidsx, pids );    }
     else if ( R >= 0 ) pids.push_back(R/2);
     else if ( P >= 0 ) pids.push_back(P);
     int64_t npids = pids.size( );

     // Fetch paths.

     vec<int64_t> ids;
     for ( int i = 0; i < npids; i++ )
          ids.push_back( 2*pids[i], 2*pids[i] + 1 );
     ReadPathVec p;

     int64_t nreads = MastervecFileObjectCount( DIR + "/a.paths" );
     for ( int i = 0; i < ids.isize( ); i++ )
     {    if ( ids[i] < 0 )
          {    cout << "\nRead id " << ids[i] << " is negative." << std::endl;
               Scram(1);    }
          if ( ids[i] >= nreads )
          {    cout << "\nRead id " << ids[i] << " is >= nreads = " << nreads
                    << "." << std::endl;
               Scram(1);    }    }
     p.Read( DIR + "/a.paths", ids );

     // Display results.

     vec< int > good_ids;       // ids to write out
     vec< quad< vec<int>, vec<int>, int, Bool > > stuff;
     int bads = 0;
     for ( int i = 0; i < npids; i++ )
     {    const ReadPath &p1 = p[2*i], &p2 = p[2*i+1];
          Bool bad = False;
          for ( int i = 0; i < (int) p1.size( ); i++ )
               if ( p1[i] >= inv.isize( ) ) bad = True;
          for ( int i = 0; i < (int) p2.size( ); i++ )
               if ( p2[i] >= inv.isize( ) ) bad = True;

          if (bad) 
          {    bads++;
               continue;    } // CONTINUE
          else good_ids.push_back( 2*pids[i], 2*pids[i]+1 );

          vec<int> x1, x2, y1, y2;
          for ( int i = 0; i < (int) p1.size( ); i++ )
               x1.push_back( p1[i] );
          for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
               x2.push_back( inv[ p2[i] ] );
          if ( E == "" || Meet2( x1, e ) || Meet2( x2, e ) ) 
               stuff.push( x1, x2, pids[i], True );
          for ( int i = 0; i < (int) p2.size( ); i++ )
               y2.push_back( p2[i] );
          for ( int i = ( (int) p1.size( ) ) - 1; i >= 0; i-- )
               y1.push_back( inv[ p1[i] ] );
          if ( E == "" || Meet2( y1, e ) || Meet2( y2, e ) ) 
               stuff.push( y2, y1, pids[i], False );    }

     // Write read IDS if required

     if (WRITE_READ_IDS) {
       Ofstream(out, "a.read_ids");
       for (size_t i = 0; i < good_ids.size(); i++)
             out << good_ids[i] << std::endl;
     }

     if ( bads > 0 ) std::cout << "\n" << bads << " of " << npids << " bad pids" << std::endl;
     cout << "\n";
     Sort(stuff);
     vec<vec<String>> rows;
     vec<int> all;
     for ( int i = 0; i < stuff.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < stuff.isize( ); j++ )
          {    if ( stuff[j].first != stuff[i].first 
                    || stuff[j].second != stuff[i].second )
               {    break;    }    }
          ostringstream out1, out2, out3;
          out1 << "[" << j-i << "]";
          out2 << printSeq( stuff[i].first ) << ".." << printSeq( stuff[i].second );
          for ( int k = i; k < j; k++ )
          {    all.append( stuff[i].first );
               all.append( stuff[i].second );    }
          if (PIDS)
          {    out3 << "  [pid=";
               for ( int k = i; k < j; k++ )
               {    if ( k > i ) out3 << ",";
                    if ( stuff[k].fourth ) out3 << "+";
                    else out3 << "-";
                    out3 << stuff[k].third;    }
               out3 << "]";    }
          vec<String> row = { out1.str( ), out2.str( ) };
          if (PIDS) row.push_back( out3.str( ) );
          rows.push_back(row);
          i = j - 1;    }
     vec<String> row = { "[" + ToString( stuff.size( ) ) + "]", "TOTAL" };
     rows.push_back(row);
     PrintTabular( cout, rows, 2 );
     cout << "\n";
     Sort(all);
     vec<int> all2;
     for ( int i = 0; i < all.isize( ); i++ )
     {    int j = all.NextDiff(i);
          if ( j - i >= 2 ) all2.push_back( all[i] );
          i = j - 1;    }
     UniqueSort(all);
     cout << "all edges: " << printSeq(all) << std::endl;
     cout << "all non-solo edges: " << printSeq(all2) << std::endl;
     
     // Dump reads.

     vecbasevector bases;
     vecqualvector quals;
     if ( OUT_READ_HEAD != "" || SCRUT )
     {    bases.Read( DIR + "/../data/frag_reads_orig.fastb", ids );
          if ( IsRegularFile( DIR + "/../data/frag_reads_orig.qualp" ) )
          {    VecPQVec qualsp;
               qualsp.Read( DIR + "/../data/frag_reads_orig.qualp", ids );
               quals.resize( qualsp.size( ) );
               for ( int i = 0; i < (int) quals.size( ); i++ )
                    qualsp[i].unpack( &quals[i] );    }
          else // for backward compatibility
          {    quals.Read( DIR + "/../data/frag_reads_orig.qualb", ids );    }
          if ( OUT_READ_HEAD != "" )
          {    vecqualvector qualsx(quals);
               vec<int> nid( ids.size( ), vec<int>::IDENTITY );
               if (OUT_ORIENT)
               {    for ( int i = 0; i < ids.isize( ); i += 2 )
                    {    int id1 = ids[i], id2 = ids[i] + 1;
                         Bool fw1 = False, fw2 = False;
                         for ( int j = 0; j < (int) p[i].size( ); j++ )
                         {    if ( Member( e, p[i][j] ) ) fw1 = True;
                              if ( Member( e, inv[p[i][j]] ) ) fw2 = True;    }
                         for ( int j = 0; j < (int) p[i+1].size( ); j++ )
                         {    if ( Member( e, p[i+1][j] ) ) fw2 = True;
                              if ( Member( e, inv[p[i+1][j]] ) ) fw1 = True;    }
                         if ( fw2 && !fw1 )
                         {    swap( nid[i], nid[i+1] );
                              swap( qualsx[i], qualsx[i+1] );    }    }    }
               Ofstream( out, OUT_READ_HEAD + ".fasta" );
               for ( int i = 0; i < nid.isize( ); i++ )
                    bases[ nid[i] ].Print( out, ids[ nid[i] ] );
               qualsx.WriteAll( OUT_READ_HEAD + ".qualb" );    }    }
     
     // Scrutinize placements of reads on the edges.

     if (SCRUT)
     {    cout << "\nscrutinizing placements:\n";
          vec<int> eplus(e);
          for ( int i = 0; i < stuff.isize( ); i++ )
          {    eplus.append( stuff[i].first );
               eplus.append( stuff[i].second );    }
          UniqueSort(eplus);
          vecbasevector edges;
          edges.Read( DIR + "/a.fastb", eplus );
          for ( int i = 0; i < npids; i++ )
          {    const ReadPath &p1 = p[2*i], &p2 = p[2*i+1];
               int pid = pids[i];
               Bool problems = False;
               ostringstream out;
               out << "\n";
               PRINT_TO( out, pid );
               Bool bad = False;
               for ( int i = 0; i < (int) p1.size( ); i++ )
                    if ( p1[i] >= inv.isize( ) ) bad = True;
               for ( int i = 0; i < (int) p2.size( ); i++ )
                    if ( p2[i] >= inv.isize( ) ) bad = True;
               if (bad)
               {    out << "BAD!\n";
                    problems = True;
                    cout << out.str( );
                    continue;    }
               vec<int> x1, x2, y1, y2;
               for ( int i = 0; i < (int) p1.size( ); i++ )
                    x1.push_back( p1[i] );
               for ( int i = ( (int) p2.size( ) ) - 1; i >= 0; i-- )
                    x2.push_back( inv[ p2[i] ] );
               if ( Meet2( x1, e ) || Meet2( x2, e ) ) 
               {    if ( x1.empty( ) ) 
                    {    out << "p1 unplaced" << std::endl;
                         problems = True;    }
                    else if ( p1.getOffset( ) < 0 )
                    {    out << "p1 has negative start" << std::endl;
                         problems = True;    }
                    else
                    {    int start = p1.getOffset( );
                         out << start << ":" << printSeq(p1) << std::endl;
                         int e = p1[0];
                         ForceAssertGe( e, 0 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         int ex = BinPosition( eplus, e );
                         ForceAssertGe( ex, 0 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                         int id = 2*pid;
                         int idx = BinPosition( ids, id );
                         const basevector& b = bases[idx];
                         const qualvector& q = quals[idx];
                         int pp = 0;
                         const int K = 200;
                         for ( int pos = 0; pos < b.isize( ); pos++ )
                         {    if ( b[pos] != edges[ex][start+pos] && q[pos] > 2 )
                              {    out << "diff at rpos = " << pos << ", "
                                        << e << "." << start+pos << " --> "
                                        << as_base( b[pos] )
                                        << "[" << int(q[pos]) << "]\n";
                                   problems = True;    }
                              if ( start + pos == edges[ex].isize( ) - 1 )
                              {    pp++;
                                   if ( pp == (int) p1.size( ) )
                                   {    out << "at end of path\n";
                                        problems = True;
                                        break;    }
                                   start -= edges[ex].isize( ) - (K-1);
                                   e = p1[pp];
                                   ForceAssertGe( e, 0 ); // XXXXXXXXXXXXXXXXXXXXXXX
                                   ex = BinPosition( eplus, e );    
                                   ForceAssertGe( ex, 0 ); // XXXXXXXXXXXXXXXXXXXXXX
                                        }    }    }
                    if (problems) std::cout << out.str( );    }    }    }    }
