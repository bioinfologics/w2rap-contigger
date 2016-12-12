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
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "graphics/BasicGraphics.h"
#include "kmers/LongReadPather.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/long/Correct1Pre.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/HBVFromEdges.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/GapToyTools.h"
#include <ctime>


void FixPath( ReadPath& p, const vec< triple<int,int,int> >& merges,
     const vec< vec<int> >& merges_index, const HyperBasevector& hb_orig )
{
     // mark path edges no longer in the hbv as id=-1
     for ( int64_t l = 0; l < (int64_t) p.size( ); l++ )
          if ( p[l] >= hb_orig.EdgeObjectCount( ) ) p[l] = -1;
     while(1)
     {    Bool changed = False;
          // for each valid path edge (not -1), store and sort a
          // list of merge ids involving the edge
          vec<int> mids;
          for ( int l = 0; l < (int) p.size( ); l++ )
               if ( p[l] >= 0 ) mids.append( merges_index[ p[l] ] );
          UniqueSort(mids);
          // for each merge id, iteratively update the path to the
          // new (merged) edges
          for ( int mj = 0; mj < mids.isize( ); mj++ )
          {    int j = mids[mj];
               int e1 = merges[j].first, e2 = merges[j].second;
               ForceAssert( e1 != e2 );
               int enew = merges[j].third;
               for ( int64_t l = 0; l < (int64_t) p.size( ); l++ )
               {
                    // case 1: adjacent path edges cover a merged pair
                    //           p[l] == e1, p[l+1] == e2 -> p[l]=p[l+1]=enew
                    if ( l < ( (int) p.size( ) ) - 1 && p[l] == e1 && p[l+1] == e2 )
                    {    p[l] = enew;
                         for ( int m = l+2; m < (int) p.size( ); m++ )
                              p[m-1] = p[m];
                         p.pop_back( );
                         changed = True;
                    }
                    // case 2: last path edge was start of a merged pair
                    else if ( l == ( (int) p.size( ) ) - 1 && p[l] == e1 )
                    {
/*                         size_t lastSkip = p.getLastSkip();
                         p.setLastSkip(lastSkip + hb_orig.EdgeLengthKmers(e2)); */
                         p[l] = enew;
                         changed = True;
                    }
                    // case 3: first path edge was end of a merged pair
                    else if ( l == 0 && p[l] == e2 )
                    {    int offset = p.getOffset();
                         p.setOffset( offset + hb_orig.EdgeLengthKmers(e1));
                         p[l] = enew;
                         changed = True;    }
                    if (changed) break;    }
               if (changed) break;    }

          if ( !changed ) break;   // dump out if we made a clean pass

     }    // while(1)
}

void RemoveUnneededVertices2( HyperBasevector& hbv, vec<int>& inv, ReadPathVec& paths, Bool debug )
{
    static int debug_serial = 0;
    ++debug_serial;
    std::ostringstream debug_fnam_head;
    debug_fnam_head << "graph." << std::setprecision(3) << debug_serial;


    if ( debug ) BinaryWriter::writeFile( debug_fnam_head.str() + ".BEFORE.hbv", hbv );

    // new algorithm
    // 1. make a list of vertices to kill
    // 2. find the "boundary" edges at the front and tail of the run to remove
    // 3. walk edges between/including boundary edges, building up new edge object,
    //    and recording edge mappings.
    // 4. add new edge object between boundary vertices;
    // 5. delete all edges marked for deletion
    //
    // this all assumes that the involution can't share edges in any
    // string of edges that we're interested in.  The only way I can see
    // this happening is a single, palindromic edge, but that would not
    // qualify for

    vec<bool> vertex_kill( hbv.N(), false );
    vec<int> vertex_queue;
    vec<int> to_left, to_right;
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);
    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Step1"<< std::endl;
    // step 1: make a list of vertices to kill
    // o----o----o----o
    // v0   v1   v2   v3
    for ( int v = 0; v < hbv.N(); ++v )
        if ( hbv.FromSize(v) == 1 && hbv.ToSize(v) == 1
                && hbv.From(v)[0] != hbv.To(v)[0] 
                && hbv.Bases( hbv.IFrom( v, 0 ) ) > 0
                && hbv.Bases( hbv.ITo( v, 0 ) ) > 0 ) {
            vertex_kill[v] = true;
            vertex_queue.push_back(v);
        }


    // step 2
    //time(&rawtime);
    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Step2"<< std::endl;
    vec<std::pair<int,int>> bound;
    while ( vertex_queue.size() ) {
        int v = vertex_queue.back();
        vertex_queue.pop_back();
        if ( !vertex_kill[v] ) continue;     // already done it
        int eleft;
        int vleft = v;
        size_t runsize = 0;
        do {
            if ( debug ) std::cout << "vleft = " << vleft << std::endl;
            runsize++;
            vertex_kill[vleft] = false;
            eleft = hbv.EdgeObjectIndexByIndexTo(vleft,0);
            vleft = hbv.To(vleft)[0];
        } while ( vertex_kill[vleft] );
        int eright;
        int vright = v;
        do {
            if ( debug ) std::cout << "vright = " << vright << std::endl;
            runsize++;
            vertex_kill[vright] = false;
            eright = hbv.EdgeObjectIndexByIndexFrom(vright,0);
            vright = hbv.From(vright)[0];
        } while ( vertex_kill[vright] );
        runsize--;      // 'v' gets counted twice

        // We rely on the fact that the involution is not
        // tied up with the path here.  We decide to push on the involution
        // here, too, so that we *know* what the inv[] of the new edge is.  However,
        // this requires that we canonicalize, so we don't do this twice.  This
        // canonicalization looks odd, but is correct (I think).
        if ( eleft < inv[eright] ) {
            // WARNING: code below relies on the fact that we're pushing on a
            // run and its involution adjacent in this list.
            bound.push(eleft,eright);
            bound.push(inv[eright], inv[eleft]);

            if ( debug ) {
                std::cout << "eleft = " << eleft << ", eright = " << eright << std::endl;
                std::cout << "inv eleft = " << inv[eleft] <<
                        ", inv eright = " << inv[eright] << std::endl;
                std::cout << "runsize = " << runsize << std::endl;
                std::cout << "===" << std::endl;
            }
        }
    }
    if ( debug ) PRINT(bound.size());

    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Step3"<< std::endl;
    // steps 3 and 4
    
    //XXX Optimization START (Gonza & BJ - 2015-09-07)
    //Variables:
    //  - Function scope (updated and returned):
    //              - hbv - AKA the graph
    //              - inv - AKA the graph transformation
    //              - paths - AKA the paths :D
    //  - Function scope (received):
    //              - bound - AKA the start-end of paths to simplify (vector of pairs)
    //
    //  - Inner scope:
    //              - edge_renumber0 = mapping to new edges 
    //              - offsets (where to find the old edge in the new edge)
    //              - to_delete
    //              - new_edge_numbers (just to update inv)
    //
    vec<int> edge_renumber0( hbv.EdgeObjectCount(), vec<int>::IDENTITY );
    vec<int> offsets(hbv.EdgeObjectCount(),0);
    vec<int> new_edge_numbers;
    vec<int> to_delete;
    //XXX: bound.size() is how many new edges we'll need, so we can pre-allocate that.
    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Step4"<< std::endl;
    while ( bound.size() ) {///TODO: we can make the Adds at the end and the paralelise this easily
        auto bounds = bound.back();
        bound.pop_back();

        int new_edge_no = hbv.EdgeObjectCount();

        int off = hbv.EdgeLengthKmers( bounds.first );
        edge_renumber0[ bounds.first ] = new_edge_no;
        to_delete.push_back( bounds.first );
        //First compute total size and offsets
        //std::cout<<"Computing size and offsets for run from e-v="<<to_right[bounds.first]<<"-> to e-v="<<to_right[bounds.second]<<"->"<<std::endl;
        for ( int v = to_right[bounds.first]; v != to_right[bounds.second]; v = hbv.From(v)[0] ) {
            //std::cout<<"Retrieving edgeFrom["<<v<<"[0]"<<std::endl;
            int edge = hbv.EdgeObjectIndexByIndexFrom(v,0);
            //std::cout<<"edge #"<<edge<<std::endl;
            to_delete.push_back(edge);
            offsets[edge] = off;
            edge_renumber0[edge] = new_edge_no;
            off += hbv.EdgeLengthKmers(edge);
        }
        //off holds the total size now!
        //Now reserve memory, create a new edge and copy things into their appropriate place
        //std::cout<<"Replacing run from e-v="<<to_right[bounds.first]<<"-> to e-v="<<to_right[bounds.second]<<"->"<<std::endl;
        basevector new_edge( hbv.EdgeObject( bounds.first ) );
        new_edge.reserve(off+hbv.K()-1);
        for ( int v = to_right[bounds.first]; v != to_right[bounds.second]; v = hbv.From(v)[0] ) {
          //std::cout<<"Retrieving edgeFrom["<<v<<"[0]"<<std::endl;
          int edge = hbv.EdgeObjectIndexByIndexFrom(v,0);
          //std::cout<<"Fetching edge object"<<std::endl;
          auto edge_ob = hbv.EdgeObject(edge);
          new_edge.resize(offsets[edge]);
          new_edge.append(edge_ob.begin(),edge_ob.end());
          //std::cout<<"Edge appended"<<std::endl;
        }

        
        //std::cout<<"Adding new edge between vertices"<<to_left[bounds.first]<<" and "<<to_right[bounds.second]<<std::endl;
        hbv.AddEdge(to_left[bounds.first], to_right[bounds.second], new_edge); //XXX: this modifies the graph topology we already are using
        //std::cout << bounds.first << "-" << bounds.second <<
            //        "->" << new_edge_no << std::endl;
        if ( debug ) {
            std::cout << "run from edge " << bounds.first << " to " << bounds.second <<
                    " replaced by edge " << new_edge_no << std::endl;
        }
        new_edge_numbers.push_back(new_edge_no);
    }
    //time(&rawtime);
    //std::cout << "[GapToyTools3.cc]          RemoveUnneededVertices2 edges before deletion: "<<hbv.EdgeObjectCount()<< std::endl;
    hbv.DeleteEdges(to_delete);
    //time(&rawtime);
    //std::cout << "[GapToyTools3.cc]          RemoveUnneededVertices2 edges after deletion: "<<hbv.EdgeObjectCount()<< std::endl;

    if ( debug )
        BinaryWriter::writeFile( debug_fnam_head.str() + ".AFTER.hbv", hbv );

    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Updating inv[x]"<< std::endl;
    // for each pair of newly created edges, update mInv
    inv.resize(hbv.EdgeObjectCount() );
    for (auto itr = new_edge_numbers.begin();
            itr != new_edge_numbers.end(); advance(itr,2)) {
         inv[itr[0]] = itr[1];
         inv[itr[1]] = itr[0];
    }
    //time(&rawtime);
    //std::cout << "[GapToyTools3.cc]          RemoveUnneededVertices2 updating paths with new edge numbers: " << ctime(&rawtime) << std::endl;
    // update the read paths for the newly created edges
    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Updating paths to new edge numbers"<< std::endl;
#pragma omp parallel for
    for ( int64_t i = 0; i < paths.size(); ++i ) {
        auto& path = paths[i];
        if ( path.size() ) {
            std::vector<int> oldPath = path;
            path.clear();
            path.setOffset( path.getOffset() + offsets[oldPath[0]]);
            //XXX: is this the same scenario?
            path.push_back( edge_renumber0[ oldPath[0] ] );
            for ( auto itr = oldPath.begin()+1; itr != oldPath.end(); ++itr )
                if ( edge_renumber0[ *itr ] != path.back() )
                    path.push_back( edge_renumber0[ *itr ] );
        }
    }
    //std::cout << "[GapToyTools3.cc] Begining RemoveUnneededVertices2 Finished!"<< std::endl;
    //XXX Optimization END
    //XXX: this is PROPERTY VALIDATION IN PRODUCTION! AWESOME!
    //time(&rawtime);
    //std::cout << "[GapToyTools3.cc]          RemoveUnneededVertices2 Begining Validate: " << ctime(&rawtime) << std::endl;
    //Validate( hbv, inv, paths );
    //TestInvolution(hbv, inv);
    //std::cout << "[GapToyTools3.cc] Remove unneeded vertices toleft finishing: " << to_left.size() << std::endl;
    //std::cout << "[GapToyTools3.cc] Remove unneeded vertices toright finishing: " << to_right.size() << std::endl;
}

void RemoveUnneededVerticesLoopsOnly( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths )
{    double clock1 = WallClockTime( );
     vec< triple<int,int,int> > merges;
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     for ( int i = 0; i < hb.N( ); i++ )
     {    if ( hb.From(i).size( ) == 1 && hb.To(i).size( ) == 1 
               && hb.From(i)[0] != i )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( i, 0 );
               int e2 = hb.EdgeObjectIndexByIndexFrom( i, 0 );
               if ( hb.EdgeLengthBases(e1) == 0 || hb.EdgeLengthBases(e2) == 0 )
                    continue;
               int v = hb.To(i)[0], w = hb.From(i)[0];
               Bool loop = ( v == w && hb.From(v).solo( ) && hb.To(v).solo( ) );
               if ( !loop ) continue;
               basevector p = hb.Cat( e1, e2 );
               int enew = hb.EdgeObjectCount( );
               int re1 = inv[e1], re2 = inv[e2];

               // v --e1--> i --e2--> w

               // Test for nasty condition.  The archetypal case here is 
               // e1 = (AT)^100, e2 = (TA)^100.  It is not possible to merge the 
               // edges and maintain consistent data structures.

               if ( re1 == e1 && re2 == e2 ) continue;

               // Start the merge.

               merges.push( e1, e2, enew );
               hb.JoinEdges( i, p );
               to_left.push_back(v), to_right.push_back(w);

               // Case 1.

               if ( re1 == e2 && re2 == e1 )
               {    // * --e1=re2--> * --e2=re1--> *

                    inv.push_back(enew);    }

               // Case 3.
               else
               {    // e1, e2, re1, re2 all different
                    int renew = hb.EdgeObjectCount( );
                    basevector rp = hb.Cat( re2, re1 );
                    merges.push( re2, re1, renew );
                    int ri = to_right[re2];
                    hb.JoinEdges( ri, rp );
                    int rv = to_left[re2], rw = to_right[re1];
                    to_left.push_back(rv), to_right.push_back(rw);
                    inv.push_back(renew, enew);    }    }    }
     vec< vec<int> > merges_index( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < merges.isize( ); i++ )
     {    merges_index[ merges[i].first ].push_back(i);
          merges_index[ merges[i].second ].push_back(i);    }
     LogTime( clock1, "removing unneeded vertices 1" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for
     for ( int64_t i = 0; i < (int64_t) paths.size( ); i++ )
          FixPath( paths[i], merges, merges_index, hb );
     LogTime( clock2, "removing unneeded vertices 2" );
     double clock3 = WallClockTime( );
     hb.RemoveEdgelessVertices( );
     LogTime( clock3, "removing unneeded vertices 3" );    }

void RemoveUnneededVerticesGeneralizedLoops( HyperBasevector& hb, vec<int>& inv,
     ReadPathVec& paths )
{    double clock1 = WallClockTime( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vec<Bool> processed( hb.N( ), False );
     vec<int> dels;
     for ( int i = 0; i < hb.N( ); i++ )
     {    if ( !processed[i] && hb.From(i).size( ) == 1 && hb.To(i).size( ) == 1 
               && hb.From(i)[0] != i )
          {    vec<int> chain;
               int v = i;
               Bool fail = False;
               while(1)
               {    chain.push_back(v);
                    v = hb.From(v)[0];
                    if ( hb.From(v).size( ) != 1 || hb.To(v).size( ) != 1 
                         || hb.From(v)[0] == v )
                    {    fail = True;
                         break;    }
                    if ( Member( chain, v ) ) break;    }
               if (fail) continue;
               vec<int> echain, rechain;
               for ( int j = 0; j < chain.isize( ); j++ )
                    echain.push_back( hb.IFrom( chain[j], 0 ) );
               for ( int j = echain.isize( ) - 1; j >= 0; j-- )
                    rechain.push_back( inv[ echain[j] ] );
               if ( Meet2( echain, rechain ) ) continue;
               basevector x = hb.Cat(echain), rx = hb.Cat(rechain);
               dels.append(echain), dels.append(rechain);
               for ( int j = 0; j < chain.isize( ); j++ )
               {    processed[to_left[chain[j]]] = True;
                    processed[to_right[chain[j]]] = True;
                    processed[to_left[rechain[j]]] = True;
                    processed[to_right[rechain[j]]] = True;    }
               hb.AddEdge( to_left[echain[0]], to_left[echain[0]], x );
               hb.AddEdge( to_left[rechain[0]], to_left[rechain[0]], rx );
               inv.push_back( hb.E( ) - 1, hb.E( ) - 2 );
               for ( int pass = 1; pass <= 2; pass++ )
               {    const vec<int>& c = ( pass == 1 ? echain : rechain );
                    for ( int j = 0; j < c.isize( ); j++ )
                    {    int e = c[j];
                         for ( int64_t u = 0; u < (int64_t) paths_index[e].size( ); u++ )
                         {    ReadPath& p = paths[ paths_index[e][u] ];
                              if ( p[0] != e ) continue;
                              int offset = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   offset += hb.Kmers( c[l] );
                              p.setOffset(offset);
                              p.resize(1);
                              if ( pass == 1 ) p[0] = hb.E( ) - 2;    
                              else p[0] = hb.E( ) - 1;    }    }    }    }    }
     hb.DeleteEdges(dels);
     CleanupCore( hb, inv, paths );    }

void RemoveSmallComponents3( HyperBasevector& hb, const Bool remove_small_cycles )
{    double clock1 = WallClockTime( );
     const int max_small_comp = 1000;
     const int min_circle = 200;
     vec<int> e_to_delete;
     vec< vec<int> > comps;
     hb.Components(comps);
     LogTime( clock1, "removing small components 1" );
     double clock2 = WallClockTime( );
     #pragma omp parallel for
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];

          int max_edge = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          for ( int l = 0; l < hb.From( o[j] ).isize( ); l++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( o[j], l );
               max_edge = std::max( max_edge, hb.EdgeLengthKmers(e) );    }
          if ( max_edge > max_small_comp ) continue;

          // Remove small cycles.

          int total_kmers = 0;
          for ( int j = 0; j < o.isize( ); j++ )
          for ( int l = 0; l < hb.From( o[j] ).isize( ); l++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( o[j], l );
               total_kmers += hb.EdgeLengthKmers(e);    }
          if ( total_kmers < min_circle && remove_small_cycles )
          {    for ( size_t j = 0; j < o.size( ); j++ )
               {    int v = o[j];
                    #pragma omp critical
                    {    for ( size_t t = 0; t < hb.From(v).size( ); t++ )
                         {    e_to_delete.push_back( hb.EdgeObjectIndexByIndexFrom( 
                                   v, t ) );    }    }    }
               continue;    }
          if ( hb.HasCycle(o) ) continue;

          int no = o.size( );
          digraphE<basevector> GX( digraphE<basevector>::COMPLETE_SUBGRAPH, hb, o );
          vec<int> L, sources, sinks, p;
          for ( int i = 0; i < GX.EdgeObjectCount( ); i++ )
               L.push_back( GX.EdgeObject(i).size( ) - hb.K( ) + 1 );
          digraphE<int> G( GX, L );
          G.Sources(sources), G.Sinks(sinks);
          G.AddVertices(2);
          for ( int j = 0; j < sources.isize( ); j++ )
               G.AddEdge( no, sources[j], 0 );
          for ( int j = 0; j < sinks.isize( ); j++ )
               G.AddEdge( sinks[j], no+1, 0 );
          for ( int e = 0; e < G.EdgeObjectCount( ); e++ )
               G.EdgeObjectMutable(e) = -G.EdgeObject(e);
          G.ShortestPath( no, no+1, p );
          int max_path = 0;
          for ( int j = 0; j < p.isize( ) - 1; j++ )
          {    int v1 = p[j], v2 = p[j+1];
               int m = 1000000000;
               for ( int l = 0; l < G.From(v1).isize( ); l++ )
               {    if ( G.From(v1)[l] == v2 )
                         m = Min( m, G.EdgeObjectByIndexFrom( v1, l ) );    }
               max_path -= m;    }

          if ( max_path <= max_small_comp )
          {    for ( size_t j = 0; j < o.size( ); j++ )
               {    int v = o[j];
                    #pragma omp critical
                    {    for ( size_t t = 0; t < hb.From(v).size( ); t++ )
                         {    e_to_delete.push_back( hb.EdgeObjectIndexByIndexFrom( 
                                   v, t ) );    }    }    }    }    }
     LogTime( clock2, "removing small components 2" );
     double clock3 = WallClockTime( );
     hb.DeleteEdges(e_to_delete);
     LogTime( clock3, "removing small components 3" );    }

void Empty( const HyperBasevector& hb, const vec< std::pair<vec<int>,vec<int>> >& pairs, 
     const vec<int64_t>& pairs_pid, vec<vec<int>>& left_empty, 
     vec<vec<int>>& right_empty, const Bool EMPTY2 )
{
     int nedges = hb.EdgeObjectCount( );
     left_empty.resize(nedges), right_empty.resize(nedges);
     if (EMPTY2)
     {    for ( int l = 0; l < pairs.isize( ); l++ )
          {    if ( pairs[l].first.nonempty( ) && pairs[l].second.empty( ) )
               {    left_empty[ pairs[l].first.back( ) ].push_back( pairs_pid[l] );
                    if ( pairs[l].first.front( ) != pairs[l].first.back( ) )
                    {    left_empty[ pairs[l].first.front( ) ].
                              push_back( pairs_pid[l] );    }    }
               if ( pairs[l].second.nonempty( ) && pairs[l].first.empty( ) )
               {    right_empty[ pairs[l].second.front( ) ]
                         .push_back( pairs_pid[l] );    
                    if ( pairs[l].second.front( ) != pairs[l].second.back( ) )
                    {    right_empty[ pairs[l].second.back( ) ]
                              .push_back( pairs_pid[l] );    }    }    }    }
     else
     {    for ( int l = 0; l < pairs.isize( ); l++ )
          {    if ( pairs[l].first.nonempty( ) && pairs[l].second.empty( ) )
                    left_empty[ pairs[l].first.back( ) ].push_back( pairs_pid[l] );
               if ( pairs[l].second.nonempty( ) && pairs[l].first.empty( ) )
               {    right_empty[ pairs[l].second.front( ) ]
                         .push_back( pairs_pid[l] );    }    }    }    }

void Validate( const HyperBasevector& hb, const vec<int>& inv, 
	       const ReadPathVec& paths ) {
    if (ValidateAllReadPaths(hb, paths) == false)
	TracebackThisProcess();
}


void TestIndex( const HyperBasevector& hb,
        const ReadPathVec& paths, const VecULongVec& invPaths)
{
    // horribly inefficient, just meant to be quick to code
    bool good = true;

    // is each edge entry reflected in a path?
    for ( size_t edge = 0; edge < invPaths.size(); ++edge )
        for ( auto const readid : invPaths[edge] )
            if ( std::find( paths[readid].begin(),
                    paths[readid].end(), edge ) == paths[readid].end() ) {
                std::cout << "the index for edge " << edge << " names readid "
                        << readid << " but the readpath says "
                        << printSeq( paths[readid] ) << std::endl;
                good = false;
            }

    // is each path entry reflected in an edge index
    for ( int64_t pathid = 0; pathid < paths.size(); ++pathid )
        for ( auto const edge : paths[pathid] )
            if ( std::find( invPaths[edge].begin(),
                    invPaths[edge].end(), pathid ) == invPaths[edge].end() ) {
                std::cout << "the path for read " << pathid << " names edge " <<
                edge << " but the index for that edge says " <<
                printSeq(invPaths[edge]) << std::endl;
                good = false;
            }

    if ( !good ) TracebackThisProcess();
}


void TestInvolution( const HyperBasevector& hb, const vec<int>& inv )
{
     time_t now = time(0);
     vec<int> to_left, to_right;
     vec<Bool> used;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     hb.Used(used);
     if ( inv.isize( ) != hb.EdgeObjectCount( ) )
     {    std::cout << "\n";
          PRINT2( hb.EdgeObjectCount( ), inv.size( ) );
          std::cout << "Involution has wrong size.\n" << "Abort." << std::endl;
          TracebackThisProcess( );    }
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    if ( !used[e] ) continue;
          if ( inv[e] < 0 || inv[e] >= hb.EdgeObjectCount( ) )
          {    std::cout << "\n";
               PRINT3( e, inv[e], hb.EdgeObjectCount( ) );
               std::cout << "Illegal involution value.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }
          basevector b = hb.EdgeObject(e);
          b.ReverseComplement( );
          if ( b != hb.EdgeObject( inv[e] ) )
          {    std::cout << "\n";
               int re = inv[e];
               PRINT4( e, re, b.size( ), hb.EdgeObject(re).size( ) );
               std::cout << "Involution value not rc.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }
          if ( inv[inv[e]] != e )
          {    std::cout << "Involution is not an involution.\n" << "Abort." << std::endl;
               TracebackThisProcess( );    }    }
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i1 = 0; i1 < hb.To(v).isize( ); i1++ )
          for ( int i2 = 0; i2 < hb.From(v).isize( ); i2++ )
          {    int e1 = hb.EdgeObjectIndexByIndexTo( v, i1 );
               int e2 = hb.EdgeObjectIndexByIndexFrom( v, i2 );
               int re1 = inv[e1], re2 = inv[e2];
               if ( to_right[re2] != to_left[re1] )
               {    std::cout << "Involution does not preserve graph structure.\n";
                    PRINT3( e1, to_left[e1], to_right[e1] );
                    PRINT3( e2, to_left[e2], to_right[e2] );
                    PRINT3( re1, to_left[re1], to_right[re1] );
                    PRINT3( re2, to_left[re2], to_right[re2] );
                    std::cout << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }
     // Note sure the following can happen:
     for ( int v = 0; v < hb.N( ); v++ )
     {    for ( int i = 0; i < hb.From(v).isize( ); i++ )
          {    int w = hb.From(v)[i];
               int e = hb.EdgeObjectIndexByIndexFrom( v, i );
               int re = inv[e];
               int rv = to_left[re], rw = to_right[re];
               if ( hb.From(rv).size( ) != hb.To(w).size( )
                    || hb.To(rw).size( ) != hb.From(v).size( ) )
               {    std::cout << "Graph structure is asymmetric.\n" << "Abort." << std::endl;
                    TracebackThisProcess( );    }    }    }     
     //std::cout << "[GapToyTools3.cc] Finishing test involution : " << ctime(&now) << std::endl;
     }

void FragDist( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, const String out_file )
{    
     // Generate fragment distribution.

     const int width = 10;
     const int max_sep = 1000;
     const int min_edge = 10000;
     vec<double> count( max_sep/width, 0 );
     for ( int64_t id1 = 0; id1 < (int64_t) paths.size( ); id1 += 2 )
     {    int64_t id2 = id1 + 1;
          if ( paths[id1].size( ) == 0 || paths[id2].size( ) == 0 ) continue;
          int e1 = paths[id1][0], e2 = inv[ paths[id2][0] ];
          int epos1 = paths[id1].getOffset( );
          if ( e1 != e2 ) continue;
          if ( hb.EdgeLengthBases(e1) < min_edge ) continue;
          int epos2 = hb.EdgeLengthBases(e2) - paths[id2].getOffset( );
          int len = epos2 - epos1;
          if ( len < 0 || len >= max_sep ) continue;
          count[ len/width ]++;    }

     // Output distribution.

     double total = Sum(count);
     {    Ofstream( out, out_file );
          out << "# fragment library size distribution" << std::endl;
          out << "# bins have diameter 10" << std::endl;
          out << "# line format:\n";
          out << "# bin_center mass" << std::endl;
          for ( int j = 0; j < count.isize( ); j++ )
               out << j * width + (width/2) << " " << count[j]/total << std::endl;    }

     // Check for abject failure.

     if ( total == 0 )
     {    Ofstream( out, out_file + ".png.FAIL" );
          Remove( out_file + ".png" );
          out << "Could not generate frags.dist.png because there was not\n"
               << "enough assembly to compute the distribution." << std::endl;
          return;    }

     // Check for missing executables.

     vec<String> missing;
     vec<String> ex = { "ps2epsi", "pstopnm", "pnmcrop", "pnmpad", "pnmtopng" };
     for ( auto executable : ex )
     {    if ( System( executable + " --version > /dev/null 2>&1" ) != 0 )
               missing.push_back(executable);    }
     if ( missing.nonempty( ) )
     {    Ofstream( out, out_file + ".png.FAIL" );
          out << "Could not generate frags.dist.png because the following "
               << "executables were not found:\n" << printSeq(missing) 
               << "." << std::endl;
          Remove( out_file + ".png" );
          return;    }

     // Make plot.
          
     String TITLE = "Fragment library size distribution";
     vec<graphics_primitive> points, lines;
     lines.push_back( SetLineWidth(1.0) );
     double xm = 0.0, ym = 0.0;
     for ( int j = 0; j < count.isize( ); j++ )
          points.push_back( Point( j * width + (width/2),  count[j]/total, 1.0 ) );
     points.push_back( SetColor(black) );
     points.append(lines);
     double x = 0, X = MaxX(points), y = 0, Y = MaxY(points);
     vec<graphics_primitive> x_axis = AxisX( x, X, 1.0, True, "", 0.0, False );
     x_axis.push_front( SetLineWidth(1.0) );
     vec<graphics_primitive> y_axis = AxisY( x, X, y, Y, True, "", 0.0, False );
     vec<graphics_primitive> points2(y_axis);
     points2.append(points);
     points = points2;
     vec< vec<graphics_primitive> > stack;
     vec<double> heights;
     stack.push_back(x_axis),    heights.push_back(0);
     stack.push_back(points),    heights.push_back(200);
     vec<graphics_primitive> title;
     title.push_back( TextCenter( TITLE, (X+x)/2.0, 0, 0, TimesBold(15) ) );
     stack.push_back(title);
     heights.push_back(20);
     String fail_msg = RenderGraphics( out_file + ".png", stack, heights, 1.0, 200, 
          50, True, 1.0, True );
     Remove( out_file + ".eps" );
     if ( fail_msg != "" )
     {    Remove( out_file + ".png" );
          Ofstream( out, out_file + ".png.FAIL" );
          out << "Could not generate frags.dist.png because something went "
               << "wrong, see below:\n\n" << fail_msg;
          return;    }    }

// UnwindThreeEdgePlasmids
//
// 1. Component has two vertices v, w.
//
// 2. Three edges: two from v to w (e1, e2), one from w to v (f)
//    (and if not, swap v and w).
//
// 3. e1, e2 well-supported (at least 10 pids each).
//
// 4. coverage of e1 and e2 within 25% of each other.
//
// Then replace by e1,f,e2,f loop.
//
// Note that this drops read placements that become ambiguous.
//
// Note that this will make a mistake if there are actually two plasmids present
// at nearly equal molarity.

void UnwindThreeEdgePlasmids(HyperBasevector& hb, vec<int>& inv, ReadPathVec& paths)
{
     // Heuristics.

     const int min_cov = 10;
     const double fudge = 0.25;
     const int min_links = 2;

     // Set up indices.

     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vec<int> to_right;
     hb.ToRight(to_right);

     // Look for the special components.

     vec< vec<int> > comps;
     hb.Components(comps);
     vec<int> dels;
     for ( size_t i = 0; i < comps.size( ); i++ )
     {    const vec<int>& o = comps[i];
          if ( o.size( ) != 2 ) continue;
          int v = o[0], w = o[1];
          // PRINT2( v, w ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( hb.From(v).size( ) != 2 ) std::swap( v, w );
          if ( hb.From(v).size( ) != 2 || hb.From(w).size( ) != 1 ) continue;
          if ( hb.From(v)[0] != w || hb.From(v)[1] != w || hb.From(w)[0] != v )
               continue;
          int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 ), f = hb.IFrom( w, 0 );
          // PRINT3( e1, e2, f ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          int re1 = inv[e1], re2 = inv[e2], rf = inv[f];
          vec<int> v1 = {e1,e2,f}, v2 = {re1,re2,rf};
          if ( Meet2( v1, v2 ) || Min(v2) < Min(v1) ) continue;

          // Make sure e1 and e2 are linked.

          vec<int64_t> pids;
          for ( int64_t j = 0; j < (int64_t) paths_index[e1].size( ); j++ )
               pids.push_back( paths_index[e1][j]/2 );
          for ( int64_t j = 0; j < (int64_t) paths_index[e2].size( ); j++ )
               pids.push_back( paths_index[e2][j]/2 );
          for ( int64_t j = 0; j < (int64_t) paths_index[f].size( ); j++ )
               pids.push_back( paths_index[f][j]/2 );
          for ( int64_t j = 0; j < (int64_t) paths_index[re1].size( ); j++ )
               pids.push_back( paths_index[re1][j]/2 );
          for ( int64_t j = 0; j < (int64_t) paths_index[re2].size( ); j++ )
               pids.push_back( paths_index[re2][j]/2 );
          for ( int64_t j = 0; j < (int64_t) paths_index[rf].size( ); j++ )
               pids.push_back( paths_index[rf][j]/2 );
          UniqueSort(pids);
          int links = 0;
          for ( int64_t j = 0; j < pids.isize( ); j++ )
          {    int64_t id1 = 2*pids[j], id2 = 2*pids[j] + 1;
               vec<int> es;
               for ( int64_t l = 0; l < (int64_t) paths[id1].size( ); l++ )
                    es.push_back( paths[id1][l], inv[ paths[id1][l] ] );
               for ( int64_t l = 0; l < (int64_t) paths[id2].size( ); l++ )
                    es.push_back( paths[id2][l], inv[ paths[id2][l] ] );
               UniqueSort(es);
               if ( BinMember( es, e1 ) && BinMember( es, e2 ) ) links++;    }
          // PRINT(links); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( links < min_links ) continue;

          // Make sure that e1 and e2 have enough support.  Otherwise we would 
          // worry that one of them is an error branch.

          int ne1 = Pids( e1, hb, inv, paths, paths_index ).size( );
          int ne2 = Pids( e2, hb, inv, paths, paths_index ).size( );
          if ( ne1 < min_cov || ne2 < min_cov ) continue;

          // Calculate a bunch of interesting stuff that we're not actually using.

          /*
          int nf = Pids( f, hb, inv, paths, paths_index ).size( ); // XXXXXXXXXXXXXX
          PRINT3( ne1, ne2, nf ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          PRINT3( hb.Kmers(e1), hb.Kmers(e2), hb.Kmers(f) ); // XXXXXXXXXXXXXXXXXXXX
          double ce1 = double(ne1) / ( hb.Kmers(e1) + hb.K( ) - 1 - 60 );
          double ce2 = double(ne2) / ( hb.Kmers(e2) + hb.K( ) - 1 - 60 );
          double r = Max(ce1,ce2)/Min(ce1,ce2);
          PRINT3( ce1, ce2, r ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          // double ce1 = RawCoverage( e1, hb, inv, paths, paths_index );
          // double ce2 = RawCoverage( e2, hb, inv, paths, paths_index );
          double cf = RawCoverage( f, hb, inv, paths, paths_index );
          if ( Max(ce1,ce2)/Min(ce1,ce2) - 1 > fudge ) continue;
          double r1 = (cf/ce1)/2.0, r2 = (cf/ce2)/2.0;
          PRINT2( r1, r2 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          double low = 1 - fudge, high = 1 + fudge;
          */

          {    
               // Start assembly edit.

               // std::cout << "editing" << std::endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               vec<int> x = {e1,f,e2,f}, rx = {rf,re2,rf,re1};
               basevector b = hb.Cat(x), rb = hb.Cat(rx);
               dels.push_back( e1, e2, f, re1, re2, rf );
               int m = hb.AddEdge( v, v, b );
               int rv = to_right[re1];
               int rm = hb.AddEdge( rv, rv, rb );    
               inv.push_back( rm, m );

               // Get the relevant read pairs.

               vec<int64_t> pids;
               for ( int e : x ) 
                    pids.append( Pids( e, hb, inv, paths, paths_index ) );
               UniqueSort(pids);

               // Remap the reads.  It would be smarter to take into account
               // pairing.

               vec<int64_t> ids;
               for ( int j = 0; j < pids.isize( ); j++ )
                    ids.push_back( 2*pids[j], 2*pids[j] + 1 );
               for ( int j = 0; j < ids.isize( ); j++ )
               {    ReadPath& p = paths[ ids[j] ];
                    if ( p.size( ) == 0 ) continue;
                    Bool fixed = False;
                    for ( int l = 0; l < (int) p.size( ); l++ )
                    {    int x = p[l];
                         int pre = 0;
                         for ( int r = 0; r < l; r++ )
                              pre += hb.Kmers( p[r] );
                         fixed = True;
                         if ( x == e1 )
                         {    p[0] = m;
                              p.addOffset(-pre);    }
                         else if ( x == e2 )
                         {    p[0] = m;
                              p.addOffset( -pre + hb.Kmers(e1) + hb.Kmers(f) );    }
                         else if ( x == re2 )
                         {    p[0] = rm;
                              p.addOffset( -pre + hb.Kmers(f) );
                              fixed = True;    }
                         else if ( x == re1 )
                         {    p[0] = rm;
                              p.addOffset( -pre + 2*hb.Kmers(f) + hb.Kmers(e2) );   }
                         else fixed = False;
                         if (fixed) 
                         {    p.resize(1);
                              break;    }    }
                    if ( !fixed ) p.resize(0);    }    }    }

     // Finish assembly edit.

     hb.DeleteEdges(dels);
     CleanupCore( hb, inv, paths );    }
