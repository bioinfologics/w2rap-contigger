///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Apr 2, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef PULLAPARTER_H_
#define PULLAPARTER_H_
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

#include "paths/long/large/GapToyTools.h"

#include "Vec.h"

class PullAparter {
     public:
          PullAparter( HyperBasevector& hbv, vec<int>& inv,
                    ReadPathVec& paths, VecULongVec& invPaths,
		    vec<int> trace_edges = vec<int>{},
                    bool debug = false,
                    int min_reads = 5, float min_mult = 5.0) :
                    mHBV(hbv), mInv(inv), mPaths(paths), mEdgeToPathIds(invPaths),
                    mTraceEdges(trace_edges), mDebug(debug), mMinReads(min_reads),
                    mMinMult(min_mult), mRemovedReadPaths(0u)
                    { hbv.ToLeft(mToLeft); hbv.ToRight(mToRight); }

          size_t getRemovedReadPaths() const { return mRemovedReadPaths; }

          bool isCanonicalRepeatEdge(int edge) const {
               ForceAssertGe(edge,0);
               ForceAssertLt(edge,mHBV.EdgeObjectCount());
               int vleft = mToLeft[edge];
               int vright = mToRight[edge];

               return ( mHBV.FromSize(vleft) == 1 && mHBV.ToSize(vleft) == 2 &&
                        mHBV.ToSize(vright) == 1 && mHBV.FromSize(vright) == 2 &&
                        mHBV.To(vleft)[0] != vright && mHBV.To(vleft)[1] != vright );

          }

          template <typename TVec>
          TVec inversePath( TVec const& p ) const {
               TVec rp;

               for ( auto itr = p.crbegin(); itr != p.crend(); ++itr )
                    rp.push_back( mInv[*itr] );

               return rp;
          }

          void scorePathSupportEnds( vec<vec<int>> const& paths,
                    vec<int>& scores, VecULongVec* readids_p = nullptr) const {

               vec<vec<int>> rpaths;
               vec<int> ends;

               for ( auto const& p : paths ) {
                    ForceAssertEq(p.size(),3u);
                    rpaths.push_back( inversePath(p) );
                    ends.push_back(p[0]);
                    ends.push_back(p[2]);
                    ends.push_back(rpaths.back()[0]);
                    ends.push_back(rpaths.back()[2]);
               }
               UniqueSort(ends);

               vec<size_t> readPathIndices;
               for ( auto const end : ends ) {
                       for ( auto const p : mEdgeToPathIds[end] ) {
                           // ensure that we're pushing on pairs
                           int offset = ( p % 2 == 0 ) ? +1 : -1 ;
                           readPathIndices.push_back(p);
                           readPathIndices.push_back(p+offset);
                       }
               }
               UniqueSort(readPathIndices);

               if ( mDebug ) std::cout << "scoring " << readPathIndices.size()
                             << " reads" << std::endl;

               scores.resize_and_set( paths.size(), 0 );
               if ( readids_p ) {
                   readids_p->clear();
                   readids_p->resize(paths.size());
               }
               ForceAssertEq(readPathIndices.size()%2,0u);
               for ( auto itr = readPathIndices.cbegin();
                         itr != readPathIndices.cend(); advance(itr,2) ) {
                    std::vector<int> readPath = mPaths[itr[0]];
                    std::vector<int> readPath1 = inversePath(mPaths[itr[1]]);

                    OverlapAppend(readPath, readPath1);

                    for ( size_t ipath=0; ipath < paths.size(); ++ipath ) {
                         auto f0 = readPath.end();
                         auto f2 = f0, r0 = f0, r2 = f0;
                         for ( auto rpItr = readPath.begin(); rpItr != readPath.end(); ++rpItr ) {
                              if ( *rpItr == paths[ipath][0] ) f0 = rpItr;
                              else if ( *rpItr == paths[ipath][2] ) f2 = rpItr;
                              else if ( *rpItr == rpaths[ipath][0] ) r0 = rpItr;
                              else if ( *rpItr == rpaths[ipath][2] ) r2 = rpItr;
                         }

                         if ( ( f0 < f2 && f2 != readPath.end() ) ||
                              ( r0 < r2 && r2 != readPath.end() ) ) {
                              scores[ipath]++;
                              if ( readids_p ) {
                                  // maybe someday we could distinguish between paths seen
                                  // entirely on one read or another
                                  (*readids_p)[ipath].push_back(itr[0]);
                                  (*readids_p)[ipath].push_back(itr[1]);
                              }
                         }
                    }
            }

            if ( mDebug ) {
                 for ( size_t i = 0; i < paths.size(); ++i ) {
                      std::cout << "score=" << scores[i] << ", path=" << printSeq(paths[i])
                                                       << std::endl;
                 }
            }

          }

          void dumpReadPath(ReadPath const& rp )
          {
              for ( auto const edge : rp )  {
                  std::cout << "(" << mToLeft[edge] << ") " << edge << " ";
                  std::cout << "(" << mToRight[edge] << ")";
              }
              std::cout << std::endl;
          }

          void nukeReadPaths( ULongVec const& rp1 )
          {
              if ( mDebug ) {
                  std::cout << "removing the following readPaths:" << std::endl;
                  std::cout << printSeq(rp1) << std::endl;
              }
              for ( auto const readid : rp1 ) {
                  for ( auto const edge : mPaths[readid] ) {
                      auto& inv = mEdgeToPathIds[edge];
                      auto new_end = std::remove( inv.begin(), inv.end(), readid );
                      inv.erase( new_end, inv.end() );
                  }
                  mPaths[readid].clear();
              }
              mRemovedReadPaths += rp1.size();
          }


          bool isSeparable(int edge, vec<vec<int>>* sep_paths, bool nukeBadReadPaths = false) {

               if ( ! isCanonicalRepeatEdge(edge) ) return false;

               int vleft = mToLeft[edge];
               int vright = mToRight[edge];

               // first find the "side" edges
               std::pair<int,int> ledges,redges;
               ForceAssertEq(mHBV.ToSize(vleft),2);
               ledges = std::make_pair( mHBV.EdgeObjectIndexByIndexTo(vleft,0),
                                mHBV.EdgeObjectIndexByIndexTo(vleft,1)  );
               ForceAssertEq(mHBV.FromSize(vright),2);
               redges = std::make_pair(  mHBV.EdgeObjectIndexByIndexFrom(vright,0),
                                mHBV.EdgeObjectIndexByIndexFrom(vright,1) );

               // determine if the inversion of any edge is here
               vec<int> alledges = { edge, ledges.first, ledges.second, redges.first, redges.second };
               for ( int e : alledges )
                    if ( Member(alledges, mInv[e]) )
                         return false;      // EARLY EXIT


               // We're hard-coded for now for the canonical repeat
               // so we have four paths through.  Let's score them
               vec<vec<int>> paths(4);
               vec<int> scores;
               paths[0] = { ledges.first, edge, redges.first };
               paths[1] = { ledges.first, edge, redges.second };
               paths[2] = { ledges.second, edge, redges.first };
               paths[3] = { ledges.second, edge, redges.second };

               if ( mDebug ) std::cout << std::endl << std::endl << "scoring edge id#" << edge << std::endl;
               VecULongVec scoreReads;

               scorePathSupportEnds(paths, scores, &scoreReads );

               vec<int> index(scores.size(), vec<int>::IDENTITY);
               ReverseSortSync(scores,index);
               int sum1 = scores[0] + scores[1];

               if ( sum1 < mMinReads ||
                         sum1 < mMinMult*scores[2] ||
                         sum1 < mMinMult*scores[3] ) {
                    if ( mDebug ) {
                         std::cout << "scoring failed: " << std::endl;
                         std::cout << "scores:    " << printSeq(scores) << std::endl;
                         std::cout << "index:     " << printSeq(index) << std::endl;
                         std::cout << "sum1:      " << sum1 << std::endl;
                         std::cout << "mMinReads: " << mMinReads << std::endl;
                         std::cout << "mMinMult:  " << mMinMult << std::endl;
                    }
                    return false;
               }

               unsigned int pathMask = ( 1 << index[0] ) | ( 1 << index[1] );

               if ( mDebug ) std::cout << "pathMask=" << pathMask << std::endl;

               if ( pathMask == 0b1001 ) {
                    if ( sep_paths ) sep_paths->push_back( paths[0], paths[3] );
                    if ( nukeBadReadPaths ) {
                        nukeReadPaths( scoreReads[1] ); nukeReadPaths(scoreReads[2]);
                    }
               } else if ( pathMask == 0b0110 ) {
                    if ( sep_paths ) sep_paths->push_back( paths[1], paths[2] );
                    if ( nukeBadReadPaths ) {
                        nukeReadPaths( scoreReads[0] ); nukeReadPaths(scoreReads[3]);
                    }
               } else {
                    if ( mDebug ) {
                            std::cout << "failing because we found a 'cross' path" << std::endl;
                            std::cout << "pathMask: " << pathMask << std::endl;
                    }
                    return false;
               }

               return true;
          }

          // pull apart a pair of paths through a canonical repeat
          // path1 essentially gets a new center edge
          // path2 only used for sanity check at the moment
          //
          // returns the new center edge given to path1
          int Separate( vec<int> const& path1, vec<int> const& path2 )
          {
               //          \______/      XX
               //          /      \      XX

               // assumes pairs of paths
               ForceAssertEq(path1.size(),3u);
               ForceAssertEq(path2.size(),3u);
               ForceAssertEq(path1[1],path2[1]);
               ForceAssertNe(path1[0],path2[0]);
               ForceAssertNe(path1[2],path2[2]);

               // repeat vertices
               int center = path1[1];
               int v1 = mToLeft[center];
               int v2 = mToRight[center];

               int newV1 = mHBV.N();
               int newV2 = newV1+1;
               mHBV.AddVertices(2);

               // repeat center edge
               ForceAssertEq( mToLeft.isize(), mHBV.EdgeObjectCount() );
               ForceAssertEq( mToRight.isize(), mHBV.EdgeObjectCount() );

               const int newCenter = mHBV.AddEdge(newV1,newV2, mHBV.EdgeObject(center));
               if (mDebug) std::cout << "new center = " << newCenter << std::endl;

               mToLeft.push_back( newV1 );
               mToRight.push_back( newV2 );

               // transfer (one set of) side edges
               mHBV.GiveEdgeNewToVx(path1[0], v1, newV1);
               mHBV.GiveEdgeNewFromVx(path1[2], v2, newV2);
               mToRight[path1[0]] = newV1;
               mToLeft[path1[2]] = newV2;

              return newCenter;
          }

          size_t SeparateAll()
          {
               vec<vec<int>> to_separate;
               for ( int i = 0; i < mHBV.EdgeObjectCount(); ++i ) {
                    if ( i < mInv[i] ) {
                        // check if separable
                        bool sep = this->isSeparable(i, &to_separate, true);

                        // edge trace info
                        if ( mTraceEdges.size() &&
                                ( Member(mTraceEdges, i) ||
                                  Member(mTraceEdges, mInv[i]) ) ) {
                            std::cout << "testing edge " << i << ", inv " <<
                                    mInv[i] << ": ";
                            if ( sep ) std::cout << "separable";
                            else std::cout << "NOT separable";
                            std::cout << std::endl;
                        }
                    }
               }
               SeparateAll(to_separate);
               ForceAssertEq(to_separate.size()%2,0u);
               return to_separate.size() / 2;
          }


          // given pairs of paths through canonical repeats,
          // pull apart each of them AND their counterparts on the
          // inverse.
          void SeparateAll( vec<vec<int>>& paths )
          {
               ForceAssertEq(paths.size() % 2, 0u);
               mEdgeToPathIds.reserve(mEdgeToPathIds.size() + paths.size() );
               mInv.reserve(mInv.size() + paths.size() );

               double clock = WallClockTime();
               if ( mTraceEdges.size() ) std::cout << "TRACE EDGES: " << printSeq(mTraceEdges) << std::endl;
               for ( auto itr = paths.begin(); itr != paths.end(); advance(itr,2) ) {
                 auto inv0 = inversePath( itr[0] );
                 auto inv1 = inversePath( itr[1] );
                 int center = itr[0][1];
                 int invcenter = inv0[1];
                 int p1center = Separate(itr[0],itr[1]);
                 int p1centerInv = Separate(inv0, inv1);
                 if ( Member(mTraceEdges,center) ) {
                     std::cout << "trace edge " << center << " new center " << p1center << std::endl;
                 }
                 ForceAssertEq( mInv.isize(), p1center );
                 mInv.push_back( p1centerInv );        // inv of center
                 mInv.push_back( p1center );           // inv of centerInv
                 MigrateReadPaths( itr[0], itr[1], p1center );
                 MigrateReadPaths( inv0, inv1, p1centerInv );

                 itr[0][1] = p1center;  // update the first path with the new center edge
               }
               std::cout << TimeSince(clock) << " used separating paths 1" << std::endl;


               //////////////////////////////////////////////////////
               // at this point we have:             now we work on this:
               //
               //      o----o----o----o              o---------o
               //
               //      o----o----o----o              o---------o
               //
               // and are fully consistent
               // in all data structures.
               mHBV.ToLeft(mToLeft);
               mHBV.ToRight(mToRight);
               // std::cout << mHBV.EdgeObjectCount() <<  "/" << mHBV.N() <<
               //  " edges/vertices before removing unneeded vertices" << std::endl;
               RemoveUnneededVertices2(mHBV, mInv, mPaths);

               double clock2 = WallClockTime();

               // remove dead edge objects
               vec<int> renumber_edges = mHBV.RemoveDeadEdgeObjects();

               // fix mToLeft and mToRight
               mHBV.ToLeft(mToLeft);
               mHBV.ToRight(mToRight);

               // fix mInv - only due to dead object removal
               for ( int oldc = 0; oldc < renumber_edges.isize(); ++oldc ) {
                    // this works, because we'll only ever be looking forward
                    // in the lists.  The equality test is obviously not necessary.
                    // the invariant is:  newc <= oldc
                    int newc = renumber_edges[oldc];
                    ForceAssertLe(newc, oldc);
                    if ( newc != -1) mInv[newc] = mInv[oldc];
               }
               mInv.resize(mHBV.EdgeObjectCount());

               // Workaround for temporary possibility that a member of inv is -1.
               /*
               for (auto const& inv : mInv ) ForceAssertGe(inv,0);
               for (auto& inv : mInv ) inv = renumber_edges[inv];
               */
               for (auto& inv : mInv ) 
               {    if ( inv >= 0 ) inv = renumber_edges[inv];    }

               // fix ReadPaths and ReadPaths index
               // ReadPaths are broken due to removal of dead edge objects
               // ReadPaths index was broken by removing unneeded vertices
               mEdgeToPathIds.clear();
               mEdgeToPathIds.resize(mHBV.EdgeObjectCount());
               for ( size_t readid = 0; readid < mPaths.size(); ++readid ) {
                    ReadPath& path = mPaths[readid];
                    for ( auto& edge : path ) {
                        if ( renumber_edges[edge] != -1 ) {
                            edge = renumber_edges[edge];
                            mEdgeToPathIds[edge].push_back(readid);
                        } else {
                            std::cout << "WARNING: read " << readid <<
                                    " contains a dead edge " << edge << std::endl;
                            path.clear(); /// NIW MUST ALSO REMOVE FROM INDEX
                        }
                    }
               }
               std::cout << TimeSince(clock2) << " used in fixing mToLeft, mToRight, "
                         << "and mEdgeToPathIds" << std::endl;
               // std::cout << mHBV.EdgeObjectCount() <<  "/" << mHBV.N() <<
               //   " edges/vertices after removing unneeded vertices and dead edges"
               //           << std::endl;

          }


          void MigrateReadPaths( vec<int> const& path1, vec<int> const& path2,
                    int newP1Center )
          {
               bool ldebug = false;
               ForceAssertEq(path1.size(), 3u);
               ForceAssertEq(path2.size(), 3u);

               // use index to find all readpaths touching the old center
               int center = path1[1];
               auto invPaths = mEdgeToPathIds[center];
               ULongVec newInv, oldInv;
               newInv.reserve(invPaths.size());
               oldInv.reserve(invPaths.size());
               bool trace =
                       Member( mTraceEdges, center )
                       || Member( mTraceEdges, newP1Center ) ;

               if ( ldebug ) {
                       std::cout << "MigrateReadPaths: center=" << center << ", paths:" << std::endl;
                       std::cout << "path1: " << printSeq(path1) << std::endl;
                       std::cout << "path2: " << printSeq(path2) << std::endl;
                       std::cout << "associated readPaths: " << printSeq(invPaths) << std::endl;
               }

               for ( auto itr = invPaths.begin(); itr != invPaths.end(); ++itr ) {
                    if ( ldebug ) std::cout << "*itr = " << *itr << std::endl;
                    // if readpath contains path1[0] or path1[2], then update with
                    // the new center edge
                    auto& readPath = mPaths[*itr];
                    if ( readPath.size() == 0 ) {
                        if ( ldebug ) std::cout << "... readpath already nulled" << std::endl;
                        continue;       // may have been nulled earlier
                    }
                    std::vector<int> extReadPath = readPath;
                    int pair_offset = ( *itr % 2 == 0 ) ? +1 : -1 ;
                    std::vector<int> extReadPath1;
                    extReadPath1 = inversePath(mPaths[*itr + pair_offset]);

                    // make the read path contain the paired read's edges, too.
                    OverlapAppend(extReadPath, extReadPath1);

                    vec<unsigned long> path_ids = { *itr, (*itr + pair_offset) } ;

                    if ( trace ) std::cout << "readPath (" << path_ids[0] << "," << path_ids[1] << "): ";

                    bool path1_support =
                            std::find(extReadPath.begin(), extReadPath.end(), path1[0])
                                != extReadPath.end() ||
                            std::find(extReadPath.begin(), extReadPath.end(), path1[2])
                                != extReadPath.end();

                    bool path2_support =
                            std::find(extReadPath.begin(), extReadPath.end(), path2[0])
                                != extReadPath.end() ||
                            std::find(extReadPath.begin(), extReadPath.end(), path2[2])
                                != extReadPath.end();

                   if ( ldebug ) { std::cout << "..."; PRINT2(path1_support, path2_support); }

                    if ( path1_support && !path2_support ) {
                         if ( ldebug || trace ) {
                             std::cout << printSeq(readPath) << "," << printSeq(extReadPath1)<< std::endl;
                             std::cout << "   gets new center " << newP1Center << std::endl;
                             std::cout << "   extended read path " << printSeq(extReadPath) << std::endl;
                         }
                         // update readPath with new center edge
                         std::replace( readPath.begin(), readPath.end(), center, newP1Center );
                         // add to new inverse list for new center edge
                         newInv.push_back(*itr);
                    } else if ( path2_support && !path1_support ) {
                         // add to new inverse list for old center edge
                         oldInv.push_back(*itr);
                         if ( ldebug || trace ) {
                             std::cout << printSeq(readPath) << std::endl;
                             std::cout << "   KEEPS OLD center " << std::endl;
                         }
                    } else {
                         // if there's no evidence linking it to an adjacent edge then
                         // we just remove the path.  At this point it better just be a
                         // single edge path or the path is no longer valid in the graph.
                         if ( ldebug || trace ) {
                             std::cout << printSeq(readPath)  << std::endl;
                             std::cout << "   has no (or conflicting) evidence " << std::endl;
                         }

                         for ( auto const edge : mPaths[*itr] ) {
                             auto& inv = mEdgeToPathIds[edge];
                             auto new_end = std::remove( inv.begin(), inv.end(), *itr );
                             inv.erase( new_end, inv.end() );
                         }
                         mPaths[*itr].clear();
                         mRemovedReadPaths++;
                         if ( mDebug && path1_support && path2_support ) {
                             std::cout << "WARNING: CONFLICTING ReadPaths:" << std::endl;
                             std::cout << "readPath1: " << printSeq(readPath) << std::endl;
                             std::cout << "inv(readPath2): " << printSeq(extReadPath1) << std::endl;
                             std::cout << "graph path1: " << printSeq(path1) << std::endl;
                             std::cout << "graph path2: " << printSeq(path2) << std::endl;
                             std::cout << "newP1Center: " << newP1Center << std::endl;
                         }

                    }
                    if (ldebug) std::cout << "XXX final *itr=" << *itr
                        << ", mPaths[*itr]=" << printSeq(mPaths[*itr])
                        << " // " << printSeq(readPath) << std::endl;
               }
               mEdgeToPathIds[center] = oldInv;
               ForceAssertEq( mEdgeToPathIds.size(), static_cast<unsigned>(newP1Center) );
               mEdgeToPathIds.push_back( newInv );
          }


     private:
          HyperBasevector& mHBV;
          vec<int>& mInv;
          ReadPathVec& mPaths;
          VecULongVec& mEdgeToPathIds;
          vec<int> mToLeft;
          vec<int> mToRight;
          vec<int> mTraceEdges;
          bool mDebug;
          int mMinReads;
          float mMinMult;
          size_t mRemovedReadPaths;
     };



#endif /* PULLAPARTER_H_ */
