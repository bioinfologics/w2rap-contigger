///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "VecUtilities.h"
#include "efasta/AmbiguityScore.h"
#include "efasta/EfastaTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/ProcessGap.h"
#include "paths/LongReadTools.h"
#include "paths/Ulink.h"
#include "paths/Uniseq.h"
#include "paths/UnipathScaffold.h"

void GetWalks0( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc,
     const vec< vec< std::pair<int,int> > >& nextsx, vec<int> use,
     vec< vec< std::pair<int,int> > >& walks1, int& bad, const double max_dev_diff )
{    
     // Heuristics.

     const int max_devs = 10;
     const int max_iterations = 10000000;
     const int max_excess_charm = K;

     // Do two passes.  In the first pass, for now just exploratory, allow a given
     // vertex to be used at most twice.
     // FOR NOW TURNED OFF AS WE START WITH PASS 2.

     vec<double> offbys;
     for ( int pass = 1; pass <= 2; pass++ )
     {    
          double max_devs_this = max_devs;
          if ( pass == 2 )
          {    if ( bad == 0 )
               {    use.clear( );
                    for ( int i = 0; i < walks1.isize( ); i++ )
                    {    for ( int j = 0; j < walks1[i].isize( ); j++ )
                              use.push_back( walks1[i][j].first );    }
                    UniqueSort(use);    }    }

          walks1.clear( ), offbys.clear( );
          vec< vec< std::pair<int,int> > > walks0;
          vec<int> charm0, charm1; // misnomer!
          vec< std::pair<int,int> > w;
          w.push( u1, 0 );
          walks0.push_back(w);
          charm0.push_back(0);
          int best_charm = 1000000000;
          bad = 0;
          int iterations = 0;
          while( walks0.nonempty( ) )
          {    if ( ++iterations > max_iterations )
               {    bad = 1;
                    break;    }
               vec< std::pair<int,int> > w = walks0.back( );
               walks0.pop_back( );
               int ch = charm0.back( );
               charm0.pop_back( );
               int v1 = w.back( ).first;
               int sepx = 0;
               for ( int j = 1; j < w.isize( ) - 1; j++ )
                    sepx += unibases[ w[j].first ].isize( ) - w[j].second;
               sepx -= w.back( ).second;
               double offby = double(sepx-sep) / double(dev);
               if ( iterations > 1 && offby > max_devs_this ) continue;
               if ( v1 == u2 && w.size( ) > 1 )
               {    if ( offby <= max_devs_this && offby >= -max_devs && 
                         ch <= best_charm + max_excess_charm )
                    {    walks1.push_back(w);
                         charm1.push_back(ch);
                         best_charm = Min( best_charm, ch );
                         offbys.push_back( Abs(offby) );    
                         max_devs_this 
                              = Min( max_devs_this, Abs(offby) + max_dev_diff );    }
                    // continue;    
                         }
               for ( int i = 0; i < nextsx[v1].isize( ); i++ )
               {    int v2 = nextsx[v1][i].first;
                    if ( use.nonempty( ) && !BinMember( use, v2 ) ) continue;
                    if ( pass == 1 )
                    {    int count = 0;
                         for ( int j = 0; j < w.isize( ); j++ )
                              if ( w[j].first == v2 ) count++;
                         if ( count > 1 ) continue;    }
                    int over = nextsx[v1][i].second;
                    int ch_new = ch + K - 1 - over;
                    if ( ch_new <= best_charm + max_excess_charm )
                    {    vec< std::pair<int,int> > wx(w);
                         wx.push( v2, over );
                         walks0.push_back(wx);    
                         charm0.push_back(ch_new);    }    }    }
          vec<Bool> to_delete( walks1.size( ), False );
          for ( int i = 0; i < walks1.isize( ); i++ )
               if ( charm1[i] > best_charm + max_excess_charm ) to_delete[i] = True;
          EraseIf( walks1, to_delete ), EraseIf( offbys, to_delete );    }

     // Select passing walks.

     SortSync( offbys, walks1 );
     if ( walks1.nonempty( ) )
     {    int x;
          for ( x = 1; x < walks1.isize( ); x++ )
               if ( offbys[x] > offbys[0] + max_dev_diff ) break;
          walks1.resize(x);    }    }

void GetWalks( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc,
     const vec< vec< std::pair<int,int> > >& nextsx, const vec<int>& use,
     vec< vec< std::pair<int,int> > >& walks1, int& bad, const double max_dev_diff )
{
     GetWalks0( u1, u2, sep, dev, unibases, K, to_rc, nextsx, use, walks1, bad,
          max_dev_diff );
     if ( bad == 0 ) return;
     vec<int> use_rc;
     for ( int j = 0; j < use.isize( ); j++ )
          use_rc.push_back( to_rc[ use[j] ] );
     Sort(use_rc);
     GetWalks0( to_rc[u2], to_rc[u1], sep, dev, unibases, K, to_rc, nextsx, use_rc,
          walks1, bad, max_dev_diff );
     {    for ( int j = 0; j < walks1.isize( ); j++ )
          {    walks1[j].ReverseMe( );
               for ( int l = 0; l < walks1[j].isize( ); l++ )
                    walks1[j][l].first = to_rc[ walks1[j][l].first ];
               for ( int l = walks1[j].isize( ) - 1; l >= 1; l-- )
                    walks1[j][l].second = walks1[j][l-1].second;
               walks1[j][0].second = 0;     }    }    }

