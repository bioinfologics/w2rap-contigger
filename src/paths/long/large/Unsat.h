#ifndef UNSAT_H
#define UNSAT_H

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

void Unsat2( const HyperBasevector& hb, const vec<int>& inv,
     const ReadPathVec& paths, vec< vec< std::pair<int,int> > >& xs,
     const String& work_dir, const int A2V );

void Unsat( const HyperBasevector& hb, const vec<int>& inv,
            const ReadPathVec& paths, vec< vec< std::pair<int,int> > >& xs,
            const String& work_dir, const int A2V );

void MergeClusters( const vec< vec< std::pair<int,int> > >& x,
                    vec< vec< std::pair<int,int> > >& y, const vec< vec<int> >& n, const int N );


void MergeClusters2( const vec< vec< std::pair<int,int> > >& x,
                    vec< vec< std::pair<int,int> > >& y, const vec< vec<int> >& n, const int N );

void MergeClusters3( const vec< vec< std::pair<int,int> > >& x,
                     vec< vec< std::pair<int,int> > >& y, const vec< vec<int> >& n, const int N );
#endif
