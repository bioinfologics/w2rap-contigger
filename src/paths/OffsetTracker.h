/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_OFFSETTRACKER_H
#define PATHS_OFFSETTRACKER_H

#include "Vec.h"

#include "paths/UnipathSeq.h"
#include "paths/MuxGraph.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "Set.h"

// Offset where both superseqs are explicit.
class ExplicitOffset {
  public:
    ExplicitOffset() {}

    ExplicitOffset( int from, int to, int amount, int source = 0 )
        : m_from( from ), m_to( to ), m_amt( amount ), m_source( source ) {
    }

    int GetFrom() const {
        return m_from;
    }
    int GetTo() const {
        return m_to;
    }
    int GetAmount() const {
        return m_amt;
    }
    int GetSource() const {
        return m_source;
    }

    bool operator< ( const ExplicitOffset& other ) const {
        return ( m_to < other.m_to ||
                 m_to == other.m_to && ( m_from < other.m_from ||
                                         m_from == other.m_from && m_amt < other.m_amt ) );
    }

    bool operator== ( const ExplicitOffset& other ) const {
        return ( m_from == other.m_from &&
                 m_to == other.m_to &&
                 m_amt == other.m_amt );
    }

    struct OrderByDecreasingAmount : public std::binary_function<ExplicitOffset,ExplicitOffset,bool> {
        bool operator() ( const ExplicitOffset& lhs, const ExplicitOffset& rhs ) const {
            return ( lhs.GetAmount() > rhs.GetAmount() );
        }
    };

    struct OrderByIncreasingAmount : public std::binary_function<ExplicitOffset,ExplicitOffset,bool> {
        bool operator() ( const ExplicitOffset& lhs, const ExplicitOffset& rhs ) const {
            return ( lhs.GetAmount() < rhs.GetAmount() );
        }
    };

  private:
    // TODO: potentially dangerous truncation of index
    int m_from;
    int m_to;
    int m_amt;
    int m_source;
};


// Offset where "to" superseq is implicit.
class ImplicitOffset {
  public:
    ImplicitOffset() {}
    ImplicitOffset( int from, int amount )
        : m_from( from ), m_amt( amount ) {}

    int GetFrom() const {
        return m_from;
    }
    int GetAmount() const {
        return m_amt;
    }

    bool operator< ( const ImplicitOffset& other ) const {
        return ( m_from < other.m_from ||
                 m_from == other.m_from && m_amt < other.m_amt );
    }

    struct OrderByFrom : public std::binary_function<ImplicitOffset,ImplicitOffset,bool> {
        bool operator() ( const ImplicitOffset& lhs, const ImplicitOffset& rhs ) const {
            return ( lhs.GetFrom() < rhs.GetFrom() );
        }
    };

  private:
    // TODO: potentially dangerous truncation of index
    int m_from;
    int m_amt;
};

TRIVIALLY_SERIALIZABLE(ImplicitOffset);
typedef SerfVec<ImplicitOffset> ImplicitOffsetVec;
typedef MasterVec<ImplicitOffsetVec> VecImplicitOffsetVec;
extern template class SmallVec< ImplicitOffset, MempoolAllocator<ImplicitOffset> >;
extern template class OuterVec<ImplicitOffsetVec>;

class MutableOffsetTracker {
  public:
    MutableOffsetTracker() {}

    MutableOffsetTracker( const int size )
        : m_offsets( size ) {}

    MutableOffsetTracker( const vecUnipathSeq& unipathSeqs,
                          const MuxGraph& inverseMuxGraph,
                          const int firstSuperSeq,
                          const vec<int>& maxOffsets );

    int size() const {
        return m_offsets.size();
    }

    void resize( const int size ) {
        m_offsets.resize( size );
    }

    bool Add( int from, int to, int offsetAmount ) {
        std::pair<std::set<ImplicitOffset>::iterator,bool> result =
            m_offsets[ to ].insert( ImplicitOffset( from, offsetAmount ) );
        return result.second;
    }

    bool Add( const ExplicitOffset& offset ) {
        std::pair<std::set<ImplicitOffset>::iterator,bool> result =
            m_offsets[ offset.GetTo() ].insert( ImplicitOffset( offset.GetFrom(), offset.GetAmount() ) );
        return result.second;
    }

    const StdSet<ImplicitOffset>& GetOffsetsTo( const int to ) const {
        return m_offsets[ to ];
    }

  private:
    vec< StdSet<ImplicitOffset> > m_offsets;
};


class OffsetTracker {
  public:
    OffsetTracker() {}

    OffsetTracker( const int size )
        : m_offsets( size ) {}

    OffsetTracker( const MutableOffsetTracker& mutableTracker );

    void ConvertFrom( const MutableOffsetTracker& mutableTracker );

    int size() const {
        return m_offsets.size();
    }

    void resize( const int size ) {
        m_offsets.resize( size );
    }

    const ImplicitOffsetVec& GetOffsetsTo( const int to ) const {
        return m_offsets[ to ];
    }

    typedef std::pair<ImplicitOffsetVec::const_iterator,
            ImplicitOffsetVec::const_iterator> Range;

    Range
    GetOffsetsToFrom( const int to, const int from ) const {
        return equal_range( m_offsets[to].begin(), m_offsets[to].end(),
                            ImplicitOffset( from, 0 ), ImplicitOffset::OrderByFrom() );
    }

    void Write( const String& filename ) const {
        m_offsets.WriteAll( filename );
    }

    void Read( const String& filename ) {
        m_offsets.ReadAll( filename );
    }

    friend class OffsetTrackerBuilder;

  private:
    VecImplicitOffsetVec m_offsets;
};

#endif
