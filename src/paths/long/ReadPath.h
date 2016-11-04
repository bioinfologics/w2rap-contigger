///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ReadPath.h
 *
 *  Created on: Dec 11, 2013
 *      Author: tsharpe
 */

#ifndef READPATH_H_
#define READPATH_H_

#include <vector>
#include <fstream>

// A description of a graph traversal by some sequence (a read, let's say).
// It's just a vector of edge IDs, but it also tells you how many bases at the
// start of the first edge to skip, and how many bases at the end of the last
// edge to skip.
class ReadPath : public std::vector<int>
{
public:
    ReadPath() :  mOffset(0) {}
    ReadPath( int offset )
    : mOffset(offset) {}

    ReadPath( int offset, const std::vector<int>& edge_list )
    : mOffset(offset) {
	this->assign(edge_list.begin(), edge_list.end());
    }

    int getOffset() const { return mOffset; }
    void setOffset( int offset ) { mOffset = offset; }
    void addOffset( int add ) { mOffset += add; }

    // FirstSkip was replaced by mOffset, which could be negative.  If mOffset is <= 0,
    // then FirstSkip is zero and the read completely covers the start of the edge.
    // If mOffset is positive, then the start of the read is past the start of the edge.
    unsigned getFirstSkip() const { return (mOffset < 0 ? 0u : static_cast<unsigned>(mOffset)); }
    void setFirstSkip( unsigned firstSkip ) { mOffset = firstSkip; }

    void push_front(int const i){
        this->insert(this->begin(),i);
    }
    bool same_read(ReadPath const& rp){
        return (mOffset==rp.getOffset() and *this==rp);
    }


private:
    int mOffset;
};
typedef std::vector<ReadPath> ReadPathVec;

void WriteReadPathVec(const ReadPathVec &rpv, const char * filename);
void LoadReadPathVec(ReadPathVec &rpv, const char * filename);

#endif /* READPATH_H_ */
