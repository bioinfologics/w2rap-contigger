// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <paths/PathFinder.h>
#include <kmers/kmatch/KMatch.h>
#include <util/OutputLog.h>
#include "CoreTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/ImprovePath.h"
#include "paths/long/large/PullAparter.h"
#include "paths/long/large/Simplify.h"
#include "kmers/kmatch/KMatch.h"
#include "GFADump.h"
#include "OverlapValidator.h"


void graph_status(const HyperBasevector &hb) {
    uint64_t total=0;
    uint64_t null_sized=0;
    std::vector<uint64_t> sizes;
    sizes.reserve(hb.EdgeObjectCount());
    for (auto i=0; i<hb.EdgeObjectCount(); ++i) {
        sizes.push_back(hb.EdgeObject(i).size()-hb.K()+1);
        if (hb.EdgeObject(i).size()>0)
        total+=hb.EdgeObject(i).size()-hb.K()+1;
        else ++null_sized;

    }
    OutputLog(2) << "GRAPH:  " << total << " " <<hb.K()<<"-mers in " <<hb.EdgeObjectCount()<<" edges";
    if (null_sized>0) OutputLog(2,false) << " ("<<null_sized<<" gap edges)";
    OutputLog(2,false) << std::endl;

    std::sort(sizes.begin(),sizes.end(),std::greater<uint64_t>());
    uint64_t n20=0,n50=0,n80=0,t20=20*total/100,t50=50*total/100,t80=80*total/100,t=0;
    for (auto s:sizes) {
        t+=s;
        if (0==n20 and t>=t20) n20=s;
        if (0==n50 and t>=t50) n50=s;
        if (0==n80 and t>=t80) { n80=s; break;};
    }
    OutputLog(2) << "GRAPH EDGES:  Nk20=" << n20 << "  Nk50=" << n50 << "  Nk80=" << n80 <<std::endl;

}

void path_status(const ReadPathVec &paths){
    uint64_t u=0,ps=0,pm=0,none=0,single=0,both=0;
    bool first=true;
    uint8_t uends=0;
    for (auto &p:paths) {
        if (p.size() == 0) {
            ++u;
            ++uends;
        }
        else if (p.size() == 1) ++ps;
        else ++pm;

        if (not first){
            if (uends==2) ++none;
            else if (uends==1) ++single;
            else ++both;
            uends=0;
        }
        first=!first;
    }
    OutputLog(2) << "PATHS:  empty: " << u << "  1: " << ps <<"  2+: " << pm << std::endl;
    OutputLog(2) << "PAIR ENDS: none: " << none << "  single: " << single << "  both: " << both << std::endl;

}

void graph_path_pairs_status(HyperBasevector &hbv, ReadPathVec & paths){
    graph_status(hbv);
    path_status(paths);
}


//Needed by old versions
void AnalyzeBranches(HyperBasevector &hb, vec<int> &to_right, const vec<int> &inv2,
                     ReadPathVec &paths2, const Bool ANALYZE_BRANCHES_REV,
                     const int min_ratio2, const Bool ANALYZE_BRANCHES_VERBOSE) {
    double clock0 = WallClockTime();
    vec<int> to_left;
    hb.ToLeft(to_left);
    for (int64_t i = 0; i < (int64_t) paths2.size(); i++) {
        ReadPath &p = paths2[i];
        for (int64_t j = 0; j < (int64_t) p.size(); j++) {
            if (p[j] >= hb.EdgeObjectCount()) p[j] = -1;
            if (j > 0 && p[j - 1] >= 0 && p[j] >= 0
                && to_right[p[j - 1]] != to_left[p[j]]) { p[j] = -1; }
        }
    }

    // Heuristics.

    const int max_dist = 4;
    const int min_ratio = 5;
    const int max_kill = 2;

    vec<std::pair<int, int> > breaks;
    vec<vec<int> > froms(hb.EdgeObjectCount()), tos(hb.EdgeObjectCount());
    LogTime(clock0, "analyzing branches 0");
    double clock1 = WallClockTime();
    for (int pass = 1; pass <= 2; pass++) {
        const int batch = 10000;
        int64_t npids = paths2.size() / 2;
#pragma omp parallel for
        for (int64_t bi = 0; bi < npids; bi += batch) {
            vec<std::pair<int, int> > PP;
            for (int64_t pid = bi; pid < Min(bi + batch, npids); pid++) {
                vec<int> x, y;
                for (int64_t j = 0; j < (int64_t) paths2[2 * pid].size(); j++)
                    x.push_back(paths2[2 * pid][j]);
                for (int64_t j = 0; j < (int64_t) paths2[2 * pid + 1].size(); j++)
                    y.push_back(paths2[2 * pid + 1][j]);
                y.ReverseMe();
                for (int j = 0; j < y.isize(); j++)
                    if (y[j] >= 0) y[j] = inv2[y[j]];
                if (pass == 2) {
                    swap(x, y);
                    x.ReverseMe(), y.ReverseMe();
                    for (int j = 0; j < x.isize(); j++)
                        if (x[j] >= 0) x[j] = inv2[x[j]];
                    for (int j = 0; j < y.isize(); j++)
                        if (y[j] >= 0) y[j] = inv2[y[j]];
                }
                std::pair<vec<int>, vec<int> > p = std::make_pair(x, y);
                vec<std::pair<int, int> > P;
                for (int j1 = 0; j1 < p.first.isize() - 1; j1++) {
                    if (p.first[j1] >= 0 && p.first[j1 + 1] >= 0)
                        P.push(p.first[j1], p.first[j1 + 1]);
                }
                for (int j1 = 0; j1 < p.second.isize() - 1; j1++) {
                    if (p.second[j1] >= 0 && p.second[j1 + 1] >= 0)
                        P.push(p.second[j1], p.second[j1 + 1]);
                }
                for (int j1 = 0; j1 < p.first.isize(); j1++) {
                    int x1 = p.first[j1];
                    if (x1 >= 0) {
                        int m = Position(p.second, x1);
                        if (m < 0 && p.second.nonempty()
                            && p.second[0] >= 0) { P.push(x1, p.second[0]); }
                    }
                }
                UniqueSort(P);
                PP.append(P);
            }
#pragma omp critical
            {
                for (int j = 0; j < PP.isize(); j++) {
                    froms[PP[j].first].
                            push_back(PP[j].second);
                    tos[PP[j].second].
                            push_back(PP[j].first);
                }
            }
        }
    }
    LogTime(clock1, "analyzing branches 1");
    double clock1b = WallClockTime();
#pragma omp parallel for
    for (int e = 0; e < hb.EdgeObjectCount(); e++) {
        Sort(froms[e]);
        Sort(tos[e]);
    }
    if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\nforward reach:\n";
    LogTime(clock1b, "analyzing branches 1b");
    double clock2 = WallClockTime();

    /*
    for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
    {    int v = to_right[e];
         if ( hb.From(v).size( ) != 2 || hb.To(v).size( ) > 1 ) continue;
         int n1 = 0, n2 = 0;
         for ( int j = 0; j < froms[e].isize( ); j++ )
         {    if ( froms[e][j] = hb.IFrom( v, 0 ) ) n1++;
              if ( froms[e][j] = hb.IFrom( v, 1 ) ) n2++;    }
         if ( n1 >= 10 && n2 <= 1 ) breaks.push( e, hb.IFrom( v, 1 ) );
         if ( n2 >= 10 && n1 <= 1 ) breaks.push( e, hb.IFrom( v, 0 ) );    }
    */

    for (int e = 0; e < hb.EdgeObjectCount(); e++) {
        int v = to_right[e];
        if (hb.From(v).size() <= 1) continue;
        if (hb.To(v).size() > 1) continue;
        vec<vec<int> > follow(hb.From(v).size());
        vec<int> branches;
        for (int j = 0; j < hb.From(v).isize(); j++) {
            int f = hb.EdgeObjectIndexByIndexFrom(v, j);
            branches.push_back(f);
        }
        int nbranches = branches.size();

        for (int j = 0; j < hb.From(v).isize(); j++) {
            int f = hb.EdgeObjectIndexByIndexFrom(v, j);
            int w = to_right[f];
            for (int l = 0; l < hb.From(w).isize(); l++)
                follow[j].push_back(hb.EdgeObjectIndexByIndexFrom(w, l));
        }

        for (int dpass = 1; dpass < max_dist; dpass++) {
            for (int i = 0; i < nbranches; i++) {
                int n = follow[i].size();
                for (int j = 0; j < n; j++) {
                    int w = to_right[follow[i][j]];
                    follow[i].append(hb.FromEdgeObj(w));
                }
                UniqueSort(follow[i]);
            }
        }

        vec<int> fr, count;
        for (int i = 0; i < froms[e].isize(); i++) {
            int j = froms[e].NextDiff(i);
            int c = j - i, f = froms[e][i];
            fr.push_back(f), count.push_back(c);
            i = j - 1;
        }
        vec<Bool> to_delete(fr.size(), False);
        for (int i = 0; i < fr.isize(); i++) {
            vec<int> homes;
            for (int j = 0; j < follow.isize(); j++)
                if (Member(follow[j], fr[i])) homes.push_back(j);
            if (homes.size() == follow.size()) count[i] = 0;
            if (homes.solo()) {
                for (int j = 0; j < fr.isize(); j++) {
                    if (fr[j] == hb.EdgeObjectIndexByIndexFrom(v, homes[0])) {
                        count[j] += count[i];
                        count[i] = 0;
                    }
                }
            }
        }
        for (int i = 0; i < fr.isize(); i++)
            if (count[i] == 0) to_delete[i] = True;
        EraseIf(fr, to_delete), EraseIf(count, to_delete);
        vec<int> s1 = fr, s2 = branches;
        Sort(s1), Sort(s2);
        if (s1 == s2 && s1.size() == 2) {
            if (count[0] < min_ratio * count[1]
                && count[1] < min_ratio * count[0]) { continue; }
        }
        ReverseSortSync(count, fr);
        if (ANALYZE_BRANCHES_VERBOSE) {
            std::cout << e << " -->";
            for (int i = 0; i < fr.isize(); i++)
                std::cout << " " << fr[i] << "[" << count[i] << "]";
        }
        if (count.size() >= 2 && count[0] >= min_ratio2 * Max(1, count[1])
            && count[1] <= max_kill && Member(branches, fr[0])) {
            if (ANALYZE_BRANCHES_VERBOSE) std::cout << " -- RECOMMEND PRUNING";
            for (int j = 0; j < branches.isize(); j++)
                if (branches[j] != fr[0]) breaks.push(e, branches[j]);
        }
        if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";
    }
    UniqueSort(breaks);
    for (int i = 0; i < breaks.isize(); i++) {
        int e = breaks[i].first, f = breaks[i].second;
        int n = hb.N();
        hb.AddVertices(2);
        hb.GiveEdgeNewFromVx(f, to_right[e], n);
        to_left[f] = n;
        int re = inv2[e], rf = inv2[f];
        if (re >= 0 && rf >= 0) {
            hb.GiveEdgeNewToVx(rf, to_right[rf], n + 1);
            to_right[rf] = n + 1;
        }
    }

    if (ANALYZE_BRANCHES_REV) {
        vec<std::pair<int, int> > breaksr;
        if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\nbackward reach:\n";

        /*
        for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
        {    int v = to_left[e];
             if ( hb.To(v).size( ) != 2 || hb.From(v).size( ) > 1 ) continue;
             int n1 = 0, n2 = 0;
             for ( int j = 0; j < tos[e].isize( ); j++ )
             {    if ( tos[e][j] = hb.ITo( v, 0 ) ) n1++;
                  if ( tos[e][j] = hb.ITo( v, 1 ) ) n2++;    }
             if ( n1 >= 10 && n2 <= 1 ) breaksr.push( hb.ITo( v, 1 ), e );
             if ( n2 >= 10 && n1 <= 1 ) breaksr.push( hb.ITo( v, 0 ), e );    }
        */

        for (int e = 0; e < hb.EdgeObjectCount(); e++) {
            int v = to_left[e];
            if (hb.To(v).size() <= 1) continue;
            if (hb.From(v).size() > 1) continue;
            vec<vec<int> > preceed(hb.To(v).size());
            vec<int> branches;
            for (int j = 0; j < hb.To(v).isize(); j++) {
                int f = hb.EdgeObjectIndexByIndexTo(v, j);
                branches.push_back(f);
            }
            int nbranches = branches.size();

            for (int j = 0; j < hb.To(v).isize(); j++) {
                int f = hb.EdgeObjectIndexByIndexTo(v, j);
                int w = to_left[f];
                for (int l = 0; l < hb.To(w).isize(); l++)
                    preceed[j].push_back(hb.EdgeObjectIndexByIndexTo(w, l));
            }

            for (int dpass = 1; dpass < max_dist; dpass++) {
                for (int i = 0; i < nbranches; i++) {
                    int n = preceed[i].size();
                    for (int j = 0; j < n; j++) {
                        int w = to_left[preceed[i][j]];
                        preceed[i].append(hb.ToEdgeObj(w));
                    }
                    UniqueSort(preceed[i]);
                }
            }

            vec<int> fr, count;
            for (int i = 0; i < tos[e].isize(); i++) {
                int j = tos[e].NextDiff(i);
                int c = j - i, f = tos[e][i];
                if (to_right[f] == to_left[e]) { fr.push_back(f), count.push_back(c); }
                i = j - 1;
            }
            vec<Bool> to_delete(fr.size(), False);
            for (int i = 0; i < fr.isize(); i++) {
                vec<int> homes;
                for (int j = 0; j < preceed.isize(); j++)
                    if (Member(preceed[j], fr[i])) homes.push_back(j);
                if (homes.size() == preceed.size()) count[i] = 0;
                if (homes.solo()) {
                    for (int j = 0; j < fr.isize(); j++) {
                        if (fr[j] == hb.EdgeObjectIndexByIndexTo(v, homes[0])) {
                            count[j] += count[i];
                            count[i] = 0;
                        }
                    }
                }
            }
            for (int i = 0; i < fr.isize(); i++)
                if (count[i] == 0) to_delete[i] = True;
            EraseIf(fr, to_delete), EraseIf(count, to_delete);
            vec<int> s1 = fr, s2 = branches;
            Sort(s1), Sort(s2);
            if (s1 == s2 && s1.size() == 2) {
                if (count[0] < min_ratio * count[1]
                    && count[1] < min_ratio * count[0]) { continue; }
            }
            ReverseSortSync(count, fr);
            if (ANALYZE_BRANCHES_VERBOSE) {
                std::cout << e << " <--";
                for (int i = 0; i < fr.isize(); i++)
                    std::cout << " " << fr[i] << "[" << count[i] << "]";
            }
            if (count.size() >= 2 && count[0] >= min_ratio2 * Max(1, count[1])
                && count[1] <= max_kill && Member(branches, fr[0])) {
                if (ANALYZE_BRANCHES_VERBOSE) std::cout << " -- RECOMMEND PRUNING";
                for (int j = 0; j < branches.isize(); j++)
                    if (branches[j] != fr[0]) breaksr.push(branches[j], e);
            }
            if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";
        }
        UniqueSort(breaksr);
        for (int i = 0; i < breaksr.isize(); i++) {
            int e = breaksr[i].first, f = breaksr[i].second;
            int n = hb.N();
            hb.AddVertices(2);
            hb.GiveEdgeNewToVx(e, to_left[f], n);
            to_right[e] = n;

            int re = inv2[e], rf = inv2[f];
            if (re >= 0 && rf >= 0) {
                hb.GiveEdgeNewFromVx(re, to_left[re], n + 1);
                to_left[re] = n + 1;
            }
        }
        breaks.append(breaksr);
    }

    int nb = breaks.size();
    for (int i = 0; i < nb; i++)
        breaks.push(inv2[breaks[i].second], inv2[breaks[i].first]);
    UniqueSort(breaks);
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t) paths2.size(); i++) {
        ReadPath &p = paths2[i];
        Bool bad = False;
        for (int j = 0; j < ((int) p.size()) - 1; j++) {
            std::pair<int, int> x = std::make_pair(p[j], p[j + 1]);
            if (BinMember(breaks, x)) bad = True;
        }
        if (bad) p.resize(0);
    }
    if (ANALYZE_BRANCHES_VERBOSE) std::cout << "\n";
    LogTime(clock2, "analyzing branches 2");
}



void update_read_placements(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals){
    path_improver pimp;
    vec<int64_t> ids;
    ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,
                 /*IMPROVE_PATHS_LARGE*/False, False);
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

}

void update_read_placements_kmatch(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals) {
    std::cout<<Date()<<": creating kmatch object"<<std::endl;
    KMatch km(31);
    km.Hbv2Index(hb,31);
    //uint64_t mapped=0,c=0;
    /*km.Hbv2Map(hb);
    uint64_t TEST_READ_COUNT=100000;
    auto wc1=WallClockTime();
    std::cout<<Date()<<": Counting matches on "<<TEST_READ_COUNT<<" reads"<<std::endl;

    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.countReadMatches(bases[i].ToString())>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in countReadMatches "<<std::endl;
    wc1=WallClockTime();
    std::cout<<Date()<<": Looking up "<<TEST_READ_COUNT<<" reads"<<std::endl;
    mapped=0;
    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.lookupRead(bases[i].ToString()).size()>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in lookupRead "<<std::endl;
    wc1=WallClockTime();
    std::cout<<Date()<<": Looking up "<<TEST_READ_COUNT<<" reads in MAP"<<std::endl;
    mapped=0;
    for (auto i=0;i<TEST_READ_COUNT;++i){
        if (km.lookupReadInMap(bases[i].ToString()).size()>0) ++mapped;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<TEST_READ_COUNT<<" reads have hits to edges, "<<TimeSince(wc1)<<" spent in lookupReadInMap "<<std::endl;
    */
    std::atomic_uint_fast64_t mapped(0),mapped50(0),mappedUniquely(0),c(0);
    #pragma omp parallel for
    for (auto i=0;i<bases.size();++i){
        auto hits=km.lookupRead(bases[i].ToString());;
        if (hits.size()>0) ++mapped;
        if (hits.size()>(bases[i].size() - 31 + 1)*.5) {
            ++mapped50;
            bool unique=true;
            int64_t e=-1;
            for (auto & h:hits){
                if (e==-1) e=h.edge_id;
                if (e!=h.edge_id) {
                    unique =false;
                    break;
                }
            }
            if (unique) ++mappedUniquely;
        }

        ++c;
        if (c%100000==0) std::cout<<Date()<<": "<<mapped<<" / "<<c<<" reads have hits to edges"<<std::endl;
    }
    std::cout<<Date()<<": "<<mapped<<" / "<<bases.size()<<" reads have hits to edges"<<std::endl;
    std::cout<<Date()<<": "<<mapped50<<" / "<<bases.size()<<" reads have 50%+ kmers hitting to edges"<<std::endl;
    std::cout<<Date()<<": "<<mappedUniquely<<" / "<<bases.size()<<" reads have 50%+ kmers hitting same edge and no kmers hitting others"<<std::endl;
}

#include <iterator>
#include <SpectraCn.hpp>
#include <unordered_set>

template<class ForwardIt, class T>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value)
{
    // Note: BOTH type T and the type after ForwardIt is dereferenced
    // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
    // This is stricter then lower_bound requirement (see above)

    first = std::lower_bound(first, last, value);
    return (first != last && value == *first) ? first : last;
}



void remove_unsupported_edges(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals, const int MAX_SUPP_DEL, const int min_mult){
    uint64_t delcount=1;
    uint64_t pass=0;
    //while(delcount) {
        vec<int> toLeft,toRight;
        hb.ToLeft(toLeft);
        hb.ToRight(toRight);
        vec<int> dels;
        {
            std::vector<int> support_in(hb.EdgeObjectCount(), 0);
            std::vector<int> support_out(hb.EdgeObjectCount(), 0);
            //this is support of the edge as exit from its input vertex.

            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j >= 1) {
                        support_in[e]++;
                        if (inv[e] >= 0) support_out[inv[e]]++;
                    }
                    if (j < (int64_t) paths[id].size() - 1) {
                        support_out[e]++;
                        if (inv[e] >= 0) support_in[inv[e]]++;
                    }
                }
            }

            //First (trivial) decisions: edges with no in and no out support
#pragma omp parallel for
            for (auto e=0;e<hb.EdgeObjectCount();e++){
                if (hb.EdgeObject(e).size()>1000) continue;
                if (support_in[e]<=MAX_SUPP_DEL && support_out[e]<=MAX_SUPP_DEL){
                    auto vfrom=toLeft[e];
                    auto vto=toRight[e];
                    bool alternative_from=hb.To(vfrom).size()==0;
                    bool alternative_to=hb.From(vto).size()==0;

                    for (auto fi=0;fi<hb.From(vfrom).isize();++fi) {
                        auto other = hb.EdgeObjectIndexByIndexFrom(vfrom, fi);
                        if (other != e and (inv[e] < 0 or other != inv[e]) and
                            support_in[other] >= min_mult * support_in[e] and
                            support_in[other] > min_mult and support_in[other] > MAX_SUPP_DEL) {
                            alternative_from = true;
                            break;
                        }
                    }

                    for (auto ti=0;ti<hb.To(vto).isize();++ti) {
                        auto other = hb.EdgeObjectIndexByIndexTo(vto, ti);
                        if (other != e and (inv[e] < 0 or other != inv[e]) and
                            support_out[other] >= min_mult * support_out[e] and
                            support_out[other] > min_mult and support_out[other] > MAX_SUPP_DEL) {
                            alternative_to = true;
                            break;
                        }
                    }
#pragma omp critical
                    if (alternative_from and alternative_to) {
                        dels.push_back(e);
                        if (inv[e]>0) dels.push_back(inv[e]);
                    }
                }
            }
        }
        auto before = hb.EdgeObjectCount();
        delcount = dels.size();

        /*
        //Update paths first, just to be sure
        for (int64_t i = 0; i < (int64_t) paths.size(); i++) {
            for (int64_t j = 0; j < (int64_t) paths[i].size(); j++) {
                if (binary_find(dels.begin(), dels.end(), paths[i][j]) != dels.end()) {
                    paths[i].resize(j);
                    break;
                }
            }
        }*/

        //GFADumpDetail("unsupported_paths_marked_detail"+std::to_string(pass),hb,inv,dels);
        hb.DeleteEdges(dels);
        Cleanup(hb, inv, paths);
        OutputLog(2) << delcount << " / " << before << " unsupported edges removed, "
                  << hb.EdgeObjectCount() << " edges after cleanup" << std::endl;
        // Improve read placements and delete funky pairs.
        OutputLog(2) << "rerouting paths and cleaning pairs" << std::endl;
        /*ReroutePaths(hb, inv, paths, bases, quals);
        DeleteFunkyPathPairs(hb, inv, bases, paths, False);
        std::cout << Date() << ": improving paths" << std::endl;
        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp, false, False);*/
        pass++;
    //}
}

void full_cleanup(HyperBasevector &hb, vec<int> &inv, ReadPathVec &paths, const int tampsize, const int hangssize){
    Tamp(hb, inv, paths, tampsize);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);
}



void SimplifyEXP(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
              const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
              const int MAX_DEL2,
              const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
              const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
              const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
              const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
              const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
              const Bool FINAL_TINY, const Bool UNWIND3, const bool RUN_PATHFINDER, const bool dump_pf_files, const bool VERBOSE_PATHFINDER) {


    // Improve read placements and delete funky pairs.
    OutputLog(2) << "rerouting paths and cleaning pairs" << std::endl;
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);
    /*if (IMPROVE_PATHS) {
        std::cout << Date() << ": improving paths" << std::endl;

        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp, IMPROVE_PATHS_LARGE, False);
    }*/
    path_status(paths);
    /*OverlapValidator oval(hb,inv,paths);
    GFADumpDetail("before_ovlpval_detail",hb,inv);
    oval.compute_overlap_support();
    oval.analyse_complex_overlaps();
    std::vector<int> paint;
    paint.reserve(hb.EdgeObjectCount());
    for (int e:oval.find_perfect_tips(1000,5)) paint.push_back(e);
    std::cout<<Date()<<": "<<paint.size()<<" perfect tips found"<<std::endl;
    GFADumpDetail("ovlpval_perfect_tips_detail",hb,inv,paint);
    hb.DeleteEdges(paint);
    Cleanup(hb,inv,paths);
    graph_status(hb);
    path_status(paths);*/


    //Remove unsupported edges in certain situations.
    const int min_mult=5;
    OutputLog(2) << "removing alternative edges with input support <="<<MAX_SUPP_DEL << std::endl;
    remove_unsupported_edges(hb,inv,paths,bases,quals,MAX_SUPP_DEL,min_mult);




    // Clean up assembly.
    OutputLog(2) << "removing small components"<<std::endl;
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);


    /*
    if (TAMP_EARLY) {
        std::cout << Date() << ": Tamping" << std::endl;
        Tamp(hb, inv, paths, 0);
    }

    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);*/
    graph_status(hb);
    path_status(paths);

    /*std::cout << Date() << ": analysing branches" << std::endl;


    AnalyzeBranches(hb, to_right, inv, paths, True, MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE2);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, MAX_DEL2);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    graph_status(hb);
    path_status(paths);*/

    OutputLog(2) << "popping bubbles" << std::endl;
    PopBubbles(hb, inv, bases, quals, paths);

    Cleanup(hb, inv, paths);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);
    Tamp(hb, inv, paths, 10);
    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);
    graph_status(hb);
    path_status(paths);

    // Pull apart.



    if (RUN_PATHFINDER) {
        //TODO: remove pull aparter once the pathfinder solves all repeats (just using distance 2 on OverlapValidator should do)
        VecULongVec invPaths;
        OutputLog(2) << "pulling apart canonical repeats" << std::endl;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        OutputLog(2) << count << " repeats separated, " << pa.getRemovedReadPaths()
                  << " read paths removed" << std::endl;
        graph_status(hb);
        path_status(paths);



        OutputLog(2) << "running pathfinder" << std::endl;
        invert(paths, invPaths, hb.EdgeObjectCount());
        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_start.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_start.paths").c_str());
        }

        PathFinder pf(hb, inv, paths, invPaths, bases, quals, 5,VERBOSE_PATHFINDER);
        pf.unroll_loops(800);
        RemoveUnneededVertices2(hb, inv, paths);
        Cleanup(hb, inv, paths);

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_unrolled_loops.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_unrolled_loops.paths").c_str());
        }
        invPaths.clear();
        invert( paths, invPaths, hb.EdgeObjectCount( ) );
        pf.untangle_complex_in_out_choices(700);
        RemoveUnneededVertices2(hb, inv, paths);
        Cleanup(hb, inv, paths);
        DeleteFunkyPathPairs(hb, inv, bases, paths, False);
        Tamp(hb, inv, paths, 10);
        RemoveHangs(hb, inv, paths, 700);
        Cleanup(hb, inv, paths);
        RemoveSmallComponents3(hb);
        graph_status(hb);
        path_status(paths);

        if (dump_pf_files) {
            BinaryWriter::writeFile(fin_dir + "/pf_end.hbv", hb);
            WriteReadPathVec(paths,(fin_dir + "/pf_end.paths").c_str());
        }
    } else {
        OutputLog(2) << "pulling apart canonical repeats" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        OutputLog(2) << count << " repeats separated, " << pa.getRemovedReadPaths() << " read paths removed" << std::endl;
        graph_status(hb);
        path_status(paths);
    }
    // Improve paths.

    if (IMPROVE_PATHS) {
        path_improver pimp;
        vec<int64_t> ids;
        OutputLog(2) << "improving paths" << std::endl;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,IMPROVE_PATHS_LARGE, False);
        path_status(paths);
    }

    // Extend paths.

    if (EXT_FINAL) {
        OutputLog(2) << "extending paths" << std::endl;
        vec<int> to_left;
        vec<int> to_right;

        hb.ToLeft(to_left), hb.ToRight(to_right);
        int ext = 0;
        auto qvItr = quals.begin();
        for (int64_t id = 0; id < (int64_t) paths.size(); id++, ++qvItr) {
            Bool verbose = False;
            const int min_gain = 20;
            ReadPath p = paths[id];
            ExtendPath2(paths[id], id, hb, to_left, to_right, bases[id], *qvItr, min_gain, verbose, EXT_FINAL_MODE);
            if (p != paths[id]) ext++;
        }
        OutputLog(2) << ext << " paths extended" << std::endl;
        path_status(paths);
    }

    // Degloop.

    if (DEGLOOP) {
        OutputLog(2) << "deglooping" << std::endl;
        Degloop(DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST);
        RemoveHangs(hb, inv, paths, 700);
        Cleanup(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }

    // Unwind three-edge plasmids.

    if (UNWIND3) {
        OutputLog(2) << "unwinding 3-edge plasmids" << std::endl;
        UnwindThreeEdgePlasmids(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }

    // Remove tiny stuff.

    if (FINAL_TINY) {
        OutputLog(2) << "removing small components" << std::endl;
        RemoveSmallComponents3(hb, True);
        Cleanup(hb, inv, paths);
        CleanupLoops(hb, inv, paths);
        RemoveUnneededVerticesGeneralizedLoops(hb, inv, paths);
        graph_status(hb);
        path_status(paths);
    }
}

void Simplify(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals,
              const int MAX_SUPP_DEL, const Bool TAMP_EARLY, const int MIN_RATIO2,
              const int MAX_DEL2,
              const Bool ANALYZE_BRANCHES_VERBOSE2, const String &TRACE_SEQ,
              const Bool DEGLOOP, const Bool EXT_FINAL, const int EXT_FINAL_MODE,
              const Bool PULL_APART_VERBOSE, const vec<int> &PULL_APART_TRACE,
              const int DEGLOOP_MODE, const double DEGLOOP_MIN_DIST,
              const Bool IMPROVE_PATHS, const Bool IMPROVE_PATHS_LARGE,
              const Bool FINAL_TINY, const Bool UNWIND3) {
    std::string name;
    // Improve read placements and delete funky pairs.
    OutputLog(2) << "rerouting paths" << std::endl;

    ReroutePaths(hb, inv, paths, bases, quals);
//    support = get_edges_support(hb, inv, paths);
//    std::cout << "After ReroutePaths, support for edge10797697: " << support[10797697] << std::endl;

    DeleteFunkyPathPairs(hb, inv, bases, paths, False);
//    support = get_edges_support(hb, inv, paths);
//    std::cout << "After DeleteFunkyPathPairs, support for edge10797697: " << support[10797697] << std::endl;


    if (IMPROVE_PATHS) {
        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,
                     IMPROVE_PATHS_LARGE, False);
        graph_path_pairs_status(hb,paths);
    }
//    support = get_edges_support(hb, inv, paths);
//    std::cout << "After ImprovePaths, support for edge10797697: " << support[10797697] << std::endl;

    // Tip clipping in situations where a small tip ( size < 2*K ) is adjacent to a large contig ( size >= 5*K )
    // TODO: Explore doing this before/after updating the paths, this should confirm suspicions of poorly rerouted paths before this step!
    {
        vec<int> dels;

        // Cleanup of "From" edges
        for (int v = 0; v < hb.N(); v++) {
            bool to_print(true);
            if (hb.From(v).size() == 2) {
                int e1 = hb.EdgeObjectIndexByIndexFrom(v, 0);
                int e2 = hb.EdgeObjectIndexByIndexFrom(v, 1);
                if (to_print) {
                    std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                    std::cout << "Length(e1) = " << hb.EdgeObject(e1).size() << " Length(e2) = " << hb.EdgeObject(e2).size()
                              << std::endl;
                }
                if (hb.EdgeObject(e1).size() > hb.EdgeObject(e2).size()) std::swap(e1, e2);
                if (hb.EdgeObject(e1).size() <= 2*hb.K() && hb.EdgeObject(e2).size() >= 5*hb.K() /*&& hb.EdgeObject(e1).size() < 2*hb.K() */) {
                    dels.push_back(e1);
                }
            }
        }

        // Cleanup of "To" edges
        for (int v = 0; v < hb.N(); v++) {
            bool to_print(true);
            if (hb.From(v).size() == 2) {
                int e1 = hb.EdgeObjectIndexByIndexFrom(v, 0);
                int e2 = hb.EdgeObjectIndexByIndexFrom(v, 1);
                if (to_print) {
                    std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                    std::cout << "Length(e1) = " << hb.EdgeObject(e1).size() << " Length(e2) = " << hb.EdgeObject(e2).size()
                              << std::endl;
                }
                if (hb.EdgeObject(e1).size() > hb.EdgeObject(e2).size()) std::swap(e1, e2);
                if (hb.EdgeObject(e1).size() <= 2*hb.K() && hb.EdgeObject(e2).size() >= 5*hb.K() /*&& hb.EdgeObject(e1).size() < 2*hb.K() */) {
                    dels.push_back(e1);
                }
            }
        }

        hb.DeleteEdges(dels);
        Cleanup(hb, inv, paths);
    }


    // Remove unsupported edges in certain situations.
    OutputLog(2) << "removing unsupported edges" << std::endl;
    {
        const int min_mult = 10;
        vec<int> dels;
        std::ofstream unsupported_edges_sizes("unsupported.sizes");
        {
            {
                vec<int> support(hb.EdgeObjectCount(), 0);
                for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                    for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                        int e = paths[id][j];
                        if (j >= 1) support[e]++;
                        if (inv[e] >= 0 && j < (int64_t) paths[id].size() - 1)
                            support[inv[e]]++;
                    }
                }

                for (int v = 0; v < hb.N(); v++) {
                    bool to_print(true);
                    if (hb.From(v).size() == 2) {
                        int e1 = hb.EdgeObjectIndexByIndexFrom(v, 0);
                        int e2 = hb.EdgeObjectIndexByIndexFrom(v, 1);
                        if (to_print) {
                            std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                            std::cout << "Support(e1) = " << support[e1] << " Support(e2) = " << support[e2]
                                      << std::endl;
                        }
                        if (support[e1] > support[e2]) std::swap(e1, e2);
                        int s1 = support[e1], s2 = support[e2];
                        if (to_print) {
                            std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                            std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;

                            std::cout << "s1(" << s1 << ") <= MAX_SUPP_DEL(" << MAX_SUPP_DEL << ") && s2(" << s2
                                      << ") >= " << min_mult << " * Max(1, " << s1 << ")" << std::endl;
                        }
                        if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1) /*&& hb.EdgeObject(e1).size() < 2*hb.K() */) {
                            if (to_print) std::cout << "dels.push_back(" << e1 << ")" << std::endl;
                            dels.push_back(e1);
                            unsupported_edges_sizes << e1 << "\t" << hb.EdgeObject(e1).size() << "\n";
                        }
                    }
                }
            }
            {
                vec<int> support(hb.EdgeObjectCount(), 0);
                for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                    for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                        int e = paths[id][j];
                        if (j < (int64_t) paths[id].size() - 1) support[e]++;
                        if (inv[e] >= 0 && j >= 1) support[inv[e]]++;
                    }
                }
                for (int v = 0; v < hb.N(); v++) {
                    bool to_print(true);
                    if (hb.To(v).size() == 2) {
                        int e1 = hb.EdgeObjectIndexByIndexTo(v, 0);
                        int e2 = hb.EdgeObjectIndexByIndexTo(v, 1);
                        if (to_print) {
                            std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                            std::cout << "Support(e1) = " << support[e1] << " Support(e2) = " << support[e2]
                                      << std::endl;
                        }
                        if (support[e1] > support[e2]) std::swap(e1, e2);
                        int s1 = support[e1], s2 = support[e2];
                        if (to_print) {
                            std::cout << "e1 = " << e1 << " e2 = " << e2 << std::endl;
                            std::cout << "s1 = " << s1 << " s2 = " << s2 << std::endl;

                            std::cout << "s1(" << s1 << ") <= MAX_SUPP_DEL(" << MAX_SUPP_DEL << ") && s2(" << s2
                                      << ") >= " << min_mult << " * Max(1, " << s1 << ")" << std::endl;
                        }
                        if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1) /* && hb.EdgeObject(e1).size() < 2*hb.K() */) {
                            if (to_print) std::cout << "dels.push_back(" << e1 << ")" << std::endl;
                            dels.push_back(e1);
                            unsupported_edges_sizes << e1 << "\t" << hb.EdgeObject(e1).size() << "\n";
                        }
                    }
                }
            }
        }
        // TODO: Print all dels here to check if all these edges are bogus
        auto before=hb.EdgeObjectCount();
        auto delcount=dels.size();

        {
            std::ofstream to_delete("to_delete.edges");
            int i = 0;
            for (const auto d:dels) {
                if (d % 2 == 0) { // Only output the canonical edges, the u
                    to_delete << "edge" << d << "\n";
                }
            }
            to_delete << std::endl;
            std::cout << "There were " << i << " deleted edges" << std::endl;
        }

        name = "step7_before_DeleteEdges";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);
        hb.DeleteEdges(dels);

        name = "step7_before_Cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        { // Cleanup (inlined)
            {
                vec<Bool> used;
                hb.Used(used);
                for (int64_t i = 0; i < (int64_t) paths.size(); i++) {
                    for (int64_t j = 0; j < (int64_t) paths[i].size(); j++) {
                        if (paths[i][j] < 0 || paths[i][j] >= hb.EdgeObjectCount() || !used[paths[i][j]]) {
                            paths[i].resize(j);
                            break;
                        }
                    }
                }
            }
            RemoveUnneededVertices2(hb, inv, paths, true);
            CleanupCore(hb, inv, paths);
        }

        OutputLog(2) << delcount << " / " <<before<<" edges removed, "<<hb.EdgeObjectCount()<<" edges after cleanup"<<std::endl;
        name = "step7_edge_cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);
        graph_path_pairs_status(hb,paths);
    }



    // Clean up assembly.

    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);
    name = "step7_removesmallcomponents1_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);
    OutputLog(2) << hb.EdgeObjectCount()<<" edges after cleanup"<<std::endl;
    graph_path_pairs_status(hb,paths);



    if (TAMP_EARLY) {
        OutputLog(2) << "early tamping" << std::endl;
        Tamp(hb, inv, paths, 0);

        name = "step7_tamp0_cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        graph_path_pairs_status(hb,paths);
    }

    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);

    name = "step7_removehangs_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    OutputLog(2) <<hb.EdgeObjectCount()<<" edges after removing hangs"<<std::endl;
    graph_path_pairs_status(hb,paths);
    OutputLog(2) << "analysing branches" << std::endl;
    vec<int> to_right;
    hb.ToRight(to_right);

    AnalyzeBranches(hb, to_right, inv, paths, True, MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE2);
    Cleanup(hb, inv, paths);

    name = "step7_analysebranches_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    RemoveHangs(hb, inv, paths, MAX_DEL2);
    Cleanup(hb, inv, paths);

    name = "step7_removehangs_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    name = "step7_removesmallcomponents2_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    OutputLog(2) << hb.EdgeObjectCount()<<" edges after branch analysis and cleanup"<<std::endl;
    graph_path_pairs_status(hb,paths);
    OutputLog(2) << "popping bubbles" << std::endl;
    PopBubbles(hb, inv, bases, quals, paths);
    Cleanup(hb, inv, paths);

    name = "step7_popbubbles_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    OutputLog(2) << hb.EdgeObjectCount()<<" edges after bubble popping and cleanup"<<std::endl;
    graph_path_pairs_status(hb,paths);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

    OutputLog(2) << "tamping (700)" << std::endl;

    Tamp(hb, inv, paths, 10);

    name = "step7_tamp10_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);

    name = "step7_removehangs2_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    name = "step7_removesmallcomponents3_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    OutputLog(2) << hb.EdgeObjectCount()<<" edges after tamping, re-removing small components and cleanup"<<std::endl;
    // Pull apart.

    {
        OutputLog(2) << "pulling apart repeats" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();

        name = "step7_pull_apart";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        OutputLog(2) << count << " repeats pulled apart." << std::endl;
        OutputLog(2) << ": there were " << pa.getRemovedReadPaths() << " read paths removed during separation."<< std::endl;
    }


    OutputLog(2) << "making paths index for PathFinder" << std::endl;
    VecULongVec invPaths;
    invert(paths, invPaths, hb.EdgeObjectCount());

    OutputLog(2) << "PathFinder: unrolling loops" << std::endl;
    PathFinder(hb, inv, paths, invPaths, bases, quals).unroll_loops(800);
    OutputLog(2) << "Removing unneeded Vertices" << std::endl;
    RemoveUnneededVertices2(hb, inv, paths);
    Cleanup(hb, inv, paths);

    name = "step7_pathfinder_unroll_loops_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    invPaths.clear();
    invert( paths, invPaths, hb.EdgeObjectCount( ) );
    OutputLog(2) << "PathFinder: analysing single-direction repeats" << std::endl;
    PathFinder(hb, inv, paths, invPaths, bases, quals).untangle_complex_in_out_choices(700);
    OutputLog(2) << "Removing unneeded Vertices" << std::endl;
    RemoveUnneededVertices2(hb, inv, paths);
    Cleanup(hb, inv, paths);

    name = "step7_pathfinder_complex_in_out_cleanup";
    BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
    GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
    SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

    graph_path_pairs_status(hb,paths);


    // Improve paths.

    if (IMPROVE_PATHS) {
        path_improver pimp;
        vec<int64_t> ids;
        ImprovePaths(paths, hb, inv, bases, quals, ids, pimp,
                     IMPROVE_PATHS_LARGE, False);
        graph_path_pairs_status(hb,paths);
    }

    // Extend paths.

    if (EXT_FINAL) {
        vec<int> to_left;
        hb.ToLeft(to_left), hb.ToRight(to_right);
        int ext = 0;
        auto qvItr = quals.begin();
        for (int64_t id = 0; id < (int64_t) paths.size(); id++, ++qvItr) {
            Bool verbose = False;
            const int min_gain = 20;
            ReadPath p = paths[id];
            ExtendPath2(paths[id], id, hb, to_left, to_right, bases[id], *qvItr,
                        min_gain, verbose, EXT_FINAL_MODE);
            if (p != paths[id]) ext++;
        }
        OutputLog(2) << ext << " paths extended" << std::endl;
        graph_path_pairs_status(hb,paths);
    }

    // Degloop.

    if (DEGLOOP) {
        OutputLog(2) << "deglooping" << std::endl;
        Degloop(DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST);
        RemoveHangs(hb, inv, paths, 700);
        Cleanup(hb, inv, paths);

        name = "step7_degloop_cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        graph_path_pairs_status(hb,paths);

    }

    // Unwind three-edge plasmids.

    if (UNWIND3) UnwindThreeEdgePlasmids(hb, inv, paths);

    // Remove tiny stuff.

    if (FINAL_TINY) {
        OutputLog(2) << "removing small components" << std::endl;
        RemoveSmallComponents3(hb, True);
        Cleanup(hb, inv, paths);

        name = "step7_removesmallcomponents_final_cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        CleanupLoops(hb, inv, paths);
        RemoveUnneededVerticesGeneralizedLoops(hb, inv, paths);

        name = "unneded_generalised_cleanup";
        BinaryWriter::writeFile(fin_dir + "/" + name + ".hbv", hb);
        GFADump(std::string(fin_dir + "/" + name), hb, inv, paths, 0, 0, false);
        SpectraCN::DumpSpectraCN(hb, inv, fin_dir,  name);

        graph_path_pairs_status(hb,paths);
    }
}


void SimplifyDV(const String &fin_dir, HyperBasevector &hb, vec<int> &inv,
              ReadPathVec &paths, const vecbasevector &bases, const VecPQVec &quals) {
    //Horrible hardcoded stuff, will be contained here
    int MAX_SUPP_DEL = 0;
    int MIN_RATIO2 = 8;
    int MAX_DEL2 = 200;
    bool ANALYZE_BRANCHES_VERBOSE2 = False;
    const String TRACE_SEQ = "";
    int EXT_FINAL_MODE = 1;
    bool PULL_APART_VERBOSE = False;
    const vec<int> PULL_APART_TRACE;
    int DEGLOOP_MODE = 1;
    float DEGLOOP_MIN_DIST = 2.5;
    bool IMPROVE_PATHS_LARGE = False;

    // Improve read placements and delete funky pairs.
    OutputLog(2) << "rerouting paths" << std::endl;
    ReroutePaths(hb, inv, paths, bases, quals);
    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

    // Remove unsupported edges in certain situations.
    OutputLog(2) << "removing unsupported edges" << std::endl;
    {
        const int min_mult = 10;
        vec<int> dels;
        {
            vec<int> support(hb.EdgeObjectCount(), 0);
            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j >= 1) support[e]++;
                    if (inv[e] >= 0 && j < (int64_t) paths[id].size() - 1)
                        support[inv[e]]++;
                }
            }
#pragma omp parallel for
            for (int v = 0; v < hb.N(); v++) {
                if (hb.From(v).size() == 2) {
                    int e1 = hb.EdgeObjectIndexByIndexFrom(v, 0);
                    int e2 = hb.EdgeObjectIndexByIndexFrom(v, 1);
                    if (support[e1] > support[e2]) std::swap(e1, e2);
                    int s1 = support[e1], s2 = support[e2];
                    if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1)) {
#pragma omp critical
                        { dels.push_back(e1); }
                    }
                }
            }
        }
        {
            vec<int> support(hb.EdgeObjectCount(), 0);
            for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
                for (int64_t j = 0; j < (int64_t) paths[id].size(); j++) {
                    int e = paths[id][j];
                    if (j < (int64_t) paths[id].size() - 1) support[e]++;
                    if (inv[e] >= 0 && j >= 1) support[inv[e]]++;
                }
            }
#pragma omp parallel for
            for (int v = 0; v < hb.N(); v++) {
                if (hb.To(v).size() == 2) {
                    int e1 = hb.EdgeObjectIndexByIndexTo(v, 0);
                    int e2 = hb.EdgeObjectIndexByIndexTo(v, 1);
                    if (support[e1] > support[e2]) std::swap(e1, e2);
                    int s1 = support[e1], s2 = support[e2];
                    if (s1 <= MAX_SUPP_DEL && s2 >= min_mult * Max(1, s1)) {
#pragma omp critical
                        { dels.push_back(e1); }
                    }
                }
            }
        }
        hb.DeleteEdges(dels);
        Cleanup(hb, inv, paths);
    }

    OutputLog(2) << "removing small components" << std::endl;
    // Clean up assembly.

    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    OutputLog(2) << "early tamping" << std::endl;
    Tamp(hb, inv, paths, 0);


    RemoveHangs(hb, inv, paths, 100);
    Cleanup(hb, inv, paths);

    OutputLog(2) << "analysing branches" << std::endl;
    vec<int> to_right;
    hb.ToRight(to_right);

    AnalyzeBranches(hb, to_right, inv, paths, True, MIN_RATIO2, ANALYZE_BRANCHES_VERBOSE2);
    Cleanup(hb, inv, paths);
    RemoveHangs(hb, inv, paths, MAX_DEL2);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    OutputLog(2) << "popping bubbles" << std::endl;
    PopBubbles(hb, inv, bases, quals, paths);
    Cleanup(hb, inv, paths);

    DeleteFunkyPathPairs(hb, inv, bases, paths, False);

    OutputLog(2) << "tamping (700)" << std::endl;

    Tamp(hb, inv, paths, 10);
    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);
    RemoveSmallComponents3(hb);
    Cleanup(hb, inv, paths);

    // Pull apart.

    {
        OutputLog(2) << "pulling apart repeats" << std::endl;
        VecULongVec invPaths;
        invert(paths, invPaths, hb.EdgeObjectCount());
        PullAparter pa(hb, inv, paths, invPaths, PULL_APART_TRACE, PULL_APART_VERBOSE, 5, 5.0);
        size_t count = pa.SeparateAll();
        OutputLog(2) << count << " repeats pulled apart" << std::endl;
        OutputLog(2) << pa.getRemovedReadPaths() << " read paths removed during separation" << std::endl;
    }

    // Improve paths.
    OutputLog(2) << "improving paths" << std::endl;
    path_improver pimp;
    vec<int64_t> ids;
    ImprovePaths(paths, hb, inv, bases, quals, ids, pimp, IMPROVE_PATHS_LARGE, False);

    // Extend paths.

    OutputLog(2) << "final path extension" << std::endl;
    vec<int> to_left;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    int ext = 0;
    auto qvItr = quals.begin();
    for (int64_t id = 0; id < (int64_t) paths.size(); id++, ++qvItr) {
        Bool verbose = False;
        const int min_gain = 20;
        ReadPath p = paths[id];
        ExtendPath2(paths[id], id, hb, to_left, to_right, bases[id], *qvItr,
                    min_gain, verbose, EXT_FINAL_MODE);
        if (p != paths[id]) ext++;
    }
    OutputLog(2) << ext << " paths extended" << std::endl;


    // Degloop.

    OutputLog(2) << "deglooping" << std::endl;
    Degloop(DEGLOOP_MODE, hb, inv, paths, bases, quals, DEGLOOP_MIN_DIST);

    RemoveHangs(hb, inv, paths, 700);
    Cleanup(hb, inv, paths);


    // Unwind three-edge plasmids.

    OutputLog(2) << "unwinding 3-edge circles" << std::endl;
    UnwindThreeEdgePlasmids(hb, inv, paths);

    // Remove tiny stuff.

    OutputLog(2) << "removing small components" << std::endl;
    RemoveSmallComponents3(hb, True);
    Cleanup(hb, inv, paths);
    CleanupLoops(hb, inv, paths);
    RemoveUnneededVerticesGeneralizedLoops(hb, inv, paths);
}