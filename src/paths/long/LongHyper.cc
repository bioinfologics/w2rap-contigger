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
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "efasta/EfastaTools.h"
#include "math/HoInterval.h"
#include "math/IntDistribution.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/long/CreateGenome.h"
//#include "paths/long/EvalByReads.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/PairInfo.h"
#include "paths/long/SupportedHyperBasevector.h"

Bool LongHyper(const VecEFasta &correctede, const vec<pairing_info> &cpartner, SupportedHyperBasevector &shb,
               const long_heuristics &heur, const long_logging_control &log_control, const long_logging &logc,
               bool useOldLRPMethod) {

    // Choose K2.

    double K2frac_mult = 1.0;
    double K2frac_to_use = heur.K2frac * K2frac_mult;//.22
    int K2 = SelectK2(correctede, K2frac_to_use, logc, heur);
    K2 = Max(K2, heur.K2_FLOOR);
    if (heur.K2_FORCE >= 0) K2 = heur.K2_FORCE;
    double fudge_mult = K2frac_mult;

    // Setup for PERTURB_TRANSLATIONS.


    const bool bUseOrgReads = false;
    vecbasevector const &Bases = vecbasevector();
    QualVecVec const &Quals = QualVecVec();
    PairsManager const &Pairs = PairsManager();

    // Expand correctede.  Note that this is exponential and thus totally unsound.

    double eclock = WallClockTime();
    if (logc.STATUS_LOGGING) ReportPeakMem();
    if (logc.STATUS_LOGGING) std::cout << Date() << ": making hyper" << std::endl;
    vecbasevector correctedv;
    vec<triple<int, int, int> > origin;
    for (size_t id = 0; id < correctede.size(); id++) {
        vec<basevector> b;
        correctede[id].ExpandTo(b);
        for (int j = 0; j < b.isize(); j++) {
            correctedv.push_back_reserve(b[j]);
            origin.push(id, j, b.size());
        }
    }
    REPORT_TIME(eclock, "used in expansion");
    int64_t read_count = correctedv.size();
    if (heur.INJECT_REF && log_control.G != 0)
        correctedv.Append(*(log_control.G));

    // Make paths.

    double pclock = WallClockTime();
    HyperBasevector hb;
    HyperKmerPath h;
    vecKmerPath paths, paths_rc;
    if (logc.STATUS_LOGGING)
        std::cout << Date() << ": calling LongReadsToPaths" << std::endl;
    unsigned const COVERAGE = 50u;
    LongReadsToPaths(correctedv, K2, COVERAGE, &hb, &h, &paths, &paths_rc);
    REPORT_TIME(pclock, "used in pathing");

    // Trace reads through h.  Ignore reads that lie entirely on one edge.
    // Note that in fact each bona fide read can give rise to several 'reads'
    // via EFASTA and that each of these are traced here.  Note also that we
    // discard reads that do not have a perfect mapping to the graph, and this
    // could have bad consequences.

    double tclock = WallClockTime();
    vec<vec<int> > usu;
    vec<fix64_6> count_fw, count_rc;
    vec<std::pair<vec<int>, vec<int> > > pairs;
    vec<vec<pair_point> > pair_data;
    vec<vec<std::pair<fix64_6, int64_t> > > weights_fw_origin, weights_rc_origin;
    {
        if (logc.STATUS_LOGGING) std::cout << Date() << ": start tracing reads" << std::endl;
        vecKmerPath hpaths;
        vec<big_tagged_rpint> hpathsdb;
        for (int e = 0; e < h.EdgeObjectCount(); e++)
            hpaths.push_back_reserve(h.EdgeObject(e));
        CreateDatabase(hpaths, hpathsdb);

        // Track in parallel:
        // us = sequences in unipath graph
        // weight_fw = weight of sequence (forward orientation)
        // weight_rc = weight of sequence (reverse orientation)
        // denom_fw = reciprocal of forward weight (up to rounding)
        // denom_rc = reciprocal of reverse weight (up to rounding)
        // left = left extension of read by sequence
        // right = right extension of read by sequence
        // ids = original read id
        // fwplace = True if forward orientation

        vec<vec<int> > us;
        vec<fix64_6> weight_fw, weight_rc;
        vec<int> denom_fw, denom_rc;
        vec<int> left, right;
        vec<int> ids;
        vec<Bool> fwplace;

        // Go through the reads.

        vec<int> starts;
        for (int id = 0; id < read_count; id++) {
            starts.push_back(us.size());
            for (int pass = 1; pass <= 2; pass++) {
                const KmerPath &p = (pass == 1 ? paths[id] : paths_rc[id]);
                vec<triple<ho_interval, int, ho_interval> > M, M2;
                vec<int> u;
                int rpos = 0;
                for (int j = 0; j < p.NSegments(); j++) {
                    const KmerPathInterval &I = p.Segment(j);
                    vec<longlong> locs;
                    Contains(hpathsdb, I, locs);
                    for (int l = 0; l < locs.isize(); l++) {
                        const big_tagged_rpint &t = hpathsdb[locs[l]];
                        int hid = t.PathId();
                        if (hid < 0) continue;
                        longlong hpos = I.Start() - t.Start();
                        longlong start = Max(I.Start(), t.Start());
                        longlong stop = Min(I.Stop(), t.Stop());
                        longlong hstart = start - t.Start();
                        for (int r = 0; r < t.PathPos(); r++)
                            hstart += hpaths[hid].Segment(r).Length();
                        longlong hstop = hstart + stop - start;
                        longlong rstart = rpos + start - I.Start();
                        longlong rstop = rstart + stop - start;
                        M.push(ho_interval(rstart, rstop), hid,
                               ho_interval(hstart, hstop));
                    }
                    rpos += I.Length();
                }
                Bool bad = False;
                for (int i = 0; i < M.isize(); i++) {
                    int j;
                    for (j = i + 1; j < M.isize(); j++) {
                        if (M[j].first.Start() != M[j - 1].first.Stop() + 1)
                            break;
                        if (M[j].second != M[j - 1].second) break;
                        if (M[j].third.Start() != M[j - 1].third.Stop() + 1)
                            break;
                    }
                    u.push_back(M[i].second);
                    Bool incomplete = False;
                    if (i > 0 && M[i].third.Start() > 0) incomplete = True;
                    if (j < M.isize() && M[j - 1].third.Stop()
                                         != hpaths[M[i].second].KmerCount() - 1) {
                        incomplete = True;
                        bad = True;
                    }
                    if (logc.TRACE_READS && i == 0) {
                        std::cout << "\n" << (pass == 1 ? "+" : "-")
                                  << origin[id].first << "."
                                  << origin[id].second + 1 << "\n";
                    }
                    if (i == 0 && j == M.isize() && !incomplete) {
                        i = j - 1;
                        continue;
                    }
                    int last = (i == 0 ? -1 : M2.back().first.Stop());
                    if (M[i].first.Start() > last + 1) {
                        if (logc.TRACE_READS) {
                            std::cout << last + 1 << "-" << M[i].first.Start() - 1
                                      << " --> MISSING\n";
                        }
                        bad = True;
                    }
                    M2.push(ho_interval(M[i].first.Start(),
                                        M[j - 1].first.Stop()), M[i].second,
                            ho_interval(M[i].third.Start(),
                                        M[j - 1].third.Stop()));
                    if (logc.TRACE_READS) {
                        std::cout << M[i].first.Start() << "-"
                                  << M[j - 1].first.Stop() << " --> "
                                  << M[i].second << "." << M[i].third.Start()
                                  << "-" << M[j - 1].third.Stop()
                                  << (incomplete ? " [INCOMPLETE]" : "")
                                  << "\n";
                    }
                    if (j == M.isize() && M[j - 1].first.Stop()
                                          < p.KmerCount() - 1) {
                        if (logc.TRACE_READS) {
                            std::cout << M[j - 1].first.Stop() + 1 << "-"
                                      << p.KmerCount() - 1
                                      << " --> MISSING\n";
                        }
                        bad = True;
                    }
                    i = j - 1;
                }
                if (!bad && u.nonempty()) {
                    if (logc.TRACE_READS) {
                        std::cout << "u =";
                        for (int j = 0; j < u.isize(); j++)
                            std::cout << " " << u[j];
                        std::cout << "\n";
                    }
                    us.push_back(u);
                    if (pass == 1) {
                        weight_fw.push_back(fix64_6(1, origin[id].third));
                        weight_rc.push_back(0);
                        denom_fw.push_back(origin[id].third);
                        denom_rc.push_back(0);
                    } else {
                        weight_rc.push_back(fix64_6(1, origin[id].third));
                        weight_fw.push_back(0);
                        denom_rc.push_back(origin[id].third);
                        denom_fw.push_back(0);
                    }
                    left.push_back(M.front().third.Start());
                    right.push_back(hb.EdgeObject(u.back()).isize()
                                    - (M.back().third.Stop() + K2));
                    ids.push_back(origin[id].first);
                    fwplace.push_back(pass == 1);
                }
            }
        }
        starts.push_back(us.size());


        // Look up pairs.

        for (int id1 = 0; id1 < cpartner.isize(); id1++) {
            const pairing_info &p = cpartner[id1];
            if (p.Status() != 1) continue;
            // std::cout << "\npair\n";
            int id2 = p.Partner(), lib = p.LibId();
            for (int m1 = starts[id1]; m1 < starts[id1 + 1]; m1++)
                for (int m2 = starts[id2]; m2 < starts[id2 + 1]; m2++) {
                    if (fwplace[m1] == fwplace[m2]) continue;
                    int n1(m1), n2(m2);
                    if (!fwplace[m1]) std::swap(n1, n2);
                    fix64_6 w(1, (denom_fw[n1] + denom_rc[n1])
                                 * (denom_fw[n2] + denom_rc[n2]));
                    int trim = right[n1] + left[n2];
                    if (logc.PRINT_INITIAL_PAIRS) {
                        std::cout << printSeq(us[n1]) << " [w=" << weight_fw[n1]
                                  << "+" << weight_rc[n1] << ",trim=" << right[n1] << "]"
                                  << " --" << lib << "--> " << printSeq(us[n2])
                                  << " [w=" << weight_fw[n2] << "+" << weight_rc[n2]
                                  << ",trim=" << left[n2] << "]"
                                  << " (trim=" << trim << ", weight=" << w << ")\n";
                    }
                    pairs.push(us[n1], us[n2]);
                    vec<pair_point> p;
                    p.push(trim, w, lib);
                    pair_data.push(p);
                }
        }
        SortSync(pairs, pair_data);
        vec<std::pair<vec<int>, vec<int> > > pairs2;
        vec<vec<pair_point> > pair_data2;
        for (int i = 0; i < pairs.isize(); i++) {
            int j = pairs.NextDiff(i);
            pairs2.push_back(pairs[i]);
            vec<pair_point> p;
            for (int k = i; k < j; k++)
                p.push_back(pair_data[k][0]);
            Sort(p);
            pair_data2.push_back(p);
            i = j - 1;
        }
        pairs = pairs2;
        pair_data = pair_data2;

        // Generate main words.

        SortSync(us, weight_fw, weight_rc, ids);
        for (int i = 0; i < us.isize(); i++) {
            int j = us.NextDiff(i);
            usu.push_back(us[i]);
            fix64_6 cfw = 0, crc = 0;
            for (int k = i; k < j; k++) {
                cfw += weight_fw[k];
                crc += weight_rc[k];
            }
            count_fw.push_back(cfw);
            count_rc.push_back(crc);
            vec<std::pair<fix64_6, int64_t> > wfw, wrc;
            for (int k = i; k < j; k++) {
                if (weight_fw[k] > 0) wfw.push(weight_fw[k], ids[k]);
                if (weight_rc[k] > 0) wrc.push(weight_rc[k], ids[k]);
            }
            weights_fw_origin.push_back(wfw);
            weights_rc_origin.push_back(wrc);
            i = j - 1;
        }
        if (logc.TRACE_READS) {
            std::cout << "\ntraces:\n";
            for (int i = 0; i < usu.isize(); i++) {
                std::cout << "[" << i + 1 << "," << std::setiosflags(std::ios::fixed)
                          << std::setprecision(1) << count_fw[i] + count_rc[i]
                          << std::resetiosflags(std::ios::fixed) << "x]";
                for (int j = 0; j < usu[i].isize(); j++)
                    std::cout << " " << usu[i][j];
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        if (logc.STATUS_LOGGING)
            std::cout << Date() << ": read tracing complete" << std::endl;
    }

    // Get median read length and create SupportedHyperBasevector.

    vec<int> len;
    for (size_t i = 0; i < correctede.size(); i++) {
        int n = correctede[i].Length1();
        if (n == 0) continue;
        len.push_back(n);
    }
    Sort(len);
    IntDistribution read_length_dist;
    read_length_dist.from_hits(len);
    vec<int> inv;
    vecbasevector hbx;
    for (int e = 0; e < hb.EdgeObjectCount(); e++)
        hbx.push_back_reserve(hb.EdgeObject(e));
    REPORT_TIME(tclock, "used in tracing");
    double iclock = WallClockTime();
    UnibaseInvolution(hbx, inv);
    REPORT_TIME(iclock, "used finding involution");
    double z1clock = WallClockTime();
    if (usu.empty()) return False;
    shb = SupportedHyperBasevector(hb, inv, usu, count_fw, count_rc,
                                   weights_fw_origin, weights_rc_origin, pairs, pair_data,
                                   correctede.size(), read_length_dist, fudge_mult);
    REPORT_TIME(z1clock, "used after tracing");
    shb.FixWeights(logc);
    shb.TestValid(logc);

    return True;
}
