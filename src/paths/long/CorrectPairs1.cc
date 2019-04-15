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
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "efasta/EfastaTools.h"
#include "paths/long/CorrectPairs1.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/ReadStack.h"
#include "util/w2rap_timers.h"
TIMELOG_CREATE_GLOBAL(CP1_Align);
TIMELOG_CREATE_GLOBAL(CP1_MakeStacks);
TIMELOG_CREATE_GLOBAL(CP1_Correct);


namespace { // open anonymous namespace

Bool cmp_ho_start_stop( const ho_interval& h1, const ho_interval& h2 )
{    if ( h1.Start( ) < h2.Start( ) ) return True;
     if ( h1.Start( ) > h2.Start( ) ) return False;
     return h1.Stop( ) > h2.Stop( );    }

} // close anonymous namespace

void CorrectPairs1( const int K, const int max_freq, vecbasevector& bases,
     QualVecVec& quals, const PairsManager& pairs, const vec<Bool>& to_edit,
     const vec<int>& trace_ids, const long_heuristics& heur, 
     //const long_logging_control& log_control, const long_logging& logc,
     VecEFasta& corrected )
{
     // Build alignments.
     TIMELOG_DECLARE_LOCAL(CP1_Align, Loop);

     TIMELOG_START_LOCAL(CP1_Align, Loop);
     FriendAligner faligns(bases,
                           heur.FF_MAKE_ALIGN_IMPL, K,
                           heur.FF_MAX_FREQ,
                           heur.FF_DOWN_SAMPLE, heur.FF_VERBOSITY);
     TIMELOG_STOP_LOCAL(CP1_Align, Loop);

     // Go through the reads.

     vec<int64_t> use;
     for (int64_t id1 = 0; id1 < (int64_t) bases.size(); id1++) {
          int64_t id2 = pairs.getPartnerID(id1);
          if (to_edit[id1] && to_edit[id2] && bases[id1].size() > 0 && id2 < id1)
               use.push_back(id1);
     }


     //int batch = (int64_t) use.size( ) / omp_get_max_threads( );
     //batch = Min( 100, Max( 1, batch ) );
     //#pragma omp parallel for schedule(dynamic, batch)
     for (int64_t id1x = 0; id1x < (int64_t) use.size(); id1x++) {
          TIMELOG_DECLARE_LOCAL(CP1_MakeStacks, Loop);
          TIMELOG_DECLARE_LOCAL(CP1_Correct, Loop);

          TIMELOG_START_LOCAL(CP1_MakeStacks, Loop);
          int64_t id1 = use[id1x];
          // Build stacks.
          std::ostringstream out;
          int64_t id1p = pairs.getPartnerID(id1);
          Bool trace = BinMember(trace_ids, id1) || BinMember(trace_ids, id1p);
          // Get aligns for this read
          Friends aligns;
          faligns.getAligns(id1, &aligns);
          readstack stack1(id1, aligns, 0, aligns.size(),
                           readstack::right_extended, bases, quals, pairs);
          Friends aligns_p;
          faligns.getAligns(id1p, &aligns_p);
          readstack stack2(id1p, aligns_p, 0, aligns_p.size(),
                           readstack::right_extended, bases, quals, pairs);
          if (stack1.Rows() > heur.MAX_STACK || stack2.Rows() > heur.MAX_STACK) {
               TIMELOG_STOP_LOCAL(CP1_MakeStacks, Loop);
               continue;
          }


          // Filter out low-quality reads.

          int total_bases = 0, total_qual = 0;
          for (int i = 0; i < stack1.Cols(); i++) {
               if (stack1.Def(0, i)) total_bases++;
               if (stack1.Qual(0, i) >= 2) total_qual += stack1.Qual(0, i);
          }
          for (int i = 0; i < stack2.Cols(); i++) {
               if (stack2.Def(0, i)) total_bases++;
               if (stack2.Qual(0, i) >= 2) total_qual += stack2.Qual(0, i);
          }
          double this_qual = double(total_qual) / double(total_bases);

          int bases_all = 0, total_all = 0;
          vec<int> ids_all;
          for (int pass = 1; pass <= 2; pass++) {
               const readstack &s = (pass == 1 ? stack1 : stack2);
               for (int j = 0; j < s.Rows(); j++)
                    ids_all.push_back(s.Id(j));
          }
          UniqueSort(ids_all);
          for (int m = 0; m < ids_all.isize(); m++) {
               int id = ids_all[m];
               for (int j = 0; j < (int) quals[id].size(); j++) {
                    bases_all++;
                    if (quals[id][j] >= 2) total_all += quals[id][j];
               }
          }
          if (bases_all == 0) bases_all++;
          double all_qual = double(total_all) / double(bases_all);


          if (all_qual - this_qual > heur.CP_MAX_QDIFF) {
               if (omp_get_thread_num() == 0) std::cout << out.str();
               TIMELOG_STOP_LOCAL(CP1_MakeStacks, Loop);
               continue;
          }

          // Remove friends having inadequate glue to the founder.

          vec<Bool> suspect;
          stack1.FlagNoise(suspect);
          stack1.Erase(suspect);
          stack2.FlagNoise(suspect);
          stack2.Erase(suspect);

          // Proceed.

          const int q_solid = 30;
          stack1.Raise1(0);
          stack1.MotifDiff(1, suspect);
          stack1.Erase(suspect);


          stack1.HighQualDiff(q_solid, 1, suspect);
          stack1.Erase(suspect);


          stack2.Raise1(0);
          stack2.MotifDiff(1, suspect);
          stack2.Erase(suspect);
          stack2.HighQualDiff(q_solid, 1, suspect);
          stack2.Erase(suspect);




          // Try to recruit other reads.  Computationally unsound.


          stack2.Reverse();
          basevector con1, con2;
          QualVec conq1, conq2;
          stack1.Consensus1(con1, conq1);
          stack2.Consensus1(con2, conq2);
          const int L = 20;

          // Check for single-base indel (second method).



          vec<int> offsets = GetOffsets1(stack1, stack2, 0, heur.DELTA_MIS);
          vec<int> final_offsets = offsets;
          TIMELOG_STOP_LOCAL(CP1_MakeStacks, Loop);

          // For each offset, create the merged stack associated to it.
          TIMELOG_START_LOCAL(CP1_Correct, Loop);
          vec<basevector> closures;
          vec<QualVec> closuresq;
          vec<int> closureso;

          for (int oj = 0; oj < offsets.isize(); oj++) {
               int minq_floor = (offsets.size() > 1 ? heur.CP_MINQ_FLOOR
                                                    : 5 /* Min( heur.CP_MINQ_FLOOR, 5 ) */ );
               int min_glue_floor = (offsets.size() > 1 ? heur.CP_MIN_GLUE
                                                        : Min(heur.CP_MIN_GLUE, 20));
               readstack stack(stack1);
               stack.Merge(stack2, offsets[oj]);
               stack.SortByPid(pairs.getPairID(id1), 0, stack1.Rows());
               stack.Unique();
               stack.Raise1(0), stack.Raise1(1);



               /*
               stack.MotifDiff(2,suspect);
               if ( suspect[0] || suspect[1] ) continue;
               stack.Erase(suspect);
               */

               stack.HighQualDiff(q_solid, 2, suspect);
               if (suspect[0] || suspect[1]) continue;
               stack.Erase(suspect);


               // stack.AddPartners( 32, 2, bases, quals, pairs );

               stack.PairWeak1(suspect);
               if (suspect[0] || suspect[1]) continue;
               stack.Erase(suspect);

               // stack.CleanColumns(2,suspect);
               // if ( !suspect[0] && !suspect[1] ) stack.Erase(suspect);

               int start, stop;
               for (start = 0; start < stack.Cols(); start++)
                    if (stack.Def(0, start)) break;
               for (stop = stack.Cols() - 1; stop >= 0; stop--)
                    if (stack.Def(1, stop)) break;
               stop++;
               Bool weird = (start >= stop);
               if (!weird) stack.Trim(start, stop);

               // Create and edit consensus.

               /*basevector con,con2;
               QualVec conq,conq2;
               stack.StrongConsensus1( con, conq, heur.CP_RAISE_ZERO );
               stack.StrongConsensus2( con2, conq2, heur.CP_RAISE_ZERO );
               if (con!=con2 || conq!=conq2) (FatalErr("StrongConsensus2 is WRONG!!!"));*/
               basevector con;
               QualVec conq;
               stack.StrongConsensus2(con, conq, heur.CP_RAISE_ZERO);


               const int protected_bases = 10;
               const int q_to_protect = 20;
               for (int j = 0; j < protected_bases; j++) {
                    if (j >= stack.Cols()) break;
                    if (stack.Def(0, j) && stack.Base(0, j) != con[j]
                        && stack.Qual(0, j) >= q_to_protect) {
                         con.Set(j, stack.Base(0, j));
                         conq[j] = stack.Qual(0, j);
                    }
               }
               for (int j = 0; j < protected_bases; j++) {
                    int jr = stack.Cols() - j - 1;
                    if (jr < 0) break;
                    if (stack.Def(1, jr) && stack.Base(1, jr) != con[jr]
                        && stack.Qual(1, jr) >= q_to_protect) {
                         con.Set(jr, stack.Base(1, jr));
                         conq[jr] = stack.Qual(1, jr);
                    }
               }

               for (int j = 0; j < con.isize(); j++) {
                    if (stack.Qual(0, j) >= 30 && stack.Base(0, j) != con[j])
                         conq[j] = 0;
                    if (stack.Qual(1, j) >= 30 && stack.Base(1, j) != con[j])
                         conq[j] = 0;
               }

               // Test for suspicious inconsistencies between the founder
               // and the consensus.

               for (int j = 0; j < con.isize(); j++)
                    for (int m = 0; m < 2; m++) {
                         if (!stack.Def(m, j) || stack.Base(m, j) == con[j])
                              continue;
                         const int flank = 5;
                         const int min_mult = 3;
                         if (j < flank || j + flank >= con.isize()) continue;
                         Bool mismatch = False;
                         for (int l = 0; l < flank; l++) {
                              if (stack.Base(m, j - l - 1) != con[j - l - 1])
                                   mismatch = True;
                              if (stack.Base(m, j + l + 1) != con[j + l + 1])
                                   mismatch = True;
                         }
                         if (mismatch) continue;
                         int mult = 0;
                         for (int r = 2; r < stack.Rows(); r++) {
                              Bool mismatch = False;
                              for (int p = j - flank; p <= j + flank; p++) {
                                   if (stack.Base(r, p) != stack.Base(m, p)) {
                                        mismatch = True;
                                        break;
                                   }
                              }
                              if (mismatch) continue;
                              if (++mult == min_mult) break;
                         }
                         if (mult == min_mult) conq[j] = 0;
                    }

               // Attempt to recover conflicted columns.

               vec<Bool> to_del(stack.Rows(), False);
               for (int j = 0; j < stack.Cols(); j++) {
                    if (conq[j] < minq_floor) {
                         const int qmin = 2;
                         const int qdelta = 10;
                         if (stack.Qual(0, j) < qmin && stack.Qual(1, j) < qmin)
                              continue;
                         if (stack.Qual(0, j) >= qmin && stack.Qual(1, j) >= qmin
                             && stack.Base(0, j) != stack.Base(1, j)
                             && Abs(stack.Qual(0, j) - stack.Qual(1, j)) < qdelta) { continue; }
                         char b;
                         if (stack.Qual(0, j) >= qmin
                             && stack.Qual(0, j) >= stack.Qual(1, j)) { b = stack.Base(0, j); }
                         else b = stack.Base(1, j);
                         for (int i = 2; i < stack.Rows(); i++) {
                              if (stack.Qual(i, j) >= qmin && stack.Base(i, j) != b) {
                                   to_del[i] = True;
                              }
                         }
                    }
               }
               stack.Erase(to_del);
               /*stack.StrongConsensus1( con, conq, heur.CP_RAISE_ZERO );
               stack.StrongConsensus2( con2, conq2, heur.CP_RAISE_ZERO );
               if (con!=con2 || conq!=conq2) (FatalErr("StrongConsensus2 is WRONG!!!"));*/
               stack.StrongConsensus2(con, conq, heur.CP_RAISE_ZERO);

               for (int j = 0; j < protected_bases; j++) {
                    if (j >= stack.Cols()) break;
                    if (stack.Def(0, j) && stack.Base(0, j) != con[j]
                        && stack.Qual(0, j) >= q_to_protect) {
                         con.Set(j, stack.Base(0, j));
                         conq[j] = stack.Qual(0, j);
                    }
               }
               for (int j = 0; j < protected_bases; j++) {
                    int jr = stack.Cols() - j - 1;
                    if (jr < 0) break;
                    if (stack.Def(1, jr) && stack.Base(1, jr) != con[jr]
                        && stack.Qual(1, jr) >= q_to_protect) {
                         con.Set(jr, stack.Base(1, jr));
                         conq[jr] = stack.Qual(1, jr);
                    }
               }

               const int qfloor = 20;

               // Probably doesn't do anything:

               Bool yes1 = False, yes2 = False;
               for (int j = 0; j < stack.Cols(); j++) {
                    if (stack.Def(0, j)) yes1 = True;
                    if (stack.Def(1, j)) yes2 = True;
               }
               if (!yes1 || !yes2) continue;

               // Determine minimum consensus quality.

               int minq = 1000000000;
               for (int j = 0; j < con.isize(); j++)
                    minq = Min(minq, (int) conq[j]);

               // Check for glue.

               vec<ho_interval> agree;
               for (int i = 0; i < stack.Rows(); i++)
                    for (int j = 0; j < stack.Cols(); j++) {
                         if (stack.Base(i, j) != con[j]) continue;
                         int k;
                         for (k = j + 1; k < stack.Cols(); k++)
                              if (stack.Base(i, k) != con[k]) break;
                         if (k - j >= 40) agree.push(j, k);
                         j = k - 1;
                    }
               sort(agree.begin(), agree.end(), cmp_ho_start_stop);
               vec<Bool> to_delete(agree.size(), False);
               for (int i = 0; i < agree.isize(); i++) {
                    int j;
                    for (j = i + 1; j < agree.isize(); j++) { if (agree[j].Stop() > agree[i].Stop()) break; }
                    for (int k = i + 1; k < j; k++)
                         to_delete[k] = True;
                    i = j - 1;
               }
               EraseIf(agree, to_delete);
               int min_glue;
               if (agree.empty() || agree[0].Start() > 0) min_glue = 0;
               else {
                    min_glue = agree[0].Length();
                    int stop = agree[0].Stop();
                    for (int i = 1; i < agree.isize(); i++) {
                         if (agree[i].Stop() > stop) {
                              min_glue = Min(min_glue, stop - agree[i].Start());
                              stop = agree[i].Stop();
                         }
                    }
                    if (stop < con.isize()) min_glue = 0;
               }
               if (minq >= minq_floor && min_glue >= min_glue_floor) {
                    closures.push_back(con), closuresq.push_back(conq);
                    closureso.push_back(offsets[oj]);
               }
          }

          // Save closures.  Currently we only save a longest one.
          // (NO, CHANGED.)

          if (closures.nonempty()) {
               if (heur.CP_CONDENSE_HOMOPOLYMERS) {
                    efasta eclosure(closures);
                    if (eclosure.AmbEventCount() == 1) {
                         String mid = eclosure.Between("{", "}");
                         mid.GlobalReplaceBy(",", "");
                         Bool homopolymer = True;
                         for (int j = 1; j < mid.isize(); j++)
                              if (mid[j] != mid[0]) homopolymer = False;
                         if (homopolymer) {
                              corrected[id1] = eclosure;
                              eclosure.ReverseComplement();
                              corrected[id1p] = eclosure;
                              continue;
                         }
                    }
               }

               int mc = 1000000000;
               for (int i = 0; i < closures.isize(); i++)
                    mc = Min(mc, closures[i].isize());
               basevector left(mc), right(mc);

               for (int j = 0; j < mc; j++) {
                    vec<char> seen;
                    for (int i = 0; i < closures.isize(); i++)
                         seen.push_back(closures[i][j]);
                    UniqueSort(seen);
                    if (seen.solo()) left.Set(j, seen[0]);
                    else {
                         left.resize(j);
                         break;
                    }
               }
               for (int j = 0; j < mc; j++) {
                    vec<char> seen;
                    for (int i = 0; i < closures.isize(); i++) {
                         int jb = closures[i].isize() - j - 1;
                         seen.push_back(closures[i][jb]);
                    }
                    UniqueSort(seen);
                    if (seen.solo()) right.Set(mc - j - 1, seen[0]);
                    else {
                         right.SetToSubOf(right, mc - j, j);
                         break;
                    }
               }

               corrected[id1] = left;
               if (left != right) {
                    right.ReverseComplement();
                    corrected[id1p] = right;
               }
          }
          TIMELOG_STOP_LOCAL(CP1_Correct, Loop);
     }

}
