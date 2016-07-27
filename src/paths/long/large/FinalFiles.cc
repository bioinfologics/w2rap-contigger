///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
//#include "paths/long/large/DiscoStats.h"
#include "paths/long/large/FinalFiles.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"

// Build final assembly files, starting from the results of scaffolding.

void FinalFiles(const HyperBasevector &hb, const vec<int> &inv, const
ReadPathVec &paths, const vec<String> &subsam_names,
                const vec<int64_t> &subsam_starts, const String &work_dir, const String &prefix,
                const int MAX_CELL_PATHS, const int MAX_DEPTH,
                const vecbasevector &G) {
     // Write some assembly files.

     TestInvolution(hb, inv);

     // Make scaffold lines.

     vec<vec<vec<vec<int>>>> linesx;
     // std::cout << "finding scaffold lines, mem usage = "
     //      << ToStringAddCommas( MemUsageBytes( ) ) << std::endl;
     FindLines(hb, inv, linesx, MAX_CELL_PATHS, MAX_DEPTH);
     SortLines(linesx, hb, inv);
     BinaryWriter::writeFile(work_dir + "/" + prefix + ".lines", linesx);
     DumpLineFiles(linesx, hb, inv, paths, work_dir);
     {
          vec<vec<covcount>> covsx;
          ComputeCoverage(hb, inv, paths, linesx, subsam_starts, covsx);
          BinaryWriter::writeFile(work_dir + "/" + prefix + ".covs", covsx);
          vec<int> npairsx;
          vec<int> llensx;
          GetLineLengths(hb, linesx, llensx);
          GetLineNpairs(hb, inv, paths, linesx, npairsx);
          BinaryWriter::writeFile(work_dir + "/" + prefix + ".lines.npairs", npairsx);
          WriteLineStats(work_dir + "/" + prefix, linesx, llensx, npairsx, covsx);

          // Generate scaffolded assembly fasta file and other output files.
          //MakeFinalFasta(hb, inv, linesx, npairsx, covsx, work_dir, prefix);
     }

     // Report edge stats.

     String edges_kb_report;
     int64_t genome_size = 0, genome_size0 = 0;
     {
          vec<int> lens;
          for (int e = 0; e < hb.EdgeObjectCount(); e++) {
               if (hb.EdgeLengthBases(e) > 0)
                    lens.push_back(hb.EdgeLengthKmers(e));
               else lens.push_back(0);
          }
          Sort(lens);
          std::ostringstream out1;
          out1 << Mean(lens);
          std::cout << "\nassembly has " << lens.size() << " edges of mean length "
          << out1.str() << edges_kb_report << std::endl;

     }
     const int min_len = 1000;
     int64_t scaffoldN50 = LineN50(hb, linesx, min_len);
     vec<int> lens;
     GetLineLengths(hb, linesx, lens);
     int64_t total1 = 0, total10 = 0, total100 = 0;
     for (int i = 0; i < lens.isize(); i++) {
          if (lens[i] >= 1000) total1 += lens[i];
          if (lens[i] >= 10000) total10 += lens[i];
          if (lens[i] >= 100000) total100 += lens[i];
     }
     total1 /= 2;
     total10 /= 2;
     total100 /= 2;


     Ofstream(sout, work_dir + "/stats");
     {
          sout << "# " << prefix << " assembly statistics" << std::endl;
          sout << std::endl;
          sout << "N50: " << ToStringAddCommas(scaffoldN50) << std::endl;
          sout << "total bases in 1 kb+ sequences: " << ToStringAddCommas(total1) << std::endl;
          sout << "total bases in 10 kb+ sequences: " << ToStringAddCommas(total10) << std::endl;
          sout << "total bases in 100 kb+ sequences: " << ToStringAddCommas(total100) << std::endl;
          std::cout << "# " << prefix << " assembly statistics" << std::endl;
          std::cout << std::endl;
          std::cout << "total N50: " << ToStringAddCommas(scaffoldN50) << std::endl;
          std::cout << "total bases in 1 kb+ sequences: " << ToStringAddCommas(total1) << std::endl;
          std::cout << "total bases in 10 kb+ sequences: " << ToStringAddCommas(total10) << std::endl;
          std::cout << "total bases in 100 kb+ sequences: " << ToStringAddCommas(total100) << std::endl;


     }


     // Compute edge coverage.
     /*
     {
          Ofstream(xout, work_dir + "/a.counts");
          int ns = subsam_names.size();
          xout << "#";
          for (int j = 0; j < ns; j++)
               xout << " " << subsam_names[j];
          xout << "\n";
          vec<vec<int>> count(ns, vec<int>(hb.E(), 0));
          for (int64_t id = 0; id < (int64_t) paths.size(); id++) {
               int ss;
               for (ss = 0; ss < ns; ss++)
                    if (ss == ns - 1 || id < subsam_starts[ss + 1]) break;
               const ReadPath &p = paths[id];
               for (int j = 0; j < (int) p.size(); j++) {
                    count[ss][p[j]]++;
                    if (inv[p[j]] != p[j])
                         count[ss][inv[p[j]]]++;
               }
          }
          for (int e = 0; e < hb.E(); e++) {
               xout << e;
               for (int j = 0; j < ns; j++)
                    xout << " " << count[j][e];
               xout << "\n";
          }
     }*/
}
