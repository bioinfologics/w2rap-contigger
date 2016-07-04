///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "TokenizeString.h"
#include "bam/ReadBAM.h"
#include "feudal/ObjectManager.h"
#include "feudal/PQVec.h"
#include "math/HoInterval.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/large/ExtractReads.h"
#include "paths/long/large/ReadNameLookup.h"

class rs_meta { // read set meta info
     public:
     String type;
     String sample;
     String lib;
     double frac;
     rs_meta( ) : type("frag"), sample("C"), lib("1"), frac(1) { }
     friend std::ostream& operator<<( std::ostream& out, const rs_meta& m )
     {    return out << "type=" << m.type << ",sample=" << m.sample << ",lib=" 
               << m.lib << ",frac=" << m.frac;    }
};

void GetAmbInt( const vecbitvector& amb, vec< std::pair<int,ho_interval> >& amb_int );

void GetCannedReferenceSequences( const String& sample, const String& species,
     const String& work_dir );

void ExtractReads( String reads, const String& work_dir, vec<String>& subsam_names,
     vec<int64_t>& subsam_starts, vecbvec* pReads, VecPQVec* quals )
{
     double lclock = WallClockTime( );


     // Parse, check and load files.  This is the main extraction path.

     std::cout << Date() << ": finding input files" << std::endl;
     reads.GlobalReplaceBy(" ", "");
     vec<vec<String> > infiles;
     vec<vec<String> > infiles_rn;
     vecbasevector &xbases = (*pReads);
     VecPQVec &xquals = (*quals);
     vecString xnames;
     vec<vec<std::pair<int, int> > > infiles_pairs;
     vec<rs_meta> infiles_meta;

     // Define read groups.
     vec<String> groups;
     String line;

     Tokenize(reads, '+', groups);

     // Extract metainfo.
     //[GONZA] This is useless for us, never used metainformation
     infiles_meta.resize(groups.size());
     for (int g = 0; g < groups.isize(); g++) {
          String meta;
          if (groups[g].Contains("::")) {
               meta = groups[g].Before("::");
               groups[g] = groups[g].After("::");
          }
          vec<String> parts;
          Tokenize(meta, ',', parts);
          rs_meta &x = infiles_meta[g];
          for (int i = 0; i < parts.isize(); i++) {
               if (!parts[i].Contains(":")) {
                    std::cout << "\nIllegal metainfo specification " << meta
                    << ".\nEach metainfo field must be of the form "
                    << "arg:value.\n" << std::endl;
                    Scram(1);
               }
               String arg = parts[i].Before(":"), val = parts[i].After(":");
               if (arg == "type") {
                    if (val == "frag" || val == "long" || val == "jump")
                         x.type = val;
                    else {
                         std::cout << "\nUnrecognized type " << val
                         << " in metainfo specification " << meta
                         << ".\n" << std::endl;
                         Scram(1);
                    }
                    if (val != "frag") {
                         std::cout << "\nCurrently only type frag is implemented."
                         << "\n" << std::endl;
                         Scram(1);
                    }
               }
               else if (arg == "sample") x.sample = val;
               else if (arg == "lib") x.lib = val;
               else if (arg == "frac") {
                    if (!val.IsDouble() || val.Double() <= 0
                        || val.Double() > 1) {
                         std::cout << "\nIllegal value " << val << " for frac in "
                         << "metainfo specification " << meta
                         << ".\n" << std::endl;
                         Scram(1);
                    }
                    x.frac = val.Double();
               }
               else {
                    std::cout << "\nIllegal argument " << arg << " in "
                    << "metainfo specification " << meta
                    << ".\n" << std::endl;
                    Scram(1);
               }
          }
     }

     // Sort groups by sample.
     //[GONZA] Useless
     subsam_names.clear();
     vec<String> ssn;
     for (int g = 0; g < groups.isize(); g++) {
          subsam_names.push_back(infiles_meta[g].sample);
          if (!Member(ssn, subsam_names.back()))
               ssn.push_back(subsam_names.back());
     }
     vec<int> ssi(subsam_names.size());
     for (int i = 0; i < subsam_names.isize(); i++)
          ssi[i] = Position(ssn, subsam_names[i]);
     SortSync(ssi, subsam_names, groups, infiles_meta);
     subsam_starts.resize_and_set(subsam_names.size(), 0);

     // Create file list.

     infiles.resize(groups.size());
     infiles_rn.resize(groups.size());
     infiles_pairs.resize(groups.size());
     for (int g = 0; g < groups.isize(); g++) {
          String s = groups[g];

          // Parse data.

          vec<String> fns;
          int bcount = 0;
          String b;
          for (int i = 0; i < s.isize(); i++) {
               if (s[i] == '{') {
                    b.push_back(s[i]);
                    bcount++;
               }
               else if (s[i] == '}') {
                    b.push_back(s[i]);
                    bcount--;
               }
               else if (s[i] == ',' && bcount == 0) {
                    fns.push_back(b);
                    b.clear();
               }
               else b.push_back(s[i]);
          }
          fns.push_back(b);
          for (int i = 0; i < fns.isize(); i++) {
               String f = fns[i];
               vec<String> fs;
               int ok = Glob(f, fs);
               if (ok != 0) {
                    std::cout << "\nFailed to glob " << f << ".\n"
                    << "This means that it does not correspond to a "
                    << "file or files according to the\n"
                    << "rules for globbing.  "
                    << "Please see the DISCOVAR de novo manual.\n"
                    << std::endl;
                    Scram(1);
               }
               infiles[g].append(fs);
          }
     }

     // Check that files are OK.

     for (int g = 0; g < groups.isize(); g++)
          for (int j = 0; j < infiles[g].isize(); j++) {
               String fn = infiles[g][j];
               if (!IsRegularFile(fn)) {
                    std::cout << "\nFailed to find file " << fn << ".\n" << std::endl;
                    Scram(1);
               }
               vec<String> suf = {".bam", ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fastb"};
               Bool ok = False;
               for (auto s : suf)
                    if (fn.Contains(s, -1)) ok = True;
               if (!ok) {
                    std::cout << "\nInput file " << fn << " has unsupported type.\n"
                    << std::endl;
                    Scram(1);
               }
          }
     vec<String> flat;
     for (int g = 0; g < groups.isize(); g++)
          flat.append(infiles[g]);
     Sort(flat);
     for (int i = 1; i < flat.isize(); i++) {
          if (flat[i] == flat[i - 1] || flat[i] == flat[i - 1] + ".gz"
              || flat[i] + ".gz" == flat[i - 1]) {
               std::cout << "\nFile " << flat[i]
               << " appears more than once in your "
               << "READS specification.\n" << std::endl;
               Scram(1);
          }
     }

     // Get first readnames for fastq files, and sort by them.

     for (int g = 0; g < groups.isize(); g++) {
          infiles_rn[g].resize(infiles[g].size());
          for (int j = 0; j < infiles[g].isize(); j++) {
               String fn = infiles[g][j];
               vec<String> suf = {".fastq", ".fastq.gz", ".fq", ".fq.gz"};
               Bool fq = False;
               for (auto s : suf)
                    if (fn.Contains(s, -1)) fq = True;
               if (!fq) continue;
               String command = "cat " + fn + " | head -1";
               if (fn.Contains(".gz", -1)) command = "z" + command;
               fast_pipe_ifstream in(command);
               getline(in, line);
               if (!line.Contains("@", 0) || line.size() == 1
                   || (line[1] == ' ' || line[1] == '/')) {
                    std::cout << "\nSomething is wrong with the first line of your "
                    << "fastq file " << fn << "\n" << std::endl;
                    Scram(1);
               }
               int p = 0;
               for (p = 0; p < line.isize(); p++)
                    if (line[p] == ' ' || line[p] == '/') break;
               infiles_rn[g][j] = line.substr(1, p - 1);
          }
          SortSync(infiles_rn[g], infiles[g]);
          for (int j = 0; j < infiles_rn[g].isize(); j++) {
               int k = infiles_rn[g].NextDiff(j);
               if (k - j > 2 && infiles_rn[g][j] != "") {
                    std::cout << "\nThere are more than two fastq files that "
                    << "start with the read name " << infiles_rn[g][j]
                    << ":\n";
                    for (int l = j; l < k; l++)
                         std::cout << "[" << l - j + 1 << "] " << infiles[g][l] << std::endl;
                    std::cout << "Therefore it's not clear how to pair the "
                    << "files.\n\n";
                    Scram(1);
               }
          }
     }

     // Precheck fastb files.

     for (int g = 0; g < groups.isize(); g++) {
          if (infiles_meta[g].type != "long") {
               for (int j = 0; j < infiles[g].isize(); j++) {
                    String fn = infiles[g][j];
                    if (!fn.Contains(".fastb", -1)) continue;
                    int64_t B = MastervecFileObjectCount(fn), Q;
                    if (B % 2 != 0) {
                         std::cout << "\nThe file " << fn << " should be interlaced "
                         << "and hence have an even number\n"
                         << "of entries.  It does not.\n" << std::endl;
                         Scram(1);
                    }
                    String fn2b = fn.RevBefore(".fastb") + ".qualb";
                    String fn2p = fn.RevBefore(".fastb") + ".qualp";
                    String both;
                    if (IsRegularFile(fn2b)) {
                         Q = MastervecFileObjectCount(fn2b);
                         both = fn.Before(".fastb") + ".{fastb,qualb}";
                    }
                    else if (IsRegularFile(fn2p)) {
                         Q = MastervecFileObjectCount(fn2p);
                         both = fn.Before(".fastb") + ".{fastb,qualp}";
                    }
                    else {
                         std::cout << "\nFor the file " << fn << ", there is no "
                         << "matching qualb or qualp file.\n" << std::endl;
                         Scram(1);
                    }
                    if (B != Q) {
                         std::cout << "\nSee unequal file size for " << both
                         << ".\n" << std::endl;
                         Scram(1);
                    }
               }
          }
     }

     // Read the files.

     int nfiles = 0;
     for (int g = 0; g < groups.isize(); g++)
          nfiles += infiles[g].size();
     std::cout << Date() << ": reading " << nfiles
     << " files (which may take a while)" << std::endl;
     for (int g = 0; g < groups.isize(); g++) {
          if (g > 0 && subsam_names[g] != subsam_names[g - 1])
               subsam_starts[g] = xbases.size();
          for (int j = 0; j < infiles[g].isize(); j++) {
               String fn = infiles[g][j];
               int64_t N0 = xbases.size();

               // Parse bam files.

               if (fn.Contains(".bam", -1)) {
                    bool const UNIQUIFY_NAMES = true;
                    vecString *pxnames = 0;
                    BAMReader bamReader(False /*USE_PF_ONLY*/, UNIQUIFY_NAMES,
                                        infiles_meta[g].frac, long(-1/*READS_TO_USE*/));
                    bamReader.readBAM(
                            fn, &xbases, &xquals, pxnames);
               }

                    // Parse fastb/qualb/qualp files.

               else if (fn.Contains(".fastb", -1)) {
                    String fn2b = fn.RevBefore(".fastb") + ".qualb";
                    String fn2p = fn.RevBefore(".fastb") + ".qualp";
                    if (IsRegularFile(fn2b)) {
                         xbases.ReadAll(fn, True);
                         vecqualvector q;
                         q.ReadAll(fn2b);
                         convertAppendParallel(q.begin(), q.end(), xquals);
                         infiles[g][j]
                                 = fn.Before(".fastb") + ".{fastb,qualb}";
                    }
                    else if (IsRegularFile(fn2p)) {
                         xbases.ReadAll(fn, True);
                         xquals.ReadAll(fn2p, True);
                         infiles[g][j] = fn.Before(".fastb")
                                         + ".{fastb,qualp}";
                    }
                    double frac = infiles_meta[g].frac;
                    if (frac < 1) {
                         int64_t total = 0, taken = 0;
                         Bool skip_next = False;
                         int64_t pos = N0;
                         for (int64_t i = N0; i < (int64_t) xbases.size(); i++) {
                              total++;
                              if (skip_next) {
                                   skip_next = False;
                                   continue;
                              }
                              if (total % 2 == 1
                                  && double(taken) / double(total) > frac) {
                                   skip_next = True;
                                   continue;
                              }
                              taken++;
                              if (pos < i) {
                                   xbases[pos] = xbases[i];
                                   xquals[pos] = xquals[i];
                              }
                              pos++;
                         }
                         xbases.resize(pos), xquals.resize(pos);
                    }
               }

                    // Parse paired fastq files.

               else if (infiles_rn[g][j] != ""
                        && j < infiles_rn[g].isize() - 1
                        && infiles_rn[g][j] == infiles_rn[g][j + 1]) {
                    infiles_pairs[g].push(j, j + 1);
                    const String &fn1 = infiles[g][j], &fn2 = infiles[g][j + 1];
                    String command1 = "cat " + fn1, command2 = "cat " + fn2;
                    if (fn1.Contains(".gz", -1)) command1 = "z" + command1;
                    if (fn2.Contains(".gz", -1)) command2 = "z" + command2;
                    fast_pipe_ifstream in1(command1), in2(command2);
                    String line1, line2;
                    int64_t total = 0, taken = 0;
                    double frac = infiles_meta[g].frac;

                    // Buffer for quality score compression in batches.

                    const int qbmax = 10000000;
                    vec<qvec> qualsbuf;
                    MempoolOwner<char> alloc;
                    for (int i = 0; i < qbmax; i++)
                         qualsbuf.emplace_back(alloc);
                    int qbcount = 0;

                    // Go through the input files.

                    basevector b1, b2;
                    while (1) {
                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() && in2.fail()) break;
                         if ((in1.fail() && !in2.fail())
                             || (in2.fail() && !in1.fail())) {
                              std::cout << "\nThe files " << fn1 << " and " << fn2
                              << " appear to be paired, yet have "
                              << "different numbers of records.\n" << std::endl;
                              Scram(1);
                         }

                         // Fetch bases.  Turn Ns into As.

                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() || in2.fail()) {
                              std::cout << "\nSee incomplete record in " << fn1
                              << " or " << fn2 << ".\n" << std::endl;
                              Scram(1);
                         }
                         for (int i = 0; i < line1.isize(); i++)
                              if (line1[i] == 'N') line1[i] = 'A';
                         for (int i = 0; i < line2.isize(); i++)
                              if (line2[i] == 'N') line2[i] = 'A';
                         b1.SetFromString(line1);
                         b2.SetFromString(line2);

                         // Skip line.

                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() || in2.fail()) {
                              std::cout << "\nSee incomplete record in " << fn1
                              << " or " << fn2 << ".\n" << std::endl;
                              Scram(1);
                         }

                         // Fetch quals.

                         getline(in1, line1), getline(in2, line2);
                         if (in1.fail() || in2.fail()) {
                              std::cout << "\nSee incomplete record in " << fn1
                              << " or " << fn2 << ".\n" << std::endl;
                              Scram(1);
                         }
                         if (b1.size() != line1.size()
                             || b2.size() != line2.size()) {
                              std::cout << "\n1: " << b1.size() << " bases "
                              << ", " << line1.size() << " quals" << std::endl;
                              std::cout << "2: " << b2.size() << " bases "
                              << ", " << line2.size() << " quals" << std::endl;
                              std::cout << "See inconsistent base/quality lengths "
                              << "in " << fn1 << " or " << fn2 << std::endl;
                              Scram(1);
                         }

                         // Check frac.

                         if (frac < 1) {
                              total++;
                              if (double(taken) / double(total) > frac)
                                   continue;
                              taken++;
                         }

                         // Save.

                         qvec &q1 = qualsbuf[qbcount++];
                         qvec &q2 = qualsbuf[qbcount++];
                         q1.resize(line1.size()), q2.resize(line2.size());
                         if (qbcount == qbmax) {
                              convertAppendParallel(qualsbuf.begin(),
                                                    qualsbuf.begin() + qbcount, xquals);
                              qbcount = 0;
                         }
                         for (int i = 0; i < line1.isize(); i++)
                              q1[i] = line1[i] - 33;
                         for (int i = 0; i < line2.isize(); i++)
                              q2[i] = line2[i] - 33;
                         xbases.push_back(b1), xbases.push_back(b2);
                    }
                    convertAppendParallel(qualsbuf.begin(),
                                          qualsbuf.begin() + qbcount, xquals);
                    j++;
               }

                    // Parse unpaired fastq files.

               else if (infiles_rn[g][j] != "") {
                    vecqualvector Q;
                    const String &fn = infiles[g][j];
                    String command = "cat " + fn;
                    if (fn.Contains(".gz", -1)) command = "z" + command;
                    fast_pipe_ifstream in(command);
                    int64_t total = 0, taken = 0;
                    double frac = infiles_meta[g].frac;
                    Bool skip_next = False;
                    while (1) {
                         getline(in, line);
                         if (in.fail()) break;

                         // Fetch bases.  Turn Ns into As.

                         getline(in, line);
                         if (in.fail()) {
                              std::cout << "\nSee incomplete record in " << fn
                              << ".\n" << std::endl;
                              Scram(1);
                         }
                         for (int i = 0; i < line.isize(); i++)
                              if (line[i] == 'N') line[i] = 'A';
                         basevector b(line);

                         // Skip line.

                         getline(in, line);
                         if (in.fail()) {
                              std::cout << "\nSee incomplete record in " << fn
                              << ".\n" << std::endl;
                              Scram(1);
                         }

                         // Fetch quals.

                         getline(in, line);
                         if (in.fail()) {
                              std::cout << "\nSee incomplete record in " << fn
                              << ".\n" << std::endl;
                              Scram(1);
                         }
                         if (b.size() != line.size()) {
                              std::cout << "\nSee " << b.size() << " bases "
                              << ", " << line.size() << " quals" << std::endl;
                              std::cout << "See inconsistent base/quality lengths "
                              << "in " << fn << ".\n" << std::endl;
                              Scram(1);
                         }

                         // Check frac.

                         if (frac < 1) {
                              total++;
                              if (skip_next) {
                                   skip_next = False;
                                   continue;
                              }
                              if (total % 2 == 1
                                  && double(taken) / double(total) > frac) {
                                   skip_next = True;
                                   continue;
                              }
                              taken++;
                         }

                         // Save.

                         qualvector q(line.size());
                         for (int i = 0; i < line.isize(); i++)
                              q[i] = line[i] - 33;
                         xbases.push_back(b);
                         Q.push_back(q);
                    }

                    // Check sanity and compress.

                    if (infiles_meta[g].type != "long" && Q.size() % 2 != 0) {
                         std::cout << "\nThe file\n" << fn
                         << "\nshould be interlaced "
                         << "and hence have an even number of entries."
                         << "  It does not.\n" << std::endl;
                         Scram(1);
                    }
                    convertAppendParallel(
                            Q.begin(), Q.end(), xquals);
               }
          }
     }

     // Generate file list.

     vec<String> f1;
     vec<std::pair<String, String> > f2;
     vec<rs_meta> f1_meta, f2_meta;
     for (int g = 0; g < infiles.isize(); g++) {
          vec<Bool> used(infiles[g].size(), False);
          for (int j = 0; j < infiles_pairs[g].isize(); j++) {
               int x1 = infiles_pairs[g][j].first;
               int x2 = infiles_pairs[g][j].second;
               f2.push(infiles[g][x1], infiles[g][x2]);
               f2_meta.push_back(infiles_meta[g]);
               used[x1] = used[x2] = True;
          }
          for (int j = 0; j < infiles[g].isize(); j++) {
               if (!used[j]) {
                    f1.push_back(infiles[g][j]);
                    f1_meta.push_back(infiles_meta[g]);
               }
          }
     }
     std::cout << "\nINPUT FILES:\n";
     std::ostringstream iout;
     for (int j = 0; j < f1.isize(); j++)
          iout << "[" << j + 1 << "," << f1_meta[j] << "]  " << f1[j] << std::endl;
     for (int j = 0; j < f2.isize(); j++) {
          iout << "[" << f1.isize() + j + 1 << "a"
          << "," << f2_meta[j] << "] " << f2[j].first << std::endl;
          iout << "[" << f1.isize() + j + 1 << "b"
          << "," << f2_meta[j] << "] " << f2[j].second << std::endl;
     }
     {
          Ofstream(ioutx, work_dir + "/input_files");
          ioutx << iout.str();
     }
     std::cout << iout.str() << "\n";

     // Fix subsam.

     int scount = 0;
     for (int i = 0; i < subsam_names.isize(); i++) {
          if (i > 0) {
               if (subsam_names[i] != subsam_names[i - 1]) {
                    subsam_names[++scount] = subsam_names[i];
                    subsam_starts[scount] = subsam_starts[i];
               }
          }
     }
     subsam_names.resize(scount + 1), subsam_starts.resize(scount + 1);
     std::cout << Date() << ": found " << subsam_names.size()
     << " samples" << std::endl;
     std::cout << Date() << ": starts = " << printSeq(subsam_starts) << std::endl;
     for (int i = 0; i < subsam_starts.isize(); i++) {
          if ((i < subsam_starts.isize() - 1
               && subsam_starts[i] == subsam_starts[i + 1])
              || (i == subsam_starts.isize() - 1
                  && subsam_starts[i] == (int64_t) xbases.size())) {
               std::cout << "\nWARNING: It looks like you've got zero reads "
               << "for sample " << subsam_names[i] << ".\n";
               std::cout << "One way this could happen, for example, is if you "
               << "provided a BAM file\nin which none of the reads were "
               << "paired." << std::endl;
          }
     }

     // Work around bug.

     /*
     for ( int64_t i = 0; i < (int64_t) xbases.size( ); i++ )
     {    if ( xbases[i].size( ) < 60 )
          {    std::cout << "See read of length " << xbases[i].size( ) << ".  ";
               std::cout << "For now, all reads must have length >= 60.\n" << std::endl;
               Scram(1);    }    }
     */

     // Check for no reads.

     if (xbases.size() == 0) {
          std::cout << "\nThere are no reads." << std::endl;
          std::cout << "Assembly cannot proceed, but please have an A1 day!"
          << std::endl << std::endl;
          Scram(1);
     }

     // Save files.


     if (xnames.size() > 0) {
          xnames.WriteAll(work_dir + "/data/frag_reads_orig.names");
          readname_lookup look(xnames);
          BinaryWriter::writeFile(
                  work_dir + "/data/frag_reads_orig.names.idx", look);
     }

     // Report stats.

     std::cout << Date() << ": data extraction complete"
     #ifdef __linux
     << ", peak = " << PeakMemUsageGBString( )
     #endif
     << std::endl;
     std::cout << TimeSince(lclock) << " used extracting reads" << std::endl;



}

void GetAmbInt( const vecbitvector& amb, vec< std::pair<int,ho_interval> >& amb_int )
{    for ( int g = 0; g < (int) amb.size( ); g++ )
     {    for ( int i = 0; i < (int) amb[g].size( ); i++ )
          {    if ( !amb[g][i] ) continue;
               int j;
               for ( j = i + 1; j < (int) amb[g].size( ); j++ )
                    if ( !amb[g][j] ) break;
               amb_int.push( g, ho_interval( i, j ) );
               i = j - 1;    }    }    }
