///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <fenv.h>
#include "paths/long/EMEC3.h"
#include "VecUtilities.h"


namespace {		// one-off local utilities


/*
 * log_domain_norm
 *
 * normalize a vector of log-probabilities to sum to zero, so that their
 * associated probabilities sum to one.  This is useful for normalizing *before*
 * exponentiating.
 *
 * requires a vector-like object (with value_type, size(), and operator[] )
 *
 */
template <class T>
void log_domain_norm( T& p ) {

    typename T::value_type maxv = p[0];		// calculate the max value
    for ( size_t i = 1; i < p.size(); ++i )
        if ( p[i] > maxv )
            maxv = p[i];

    typename T::value_type sum = 0;		// calculate a sum by dividing out the max value and then
    for ( size_t i = 0; i < p.size(); ++i )	// multiplying it back in (really subtracting and adding )
        sum += exp( p[i] - maxv );
    sum = maxv + log( sum );

    for ( size_t i = 0; i < p.size(); ++i )	// subtract to "divide" by the sum and normalize
        p[i] = p[i] - sum;
}

/*
 * print our base vector which might contain spaces
 */
std::ostream& operator<<(std::ostream& os, const _BaseVec& v ) {
    for (size_t i = 0; i < v.size(); ++i ) {
        os << ( v[i] == 32 ? ' ' : Base::val2Char(v[i]) );
        if ( (i+1) % 80 == 0 )
            os << std::endl;
    }
    return os;
}


/*
 * printq - prints (extended) quality scores with dashes for negative numbers,
 * digits for 0-9, a-z for 10-35, and A-Z for 36-61
 */
std::ostream& printq(std::ostream& os, const _QualVec& v ) {
    for (size_t i = 0; i < v.size(); ++i ) {
        if ( v[i] < 0 )
            os << "-";
        else if ( v[i] < 10 )
            os << v[i];
        else if ( v[i] < 10+26 )
            os << static_cast<char>( v[i] + 'a' - 10 );
        else if ( v[i] < 10+26+26 )
            os << static_cast<char>( v[i] + 'A' - 10 - 26);
        else
            os << "X";
        if ( (i+1 ) % 80 == 0 )
            os << std::endl;
    }
    return os;
}


/*
 * Precomputed - table of log-probabilities from quality scores:
 *
 * Given quality score Q, return:
 *  q = Pr{Error} = 10^(-Q/10) 		Q2p
 *  p = 1-q = Pr{!Error}		Q2q
 *  log q and log p			Q2logq and Q2logp
 *
 *  This drastically speeds up things.
 */

class Precomputed {
  private:
    static const int max_Q = 55;

    Precomputed() : _Q2q(max_Q+1), _Q2p(max_Q+1), _Q2logq(max_Q+1), _Q2logp(max_Q+1)  {
        for ( int i = 0; i <= max_Q; ++i ) {
            double Q = static_cast<double>(i);
            double q = pow(10, -1* static_cast<double>(Q) / 10);
            _Q2q[i] = q;
            _Q2p[i] = 1.0 - q;
            _Q2logq[i] = log(_Q2q[i]);
            _Q2logp[i] = log(_Q2p[i]);	// there will be a -Inf in zero here... good.
        }
    };

    std::vector<double> _Q2q;
    std::vector<double> _Q2p;
    std::vector<double> _Q2logq;
    std::vector<double> _Q2logp;


  public:

    static double Q2p(int Q) {
        ForceAssertGe(Q, 0);
        return me._Q2p[Q];
    };
    static double Q2q(int Q) {
        ForceAssertGe(Q, 0);
        return me._Q2q[Q];
    };
    static double Q2logp(int Q) {
        ForceAssertGe(Q, 0);
        return me._Q2logp[Q];
    };
    static double Q2logq(int Q) {
        ForceAssertGe(Q, 0);
        return me._Q2logq[Q];
    };
    static Precomputed me;
};

Precomputed Precomputed::me;		// static member houses the table



/*
 * imax - find the index of the maximum element of a vector-like thing
 */
template <class T>
size_t imax( T& v ) {
    size_t max_i = 0;
    typename T::value_type max_v = v[0];
    for ( size_t i = 1; i < v.size(); ++i )
        if ( v[i] > max_v ) {
            max_v = v[i];
            max_i = i;
        }
    return max_i;
}



/*
 * init_prior - intialize the prior based on a given read/quals pairs.
 *
 * This is usually used to generate a probabilitistic (one probability per base) prior based on a
 * read and its associated quality scores.
 *
 */
void init_prior( const _BaseVec& read, const _QualVec& quals, std::vector<BaseProb>& prior, bool debug = false ) {
    size_t read_len = read.size();
    prior.resize(read_len, BaseProb(4,0.));

    if ( debug )
        for ( size_t i = 0; i < prior.size(); ++i )
            std::cout << "read,qual[" << i << "]=" << Base::val2Char(read[i]) << "," << quals[i] << std::endl;

    for ( size_t i = 0; i < read_len; ++i ) {
        unsigned char call = read[i];
        unsigned char qual = quals[i];	// NIW -- attention -- why is this u_char?

        if ( qual == 0 ) qual = 30;		// special case: Q0 has been previously corrected

        ForceAssertGt(qual, 0);
        ForceAssertLt(qual, 55);

        double perror = ( qual == 1 || qual == 2 ) ? 0.66 : pow(10, -1.* static_cast<double>(qual) / 10.);  // special case: Q1 and Q2 are unknown -- for now we use Pr{Error}=0.66

        prior[i][0] = prior[i][1] = prior[i][2] = prior[i][3] = perror / 3;		// 3 bases get perror/3, 1 base gets 1-perror.
        prior[i][call] = 1.0 - perror;

    }

    if ( debug ) {
        for ( size_t i = 0; i < prior.size(); ++i )
            std::cout << "base,qual=" << static_cast<int>(read[i]) << "," << quals[i] << "=== prior[" << i << "]=" << prior[i][0] << "," << prior[i][1] << "," << prior[i][2] << "," << prior[i][3] << std::endl;
        std::cout  << std::endl;
    }
}



/*
 * init_truth - initialize the true, underlying sequence, parameter T in our model.
 */
void init_truth( const _BaseVecVec& reads, const _QualVecVec& quals, _BaseVec& t ) {
    // okay, so in many situations we might initialize with a simple majority consensus
    // but in the case where there may be many non-friends, this will probably just initialize
    // with a mess.  It also depends on what gets used first - the friends or the
    // truth.  It's probably best to start off with the founder read for now.

    t = reads[0];
    return;
}


/*
 * phi_f - distribution around the true sequence T
 */
double phi_f( const _BaseVec& D_j, const _QualVec& Q_j, const _BaseVec& T ) {
    long double result = 0.0;
    double log3 = log(3.0);

    for ( size_t i = 0; i < D_j.size(); ++i ) {
        if ( D_j[i] != ' ' && Q_j[i] >= 0 ) {
            int qual = Q_j[i];
            if ( qual == 0 ) qual = 30;
            result += ( D_j[i] == T[i] ? Precomputed::Q2logp(qual) : (Precomputed::Q2logq(qual) - log3) );
//      result += ( D_j[i] == T[i] ? Precomputed::Q2logp(Q_j[i]) : (Precomputed::Q2logq(Q_j[i]) - log3) );
        }
    }

    long double exp_result = exp(result);

    return exp_result;
}


/*
 * phi_NOTf - alternate model for random bases - 25% for anything
 */
double phi_NOTf( const _BaseVec& D_j, const _QualVec& Q_j, const _BaseVec& T ) {
    double result = 0.0;
    double logp = log(0.25);
    double logq = log(0.75 / 3.0);

    for ( size_t i = 0; i < D_j.size(); ++i ) {
        if ( D_j[i] != ' ' && Q_j[i] >= 0 ) {
            result += ( D_j[i] == T[i] ? logp : logq );
        }
    }

    result = exp(result);

    return result;
}



/*
 * estimate_friends -
 *
 * given a read stack, quality scores, the truth, and a prior on friendship A, estimate the probability
 * that each read in the stack is a true friend of the "truth" sequence.
 */
double estimate_friends( const _BaseVecVec& reads, const _QualVecVec& quals, const _BaseVec& truth, double a, std::vector<double>& pfriend, bool debug = false ) {
    double rmsdiff = 0.0;

    for ( size_t j = 0; j < reads.size(); ++j ) {		// j == 0 is the founder read, we only evaluate for debugging purposes

        double p1 = phi_f(reads[j], quals[j], truth);
        double p2 = phi_NOTf(reads[j], quals[j], truth);

        double new_friend = a*p1 / (a*p1 + (1.0-a)*p2);

        if ( j == 0 && new_friend < 0.5 ) {
            #pragma omp critical
            {
                std::cout << "-------" << std::endl;
                std::cout << "A=" << a << std::endl;
                std::cout << "P1=" << p1 << std::endl;
                std::cout << "P2=" << p2 << std::endl;
                std::cout << "READ0:" << std::endl << reads[0] << std::endl;;
                std::cout << "TRUTH:" << std::endl << truth << std::endl;
                std::cout << "READ0_Q:" << std::endl;
                printq(std::cout, quals[j]) << std::endl;
                double not_friend = (1.0-a)*p2 / ( a*p1 + (1.0-a)*p2 );
                std::cout << "friend=" << new_friend << ", not_friend=" << not_friend << ", sum=" << new_friend+not_friend << std::endl;
            }
        }

        if ( debug )  {
            std::cout << "READ " << j << std::endl;
            std::cout << "TRUTH: " << truth << std::endl;
            std::cout << "READ : " << reads[j] << std::endl;
            std::cout << "QUAL : ";
            printq(std::cout,quals[j]) << std::endl;

            double not_friend = (1.0-a)*p2 / ( a*p1 + (1.0-a)*p2 );
            std::cout << "friend=" << new_friend << ", not_friend=" << not_friend << ", sum=" << new_friend+not_friend << std::endl;
        }

        if ( j > 0 ) {	// we might run to zero for debugging -- check to see if the founder read is itself a friend, but only store values for j > 0
            rmsdiff += ( pfriend[j] - new_friend) * ( pfriend[j] - new_friend );
        }
        pfriend[j] = std::min(std::max(new_friend, 0.0001),1.0-0.0001);

    }


    return rmsdiff;
}


/*
 * estimate_A - update the prior on friendship
 *
 * This turns out to be the mean of our friend estimates.
 *
 */
void estimate_A( std::vector<double>& pfriend, double& a_return ) {
    long double a = 0.0;

    for ( size_t j = 1; j < pfriend.size(); ++j )
        a += pfriend[j];

    a_return= static_cast<double>( a / static_cast<long double>( pfriend.size() ) );
}



/*
 * estimate_truth - given friendships, determine the ML base at each
 * position
 *
 * reads and quals are the usual readstack and associated quality
 * scores.
 *
 * pfriend is a vector of friendship probabilities
 *
 * t_prior is the probabilistic prior - usually based on the founder
 * read and its quality scores
 *
 * truth and truthq are the output true estimated sequence and
 * associated quality (the q is currently unused)
 */
void estimate_truth( const _BaseVecVec& reads, const _QualVecVec& quals, const std::vector<double>& pfriend, const std::vector<BaseProb>& t_prior, _BaseVec& truth, _QualVec& truthq, bool final, bool debug = false, const size_t serialno = 0 ) {
    long double log25 = log(0.25);
    long double log75 = log(0.75);
    long double log3 = log(3.0);
    long double eps = 1e-10;

    std::vector<int> quality_sum( truth.size(), 0);

    for ( size_t i = 0; i < truth.size(); ++i ) {	// for each base position
        vec<long double> p(4,0.);
        size_t count = 0;

        if ( 0 && quals[0][i] == 0 ) {
            truth[i] = reads[0][i];
            truthq[i] = quals[0][i];
            continue;
        }

        for ( size_t j = 1; j < reads.size(); ++j ) {	// for each friend read (past the founder)

            if ( reads[j].size() > i && reads[j][i] != ' ' && quals[j][i] > 0 ) {	// NIW -- how to handle Q1 and Q2...
                count++;
                quality_sum[i] += quals[j][i];
                ForceAssertGe(reads[j][i],0);
                ForceAssertLe(reads[j][i],4);

                int qual = quals[j][i];
                if ( qual == 0) qual = 30;					// this is a hack, but we'd like to include previously corrected bases
                const size_t Dij = static_cast<size_t>(reads[j][i]);
                const double logqij = Precomputed::Q2logq(qual);
                const double logpij = Precomputed::Q2logp(qual);

                double mis, hit;

                mis = (1.0 - pfriend[j]) * (log75 - log3) + pfriend[j] * (logqij - log3);
                hit = (1.0 - pfriend[j]) * log25 + pfriend[j] * logpij;

                if ( debug && ( i == 24 || i == 25)  ) { // NIW TEMPORARY
                    #pragma omp critical
                    std::cout << "pfriend[j]=" << pfriend[j] << ", logqij=" << logqij << ", logpij=" << logpij << ", mis=" << mis << ", hit=" << hit << std::endl;
                }

                p[0] += mis;		// all bases get credit for a miss
                p[1] += mis;
                p[2] += mis;
                p[3] += mis;
                p[Dij] -= mis;		// cheaper to subtract it back out and add in the hit, than to branch here
                p[Dij] += hit;
            }
        }

        if ( count < 3 ) { 		// NIW TEMPORARY - and quality_sum < X ?

            if ( debug && ( i == 24 || i == 25 ) )
                std::cout <<  "TOOK BASE FROM READ due to lack of evidence" << std::endl;

            truth[i] = reads[0][i];
            truthq[i] = quals[0][i];

        } else {

            if ( debug && ( i == 24 || i == 25 ) ) { 	// NIW TEMPORARY
                std::cout << "p[]=" << exp(p[0]) << "," << exp(p[1]) << "," << exp(p[2]) << "," << exp(p[3]) << std::endl;
                std::cout << "prior=" << t_prior[i][0] << "," << t_prior[i][1] << "," << t_prior[i][2] << "," << t_prior[i][3] << std::endl;
            }


            long double gamma=std::max(1.0, 1.0*count);		// gamma is a weight on the prior
            for (size_t idx = 0; idx < 4; ++idx ) {
                p[idx] =  p[idx] + gamma*log(t_prior[i][idx]) ;
            }

            log_domain_norm(p);	// normalize so that probabilities (after exponentiation) sum to one.

            for ( size_t idx = 0; idx < 4; ++idx )
                p[idx] = exp(p[idx]);

            if ( debug && ( i == 24 || i == 25 ) ) { // NIW TEMPORARY
                std::cout << "with prior p[]=" << p[0] << "," << (p[1]) << "," << (p[2]) << "," << (p[3]) << std::endl;
                std::cout << "prior=" << t_prior[i][0] << "," << t_prior[i][1] << "," << t_prior[i][2] << "," << t_prior[i][3] << std::endl;
            }

            long double sum = 0.0;
            for ( size_t idx = 0; idx < 4; ++idx )
                sum += p[idx];

#if 0
            if ( ! ( sum > 0.0 ) ) {
                #pragma omp critical
                {
                    std::cout << "serialno=" << serialno << std::endl;
                    std::cout << "withOUT prior p[])=" << exp(p_tmp[0]) << "," << (exp(p_tmp[1])) << "," << (exp(p_tmp[2])) << "," << (exp(p_tmp[3])) << std::endl;
                    std::cout << "with prior p[]=" << p[0] << "," << (p[1]) << "," << (p[2]) << "," << (p[3]) << std::endl;
                    std::cout << "with prior before exp p[]=" << p_tmp2[0] << "," << (p_tmp2[1]) << "," << (p_tmp2[2]) << "," << (p_tmp2[3]) << std::endl;
                    log_domain_norm( p_tmp2 );
                    std::cout << "with prior after lognorm-exp p[]=" << exp(p_tmp2[0]) << "," << exp(p_tmp2[1]) << "," << exp(p_tmp2[2]) << "," << exp(p_tmp2[3]) << std::endl;
                    std::cout << "prior=" << t_prior[i][0] << "," << t_prior[i][1] << "," << t_prior[i][2] << "," << t_prior[i][3] << std::endl;
                    std::cout << "gamma=" << gamma << std::endl;

                    for ( size_t ji = 0; ji < pfriend.size(); ji++ )
                        std::cout << "pfriend[" << ji << "]=" << pfriend[ji] << std::endl;
                }
            }
#endif
            ForceAssertGt(sum, 0.0);

            for ( size_t idx = 0; idx < 4; ++idx ) {	// normalize and clamp -- this is primarily for Q scores which are going to be fudged later
                p[idx] /= sum;
                p[idx] = std::max<long double>(1e-5, std::min<long double>(1.0-1e-5, p[idx]));
            }




            // index and probability associated with the most probable base

            long double vmax = p[0];
            size_t imax = 0;
            for ( size_t idx = 1; idx < 4; ++idx )
                if ( p[idx] > vmax ) {
                    vmax = p[idx];
                    imax = idx;
                }


#if 0
            truth[i] = imax;
#else
            vec<size_t> index(4,0);
            index[1]=1;
            index[2]=2;
            index[3]=3;
            ReverseSortSync( p, index );
            ForceAssertEq(imax, index[0]);
//#pragma omp critical
//      std::cout << "imax=" << imax << ", index[0]=" << index[0] << ", p[0]=" << p[0] << ", p[1]=" << p[1] << ", p[2]=" << p[2] << ", p[3]=" << p[3] << std::endl;

            if ( final ) {
                if ( p[0] > 0.9 && p[1] < 0.1 ) {
                    truth[i] = index[0];
                } else {
                    truth[i] = reads[0][i];
                }
                truthq[i] = 0;
            } else {
                truth[i] = index[0];
                truthq[i] = static_cast<int>(-10.0*log10(1.0-p[index[0]]));	// this value may be fudged in run_EMEC3
            }
#endif

            //    std::cout << "p[imax]=" << p[imax] << ", truthq[" << i << "]=" << truthq[i] << std::endl;

            if ( debug && ( i == 24 || i == 25 ) )
                std::cout << "WE DECIDED=" << Base::val2Char(truth[i]) << std::endl;
        }
    }
}


/*
 * hack_q - updated quality scores for the new read.
 *
 * For now we use the quality scores from the original read and ZERO if we've updated a position.
 * This is consistent with the previous code in LongProto.
 *
 */
void hack_q( const _BaseVec& call0, const _QualVec& call0q, const _BaseVec& t, _QualVec& q ) {
    ForceAssert( call0.size() == t.size() );
    ForceAssert( call0q.size() == t.size() );
    ForceAssert( q.size() == t.size() );

    for ( size_t i = 0; i < q.size(); ++i ) {
        if ( call0[i] == t[i] )
            q[i] = call0q[i];
        else
            q[i] = 0;
    }
}

void hack_q3( const _BaseVec& call0, const _QualVec& call0q, const _BaseVec& t, _QualVec& q ) {
    ForceAssert( call0.size() == t.size() );
    ForceAssert( call0q.size() == t.size() );
    ForceAssert( q.size() == t.size() );

    for ( size_t i = 0; i < q.size(); ++i ) {
        q[i] = call0q[i];
    }
}

void hack_q2( const _BaseVecVec& call, const _QualVecVec& callq, const std::vector<double>& pfriend, const _BaseVec& t, _QualVec& q) {
    ForceAssert( call[0].size() == t.size() );
    ForceAssert( q.size() == t.size() );
    ForceAssert( call.size() == pfriend.size() );

    for ( size_t i = 0; i < t.size(); ++i ) {
        std::vector<long double> credit(4,0.);
        std::vector<long double> credit2(4,0.);
        if ( t[i] != call[0][i] ) {
            credit[call[0][i]] = callq[0][i];
            credit2[call[0][i]] = 1.0;
            for ( size_t j = 1; j < call.size(); ++j ) {
                unsigned char ci = call[j][i];
                if ( ci != 32 && callq[j][i] > 0 ) {
                    credit[ci] += callq[j][i]*pfriend[j];
                    credit2[ci] += pfriend[j];
                }
            }

            long double sum = 0.;
            for ( size_t j = 0; j < 4; ++j )
                sum+=credit[j];
            for ( size_t j = 0; j < 4; ++j )
                credit[j] /= sum;

            sum = 0.;
            for ( size_t j = 0; j < 4; ++j )
                sum+=credit2[j];
            for ( size_t j = 0; j < 4; ++j )
                credit2[j] /= sum;


            ForceAssertGt( sum, 0. );

            double perror = std::max<long double>( 1e-4, std::min<long double>( 0.75 , (1.0 - credit[t[i]] ) ) );

            if ( perror > 0.74 ) {
                std::cout << "high error: credit=" << credit[0] << "," << credit[1] << "," << credit[2] << "," << credit[3] << std::endl;
                std::cout << "___________ credit=" << credit2[0] << "," << credit2[1] << "," << credit2[2] << "," << credit2[3] << std::endl;
                double alterror = std::max<long double>( 1e-4, std::min<long double>( 0.75 , (1.0 - credit2[t[i]] ) ) );
                std::cout << "alterror: " <<  alterror << std::endl;
            }

            q[i] = static_cast<int>( -10.0 * log10( perror ) );

            ForceAssertLt(q[i], 55);	// somewhat arbitrary, but we're looking for anomalies atm.
        } else {
            q[i] = callq[0][i];
        }
    }
}

} // namespace




/*
 * run_EMEC3 - main entrypoint
 *
 * call  - stack comprising the "founder" read (always position 0) and friends
 * callq - quality scores for reads in call
 * t, q  - output "true" read and updated quality score (zeros for edited positions)
 * pfriend - vector of friendship probabilities
 * trim_to - currently unused
 * debug_level - currently >0 for debugging
 *
 */
void run_EMEC3( const _BaseVecVec& call, const _QualVecVec& callq, _BaseVec& t, _QualVec& q, std::vector<double>& pfriend, int& trim_to, const unsigned debug, const size_t serialno ) {

    trim_to = call[0].size();

    if ( debug > 0 )
        std::cout << "********* START OF RUN ***********" << std::endl;
    size_t nfriends = call.size()-1;

    if ( nfriends < 3 ) {
        // no work to do
        t = call[0];
        q = callq[0];
        return;
    }

    double a = 0.5;			// A - the prior on friendship

    pfriend.resize(nfriends+1, 0.5);	// responsibilities p(F_j = 1 | \theta, \D)

    std::vector<BaseProb> t_prior;
    init_prior( call[0], callq[0], t_prior, false );	// prior on the parameter T - derived from the founder read
    init_truth( call, callq, t );					// initialize t -- currently based on the founder
    q.resize(t.size());


    if ( 0 && serialno == 3152 ) {	// NIW TEMPORARY
        std::cout << "orig0_pos=" << 66 << ": " << Base::val2Char(call[0][66]) << "(" << callq[0][66] << ")" << std::endl;
        std::cout << "prior=" << t_prior[66] << std::endl;
        std::cout << "truth=" << Base::val2Char(t[66]) << std::endl;
    }

    double olddiff = std::numeric_limits<double>::max();

    const unsigned int maxiter = 20;
    unsigned int iter = 0;
    for ( ; iter < maxiter; ++iter ) {
        double thisdiff;



        // E-step
        thisdiff = estimate_friends( call, callq, t, a, pfriend, debug > 0);

        if ( debug > 0 ) {
            for ( size_t i = 0; i < pfriend.size(); ++i )
                std::cout << "pfriend[" << i << "]=" << pfriend[i] << std::endl;
        }

        double fcount = 0.0;
        for ( size_t i = 1; i < pfriend.size(); ++i )
            if ( pfriend[i] > 0.5 )
                fcount+=1.0;
        if ( fcount < 3.0 || fcount / pfriend.size() < 0.05 )  {
            if ( debug > 0 )
                std::cout << "short circuting due to lack of friends" << std::endl;
            t = call[0];
            q = callq[0];

            return;
        }




        estimate_A( pfriend, a );

        if (debug > 0)
            std::cout << "A=" << a << std::endl;

        estimate_truth( call, callq, pfriend, t_prior, t, q, false, debug > 0, serialno);

#if 0
        if ( 0 && serialno == 3152 )  {			// NIW TEMPORARY
            std::vector<long long> evidence(4,0);

            for ( size_t j = 1; j < call.size(); ++j ) {
                if ( call[j][66] == ' ' ) continue;
                std::cout << ( ( pfriend[j-1]>0.5) ? "FRIEND" : "friend" );
                std::cout << " " << j << "(" << pfriend[j-1] << ")" << " " << ( call[j][66] == 32 ? ' ' : Base::val2Char(call[j][66]) ) << "(" << callq[j][66] << ")" <<std::endl;
                if ( call[j][66] != 32 )
                    evidence[call[j][66]]+=callq[j][66];
            }

            std::cout << "evidence[A]=" << evidence[0] << ", [C]=" << evidence[1] << ", [G]=" << evidence[2] << ", [T]=" << evidence[3] << std::endl;
        }
#endif

//    hack_q(call[0], callq[0], t, q);

        if (debug > 0)
            std::cout << "T=" << t <<std::endl;

        if ( debug > 0)
            std::cout << "olddiff=" << olddiff << ", thisdiff=" << thisdiff << std::endl;

        if ( std::abs(olddiff - thisdiff) < 1e-3 )
            break;

        olddiff = thisdiff;
    }

    if ( pfriend[0] < 0.5 ) {
        #pragma omp critical
        {
            std::cout << Date() << ": EM converged to wrong mode: id=" << serialno << ", pfriend[0]=" << pfriend[0] << std::endl;
        }
        t = call[0];
        q = callq[0];

        return;
    }

    if ( iter >= maxiter ) {
        #pragma omp critical
        {
            std::cout << Date() << ": EM didn't converge: id=" << serialno << ", pfriend[0]=" << pfriend[0] << ", olddiff=" << olddiff << std::endl;
        }
#if 0
        t = call[0];
        q = callq[0];
        return;
#endif
        for ( size_t i = 1; i < pfriend.size(); ++i )
            pfriend[i] = 1.;
    }


    //  various strategies for recomputing Q.  For now we just zero out
    //  bases we change below after the consensus.
//  hack_q(call[0], callq[0], t, q);
//  hack_q3(call[0], callq[0], t, q);
//  hack_q2(call, callq, pfriend, t, q);


    // check how many friends are available...
    nfriends = 0;
    for ( size_t i = 1; i < pfriend.size(); ++i )
        if ( pfriend[i] > 0.5 )
            nfriends++;

    if ( nfriends < 3 || nfriends < 0.05 * pfriend.size() ) {
        // estimates are unreliable with too little data
        t = call[0];
        q = callq[0];
    }


    /*
     * We've clustered code into friends and non-friends and T should be
     * the maximum-likelihood estimate of the "center" of the friend
     * cluster.  That's not necessarily very conservative w.r.t.
     */
    else {
#if 0
        estimate_truth( call, callq, pfriend, t_prior, t, q, true, debug > 0, serialno);
#else

        for ( size_t i = 0; i < t.size(); ++i ) {
            t[i] = call[0][i];
            q[i] = callq[0][i];

            BaseProb ptruth(4, 0.);
            long double sum = 0.0;
            size_t friend_count = 0;
            for ( size_t j = 1; j < call.size(); ++j ) {
                size_t c = call[j][i];
                int q = callq[j][i] == 0 ? 30 : callq[j][i];
                if ( c != 32 && pfriend[j] > 0.5 ) {
                    ptruth[c] += pfriend[j];
                    sum += pfriend[j];
                    friend_count++;
                }
            }

            vec<size_t> ids(ptruth.size());
            for ( size_t k = 0; k < ptruth.size(); ++k ) {
                ids[k] = k;
            }

            ReverseSortSync( ptruth, ids );

            if ( sum > 0 && friend_count > 3 ) {		// second case obviates the first, but we'll leave it for the moment
                for ( size_t k = 0; k < ptruth.size(); ++k ) {
                    ptruth[k] /= sum;
                }
                if ( ptruth[0] > 0.70 && ptruth[1] < 0.10 ) {

                    t[i] = ids[0];
                    if ( t[i] != call[0][i] )
                        q[i] = 0;
                }
            }

        }

#endif
    }



}

