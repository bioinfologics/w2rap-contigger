#ifndef MATCHRESULT_INCLUDED
#define MATCHRESULT_INCLUDED 1
typedef struct {
  unsigned long int query_id;
  unsigned long int target_id;
  unsigned long int query_size;
  unsigned long int qstart;
  unsigned long int tstart;
  unsigned long int len;
  bool reverse;
  double prob;
  double ident;
} t_result;
#endif //MATCHRESULT_INCLUDED
