// Fixed-storage, integer-only hot-path observation; one ledger per native run.
#ifndef OBSERVE85_H
#define OBSERVE85_H
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include <cstdio>
namespace Observe85 {
enum Source { SETUP, CORE, ACCEPT, SIB_MID, FRESH, NSOURCE };
inline unsigned long long entered[NSOURCE]={}, completed[NSOURCE]={};
inline unsigned long long forwards[NSOURCE]={}, reverses[NSOURCE]={}, builds[NSOURCE]={};
inline unsigned long long virial_requested[NSOURCE]={};
inline void reset() {
  for(int i=0;i<NSOURCE;i++) entered[i]=completed[i]=forwards[i]=reverses[i]=builds[i]=virial_requested[i]=0;
}
inline void pair(LAMMPS_NS::Pair *p,int e,int v,Source s) {
  ++entered[s]; if(v & 3) ++virial_requested[s]; p->compute(e,v); ++completed[s];
}
inline void forward(LAMMPS_NS::Comm *c,Source s) {++forwards[s]; c->forward_comm();}
inline void reverse(LAMMPS_NS::Comm *c,Source s) {++reverses[s]; c->reverse_comm();}
inline void build(LAMMPS_NS::Neighbor *n,int flag,Source s) {++builds[s]; n->build(flag);}
inline void emit(FILE *f) {
  if(!f) return;
  const char *names[NSOURCE]={"setup","core","accept","sib_mid","fresh"};
  std::fprintf(f,"OBSERVE85 {");
  for(int i=0;i<NSOURCE;i++) std::fprintf(f,
    "%s\"%s\":{\"entered\":%llu,\"completed\":%llu,\"forward\":%llu,\"reverse\":%llu,\"build\":%llu,\"virial_requested\":%llu}",
    i?",":"",names[i],entered[i],completed[i],forwards[i],reverses[i],builds[i],virial_requested[i]);
  std::fprintf(f,"}\n");
}
}
#endif
