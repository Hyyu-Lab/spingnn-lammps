// LP-SIB: lagged-predictor variant of the repository spin/sib (normalized
// predicted midpoint).  One pair evaluation per step, at the midpoint.
#ifdef FIX_CLASS
FixStyle(spin/lpsib,FixSpinLPSIB);
#else
#ifndef LMP_FIX_SPIN_LPSIB_H
#define LMP_FIX_SPIN_LPSIB_H
#include "fix.h"
namespace LAMMPS_NS {
class FixSpinLPSIB : public Fix {
 public:
  FixSpinLPSIB(class LAMMPS *,int,char **);
  ~FixSpinLPSIB() override;
  int setmask() override;
  void *extract(const char *,int &) override;
  void init() override;
  int modify_param(int,char **) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;
  void grow_arrays(int) override;
  void copy_arrays(int,int,int) override;
  void set_arrays(int) override;
  double memory_usage() override;
  double compute_vector(int) override;
 private:
  double temp,alpha,hbar,gil_factor,D,C_spin,dts;
  int seed,accept,cache_valid;
  bigint nsteps,naccept;
  double maxnormdev;
  class RanMars *random;
  double **store;             // [0..2] cached midpoint fm, [3..5] s_n, [6..8] R
  int npairspin;
  class PairSpin **spin_pairs;
  int lattice_moving,predictor_fresh,ou_enabled,ou_seed;
  double ou_temp,ou_damp;
  class RanMars *ou_random;
  int nprecession;
  class FixPrecessionSpin **precessions;
  int box_iso,box_initialized;
  int box_ou_enabled,box_ou_seed,npt_test_negative;
  double box_ou_temp,box_ou_damp,box_noise[2];
  class RanMars *box_random;
  double box_V,box_pV,box_W,box_Pbar,box_P,box_shape[3],box_origin[3],box_virial[6];
  void box_init();
  void box_T(double);
  void box_D(double);
  void box_sync();
  void box_check_pair(class Pair *,const char *);
  void box_accept_comm();
  void box_bath_half(int);
  void ou_half(int);
  void field(int);            // replace fm[i] with the pair-spin field at current sp
  void cayley(const double *s,const double *F,const double *R,double dtspin,const double *x,double *out);
};
}
#endif
#endif
