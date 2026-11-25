// Step85: unified frozen/moving LP-SIB. Frozen spin arithmetic retained.
//
// Identical to the repository spin/sib (fix_spin_sib.cpp, the normalized
// predicted-midpoint variant) except for the predictor field: spin/sib
// evaluates the field at s_n, LP-SIB reuses the field evaluated at the previous
// step's normalized predicted midpoint.  The single pair evaluation per step is
// the ordinary Verlet force call, which therefore happens at the midpoint:
//   initial_integrate: predictor with cached field, sp <- unit((s_n+s~)/2)
//   [Verlet: forward comm, force clear, pair compute at the midpoint]
//   final_integrate:   field at midpoint -> cache; corrector from s_n -> s_{n+1}
//   end_of_step (every `accept` steps): one extra pair evaluation at the
//     accepted state so that dumps/thermo/compute spin see fm and pe of s_{n+1}.
// The extra evaluation does not touch the integrator state, so the trajectory
// does not depend on `accept`.  The first step of every run uses the field at
// s_n from the Verlet setup force call (the cache is invalidated in init()).
// Noise, Cayley matrices and units follow spin/sib exactly.
#include "fix_spin_lpsib.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "pair_spin.h"
#include "fix_precession_spin.h"
#include "domain.h"
#include "neighbor.h"
#include "pair_hybrid.h"
#include "observe85.h"
#include <vector>
#include <algorithm>
#include <cstdio>
#include "random_mars.h"
#include "update.h"
#include <cmath>
#include <cstring>
using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
void vec2mat(const double v[3],double m[3][3]){
  m[0][0]=0.0;  m[0][1]=v[2]; m[0][2]=-v[1];
  m[1][0]=-v[2];m[1][1]=0.0;  m[1][2]=v[0];
  m[2][0]=v[1]; m[2][1]=-v[0];m[2][2]=0.0;
}
}

FixSpinLPSIB::FixSpinLPSIB(LAMMPS *lmp,int narg,char **arg):Fix(lmp,narg,arg),
  random(nullptr),store(nullptr),npairspin(0),spin_pairs(nullptr),
  lattice_moving(0),predictor_fresh(0),ou_enabled(0),ou_random(nullptr),nprecession(0),precessions(nullptr),
  box_iso(0),box_initialized(0),box_ou_enabled(0),box_ou_seed(0),npt_test_negative(0),
  box_ou_temp(0),box_ou_damp(0),box_noise{},box_random(nullptr),
  box_V(0),box_pV(0),box_W(0),box_Pbar(0),box_P(0),box_virial{} {
  // fix ID GROUP spin/lpsib T alpha seed lattice frozen accept N
  if(narg<10) error->all(FLERR,"Illegal fix spin/lpsib command");
  temp=utils::numeric(FLERR,arg[3],false,lmp);
  alpha=utils::numeric(FLERR,arg[4],false,lmp);
  seed=utils::inumeric(FLERR,arg[5],false,lmp);
  if(std::strcmp(arg[6],"lattice")) error->all(FLERR,"spin/lpsib expected lattice frozen|moving");
  if(!std::strcmp(arg[7],"moving")) lattice_moving=1;
  else if(std::strcmp(arg[7],"frozen")) error->all(FLERR,"spin/lpsib expected lattice frozen|moving");
  if(std::strcmp(arg[8],"accept")) error->all(FLERR,"Expected: accept N");
  accept=utils::inumeric(FLERR,arg[9],false,lmp);
  for(int j=10;j<narg;){
    if(!std::strcmp(arg[j],"predictor") && j+1<narg){
      if(!std::strcmp(arg[j+1],"lagged")) predictor_fresh=0;
      else if(!std::strcmp(arg[j+1],"fresh")) predictor_fresh=1;
      else error->all(FLERR,"spin/lpsib predictor must be lagged|fresh");
      j+=2;
    } else if(!std::strcmp(arg[j],"ou") && j+3<narg){
      if(ou_enabled) error->all(FLERR,"spin/lpsib duplicate OU bath");
      ou_temp=utils::numeric(FLERR,arg[j+1],false,lmp);
      ou_damp=utils::numeric(FLERR,arg[j+2],false,lmp);
      ou_seed=utils::inumeric(FLERR,arg[j+3],false,lmp);
      if(!lattice_moving || !std::isfinite(ou_temp+ou_damp) || ou_temp<0 || ou_damp<=0 ||
         ou_seed<=0 || ou_seed>900000000 || ou_seed==seed)
        error->all(FLERR,"spin/lpsib invalid independent OU bath");
      ou_enabled=1; j+=4;
    } else if(!std::strcmp(arg[j],"box_ou") && j+3<narg){
      if(box_ou_enabled) error->all(FLERR,"S85B_DUPLICATE_BOX_OU");
      box_ou_temp=utils::numeric(FLERR,arg[j+1],false,lmp);
      box_ou_damp=utils::numeric(FLERR,arg[j+2],false,lmp);
      box_ou_seed=utils::inumeric(FLERR,arg[j+3],false,lmp);
      if(!std::isfinite(box_ou_temp+box_ou_damp) || box_ou_temp<0 || box_ou_damp<=0 ||
         box_ou_seed<=0 || box_ou_seed>900000000)
        error->all(FLERR,"S85B_INVALID_BOX_OU");
      box_ou_enabled=1; j+=4;
    } else if(!std::strcmp(arg[j],"box") && j+3<narg){
      if(box_iso || std::strcmp(arg[j+1],"iso")) error->all(FLERR,"S85N_BOX_SYNTAX");
      box_Pbar=utils::numeric(FLERR,arg[j+2],false,lmp);
      box_W=utils::numeric(FLERR,arg[j+3],false,lmp);
      if(!std::isfinite(box_Pbar) || !std::isfinite(box_W) || box_W<=0)
        error->all(FLERR,"S85N_BOX_PARAMETERS");
      box_iso=1; j+=4;
    } else error->all(FLERR,"spin/lpsib unknown optional keyword");
  }
  if(box_iso && !lattice_moving) error->all(FLERR,"S85N_BOX_REQUIRES_MOVING");
  // Step91: fresh uses a second midpoint force/virial call for canonical-V NPT.
  if(box_ou_enabled && (!box_iso || box_ou_seed==seed || (ou_enabled && box_ou_seed==ou_seed)))
    error->all(FLERR,"S85B_BOX_OU_REQUIRES_INDEPENDENT_STREAM");
  if(box_iso) box_change=BOX_CHANGE_SIZE; // reject conflicting fixes explicitly in init
  if(temp<0 || alpha<0 || seed<=0 || accept<0 || !std::isfinite(temp+alpha))
    error->all(FLERR,"Invalid spin/lpsib parameters");
  if(!atom->sp_flag || std::strcmp(update->unit_style,"metal"))
    error->all(FLERR,"spin/lpsib requires atom/spin style and metal units");
  if(comm->nprocs!=1) error->all(FLERR,"spin/lpsib currently supports exactly one MPI rank");
  random=new RanMars(lmp,seed+comm->me);
  if(ou_enabled) ou_random=new RanMars(lmp,ou_seed);
  if(box_ou_enabled) box_random=new RanMars(lmp,box_ou_seed);
  time_integrate=1;
  nevery=accept>0?accept:1;
  vector_flag=1; size_vector=box_iso?16:3; global_freq=1; extvector=0;
  peratom_flag=1; size_peratom_cols=box_iso?21:(lattice_moving?15:9); peratom_freq=1;
  grow_arrays(atom->nmax); atom->add_callback(0);
  for(int i=0;i<atom->nmax;i++) set_arrays(i);
  cache_valid=0; nsteps=naccept=0; maxnormdev=0.0;
}

FixSpinLPSIB::~FixSpinLPSIB(){
  atom->delete_callback(id,0);
  memory->destroy(store);
  delete random;
  delete ou_random;
  delete box_random;
  delete [] precessions;
  delete [] spin_pairs;
}

int FixSpinLPSIB::setmask(){
  int mask=INITIAL_INTEGRATE|FINAL_INTEGRATE;
  if(accept>0) mask|=END_OF_STEP;
  return mask;
}

void FixSpinLPSIB::grow_arrays(int n){memory->grow(store,n,box_iso?21:15,"spin/lpsib:store"); array_atom=store;}
void FixSpinLPSIB::copy_arrays(int i,int j,int){for(int k=0;k<(box_iso?21:15);k++) store[j][k]=store[i][k];}
void FixSpinLPSIB::set_arrays(int i){for(int k=0;k<(box_iso?21:15);k++) store[i][k]=0.0;}
double FixSpinLPSIB::memory_usage(){return (double)atom->nmax*(box_iso?21:15)*sizeof(double);}

void FixSpinLPSIB::init(){
  if(std::strcmp(update->integrate_style,"verlet")) error->all(FLERR,"spin/lpsib requires run_style verlet");
  if(box_iso){
    if(domain->triclinic || domain->dimension!=3 || !domain->xperiodic || !domain->yperiodic || !domain->zperiodic)
      error->all(FLERR,"S85N_BOX_REQUIRES_ORTHOGONAL_PERIODIC_3D");
    if(neighbor->every!=1 || neighbor->delay!=0 || !neighbor->dist_check || neighbor->build_once)
      error->all(FLERR,"S85N_BOX_NEIGHBOR_POLICY");
    box_check_pair(force->pair,force->pair_style);
  }
  if(lattice_moving){
    if(igroup!=0) error->all(FLERR,"spin/lpsib moving requires group all");
    if(!force->newton_pair) error->all(FLERR,"spin/lpsib moving requires newton on");
    if(atom->molecular!=Atom::ATOMIC || force->kspace || domain->triclinic || atom->bec_flag)
      error->all(FLERR,"spin/lpsib moving unsupported molecular/kspace/triclinic/BEC configuration");
  }
  std::vector<FixPrecessionSpin*> prec;
  for(int j=0;j<modify->nfix;j++){
    Fix *other=modify->fix[j]; if(other==this) continue;
    if(box_iso && (other->box_change || utils::strmatch(other->style,"^npt") ||
        utils::strmatch(other->style,"^nph") || utils::strmatch(other->style,"^press/berendsen") ||
        utils::strmatch(other->style,"^box/relax") || utils::strmatch(other->style,"^deform")))
      error->all(FLERR,"S85N_BOX_CONFLICTING_FIX");
    const bool overlap=(other->groupbit & groupbit)!=0 || lattice_moving;
    if(utils::strmatch(other->style,"^precession/spin")){
      auto *q=dynamic_cast<FixPrecessionSpin*>(other);
      if(!q || q->stt_flag) error->all(FLERR,"spin/lpsib requires conservative precession: STT forbidden");
      prec.push_back(q); continue;
    }
    if(utils::strmatch(other->style,"^setforce/spin") || utils::strmatch(other->style,"^langevin/spin"))
      error->all(FLERR,"spin/lpsib does not support setforce/spin or langevin/spin");
    if(lattice_moving && overlap && (utils::strmatch(other->style,"^langevin") ||
        utils::strmatch(other->style,"^temp/") || utils::strmatch(other->style,"^gle") ||
        utils::strmatch(other->style,"^gld") || utils::strmatch(other->style,"^viscous") ||
        utils::strmatch(other->style,"^brownian")))
      error->all(FLERR,"spin/lpsib moving forbids external lattice thermostat; use in-fix OU");
    if(other->time_integrate) error->all(FLERR,"spin/lpsib must be the only time integrator");
    if(lattice_moving && overlap && (modify->fmask[j] & (PRE_FORCE|POST_FORCE|POST_INTEGRATE)))
      error->all(FLERR,"spin/lpsib moving unsupported external force/integration modifier");
  }
  delete [] precessions; nprecession=(int)prec.size();
  precessions=nprecession?new FixPrecessionSpin*[nprecession]:nullptr;
  for(int j=0;j<nprecession;j++) precessions[j]=prec[j];
  delete [] spin_pairs; spin_pairs=nullptr; npairspin=0;
  std::vector<PairSpin*> pairs;
  for(int n=1;;n++){
    Pair *raw=force->pair_match("^spin",0,n);
    if(!raw) break;
    auto *q=dynamic_cast<PairSpin*>(raw);
    if(!q) error->all(FLERR,"Spin pair style is not a PairSpin");
    if(std::find(pairs.begin(),pairs.end(),q)!=pairs.end()) break;
    pairs.push_back(q);
  }
  npairspin=(int)pairs.size();
  spin_pairs=npairspin?new PairSpin*[npairspin]:nullptr;
  for(int j=0;j<npairspin;j++) spin_pairs[j]=pairs[j];
  if(comm->me==0 && screen) std::fprintf(screen,
    "S85_NATIVE_CONSTANTS {\"ftm2v\":%.17g,\"hplanck\":%.17g,\"MY_2PI\":%.17g,\"hbar\":%.17g,\"boltz\":%.17g,\"mvv2e\":%.17g}\n",
    force->ftm2v,force->hplanck,MathConst::MY_2PI,force->hplanck/MathConst::MY_2PI,force->boltz,force->mvv2e);
  gil_factor=1.0/(1.0+alpha*alpha);
  hbar=force->hplanck/MathConst::MY_2PI;
  D=alpha*gil_factor*force->boltz*temp;
  C_spin=std::sqrt(2.0*D);
  dts=2.0*update->dt*gil_factor/hbar;
  cache_valid=0;
  if(box_iso) box_init();
}

// Same as FixNVESpinSIB::ComputeInteractionsSpin with only pair/spin terms.
void FixSpinLPSIB::field(int i){
  double fmi[3]={0.0,0.0,0.0};
  for(int k=0;k<npairspin;k++) spin_pairs[k]->compute_single_pair(i,fmi);
  for(int k=0;k<nprecession;k++) precessions[k]->compute_single_precession(i,atom->sp[i],fmi);
  for(int k=0;k<3;k++) atom->fm[i][k]=fmi[k];
}

// out = (I - dtspin/2 A - sqrt(dtspin)/2 B)^-1 (I + dtspin/2 A + sqrt(dtspin)/2 B) x,
// A = mat(-F - alpha s x F), B = mat(C(-R - alpha s x R)), as in spin/sib.
void FixSpinLPSIB::cayley(const double *s,const double *F,const double *R,double dtspin,const double *x,double *out){
  double t[3],a[3],b[3],Am[3][3],Bm[3][3],P[3][3],Q[3][3],Pi[3][3],y[3];
  MathExtra::cross3(s,F,t); for(int k=0;k<3;k++) a[k]=-F[k]-alpha*t[k];
  MathExtra::cross3(s,R,t); for(int k=0;k<3;k++) b[k]=C_spin*(-R[k]-alpha*t[k]);
  vec2mat(a,Am); vec2mat(b,Bm);
  const double sq=std::sqrt(dtspin);
  for(int r=0;r<3;r++) for(int c=0;c<3;c++){
    const double I=(r==c)?1.0:0.0;
    P[r][c]=I-0.5*dtspin*Am[r][c]-0.5*sq*Bm[r][c];
    Q[r][c]=I+0.5*dtspin*Am[r][c]+0.5*sq*Bm[r][c];
  }
  const double det=P[0][0]*P[1][1]*P[2][2]+P[0][1]*P[1][2]*P[2][0]+P[0][2]*P[1][0]*P[2][1]
    -P[0][0]*P[1][2]*P[2][1]-P[0][1]*P[1][0]*P[2][2]-P[0][2]*P[1][1]*P[2][0];
  if(det==0.0) error->one(FLERR,"spin/lpsib singular Cayley matrix");
  const double di=1.0/det;
  Pi[0][0]= di*(P[1][1]*P[2][2]-P[1][2]*P[2][1]);
  Pi[0][1]=-di*(P[0][1]*P[2][2]-P[0][2]*P[2][1]);
  Pi[0][2]= di*(P[0][1]*P[1][2]-P[0][2]*P[1][1]);
  Pi[1][0]=-di*(P[1][0]*P[2][2]-P[1][2]*P[2][0]);
  Pi[1][1]= di*(P[0][0]*P[2][2]-P[0][2]*P[2][0]);
  Pi[1][2]=-di*(P[0][0]*P[1][2]-P[0][2]*P[1][0]);
  Pi[2][0]= di*(P[1][0]*P[2][1]-P[1][1]*P[2][0]);
  Pi[2][1]=-di*(P[0][0]*P[2][1]-P[0][1]*P[2][0]);
  Pi[2][2]= di*(P[0][0]*P[1][1]-P[0][1]*P[1][0]);
  MathExtra::matvec(Q,x,y);
  MathExtra::matvec(Pi,y,out);
}

void FixSpinLPSIB::initial_integrate(int){
  double **sp=atom->sp, **fm=atom->fm; int *mask=atom->mask; const int nlocal=atom->nlocal;
  if(!cache_valid){
    // Field at s_n: Verlet setup has just evaluated the pair style at s_n.
    for(int i=0;i<nlocal;i++) if(mask[i]&groupbit) field(i);
    for(int i=0;i<nlocal;i++) if((mask[i]&groupbit)&&sp[i][3]>0) for(int k=0;k<3;k++) store[i][k]=fm[i][k];
    cache_valid=1;
  }
  if(box_iso){
    if(ou_enabled || box_ou_enabled) box_bath_half(0);
    box_D(0.5*update->dt);
  }
  else if(lattice_moving){
    if(ou_enabled) ou_half(9);
    for(int i=0;i<atom->nlocal;i++) if(atom->mask[i]&groupbit)
      for(int k=0;k<3;k++) atom->x[i][k]+=0.5*update->dt*atom->v[i][k];
  }
  for(int i=0;i<nlocal;i++){
    if(!((mask[i]&groupbit)&&sp[i][3]>0)) continue;
    double F[3],pred[3];
    for(int k=0;k<3;k++){
      F[k]=store[i][k]*hbar;
      store[i][3+k]=sp[i][k];
      store[i][6+k]=random->gaussian();
    }
    if(predictor_fresh) continue;
    const double dtspin=dts/sp[i][3];
    cayley(&store[i][3],F,&store[i][6],dtspin,&store[i][3],pred);
    for(int k=0;k<3;k++) sp[i][k]=0.5*(store[i][3+k]+pred[k]);
    if(lattice_moving){
      double norm=std::sqrt(MathExtra::dot3(sp[i],sp[i]));
      if(!std::isfinite(norm) || 2.0*norm<1e-6) error->one(FLERR,"S85_MIDPOINT_DENOMINATOR");
    }
    MathExtra::norm3(sp[i]);
  }
}

void FixSpinLPSIB::final_integrate(){
  double **sp=atom->sp, **fm=atom->fm; int *mask=atom->mask; const int nlocal=atom->nlocal;
  if(predictor_fresh){
    // Ordinary Verlet evaluated (xhalf,sold). Finish fresh prediction now.
    // All fields are gathered before changing any spin.
    for(int i=0;i<nlocal;i++) if(mask[i]&groupbit) field(i);
    for(int i=0;i<nlocal;i++){
      if(!((mask[i]&groupbit)&&sp[i][3]>0)) continue;
      double F[3],pred[3];
      for(int k=0;k<3;k++) F[k]=fm[i][k]*hbar;
      const double dtspin=dts/sp[i][3];
      cayley(&store[i][3],F,&store[i][6],dtspin,&store[i][3],pred);
      for(int k=0;k<3;k++) sp[i][k]=0.5*(store[i][3+k]+pred[k]);
      const double norm=std::sqrt(MathExtra::dot3(sp[i],sp[i]));
      if(!std::isfinite(norm) || 2.0*norm<1e-6) error->one(FLERR,"S85_MIDPOINT_DENOMINATOR");
      MathExtra::norm3(sp[i]);
    }
    Observe85::forward(comm,Observe85::FRESH);
    size_t nbytes=sizeof(double)*atom->nlocal;
    if(force->newton) nbytes+=sizeof(double)*atom->nghost;
    if(nbytes){
      std::memset(&atom->f[0][0],0,3*nbytes);
      atom->avec->force_clear(0,nbytes);
    }
    // The fresh predictor used the ordinary Verlet call; this second call
    // is at the spin midpoint. NPT must request and use its virial, never
    // the predictor-state virial. No-box arithmetic/flags remain unchanged.
    Observe85::pair(force->pair,1,box_iso?1:0,Observe85::FRESH);
    if(force->newton) Observe85::reverse(comm,Observe85::FRESH);
    for(int j=0;j<nprecession;j++) precessions[j]->post_force(0);
  }
  // Pair/field now at (xhalf,mid), either ordinary call(A) or fresh second call(C).
  for(int i=0;i<nlocal;i++) if(mask[i]&groupbit) field(i);
  for(int i=0;i<nlocal;i++){
    if(!((mask[i]&groupbit)&&sp[i][3]>0)) continue;
    double F[3],mid[3],next[3];
    for(int k=0;k<3;k++){ store[i][k]=fm[i][k]; F[k]=fm[i][k]*hbar; mid[k]=sp[i][k]; }
    const double dtspin=dts/sp[i][3];
    cayley(mid,F,&store[i][6],dtspin,&store[i][3],next);
    for(int k=0;k<3;k++) sp[i][k]=next[k];
    MathExtra::norm3(sp[i]);
    const double n=std::sqrt(sp[i][0]*sp[i][0]+sp[i][1]*sp[i][1]+sp[i][2]*sp[i][2]);
    if(std::fabs(n-1.0)>maxnormdev) maxnormdev=std::fabs(n-1.0);
  }
  if(box_iso){
    if(!force->pair->vflag_global) error->all(FLERR,"S85N_MISSING_CORE_VIRIAL");
    const double L=std::pow(box_V,1.0/3.0);
    for(int k=0;k<6;k++) box_virial[k]=force->pair->virial[k];
    for(int i=0;i<nlocal;i++) for(int k=0;k<3;k++)
      store[i][18+k]+=update->dt*L*atom->f[i][k];
    box_pV+=update->dt*((box_virial[0]+box_virial[1]+box_virial[2])/(3.0*box_V)-box_P);
    if(npt_test_negative==2) box_pV-=update->dt*force->boltz*temp/box_V;
    box_D(0.5*update->dt);
    if(ou_enabled || box_ou_enabled) box_bath_half(1);
  } else if(lattice_moving){
    for(int i=0;i<nlocal;i++) if(mask[i]&groupbit){
      const double mass=atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]];
      const double kick=update->dt*force->ftm2v/mass;
      for(int k=0;k<3;k++){
        atom->v[i][k]+=kick*atom->f[i][k];
        atom->x[i][k]+=0.5*update->dt*atom->v[i][k];
      }
    }
    if(ou_enabled) ou_half(12);
  }
  nsteps++;
}

// Accepted-state pair evaluation for output only; does not alter the cache.
void FixSpinLPSIB::end_of_step(){
  if(box_iso) box_accept_comm();
  else Observe85::forward(comm,Observe85::ACCEPT);
  size_t nbytes=sizeof(double)*atom->nlocal;
  if(force->newton) nbytes+=sizeof(double)*atom->nghost;
  if(nbytes){
    std::memset(&atom->f[0][0],0,3*nbytes);
    atom->avec->force_clear(0,nbytes);
  }
  Observe85::pair(force->pair,1,box_iso?1:0,Observe85::ACCEPT);
  if(force->newton) Observe85::reverse(comm,Observe85::ACCEPT);
  for(int j=0;j<nprecession;j++) precessions[j]->post_force(0);
  naccept++;
}

double FixSpinLPSIB::compute_vector(int n){
  if(n==0) return (double)nsteps;
  if(n==1) return (double)naccept;
  if(n==2) return maxnormdev;
  if(n==3) return box_V;
  if(n==4) return box_pV;
  if(n==5) return box_pV/box_W/box_V;
  if(n==6) return box_P;
  if(n==7) return force->nktv2p;
  if(n==14 || n==15) return box_noise[n-14];
  return box_virial[n-8];
}

void *FixSpinLPSIB::extract(const char *name,int &dim){
  dim=0;
  if(!std::strcmp(name,"lpsib_box_iso")) return &box_iso;
  if(!std::strcmp(name,"spin_lattice_moving")) return &lattice_moving;
  return nullptr;
}

void FixSpinLPSIB::ou_half(int offset){
  const double decay=std::exp(-0.5*update->dt/ou_damp);
  const double variance=-std::expm1(-update->dt/ou_damp)*force->boltz*ou_temp*force->ftm2v;
  for(int i=0;i<atom->nlocal;i++) if(atom->mask[i]&groupbit){
    const double mass=atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]];
    const double sigma=std::sqrt(variance/mass);
    for(int k=0;k<3;k++){
      store[i][offset+k]=ou_random->gaussian();
      atom->v[i][k]=decay*atom->v[i][k]+sigma*store[i][offset+k];
    }
  }
}


int FixSpinLPSIB::modify_param(int narg,char **arg){
  if(!std::strcmp(arg[0],"npt_test_negative")){
    if(narg<2 || !box_iso || box_initialized || !ou_enabled || !box_ou_enabled || alpha<=0 ||
       temp<=0 || temp!=ou_temp || temp!=box_ou_temp)
      error->all(FLERR,"S85B_TEST_NEGATIVE_REQUIRES_THREE_EQUAL_BATHS");
    if(!std::strcmp(arg[1],"no_box_noise")) npt_test_negative=1;
    else if(!std::strcmp(arg[1],"wrong_measure")) npt_test_negative=2;
    else error->all(FLERR,"S85B_UNKNOWN_TEST_NEGATIVE");
    if(comm->me==0 && screen) std::fprintf(screen,"S85B_TEST_NEGATIVE %s\n",arg[1]);
    return 2;
  }
  if(!std::strcmp(arg[0],"box_momentum")){
    if(narg<2 || !box_iso || box_initialized) error->all(FLERR,"S85N_BOX_INITIAL_MOMENTUM_ONLY");
    box_pV=utils::numeric(FLERR,arg[1],false,lmp);
    if(!std::isfinite(box_pV)) error->all(FLERR,"S85N_BOX_INITIAL_MOMENTUM_ONLY");
    return 2;
  }
  return 0;
}

void FixSpinLPSIB::box_check_pair(Pair *p,const char *name){
  if(!p || !name || !p->compute_flag || p->tail_flag) error->all(FLERR,"S85N_BOX_PAIR_VIRIAL_UNVALIDATED");
  if(!std::strcmp(name,"hybrid") || !std::strcmp(name,"hybrid/overlay")){
    auto *hybrid=dynamic_cast<PairHybrid*>(p);
    if(!hybrid) error->all(FLERR,"S85N_BOX_PAIR_VIRIAL_UNVALIDATED");
    for(int j=0;j<hybrid->nstyles;j++) box_check_pair(hybrid->styles[j],hybrid->keywords[j]);
    return;
  }
  if(std::strcmp(name,"lj/cut") && std::strcmp(name,"spin/exchange") && std::strcmp(name,"zero"))
    error->all(FLERR,"S85N_BOX_PAIR_VIRIAL_UNVALIDATED");
}

// Exact conditional OU in canonical momenta.  It adds no force call or
// communication.  The legacy no-box OU and the all-baths-off arithmetic stay
// untouched.  The last half-step only needs to synchronize physical velocity.
void FixSpinLPSIB::box_bath_half(int half){
  const double L=std::pow(box_V,1.0/3.0);
  if(ou_enabled){
    const double decay=std::exp(-0.5*update->dt/ou_damp);
    const double variance=-std::expm1(-update->dt/ou_damp)*force->boltz*ou_temp;
    for(int i=0;i<atom->nlocal;i++){
      const double m=(atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]])*force->mvv2e;
      const double sigma=std::sqrt(variance*m)*L;
      for(int k=0;k<3;k++){
        const double noise=ou_random->gaussian();
        store[i][9+3*half+k]=noise;
        store[i][18+k]=decay*store[i][18+k]+sigma*noise;
        atom->v[i][k]=store[i][18+k]/(m*L);
      }
    }
  }
  if(box_ou_enabled){
    const double decay=std::exp(-0.5*update->dt/box_ou_damp);
    const double variance=-std::expm1(-update->dt/box_ou_damp)*box_W*force->boltz*box_ou_temp;
    box_noise[half]=box_random->gaussian();
    const double sigma=npt_test_negative==1?0.0:std::sqrt(variance);
    box_pV=decay*box_pV+sigma*box_noise[half];
  }
}

void FixSpinLPSIB::box_init(){
  box_P=box_Pbar/force->nktv2p;
  if(!box_initialized){
    box_V=domain->xprd*domain->yprd*domain->zprd;
    const double L=std::pow(box_V,1.0/3.0);
    for(int k=0;k<3;k++){
      box_origin[k]=domain->boxlo[k];box_shape[k]=domain->prd[k]/L;
    }
    for(int i=0;i<atom->nlocal;i++){
      double xu[3];domain->unmap(atom->x[i],atom->image[i],xu);
      const double m=(atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]])*force->mvv2e;
      if(!std::isfinite(m) || m<=0) error->one(FLERR,"S85N_POSITIVE_MASS");
      for(int k=0;k<3;k++){
        store[i][15+k]=(xu[k]-box_origin[k])/L;
        store[i][18+k]=m*L*atom->v[i][k];
      }
    }
    box_initialized=1;
  }
  if(comm->me==0 && screen) std::fprintf(screen,
    "S85N_UNITS {\"nktv2p\":%.17g,\"Pbar\":%.17g,\"P_eV_A3\":%.17g,\"Wv\":%.17g,\"V\":%.17g,\"pV\":%.17g}\n",
    force->nktv2p,box_Pbar,box_P,box_W,box_V,box_pV);
}

void FixSpinLPSIB::box_T(double a){
  const double L=std::pow(box_V,1.0/3.0);double twiceK=0.;
  for(int i=0;i<atom->nlocal;i++){
    const double m=(atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]])*force->mvv2e;
    for(int k=0;k<3;k++){
      const double p=store[i][18+k],denom=m*L*L;
      twiceK+=p*p/denom;
      store[i][15+k]+=a*p/denom;
    }
  }
  box_pV+=a*twiceK/(3.0*box_V);
}

void FixSpinLPSIB::box_D(double a){
  box_T(0.5*a);
  box_V+=a*box_pV/box_W;
  if(!std::isfinite(box_V) || box_V<=0) error->all(FLERR,"S85N_NONPOSITIVE_VOLUME");
  box_T(0.5*a);box_sync();
}

void FixSpinLPSIB::box_sync(){
  const double L=std::pow(box_V,1.0/3.0);
  for(int k=0;k<3;k++) domain->boxhi[k]=box_origin[k]+box_shape[k]*L;
  domain->set_global_box();domain->set_local_box();
  for(int i=0;i<atom->nlocal;i++){
    const double m=(atom->rmass?atom->rmass[i]:atom->mass[atom->type[i]])*force->mvv2e;
    for(int k=0;k<3;k++){
      atom->x[i][k]=box_origin[k]+L*store[i][15+k];
      atom->v[i][k]=store[i][18+k]/(m*L);
    }
    atom->image[i]=((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
    domain->remap(atom->x[i],atom->image[i]);
  }
}

void FixSpinLPSIB::box_accept_comm(){
  if(neighbor->check_distance()){
    domain->pbc();domain->reset_box();comm->setup();
    if(neighbor->style) neighbor->setup_bins();
    comm->exchange();comm->borders();
    Observe85::build(neighbor,1,Observe85::ACCEPT);
  } else Observe85::forward(comm,Observe85::ACCEPT);
}
