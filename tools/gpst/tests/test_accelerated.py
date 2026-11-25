"""Deterministic contracts and algebra; canonical calibration is run separately."""
import dataclasses
import inspect
from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.special import eval_legendre

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'accelerated'))
import tspin_v1_1 as tspin
from tspin_v1_1.geometry import FIXED, legendre, minimum_image, moment_arrays
from tspin_v1_1.basis import learn_basis, evaluate_basis, upper_bounds
from tspin_v1_1.inference import crossfit, pooled_pair
from source_guard import errors as source_errors

POS = np.array([[0.,0.,0.],[1.,0.,0.],[.2,1.1,0.],[.1,.3,1.3]])
BOX = np.eye(3)*20
RC = 3.4


def unit(x):
    return x/np.linalg.norm(x, axis=-1, keepdims=True)


def arrays(n=4, frames=12, seed=80):
    rng = np.random.default_rng(seed)
    return dict(atom_ids=np.arange(1,n+1), magnetic_mask=np.ones(n,bool),
                positions=np.broadcast_to(POS,(frames,n,3)).copy(),
                cell=np.broadcast_to(BOX,(frames,3,3)).copy(), pbc=np.zeros(3,bool),
                spins=unit(rng.normal(size=(frames,n,3))),
                spin_forces=rng.normal(size=(frames,n,3)), frame_index=np.arange(frames),
                force_contract_confirmed=True)


def phi_grad(s, name, graph):
    i,j,w,u = [graph[k] for k in ('i','j','w','direction')]
    grad = np.zeros_like(s)
    if name.startswith('pair_'):
        order = int(name[6:]); p,dp = legendre(np.sum(s[i]*s[j],axis=1),order)
        phi = w@p[order]
        np.add.at(grad,i,(w*dp[order])[:,None]*s[j])
        np.add.at(grad,j,(w*dp[order])[:,None]*s[i])
        d = 2*order*(order+1)*phi
    elif name == 'geometry_Q2':
        grad = np.einsum('nij,nj->ni',graph['Q'],s)
        phi = .5*np.sum(grad*s)
        d = 6*phi-np.trace(graph['Q'],axis1=1,axis2=2).sum()
    elif name == 'geometry_zonal4':
        phi = 0.
        for idx in (i,j):
            p,dp = legendre(np.sum(s[idx]*u,axis=1),4)
            phi += w@p[4]
            np.add.at(grad,idx,(w*dp[4])[:,None]*u)
        d = 20*phi
    else:
        phi = np.sum(w*np.sum(u*np.cross(s[i],s[j]),axis=1))
        np.add.at(grad,i,w[:,None]*np.cross(s[j],u))
        np.add.at(grad,j,w[:,None]*np.cross(u,s[i]))
        d = 4*phi
    return phi, grad-s*np.sum(s*grad,axis=1,keepdims=True), d


@pytest.mark.parametrize('name',FIXED)
def test_gradient_divergence_h2(name):
    rng=np.random.default_rng(801);s=unit(rng.normal(size=(4,3)));f=rng.normal(size=s.shape)
    g=tspin.build_graph(POS,BOX,np.zeros(3,bool),RC)
    phi,x,d=phi_grad(s,name,g);r=moment_arrays(s[None],f[None],g);idx=FIXED.index(name)
    assert_allclose(r['N'][0,idx],np.sum(x*f),rtol=1e-12,atol=1e-12)
    assert_allclose(r['D'][0,idx],d,rtol=1e-12,atol=1e-12)
    errs=[]
    for h in (2e-3,1e-3,5e-4):
        div=0.;grad_error=0.
        for site in range(4):
            e=unit(np.cross(s[site],np.eye(3)[np.argmin(abs(s[site]))]))
            for tangent in (e,np.cross(s[site],e)):
                sp=s.copy();sm=s.copy()
                sp[site]=np.cos(h)*s[site]+np.sin(h)*tangent
                sm[site]=np.cos(h)*s[site]-np.sin(h)*tangent
                pp,xp,_=phi_grad(sp,name,g);pm,xm,_=phi_grad(sm,name,g)
                div += tangent@(xp[site]-xm[site])/(2*h)
                grad_error=max(grad_error,abs((pp-pm)/(2*h)-x[site]@tangent))
        errs.append(abs(div+d)/max(1,abs(d)))
        assert grad_error < .002
    assert errs[-1] < 1e-4
    assert all(a/b > 3.8 for a,b in zip(errs,errs[1:]))


@pytest.mark.parametrize('degree',[1,2,4,8])
def test_legendre(degree):
    x=np.linspace(-1,1,201);p,dp=legendre(x)
    assert_allclose(p[degree],eval_legendre(degree,x),rtol=1e-12,atol=1e-12)
    c=np.zeros(degree+1);c[-1]=1
    assert_allclose(dp[degree],np.polynomial.legendre.legval(x,np.polynomial.legendre.legder(c)),rtol=1e-12,atol=1e-12)


@pytest.mark.parametrize('operation',['rotation','translation','id_permutation','type_rename','triclinic','radial'])
def test_invariance(operation):
    a=arrays();base=tspin.compute_moments(tspin.validate_frame(**a),RC)
    b={k:v.copy() if isinstance(v,np.ndarray) else v for k,v in a.items()}
    if operation=='rotation':
        q,_=np.linalg.qr(np.random.default_rng(80).normal(size=(3,3)))
        for key in ('positions','spins','spin_forces','cell'):b[key]=b[key]@q
    elif operation=='translation':b['positions']+=np.array([1.25,-2.75,.83])
    elif operation=='id_permutation':
        order=np.array([2,0,3,1])
        for key in ('atom_ids','magnetic_mask'):b[key]=b[key][order]
        for key in ('positions','spins','spin_forces'):b[key]=b[key][:,order]
    elif operation=='type_rename':
        # Types intentionally cannot enter the API. External labels can change.
        labels=np.array([17,29,17,29]);renamed=np.array([113,7,113,7]);assert not np.array_equal(labels,renamed)
    elif operation=='triclinic':b['cell'][:]=np.array([[20.,0,0],[8.,20.,0],[3.,2.,20.]])
    else:b['spin_forces']+=3*b['spins']
    got=tspin.compute_moments(tspin.validate_frame(**b),RC)
    for key in ('N','D'):
        assert np.max(abs(got[key]-base[key])/np.maximum(1,abs(base[key]))) <= 1e-12


def test_periodic_bruteforce_and_supercell():
    rng=np.random.default_rng(81);box=np.array([[7.,0,0],[3.5,7.,0],[1.,2.,7.]])
    pos=rng.random((24,3))@box;rc=2.1;periodic=np.ones(3,bool)
    g=tspin.build_graph(pos,box,periodic,rc)
    shifts=np.array(list(product(range(-2,3),repeat=3)))@box
    expected={}
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            candidates=pos[j]-pos[i]-shifts;v=candidates[np.argmin(np.sum(candidates*candidates,axis=1))]
            if np.linalg.norm(v)<rc:expected[i,j]=v
    assert set(expected)==set(zip(g['i'],g['j']))
    assert_allclose(g['distance'],[np.linalg.norm(expected[i,j]) for i,j in zip(g['i'],g['j'])],atol=1e-12)
    s=unit(rng.normal(size=(3,24,3)));f=rng.normal(size=s.shape);r=moment_arrays(s,f,g)
    pos2=np.concatenate([pos,pos+box[0]]);box2=box.copy();box2[0]*=2
    g2=tspin.build_graph(pos2,box2,periodic,rc,shell_cuts=g['shell_cuts'])
    r2=moment_arrays(np.tile(s,(1,2,1)),np.tile(f,(1,2,1)),g2)
    for key in ('N','D'):assert_allclose(r2[key],2*r[key],rtol=1e-12,atol=1e-12)
    assert_allclose(r2['N'].sum(0)/r2['D'].sum(0),r['N'].sum(0)/r['D'].sum(0),rtol=1e-12)
    dr=rng.normal(size=(25,3))*8
    exact=minimum_image(dr,box,periodic)
    large=np.array(list(product(range(-6,7),repeat=3)))@box
    brute=np.min(np.sum((dr[:,None]-large)**2,axis=2),axis=1)
    assert_allclose(np.sum(exact*exact,axis=1),brute,rtol=1e-12,atol=1e-12)


def test_multispecies_nonmagnetic_variable_cell():
    a=arrays();a['atom_ids']=np.array([4,2,7,1,9,8]);a['magnetic_mask']=np.array([1,0,1,1,0,1],bool)
    a['positions']=np.insert(a['positions'],[1,3],[[8.,8.,8.],[9.,9.,9.]],axis=1)
    a['cell']*=np.linspace(1,1.02,len(a['cell']))[:,None,None]
    # Two externally named magnetic species share the same directional contract.
    frame=tspin.validate_frame(**a);r=tspin.compute_moments(frame,RC)
    assert r['N'].shape==(12,7) and np.isfinite(r['N']).all()
    for k in range(12):
        g=tspin.build_graph(frame.positions[k,frame.magnetic_mask],frame.cell[k],frame.pbc,RC,shell_cuts=r['shell_cuts'])
        single=moment_arrays(frame.spins[k:k+1],frame.spin_forces[k:k+1],g)
        assert_allclose(r['N'][k],single['N'][0],atol=1e-12)


@pytest.mark.parametrize('fault,code',[('unit','NONUNIT_SPIN'),('nan','BAD_INPUT'),('inf','BAD_INPUT'),('id','ID_MISMATCH'),('shape','ID_MISMATCH'),('force','FORCE_CONTRACT_MISSING')])
def test_rejections(fault,code):
    a=arrays()
    if fault=='unit':a['spins']*=2
    elif fault=='nan':a['spin_forces'][0,0,0]=np.nan
    elif fault=='inf':a['positions'][0,0,0]=np.inf
    elif fault=='id':a['atom_ids'][1]=a['atom_ids'][0]
    elif fault=='shape':a['cell']=a['cell'][0]
    else:a['force_contract_confirmed']=False
    with pytest.raises(ValueError,match=code):tspin.validate_frame(**a)


def test_direct_dataclass_cannot_bypass_checks():
    a=arrays();a['spins']*=2
    with pytest.raises(ValueError,match='NONUNIT_SPIN'):tspin.compute_moments(tspin.FrameInput(**a),RC)


def test_no_reference_temperature_or_default_force(tmp_path):
    a=arrays();a.pop('force_contract_confirmed');np.savez(tmp_path/'frames.npz',**a)
    (tmp_path/'meta.json').write_text('malformed: must not be read')
    frame=tspin.load_frame_npz(tmp_path/'frames.npz',force_contract_confirmed=True)
    assert len(frame.spins)==12
    for func in (tspin.validate_frame,tspin.compute_moments,tspin.analyze,tspin.load_frame_npz):
        assert not {'temperature','temperature_K','T','theta'}.intersection(inspect.signature(func).parameters)
    with pytest.raises(TypeError):tspin.from_lammps_fm(np.ones(3))
    assert_allclose(tspin.from_lammps_fm(np.ones(3),.123),.123)


def test_pooled_and_fieller():
    n=np.array([2.,8.]);d=np.array([1.,4.]);assert tspin.pooled_ratio(n,d,kb=1)['temperature_K']==2
    with pytest.raises(ValueError,match='only pooled'):tspin.pooled_ratio(n,d,kb=1,method='frame_average')
    a,b=pooled_pair([np.array([[2.]]),np.full((3,1),8.)],[np.array([[1.]]),np.full((3,1),2.)])
    assert_allclose(a/b,26/7)
    assert tspin.fieller(0,0,np.zeros((2,2)),23)['status']=='ZERO_SIGNAL'
    assert tspin.fieller(1,0,np.zeros((2,2)),23)['status']=='EMPTY_CI'
    assert tspin.fieller(1,0,np.eye(2),23)['status']=='WEAK_DENOMINATOR'
    assert len(tspin.fieller(4,.1,np.diag([.01,1.]),23)['confidence_set'])==2


def test_zero_signal_and_bounds():
    a=arrays();a['spin_forces'][:]=0;r=tspin.compute_moments(tspin.validate_frame(**a),RC)
    out=tspin.analyze([r]*6,list(range(6)),kb=1)
    assert all(c['status']=='ZERO_SIGNAL' and c['temperature_K'] is None for c in out['channels'].values())
    basis=learn_basis([r]);_,d=evaluate_basis(r,basis);u=upper_bounds(r,basis)
    assert np.max(abs(d)-u)<=1e-12
    # Cubic bond geometry gives an isotropic Q on each site, hence exact null.
    pos=np.array(list(product(range(4),repeat=3)),float);g=tspin.build_graph(pos,np.eye(3)*4,np.ones(3,bool),1.1)
    rng=np.random.default_rng(82);s=unit(rng.normal(size=(16,64,3)));f=rng.normal(size=s.shape)
    rr=moment_arrays(s,f,g);assert np.max(abs(rr['N'][:,4]))<1e-12 and np.max(abs(rr['D'][:,4]))<1e-12


def test_crossfit_training_isolation():
    rng=np.random.default_rng(83);g=tspin.build_graph(POS,BOX,np.zeros(3,bool),RC)
    rows=[]
    for _ in range(6):
        s=unit(rng.normal(size=(512,4,3)));s[:,:,2]+=1;s=unit(s);f=np.zeros_like(s);f[:,:,2]=1
        rows.append(moment_arrays(s,f,g))
    _,_,before=crossfit(rows,np.arange(6)%2)
    changed=[{k:v.copy() for k,v in r.items()} for r in rows]
    for i in (0,2,4):changed[i]['FS']*=19
    _,_,after=crossfit(changed,np.arange(6)%2)
    assert before[0]==after[0]


def test_source_guard():
    assert source_errors(ROOT/'accelerated/tspin_v1_1')==[]


def test_uniform_incompatible_force():
    # Exact degree-two spherical design, with an attractive force but uniform
    # independent directions: positive N, exactly zero population D.
    axes=np.concatenate([np.eye(3),-np.eye(3)])
    s=np.array(list(product(axes,repeat=4)))
    f=s.sum(1)[:,None]-s
    g=tspin.build_graph(POS,BOX,np.zeros(3,bool),RC)
    r=moment_arrays(s,f,g)
    assert r['N'][:,0].mean()>0
    assert abs(r['D'][:,0].mean())<1e-12
    rejection=tspin.fieller(r['N'][:,0].mean(),r['D'][:,0].mean(),np.zeros((2,2)),23)
    assert rejection['status']=='EMPTY_CI'
    out=tspin.analyze([r]*6,list(range(6)),kb=1)
    assert out['channels']['CF_GMM']['status']=='MOMENT_MODEL_MISFIT'


def test_relative_cutoff_and_stale_graph():
    a=arrays();frame=tspin.validate_frame(**a)
    g=tspin.build_graph(POS,BOX,np.zeros(3,bool),{'nearest_multiplier':3.})
    assert_allclose(g['rc'],3.)
    bad=dataclasses.replace(frame,positions=frame.positions+1)
    with pytest.raises(ValueError,match='ID_MISMATCH'):tspin.compute_moments(bad,3.,graph=g)
    with pytest.raises(ValueError,match='unique image'):tspin.build_graph(POS,np.eye(3)*4,np.ones(3,bool),2.)


def test_same_lattice_triclinic_representation():
    rng=np.random.default_rng(8001);box=np.eye(3)*10
    shear=np.array([[1,0,0],[1,1,0],[0,0,1]])
    pos=rng.random((80,3))@box
    a=tspin.build_graph(pos,box,np.ones(3,bool),2.)
    b=tspin.build_graph(pos,shear@box,np.ones(3,bool),2.,shell_cuts=a['shell_cuts'])
    assert np.array_equal(a['i'],b['i']) and np.array_equal(a['j'],b['j'])
    s=unit(rng.normal(size=(8,80,3)));f=rng.normal(size=s.shape)
    x=moment_arrays(s,f,a);y=moment_arrays(s,f,b)
    for key in ('N','D'):assert np.max(abs(x[key]-y[key])/np.maximum(1,abs(x[key])))<=1e-12


@pytest.mark.parametrize('dtype',[np.int64,np.uint64])
def test_descending_frame_indices_rejected(dtype):
    a=arrays()
    a['frame_index']=np.arange(12,0,-1,dtype=dtype)
    with pytest.raises(ValueError,match='ID_MISMATCH'):
        tspin.validate_frame(**a)


@pytest.mark.parametrize('cutoff',[2.,{'nearest_multiplier':2.}])
def test_prepared_graph_cutoff_mismatch_rejected(cutoff):
    frame=tspin.validate_frame(**arrays())
    graph=tspin.build_graph(POS,BOX,np.zeros(3,bool),RC)
    with pytest.raises(ValueError,match='ID_MISMATCH'):
        tspin.compute_moments(frame,cutoff,graph=graph)


@pytest.mark.parametrize('dtype',[np.int64,np.uint64])
def test_increasing_frame_indices_and_equivalent_cached_cutoffs(dtype):
    a=arrays();a['frame_index']=np.arange(12,dtype=dtype)
    frame=tspin.validate_frame(**a)
    graph=tspin.build_graph(POS,BOX,np.zeros(3,bool),RC)
    expected=tspin.compute_moments(frame,RC)
    actual=tspin.compute_moments(frame,{'nearest_multiplier':RC},graph=graph)
    for key in ('N','D','shell_N','shell_D','SS','FS','shell_SS','shell_AN'):
        assert_allclose(actual[key],expected[key],rtol=0,atol=0)
