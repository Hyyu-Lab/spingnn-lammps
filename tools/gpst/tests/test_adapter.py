"""Small fixed interface fixtures and exact rejection-text differences."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pytest
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dump_adapter', ROOT/'dump_to_moments.py')
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


def fixed_dump(spin_length=1., step=0):
    header = f'ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n4\nITEM: BOX BOUNDS pp pp pp\n0 20\n0 20\n0 20\nITEM: ATOMS ' + ' '.join(adapter.FIELDS)+'\n'
    # Explicit unit directions; fm chosen to give force components (2,4,6)
    # under the synthetic hbar=2 conversion. This is not a physical constant.
    rows = [[1,0,0,0,spin_length,0,0,2.2,1,2,3],
            [2,1,0,0,0,1,0,2.2,1,2,3],
            [3,0,1,0,0,0,1,2.2,1,2,3],
            [4,0,0,1,-1,0,0,2.2,1,2,3]]
    return header+'\n'.join(' '.join(map(str,row)) for row in rows)+'\n'


def read(path):
    return adapter.read_dump(path,hbar=2.,accept_every=10,force_contract_confirmed=True)


def test_fixed_fixture_and_conversion(tmp_path):
    path=tmp_path/'tiny.dump';path.write_text(fixed_dump());frame=read(path)
    assert np.array_equal(frame.atom_ids,[1,2,3,4])
    assert np.array_equal(frame.spin_forces,np.broadcast_to([2.,4.,6.],(1,4,3)))
    moments=adapter.tspin.compute_moments(frame,3.)
    assert np.isfinite(moments['N']).all() and np.isfinite(moments['D']).all()
    # Analytic independent pooled-ratio fixture, including exact Fieller CI.
    ratio=adapter.tspin.pooled_ratio([2.,8.],[1.,4.],kb=1.)
    assert ratio['temperature_K']==2.
    ci=adapter.tspin.fieller(2.,1.,np.zeros((2,2)),23)
    assert ci['confidence_set']==[[2.,2.]]


@pytest.mark.parametrize('raw,expected',[(fixed_dump(2.),'NONUNIT_SPIN'),
    (fixed_dump(step=1),'DUMP_NOT_ACCEPTED_STATE'),
    (fixed_dump().replace('c_gp[7]','fmx'),'DUMP_FIELD_CONTRACT')])
def test_negative_difference(tmp_path,raw,expected):
    path=tmp_path/'tiny.dump';path.write_text(fixed_dump());baseline=[]
    try:read(path)
    except ValueError as e:baseline=[str(e)]
    assert baseline==[]
    path.write_text(raw)
    with pytest.raises(ValueError) as exc:read(path)
    assert expected in set([str(exc.value)])-set(baseline)


def test_explicit_contract_and_hbar(tmp_path):
    p=tmp_path/'tiny.dump';p.write_text(fixed_dump());read(p)
    with pytest.raises(ValueError,match='FORCE_CONTRACT_MISSING'):
        adapter.read_dump(p,hbar=2.,accept_every=10,force_contract_confirmed=False)
    with pytest.raises(ValueError,match='BAD_INPUT'):
        adapter.read_dump(p,hbar=0.,accept_every=10,force_contract_confirmed=True)
