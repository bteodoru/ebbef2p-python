import pytest

from ebbef2p.structure import Structure

def test_add_beam():

    s = Structure('test')
    
    s.add_beam(coord=[0, 10],  E=29e6, I=350)
    
    # check parameters of the beam to ensure they match the input
    assert s.beams[0].coord == [0, 10], "beam coords does not match input"
    assert s.beams[0].length == 10, "beam length does not match input"
    assert s.beams[0].E == 29e6, "Young's modulus does not match input"
    assert s.beams[0].I == 350, "area moment of inertia does not match input"

def test_add_nodal_load():

    s = Structure('test')

    s.add_nodal_load(value=10000, position=3, type='fz')

    assert s.nodal_loads[0].value == 10000
    assert s.nodal_loads[0].position == 3
    assert s.nodal_loads[0].type == 'fz'