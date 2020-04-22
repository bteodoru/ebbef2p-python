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