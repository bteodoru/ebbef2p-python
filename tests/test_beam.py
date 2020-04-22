import pytest

from ebbef2p.beam import Beam

def test_beam_params():
    
    beam = Beam(coord=[0, 10],  E=29e6, I=350)
    
    # check parameters of the beam to ensure they match the input
    assert beam.coord == [0, 10], "beam coords does not match input"
    assert beam.length == 10, "beam length does not match input"
    assert beam.E == 29e6, "Young's modulus does not match input"
    assert beam.I == 350, "area moment of inertia does not match input"