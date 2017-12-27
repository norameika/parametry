import numpy as np


def gen_square():
    csecs = [[[np.cos(t), np.sin(t), - 0.5] for t in np.linspace(0, np.pi * 2, 5)[:-1]],
             [[np.cos(t), np.sin(t), + 0.5] for t in np.linspace(0, np.pi * 2, 5)[:-1]]]
    for csec in csecs:
        yield csec


def gen_pillar():
    csec = [[np.cos(t), np.sin(t), 0] for t in np.linspace(0, np.pi * 2, 181)[:-1]]
    csecs = [np.array(csec) + np.array([0, 0, z]) for z in np.linspace(-2, 2, 2)]

    for csec in csecs:
        yield csec
