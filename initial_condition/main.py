import os
import sys
import numpy as np


def init_time(dest):
    # iterator and time
    step = np.array(0, dtype=np.uint64)
    time = np.array(0, dtype=np.float64)
    np.save(f"{dest}/step.npy", step)
    np.save(f"{dest}/time.npy", time)
    return


def init_domain(lengths, glsizes, dest):
    xf = np.linspace(0., lengths[0], glsizes[0] + 1, endpoint=True)
    xc = 0.
    xc = np.append(xc, 0.5 * xf[:-1] + 0.5 * xf[1:])
    xc = np.append(xc, lengths[0])
    np.save(f"{dest}/xf.npy", np.array(xf, dtype=np.float64))
    np.save(f"{dest}/xc.npy", np.array(xc, dtype=np.float64))
    np.save(f"{dest}/glsizes.npy", np.array(glsizes, dtype=np.uint64))
    np.save(f"{dest}/lengths.npy", np.array(lengths, dtype=np.float64))
    return xc


def init_interface(lengths, glsizes, xc, dest):
    d = 1.0
    r = 0.5 * d
    cx = d
    cy = d
    lx = lengths[0]
    ly = lengths[1]
    nx = glsizes[0]
    ny = glsizes[1]
    dx = lx / nx
    dy = ly / ny
    yc = np.linspace(0.5 * dy, ly - 0.5 * dy, ny)
    xc, yc = np.meshgrid(xc, yc)
    eta = r - np.sqrt(np.power(xc - cx, 2.) + np.power(yc - cy, 2.))
    vof = 0.5 * (1. + np.tanh(0.5 * glsizes[0] * eta))
    np.save(f"{dest}/vof.npy", vof)
    return


def init_fluid(lengths, glsizes, dest):
    lx = lengths[0]
    ly = lengths[1]
    dx = lx / glsizes[0]
    dy = ly / glsizes[1]
    shape0 = (glsizes[1], glsizes[0] + 1)
    shape1 = (glsizes[1], glsizes[0] + 2)
    ux  = np.zeros(shape0, dtype=np.float64)
    uy  = np.zeros(shape1, dtype=np.float64)
    p   = np.zeros(shape1, dtype=np.float64)
    psi = np.zeros(shape1, dtype=np.float64)
    np.save(f"{dest}/ux.npy", ux)
    np.save(f"{dest}/uy.npy", uy)
    np.save(f"{dest}/p.npy", p)
    np.save(f"{dest}/psi.npy", psi)
    return


def main():
    lengths = list()
    lengths.append(float(os.environ["lx"]))
    lengths.append(float(os.environ["ly"]))
    glsizes = list()
    glsizes.append(int(os.environ["glisize"]))
    glsizes.append(int(os.environ["gljsize"]))
    dest = sys.argv[1]
    # sanitise
    ndims = len(lengths)
    assert 2 == ndims
    # init and save
    init_time(dest)
    xc = init_domain(lengths, glsizes, dest)
    init_interface(lengths, glsizes, xc, dest)
    init_fluid(lengths, glsizes, dest)
    return


main()
