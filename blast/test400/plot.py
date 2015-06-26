#!/usr/bin/env python

"""Hack to plot cross-sections of 3D run in python"""

import os

import matplotlib.pyplot as plt

import clawpack.pyclaw.solution
import clawpack.clawutil.data

# Paths
output_path = "_output"
slice_index = [10, 10, 10]

# Load data objects
claw = clawpack.clawutil.data.ClawInputData(3)
claw.read(os.path.join(output_path, "claw.data"))

for n in xrange(claw.num_output_times):
    sol = clawpack.pyclaw.solution.Solution(n, path=output_path)

    fig = plt.figure(figsize=(8*3, 6))
    axes = [fig.add_subplot(1, 3, i + 1) for i in xrange(3)] 
    # for state in sol.states:
    state = sol.states[0]
    X = state.grid.p_centers[0]
    Y = state.grid.p_centers[1]
    Z = state.grid.p_centers[2]

    axes[0].pcolor(X[:, :, slice_index[2]], Y[:, :, slice_index[2]],
                state.q[0, :, :, slice_index[2]], vmin=0.0, vmax=1.0)
    axes[0].set_title("Z-Slice")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")

    axes[1].pcolor(X[:, slice_index[1], :], Z[:, slice_index[1], :],
                state.q[0, :, slice_index[1], :], vmin=0.0, vmax=1.0)
    axes[1].set_title("Y-Slice")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("z")

    axes[2].pcolor(Y[slice_index[1], :, :], Z[slice_index[1], :, :], 
                state.q[0, slice_index[0], :, :], vmin=0.0, vmax=1.0)
    axes[2].set_title("X-Slice")
    axes[2].set_xlabel("y")
    axes[2].set_ylabel("z")

    fig.suptitle("Solution at t = %s" % sol.t)

plt.show()
