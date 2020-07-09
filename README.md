# Pigeon

A MPC trajectory tracking controller with ROS integration for a full-scale vehicle (used on Stanford's X1 platform).

This controller is based off ["Safe Driving Envelopes for Path Tracking in Autonomous Vehicles"](https://dynamicdesignlab.sites.stanford.edu/sites/g/files/sbiybj9456/f/brown_safe_driving_envelopes_2016.pdf) by Matthew Brown, Joseph Funke, Stephen Erlien, J. Christian Gerdes. This paper focuses on _path_ tracking, where the longitudinal tracking error is not a major concern, and is decoupled from the lateral dynamics. A PD controller used to track a desired velocity profile. In our work, we care about longitudinal tracking error and have developed a "coupled" MPC tracking controller where the longitudinal and lateral dynamics are considered simulataneously.
Pigeon.jl is implemented in Julia. Unfortunately, this was developed on Julia 1.01, and there are some issues in some of the packages we use that prevent us from upgrading it. 

## DON'T USE THE MASTER BRANCH (will fix it soon.)
The `master` branch has not been updated, but we plan to do that in the future. Instead, please look at the following branches. There are multiple branches because there were different use cases (see the `README.md` from that branch). We plan on making the `master` branch contain the simplest form of Pigeon.jl in the future.

### Branch description
In reference to our IJRR paper ["On Infusing Reachability-Based Safety Assurance within Probabilistic Planning Frameworks for Human-Robot Vehicle Interactions"](http://asl.stanford.edu/wp-content/papercite-data/pdf/Leung.Schmerling.ea.IJRR19.pdf), the following branches were used to perform the analysis from the paper.

Note that in all these branches, they are capable of performing just ordinary trajectory tracking (i.e., without HJI, without walls, etc)


`both_walls`: This branch was used to perform the wall-collision avoidance analysis in the paper. Instructions are found on the README in this branch.

`dynamic_obstacle`: This branch was used to perform the the baseline comparison method of assuming the human car was a dynamic obstacle traveling at the current velocity over the trajectory planning horizon. Instructions are found on the README in this branch.

[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)
