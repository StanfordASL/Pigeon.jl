# Pigeon

A MPC trajectory tracking controller with ROS integration for a full-scale vehicle (used on Stanford's X1 platform).

This controller is based off the paper ["Safe Driving Envelopes for Path Tracking in Autonomous Vehicles"](https://dynamicdesignlab.sites.stanford.edu/sites/g/files/sbiybj9456/f/brown_safe_driving_envelopes_2016.pdf) by Matthew Brown, Joseph Funke, Stephen Erlien, J. Christian Gerdes. That paper focuses on _path_ tracking, where the longitudinal tracking error is not a major concern, and is decoupled from the lateral dynamics. A PD controller is used to track a desired velocity profile. 

In our work, we care about longitudinal _and_ lateral tracking error and have developed a "coupled" MPC tracking controller where the longitudinal and lateral dynamics are considered simulataneously.
Pigeon.jl is implemented in Julia. Unfortunately, this was developed on Julia 1.01, and there are some issues in some of the packages we use that prevent us from upgrading it. 

## DON'T USE THE MASTER BRANCH (will fix it soon.)
The `master` branch has not been updated, but we plan to do that in the future. Instead, please look at the following branches. There are multiple branches because there were different use cases (see the `README.md` from that branch). We plan on making the `master` branch contain the simplest form of Pigeon.jl in the future.


## Installing Pigeon

This has been tested on Julia 1.0.1, and 1.0.5. For some reason, it's slower on later versions of Julia.

1. download [Julia 1.0.1, or 1.0.5](https://github.com/JuliaLang/julia/)
2. You need the Pigeon `Project.toml` and `Manifest.toml` files in the env folder.
   Make a folder anywhere you want (e.g., `~/.julia/environmentments/pigeon` or `~/projects/pigeon`)
   Add the `Project.toml` and `Manifest.toml` files from the repo into that folder. You will need to start Julia from the same folder where the `.toml` files live.
3. Open a Julia repl, enter package manager, and add the Pigeon package.
   > ]

   > activate .

   > instantiate

   > add git@github.com:StanfordASL/Pigeon.jl.git
   or for a specific branch 
   > add git@github.com:StanfordASL/Pigeon.jl.git#<branch_name>

4. Back in the repl

   > using Pigeon

5. Pigeon should be up and running.

If you want to develop, or make some changes to Pigeon, you can type
> ]

> develop Pigeon

and you can find the repo now in `~/.julia/dev/Pigeon`. We recommend doing this because you can access the Pigeon repo here, and change branches. You can also perform commits and other git operations.


## Branch description
In reference to our IJRR paper ["On Infusing Reachability-Based Safety Assurance within Probabilistic Planning Frameworks for Human-Robot Vehicle Interactions"](http://asl.stanford.edu/wp-content/papercite-data/pdf/Leung.Schmerling.ea.IJRR19.pdf), the following branches were used to perform the analysis from the paper.

Note that in all these branches, they are capable of performing just ordinary trajectory tracking (i.e., without HJI, without walls, etc)

`vanilla_mpc`: This is a vanilla MPC tracking controller without any extra functionalities such as the HJI safety controller, or wall collision avoidance. If you just want a trajectory tracking controller with the vehicle model described in our [IJRR paper]((http://asl.stanford.edu/wp-content/papercite-data/pdf/Leung.Schmerling.ea.IJRR19.pdf)), this this is the package for you. Instructions are found on the README in this branch. We plan on making a tidier version of this branch the `master` branch.

`both_walls`: This branch was used to perform the wall-collision avoidance analysis in the paper. Instructions are found on the README in this branch.

`dynamic_obstacle`: This branch was used to perform the the baseline comparison method of assuming the human car was a dynamic obstacle traveling at the current velocity over the trajectory planning horizon. Instructions are found on the README in this branch.


## Future work and TODOs
We have primarily used Pigeon.jl as a ROS node and have been running it on X1, a full-scale test vehicle. Pigeon can be used in simulation without ROS. We plan on providing more documentation on this in the future. In the meanwhile, see `model_predictive_control.jl` for an idea.




[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)
