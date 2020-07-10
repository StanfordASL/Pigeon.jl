# Pigeon

This is a "vanilla" Pigeon MPC tracking controller. The HJI and wall collision avoidance functionalities from the other branches are removed. Given a trajectory, this will compute controls to track that trajectory.


## To run:
In the julia terminal:
To activate the julia project environment.
> ]
> activate <path to env\>

To load Pigeon
> using Pigeon

To start running the Pigeon ROS node
> Pigeon.start_ROS_node()



You can change the settings to the Pigeon mpc by typing things like:

> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_r=50.0)

But the current weights in the objective seems to work pretty well, so we don't recommend changing them unless you feel it is necessary.





[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)