# Pigeon

Code used for the wall collision avoidance analysis in the paper.

HJI safety-preserving controller with the option to avoid a wall using HJI_WALL value function, or with state (lateral error) constraints.


## To run:
In the julia terminal:
To activate the julia project environment.
> ]
> activate <path to env>

To load Pigeon
> using Pigeon

To start running the Pigeon ROS node
> Pigeon.start_ROS_node()



You can change the settings to the Pigeon mpc by typing things like:
To switch to the HJI optimal policy (switching controller)
> Pigeon.use_HJI_policy[] = true  

To change the weights on the objective function using something like this:
> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_HJI=500.0, W_erbd=500.0)

For no wall collision avoidance and with HJI safety-preserving controller:
> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_HJI=500.0, W_erbd=0.0, W_WALL=0.0)

For HJI safety preserving wall collision avoidance and with HJI safety-preserving controller:
> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_HJI=500.0, W_erbd=0.0, W_WALL=500.0)

For stte constraint wall collision avoidance and with HJI safety-preserving controller:
> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_HJI=500.0, W_erbd=500.0, W_WALL=0.0)


Note that the non-zero values can be changed to different numbers. These are just weights in the cost function.

To change the buffer around the avoid set, affecting when the safety controller kicks in.
> Pigeon.X1CMPC.HJI_ϵ = 0.05






[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)


[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)
fin