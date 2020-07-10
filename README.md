# Pigeon

Code used for the baseline comparison analysis in the paper.

The baseline comparison technique with HJI-based safety-preserving control constraints.
The other (human) car is treated as a dynamic obstacle that will continue moving at it's current velocity. The robot car applies state (i.e. lateral error) constraints along the trajectory to stay away from this projection.


## To run:
In the julia terminal:
To activate the julia project environment.
> ]\
> activate <path to env\>

To load Pigeon
> using Pigeon

To start running the Pigeon ROS node
> Pigeon.start_ROS_node()



You can change the settings to the Pigeon mpc by typing things like:

To switch to the HJI optimal policy (switching controller)
> Pigeon.use_HJI_policy[] = true  

To change the weights on the objective function using something like this:
> Pigeon.X1CMPC.control_params = Pigeon.CoupledControlParams(W_HJI=0.0, W_erbd=500.0)

By default, `W_HJI = 500.0` and `W_erbd = 0.0` which corresponds to using the HJI safety controller but with no wall collision avoidance.

`W_HJI = 500.0` and `W_erbd = 500.0` which corresponds to using the HJI safety controller and using state constraint for wall collision avoidance.

`W_HJI = 0.0` and `W_erbd = 500.0` which corresponds to _not_ using the HJI safety controller, but using state constraints to avoid the other car and wall by making a constant velocity dynamical obstacle assumption (baseline method).

To change the buffer around the avoid set, affecting when the safety controller kicks in.
> Pigeon.X1CMPC.HJI_ϵ = 0.05






[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)
