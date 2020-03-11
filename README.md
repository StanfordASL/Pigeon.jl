# Pigeon

This has been tested on Julia 1.0.1, and 1.0.5. For some reason, it's slower on later versions of Julia.

1. download [Julia 1.0.1, or 1.0.5](https://github.com/JuliaLang/julia/)
2. You need the Pigeon `Project.toml` and `Manifest.toml` files.
   Make a folder anywhere you want (e.g., `~/.julia/environmentments/pigeon` or `~/projects/pigeon`)
   Add the `Project.toml` and `Manifest.toml` files from the repo into that folder. You will need to start Julia from the same folder where the `.toml` files live, and 
3. Open a Julia repl, enter package manager, and add the Pigeon package.
   > ]
   
   > activate .
   
   > instantiate
   
   > add git@github.com:StanfordASL/Pigeon.jl.git#vanilla_mpc
   
4. Back in the repl
   > using Pigeon
   
5. Pigeon should be up and running.

If you want to develop Pigeon, you can type
> ]

> develop Pigeon

and you can find the repo now in `~/.julia/dev/Pigeon`.



[![Build Status](https://travis-ci.org/schmrlng/Pigeon.jl.svg?branch=master)](https://travis-ci.org/schmrlng/Pigeon.jl)
[![codecov.io](http://codecov.io/github/schmrlng/Pigeon.jl/coverage.svg?branch=master)](http://codecov.io/github/schmrlng/Pigeon.jl?branch=master)
