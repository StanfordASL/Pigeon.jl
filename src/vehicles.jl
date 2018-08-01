function X1()
    X1 = Dict{Symbol,Float64}()

    # Mass and Moments of Inertia
    X1[:mfl] = 484                                                     # weight at front left wheel (kg)
    X1[:mfr] = 455                                                     # weight at front right wheel (kg)
    X1[:mrl] = 521                                                     # weight at rear left wheel (kg)
    X1[:mrr] = 504                                                     # weight at rear right wheel (kg)
    X1[:m] = X1[:mfl] + X1[:mfr] + X1[:mrl] + X1[:mrr]                 # total vehicle mass (kg)
    X1[:Ixx] = 175                                                     # pitch moment of inertia (kg*m^2)
    X1[:Iyy] = 1000                                                    # roll moment of inertia (kg*m^2)
    X1[:Izz] = 2900                                                    # yaw moment of inertia (kg*m^2)

    # Dimensions
    X1[:L] = 2.87                                                      # wheelbase (m)
    X1[:d]  = 1.63                                                     # track width (m)
    X1[:a] = (X1[:mrl] + X1[:mrr])/X1[:m]*X1[:L]                       # distance from CG to front axle (m)
    X1[:b] = (X1[:mfl] + X1[:mfr])/X1[:m]*X1[:L]                       # distance from CG to rear axle (m)
    X1[:hf] = 0.1                                                      # front roll center height (m)
    X1[:hr] = 0.1                                                      # rear roll center height (m)
    X1[:h1] = 0.37                                                     # CG height w.r.t. line connecting front and rear roll centers (m)
    X1[:h] = X1[:hf]*X1[:b]/X1[:L] + X1[:hr]*X1[:a]/X1[:L] + X1[:h1]   # CG height (m)
    X1[:ab] = X1[:a] + 0.4953                                          # distance from CG to front bumper (m)
    X1[:bb] = X1[:b] + 0.5715                                          # distance from CG to rear bumper (m)
    X1[:w] = 1.87                                                      # physical width (m)

    # Tire Model Parameters
    X1[:μ] = 0.92                                                      # coefficient of friction
    X1[:Cαf] = 140e3                                                   # front tire (pair) cornering stiffness (N/rad)
    X1[:Cαr] = 190e3                                                   # rear tire (pair) cornering stiffness (N/rad)

    # Longitudinal Actuation Parameters
    X1[:maxFx] = 5600                                                  # max positive longitudinal force (N)
    X1[:maxPower] = 75e3                                               # max motor power output (W)

    # Longitudinal Drag Force Parameters (FxDrag = Cd0 + Cd1*Ux + Cd2*Ux^2)
    X1[:Cd0] = 241.0                                                   # rolling resistance (N)
    X1[:Cd1] = 25.1                                                    # linear drag coefficint (N/(m/s))
    X1[:Cd2] = 0.0                                                     # quadratic "aero" drag coefficint (N/(m/s)^2)

    # Drive and Brake Distribution
    X1[:fwd_frac] = 0.0                                                # front wheel drive fraction for implementing desired Fx
    X1[:rwd_frac] = 1 - X1[:fwd_frac]                                  # rear wheel drive fraction for implementing desired Fx
    X1[:fwb_frac] = 0.6                                                # front wheel brake fraction for implementing desired Fx
    X1[:rwb_frac] = 1 - X1[:fwb_frac]                                  # rear wheel brake fraction for implementing desired Fx


    ## Janky HJI stuff
    # Maximum Brake Force
    X1[:minFx] = max(-X1[:m]*GRAVITY*X1[:a]*X1[:μ]/(X1[:L]*X1[:rwb_frac] + X1[:μ]*X1[:h]),
                     -X1[:m]*GRAVITY*X1[:b]*X1[:μ]/(X1[:L]*X1[:fwb_frac] - X1[:μ]*X1[:h]))
    X1[:maxAx] = X1[:maxFx]/X1[:m]
    X1[:minAx] = X1[:minFx]/X1[:m]
    X1[:maxP2mx] = X1[:maxPower]/X1[:m]

    # Maximum Steering
    X1[:d_max] = 18*pi/180;
    X1[:w_per_v_max_lowspeed] = tan(X1[:d_max])/X1[:L]

    X1[:maxA_approx] = 0.9*X1[:μ]*GRAVITY

    X1
end
