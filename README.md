# StokesVEP_2D

This is a 2D finite difference, staggered grid code, solving the conservation of momentum equation using an iterative pseudo-transient solver for different kinematic boundary conditions (so far: pure and simple shear, i.e. pureshear=true or pureshear=false, respectively) with a focus on the effect of different strain-dependent weakening and hardening (SDWH) rheologies on the localization of deformation and the formation of shear bands. 

On can choose either a plastic-strain softening (PSS, i.e., PSS=true) or a viscous-strain softening (VSS, i.e., PSS=false) rheology. 

The initial condition can be defined as 
    1) a gaussian strain distribution with a maximum at the center (0,0) of the model domain or
    2) a randomly distributed strain field with a maximum amplitude \gamma_0.
