# Initialization
using Plots, LazyGrids, Printf, LinearAlgebra, LoopVectorization, MAT
import Statistics: mean 

@views function av(Aout,Ain) 
        Aout .= 0.25.*(Ain[1:end-1,1:end-1].+Ain[2:end,1:end-1].+
                    Ain[1:end-1,2:end].+Ain[2:end,2:end])                         
end
    
@views function main(N,η_reg,γ0)
do_save     =   true
pureshear   =   true
PSS         =   false
γinc        =   false
Ddir        =   "D:/Users/lukas/Numerics/BACKUP/progs/src/MATLAB/Projects/2D_VEP_SDW_reg"
if PSS
    Type        =   "PSS"
else
    Type        =   "VSS"
end
if pureshear
    Def     =   "PS"
else
    Def     =   "SS"
end
# Model constants ====================================================
# Model dimensions ---------------------------------------------------
xmin        =   0.0
xmax        =   1.0
ymin        =   .0
ymax        =   1.0
# Viscous inclusion --------------------------------------------------
# Inclusion radius ---------------------------------------------------
ri          =   .3
# Viscosity ----------------------------------------------------------
ηi          =   1.0
ηm          =   1.0
## Plastic constants --------------------------------------------------
τyield0     =   175        # Background yield stress
η_reg       =   1.2e-2      # Regularization viscosity
# Elastic constants --------------------------------------------------
G           =   10           # Elastic shear module
# Kinematic Boundary conditions --------------------------------------
εbg         =   -1.0        # Background strain rate
# ================================================================== #
# Strain dependent weakening =========================================
Dmax        =   0.9
γcr         =   10
# γ0          =   0.0
# Strain hardening ================================================= #
T           =   1.0;          # Non-dimensional temperature
ηγ          =   82.8931     # Temperature-dependent healing constant
B           =   2.44e9      # Healing time scale
H           =   0 # B * exp( -ηγ/2 * ( 1/(T+1) - 1/2 ) )
@printf("H(T) = %2.2e\n",H)
# ================================================================== #
# Time constants =====================================================
nt          =   400                # Number of iterations
t           =   .0                 # time
# ================================================================== #
# Numerical parameters ===============================================
ncx         =   N 
ncy         =   N 
Δx, Δy      =   (xmax-xmin)/(ncx+1), (ymax-ymin)/(ncy+1)
success     =   false
# Vertices coordinates -----------------------------------------------
xv          =   LinRange(xmin,xmax,ncx+1)
yv          =   LinRange(ymin,ymax,ncy+1)
yvx         =   LinRange(ymin - Δy/2.0,ymax + Δy/2.0,ncy+2)
xvy         =   LinRange(xmin - Δx/2.0,xmax + Δx/2.0,ncx+2)
# Center coordinates -------------------------------------------------
xc          =   0.5*(xv[1:end-1].+xv[2:end])
yc          =   0.5*(yv[1:end-1].+yv[2:end])
# Mesh grid ----------------------------------------------------------
(xc2,yc2)   =   ndgrid(xc,yc)
(xv2,yv2)   =   ndgrid(xv,yv)
(xvx2,yvx2) =   ndgrid(xv,yvx)
(xvy2,yvy2) =   ndgrid(xvy,yv)
# ================================================================== #
# Memory allocation ==================================================
# Kinematic ----------------------------------------------------------
P       =   zeros(ncx+2,ncy)        # Pressure
Vx      =   zeros(ncx+1,ncy+2)      # Horizontal velocity
Vy      =   zeros(ncx+2,ncy+1)      # Vertical velocity
Vxc     =   zeros(ncx,ncy)          # Horizontal velocity; centroids
Vyc     =   zeros(ncx,ncy)          # Vertical velocity; centroids
Vc      =   zeros(ncx,ncy)          # absolut velocity; centroids
dvxdx   =   zeros(ncx,ncy)          # Horizontal velocity gradient
dvydy   =   zeros(ncx,ncy)          # Vertical velocity gradient
dvxdy   =   zeros(ncx+1,ncy+1)      # Horizontal shear velocity gradient
dvydx   =   zeros(ncx+1,ncy+1)      # Vertical shear velocity gradient
# --------------------------------------------------------------------
εxx     =   zeros(ncx,ncy)          # Normal horizontal strain rate; centroids
εyy     =   zeros(ncx,ncy)          # Normal vertical strain rate; centroids
εxy     =   zeros(ncx,ncy)          # Shear strain rate; centroids
εII     =   zeros(ncx,ncy)          # Second invariant; centroids
εxxv    =   zeros(ncx+1,ncy+1)      # Normal strain rate; vertices
εyyv    =   zeros(ncx+1,ncy+1)      # Normal strain rate; vertices
εxyv    =   zeros(ncx+1,ncy+1)      # Shear strain rate; vertices
εIIv    =   zeros(ncx+1,ncy+1)      # Second invariant; vertices

εxxeff  =   zeros(ncx,ncy)          # Effective horizontal strain rate; centroids
εyyeff  =   zeros(ncx,ncy)          # Effective vertical strain rate; centroids
εxyeff  =   zeros(ncx,ncy)          # Effective shear strain rate; centroids
εIIeff  =   zeros(ncx,ncy)          # Effective shear strain rate; centroids
εxyeffv =   zeros(ncx+1,ncy+1)      # Effective shear strain rate; vertices 
εxxpl   =   zeros(ncx,ncy)          # Plastic strain rate; centroids
εyypl   =   zeros(ncx,ncy)          # Plastic strain rate; centroids
εxypl   =   zeros(ncx,ncy)          # Plastic strain rate; centroids
εxxplv  =   zeros(ncx+1,ncy+1)      # Plastic strain rate; vertices
εyyplv  =   zeros(ncx+1,ncy+1)      # Plastic strain rate; vertices
εxyplv  =   zeros(ncx+1,ncy+1)      # Plastic strain rate; vertices
εIIpl   =   zeros(ncx,ncy)          # Plastic strain rate; centroids
εIIplv  =   zeros(ncx+1,ncy+1)      # Plastic strain rate; vertices

εxxel   =   zeros(ncx+1,ncy+1)      # Elastic horizontal strain rate; vertices
εyyel   =   zeros(ncx+1,ncy+1)      # Elastic vertical strain rate; vertices
εxyel   =   zeros(ncx+1,ncy+1)      # Elastic shear strain rate; vertices 
εIIel   =   zeros(ncx+1,ncy+1)      # Elastic strain rate; vertices
εxxvi   =   zeros(ncx+1,ncy+1)      # Viscous horizontal strain rate; vertices
εyyvi   =   zeros(ncx+1,ncy+1)      # Viscous vertical strain rate; vertices
εxyvi   =   zeros(ncx+1,ncy+1)      # Viscous shear strain rate; vertices 
εIIvi   =   zeros(ncx+1,ncy+1)      # Viscous strain rate; vertices
εIIvic  =   zeros(ncx,ncy)          # Viscous strain rate; centroids
εcheck  =   zeros(ncx+1,ncy+1)      # 

τxx     =   zeros(ncx+2,ncy)        # Normal horizontal stress; centroids
τyy     =   zeros(ncx,ncy)          # Normal vertical stress; centroids
τxy     =   zeros(ncx,ncy)          # Shear stress; centroids
τII     =   zeros(ncx,ncy)          # Second invariant; centroids
τIIeff  =   zeros(ncx,ncy)          # Effective strain rate; vertices
τxxv    =   zeros(ncx+1,ncy+1)      # Normal horizontal stress; vertices
τyyv    =   zeros(ncx+1,ncy+1)      # Normal vertical stress; vertices
τxyv    =   zeros(ncx+1,ncy+1)      # Shear stress; vertices

τxxo    =   zeros(ncx+2,ncy)        # Old normal stress - x 
τyyo    =   zeros(ncx,ncy)          # Old normal stress - y
τxyo    =   zeros(ncx,ncy)          # Old shear stress
τxxov   =   zeros(ncx+1,ncy+1)      # Old normal stress - x 
τyyov   =   zeros(ncx+1,ncy+1)      # Old normal stress - y
τxyov   =   zeros(ncx+1,ncy+1)      # Old shear stress

dQdτxx  =   zeros(ncx,ncy)          #
dQdτyy  =   zeros(ncx,ncy)          #
dQdτxy  =   zeros(ncx,ncy)          #
# Rheological---------------------------------------------------------
ηc      =   zeros(ncx,ncy)          # Effective viscosity; centroids
ηv      =   zeros(ncx+1,ncy+1)      # Effective viscosity; vertices
ηBC     =   zeros(ncy-1)            # For BC
ηvi     =   zeros(ncx+1,ncy+1)      # Viscous viscosity; vertices
F       =   zeros(ncx,ncy)          # Yield function; vertices
Fchk    =   zeros(ncx,ncy)          # 
λp      =   zeros(ncx,ncy)          # 
ispl    =   zeros(ncx,ncy)          # Plasticity index
τyield  =   zeros(ncx,ncy)          # Yield stress
# Iterative ----------------------------------------------------------
∇V      =   zeros(ncx,ncy)          #
Fp      =   zeros(ncx,ncy)          # 
Fx      =   zeros(ncx+1,ncy)        # 
Fy      =   zeros(ncx,ncy+1)        # 
dvxdτ   =   zeros(ncx+1,ncy)        # 
dvydτ   =   zeros(ncx,ncy+1)        #
Δτv     =   zeros(ncx+1,ncy+1)      # 
Δτvx    =   zeros(ncx+1,ncy)        # 
Δτvy    =   zeros(ncx,ncy+1)        # 
κΔτp    =   zeros(ncx,ncy)          # 
# Damage -------------------------------------------------------------
γc      =   zeros(ncx,ncy)          # Apparent strain; centroids
Dc      =   zeros(ncx,ncy)          # Damage
Dv      =   zeros(ncx+1,ncy+1)      # Damage; vertices
# Time parameter -----------------------------------------------------
τmean   =   zeros(nt)               # Mean stress
τmax    =   zeros(nt)               # Maximum stress
time    =   zeros(nt)               # Time
γmean   =   zeros(nt)               # Mean strain
γmax    =   zeros(nt)               # Maximum strain
# ================================================================== #
# Initialization =====================================================
# velocity boundary conditions ---------------------------------------
if pureshear 
    Vx      .=  -εbg.*xvx2
    Vy      .=  εbg.*yvy2
else
    Vx      .=  2 .* εbg.*yvx2  
    VBCN    =   2 .* εbg*ymax
    VBCS    =   2 .* εbg*ymin
end
# Viscous inclusion --------------------------------------------------
α           =   0.0
a_ell,b_ell =   1.0, 1.0
x_ell       =   xv2 .* cos(α) .+ yv2 .* sin(α) 
y_ell       =   -xv2 .* sin(α) .+ yv2 .* cos(α) 
Elli        =   (x_ell ./ a_ell).^2 .+ (y_ell ./ b_ell).^2
ηvi                 .=  ηm
ηvi[Elli .< ri.^2]  .=  ηi
# --------------------------------------------------------------------
# Initial strain -----------------------------------------------------
if γinc
    γc      .=  γ0 .* exp.( .- (xc2.^2 .+ yc2.^2) ./ ri .^ 2 )
else
    γc      .=  rand(Float64,ncx,ncy) .* γ0
end

# ================================================================== #
if do_save
    anim = Animation()
    anim2 = Animation()
    anim3 = Animation()
    outdir  =   "$(Ddir)/$(Type)/$(Def)_gam0_$(γ0)_Dmax_$(Dmax)_etar_$(η_reg)_ncx_$(ncx)_gaminc_$(γinc)"
    chdir   =   isdir(outdir)
    if !chdir
        mkdir(outdir)
    end
end
# Start time loop ====================================================
for it = 1:nt
    # Time step length -----------------------------------------------
    Vxc     .= (Vx[1:end-1,2:end-1].+Vx[2:end,2:end-1])./2.0
    Vyc     .= (Vy[2:end-1,1:end-1].+Vy[2:end-1,2:end])./2.0
    Δt      =   0.5 * min(Δx,Δy) / max(maximum(abs.(Vxc)),maximum(abs.(Vyc)));
    # Elastic viscosity ----------------------------------------------
    ηel     =   G*Δt
    # Time -----------------------------------------------------------
    t           += Δt
    time[it]    = t
    # Damage update --------------------------------------------------
    Dc      .=  Dmax .*γc ./ γcr        
    # Viscosity initialization ---------------------------------------
    if !PSS
        av(Dv[2:end-1,2:end-1],Dc)
        Dv[1,:] = Dv[2,:]; Dv[end,:] = Dv[end-1,:]; Dv[:,1] = Dv[:,2]; Dv[:,end] = Dv[:,end-1]
        ηvi                 .=  ηm  .* (1 .- Dv)
        ηvi[Elli .< ri.^2]  .=  ηi .* (1 .- Dv[Elli .< ri.^2])
    end
    ηv      .=  (1.0 ./ ηvi .+ 1.0 ./ ηel ).^(-1.0)
    av(ηc,ηv)
    av(ηv[2:end-1,2:end-1],ηc)
    @printf("Step %05d --- Δηve: %2.2e, min(ηve): %2.2e, max(ηve): %2.2e\n", 
            it, log10(maximum(ηc)/minimum(ηc)),minimum(ηc), maximum(ηc))
    λp      .= 0.0    
    # Update stress field --------------------------------------------
    τxxo    .=  τxx
    τyyo    .=  τyy
    τxyo    .=  τxy
    τxyov   .=  τxyv
    # Stokes residual evaluation =====================================
    # Iterative parameters -------------------------------------------
    Reopt   =   5*pi
    cfl     =   0.5
    ρ       =   cfl*Reopt/ncx
    @time @views for iter = 1:50000    
        # Viscosity update -------------------------------------------           
        @tturbo @. ηv   =  (1.0 / ηvi + 1.0 / ηel )^(-1.0) 
        @tturbo @. ηc   = (ηv[1:end-1,1:end-1] + ηv[2:end,1:end-1] + 
                    ηv[1:end-1,2:end] + ηv[2:end,2:end]) / 4.0    
        # Boundary velocity ------------------------------------------
        if pureshear
            @tturbo @. Vx[:,1]   = Vx[:,2]                  # bottom
            @tturbo @. Vx[:,end] = Vx[:,end-1]              # top
            @tturbo @. Vy[1,:]      = Vy[2,:]               # left
            @tturbo @. Vy[end,:]    = Vy[end-1,:]           # right
        else
            @tturbo @. Vx[:,1]   = 2 * VBCS - Vx[:,2]       # bottom
            @tturbo @. Vx[:,end] = 2 * VBCN - Vx[:,end-1]   # top
            @tturbo @. Vy[1,:]      = Vy[end-1,:]           # left
            @tturbo @. Vy[end,:]    = Vy[2,:]               # right
        end        
        # Calculate velocity gradient tensor ------------------------- 
        @tturbo @. dvxdx = (Vx[2:end,2:end-1]-Vx[1:end-1,2:end-1])/Δx
        @tturbo @. dvydy = (Vy[2:end-1,2:end]-Vy[2:end-1,1:end-1])/Δy
        @tturbo @. dvxdy = (Vx[:,2:end]-Vx[:,1:end-1])/Δy
        @tturbo @. dvydx = (Vy[2:end,:]-Vy[1:end-1,:])/Δx
        # Calculate divergence of velocity ---------------------------
        @tturbo @. ∇V   =  dvxdx + dvydy
        # Calculate strain rate tensor -------------------------------
        @tturbo @. εxx  =  dvxdx - 1/3 *∇V
        @tturbo @. εyy  =  dvydy - 1/3 *∇V
        @tturbo @. εxyv =  1/2 *( dvxdy + dvydx )
        @tturbo @. εxy  = (εxyv[1:end-1,1:end-1] + εxyv[2:end,1:end-1] + 
                    εxyv[1:end-1,2:end] + εxyv[2:end,2:end]) / 4.0                    
        @tturbo @. εII = sqrt(0.5*(εxx^2 + εyy^2) + εxy^2)        
        # Visco-elastic strain rates ---------------------------------
        @tturbo @. εxxeff  = εxx + (τxxo[2:end-1,:] / 2.0 / ηel)
        @tturbo @. εyyeff  = εyy + (τyyo / 2.0 / ηel)
        @tturbo @. εxyeffv = εxyv + (τxyov / 2.0 / ηel)
        @tturbo @. εxyeff = (εxyeffv[1:end-1,1:end-1] + εxyeffv[2:end,1:end-1] + 
                    εxyeffv[1:end-1,2:end] + εxyeffv[2:end,2:end]) / 4.0    
        @tturbo @. εIIeff = sqrt(0.5*(εxxeff^2 + εyyeff^2) + εxyeff^2)
        ## Trial stress -----------------------------------------------
        @tturbo @. τxx[2:end-1,:]  =  2.0*ηc*εxxeff
        @tturbo @. τyy             =  2.0*ηc*εyyeff
        @tturbo @. τxy             =  2.0*ηc*εxyeff        
        @tturbo @. τII = sqrt(0.5*(τxx[2:end-1,:]^2 + τyy^2) + τxy^2 )  
        # Plasticity -------------------------------------------------
        if PSS
            @tturbo @. τyield   =   τyield0 * (1 - Dc)
        else
            @tturbo @. τyield   =   τyield0
        end
        @tturbo @. F        =   τII - τyield
        @tturbo @. ispl     =   0.0
        @tturbo @. ispl     =   F >= 0.0
        @tturbo @. λp       =   F*ispl / ( ηc + η_reg )
        @tturbo @. dQdτxx   =   0.5*τxx[2:end-1,:]/τII
        @tturbo @. dQdτyy   =   0.5*τyy/τII
        @tturbo @. dQdτxy   =       τxy/τII
        # Plastic correction 
        @tturbo @. τxx[2:end-1,:] =   
                                2.0 * ηc*(εxxeff -     λp*dQdτxx)
        @tturbo @. τyy      =   2.0 * ηc*(εyyeff -     λp*dQdτyy)
        @tturbo @. τxy      =   2.0 * ηc*(εxyeff - 0.5*λp*dQdτxy)
        @tturbo @. τII      = 	sqrt(0.5*(τxx[2:end-1,:]^2 + τyy^2) + τxy^2)
        @tturbo @. Fchk     =   τII - τyield - λp * η_reg
        # Effective viscosity ----------------------------------------
        @tturbo @. ηc = τII / 2.0 / εIIeff
        @tturbo @. ηv[2:end-1,2:end-1] = 
                    (ηc[1:end-1,1:end-1] + ηc[2:end,1:end-1] + 
                     ηc[1:end-1,2:end]   + ηc[2:end,2:end]) ./ 4.0
        if pureshear
            @tturbo @. ηv[1,:]   = ηv[2,:]
            @tturbo @. ηv[end,:] = ηv[end-1,:]
        else
            @tturbo @. ηBC = (ηc[1,1:end-1] + ηc[end,1:end-1] + ηc[1,2:end]   + ηc[end,2:end]) ./ 4.0
            @tturbo @. ηv[1,2:end-1] = ηBC
            @tturbo @. ηv[end,2:end-1] = ηBC
        end
        @tturbo @. ηv[:,1]   = ηv[:,2]
        @tturbo @. ηv[:,end] = ηv[:,end-1]
        @tturbo @. τxyv = 2.0 * ηv * εxyeffv
        # PT time step -----------------------------------------------
        @tturbo @. Δτv  = ρ*min(Δx,Δy)^2 / ηv / 4.1 * cfl
        @tturbo @. Δτvx = (Δτv[:,1:end-1] + Δτv[:,2:end]) / 2
        @tturbo @. Δτvy = (Δτv[1:end-1,:] + Δτv[2:end,:]) / 2
        @tturbo @. κΔτp = cfl * ηc * Δx / (xmax-xmin)
        # Define residuals ===========================================
        # Continuity equation ----------------------------------------
        @tturbo @. Fp    =  -∇V
        if pureshear==false Fp  .-= mean(Fp) end
        # x - stokes equation ----------------------------------------
        if pureshear # pure shear boundary conditions
            @tturbo @. Fx[2:end-1,:]    =  
                -(P[3:end-1,:] - P[2:end-2,:])/Δx + 
                (τxx[3:end-1,:] - τxx[2:end-2,:])/Δx + 
                (τxyv[2:end-1,2:end] - τxyv[2:end-1,1:end-1])/Δy
        else # simple shear boundary conditions; periodic
            @tturbo @. P[1,:]       = P[end-1,:]
            @tturbo @. P[end,:]     = P[2,:]
            @tturbo @. τxx[1,:]     = τxx[end-1,:]
            @tturbo @. τxx[end,:]   = τxx[2,:]
            @tturbo @. Fx   =  
                -(P[2:end,:] - P[1:end-1,:])/Δx + 
                (τxx[2:end,:] - τxx[1:end-1,:])/Δx + 
                (τxyv[:,2:end] - τxyv[:,1:end-1])/Δy
        end
        # y - stokes equation ----------------------------------------
        @tturbo @. Fy[:,2:end-1]    =  
                    -(P[2:end-1,2:end] - P[2:end-1,1:end-1])/Δy + 
                    (τyy[:,2:end] - τyy[:,1:end-1])/Δy + 
                    (τxyv[2:end,2:end-1] - τxyv[1:end-1,2:end-1])/Δx
        # Calculate rate update --------------------------------------
        @tturbo @. dvxdτ   =  (1-ρ) * dvxdτ + Fx
        @tturbo @. dvydτ   =  (1-ρ) * dvydτ + Fy
        # Update velocity and pressure -------------------------------
        @tturbo @. Vx[:,2:end-1]   = Vx[:,2:end-1] + Δτvx / ρ * dvxdτ
        @tturbo @. Vy[2:end-1,:]   = Vy[2:end-1,:] + Δτvy / ρ * dvydτ
        @tturbo @. P[2:end-1,:]    = P[2:end-1,:] + κΔτp * Fp
        # Check converge ---------------------------------------------
        if iter%2000==0 || iter==1
            @printf("   Iter %05d --- Fx: %2.2e Fy: %2.2e Fp: %2.2e\n", 
            iter, norm(Fx)/length(Fx), norm(Fy)/length(Fy), norm(Fp)/length(Fp))
            converge = 
                norm(Fx)/length(Fx)<1e-11 && 
                norm(Fy)/length(Fy)<1e-11 && 
                norm(Fp)/length(Fp)<1e-11
            diverge  = 
                norm(Fx)/length(Fx)>1e2  && 
                norm(Fy)/length(Fy)>1e2  && 
                norm(Fp)/length(Fp)>1e2
            if converge success=true; break end
            if diverge success=false; break end
            if isnan(norm(Fx)) success=false; break end
        end
    end # End Rheological iterations
    @printf("Step %05d --- Δηvep: %2.2e, min(ηvep): %2.2e, max(ηvep): %2.2e\n\n", 
        it, log10(maximum(ηc)/minimum(ηc)),minimum(ηc), maximum(ηc))
    if success==false break end
    # Calculate properties on the vertices ===========================   
    av(εxxv[2:end-1,2:end-1],εxx)
    εxxv[1,:] = εxxv[2,:]; εxxv[end,:] = εxxv[end-1,:]; εxxv[:,1] = εxxv[:,2]; εxxv[:,end] = εxxv[:,end-1]
    av(εyyv[2:end-1,2:end-1],εyy)
    εyyv[1,:] = εyyv[2,:]; εyyv[end,:] = εyyv[end-1,:]; εyyv[:,1] = εyyv[:,2]; εyyv[:,end] = εyyv[:,end-1]
    @. εIIv = sqrt(0.5*(εxxv^2 + εyyv^2) + εxyv^2)
    @. τxxv[2:end-1,2:end-1] = (τxx[2:end-2,1:end-1] + τxx[3:end-1,1:end-1] + τxx[2:end-2,2:end] + τxx[3:end-1,2:end]) / 4.0
    @. τxxv[1,:] = τxxv[2,:]; τxxv[end,:] = τxxv[end-1,:]; τxxv[:,1] = τxxv[:,2]; τxxv[:,end] = τxxv[:,end-1]
    av(τyyv[2:end-1,2:end-1],τyy)
    @. τyyv[1,:] = τyyv[2,:]; τyyv[end,:] = τyyv[end-1,:]; τyyv[:,1] = τyyv[:,2]; τyyv[:,end] = τyyv[:,end-1]
    @. τxxov[2:end-1,2:end-1] = (τxxo[2:end-2,1:end-1] + τxxo[3:end-1,1:end-1] + τxxo[2:end-2,2:end] + τxxo[3:end-1,2:end]) / 4.0
    @. τxxov[1,:] = τxxov[2,:]; τxxov[end,:] = τxxov[end-1,:]; τxxov[:,1] = τxxov[:,2]; τxxov[:,end] = τxxov[:,end-1]
    av(τyyov[2:end-1,2:end-1],τyyo)
    @. τyyov[1,:] = τyyov[2,:]; τyyov[end,:] = τyyov[end-1,:]; τyyov[:,1] = τyyov[:,2]; τyyov[:,end] = τyyov[:,end-1]
    # Viscous ========================================================
    @. εxxvi = τxxv / (2.0 * ηvi)
    @. εyyvi = τyyv / (2.0 * ηvi)
    @. εxyvi = τxyv / (2.0 * ηvi)
    @. εIIvi = sqrt( 0.5*(εxxvi^2 + εyyvi^2) + εxyvi^2)  
    av(εIIvic,εIIvi)
    # Plasticity =====================================================
    @. εxxpl = λp*dQdτxx
    av(εxxplv[2:end-1,2:end-1],εxxpl)
    @. εxxplv[1,:] = εxxplv[2,:]; εxxplv[end,:] = εxxplv[end-1,:]; εxxplv[:,1] = εxxplv[:,2]; εxxplv[:,end] = εxxplv[:,end-1]
    @. εyypl = λp*dQdτyy
    av(εyyplv[2:end-1,2:end-1],εyypl)
    @. εyyplv[1,:] = εyyplv[2,:]; εyyplv[end,:] = εyyplv[end-1,:]; εyyplv[:,1] = εyyplv[:,2]; εyyplv[:,end] = εyyplv[:,end-1]
    @. εxypl = λp*dQdτxy
    av(εxyplv[2:end-1,2:end-1],εxypl)
    @. εxyplv[1,:] = εxyplv[2,:]; εxyplv[end,:] = εxyplv[end-1,:]; εxyplv[:,1] = εxyplv[:,2]; εxyplv[:,end] = εxyplv[:,end-1]
    @. εIIpl[ispl==1] = λp[ispl==1] / 2.0  
    av(εIIplv[2:end-1,2:end-1],εIIpl)
    @. εIIplv[1,:] = εIIplv[2,:]; εIIplv[end,:] = εIIplv[end-1,:]; εIIplv[:,1] = εIIplv[:,2]; εIIplv[:,end] = εIIplv[:,end-1]
    # Elasticity =====================================================
    @. εxxel = (τxxv - τxxov) / (2.0 * ηel)
    @. εyyel = (τyyv - τyyov) / (2.0 * ηel)
    @. εxyel = (τxyv - τxyov) / (2.0 * ηel)
    @. εIIel = sqrt( 0.5*(εxxel^2 + εyyel^2) + εxyel^2)
    # Effective Rheology =============================================
    τIIeff  .=  2.0 .* ηc .* εIIeff
    ## Check Strain Rate Calculations ================================    
    # εcheck  .=  εIIv .- εIIvi .- εIIplv .- εIIel
    # εcheck  .=  εxxv .- εxxvi .- εxxel .- εxxplv
    εcheck  .=  εIIplv ./ ( εIIel + εIIvi )
    # Calculate von Mises strain ===================================== 
    # γc          .+= Δt .* εII    
    if PSS
        γc       .= Δt .* εIIpl .+ γc .* ( 1 - H * Δt )
    else
        γc       .= Δt .* εIIvic .+ γc .* ( 1 - H * Δt )
    end
    γc           .=  min.(γc,γcr)    
    # Calculate time parameters ======================================
    τmean[it]   =   norm(τII)/length(τII)
    τmax[it]    =   maximum(τII)
    γmean[it]   =   norm(γc)/length(γc)
    γmax[it]    =   maximum(γc)
    # Total velocity =================================================
    Vxc     .= (Vx[1:end-1,2:end-1].+Vx[2:end,2:end-1])./2.0
    Vyc     .= (Vy[2:end-1,1:end-1].+Vy[2:end-1,2:end])./2.0
    Vc      .= sqrt.(Vxc.^2 + Vyc.^2)    
    if do_save
        # Create Animation ============================================
        p1 = heatmap(xc,yc,Vc',title="V")
        p2 = heatmap(xc,yc,P[2:end-1,:]',title="P")    
        p3 = heatmap(xc,yc,γc',title="γc")
        p4 = heatmap(xc,yc,log10.(ηc'),title="log10(ηc)")
        plot(p1,p2,p3,p4,aspect_ratio=1,xlim=(xmin,xmax),
            ylim=(ymin,ymax),xlabel="x",ylabel="y"); frame(anim)
        p5 = heatmap(xv,yv,εIIv',title="εII")
        p6 = heatmap(xv,yv,εIIvi',title="εIIvi")
        p7 = heatmap(xv,yv,εIIel',title="εIIel")
        p8 = heatmap(xv,yv,εIIplv',title="εIIpl")
        plot(p5,p6,p7,p8,aspect_ratio=1,xlim=(xmin,xmax),
            ylim=(ymin,ymax),xlabel="x",ylabel="y"); frame(anim2)
        p9 = heatmap(xc,yc,τyield',title="τyield")
        p10 = heatmap(xc,yc,τII',title="τII")
        p11 = heatmap(xv,yv,log10.(εcheck'),title="log( εIIpl / (εIIel + εIIvi) )")
        plot(p9,p10,p11,aspect_ratio=1,xlim=(xmin,xmax),
            ylim=(ymin,ymax),xlabel="x",ylabel="y"); frame(anim3)    
        file = matopen(string("$(outdir)/Data", lpad(it,4,"0"),   ".mat"), "w");
        write(file, "Pt",   Array(P));
        write(file, "Vx",   Array(Vxc));
        write(file, "Vy",   Array(Vyc)); 
        write(file, "Vc",   Array(Vc)); 
        write(file, "Tii",  Array(τII));
        write(file, "Eiiv", Array(εIIv));
        write(file, "Eiivi", Array(εIIvi));
        write(file, "Eiiel", Array(εIIel));
        write(file, "Eiipl", Array(εIIpl));
        write(file, "Gamma", Array(γc));
        write(file, "eta", Array(ηc));
        close(file)    
    end
end     # end time loop
# Plot ===============================================================
if success==true
    if do_save    
        gif(anim, "$(outdir)/Kinamatics_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).gif", fps = 15)
        gif(anim2, "$(outdir)/StrainRates_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).gif", fps = 15)
        gif(anim3, "$(outdir)/Stress_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).gif", fps = 15)
    end
    p1 = heatmap(xc,yc,Vc',title="V")
    p2 = heatmap(xc,yc,P[2:end-1,:]',title="P")
    p3 = heatmap(xc,yc,γc',title="γc")
    p4 = heatmap(xc,yc,log10.(ηc'),title="log10(ηc)")

    p5 = heatmap(xv,yv,εIIv',title="εII")
    p6 = heatmap(xv,yv,εIIvi',title="εIIvi")
    p7 = heatmap(xv,yv,εIIel',title="εIIel")
    p8 = heatmap(xv,yv,εIIplv',title="εIIpl")

    p9 = heatmap(xc,yc,τyield',title="τyield")
    p10 = heatmap(xc,yc,τII',title="τII")
    #p11 = heatmap(xv,yv,εcheck',title="εcheck")
    p11 = heatmap(xv,yv,log10.(εcheck'),title="log( εIIpl / (εIIel + εIIvi) )")

    display(plot(p1,p2,p3,p4,aspect_ratio=1,xlim=(xmin,xmax),
        ylim=(ymin,ymax),xlabel="x",ylabel="y"))
    if do_save        
        savefig("$(outdir)/Kinamatics_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).png")
    end
    display(plot(p5,p6,p7,p8,aspect_ratio=1,xlim=(xmin,xmax),
        ylim=(ymin,ymax),xlabel="x",ylabel="y"))
    if do_save    
        savefig("$(outdir)/StrainRates_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).png")
    end
        display(plot(p9,p10,p11,aspect_ratio=1,xlim=(xmin,xmax),
        ylim=(ymin,ymax),xlabel="x",ylabel="y"))
    if do_save    
        savefig("$(outdir)/Stress_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).png")
    end

    p11 = plot(time,τmean,xlabel="time",ylabel="τmean")
    p12 = plot(time,τmax,xlabel="time",ylabel="τmax")
    p13 = plot(time,γmean,xlabel="time",ylabel="γmean")
    p14 = plot(time,γmax,xlabel="time",ylabel="γmax")
    display(plot(p11,p12,p13,p14))    
    if do_save    
        savefig("$(outdir)/TimeSeries_$(Def)_gam0_$(γ0)_Dmax_$(Dmax).png")
    end

    @printf("ΔF: %2.2e, min(Fchk): %2.2e, max(Fchk): %2.2e\n", 
        norm(Fchk)/length(Fchk),
        minimum(Fchk), maximum(Fchk))
    @printf("Δε: %2.2e, min(εcheck): %2.2e, max(εcheck): %2.2e\n",
        norm(εcheck)/length(εcheck),
        minimum(εcheck), maximum(εcheck))
end

if do_save
    file = matopen(string("$(outdir)/TimeSeries$(Def)_gam0_$(γ0)_Dmax_$(Dmax)",   ".mat"), "w");
    write(file, "time",    Array(time));
    write(file, "Tmean",    Array(τmean));
    write(file, "Tmax",     Array(τmax));
    write(file, "Gmean",    Array(γmean));
    write(file, "Gmax",     Array(γmax));
    close(file)
end

end     # end main function

# ==================================================================== #
# ==================================================================== #
# ==================================================================== #
N       = 50
η_reg   =   1.2e-2
γ0      = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0] # 0.0, 1.0

for igamma=1:length(γ0)
    main(N,η_reg,γ0[igamma])
end
