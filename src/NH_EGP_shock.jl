#Endogenous grid point method Non homothetic preferences
#EGP

#BEFORE SHOCK
#Equm iter  58 , r =  0.04179245507908175 , KL ratio new:  5.637648612596845 , KL ratio: 5.63220060883385  KL diff: 0.09672957590411979%
#181.810161 seconds (742.60 M allocations: 81.140 GiB, 7.38% gc time, 0.51% compilation time)
#Mean assets: 5.637648612596844
#Total assets: 56376.48612596845
#Total labor: 10029.160753838052
#Mean assets: 5.637648612596844
#Fraction borrowing constrained: 0.04%
#10th Percentile: 1.962969437643963
#50th Percentile: 5.164911119689187
#90th Percentile: 9.96729777752071
#99th Percentile: 14.335911091105434

#POST SHOCK
#Equm iter  91 , r =  0.04153646520021015 , KL ratio new:  6.133861305325042 , KL ratio: 6.127757574335207  KL diff: 0.09960790575984646%
#367.075040 seconds (1.41 G allocations: 135.060 GiB, 9.03% gc time, 0.34% compilation time)
#Mean assets: 6.133861305325042
#Total assets: 61338.61305325042
#Total labor: 10029.160753838052
#Mean assets: 6.133861305325042
###Fraction borrowing constrained: 0.04%
#10th Percentile: 2.1819423467922543
#50th Percentile: 5.604769190346959
#90th Percentile: 10.819832645248837
#99th Percentile: 15.61411847245358

using NLsolve, Plots, Parameters, Distributions, Random, Statistics, StatsPlots, Interpolations
cd("/Users/antoineding/Documents/GitHub/IIRUHH_master/src")
include("discrete_normal.jl")
include("lininterp1.jl")

# PARAMETERS

@with_kw struct Calibration
    #Households
    σ::Float64=0.5                      # elasticity of relative demand with respect to price sigma=0.5 completementary goods
    ζ::Float64=2.0                      # Intertemporal elasticity of substitution
    θ::Float64=1/ζ                      # Inverse of intertemporal elasticity of substitution 
    γ::Array{Float64}=[1/3, 1/3, 1/3]   # intensity in each good
    ϵ::Array{Float64}=[0.8, 1.0, 1.25]  # elasticity of relative demand with respect to income in luxury good sector
    ε::Array{Float64}=[ϵ[1]/ϵ[2], ϵ[2]/ϵ[2], ϵ[3]/ϵ[2] ]
    ρ::Float64 =(σ-1)/σ                 
    β::Float64 = 0.96                  # Discount factor

    #Production
    α::Float64=0.4                      # Capital share
    δ::Float64=0.1                      # Capital depreciation
    Z::Array{Float64}=[0.9, 1.05, 1.20]   # Sector productivity
end
cal = Calibration()

## Intertemporal Utility
@unpack σ, ζ, θ, γ, ϵ, ρ, β, α, δ, Z = cal
if θ==1
    u(c) = log.(c)
else
    u(c) = (c.^(1-θ).-1)./(1-θ)
end 

u1(u) = u.^(-θ)
u1inv(u) = u.^(-1 ./(θ))


##NH preferences Intratemporal utility
function solvingNHEGP(u, cash::Float64, p_c::Vector{Float64}; cal=cal)
    @unpack σ, ζ, γ, ϵ, ρ, β, α, δ, Z = cal
    out=1-sum(γ[i]^(1/σ)*(((p_c[i]/cash)^(-σ)*γ[i])/u^((1-σ)*σ*ϵ[i]))^ρ for i=1:length(p_c))
    return out
end
function NHUtilityEGP(cash, p_c::Vector{Float64}; cal=cal)
    @unpack σ, ζ, γ, ϵ, ρ, β, α, δ, Z = cal

    #Utility level for given C endowment
    res = nlsolve(u->[solvingNHEGP(u[1], cash, p_c)], [1.0], xtol=TolX)
    U=res.zero[1]

    #Wealth and Minimum expenditure for the utility at given endowment level
    E=cash
    ExpNH=sum(γ[i]*U^(ϵ[i]*(1-σ)^2) * (p_c[i]^(1-σ)) for i = 1:length(p_c))^(1/(1-σ))
    
    #Optimal Demand with income effect
    C1=γ[1]*U^(ϵ[1]*(1-σ)^2) * (p_c[1]/E)^(-σ)
    C2=γ[2]*U^(ϵ[2]*(1-σ)^2) * (p_c[2]/E)^(-σ)
    C3=γ[3]*U^(ϵ[3]*(1-σ)^2) * (p_c[3]/E)^(-σ)

    #Expenditure share with income effect
    ω_p=γ[1]*U^(ϵ[1]*(1-σ)^2) * (p_c[1]/E)^(1-σ)
    ω_n=γ[2]*U^(ϵ[2]*(1-σ)^2) * (p_c[2]/E)^(1-σ)
    ω_l=γ[3]*U^(ϵ[3]*(1-σ)^2) * (p_c[3]/E)^(1-σ)

    ω=[ω_p ω_n ω_l]
    return [C1, C2, C3, U, E, ω, ExpNH]#, U, E, ω]

end

## income risk: discretized N(mu,sigma^2)
mu_y = 1
sd_y = 0.2
ny = 10

## asset grids
na = 40
amax = 50
borrow_lim = 0
agrid_par = 0.5 # 1 for linear, 0 for L-shaped

## computation
max_iter = 200
tol_iter = 5.0e-3
Nsim = 10000
Tsim = 500

# computation KL
maxiter_KL = 200
tol_KL = 1.0e-3
step_KL = 0.01
rguess = 1/β-1-0.0001 # a bit lower than inverse of discount rate
KLratioguess = ((rguess + δ)/α.*Z[2])^(1/(α-1))

# OPTIONS
Display = 2
MakePlots = 1

## which function to interpolation 
InterpCon = 0
InterpEMUC = 1

## tolerance for non-linear solver
TolX=1.0e-2

# SET UP GRIDS

## assets
agrid = range(0,1,length=na)
agrid = agrid.^(1 ./ agrid_par)
agrid = borrow_lim .+ (amax.-borrow_lim).*agrid

## income: disretize normal distribution
width = nlsolve(x -> discrete_normal(ny,mu_y,sd_y,x...)[1],[2.0]).zero
temp, ygrid, ydist = discrete_normal(ny,mu_y,sd_y,width...)
ycumdist = cumsum(ydist)

# DRAW RANDOM NUMBERS
Random.seed!(2022)
# Random value for revenue for Nsim people Tsim period
yrand = rand(Nsim,Tsim)
yindsim = zeros(Int,Nsim,Tsim)

for it = 1:Tsim
    #Income realization and ranking in the distribution there is 10 rank corresponding to deciles
    yindsim[yrand[:,it].<=ycumdist[1],it] .= 1
    for iy = 2:ny
        yindsim[(yrand[:,it].>ycumdist[iy-1]) .& (yrand[:,it].<=ycumdist[iy]),it] .= iy;
    end
end

# ITERATE OVER KL RATIO
KLratio = KLratioguess

iterKL = 0
KLdiff = 1

#first guess price, capital and labor
p0=[0.5, 1.0, 2.0]
Kguess=KLratio.*Nsim
Lguess= Nsim

@time while iterKL<=maxiter_KL && abs(KLdiff)>tol_KL
    p=p0
    iterKL = iterKL + 1

    #Updated firms
    r = α.*Z[2].*KLratio^(α-1) - δ
    R = 1+r
    wage = (1-α).*Z[2].*KLratio^α

    ## initialize consumption function in first iteration only
    if iterKL==1
        cash = zeros(na,ny)
        global Uguess = zeros(na,ny)
        global cashU=zeros(na,ny)
        for iy = 1:ny
            cash[:,iy] = R.*agrid .+ wage.*ygrid[iy]
            for ia=1:na

                # Utility and verify if spending is equal to available resource E
                Uguess[ia,iy] = NHUtilityEGP(cash[ia,iy],p)[4]
                cashU[ia,iy] = NHUtilityEGP(cash[ia,iy],p)[5]

            end
        end
        U=Uguess
    end

    ## solve for policy functions with EGP
    
    iter = 0
    Udiff = 1000
    U=copy(Uguess)

    if iterKL>=2
        U=copy(Ulast)
    end

    while iter <= max_iter && Udiff>tol_iter
        iter = iter + 1
        global Ulast = copy(U)
        global Elast = copy(cashU)
        global sav = zeros(na,ny)
        global C1 = zeros(na,ny)
        global C2 = zeros(na,ny)
        global C3 = zeros(na,ny)
        global S1 = zeros(na,ny)
        global S2 = zeros(na,ny)
        global S3 = zeros(na,ny)

        E=zeros(na)
        emuc = u1(Ulast)*ydist
        muc1 = β.*R.*emuc
        con1 = u1inv(muc1)
        
        #solve for Expenditure
        for iy=1:ny
            for ia=1:na
                E[ia]=nlsolve(x -> [NHUtilityEGP(x[1], p)[4].-con1[ia]], [1.0]).zero[1]
            end
        end

        ## loop over income
        ass1 = zeros(na,ny)
        for iy = 1:ny
            ## loop over current period ssets
            for ia  = 1:na
                ass1[:,iy] = (E.+ agrid.-wage.*ygrid[iy])./R
                if agrid[ia]<ass1[1,iy] # borrowing constraint binds
                    sav[ia,iy] = borrow_lim                
                else # borrowing constraint does not bind;
                    sav[ia,iy] = lininterp1(ass1[:,iy], agrid, agrid[ia])
                end          
                
                # Utility and verify if spending is equal to available resource E
                U[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy], p)[4]
                cashU[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[5]

                # Quantity primary, normal, luxury goods
                C1[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[1]
                C2[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[2]
                C3[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[3]

                #Share in primary, normal, luxury goods
                S1[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[6][1]
                S2[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[6][2]
                S3[ia,iy] = NHUtilityEGP(wage.*ygrid[iy].+ R.*agrid[ia].- sav[ia,iy],p)[6][3]
            end
        end
     
        Udiff = maximum(abs.(U - Ulast))
        if Display>=2
            println("Iteration no. $iter, max con fn diff is  $Udiff")
        end
    end


    ## simulate: start at assets from last interation
    if iterKL>=1
        global asim = zeros(Nsim,Tsim)
    elseif iterKL>1
        # asim[:,1] = Ea.*ones(Nsim,1)
        asim[:,1] = asim[:,Tsim]
    end
    
    ## create interpolating function
    savinterp = Array{Any}(undef,ny)
    for iy = 1:ny
        savinterp[iy] = interpolate((agrid,), sav[:,iy], Gridded(Linear()))
        savinterp[iy] = extrapolate(savinterp[iy], Flat())   # gives edge
    end
    

    ## loop over time periods
    for it = 1:Tsim
        if it < Tsim
            for iy = 1:ny
                asim[yindsim[:,it].==iy,it+1] = savinterp[iy](asim[yindsim[:,it].==iy,it])
            end
        end
    end

    ysim=ygrid[yindsim]

    ## assign actual labor income values
    labincsim = wage.*ysim

    ## mean assets and efficiency units
    global Ea = mean(asim[:,Tsim])
    L = mean(ysim[:,Tsim])
    
    KLrationew = Ea./ L
    
    KLdiff = KLrationew./KLratio - 1
    if Display>=1
        println("Equm iter  $iterKL , r =  $r , KL ratio new:  $KLrationew , KL ratio: $KLratio  KL diff: " * string(KLdiff*100) * "%")
    end

    KLratio = (1-step_KL)*KLratio + step_KL*KLrationew

end

ysim=ygrid[yindsim]
# MAKE PLOTS

if MakePlots==1
## consumption policy function
    p1 = plot(agrid, [Ulast[:,1] Ulast[:,ny]], xlims=(0,amax), title="Consumption Policy Function", color=[:blue :red], label=["Lowest income state" "Highest income state"])
    display(p1)
    savefig(p1, "Figures/Utilityobtained_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")

    p2 = plot(agrid, [S1[:,1] S1[:,ny] S2[:,1] S2[:,ny] S3[:,1] S3[:,ny]], 
    xlims=(0,amax), title="Expenditure shares", color=[:blue :red], 
    label=["Lowest income state Primary" "Highest income state primary" "Lowest income state Normal" "Highest income state Normal" "Lowest income state Luxury" "Highest income state Luxury"])
    display(p2)
    savefig(p2, "Figures/Share_cons_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")


    p3 = plot(agrid, [C1[:,ny] C2[:,ny] C3[:,ny]], 
    xlims=(0,amax), title="Consumption share shares", color=[:blue :red :black], 
    label=["Highest income state primary" "Highest income state Normal"  "Highest income state Luxury"])
    display(p3)
    savefig(p3, "Figures/cons_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")

    ## savings policy function
    p4 = plot(agrid, [sav[:,1].-agrid[:,1] sav[:,ny].-agrid[:,1]], xlims=(0,amax), title="Savings Policy Function", color=[:blue :red], legend=false)
    plot!(agrid, zeros(na,1), color=:black, lw=0.5)
    display(p4)
    savefig(p4, "Figures/saving_decision_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")


    ## income distribution
    p5 = histogram(ysim[:,Tsim], bins=[2*ygrid[1]-ygrid[2];ygrid].+(ygrid[2]-ygrid[1])/2, title="Income distribution", color=RGB(0,0.5,0.5), linecolor=:blue, legend=false)
    display(p5)
    savefig(p5, "Figures/income_dist_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")

    ## asset distribution
    p6 = histogram(asim[:,Tsim], nbins=100, title="Asset distribution", color=RGB(.7,.7,.7), linecolor=:black, legend=false)
    display(p6)
    savefig(p6, "Figures/Asset_dist_EGP_shock.png") # save the most recent fig as filename_string (such as "output.png")

    ## convergence check
    p7 = plot(1:Tsim, mean(asim,dims=1)', ylims=(0,2*Ea), title="Mean Asset Convergence", xlabel="Time Period", color=:black, lw=1.5, legend=false)
    display(p7)

    ## asset distribution statistics
    aysim = asim[:,Tsim]./mean(ysim[:,Tsim])
    println("Mean assets: " * string(mean(aysim)))
    println("Total assets: " * string(cumsum(aysim)[end]))
    println("Total labor: " * string(cumsum(ysim[:,end])[end]))
    println("Mean assets: " * string(mean(aysim)))

    println("Fraction borrowing constrained: " * string(sum(aysim.==borrow_lim)/Nsim * 100) * '%')
    println("10th Percentile: " * string(quantile(aysim,.1)))
    println("50th Percentile: " * string(quantile(aysim,.5)))
    println("90th Percentile: " * string(quantile(aysim,.9)))
    println("99th Percentile: " * string(quantile(aysim,.99)))
end