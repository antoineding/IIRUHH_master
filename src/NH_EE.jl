#Euler equation Non homothetic preferences
#EE Equation
#Equm iter 11, r = 0.041823840423476505, KL ratio new: 5.62472384000676, KL ratio: 5.630123436718476 KL diff: -0.09590547653894888\%
#139.061631 seconds (377.04 M allocations: 74.893 GiB, 8.09\% gc time, 0.84\% compilation time)
#Mean assets: 5.62472384000676
#Fraction borrowing constrained: 0.036000000000000004\%
#10th Percentile: 1.9191068770885877
#50th Percentile: 5.172179798601983
#90th Percentile: 9.956860008701538
#99th Percentile: 14.510228336287929

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
    β::Float64 = 0.96                   # Discount factor

    #Production
    α::Float64=0.4                      # Capital share
    δ::Float64=0.1                      # Capital depreciation
    Z::Array{Float64}=[0.9, 1.0, 1.2]   # Sector productivity
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
function solvingNHEE(u, cash::Float64, p_c::Vector{Float64}; cal=cal)
    @unpack σ, ζ, γ, ϵ, ρ, β, α, δ, Z = cal
    out=1-sum(γ[i]^(1/σ)*(((p_c[i]/cash)^(-σ)*γ[i])/u^((1-σ)*σ*ϵ[i]))^ρ for i=1:length(p_c))
    return out
end
function NHUtilityEE(cash, p_c::Vector{Float64}; cal=cal)
    @unpack σ, ζ, γ, ϵ, ρ, β, α, δ, Z = cal

    #Utility level for given C endowment
    res = nlsolve(u->[solvingNHEE(u[1], cash, p_c)], [1.0], xtol=TolX)
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
Nsim = 50000
Tsim = 500

# computation KL
maxiter_KL = 70
tol_KL = 1.0e-3
step_KL = 0.01
rguess = 1/β-1-0.0001 # a bit lower than inverse of discount rate
KLratioguess = ((rguess + δ)/α)^(1/(α-1))

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
p=[0.5, 1.0, 2.0]


@time while iterKL<=maxiter_KL && abs(KLdiff)>tol_KL
    iterKL = iterKL + 1

    r = α.*KLratio^(α-1) - δ
    R = 1+r
    wage = (1-α).* KLratio^α
    

    ## initialize consumption function in first iteration only
    if iterKL==1
        cash = zeros(na,ny)
        global Uguess = zeros(na,ny)
        global cashU=zeros(na,ny)
        for iy = 1:ny
            for ia=1:na
                cash = R.*agrid[ia] .+ wage.*ygrid[iy]
                Uguess[ia,iy] = NHUtilityEE(cash,p)[4]
                cashU[ia,iy] = NHUtilityEE(cash,p)[5]
            end
        end
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
        global emuc = u1(Ulast)*ydist

        ## loop over income
        for iy = 1:ny
            ## loop over current period ssets
            for ia  = 1:na
                 cash = R.*agrid[ia] .+ wage.*ygrid[iy]
                    if u1(cash-borrow_lim) >= β.*R.*lininterp1(agrid,emuc,borrow_lim) # check if borrowing constrained
                        sav[ia,iy] = borrow_lim
                    else
                        sav[ia,iy] = nlsolve(x -> [u1(NHUtilityEE(cash.-x[1], p)[4]).-β.*R.*lininterp1(agrid,emuc,x[1])], [cash.-Elast[ia,iy]]).zero[1] # EE equation iteration
                    end           
                # Quantity primary, normal, luxury goods
                C1[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[1]
                C2[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[2]
                C3[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[3]

                # Utility and verify if spending is equal to available resource E
                U[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[4]
                cashU[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[5]

                #Share in primary, normal, luxury goods
                S1[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[6][1]
                S2[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[6][2]
                S3[ia,iy] = NHUtilityEE(cash .- sav[ia,iy],p)[6][3]       
            end
        end
     
        Udiff = maximum(abs.(U - Ulast))
        if Display>=2
            println("Iteration no. " * string(iter), " max con fn diff is " * string(Udiff))
        end
    end

    for ia=1:na
        for iy=1:ny
            if sav[ia,iy]<=0
                sav[ia,iy]=0
            end
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
    end
    

    ## loop over time periods
    for it = 1:Tsim

        ## asset choice
        if it<Tsim
            for iy = 1:ny
                asim[yindsim[:,it].==iy,it+1] = savinterp[iy](asim[yindsim[:,it].==iy,it])
            end
        end
    end

    ysim = ygrid[yindsim]

    ## assign actual labor income values
    labincsim = wage.*ysim

    ## mean assets and efficiency units
    global Ea = mean(asim[:,Tsim])
    L = mean(ysim[:,Tsim])
    
    KLrationew = Ea/ L
    
    KLdiff = KLrationew/KLratio - 1
    if Display>=1
        println("Equm iter " * string(iterKL) * ", r = " * string(r) * ", KL ratio new: " * string(KLrationew) * ", KL ratio: " * string(KLratio) *" KL diff: " * string(KLdiff*100) * "%")
    end

    KLratio = (1-step_KL)*KLratio + step_KL*KLrationew
end

ysim = ygrid[yindsim]

# MAKE PLOTS
if MakePlots==1
    
    ## consumption policy function
    p1 = plot(agrid, [Ulast[:,1] Ulast[:,ny]], xlims=(0,amax), title="Consumption Policy Function", color=[:blue :red], label=["Lowest income state" "Highest income state"])
    display(p1)
    p1C1 = plot(agrid, [S1[:,1] S1[:,ny] S2[:,1] S2[:,ny] S3[:,1] S3[:,ny]], xlims=(0,amax), title="Consumption Policy Function", color=[:blue :red], label=["Lowest income state" "Highest income state"])
    display(p1C1)
    
    ## savings policy function
    p2 = plot(agrid, [sav[:,1].-agrid[:,1] sav[:,ny].-agrid[:,1]], xlims=(0,amax), title="Savings Policy Function", color=[:blue :red], legend=false)
    plot!(agrid, zeros(na,1), color=:black, lw=0.5)
    display(p2)
    
    ## nice zoom
    xlimits = (0,1)
    xlimind = trues(na)
    if minimum(agrid) < xlimits[1]
        xlimind = xlimind .& (agrid.>=maximum(agrid[agrid<xlimits[1]]))
    elseif minimum(agrid) > xlimits[2]
        xlimind .= 0
    end
    if maximum(agrid) > xlimits[2]
        xlimind = xlimind .& (agrid.<=minimum(agrid[agrid.>xlimits[2]]))
    elseif maximum(agrid) < xlimits[1]
        xlimind .= 0
    end

    ## consumption policy function: zoomed in
    p3 = plot(agrid[xlimind], [Ulast[xlimind,1] Ulast[xlimind,ny]], xlims=xlimits, title="Consumption: Zoomed", marker=:circle, color=[:blue :red], linewidth=2, legend=false)
    plot!(show=true)
    display(p3)

    ## savings policy function: zoomed in
    p4 = plot(agrid[xlimind], [sav[xlimind,1].-agrid[xlimind] sav[xlimind,ny].-agrid[xlimind]], xlims=xlimits, title="Savings: Zoomed (a'-a)", marker=:circle, color=[:blue :red], linewidth=2, legend=false)
    plot!(agrid, zeros(na,1), color=:black, lw=0.5)
    display(p4)
    
    ## income distribution
    p5 = histogram(ysim[:,Tsim], bins=[2*ygrid[1]-ygrid[2];ygrid].+(ygrid[2]-ygrid[1])/2, title="Income distribution", color=RGB(0,0.5,0.5), linecolor=:blue, legend=false)
    display(p5)

    ## asset distribution
    p6 = histogram(asim[:,Tsim], nbins=100, title="Asset distribution", color=RGB(.7,.7,.7), linecolor=:black, legend=false)
    display(p6)
    
    ## convergence check
    p7 = plot(1:Tsim, mean(asim,dims=1)', ylims=(0,2*Ea), title="Mean Asset Convergence", xlabel="Time Period", color=:black, lw=1.5, legend=false)
    display(p7)

    ## asset distribution statistics
    aysim = asim[:,Tsim]./mean(ysim[:,Tsim])
    println("Mean assets: " * string(mean(aysim)))
    println("Fraction borrowing constrained: " * string(sum(aysim.==borrow_lim)/Nsim * 100) * '%')
    println("10th Percentile: " * string(quantile(aysim,.1)))
    println("50th Percentile: " * string(quantile(aysim,.5)))
    println("90th Percentile: " * string(quantile(aysim,.9)))
    println("99th Percentile: " * string(quantile(aysim,.99)))
end