using DifferentialEquations, ModelingToolkit, JLD, QuasiMonteCarlo, GlobalSensitivity, Statistics
include("/home/harry/Desktop/PhD/Year 3/Measurement Project/Elements.jl")

begin
    ## Input Parameters ## 
    τ = 0.81
    # Left Ventricle 
    Emax_lv = 2.8
    Emin_lv = 0.07
    V0_lv = 20.0
    τes_lv = 0.269*τ
    τep_lv = 0.452*τ
    # Right Ventricle 
    Emax_rv = 0.45
    Emin_rv = 0.035
    V0_rv = 30.0
    τes_rv = 0.269*τ
    τep_rv = 0.452*τ
    # Left Atrium
    Emax_la = 0.13
    Emin_la = 0.09
    V0_la = 3.0
    τes_la = 0.110*τ
    τep_la = 0.18*τ
    Eshift_la = 0.85
    # Right Atrium
    Emax_ra = 0.09
    Emin_ra = 0.045
    V0_ra = 7.0
    τes_ra = 0.110*τ
    τep_ra = 0.18*τ
    Eshift_ra = 0.85 
    # Valve parameters
    Rav = 0.01
    Rmv = 0.005
    Rpv = 0.01
    Rtv = 0.005
    ## Model Parameters 
    # Systemic arteries
    Rsa = 0.0448
    Csa = 0.983
    # Systemic vascular bed
    Rsvb = 0.824
    # Systemic Veins
    Rsv = 0.0269
    Csv = 29.499
    # Pulmonary Arteries
    Rpa = 0.003
    Cpa = 6.7
    # Pulmonary Vascular bed 
    Rpvb = 0.0552
    # Pulmonary Veins
    RPV = 0.0018
    Cpv = 15.8
    # Inital Conditions 
    vlv_0 = 149.6
    vrv_0 = 189.2
    vla_0 = 71.0
    vra_0 = 67.0
    vsa_0 = 98.3
    vsv_0 = 117.996
    vpa_0 = 100.5
    vpv_0 = 126.4
end 


# Build the Model 
@parameters t

# Ventricle & Atrium for Shi 
@named LV = ShiChamberV(V₀=V0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τₑₛ=τes_lv, τₑₚ=τep_lv)
@named RV = ShiChamberV(V₀=V0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τₑₛ=τes_rv, τₑₚ=τep_rv)
@named LA = ShiChamberA(V₀=V0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τₑₛ=τes_la, τₑₚ=τep_la, Eshift=Eshift_la)
@named RA = ShiChamberA(V₀=V0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τₑₛ=τes_ra, τₑₚ=τep_ra, Eshift=Eshift_ra)
# Heart valves

@named AV = ResistorDiode(R = Rav)
@named PV = ResistorDiode(R = Rpv)
@named MV = ResistorDiode(R = Rmv)
@named TV = ResistorDiode(R = Rtv)

# System Model 

#SA is systemic arterys, Svb is systemic vacular bed, Sv systemic veins
# Systemic Loop
@named SA = CR(C = Csa, R = Rsa)
@named SVB = Resistor(R = Rsvb)
@named SV = CR(C = Csv, R = Rsv)
# pulmonary loop 
@named PA = CR(C = Cpa, R = Rpa)
@named PVB = Resistor(R = Rpvb)
@named PVe = CR(C = Cpv, R = Rpv)

circ_eqs = [
    #Ventricle out 
    connect(LV.out, AV.in)
    connect(AV.out, SA.in)
    connect(SA.out, SVB.in)
    connect(SVB.out, SV.in)
    connect(SV.out, RA.in)
    connect(RA.out, TV.in)
    connect(TV.out, RV.in)
    connect(RV.out, PV.in)
    connect(PV.out, PA.in)
    connect(PA.out, PVB.in)
    connect(PVB.out, PVe.in)
    connect(PVe.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODAE system
@named _circ_model = ODESystem(circ_eqs, t)

##
@named circ_model = compose(_circ_model, [LV, LA, RV, RA, AV, MV, PV, TV, SA, SVB, SV, PA, PVB, PVe])
    
## And simplify it
@time circ_sys = structural_simplify(circ_model)

u0 = [vlv_0, vla_0, vrv_0, vra_0, vsa_0, vsv_0, vsa_0, vpv_0, 0.0, 0.0] 

@time prob = ODEProblem(circ_sys,u0, (0.0, 10.0))

x = LinRange(8*τ, 9*τ, 300)
@time sol = solve(prob, Rodas5P(autodiff = false), reltol = 1e-5, abstol = 1e-5, saveat = x)

#plot(sol, idxs = [LV.p, RV.p])
#plot(sol, idxs = [SA.C.p, LV.p])
#plot(sol, idxs = [LA.p, RA.p])
#plot(sol, idxs = [AV.q, MV.q, TV.q, PV.q])


circ_time_post = function (p)
    out_global
end

###### Global Sensitiivty analysis Level 1 ####
# BP measurement Max(Psa)/Min(Psa)

circ_time_idxs_chunked = function (p)

    global parameterset = p

    Np::Int64 = size(p,2)

    println("Running ", Np, " cases.")

    chunksize::Int64 = Np/chunks

    println("Using ", chunks, " chunk(s) of size ", chunksize, ".")


    out = zeros(1,Np)

    for k in 1:chunks
        offset = Int((k - 1) * Np/chunks)
        startindx = Int(offset + 1)
        endindx = Int(k * Np/chunks)

        println("Starting chunk ", k, ", from ", startindx, " to ", endindx, " with offset ", offset, ".")

        pchunk = p[:,startindx:endindx]

        prob_func(prob,i,repeat) = remake(prob; u0 = [vlv_0, vla_0, vrv_0, vra_0, vsa_0, vsv_0, vsa_0, vpv_0, 0.0, 0.0], p=pchunk[:,i])

        ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)

        @time "    ODE Ensemble solved in " sol = solve(ensemble_prob, Rodas5P(autodiff = false), reltol = 1e-5, abstol = 1e-5, EnsembleThreads();saveat = x,trajectories=chunksize)
        @time "    Results processed in" Threads.@threads for i in 1:chunksize
          # println(i+offset)
          ## 0 - Output ##
          # Pressure
          out[1,i+offset] = (maximum(sol[i][SA.C.p])/minimum(sol[i][SA.C.p])) #  bp 

        end
    end
    global out_global = out
    out
end

chunks::Int64 = 10000

samples = 75000
lb = prob.p - 0.15*prob.p
ub = prob.p + 0.15*prob.p
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(samples,lb,ub,sampler)
@time res = gsa(circ_time_idxs_chunked,Sobol(nboot = 1000),A,B,batch = true, Ei_estimator = :Jansen1999)

#save("Level2.jld", "data", res)

res = load("/home/harry/Desktop/PhD/Year 3/Measurement Project/Discrete Levels/Level 1/1_D.jld")["data"]

res.S1_Conf_Int * 1.96
res.ST

using LaTeXStrings, CairoMakie, LinearAlgebra
CairoMakie.activate!(type = "svg")
P = [L"LV_V₀", L"LV_Eₘᵢₙ", L"LV_Eₘₐₓ", L"LV_τₑₛ", L"LV_τₑₚ", L"LA_V₀", L"LA_Eₘᵢₙ", L"LA_Eₘₐₓ", L"LA_τₑₛ", L"LA_τₑₚ", L"LA_Eshift", L"RV_V₀", L"RV_Eₘᵢₙ", L"RV_Eₘₐₓ", L"RV_τₑₛ", L"RV_τₑₚ", L"RA_V₀", L"RA_Eₘᵢₙ", L"RA_Eₘₐₓ", L"RA_τₑₛ", L"RA_τₑₚ", L"RA_Eshift", "AV_R", L"MV_R", L"PV_R", L"TV_R", L"SA_R", L"SA_C", L"SVB_R", L"SV_R", L"SV_C", L"PA_R", L"PA_C", "PVB_R", L"PVe_R", L"PVe_C"]

# Calculate overall parameter influence Ej 
F = transpose(res.ST)*res.ST
e_deomp=eigen(F)

λ = abs.(e_deomp.values)
Q = abs.(e_deomp.vectors)

e_value_sum = sum(λ)

e = Vector{Float64}(undef,36)
for i in 1:36
    for j in 1:36
    e[i] = sum(λ[j]*Q[i,j])/e_value_sum
    end 
end 
p=sortperm(e,rev=true)


f = Figure(size = (1600,1000));
ax = Axis(f[1,1],title="Sobol ST Sensitivity Matrix-Parameter Importance PCA Method", xticks = (1:36, P[p]), xlabel = "Parameters", ylabel = "Importance")
CairoMakie.scatter!( e[p])
hlines!(0.01)
f

show(e[p])
show(p)
#Parameters
@show P[28],P[29],P[4] 

### Examine Parameter Sloppiness ###
using LinearAlgebra, Plots
ST = map(x -> isinf(x) ? zero(x) : x,res.ST)
ST = map(x -> isnan(x) ? zero(x) : x, ST)
F = transpose(ST)*ST

e_deomp=eigen(F)

λ = abs.(e_deomp.values)#/maximum(abs.(e_deomp.values)) .+1e-16

scatter(λ, scale = :log10)

x = ones(36)
y = λ
scatter(x,y,yscale = :log10)
save("/home/harry/Desktop/PhD/Year 3/Measurement Project/Sloppiness Plot/Dis_levels/1.jld","data",λ)