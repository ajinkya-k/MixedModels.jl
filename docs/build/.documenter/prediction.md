
# Prediction and simulation in Mixed-Effects Models {#Prediction-and-simulation-in-Mixed-Effects-Models}

We recommend the [MixedModelsSim.jl](https://github.com/RePsychLing/MixedModelsSim.jl/) package and associated documentation for useful tools in constructing designs to simulate. For now, we&#39;ll use the sleep study data as a starting point.

```julia
using DataFrames
using MixedModels
using StatsBase
# use a DataFrame to make it easier to change things later
slp = DataFrame(MixedModels.dataset(:sleepstudy))
slpm = fit(MixedModel, @formula(reaction ~ 1 + days + (1|subj)), slp)
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -897.0393  1794.0786  1802.0786  1802.3072  1814.8505

Variance components:
            Column    Variance Std.Dev.
subj     (Intercept)  1296.8692 36.0121
Residual               954.5279 30.8954
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
──────────────────────────────────────────────────
                Coef.  Std. Error      z  Pr(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405     9.50618   26.45    <1e-99
days          10.4673    0.801735  13.06    <1e-38
──────────────────────────────────────────────────
```


## Prediction

The simplest form of prediction are the fitted values from the model: they are indeed the model&#39;s predictions for the observed data.

```julia
predict(slpm) ≈ fitted(slpm)
```


```
true
```


When generalizing to new data, we need to consider what happens if there are new, previously unobserved levels of the grouping variable(s). MixedModels.jl provides three options:
1. `:error`: error on encountering unobserved levels
  
2. `:population`: use population values (i.e. only the fixed effects) for observations with unobserved levels
  
3. `:missing`: return `missing` for observations with unobserved levels.
  

Providing either no prediction (`:error`, `:missing`) or providing the population-level values seem to be the most reasonable ways for _predicting_ new values. For _simulating_ new values based on previous estimates of the variance components, use `simulate`.

In the case where there are no new levels of the grouping variable, all three of these methods provide the same results:

```julia
predict(slpm, slp; new_re_levels=:population) ≈ fitted(slpm)
```


```
true
```


```julia
predict(slpm, slp; new_re_levels=:missing) ≈ fitted(slpm)
```


```
true
```


```julia
predict(slpm, slp; new_re_levels=:error) ≈ fitted(slpm)
```


```
true
```


In the case where there are new levels of the grouping variable, these methods differ.

```julia
# create a new level
slp2 = transform(slp, :subj => ByRow(x -> (x == "S308" ? "NEW" : x)) => :subj)
```


```
180×3 DataFrame
 Row │ subj    days  reaction
     │ String  Int8  Float64
─────┼────────────────────────
   1 │ NEW        0   249.56
   2 │ NEW        1   258.705
   3 │ NEW        2   250.801
   4 │ NEW        3   321.44
   5 │ NEW        4   356.852
   6 │ NEW        5   414.69
   7 │ NEW        6   382.204
   8 │ NEW        7   290.149
  ⋮  │   ⋮      ⋮       ⋮
 174 │ S372       3   310.632
 175 │ S372       4   287.173
 176 │ S372       5   329.608
 177 │ S372       6   334.482
 178 │ S372       7   343.22
 179 │ S372       8   369.142
 180 │ S372       9   364.124
              165 rows omitted
```


```julia
try
  predict(slpm, slp2; new_re_levels=:error)
catch e
  show(e)
end
```


```
ArgumentError("New level encountered in subj")
```


```julia
predict(slpm, slp2; new_re_levels=:missing)
```


```
180-element Vector{Union{Missing, Float64}}:
    missing
    missing
    missing
    missing
    missing
    missing
    missing
    missing
    missing
    missing
   ⋮
 279.92212396847816
 290.38940992807414
 300.8566958876701
 311.32398184726617
 321.79126780686215
 332.25855376645814
 342.7258397260541
 353.1931256856501
 363.66041164524614
```


```julia
predict(slpm, slp2; new_re_levels=:population)
```


```
180-element Vector{Float64}:
 251.40510484848508
 261.8723908080811
 272.3396767676771
 282.80696272727306
 293.2742486868691
 303.7415346464651
 314.20882060606107
 324.67610656565705
 335.14339252525303
 345.6106784848491
   ⋮
 279.92212396847816
 290.38940992807414
 300.8566958876701
 311.32398184726617
 321.79126780686215
 332.25855376645814
 342.7258397260541
 353.1931256856501
 363.66041164524614
```


::: tip Note

Currently, we do not support predicting based on a subset of the random effects.

:::

::: tip Note

`predict` is deterministic (within the constraints of floating point) and never adds noise to the result. If you want to construct prediction intervals, then `simulate` will generate new data with noise (including new values of the random effects).

:::

For generalized linear mixed models, there is an additional keyword argument to `predict`: `type` specifies whether the predictions are returned on the scale of the linear predictor (`:linpred`) or on the level of the response `(:response)` (i.e. the level at which the values were originally observed).

```julia
cbpp = DataFrame(MixedModels.dataset(:cbpp))
cbpp.rate = cbpp.incid ./ cbpp.hsz
gm = fit(MixedModel, @formula(rate ~ 1 + period + (1|herd)), cbpp, Binomial(), wts=float(cbpp.hsz))
predict(gm, cbpp; type=:response) ≈ fitted(gm)
```


```
false
```


```julia
logit(x) = log(x / (1 - x))
predict(gm, cbpp; type=:linpred) ≈ logit.(fitted(gm))
```


```
false
```


## Simulation

In contrast to `predict`, `simulate` and `simulate!` introduce randomness. This randomness occurs both at the level of the observation-level (residual) variance and at the level of the random effects, where new conditional modes are sampled based on the specified covariance parameter (θ; see [Details of the parameter estimation](/optimization#Details-of-the-parameter-estimation)), which defaults to the estimated value of the model. For reproducibility, we specify a pseudorandom generator here; if none is provided, the global PRNG is taken as the default.

The simplest example of `simulate` takes a fitted model and generates a new response vector based on the existing model matrices combined with noise.

```julia
using Random
ynew = simulate(MersenneTwister(42), slpm)
```


```
180-element Vector{Float64}:
 252.6752171388051
 223.3096173669927
 248.6842723713382
 255.64305141427212
 255.08087250538287
 286.58862952243004
 309.5056915487952
 286.88987784213936
 281.3989535401713
 282.98845110151944
   ⋮
 250.81655745160242
 271.2543649019418
 248.48273283810676
 305.9014182133385
 260.66238230893254
 298.11865200662567
 360.1816309547196
 353.9285987403432
 382.89445970316854
```


The simulated response can also be placed in a pre-allocated vector:

```julia
ynew2 = zeros(nrow(slp))
simulate!(MersenneTwister(42), ynew2, slpm)
ynew2 ≈ ynew
```


```
true
```


Or even directly replace the previous response vector in a model, at which point the model must be refit to the new values:

```julia
slpm2 = deepcopy(slpm)
refit!(simulate!(MersenneTwister(42), slpm2))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -895.6915  1791.3830  1799.3830  1799.6116  1812.1549

Variance components:
            Column   Variance Std.Dev.
subj     (Intercept)  912.8712 30.2138
Residual              972.8494 31.1905
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  248.336      8.32982   29.81    <1e-99
days           8.00599    0.809393   9.89    <1e-22
───────────────────────────────────────────────────
```


This inplace simulation actually forms the basis of [`parametricbootstrap`](/bootstrap#MixedModels.parametricbootstrap).

Finally, we can also simulate the response from entirely new data.

```julia
df = DataFrame(days = repeat(1:10, outer=20), subj=repeat(1:20, inner=10))
df[!, :subj] = string.("S", lpad.(df.subj, 2, "0"))
df[!, :reaction] .= 0
df
```


```
200×3 DataFrame
 Row │ days   subj    reaction
     │ Int64  String  Int64
─────┼─────────────────────────
   1 │     1  S01            0
   2 │     2  S01            0
   3 │     3  S01            0
   4 │     4  S01            0
   5 │     5  S01            0
   6 │     6  S01            0
   7 │     7  S01            0
   8 │     8  S01            0
  ⋮  │   ⋮      ⋮        ⋮
 194 │     4  S20            0
 195 │     5  S20            0
 196 │     6  S20            0
 197 │     7  S20            0
 198 │     8  S20            0
 199 │     9  S20            0
 200 │    10  S20            0
               185 rows omitted
```


```julia
ysim = simulate(MersenneTwister(42), slpm, df)
```


```
200-element Vector{Float64}:
 255.00432150851708
 225.63872173670472
 251.0133767410502
 257.9721557839841
 257.40997687509486
 288.91773389214205
 311.83479591850715
 289.21898221185137
 283.72805790988326
 285.31755547123146
   ⋮
 197.4028748844943
 138.22165310892973
 187.19029860225433
 215.87081467989628
 220.4557805243008
 267.6032837398289
 245.3979923049997
 304.0437000517243
 253.65615200259157
```


Note that this is a convenience method for creating a new model and then using the parameters from the old model to call `simulate` on that model. In other words, this method incurs the cost of constructing a new model and then discarding it. If you could re-use that model (e.g., fitting that model as part of a simulation study), it often makes sense to do these steps to perform these steps explicitly and avoid the unnecessary construction and discarding of an intermediate model:

```julia
msim = LinearMixedModel(@formula(reaction ~ 1 + days + (1|subj)), df)
simulate!(MersenneTwister(42), msim; θ=slpm.θ, β=slpm.β, σ=slpm.σ)
response(msim) ≈ ysim
```


```
true
```


```julia
fit!(msim)
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -996.8296  1993.6592  2001.6592  2001.8643  2014.8524

Variance components:
            Column    Variance Std.Dev.
subj     (Intercept)  1292.4331 35.9504
Residual               956.1199 30.9212
 Number of obs: 200; levels of grouping factors: 20

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  250.284      9.32369   26.84    <1e-99
days           8.08794    0.761227  10.62    <1e-25
───────────────────────────────────────────────────
```


For simulating from generalized linear mixed models, there is no `type` option because the observation-level always occurs at the level of the response and not of the linear predictor.

::: warning Warning

Simulating the model response in place may not yield the same result as simulating into a pre-allocated or new vector, depending on choice of pseudorandom number generator. Random number generation in Julia allows optimization based on type, and the internal storage type of the model response (currently a view into a matrix storing the concatenated fixed-effects model matrix and the response) may not match the type of a pre-allocated or new vector. See also [discussion here](https://discourse.julialang.org/t/weird-prng-behavior/63186).

:::

::: tip Note

All the methods that take new data as a table construct an additional `MixedModel` behind the scenes, even when the new data is exactly the same as the data that the model was fitted to. For the simulation methods in particular, these thus form a convenience wrapper for constructing a new model and calling `simulate` without new data on that model with the parameters from the original model.

:::
