
# API

In addition to its own functionality, MixedModels.jl also implements extensive support for the [`StatsAPI.StatisticalModel`](https://github.com/JuliaStats/StatsAPI.jl/blob/main/src/statisticalmodel.jl) and [`StatsAPI.RegressionModel`](https://github.com/JuliaStats/StatsAPI.jl/blob/main/src/regressionmodel.jl) API.

## Types
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.BlockDescription' href='#MixedModels.BlockDescription'><span class="jlbinding">MixedModels.BlockDescription</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BlockDescription
```


Description of blocks of `A` and `L` in a [`LinearMixedModel`](/api#MixedModels.LinearMixedModel)

**Fields**
- `blknms`: Vector{String} of block names
  
- `blkrows`: Vector{Int} of the number of rows in each block
  
- `ALtypes`: Matrix{String} of datatypes for blocks in `A` and `L`.
  

When a block in `L` is the same type as the corresponding block in `A`, it is described with a single name, such as `Dense`.  When the types differ the entry in `ALtypes` is of the form `Diag/Dense`, as determined by a `shorttype` method.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/blockdescription.jl#L1-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.BlockedSparse' href='#MixedModels.BlockedSparse'><span class="jlbinding">MixedModels.BlockedSparse</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
BlockedSparse{Tv,S,P}
```


A `SparseMatrixCSC` whose nonzeros form blocks of rows or columns or both.

**Members**
- `cscmat`: `SparseMatrixCSC{Tv, Int32}` representation for general calculations
  
- `nzasmat`: nonzeros of `cscmat` as a dense matrix
  
- `colblkptr`: pattern of blocks of columns
  

The only time these are created are as products of `ReMat`s.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/arraytypes.jl#L51-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.FeMat' href='#MixedModels.FeMat'><span class="jlbinding">MixedModels.FeMat</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FeMat{T,S}
```


A matrix and a (possibly) weighted copy of itself.

Typically, an `FeMat` represents the fixed-effects model matrix with the response (`y`) concatenated as a final column.

::: tip Note

`FeMat` is not the same as [`FeTerm`](/api#MixedModels.FeTerm).

:::

**Fields**
- `xy`: original matrix, called `xy` b/c in practice this is `hcat(fullrank(X), y)`
  
- `wtxy`: (possibly) weighted copy of `xy` (shares storage with `xy` until weights are applied)
  

Upon construction the `xy` and `wtxy` fields refer to the same matrix


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L89-L105" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.FeTerm' href='#MixedModels.FeTerm'><span class="jlbinding">MixedModels.FeTerm</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
FeTerm{T,S}
```


Term with an explicit, constant matrix representation

Typically, an `FeTerm` represents the model matrix for the fixed effects.

::: tip Note

`FeTerm` is not the same as [`FeMat`](/api#MixedModels.FeMat)!

:::

**Fields**
- `x`: full model matrix
  
- `piv`: pivot `Vector{Int}` for moving linearly dependent columns to the right
  
- `rank`: computational rank of `x`
  
- `cnames`: vector of column names
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L1-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.FeTerm-Tuple{SparseArrays.SparseMatrixCSC, AbstractVector{String}}' href='#MixedModels.FeTerm-Tuple{SparseArrays.SparseMatrixCSC, AbstractVector{String}}'><span class="jlbinding">MixedModels.FeTerm</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
FeTerm(X::SparseMatrixCSC, cnms)
```


Convenience constructor for a sparse [`FeTerm`](/api#MixedModels.FeTerm) assuming full rank, identity pivot and unit weights.

Note: automatic rank deficiency handling may be added to this method in the future, as discussed in the vignette &quot;[Rank deficiency in mixed-effects models](/rankdeficiency#Rank-deficiency-in-mixed-effects-models)&quot; for general `FeTerm`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L46-L53" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.FeTerm-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T' href='#MixedModels.FeTerm-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T'><span class="jlbinding">MixedModels.FeTerm</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
FeTerm(X::AbstractMatrix, cnms)
```


Convenience constructor for [`FeTerm`](/api#MixedModels.FeTerm) that computes the rank and pivot with unit weights.

See the vignette &quot;[Rank deficiency in mixed-effects models](/rankdeficiency#Rank-deficiency-in-mixed-effects-models)&quot; for more information on the computation of the rank and pivot.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L24-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.GaussHermiteNormalized' href='#MixedModels.GaussHermiteNormalized'><span class="jlbinding">MixedModels.GaussHermiteNormalized</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GaussHermiteNormalized{K}
```


A struct with 2 SVector{K,Float64} members
- `z`: abscissae for the K-point Gauss-Hermite quadrature rule on the Z scale
  
- `wt`: Gauss-Hermite weights normalized to sum to unity
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/gausshermite.jl#L33-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.GeneralizedLinearMixedModel' href='#MixedModels.GeneralizedLinearMixedModel'><span class="jlbinding">MixedModels.GeneralizedLinearMixedModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
GeneralizedLinearMixedModel
```


Generalized linear mixed-effects model representation

**Fields**
- `LMM`: a [`LinearMixedModel`](/api#MixedModels.LinearMixedModel) - the local approximation to the GLMM.
  
- `β`: the pivoted and possibly truncated fixed-effects vector
  
- `β₀`: similar to `β`. Used in the PIRLS algorithm if step-halving is needed.
  
- `θ`: covariance parameter vector
  
- `b`: similar to `u`, equivalent to `broadcast!(*, b, LMM.Λ, u)`
  
- `u`: a vector of matrices of random effects
  
- `u₀`: similar to `u`.  Used in the PIRLS algorithm if step-halving is needed.
  
- `resp`: a `GlmResp` object
  
- `η`: the linear predictor
  
- `wt`: vector of prior case weights, a value of `T[]` indicates equal weights.
  

The following fields are used in adaptive Gauss-Hermite quadrature, which applies only to models with a single random-effects term, in which case their lengths are the number of levels in the grouping factor for that term.  Otherwise they are zero-length vectors.
- `devc`: vector of deviance components
  
- `devc0`: vector of deviance components at offset of zero
  
- `sd`: approximate standard deviation of the conditional density
  
- `mult`: multiplier
  

**Properties**

In addition to the fieldnames, the following names are also accessible through the `.` extractor
- `theta`: synonym for `θ`
  
- `beta`: synonym for `β`
  
- `σ` or `sigma`: common scale parameter (value is `NaN` for distributions without a scale parameter)
  
- `lowerbd`: vector of lower bounds on the combined elements of `β` and `θ`
  
- `formula`, `trms`, `A`, `L`, and `optsum`: fields of the `LMM` field
  
- `X`: fixed-effects model matrix
  
- `y`: response vector
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L1-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.Grouping' href='#MixedModels.Grouping'><span class="jlbinding">MixedModels.Grouping</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
struct Grouping <: StatsModels.AbstractContrasts end
```


A placeholder type to indicate that a categorical variable is only used for grouping and not for contrasts.  When creating a `CategoricalTerm`, this skips constructing the contrasts matrix which makes it robust to large numbers of levels, while still holding onto the vector of levels and constructing the level-to-index mapping (`invindex` field of the `ContrastsMatrix`.).

Note that calling `modelcols` on a `CategoricalTerm{Grouping}` is an error.

**Examples**

```julia
julia> schema((; grp = string.(1:100_000)))
# out-of-memory error

julia> schema((; grp = string.(1:100_000)), Dict(:grp => Grouping()))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/grouping.jl#L1-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.LikelihoodRatioTest' href='#MixedModels.LikelihoodRatioTest'><span class="jlbinding">MixedModels.LikelihoodRatioTest</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LikelihoodRatioTest
```


Results of MixedModels.likelihoodratiotest

**Fields**
- `formulas`: Vector of model formulae
  
- `models`: NamedTuple of the `dof` and `deviance` of the models
  
- `tests`: NamedTuple of the sequential `dofdiff`, `deviancediff`,          and resulting `pvalues`
  

**Properties**
- `deviance` : note that this is actually -2 log likelihood for linear models              (i.e. without subtracting the constant for a saturated model)
  
- `pvalues`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/likelihoodratiotest.jl#L1-L17" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.LinearMixedModel' href='#MixedModels.LinearMixedModel'><span class="jlbinding">MixedModels.LinearMixedModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LinearMixedModel(y, Xs, form, wts=[], σ=nothing, amalgamate=true)
```


Private constructor for a LinearMixedModel.

To construct a model, you only need the response (`y`), already assembled model matrices (`Xs`), schematized formula (`form`) and weights (`wts`). Everything else in the structure can be derived from these quantities.

::: tip Note

This method is internal and experimental and so may change or disappear in a future release without being considered a breaking change.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L83-L95" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.LinearMixedModel-2' href='#MixedModels.LinearMixedModel-2'><span class="jlbinding">MixedModels.LinearMixedModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
LinearMixedModel
```


Linear mixed-effects model representation

**Fields**
- `formula`: the formula for the model
  
- `reterms`: a `Vector{AbstractReMat{T}}` of random-effects terms.
  
- `Xymat`: horizontal concatenation of a full-rank fixed-effects model matrix `X` and response `y` as an `FeMat{T}`
  
- `feterm`: the fixed-effects model matrix as an `FeTerm{T}`
  
- `sqrtwts`: vector of square roots of the case weights.  Can be empty.
  
- `parmap` : Vector{NTuple{3,Int}} of (block, row, column) mapping of θ to λ
  
- `dims` : NamedTuple{(:n, :p, :nretrms),NTuple{3,Int}} of dimensions.  `p` is the rank of `X`, which may be smaller than `size(X, 2)`.
  
- `A`: a `Vector{AbstractMatrix}` containing the row-major packed lower triangle of `hcat(Z,X,y)'hcat(Z,X,y)`
  
- `L`: the blocked lower Cholesky factor of `Λ'AΛ+I` in the same Vector representation as `A`
  
- `optsum`: an [`OptSummary`](/api#MixedModels.OptSummary) object
  

**Properties**
- `θ` or `theta`: the covariance parameter vector used to form λ
  
- `β` or `beta`: the fixed-effects coefficient vector
  
- `λ` or `lambda`: a vector of lower triangular matrices repeated on the diagonal blocks of `Λ`
  
- `σ` or `sigma`: current value of the standard deviation of the per-observation noise
  
- `b`: random effects on the original scale, as a vector of matrices
  
- `u`: random effects on the orthogonal scale, as a vector of matrices
  
- `lowerbd`: lower bounds on the elements of θ
  
- `X`: the fixed-effects model matrix
  
- `y`: the response vector
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1-L30" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.LinearMixedModel-Union{Tuple{T}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any, Any}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any, Any, Any}} where T' href='#MixedModels.LinearMixedModel-Union{Tuple{T}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any, Any}, Tuple{AbstractArray, MixedModels.FeTerm{T}, AbstractVector{<:AbstractReMat{T}}, StatsModels.FormulaTerm, Any, Any, Any}} where T'><span class="jlbinding">MixedModels.LinearMixedModel</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
LinearMixedModel(y, feterm, reterms, form, wts=[], σ=nothing; amalgamate=true)
```


Private constructor for a `LinearMixedModel` given already assembled fixed and random effects.

To construct a model, you only need a vector of `FeMat`s (the fixed-effects model matrix and response), a vector of `AbstractReMat` (the random-effects model matrices), the formula and the weights. Everything else in the structure can be derived from these quantities.

::: tip Note

This method is internal and experimental and so may change or disappear in a future release without being considered a breaking change.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L140-L153" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.MixedModel' href='#MixedModels.MixedModel'><span class="jlbinding">MixedModels.MixedModel</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MixedModel
```


Abstract type for mixed models.  MixedModels.jl implements two subtypes: `LinearMixedModel` and `GeneralizedLinearMixedModel`.  See the documentation for each for more details.

This type is primarily used for dispatch in `fit`.  Without a distribution and link function specified, a `LinearMixedModel` will be fit.  When a distribution/link function is provided, a `GeneralizedLinearModel` is fit, unless that distribution is `Normal` and the link is `IdentityLink`, in which case the resulting GLMM would be equivalent to a `LinearMixedModel` anyway and so the simpler, equivalent `LinearMixedModel` will be fit instead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/MixedModels.jl#L166-L179" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.MixedModelBootstrap' href='#MixedModels.MixedModelBootstrap'><span class="jlbinding">MixedModels.MixedModelBootstrap</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MixedModelBootstrap{T<:AbstractFloat} <: MixedModelFitCollection{T}
```


Object returned by `parametericbootstrap` with fields
- `fits`: the parameter estimates from the bootstrap replicates as a vector of named tuples.
  
- `λ`: `Vector{LowerTriangular{T,Matrix{T}}}` containing copies of the λ field from `ReMat` model terms
  
- `inds`: `Vector{Vector{Int}}` containing copies of the `inds` field from `ReMat` model terms
  
- `lowerbd`: `Vector{T}` containing the vector of lower bounds (corresponds to the identically named field of [`OptSummary`](/api#MixedModels.OptSummary))
  
- `fcnames`: NamedTuple whose keys are the grouping factor names and whose values are the column names
  

The schema of `fits` is, by default,

```
Tables.Schema:
 :objective  T
 :σ          T
 :β          NamedTuple{β_names}{NTuple{p,T}}
 :se         StaticArrays.SArray{Tuple{p},T,1,p}
 :θ          StaticArrays.SArray{Tuple{k},T,1,k}
```


where the sizes, `p` and `k`, of the `β` and `θ` elements are determined by the model.

Characteristics of the bootstrap replicates can be extracted as properties.  The `σs` and `σρs` properties unravel the `σ` and `θ` estimates into estimates of the standard deviations and correlations of the random-effects terms.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L8-L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.MixedModelFitCollection' href='#MixedModels.MixedModelFitCollection'><span class="jlbinding">MixedModels.MixedModelFitCollection</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
MixedModelFitCollection{T<:AbstractFloat}
```


Abstract supertype for [`MixedModelBootstrap`](/api#MixedModels.MixedModelBootstrap) and related functionality in other packages.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.MixedModelProfile' href='#MixedModels.MixedModelProfile'><span class="jlbinding">MixedModels.MixedModelProfile</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
 MixedModelProfile{T<:AbstractFloat}
```


Type representing a likelihood profile of a [`LinearMixedModel`](/api#MixedModels.LinearMixedModel), including associated interpolation splines.

The function [`profile`](/api#MixedModels.profile-Tuple{LinearMixedModel}) is used for computing profiles, while [`confint`](/api#StatsAPI.confint-Tuple{MixedModelProfile}) provides a useful method for constructing confidence intervals from a `MixedModelProfile`.

::: tip Note

The exact fields and their representation are considered implementation details and are **not** part of the public API.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/profile.jl#L1-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.OptSummary' href='#MixedModels.OptSummary'><span class="jlbinding">MixedModels.OptSummary</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
OptSummary
```


Summary of an optimization

**Fields**

**Tolerances, initial and final values**
- `initial`: a copy of the initial parameter values in the optimization
  
- `finitial`: the initial value of the objective
  
- `lowerbd`: lower bounds on the parameter values
  
- `final`: a copy of the final parameter values from the optimization
  
- `fmin`: the final value of the objective
  
- `feval`: the number of function evaluations  Available backends can be examined via `OPTIMIZATION_BACKENDS`.
  
- `returnvalue`: the return value, as a `Symbol`. The available return values will differ between backends.
  
- `xtol_zero_abs`: the tolerance for a near zero parameter to be considered practically zero
  
- `ftol_zero_abs`: the tolerance for change in the objective for setting a near zero parameter to zero
  
- `maxfeval`: as in NLopt (`maxeval`) and PRIMA (`maxfun`)
  

**Choice of optimizer and backend**
- `optimizer`: the name of the optimizer used, as a `Symbol`
  
- `backend`: the optimization library providing the optimizer, default is `NLoptBackend`.
  

**Backend-specific fields**
- `ftol_rel`: as in NLopt, not used in PRIMA
  
- `ftol_abs`: as in NLopt, not used in PRIMA
  
- `xtol_rel`: as in NLopt, not used in PRIMA
  
- `xtol_abs`: as in NLopt, not used in PRIMA
  
- `initial_step`: as in NLopt, not used in PRIMA
  
- `maxtime`: as in NLopt, not used in PRIMA
  
- `rhobeg`: as in PRIMA, not used in NLopt
  
- `rhoend`: as in PRIMA, not used in NLopt
  

**MixedModels-specific fields, unrelated to the optimizer**
- `fitlog`: A vector of tuples of parameter and objectives values from steps in the optimization
  
- `nAGQ`: number of adaptive Gauss-Hermite quadrature points in deviance evaluation for GLMMs
  
- `REML`: use the REML criterion for LMM fits
  
- `sigma`: a priori value for the residual standard deviation for LMM
  

::: tip Note

The internal storage of the parameter values within `fitlog` may change in the future to use a different subtype of `AbstractVector` (e.g., `StaticArrays.SVector`) for each snapshot without being considered a breaking change.

:::

::: tip Note

The exact order and number of fields may change as support for additional backends and features thereof are added. In other words: use the keyword constructor and do **not** use the positional constructor.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/optsummary.jl#L1-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.PCA' href='#MixedModels.PCA'><span class="jlbinding">MixedModels.PCA</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
PCA{T<:AbstractFloat}
```


Principal Components Analysis

**Fields**
- `covcorr` covariance or correlation matrix
  
- `sv` singular value decomposition
  
- `rnames` rownames of the original matrix
  
- `corr` is this a correlation matrix?
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/pca.jl#L1-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.RaggedArray' href='#MixedModels.RaggedArray'><span class="jlbinding">MixedModels.RaggedArray</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
RaggedArray{T,I}
```


A &quot;ragged&quot; array structure consisting of values and indices

**Fields**
- `vals`: a `Vector{T}` containing the values
  
- `inds`: a `Vector{I}` containing the indices
  

For this application a `RaggedArray` is used only in its `sum!` method.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L85-L95" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.ReMat' href='#MixedModels.ReMat'><span class="jlbinding">MixedModels.ReMat</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ReMat{T,S} <: AbstractMatrix{T}
```


A section of a model matrix generated by a random-effects term.

**Fields**
- `trm`: the grouping factor as a `StatsModels.CategoricalTerm`
  
- `refs`: indices into the levels of the grouping factor as a `Vector{Int32}`
  
- `levels`: the levels of the grouping factor
  
- `cnames`: the names of the columns of the model matrix generated by the left-hand side of the term
  
- `z`: transpose of the model matrix generated by the left-hand side of the term
  
- `wtz`: a weighted copy of `z` (`z` and `wtz` are the same object for unweighted cases)
  
- `λ`: a `LowerTriangular` or `Diagonal` matrix of size `S×S`
  
- `inds`: a `Vector{Int}` of linear indices of the potential nonzeros in `λ`
  
- `adjA`: the adjoint of the matrix as a `SparseMatrixCSC{T}`
  
- `scratch`: a `Matrix{T}`
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L3-L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.TableColumns' href='#MixedModels.TableColumns'><span class="jlbinding">MixedModels.TableColumns</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
TableColumns
```


A structure containing the column names for the numeric part of the profile table.

The struct also contains a Dict giving the column ranges for Symbols like `:σ` and `:β`. Finally it contains a scratch vector used to accumulate to values in a row of the profile table.

::: tip Note

This is an internal structure used in [`MixedModelProfile`](/api#MixedModels.MixedModelProfile). As such, it may change or disappear in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/utilities.jl#L1-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.UniformBlockDiagonal' href='#MixedModels.UniformBlockDiagonal'><span class="jlbinding">MixedModels.UniformBlockDiagonal</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
UniformBlockDiagonal{T}
```


Homogeneous block diagonal matrices.  `k` diagonal blocks each of size `m×m`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/arraytypes.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.VarCorr' href='#MixedModels.VarCorr'><span class="jlbinding">MixedModels.VarCorr</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
VarCorr
```


Information from the fitted random-effects variance-covariance matrices.

**Members**
- `σρ`: a `NamedTuple` of `NamedTuple`s as returned from `σρs`
  
- `s`: the estimate of the per-observation dispersion parameter
  

The main purpose of defining this type is to isolate the logic in the show method.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/varcorr.jl#L1-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Exported Functions {#Exported-Functions}
<details class='jldocstring custom-block' open>
<summary><a id='LinearAlgebra.cond-Tuple{MixedModel}' href='#LinearAlgebra.cond-Tuple{MixedModel}'><span class="jlbinding">LinearAlgebra.cond</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cond(m::MixedModel)
```


Return a vector of condition numbers of the λ matrices for the random-effects terms


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L22-L26" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='LinearAlgebra.logdet-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#LinearAlgebra.logdet-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">LinearAlgebra.logdet</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
logdet(m::LinearMixedModel)
```


Return the value of `log(det(Λ'Z'ZΛ + I)) + m.optsum.REML * log(det(LX*LX'))` evaluated in place.

Here LX is the diagonal term corresponding to the fixed-effects in the blocked lower Cholesky factor.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linalg/logdet.jl#L17-L25" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.GHnorm-Tuple{Int64}' href='#MixedModels.GHnorm-Tuple{Int64}'><span class="jlbinding">MixedModels.GHnorm</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
GHnorm(k::Int)
```


Return the (unique) GaussHermiteNormalized{k} object.

The function values are stored (memoized) when first evaluated.  Subsequent evaluations for the same `k` have very low overhead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/gausshermite.jl#L72-L79" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.coefpvalues-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T' href='#MixedModels.coefpvalues-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.coefpvalues</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
coefpvalues(bsamp::MixedModelFitCollection)
```


Return a rowtable with columns `(:iter, :coefname, :β, :se, :z, :p)`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L560-L564" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.condVar-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.condVar-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.condVar</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
condVar(m::LinearMixedModel)
```


Return the conditional variances matrices of the random effects.

The random effects are returned by `ranef` as a vector of length `k`, where `k` is the number of random effects terms.  The `i`th element is a matrix of size `vᵢ × ℓᵢ`  where `vᵢ` is the size of the vector-valued random effects for each of the `ℓᵢ` levels of the grouping factor.  Technically those values are the modes of the conditional distribution of the random effects given the observed data.

This function returns an array of `k` three dimensional arrays, where the `i`th array is of size `vᵢ × vᵢ × ℓᵢ`.  These are the diagonal blocks from the conditional variance-covariance matrix,

```
s² Λ(Λ'Z'ZΛ + I)⁻¹Λ'
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L295-L312" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.condVartables-Union{Tuple{MixedModel{T}}, Tuple{T}} where T' href='#MixedModels.condVartables-Union{Tuple{MixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.condVartables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
condVartables(m::LinearMixedModel)
```


Return the conditional covariance matrices of the random effects as a `NamedTuple` of columntables


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L343-L347" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fitted!-Union{Tuple{T}, Tuple{AbstractArray{T}, LinearMixedModel{T}}} where T' href='#MixedModels.fitted!-Union{Tuple{T}, Tuple{AbstractArray{T}, LinearMixedModel{T}}} where T'><span class="jlbinding">MixedModels.fitted!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fitted!(v::AbstractArray{T}, m::LinearMixedModel{T})
```


Overwrite `v` with the fitted values from `m`.

See also `fitted`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L510-L516" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fixef-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.fixef-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.fixef</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fixef(m::MixedModel)
```


Return the fixed-effects parameter vector estimate of `m`.

In the rank-deficient case the truncated parameter vector, of length `rank(m)` is returned. This is unlike `coef` which always returns a vector whose length matches the number of columns in `X`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L551-L559" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fixefnames-Tuple{LinearMixedModel}' href='#MixedModels.fixefnames-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.fixefnames</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fixefnames(m::MixedModel)
```


Return a (permuted and truncated in the rank-deficient case) vector of coefficient names.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L562-L566" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fnames-Tuple{MixedModel}' href='#MixedModels.fnames-Tuple{MixedModel}'><span class="jlbinding">MixedModels.fnames</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fnames(m::MixedModel)
```


Return the names of the grouping factors for the random-effects terms.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L572-L576" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fulldummy-Tuple{StatsModels.CategoricalTerm}' href='#MixedModels.fulldummy-Tuple{StatsModels.CategoricalTerm}'><span class="jlbinding">MixedModels.fulldummy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fulldummy(term::CategoricalTerm)
```


Assign &quot;contrasts&quot; that include all indicator columns (dummy variables) and an intercept column.

This will result in an under-determined set of contrasts, which is not a problem in the random effects because of the regularization, or &quot;shrinkage&quot;, of the conditional modes.

The interaction of `fulldummy` with complex random effects is subtle and complex with numerous potential edge cases. As we discover these edge cases, we will document and determine their behavior. Until such time, please check the model summary to verify that the expansion is working as you expected. If it is not, please report a use case on GitHub.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/randomeffectsterm.jl#L195-L207" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.issingular' href='#MixedModels.issingular'><span class="jlbinding">MixedModels.issingular</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
issingular(m::MixedModel, θ=m.θ; atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)
```


Test whether the model `m` is singular if the parameter vector is `θ`.

Equality comparisons are used b/c small non-negative θ values are replaced by 0 in `fit!`.

::: tip Note

For `GeneralizedLinearMixedModel`, the entire parameter vector (including β in the case `fast=false`) must be specified if the default is not used.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L60-L70" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.issingular-Tuple{MixedModels.MixedModelFitCollection}' href='#MixedModels.issingular-Tuple{MixedModels.MixedModelFitCollection}'><span class="jlbinding">MixedModels.issingular</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
issingular(bsamp::MixedModelFitCollection;
           atol::Real=0, rtol::Real=atol>0 ? 0 : √eps))
```


Test each bootstrap sample for singularity of the corresponding fit.

Equality comparisons are used b/c small non-negative θ values are replaced by 0 in `fit!`.

See also [`issingular(::MixedModel)`](/api#MixedModels.issingular).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L410-L419" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.lowerbd-Union{Tuple{ReMat{T}}, Tuple{T}} where T' href='#MixedModels.lowerbd-Union{Tuple{ReMat{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.lowerbd</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
lowerbd{T}(A::ReMat{T})
```


Return the vector of lower bounds on the parameters, `θ` associated with `A`

These are the elements in the lower triangle of `A.λ` in column-major ordering. Diagonals have a lower bound of `0`.  Off-diagonals have a lower-bound of `-Inf`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L177-L184" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.objective!' href='#MixedModels.objective!'><span class="jlbinding">MixedModels.objective!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
objective!(m::MixedModel, θ)
objective!(m::MixedModel)
```


Equivalent to `objective(updateL!(setθ!(m, θ)))`.

When `m` has a single, scalar random-effects term, `θ` can be a scalar.

The one-argument method curries and returns a single-argument function of `θ`.

Note that these methods modify `m`. The calling function is responsible for restoring the optimal `θ`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L120-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.objective-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.objective-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.objective</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
objective(m::LinearMixedModel)
```


Return negative twice the log-likelihood of model `m`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L784-L788" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.parametricbootstrap-Union{Tuple{T}, Tuple{Random.AbstractRNG, Integer, MixedModel{T}}, Tuple{Random.AbstractRNG, Integer, MixedModel{T}, Type{<:AbstractFloat}}} where T' href='#MixedModels.parametricbootstrap-Union{Tuple{T}, Tuple{Random.AbstractRNG, Integer, MixedModel{T}}, Tuple{Random.AbstractRNG, Integer, MixedModel{T}, Type{<:AbstractFloat}}} where T'><span class="jlbinding">MixedModels.parametricbootstrap</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
parametricbootstrap([rng::AbstractRNG], nsamp::Integer, m::MixedModel{T}, ftype=T;
    β = fixef(m), σ = m.σ, θ = m.θ, progress=true, optsum_overrides=(;))
```


Perform `nsamp` parametric bootstrap replication fits of `m`, returning a `MixedModelBootstrap`.

The default random number generator is `Random.GLOBAL_RNG`.

`ftype` can be used to store the computed bootstrap values in a lower precision. `ftype` is not a named argument because named arguments are not used in method dispatch and thus specialization. In other words, having `ftype` as a positional argument has some potential performance benefits.

**Keyword Arguments**
- `β`, `σ`, and `θ` are the values of `m`&#39;s parameters for simulating the responses.
  
- `σ` is only valid for `LinearMixedModel` and `GeneralizedLinearMixedModel` for
  

families with a dispersion parameter.
- `progress` controls whether the progress bar is shown. Note that the progress
  

bar is automatically disabled for non-interactive (i.e. logging) contexts.
- `optsum_overrides` is used to override values of [OptSummary](/api#MixedModels.OptSummary) in the models
  

fit during the bootstrapping process. For example, `optsum_overrides=(;ftol_rel=1e-08)` reduces the convergence criterion, which can greatly speed up the bootstrap fits. Taking advantage of this speed up to increase `n` can often lead to better estimates of coverage intervals.

::: tip Note

All coefficients are bootstrapped. In the rank deficient case, the inestimatable coefficients are treated as -0.0 in the simulations underlying the bootstrap, which will generally result in their estimate from the simulated data also being being inestimable and thus set to -0.0. **However this behavior may change in future releases to explicitly drop the extraneous columns before simulation and thus not include their estimates in the bootstrap result.**

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L178-L211" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.pirls!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}, Tuple{GeneralizedLinearMixedModel{T}, Any}, Tuple{GeneralizedLinearMixedModel{T}, Any, Any}} where T' href='#MixedModels.pirls!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}, Tuple{GeneralizedLinearMixedModel{T}, Any}, Tuple{GeneralizedLinearMixedModel{T}, Any, Any}} where T'><span class="jlbinding">MixedModels.pirls!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pirls!(m::GeneralizedLinearMixedModel)
```


Use Penalized Iteratively Reweighted Least Squares (PIRLS) to determine the conditional modes of the random effects.

When `varyβ` is true both `u` and `β` are optimized with PIRLS.  Otherwise only `u` is optimized and `β` is held fixed.

Passing `verbose = true` provides verbose output of the iterations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L574-L584" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.profile-Tuple{LinearMixedModel}' href='#MixedModels.profile-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.profile</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
profile(m::LinearMixedModel; threshold = 4)
```


Return a `MixedModelProfile` for the objective of `m` with respect to the fixed-effects coefficients.

`m` is `refit!` if `!isfitted(m)`.

Profiling starts at the parameter estimate and continues until reaching a parameter bound or the absolute value of ζ exceeds `threshold`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/profile.jl#L25-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.profilevc-Union{Tuple{T}, Tuple{LinearMixedModel{T}, T, AbstractVector{T}}} where T' href='#MixedModels.profilevc-Union{Tuple{T}, Tuple{LinearMixedModel{T}, T, AbstractVector{T}}} where T'><span class="jlbinding">MixedModels.profilevc</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 profilevc(m::LinearMixedModel{T}, val::T, rowj::AbstractVector{T}) where {T}
```


Profile an element of the variance components.

::: tip Note

This method is called by `profile` and currently considered internal. As such, it may change or disappear in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/vcpr.jl#L2-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.profileσ-Union{Tuple{T}, Tuple{LinearMixedModel{T}, MixedModels.TableColumns{T}}} where T' href='#MixedModels.profileσ-Union{Tuple{T}, Tuple{LinearMixedModel{T}, MixedModels.TableColumns{T}}} where T'><span class="jlbinding">MixedModels.profileσ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
profileσ(m::LinearMixedModel, tc::TableColumns; threshold=4)
```


Return a Table of the profile of `σ` for model `m`.  The profile extends to where the magnitude of ζ exceeds `threshold`.

::: tip Note

This method is called by `profile` and currently considered internal. As such, it may change or disappear in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/sigmapr.jl#L32-L39" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.pwrss-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.pwrss-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.pwrss</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pwrss(m::LinearMixedModel)
```


The penalized, weighted residual sum-of-squares.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L845-L849" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.ranef-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.ranef-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.ranef</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ranef(m::LinearMixedModel; uscale=false)
```


Return, as a `Vector{Matrix{T}}`, the conditional modes of the random effects in model `m`.

If `uscale` is `true` the random effects are on the spherical (i.e. `u`) scale, otherwise on the original scale.

For a named variant, see [`raneftables`](/constructors#MixedModels.raneftables).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L895-L904" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.raneftables-Union{Tuple{MixedModel{T}}, Tuple{T}} where T' href='#MixedModels.raneftables-Union{Tuple{MixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.raneftables</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
raneftables(m::MixedModel; uscale = false)
```


Return the conditional means of the random effects as a `NamedTuple` of Tables.jl-compliant tables.

::: tip Note

The API guarantee is only that the NamedTuple contains Tables.jl tables and not on the particular concrete type of each table.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L166-L173" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.refit!-Tuple{GeneralizedLinearMixedModel}' href='#MixedModels.refit!-Tuple{GeneralizedLinearMixedModel}'><span class="jlbinding">MixedModels.refit!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
refit!(m::GeneralizedLinearMixedModel[, y::Vector];
       fast::Bool = (length(m.θ) == length(m.optsum.final)),
       nAGQ::Integer = m.optsum.nAGQ,
       kwargs...)
```


Refit the model `m` after installing response `y`.

If `y` is omitted the current response vector is used.

If not specified, the `fast` and `nAGQ` options from the previous fit are used. `kwargs` are the same as [`fit!`](/api#StatsAPI.fit!-Union{Tuple{GeneralizedLinearMixedModel{T}},%20Tuple{T}}%20where%20T)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L646-L658" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.refit!-Tuple{LinearMixedModel}' href='#MixedModels.refit!-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.refit!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
refit!(m::LinearMixedModel[, y::Vector]; REML=m.optsum.REML, kwargs...)
```


Refit the model `m` after installing response `y`.

If `y` is omitted the current response vector is used. `kwargs` are the same as [`fit!`](/api#StatsAPI.fit!-Union{Tuple{GeneralizedLinearMixedModel{T}},%20Tuple{T}}%20where%20T).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L960-L967" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.replicate-Tuple{Function, Integer}' href='#MixedModels.replicate-Tuple{Function, Integer}'><span class="jlbinding">MixedModels.replicate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
replicate(f::Function, n::Integer; progress=true)
```


Return a vector of the values of `n` calls to `f()` - used in simulations where the value of `f` is stochastic.

`progress` controls whether the progress bar is shown. Note that the progress bar is automatically disabled for non-interactive (i.e. logging) contexts.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L126-L133" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.restoreoptsum!-Tuple{MixedModel, Any}' href='#MixedModels.restoreoptsum!-Tuple{MixedModel, Any}'><span class="jlbinding">MixedModels.restoreoptsum!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
restoreoptsum!(m::MixedModel, io::IO; atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)
restoreoptsum!(m::MixedModel, filename; atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)
```


Read, check, and restore the `optsum` field from a JSON stream or filename.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/serialization.jl#L1-L6" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.restorereplicates-Union{Tuple{T}, Tuple{Any, MixedModel{T}}, Tuple{Any, MixedModel{T}, Type{<:AbstractFloat}}} where T' href='#MixedModels.restorereplicates-Union{Tuple{T}, Tuple{Any, MixedModel{T}}, Tuple{Any, MixedModel{T}, Type{<:AbstractFloat}}} where T'><span class="jlbinding">MixedModels.restorereplicates</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
restorereplicates(f, m::MixedModel{T})
restorereplicates(f, m::MixedModel{T}, ftype::Type{<:AbstractFloat})
restorereplicates(f, m::MixedModel{T}, ctype::Type{<:MixedModelFitCollection{S}})
```


Restore replicates from `f`, using `m` to create the desired subtype of [`MixedModelFitCollection`](/api#MixedModels.MixedModelFitCollection).

`f` can be any entity supported by `Arrow.Table`. `m` does not have to be fitted, but it must have been constructed with the same structure as the source of the saved replicates.

The two-argument method constructs a [`MixedModelBootstrap`](/api#MixedModels.MixedModelBootstrap) with the same eltype as `m`. If an eltype is specified as the third argument, then a `MixedModelBootstrap` is returned. If a subtype of `MixedModelFitCollection` is specified as the third argument, then that is the return type.

See also [`savereplicates`](/api#MixedModels.savereplicates-Tuple{Any,%20MixedModels.MixedModelFitCollection}), [`restoreoptsum!`](/api#MixedModels.restoreoptsum!-Tuple{MixedModel,%20Any}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L71-L87" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.saveoptsum-Tuple{IO, MixedModel}' href='#MixedModels.saveoptsum-Tuple{IO, MixedModel}'><span class="jlbinding">MixedModels.saveoptsum</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
saveoptsum(io::IO, m::MixedModel)
saveoptsum(filename, m::MixedModel)
```


Save `m.optsum` (w/o the `lowerbd` field) in JSON format to an IO stream or a file

The reason for omitting the `lowerbd` field is because it often contains `-Inf` values that are not allowed in JSON.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/serialization.jl#L129-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.savereplicates-Tuple{Any, MixedModels.MixedModelFitCollection}' href='#MixedModels.savereplicates-Tuple{Any, MixedModels.MixedModelFitCollection}'><span class="jlbinding">MixedModels.savereplicates</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
savereplicates(f, b::MixedModelFitCollection)
```


Save the replicates associated with a [`MixedModelFitCollection`](/api#MixedModels.MixedModelFitCollection), e.g. [`MixedModelBootstrap`](/api#MixedModels.MixedModelBootstrap) as an Arrow file.

See also [`restorereplicates`](/api#MixedModels.restorereplicates-Union{Tuple{T},%20Tuple{Any,%20MixedModel{T}},%20Tuple{Any,%20MixedModel{T},%20Type{<:AbstractFloat}}}%20where%20T), [`saveoptsum`](/api#MixedModels.saveoptsum-Tuple{IO,%20MixedModel})

::: tip Note

**Only** the replicates are saved, not the entire contents of the `MixedModelFitCollection`. `restorereplicates` requires a model compatible with the bootstrap to restore the full object.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L120-L131" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.sdest-Tuple{LinearMixedModel}' href='#MixedModels.sdest-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.sdest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sdest(m::LinearMixedModel)
```


Return the estimate of σ, the standard deviation of the per-observation noise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L987-L991" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.sdest-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.sdest-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.sdest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sdest(m::GeneralizedLinearMixedModel)
```


Return the estimate of the dispersion, i.e. the standard deviation of the per-observation noise.

For models with a dispersion parameter ϕ, this is simply ϕ. For models without a dispersion parameter, this value is `missing`. This differs from `disperion`, which returns `1` for models without a dispersion parameter.

For Gaussian models, this parameter is often called σ.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L710-L720" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.setθ!-Union{Tuple{T}, Tuple{LinearMixedModel{T}, AbstractVector}} where T' href='#MixedModels.setθ!-Union{Tuple{T}, Tuple{LinearMixedModel{T}, AbstractVector}} where T'><span class="jlbinding">MixedModels.setθ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setθ!(m::LinearMixedModel, v)
```


Install `v` as the θ parameters in `m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L994-L998" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.setθ!-Union{Tuple{T}, Tuple{MixedModels.MixedModelFitCollection{T}, AbstractVector{T}}} where T' href='#MixedModels.setθ!-Union{Tuple{T}, Tuple{MixedModels.MixedModelFitCollection{T}, AbstractVector{T}}} where T'><span class="jlbinding">MixedModels.setθ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setθ!(bsamp::MixedModelFitCollection, θ::AbstractVector)
setθ!(bsamp::MixedModelFitCollection, i::Integer)
```


Install the values of the i&#39;th θ value of `bsamp.fits` in `bsamp.λ`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L449-L454" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.shortestcovint' href='#MixedModels.shortestcovint'><span class="jlbinding">MixedModels.shortestcovint</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
shortestcovint(v, level = 0.95)
```


Return the shortest interval containing `level` proportion of the values of `v`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L475-L479" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.shortestcovint-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}, Tuple{MixedModels.MixedModelFitCollection{T}, Any}} where T' href='#MixedModels.shortestcovint-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}, Tuple{MixedModels.MixedModelFitCollection{T}, Any}} where T'><span class="jlbinding">MixedModels.shortestcovint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
shortestcovint(bsamp::MixedModelFitCollection, level = 0.95)
```


Return the shortest interval containing `level` proportion for each parameter from `bsamp.allpars`.

::: warning Warning

Currently, correlations that are systematically zero are included in the the result. This may change in a future release without being considered a breaking change.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L495-L504" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.simulate' href='#MixedModels.simulate'><span class="jlbinding">MixedModels.simulate</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



See [`simulate!`](/api#MixedModels.simulate!-Tuple{Random.AbstractRNG,%20AbstractVector,%20LinearMixedModel,%20NamedTuple{names,%20T}%20where%20{N,%20names,%20T<:NTuple{N,%20AbstractVector}}})


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/simulate.jl#L1-L3" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.simulate!-Tuple{Random.AbstractRNG, AbstractVector, LinearMixedModel, NamedTuple{names, T} where {N, names, T<:NTuple{N, AbstractVector}}}' href='#MixedModels.simulate!-Tuple{Random.AbstractRNG, AbstractVector, LinearMixedModel, NamedTuple{names, T} where {N, names, T<:NTuple{N, AbstractVector}}}'><span class="jlbinding">MixedModels.simulate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulate!([rng::AbstractRNG,] y::AbstractVector, m::MixedModel{T}[, newdata];
                β = coef(m), σ = m.σ, θ = T[], wts=m.wts)
simulate([rng::AbstractRNG,] m::MixedModel{T}[, newdata];
                β = coef(m), σ = m.σ, θ = T[], wts=m.wts)
```


Simulate a new response vector, optionally overwriting a pre-allocated vector.

New data can be optionally provided in tabular format.

This simulation includes sampling new values for the random effects. Thus in contrast to `predict`, there is no distinction in between &quot;new&quot; and &quot;old&quot; / previously observed random-effects levels.

Unlike `predict`, there is no `type` parameter for `GeneralizedLinearMixedModel` because the noise term in the model and simulation is always on the response scale.

The `wts` argument is currently ignored except for `GeneralizedLinearMixedModel` models with a `Binomial` distribution.

::: tip Note

Note that `simulate!` methods with a `y::AbstractVector` as the first argument (besides the RNG) and `simulate` methods return the simulated response. This is in contrast to `simulate!` methods with a `m::MixedModel` as the first argument, which modify the model&#39;s response and return the entire modified model.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/simulate.jl#L96-L122" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.simulate!-Union{Tuple{T}, Tuple{Random.AbstractRNG, LinearMixedModel{T}}} where T' href='#MixedModels.simulate!-Union{Tuple{T}, Tuple{Random.AbstractRNG, LinearMixedModel{T}}} where T'><span class="jlbinding">MixedModels.simulate!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
simulate!(rng::AbstractRNG, m::MixedModel{T}; β=fixef(m), σ=m.σ, θ=T[])
simulate!(m::MixedModel; β=fixef(m), σ=m.σ, θ=m.θ)
```


Overwrite the response (i.e. `m.trms[end]`) with a simulated response vector from model `m`.

This simulation includes sampling new values for the random effects.

`β` can be specified either as a pivoted, full rank coefficient vector (cf. [`fixef`](/constructors#MixedModels.fixef)) or as an unpivoted full dimension coefficient vector (cf. [`coef`](/constructors#StatsAPI.coef)), where the entries corresponding to redundant columns will be ignored.

::: tip Note

Note that `simulate!` methods with a `y::AbstractVector` as the first argument (besides the RNG) and `simulate` methods return the simulated response. This is in contrast to `simulate!` methods with a `m::MixedModel` as the first argument, which modify the model&#39;s response and return the entire modified model.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/simulate.jl#L20-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.sparseL-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.sparseL-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.sparseL</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sparseL(m::LinearMixedModel; fname::Symbol=first(fnames(m)), full::Bool=false)
```


Return the lower Cholesky factor `L` as a `SparseMatrix{T,Int32}`.

`full` indicates whether the parts of `L` associated with the fixed-effects and response are to be included.

`fname` specifies the first grouping factor to include. Blocks to the left of the block corresponding  to `fname` are dropped. The default is the first, i.e., leftmost block and hence all blocks.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1113-L1123" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.stderror!-Union{Tuple{T}, Tuple{Tv}, Tuple{AbstractVector{Tv}, LinearMixedModel{T}}} where {Tv, T}' href='#MixedModels.stderror!-Union{Tuple{T}, Tuple{Tv}, Tuple{AbstractVector{Tv}, LinearMixedModel{T}}} where {Tv, T}'><span class="jlbinding">MixedModels.stderror!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
stderror!(v::AbstractVector, m::LinearMixedModel)
```


Overwrite `v` with the standard errors of the fixed-effects coefficients in `m`

The length of `v` should be the total number of coefficients (i.e. `length(coef(m))`). When the model matrix is rank-deficient the coefficients forced to `-0.0` have an undefined (i.e. `NaN`) standard error.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1179-L1187" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.updateL!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.updateL!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.updateL!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
updateL!(m::LinearMixedModel)
```


Update the blocked lower Cholesky factor, `m.L`, from `m.A` and `m.reterms` (used for λ only)

This is the crucial step in evaluating the objective, given a new parameter value.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1246-L1252" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.varest-Tuple{LinearMixedModel}' href='#MixedModels.varest-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.varest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
varest(m::LinearMixedModel)
```


Returns the estimate of σ², the variance of the conditional distribution of Y given B.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1287-L1291" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.varest-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.varest-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.varest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
varest(m::GeneralizedLinearMixedModel)
```


Returns the estimate of ϕ², the variance of the conditional distribution of Y given B.

For models with a dispersion parameter ϕ, this is simply ϕ². For models without a dispersion parameter, this value is `missing`. This differs from `disperion`, which returns `1` for models without a dispersion parameter.

For Gaussian models, this parameter is often called σ².


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L814-L824" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.zerocorr-Tuple{Any}' href='#MixedModels.zerocorr-Tuple{Any}'><span class="jlbinding">MixedModels.zerocorr</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
zerocorr(term::RandomEffectsTerm)
```


Remove correlations between random effects in `term`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/randomeffectsterm.jl#L231-L235" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Statistics.std-Tuple{LinearMixedModel}' href='#Statistics.std-Tuple{LinearMixedModel}'><span class="jlbinding">Statistics.std</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
std(m::MixedModel)
```


Return the estimated standard deviations of the random effects as a `Vector{Vector{T}}`.

FIXME: This uses an old convention of isfinite(sdest(m)).  Probably drop in favor of m.σs


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1166-L1172" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.confint-Tuple{MixedModelProfile}' href='#StatsAPI.confint-Tuple{MixedModelProfile}'><span class="jlbinding">StatsAPI.confint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
confint(pr::MixedModelProfile; level::Real=0.95)
```


Compute profile confidence intervals for coefficients and variance components, with confidence level level (by default 95%).

::: tip Note

The API guarantee is for a Tables.jl compatible table. The exact return type is an implementation detail and may change in a future minor release without being considered breaking.

:::

::: tip Note

The &quot;row names&quot; indicating the associated parameter name are guaranteed to be unambiguous, but their precise naming scheme is not yet stable and may change in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/profile.jl#L64-L78" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.confint-Union{Tuple{MixedModelBootstrap{T}}, Tuple{T}} where T' href='#StatsAPI.confint-Union{Tuple{MixedModelBootstrap{T}}, Tuple{T}} where T'><span class="jlbinding">StatsAPI.confint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
confint(pr::MixedModelBootstrap; level::Real=0.95, method=:shortest)
```


Compute bootstrap confidence intervals for coefficients and variance components, with confidence level level (by default 95%).

The keyword argument `method` determines whether the `:shortest`, i.e. highest density, interval is used or the `:equaltail`, i.e. quantile-based, interval is used. For historical reasons, the default is `:shortest`, but `:equaltail` gives the interval that is most comparable to the profile and Wald confidence intervals.

::: tip Note

The API guarantee is for a Tables.jl compatible table. The exact return type is an implementation detail and may change in a future minor release without being considered breaking.

:::

::: tip Note

The &quot;row names&quot; indicating the associated parameter name are guaranteed to be unambiguous, but their precise naming scheme is not yet stable and may change in a future release without being considered breaking.

:::

See also [`shortestcovint`](/bootstrap#MixedModels.shortestcovint).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L340-L360" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.confint-Union{Tuple{MixedModel{T}}, Tuple{T}} where T' href='#StatsAPI.confint-Union{Tuple{MixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">StatsAPI.confint</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
confint(pr::MixedModelProfile; level::Real=0.95)
```


Compute profile confidence intervals for (fixed effects) coefficients, with confidence level `level` (by default 95%).

::: tip Note

The API guarantee is for a Tables.jl compatible table. The exact return type is an implementation detail and may change in a future minor release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L352-L362" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.deviance-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}, Tuple{GeneralizedLinearMixedModel{T}, Any}} where T' href='#StatsAPI.deviance-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}, Tuple{GeneralizedLinearMixedModel{T}, Any}} where T'><span class="jlbinding">StatsAPI.deviance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
deviance(m::GeneralizedLinearMixedModel{T}, nAGQ=1)::T where {T}
```


Return the deviance of `m` evaluated by the Laplace approximation (`nAGQ=1`) or `nAGQ`-point adaptive Gauss-Hermite quadrature.

If the distribution `D` does not have a scale parameter the Laplace approximation is the squared length of the conditional modes, $u$, plus the determinant of $Λ'Z'WZΛ + I$, plus the sum of the squared deviance residuals.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L75-L84" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.dof_residual-Tuple{MixedModel}' href='#StatsAPI.dof_residual-Tuple{MixedModel}'><span class="jlbinding">StatsAPI.dof_residual</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dof_residual(m::MixedModel)
```


Return the residual degrees of freedom of the model.

::: tip Note

The residual degrees of freedom for mixed-effects models is not clearly defined due to partial pooling. The classical `nobs(m) - dof(m)` fails to capture the extra freedom granted by the random effects, but `nobs(m) - nranef(m)` would overestimate the freedom granted by the random effects. `nobs(m) - sum(leverage(m))` provides a nice balance based on the relative influence of each observation, but is computationally expensive for large models. This problem is also fundamentally related to [long-standing debates](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#why-doesnt-lme4-display-denominator-degrees-of-freedomp-values-what-other-options-do-i-have) about the appropriate treatment of the denominator degrees of freedom for $F$-tests. In the future, MixedModels.jl may provide additional methods allowing the user to choose the computation to use.

:::

::: warning Warning

Currently, the residual degrees of freedom is computed as `nobs(m) - dof(m)`, but this may change in the future without being considered a breaking change because there is no canonical definition of the residual degrees of freedom in a mixed-effects model.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L33-L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.fit!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T' href='#StatsAPI.fit!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">StatsAPI.fit!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fit!(m::GeneralizedLinearMixedModel; fast=false, nAGQ=1,
                                     verbose=false, progress=true,
                                     thin::Int=1,
                                     init_from_lmm=Set())
```


Optimize the objective function for `m`.

When `fast` is `true` a potentially much faster but slightly less accurate algorithm, in which `pirls!` optimizes both the random effects and the fixed-effects parameters, is used.

If `progress` is `true`, the default, a `ProgressMeter.ProgressUnknown` counter is displayed. during the iterations to minimize the deviance.  There is a delay before this display is initialized and it may not be shown at all for models that are optimized quickly.

If `verbose` is `true`, then both the intermediate results of both the nonlinear optimization and PIRLS are also displayed on standard output.

The `thin` argument is ignored: it had no impact on the final model fit and the logic around thinning the `fitlog` was needlessly complicated for a trivial performance gain.

By default, the starting values for model fitting are taken from a (non mixed, i.e. marginal ) GLM fit. Experience with larger datasets (many thousands of observations and/or hundreds of levels of the grouping variables) has suggested that fitting a (Gaussian) linear mixed model on the untransformed data may provide better starting values and thus overall faster fits even though an entire LMM must be fit before the GLMM can be fit. `init_from_lmm` can be used to specify which starting values from an LMM to use. Valid options are any collection (array, set, etc.) containing one or more of `:β` and `:θ`, the default is the empty set.

::: tip Note

Initializing from an LMM requires fitting the entire LMM first, so when `progress=true`, there will be two progress bars: first for the LMM, then for the GLMM.

:::

::: warning Warning

The `init_from_lmm` functionality is experimental and may change or be removed entirely without being considered a breaking change.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L204-L243" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.fit!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#StatsAPI.fit!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">StatsAPI.fit!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fit!(m::LinearMixedModel; progress::Bool=true, REML::Bool=m.optsum.REML,
                          σ::Union{Real, Nothing}=m.optsum.sigma,
                          thin::Int=typemax(Int),
                          fitlog::Bool=true)
```


Optimize the objective of a `LinearMixedModel`.  When `progress` is `true` a `ProgressMeter.ProgressUnknown` display is shown during the optimization of the objective, if the optimization takes more than one second or so.

The `thin` argument is ignored: it had no impact on the final model fit and the logic around thinning the `fitlog` was needlessly complicated for a trivial performance gain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L427-L439" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.leverage-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#StatsAPI.leverage-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">StatsAPI.leverage</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
leverage(::LinearMixedModel)
```


Return the diagonal of the hat matrix of the model.

For a linear model, the sum of the leverage values is the degrees of freedom for the model in the sense that this sum is the dimension of the span of columns of the model matrix.  With a bit of hand waving a similar argument could be made for linear mixed-effects models. The hat matrix is of the form $[ZΛ X][L L']⁻¹[ZΛ X]'$.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L708-L717" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.modelmatrix-Tuple{MixedModel}' href='#StatsAPI.modelmatrix-Tuple{MixedModel}'><span class="jlbinding">StatsAPI.modelmatrix</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
modelmatrix(m::MixedModel)
```


Returns the model matrix `X` for the fixed-effects parameters, as returned by [`coef`](/constructors#StatsAPI.coef).

This is always the full model matrix in the original column order and from a field in the model struct.  It should be copied if it is to be modified.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L108-L115" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.predict-Tuple{LinearMixedModel, NamedTuple{names, T} where {N, names, T<:NTuple{N, AbstractVector}}}' href='#StatsAPI.predict-Tuple{LinearMixedModel, NamedTuple{names, T} where {N, names, T<:NTuple{N, AbstractVector}}}'><span class="jlbinding">StatsAPI.predict</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
StatsAPI.predict(m::LinearMixedModel, newdata;
                new_re_levels=:missing)
StatsAPI.predict(m::GeneralizedLinearMixedModel, newdata;
                new_re_levels=:missing, type=:response)
```


Predict response for new data.

::: tip Note

Currently, no in-place methods are provided because these methods internally construct a new model and therefore allocate not just a response vector but also many other matrices.

:::

::: warning Warning

`newdata` should contain a column for the response (dependent variable) initialized to some numerical value (not `missing`), because this is used to construct the new model used in computing the predictions. `missing` is not valid because `missing` data are dropped before constructing the model matrices.

:::

::: warning Warning

These methods construct an entire MixedModel behind the scenes and as such may use a large amount of memory when `newdata` is large.

:::

::: warning Warning

Rank-deficiency can lead to surprising but consistent behavior. For example, if there are two perfectly collinear predictors `A` and `B` (e.g. constant multiples of each other), then it is possible that `A` will be pivoted out in the fitted model and thus the associated coefficient is set to zero. If predictions are then generated on new data where `B` has been set to zero but `A` has not, then there will no contribution from neither `A` nor `B` in the resulting predictions.

:::

The keyword argument `new_re_levels` specifies how previously unobserved values of the grouping variable are handled. Possible values are:
- `:population`: return population values for the relevant grouping variable.  In other words, treat the associated random effect as 0.  If all grouping variables have new levels, then this is equivalent to  just the fixed effects.
  
- `:missing`: return `missing`.
  
- `:error`: error on this condition. The error type is an implementation detail:  you should not rely on a particular type of error being thrown.
  

If you want simulated values for unobserved levels of the grouping variable, consider the [`simulate!`](/api#MixedModels.simulate!-Tuple{Random.AbstractRNG,%20AbstractVector,%20LinearMixedModel,%20NamedTuple{names,%20T}%20where%20{N,%20names,%20T<:NTuple{N,%20AbstractVector}}}) and `simulate` methods.

Predictions based purely on the fixed effects can be obtained by specifying previously unobserved levels of the random effects and setting `new_re_levels=:population`. Similarly, the contribution of any grouping variable can be excluded by specifying previously unobserved levels, while including previously observed levels of the other grouping variables. In the future, it may be possible to specify a subset of the grouping variables or overall random-effects structure to use, but not at this time.

::: tip Note

`new_re_levels` impacts only the behavior for previously unobserved random effects levels, i.e. new RE levels. For previously observed random effects levels, predictions take both the fixed and random effects into account.

:::

For `GeneralizedLinearMixedModel`, the `type` parameter specifies whether the predictions should be returned on the scale of linear predictor (`:linpred`) or on the response scale (`:response`). If you don&#39;t know the difference between these terms, then you probably want `type=:response`.

Regression weights are not yet supported in prediction. Similarly, offsets are also not supported for `GeneralizedLinearMixedModel`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/predict.jl#L1-L69" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.response-Tuple{MixedModel}' href='#StatsAPI.response-Tuple{MixedModel}'><span class="jlbinding">StatsAPI.response</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
response(m::MixedModel)
```


Return the response vector for the model.

For a linear mixed model this is a `view` of the last column of the `XyMat` field. For a generalized linear mixed model this is the `m.resp.y` field. In either case it should be copied if it is to be modified.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L182-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.vcov-Tuple{MixedModel}' href='#StatsAPI.vcov-Tuple{MixedModel}'><span class="jlbinding">StatsAPI.vcov</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
vcov(m::MixedModel; corr=false)
```


Returns the variance-covariance matrix of the fixed effects. If `corr` is `true`, the correlation of the fixed effects is returned instead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L222-L227" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Tables.columntable-Tuple{OptSummary}' href='#Tables.columntable-Tuple{OptSummary}'><span class="jlbinding">Tables.columntable</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
columntable(s::OptSummary, [stack::Bool=false])
```


Return `s.fitlog` as a `Tables.columntable`.

When `stack` is false (the default), there will be 3 columns in the result:
- `iter`: the iteration number
  
- `objective`: the value of the objective at that sample
  
- `θ`: the parameter vector at that sample
  

When `stack` is true, there will be 4 columns: `iter`, `objective`, `par`, and `value` where `value` is the stacked contents of the `θ` vectors (the equivalent of `vcat(θ...)`) and `par` is a vector of parameter numbers.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/optsummary.jl#L94-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Methods from `StatsAPI.jl`, `StatsBase.jl`, `StatsModels.jl` and `GLM.jl` {#Methods-from-StatsAPI.jl,-StatsBase.jl,-StatsModels.jl-and-GLM.jl}

```julia
aic
aicc
bic
coef
coefnames
coeftable
deviance
dispersion
dispersion_parameter
dof
dof_residual
fit
fit!
fitted
formula
isfitted
islinear
leverage
loglikelihood
meanresponse
modelmatrix
model_response
nobs
predict
residuals
response
responsename
StatsModels.lrtest # not exported
std
stderror
vcov
weights
```


### MixedModels.jl &quot;alternatives&quot; and extensions to StatsAPI and GLM functions {#MixedModels.jl-"alternatives"-and-extensions-to-StatsAPI-and-GLM-functions}

The following are MixedModels.jl-specific functions and not simply methods for functions defined in StatsAPI and GLM.jl.

```julia
coefpvalues
condVar
condVarTables
fitted!
fixef
fixefnames
likelihoodratiotest # not exported
pwrss
ranef
raneftables
refit!
shortestcovint
sdest
simulate
simulate!
stderrror!
varest
```


## Non-Exported Functions {#Non-Exported-Functions}

Note that unless discussed elsewhere in the online documentation, non-exported functions should be considered implementation details.
<details class='jldocstring custom-block' open>
<summary><a id='Base.copy-Union{Tuple{ReMat{T, S}}, Tuple{S}, Tuple{T}} where {T, S}' href='#Base.copy-Union{Tuple{ReMat{T, S}}, Tuple{S}, Tuple{T}} where {T, S}'><span class="jlbinding">Base.copy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Base.copy(ReMat{T,S})
```


Return a shallow copy of ReMat.

A shallow copy shares as much internal storage as possible with the original ReMat. Only the vector `λ` and the `scratch` matrix are copied.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/fixefpr.jl#L9-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.size-Tuple{MixedModel}' href='#Base.size-Tuple{MixedModel}'><span class="jlbinding">Base.size</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
size(m::MixedModel)
```


Returns the size of a mixed model as a tuple of length four: the number of observations, the number of (non-singular) fixed-effects parameters, the number of conditional modes (random effects), the number of grouping variables


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L210-L216" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GLM.wrkresp!-Union{Tuple{T}, Tuple{AbstractVector{T}, GLM.GlmResp{Vector{T}}}} where T<:AbstractFloat' href='#GLM.wrkresp!-Union{Tuple{T}, Tuple{AbstractVector{T}, GLM.GlmResp{Vector{T}}}} where T<:AbstractFloat'><span class="jlbinding">GLM.wrkresp!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
GLM.wrkresp!(v::AbstractVector{T}, resp::GLM.GlmResp{AbstractVector{T}})
```


A copy of a method from GLM that generalizes the types in the signature


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L122-L126" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.LD-Union{Tuple{LinearAlgebra.Diagonal{T, V} where V<:AbstractVector{T}}, Tuple{T}} where T<:Number' href='#MixedModels.LD-Union{Tuple{LinearAlgebra.Diagonal{T, V} where V<:AbstractVector{T}}, Tuple{T}} where T<:Number'><span class="jlbinding">MixedModels.LD</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
LD(A::Diagonal)
LD(A::HBlikDiag)
LD(A::DenseMatrix)
```


Return `log(det(tril(A)))` evaluated in place.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linalg/logdet.jl#L1-L7" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.adjA-Tuple{AbstractVector, AbstractMatrix}' href='#MixedModels.adjA-Tuple{AbstractVector, AbstractMatrix}'><span class="jlbinding">MixedModels.adjA</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
adjA(refs::AbstractVector, z::AbstractMatrix{T})
```


Returns the adjoint of an `ReMat` as a `SparseMatrixCSC{T,Int32}`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L87-L91" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.allpars-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T' href='#MixedModels.allpars-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.allpars</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
allpars(bsamp::MixedModelFitCollection)
```


Return a tidy (column)table with the parameter estimates spread into columns of `iter`, `type`, `group`, `name` and `value`.

::: warning Warning

Currently, correlations that are systematically zero are included in the the result. This may change in a future release without being considered a breaking change.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L279-L289" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.amalgamate-Union{Tuple{Vector{<:AbstractReMat{T}}}, Tuple{T}} where T' href='#MixedModels.amalgamate-Union{Tuple{Vector{<:AbstractReMat{T}}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.amalgamate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
amalgamate(reterms::Vector{AbstractReMat})
```


Combine multiple ReMat with the same grouping variable into a single object.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L33-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.average-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractFloat' href='#MixedModels.average-Union{Tuple{T}, Tuple{T, T}} where T<:AbstractFloat'><span class="jlbinding">MixedModels.average</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
average(a::T, b::T) where {T<:AbstractFloat}
```


Return the average of `a` and `b`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L46-L50" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.block-Tuple{Integer, Integer}' href='#MixedModels.block-Tuple{Integer, Integer}'><span class="jlbinding">MixedModels.block</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
block(i, j)
```


Return the linear index of the `[i,j]` position (&quot;block&quot;) in the row-major packed lower triangle.

Use the row-major ordering in this case because the result depends only on `i` and `j`, not on the overall size of the array.

When `i == j` the value is the same as `kp1choose2(i)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/blocks.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.cholUnblocked!' href='#MixedModels.cholUnblocked!'><span class="jlbinding">MixedModels.cholUnblocked!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
cholUnblocked!(A, Val{:L})
```


Overwrite the lower triangle of `A` with its lower Cholesky factor.

The name is borrowed from [https://github.com/andreasnoack/LinearAlgebra.jl] because these are part of the inner calculations in a blocked Cholesky factorization.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linalg/cholUnblocked.jl#L1-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.copyscaleinflate!' href='#MixedModels.copyscaleinflate!'><span class="jlbinding">MixedModels.copyscaleinflate!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
copyscaleinflate!(L::AbstractMatrix, A::AbstractMatrix, Λ::ReMat)
```


Overwrite L with `Λ'AΛ + I`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L561-L565" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.corrmat-Union{Tuple{ReMat{T}}, Tuple{T}} where T' href='#MixedModels.corrmat-Union{Tuple{ReMat{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.corrmat</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
corrmat(A::ReMat)
```


Return the estimated correlation matrix for `A`.  The diagonal elements are 1 and the off-diagonal elements are the correlations between those random effect terms

**Example**

Note that trailing digits may vary slightly depending on the local platform.

```julia
julia> using MixedModels

julia> mod = fit(MixedModel,
                 @formula(rt_trunc ~ 1 + spkr + prec + load + (1 + spkr + prec | subj)),
                 MixedModels.dataset(:kb07));

julia> VarCorr(mod)
Variance components:
             Column      Variance  Std.Dev.  Corr.
subj     (Intercept)     136591.782 369.583
         spkr: old        22922.871 151.403 +0.21
         prec: maintain   32348.269 179.856 -0.98 -0.03
Residual                 642324.531 801.452

julia> MixedModels.corrmat(mod.reterms[1])
3×3 LinearAlgebra.Symmetric{Float64,Array{Float64,2}}:
  1.0        0.214816   -0.982948
  0.214816   1.0        -0.0315607
 -0.982948  -0.0315607   1.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L729-L761" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.cpad-Tuple{String, Integer}' href='#MixedModels.cpad-Tuple{String, Integer}'><span class="jlbinding">MixedModels.cpad</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cpad(s::AbstractString, n::Integer)
```


Return a string of length `n` containing `s` in the center (more-or-less).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L53-L57" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.densify' href='#MixedModels.densify'><span class="jlbinding">MixedModels.densify</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
densify(S::SparseMatrix, threshold=0.1)
```


Convert sparse `S` to `Diagonal` if `S` is diagonal or to `Array(S)` if the proportion of nonzeros exceeds `threshold`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L60-L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.deviance!' href='#MixedModels.deviance!'><span class="jlbinding">MixedModels.deviance!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
deviance!(m::GeneralizedLinearMixedModel, nAGQ=1)
```


Update `m.η`, `m.μ`, etc., install the working response and working weights in `m.LMM`, update `m.LMM.A` and `m.LMM.R`, then evaluate the `deviance`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L135-L140" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.feL-Tuple{LinearMixedModel}' href='#MixedModels.feL-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.feL</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
feL(m::LinearMixedModel)
```


Return the lower Cholesky factor for the fixed-effects parameters, as an `LowerTriangular` `p × p` matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L414-L419" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fixef!-Union{Tuple{T}, Tuple{Tv}, Tuple{AbstractVector{Tv}, LinearMixedModel{T}}} where {Tv, T}' href='#MixedModels.fixef!-Union{Tuple{T}, Tuple{Tv}, Tuple{AbstractVector{Tv}, LinearMixedModel{T}}} where {Tv, T}'><span class="jlbinding">MixedModels.fixef!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fixef!(v::Vector{T}, m::MixedModel{T})
```


Overwrite `v` with the pivoted fixed-effects coefficients of model `m`

For full-rank models the length of `v` must be the rank of `X`.  For rank-deficient models the length of `v` can be the rank of `X` or the number of columns of `X`.  In the latter case the calculated coefficients are padded with -0.0 out to the number of columns.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L529-L537" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fname-Tuple{ReMat}' href='#MixedModels.fname-Tuple{ReMat}'><span class="jlbinding">MixedModels.fname</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fname(A::ReMat)
```


Return the name of the grouping factor as a `Symbol`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L121-L125" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.getθ!-Union{Tuple{T}, Tuple{AbstractVector{T}, ReMat{T}}} where T' href='#MixedModels.getθ!-Union{Tuple{T}, Tuple{AbstractVector{T}, ReMat{T}}} where T'><span class="jlbinding">MixedModels.getθ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getθ!(v::AbstractVector{T}, A::ReMat{T}) where {T}
```


Overwrite `v` with the elements of the blocks in the lower triangle of `A.Λ` (column-major ordering)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L132-L136" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.getθ-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.getθ-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.getθ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
getθ(m::LinearMixedModel)
```


Return the current covariance parameter vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L592-L596" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.indmat' href='#MixedModels.indmat'><span class="jlbinding">MixedModels.indmat</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
indmat(A::ReMat)
```


Return a `Bool` indicator matrix of the potential non-zeros in `A.λ`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L158-L162" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.isconstant-Tuple{Any}' href='#MixedModels.isconstant-Tuple{Any}'><span class="jlbinding">MixedModels.isconstant</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
isconstant(x::Array)
isconstant(x::Tuple)
```


Are all elements of the iterator the same?  That is, is it constant?


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L29-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.isfullrank-Tuple{MixedModels.FeTerm}' href='#MixedModels.isfullrank-Tuple{MixedModels.FeTerm}'><span class="jlbinding">MixedModels.isfullrank</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
isfullrank(A::FeTerm)
```


Does `A` have full column rank?


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L82-L86" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.isnested-Tuple{ReMat, ReMat}' href='#MixedModels.isnested-Tuple{ReMat, ReMat}'><span class="jlbinding">MixedModels.isnested</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
isnested(A::ReMat, B::ReMat)
```


Is the grouping factor for `A` nested in the grouping factor for `B`?

That is, does each value of `A` occur with just one value of B?


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L190-L196" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.kchoose2-Tuple{Any}' href='#MixedModels.kchoose2-Tuple{Any}'><span class="jlbinding">MixedModels.kchoose2</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
kchoose2(k)
```


The binomial coefficient `k` choose `2` which is the number of elements in the packed form of the strict lower triangle of a matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/blocks.jl#L16-L21" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.kp1choose2-Tuple{Any}' href='#MixedModels.kp1choose2-Tuple{Any}'><span class="jlbinding">MixedModels.kp1choose2</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
kp1choose2(k)
```


The binomial coefficient `k+1` choose `2` which is the number of elements in the packed form of the lower triangle of a matrix.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/blocks.jl#L26-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.likelihoodratiotest-Tuple{Vararg{MixedModel}}' href='#MixedModels.likelihoodratiotest-Tuple{Vararg{MixedModel}}'><span class="jlbinding">MixedModels.likelihoodratiotest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
likelihoodratiotest(m::MixedModel...)
likelihoodratiotest(m0::LinearModel, m::MixedModel...)
likelihoodratiotest(m0::GeneralizedLinearModel, m::MixedModel...)
likelihoodratiotest(m0::TableRegressionModel{LinearModel}, m::MixedModel...)
likelihoodratiotest(m0::TableRegressionModel{GeneralizedLinearModel}, m::MixedModel...)
```


Likeihood ratio test applied to a set of nested models.

::: tip Note

The nesting of the models is not checked.  It is incumbent on the user to check this. This differs from `StatsModels.lrtest` as nesting in mixed models, especially in the random effects specification, may be non obvious.

:::

::: tip Note

For comparisons between mixed and non-mixed models, the deviance for the non-mixed model is taken to be -2 log likelihood, i.e. omitting the additive constant for the fully saturated model. This is in line with the computation of the deviance for mixed models.

:::

This functionality may be deprecated in the future in favor of `StatsModels.lrtest`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/likelihoodratiotest.jl#L46-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.nranef-Tuple{ReMat}' href='#MixedModels.nranef-Tuple{ReMat}'><span class="jlbinding">MixedModels.nranef</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nranef(A::ReMat)
```


Return the number of random effects represented by `A`.  Zero unless `A` is an `ReMat`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L110-L114" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.nθ-Tuple{ReMat}' href='#MixedModels.nθ-Tuple{ReMat}'><span class="jlbinding">MixedModels.nθ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
nθ(A::ReMat)
```


Return the number of free parameters in the relative covariance matrix λ


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L170-L174" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.opt_params' href='#MixedModels.opt_params'><span class="jlbinding">MixedModels.opt_params</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
opt_params(::Val{backend})
```


Return a collection of the fields of [`OptSummary`](/api#MixedModels.OptSummary) used by backend.

`:xtol_zero_abs`, `:ftol_zero_abs` do not need to be specified because they are used _after_ optimization and are thus shared across backends.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/optsummary.jl#L179-L186" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.optimize!' href='#MixedModels.optimize!'><span class="jlbinding">MixedModels.optimize!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
optimize!(::LinearMixedModel, ::Val{backend}; kwargs...)
optimize!(::GeneralizedLinearMixedModel, ::Val{backend}; kwargs...)
```


Perform optimization on a mixed model, minimizing the objective.

`optimize!` set ups the call to the backend optimizer using the options contained in the model&#39;s `optsum` field. It then calls that optimizer and returns `xmin, fmin`. Providing support for a new backend involves defining appropriate `optimize!` methods with the second argument of type `::Val{:backend_name}` and adding `:backend_name` to `OPTIMIZATION_BACKENDS`. Similarly, a method `opt_params(::Val{:backend_name})` should be defined, which returns the optimization parameters (e.g. `xtol_abs` or `rho_end`) used by the backend.

Common keyword arguments are `progress` to show a progress meter as well as `nAQG` and `fast` for GLMMs.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/optsummary.jl#L160-L176" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.optsumj-Tuple{OptSummary, Integer}' href='#MixedModels.optsumj-Tuple{OptSummary, Integer}'><span class="jlbinding">MixedModels.optsumj</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
optsumj(os::OptSummary, j::Integer)
```


Return an `OptSummary` with the `j`&#39;th component of the parameter omitted.

`os.final` with its j&#39;th component omitted is used as the initial parameter.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/thetapr.jl#L2-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.parsej-Tuple{Symbol}' href='#MixedModels.parsej-Tuple{Symbol}'><span class="jlbinding">MixedModels.parsej</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
parsej(sym::Symbol)
```


Return the index from symbol names like `:θ1`, `:θ01`, etc.

::: tip Note

This method is internal.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/utilities.jl#L87-L94" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.pivot-Tuple{MixedModel}' href='#MixedModels.pivot-Tuple{MixedModel}'><span class="jlbinding">MixedModels.pivot</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pivot(m::MixedModel)
pivot(A::FeTerm)
```


Return the pivot associated with the FeTerm.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/Xymat.jl#L64-L69" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.prfit!' href='#MixedModels.prfit!'><span class="jlbinding">MixedModels.prfit!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
prfit!(m::LinearMixedModel; kwargs...)
```


Fit a mixed model using the [PRIMA](https://github.com/libprima/PRIMA.jl) implementation of the BOBYQA optimizer.

::: warning Experimental feature

This function is an experimental feature that will go away in the future. Do **not** rely on it, unless you are willing to pin the precise MixedModels.jl version. The purpose of the function is to provide the MixedModels developers a chance to explore the performance of the PRIMA implementation without the large and potentially breaking changes it would take to fully replace the current NLopt backend with a PRIMA backend or a backend supporting a range of optimizers.

:::

::: tip Package extension

In order to reduce the dependency burden, all methods of this function are implemented in a package extension and are only defined when PRIMA.jl is loaded by the user.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/prima.jl#L2-L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.profileσs!-Union{Tuple{T}, Tuple{NamedTuple, MixedModels.TableColumns{T}}} where T' href='#MixedModels.profileσs!-Union{Tuple{T}, Tuple{NamedTuple, MixedModels.TableColumns{T}}} where T'><span class="jlbinding">MixedModels.profileσs!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 profileσs!(val::NamedTuple, tc::TableColumns{T}; nzlb=1.0e-8) where {T}
```


Profile the variance components.

::: tip Note

This method is called by `profile` and currently considered internal. As such, it may change or disappear in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/vcpr.jl#L27-L35" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.ranef!-Union{Tuple{T}, Tuple{Vector, LinearMixedModel{T}, AbstractArray{T}, Bool}} where T' href='#MixedModels.ranef!-Union{Tuple{T}, Tuple{Vector, LinearMixedModel{T}, AbstractArray{T}, Bool}} where T'><span class="jlbinding">MixedModels.ranef!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ranef!(v::Vector{Matrix{T}}, m::MixedModel{T}, β, uscale::Bool) where {T}
```


Overwrite `v` with the conditional modes of the random effects for `m`.

If `uscale` is `true` the random effects are on the spherical (i.e. `u`) scale, otherwise on the original scale

`β` is the truncated, pivoted coefficient vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L852-L861" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.rankUpdate!' href='#MixedModels.rankUpdate!'><span class="jlbinding">MixedModels.rankUpdate!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
rankUpdate!(C, A)
rankUpdate!(C, A, α)
rankUpdate!(C, A, α, β)
```


A rank-k update, C := α_A&#39;A + β_C, of a Hermitian (Symmetric) matrix.

`α` and `β` both default to 1.0.  When `α` is -1.0 this is a downdate operation. The name `rankUpdate!` is borrowed from [https://github.com/andreasnoack/LinearAlgebra.jl]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linalg/rankUpdate.jl#L1-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.rePCA-Tuple{LinearMixedModel}' href='#MixedModels.rePCA-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.rePCA</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rePCA(m::LinearMixedModel; corr::Bool=true)
```


Return a named tuple of the normalized cumulative variance of a principal components analysis of the random effects covariance matrices or correlation matrices when `corr` is `true`.

The normalized cumulative variance is the proportion of the variance for the first principal component, the first two principal components, etc.  The last element is always 1.0 representing the complete proportion of the variance.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L913-L923" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.reevaluateAend!-Tuple{LinearMixedModel}' href='#MixedModels.reevaluateAend!-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.reevaluateAend!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reevaluateAend!(m::LinearMixedModel)
```


Reevaluate the last column of `m.A` from `m.Xymat`.  This function should be called after updating the response.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L940-L945" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.refitσ!-Union{Tuple{T}, Tuple{LinearMixedModel{T}, Any, MixedModels.TableColumns{T}, T, Bool}} where T' href='#MixedModels.refitσ!-Union{Tuple{T}, Tuple{LinearMixedModel{T}, Any, MixedModels.TableColumns{T}, T, Bool}} where T'><span class="jlbinding">MixedModels.refitσ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
refitσ!(m::LinearMixedModel{T}, σ::T, tc::TableColumns{T}, obj::T, neg::Bool)
```


Refit the model `m` with the given value of `σ` and return a NamedTuple of information about the fit.

`obj` and `neg` allow for conversion of the objective to the `ζ` scale and `tc` is used to return a NamedTuple

::: tip Note

This method is internal and may change or disappear in a future release without being considered breaking.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/profile/sigmapr.jl#L1-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.schematize' href='#MixedModels.schematize'><span class="jlbinding">MixedModels.schematize</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
schematize(f, tbl, contrasts::Dict{Symbol}, Mod=LinearMixedModel)
```


Find and apply the schema for f in a way that automatically uses `Grouping()` contrasts when appropriate.

::: tip Warn

This is an internal method.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/schema.jl#L53-L61" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.sdcorr-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T' href='#MixedModels.sdcorr-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.sdcorr</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
sdcorr(A::AbstractMatrix{T}) where {T}
```


Transform a square matrix `A` with positive diagonals into an `NTuple{size(A,1), T}` of standard deviations and a tuple of correlations.

`A` is assumed to be symmetric and only the lower triangle is used.  The order of the correlations is row-major ordering of the lower triangle (or, equivalently, column-major in the upper triangle).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/utilities.jl#L163-L172" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.setβθ!-Tuple{GeneralizedLinearMixedModel, Any}' href='#MixedModels.setβθ!-Tuple{GeneralizedLinearMixedModel, Any}'><span class="jlbinding">MixedModels.setβθ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
setβθ!(m::GeneralizedLinearMixedModel, v)
```


Set the parameter vector, `:βθ`, of `m` to `v`.

`βθ` is the concatenation of the fixed-effects, `β`, and the covariance parameter, `θ`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L675-L681" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.ssqdenom-Tuple{LinearMixedModel}' href='#MixedModels.ssqdenom-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.ssqdenom</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ssqdenom(m::LinearMixedModel)
```


Return the denominator for penalized sums-of-squares.

For MLE, this value is the number of observations. For REML, this value is the number of observations minus the rank of the fixed-effects matrix. The difference is analogous to the use of n or n-1 in the denominator when calculating the variance.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1151-L1160" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.statsrank-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:AbstractFloat' href='#MixedModels.statsrank-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:AbstractFloat'><span class="jlbinding">MixedModels.statsrank</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
statsrank(x::Matrix{T}, ranktol::Real=1e-8) where {T<:AbstractFloat}
```


Return the numerical column rank and a pivot vector.

The rank is determined from the absolute values of the diagonal of R from a pivoted QR decomposition, relative to the first (and, hence, largest) element of this vector.

In the full-rank case the pivot vector is `collect(axes(x, 2))`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linalg/pivot.jl#L1-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.tidyβ-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T' href='#MixedModels.tidyβ-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.tidyβ</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tidyβ(bsamp::MixedModelFitCollection)
```


Return a tidy (row)table with the parameter estimates spread into columns of `iter`, `coefname` and `β`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L541-L545" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.tidyσs-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T' href='#MixedModels.tidyσs-Union{Tuple{MixedModels.MixedModelFitCollection{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.tidyσs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tidyσs(bsamp::MixedModelFitCollection)
```


Return a tidy (row)table with the estimates of the variance components (on the standard deviation scale) spread into columns of `iter`, `group`, `column` and `σ`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L581-L586" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.unfit!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.unfit!-Union{Tuple{LinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.unfit!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
unfit!(model::MixedModel)
```


Mark a model as unfitted.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1233-L1237" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.unscaledre!' href='#MixedModels.unscaledre!'><span class="jlbinding">MixedModels.unscaledre!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
unscaledre!(y::AbstractVector{T}, M::ReMat{T}) where {T}
unscaledre!(rng::AbstractRNG, y::AbstractVector{T}, M::ReMat{T}) where {T}
```


Add unscaled random effects simulated from `M` to `y`.

These are unscaled random effects (i.e. they incorporate λ but not σ) because the scaling is done after the per-observation noise is added as a standard normal.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/simulate.jl#L296-L304" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.updateA!-Tuple{LinearMixedModel}' href='#MixedModels.updateA!-Tuple{LinearMixedModel}'><span class="jlbinding">MixedModels.updateA!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
updateA!(m::LinearMixedModel)
```


Update the cross-product array, `m.A`, from `m.reterms` and `m.Xymat`

This is usually done after a reweight! operation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1206-L1212" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.updateη!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T' href='#MixedModels.updateη!-Union{Tuple{GeneralizedLinearMixedModel{T}}, Tuple{T}} where T'><span class="jlbinding">MixedModels.updateη!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
updateη!(m::GeneralizedLinearMixedModel)
```


Update the linear predictor, `m.η`, from the offset and the `B`-scale random effects.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L796-L800" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.σvals!-Tuple{AbstractVector, ReMat, Number}' href='#MixedModels.σvals!-Tuple{AbstractVector, ReMat, Number}'><span class="jlbinding">MixedModels.σvals!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
σvals!(v::AbstractVector, A::ReMat, sc::Number)
```


Overwrite v with the standard deviations of the random effects associated with `A`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/remat.jl#L634-L638" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.σρ!-Union{Tuple{T}, Tuple{AbstractVector{<:Union{Missing, T}}, LinearAlgebra.LowerTriangular, Any}} where T' href='#MixedModels.σρ!-Union{Tuple{T}, Tuple{AbstractVector{<:Union{Missing, T}}, LinearAlgebra.LowerTriangular, Any}} where T'><span class="jlbinding">MixedModels.σρ!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
σρ!(v, t, σ)
```


push! `σ` times the row lengths (σs) and the inner products of normalized rows (ρs) of `t` onto `v`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/bootstrap.jl#L639-L643" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsModels.isnested-Tuple{MixedModel, MixedModel}' href='#StatsModels.isnested-Tuple{MixedModel, MixedModel}'><span class="jlbinding">StatsModels.isnested</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
isnested(m1::MixedModel, m2::MixedModel; atol::Real=0.0)
```


Indicate whether model `m1` is nested in model `m2`, i.e. whether `m1` can be obtained by constraining some parameters in `m2`. Both models must have been fitted on the same data. This check is conservative for `MixedModel`s and may reject nested models with different parameterizations as being non nested.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/likelihoodratiotest.jl#L255-L262" target="_blank" rel="noreferrer">source</a></Badge>

</details>

