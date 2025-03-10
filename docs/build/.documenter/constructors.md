
# Model constructors {#Model-constructors}

The `LinearMixedModel` type represents a linear mixed-effects model. Typically it is constructed from a `Formula` and an appropriate `Table` type, usually a `DataFrame`.

## Examples of linear mixed-effects model fits {#Examples-of-linear-mixed-effects-model-fits}

For illustration, several data sets from the _lme4_ package for _R_ are made available in `.arrow` format in this package. Often, for convenience, we will convert these to `DataFrame`s. These data sets include the `dyestuff` and `dyestuff2` data sets.

```julia
using DataFrames, MixedModels, StatsModels
dyestuff = MixedModels.dataset(:dyestuff)
```


```
Arrow.Table with 30 rows, 2 columns, and schema:
 :batch  String
 :yield  Int16
```


```julia
describe(DataFrame(dyestuff))
```

<div v-html="`&lt;div&gt;&lt;div style = &quot;float: left;&quot;&gt;&lt;span&gt;2×7 DataFrame&lt;/span&gt;&lt;/div&gt;&lt;div style = &quot;clear: both;&quot;&gt;&lt;/div&gt;&lt;/div&gt;&lt;div class = &quot;data-frame&quot; style = &quot;overflow-x: scroll;&quot;&gt;&lt;table class = &quot;data-frame&quot; style = &quot;margin-bottom: 6px;&quot;&gt;&lt;thead&gt;&lt;tr class = &quot;header&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;Row&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;variable&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;mean&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;min&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;median&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;max&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;nmissing&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;eltype&lt;/th&gt;&lt;/tr&gt;&lt;tr class = &quot;subheader headerLastRow&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;&lt;/th&gt;&lt;th title = &quot;Symbol&quot; style = &quot;text-align: left;&quot;&gt;Symbol&lt;/th&gt;&lt;th title = &quot;Union{Nothing, Float64}&quot; style = &quot;text-align: left;&quot;&gt;Union…&lt;/th&gt;&lt;th title = &quot;Any&quot; style = &quot;text-align: left;&quot;&gt;Any&lt;/th&gt;&lt;th title = &quot;Union{Nothing, Float64}&quot; style = &quot;text-align: left;&quot;&gt;Union…&lt;/th&gt;&lt;th title = &quot;Any&quot; style = &quot;text-align: left;&quot;&gt;Any&lt;/th&gt;&lt;th title = &quot;Int64&quot; style = &quot;text-align: left;&quot;&gt;Int64&lt;/th&gt;&lt;th title = &quot;DataType&quot; style = &quot;text-align: left;&quot;&gt;DataType&lt;/th&gt;&lt;/tr&gt;&lt;/thead&gt;&lt;tbody&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;1&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;batch&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;A&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;F&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;String&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;2&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;yield&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;1527.5&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;1440&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;1530.0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;1635&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;Int16&lt;/td&gt;&lt;/tr&gt;&lt;/tbody&gt;&lt;/table&gt;&lt;/div&gt;`"></div>

### The `@formula` language in Julia {#The-@formula-language-in-Julia}

MixedModels.jl builds on the the _Julia_ formula language provided by [StatsModels.jl](https://juliastats.org/StatsModels.jl/stable/formula/), which is similar to the formula language in _R_ and is also based on the notation from Wilkinson and Rogers ([1973](https://dx.doi.org/10.2307/2346786)). There are two ways to construct a formula in Julia.  The first way is to enclose the formula expression in the `@formula` macro:
<details class='jldocstring custom-block' open>
<summary><a id='StatsModels.@formula' href='#StatsModels.@formula'><span class="jlbinding">StatsModels.@formula</span></a> <Badge type="info" class="jlObjectType jlMacro" text="Macro" /></summary>



```julia
@formula(ex)
```


Capture and parse a formula expression as a `Formula` struct.

A formula is an abstract specification of a dependence between _left-hand_ and _right-hand_ side variables as in, e.g., a regression model.  Each side specifies at a high level how tabular data is to be converted to a numerical matrix suitable for modeling.  This specification looks something like Julia code, is represented as a Julia `Expr`, but uses special syntax.  The `@formula` macro takes an expression like `y ~ 1 + a*b`, transforms it according to the formula syntax rules into a lowered form (like `y ~ 1 + a + b + a&b`), and constructs a `Formula` struct which captures the original expression, the lowered expression, and the left- and right-hand-side.

Operators that have special interpretations in this syntax are
- `~` is the formula separator, where it is a binary operator (the first argument is the left-hand side, and the second is the right-hand side.
  
- `+` concatenates variables as columns when generating a model matrix.
  
- `&` represents an _interaction_ between two or more variables, which corresponds to a row-wise kronecker product of the individual terms (or element-wise product if all terms involved are continuous/scalar).
  
- `*` expands to all main effects and interactions: `a*b` is equivalent to `a+b+a&b`, `a*b*c` to `a+b+c+a&b+a&c+b&c+a&b&c`, etc.
  
- `1`, `0`, and `-1` indicate the presence (for `1`) or absence (for `0` and `-1`) of an intercept column.
  

The rules that are applied are
- The associative rule (un-nests nested calls to `+`, `&`, and `*`).
  
- The distributive rule (interactions `&` distribute over concatenation `+`).
  
- The `*` rule expands `a*b` to `a+b+a&b` (recursively).
  
- Subtraction is converted to addition and negation, so `x-1` becomes `x + -1` (applies only to subtraction of literal 1).
  
- Single-argument `&` calls are stripped, so `&(x)` becomes the main effect `x`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsModels.jl/blob/v0.7.4/src/formula.jl#L23-L59" target="_blank" rel="noreferrer">source</a></Badge>

</details>


The second way is to combine `Term`s with operators like `+`, `&`, `~`, and others at &quot;run time&quot;.  This is especially useful if you wish to create a formula from a list a variable names.  For instance, the following are equivalent:

```julia
@formula(y ~ 1 + a + b + a & b) == (term(:y) ~ term(1) + term(:a) + term(:b) + term(:a) & term(:b))
```


```
true
```


MixedModels.jl provides additional formula syntax for representing _random-effects terms_.  Most importantly, `|` separates random effects and their grouping factors (as in the formula extension used by the _R_ package [`lme4`](https://cran.r-project.org/web/packages/lme4/index.html).  Much like with the base formula language, `|` can be used within the `@formula` macro and to construct a formula programmatically:

```julia
@formula(y ~ 1 + a + b + (1 + a + b | g))
```


```
FormulaTerm
Response:
  y(unknown)
Predictors:
  1
  a(unknown)
  b(unknown)
  (a,b,g)->(1 + a + b) | g
```


```julia
terms = sum(term(t) for t in [1, :a, :b])
group = term(:g)
response = term(:y)
response ~ terms + (terms | group)
```


```
FormulaTerm
Response:
  y(unknown)
Predictors:
  1
  a(unknown)
  b(unknown)
  (1 + a + b | g)
```


### Models with simple, scalar random effects {#Models-with-simple,-scalar-random-effects}

A basic model with simple, scalar random effects for the levels of `batch` (the batch of an intermediate product, in this case) is declared and fit as

```julia
fm = @formula(yield ~ 1 + (1|batch))
fm1 = fit(MixedModel, fm, dyestuff)
```


```
Linear mixed model fit by maximum likelihood
 yield ~ 1 + (1 | batch)
   logLik   -2 logLik     AIC       AICc        BIC    
  -163.6635   327.3271   333.3271   334.2501   337.5307

Variance components:
            Column    Variance Std.Dev.
batch    (Intercept)  1388.3332 37.2603
Residual              2451.2501 49.5101
 Number of obs: 30; levels of grouping factors: 6

  Fixed-effects parameters:
────────────────────────────────────────────────
              Coef.  Std. Error      z  Pr(>|z|)
────────────────────────────────────────────────
(Intercept)  1527.5     17.6946  86.33    <1e-99
────────────────────────────────────────────────
```


(If you are new to Julia you may find that this first fit takes an unexpectedly long time, due to Just-In-Time (JIT) compilation of the code. The subsequent calls to such functions are much faster.)

```julia
using BenchmarkTools
dyestuff2 = MixedModels.dataset(:dyestuff2)
@benchmark fit(MixedModel, $fm, $dyestuff2)
```


```
BenchmarkTools.Trial: 10000 samples with 1 evaluation per sample.
 Range (min … max):  52.291 μs … 184.542 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     57.625 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   57.469 μs ±   3.115 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

        ▁▁ ▃▃▂▁      ▁▃▂▅▄▇▃█▆▇▂▁                               
  ▁▂▃▅▆█████████▇▆▆▅▇█████████████▇▅▄▃▃▃▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▄
  52.3 μs         Histogram: frequency by time         66.7 μs <

 Memory estimate: 49.91 KiB, allocs estimate: 1073.
```


By default, the model is fit by maximum likelihood. To use the `REML` criterion instead, add the optional named argument `REML=true` to the call to `fit`

```julia
fm1reml = fit(MixedModel, fm, dyestuff, REML=true)
```


```
Linear mixed model fit by REML
 yield ~ 1 + (1 | batch)
 REML criterion at convergence: 319.65427684225943

Variance components:
            Column    Variance Std.Dev.
batch    (Intercept)  1764.0537 42.0006
Residual              2451.2492 49.5101
 Number of obs: 30; levels of grouping factors: 6

  Fixed-effects parameters:
────────────────────────────────────────────────
              Coef.  Std. Error      z  Pr(>|z|)
────────────────────────────────────────────────
(Intercept)  1527.5     19.3834  78.80    <1e-99
────────────────────────────────────────────────
```


### Floating-point type in the model {#Floating-point-type-in-the-model}

The type of `fm1`

```julia
typeof(fm1)
```


```
LinearMixedModel{Float64}
```


includes the floating point type used internally for the various matrices, vectors, and scalars that represent the model. At present, this will always be `Float64` because the parameter estimates are optimized using the [`NLopt` package](https://github.com/JuliaOpt/NLopt.jl) which calls compiled C code that only allows for optimization with respect to a `Float64` parameter vector.

So in theory other floating point types, such as `BigFloat` or `Float32`, can be used to define a model but in practice only `Float64` works at present.
> 
> In theory, theory and practice are the same.  In practice, they aren&#39;t.  – Anon
> 


### Simple, scalar random effects {#Simple,-scalar-random-effects}

A simple, scalar random effects term in a mixed-effects model formula is of the form `(1|G)`. All random effects terms end with `|G` where `G` is the _grouping factor_ for the random effect. The name or, more generally the expression, `G`, should evaluate to a categorical array that has a distinct set of _levels_. The random effects are associated with the levels of the grouping factor.

A _scalar_ random effect is, as the name implies, one scalar value for each level of the grouping factor. A _simple, scalar_ random effects term is of the form, `(1|G)`. It corresponds to a shift in the intercept for each level of the grouping factor.

### Models with vector-valued random effects {#Models-with-vector-valued-random-effects}

The _sleepstudy_ data are observations of reaction time, `reaction`, on several subjects, `subj`, after 0 to 9 days of sleep deprivation, `days`. A model with random intercepts and random slopes for each subject, allowing for within-subject correlation of the slope and intercept, is fit as

```julia
sleepstudy = MixedModels.dataset(:sleepstudy)
fm2 = fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy)
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -875.9697  1751.9393  1763.9393  1764.4249  1783.0971

Variance components:
            Column    Variance Std.Dev.   Corr.
subj     (Intercept)  565.51066 23.78047
         days          32.68212  5.71683 +0.08
Residual              654.94145 25.59182
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
──────────────────────────────────────────────────
                Coef.  Std. Error      z  Pr(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405      6.63226  37.91    <1e-99
days          10.4673     1.50224   6.97    <1e-11
──────────────────────────────────────────────────
```


### Models with multiple, scalar random-effects terms {#Models-with-multiple,-scalar-random-effects-terms}

A model for the _Penicillin_ data incorporates random effects for the plate, and for the sample. As every sample is used on every plate these two factors are _crossed_.

```julia
penicillin = MixedModels.dataset(:penicillin)
fm3 = fit(MixedModel, @formula(diameter ~ 1 + (1|plate) + (1|sample)), penicillin)
```


```
Linear mixed model fit by maximum likelihood
 diameter ~ 1 + (1 | plate) + (1 | sample)
   logLik   -2 logLik     AIC       AICc        BIC    
  -166.0942   332.1883   340.1883   340.4761   352.0676

Variance components:
            Column   Variance Std.Dev. 
plate    (Intercept)  0.714979 0.845565
sample   (Intercept)  3.135193 1.770648
Residual              0.302426 0.549933
 Number of obs: 144; levels of grouping factors: 24, 6

  Fixed-effects parameters:
─────────────────────────────────────────────────
               Coef.  Std. Error      z  Pr(>|z|)
─────────────────────────────────────────────────
(Intercept)  22.9722    0.744596  30.85    <1e-99
─────────────────────────────────────────────────
```


In contrast, the `cask` grouping factor is _nested_ within the `batch` grouping factor in the _Pastes_ data.

```julia
pastes = DataFrame(MixedModels.dataset(:pastes))
describe(pastes)
```

<div v-html="`&lt;div&gt;&lt;div style = &quot;float: left;&quot;&gt;&lt;span&gt;3×7 DataFrame&lt;/span&gt;&lt;/div&gt;&lt;div style = &quot;clear: both;&quot;&gt;&lt;/div&gt;&lt;/div&gt;&lt;div class = &quot;data-frame&quot; style = &quot;overflow-x: scroll;&quot;&gt;&lt;table class = &quot;data-frame&quot; style = &quot;margin-bottom: 6px;&quot;&gt;&lt;thead&gt;&lt;tr class = &quot;header&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;Row&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;variable&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;mean&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;min&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;median&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;max&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;nmissing&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;eltype&lt;/th&gt;&lt;/tr&gt;&lt;tr class = &quot;subheader headerLastRow&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;&lt;/th&gt;&lt;th title = &quot;Symbol&quot; style = &quot;text-align: left;&quot;&gt;Symbol&lt;/th&gt;&lt;th title = &quot;Union{Nothing, Float64}&quot; style = &quot;text-align: left;&quot;&gt;Union…&lt;/th&gt;&lt;th title = &quot;Any&quot; style = &quot;text-align: left;&quot;&gt;Any&lt;/th&gt;&lt;th title = &quot;Union{Nothing, Float64}&quot; style = &quot;text-align: left;&quot;&gt;Union…&lt;/th&gt;&lt;th title = &quot;Any&quot; style = &quot;text-align: left;&quot;&gt;Any&lt;/th&gt;&lt;th title = &quot;Int64&quot; style = &quot;text-align: left;&quot;&gt;Int64&lt;/th&gt;&lt;th title = &quot;DataType&quot; style = &quot;text-align: left;&quot;&gt;DataType&lt;/th&gt;&lt;/tr&gt;&lt;/thead&gt;&lt;tbody&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;1&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;batch&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;A&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;J&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;String&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;2&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;cask&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;a&lt;/td&gt;&lt;td style = &quot;font-style: italic; text-align: left;&quot;&gt;&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;c&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;String&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;3&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;strength&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;60.0533&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;54.2&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;59.3&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;66.0&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;Float64&lt;/td&gt;&lt;/tr&gt;&lt;/tbody&gt;&lt;/table&gt;&lt;/div&gt;`"></div>

This can be expressed using the solidus (the &quot;`/`&quot; character) to separate grouping factors, read &quot;`cask` nested within `batch`&quot;:

```julia
fm4a = fit(MixedModel, @formula(strength ~ 1 + (1|batch/cask)), pastes)
```


```
Linear mixed model fit by maximum likelihood
 strength ~ 1 + (1 | batch) + (1 | batch & cask)
   logLik   -2 logLik     AIC       AICc        BIC    
  -123.9972   247.9945   255.9945   256.7217   264.3718

Variance components:
                Column   Variance Std.Dev. 
batch & cask (Intercept)  8.433617 2.904069
batch        (Intercept)  1.199180 1.095071
Residual                  0.678002 0.823409
 Number of obs: 60; levels of grouping factors: 30, 10

  Fixed-effects parameters:
─────────────────────────────────────────────────
               Coef.  Std. Error      z  Pr(>|z|)
─────────────────────────────────────────────────
(Intercept)  60.0533    0.642136  93.52    <1e-99
─────────────────────────────────────────────────
```


If the levels of the inner grouping factor are unique across the levels of the outer grouping factor, then this nesting does not need to expressed explicitly in the model syntax. For example, defining `sample` to be the combination of `batch` and `cask`, yields a naming scheme where the nesting is apparent from the data even if not expressed in the formula. (That is, each level of `sample` occurs in conjunction with only one level of `batch`.) As such, this model is equivalent to the previous one.

```julia
pastes.sample = (string.(pastes.cask, "&",  pastes.batch))
fm4b = fit(MixedModel, @formula(strength ~ 1 + (1|sample) + (1|batch)), pastes)
```


```
Linear mixed model fit by maximum likelihood
 strength ~ 1 + (1 | sample) + (1 | batch)
   logLik   -2 logLik     AIC       AICc        BIC    
  -123.9972   247.9945   255.9945   256.7217   264.3718

Variance components:
            Column   Variance Std.Dev. 
sample   (Intercept)  8.433617 2.904069
batch    (Intercept)  1.199178 1.095070
Residual              0.678002 0.823409
 Number of obs: 60; levels of grouping factors: 30, 10

  Fixed-effects parameters:
─────────────────────────────────────────────────
               Coef.  Std. Error      z  Pr(>|z|)
─────────────────────────────────────────────────
(Intercept)  60.0533    0.642136  93.52    <1e-99
─────────────────────────────────────────────────
```


In observational studies it is common to encounter _partially crossed_ grouping factors. For example, the _InstEval_ data are course evaluations by students, `s`, of instructors, `d`. Additional covariates include the academic department, `dept`, in which the course was given and `service`, whether or not it was a service course.

```julia
insteval = MixedModels.dataset(:insteval)
fm5 = fit(MixedModel, @formula(y ~ 1 + service * dept + (1|s) + (1|d)), insteval)
```


```
Linear mixed model fit by maximum likelihood
 y ~ 1 + service + dept + service & dept + (1 | s) + (1 | d)
    logLik     -2 logLik       AIC         AICc          BIC     
 -118792.7767  237585.5534  237647.5534  237647.5804  237932.8763

Variance components:
            Column   Variance Std.Dev. 
s        (Intercept)  0.105418 0.324681
d        (Intercept)  0.258416 0.508347
Residual              1.384728 1.176745
 Number of obs: 73421; levels of grouping factors: 2972, 1128

  Fixed-effects parameters:
────────────────────────────────────────────────────────────────
                              Coef.  Std. Error      z  Pr(>|z|)
────────────────────────────────────────────────────────────────
(Intercept)              3.27628      0.0793647  41.28    <1e-99
service: Y               0.0116044    0.0699321   0.17    0.8682
dept: D02               -0.0411091    0.120331   -0.34    0.7326
dept: D03                0.00967413   0.108411    0.09    0.9289
dept: D04                0.105017     0.0944964   1.11    0.2664
dept: D05                0.0828643    0.11148     0.74    0.4573
dept: D06               -0.01194      0.0978342  -0.12    0.9029
dept: D07                0.0992679    0.110598    0.90    0.3694
dept: D08                0.0575337    0.127935    0.45    0.6529
dept: D09               -0.00263181   0.107085   -0.02    0.9804
dept: D10               -0.223423     0.099838   -2.24    0.0252
dept: D11                0.0129816    0.110639    0.12    0.9066
dept: D12                0.00503825   0.0944243   0.05    0.9574
dept: D14                0.0050827    0.109041    0.05    0.9628
dept: D15               -0.0466719    0.101942   -0.46    0.6471
service: Y & dept: D02  -0.144352     0.0929527  -1.55    0.1204
service: Y & dept: D03   0.0174078    0.0927237   0.19    0.8511
service: Y & dept: D04  -0.0381262    0.0810901  -0.47    0.6382
service: Y & dept: D05   0.0596632    0.123952    0.48    0.6303
service: Y & dept: D06  -0.254044     0.080781   -3.14    0.0017
service: Y & dept: D07  -0.151634     0.11157    -1.36    0.1741
service: Y & dept: D08   0.0508942    0.112189    0.45    0.6501
service: Y & dept: D09  -0.259448     0.0899448  -2.88    0.0039
service: Y & dept: D10   0.25907      0.111369    2.33    0.0200
service: Y & dept: D11  -0.276577     0.0819621  -3.37    0.0007
service: Y & dept: D12  -0.0418489    0.0792928  -0.53    0.5977
service: Y & dept: D14  -0.256742     0.0931016  -2.76    0.0058
service: Y & dept: D15   0.24042      0.0982071   2.45    0.0144
────────────────────────────────────────────────────────────────
```


### Simplifying the random effect correlation structure {#Simplifying-the-random-effect-correlation-structure}

MixedModels.jl estimates not only the _variance_ of the effects for each random effect level, but also the _correlation_ between the random effects for different predictors. So, for the model of the _sleepstudy_ data above, one of the parameters that is estimated is the correlation between each subject&#39;s random intercept (i.e., their baseline reaction time) and slope (i.e., their particular change in reaction time per day of sleep deprivation). In some cases, you may wish to simplify the random effects structure by removing these correlation parameters. This often arises when there are many random effects you want to estimate (as is common in psychological experiments with many conditions and covariates), since the number of random effects parameters increases as the square of the number of predictors, making these models difficult to estimate from limited data.

The special syntax `zerocorr` can be applied to individual random effects terms inside the `@formula`:

```julia
fm2zerocorr_fm = fit(MixedModel, @formula(reaction ~ 1 + days + zerocorr(1 + days|subj)), sleepstudy)
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + zerocorr(1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -876.0016  1752.0033  1762.0033  1762.3481  1777.9680

Variance components:
            Column    Variance Std.Dev.   Corr.
subj     (Intercept)  584.25897 24.17145
         days          33.63281  5.79938   .  
Residual              653.11578 25.55613
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
──────────────────────────────────────────────────
                Coef.  Std. Error      z  Pr(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405      6.70771  37.48    <1e-99
days          10.4673     1.51931   6.89    <1e-11
──────────────────────────────────────────────────
```


Alternatively, correlations between parameters can be removed by including them as separate random effects terms:

```julia
fit(MixedModel, @formula(reaction ~ 1 + days + (1|subj) + (days|subj)), sleepstudy)
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj) + (days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -876.0016  1752.0033  1762.0033  1762.3481  1777.9680

Variance components:
            Column    Variance Std.Dev.   Corr.
subj     (Intercept)  584.25897 24.17145
         days          33.63281  5.79938   .  
Residual              653.11578 25.55613
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
──────────────────────────────────────────────────
                Coef.  Std. Error      z  Pr(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405      6.70771  37.48    <1e-99
days          10.4673     1.51931   6.89    <1e-11
──────────────────────────────────────────────────
```


Finally, for predictors that are categorical, MixedModels.jl will estimate correlations between each level. Notice the large number of correlation parameters if we treat `days` as a categorical variable by giving it contrasts:

```julia
fit(MixedModel, @formula(reaction ~ 1 + days + (1 + days|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -805.3994  1610.7987  1742.7987  1821.0642  1953.5339

Variance components:
            Column    Variance  Std.Dev.   Corr.
subj     (Intercept)   954.66067 30.89758
         days: 1       498.34749 22.32370 -0.30
         days: 2       916.26827 30.26992 -0.57 +0.75
         days: 3      1268.96763 35.62257 -0.37 +0.72 +0.87
         days: 4      1487.09850 38.56292 -0.31 +0.59 +0.67 +0.91
         days: 5      2302.20224 47.98127 -0.25 +0.46 +0.45 +0.70 +0.85
         days: 6      3851.78092 62.06272 -0.27 +0.30 +0.48 +0.70 +0.77 +0.75
         days: 7      1808.39158 42.52519 -0.16 +0.22 +0.47 +0.50 +0.63 +0.64 +0.71
         days: 8      3159.98677 56.21376 -0.20 +0.29 +0.36 +0.56 +0.73 +0.90 +0.73 +0.74
         days: 9      3086.04128 55.55215 +0.06 +0.26 +0.16 +0.38 +0.59 +0.78 +0.38 +0.53 +0.85
Residual                18.79114  4.33488
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652       7.35396  34.90    <1e-99
days: 1        7.84395     5.45654   1.44    0.1506
days: 2        8.71009     7.27954   1.20    0.2315
days: 3       26.3402      8.51975   3.09    0.0020
days: 4       31.9976      9.2035    3.48    0.0005
days: 5       51.8666     11.4012    4.55    <1e-05
days: 6       55.5265     14.6995    3.78    0.0002
days: 7       62.0988     10.1269    6.13    <1e-09
days: 8       79.9777     13.3283    6.00    <1e-08
days: 9       94.1994     13.1733    7.15    <1e-12
───────────────────────────────────────────────────
```


Separating the `1` and `days` random effects into separate terms removes the correlations between the intercept and the levels of `days`, but not between the levels themselves:

```julia
fit(MixedModel, @formula(reaction ~ 1 + days + (1|subj) + (days|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj) + (days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -827.7749  1655.5498  1769.5498  1823.7465  1951.5483

Variance components:
            Column    Variance Std.Dev.  Corr.
subj     (Intercept)   789.6792 28.1012
         days: 1         0.0000  0.0000   .  
         days: 2       357.7766 18.9150   .    +NaN
         days: 3       684.7231 26.1672   .    +NaN +0.81
         days: 4       949.6394 30.8162   .    +NaN +0.57 +0.91
         days: 5      1751.2433 41.8479   .    +NaN +0.26 +0.66 +0.87
         days: 6      3357.4474 57.9435   .    +NaN +0.45 +0.72 +0.80 +0.76
         days: 7      1539.1760 39.2323   .    +NaN +0.39 +0.42 +0.59 +0.62 +0.71
         days: 8      2735.1879 52.2990   .    +NaN +0.22 +0.52 +0.75 +0.93 +0.72 +0.75
         days: 9      2767.1610 52.6038   .    +NaN -0.05 +0.28 +0.57 +0.80 +0.34 +0.52 +0.87
Residual               135.0061 11.6192
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652       7.16738  35.81    <1e-99
days: 1        7.84395     3.87307   2.03    0.0428
days: 2        8.71009     5.90569   1.47    0.1402
days: 3       26.3402      7.28291   3.62    0.0003
days: 4       31.9976      8.23155   3.89    0.0001
days: 5       51.8667     10.5968    4.89    <1e-06
days: 6       55.5265     14.196     3.91    <1e-04
days: 7       62.0988     10.0255    6.19    <1e-09
days: 8       79.9777     12.9211    6.19    <1e-09
days: 9       94.1994     12.9897    7.25    <1e-12
───────────────────────────────────────────────────
```


(Notice that the variance component for `days: 1` is estimated as zero, so the correlations for this component are undefined and expressed as `NaN`, not a number.)

An alternative is to force all the levels of `days` as indicators using `fulldummy` encoding.
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fulldummy' href='#MixedModels.fulldummy'><span class="jlbinding">MixedModels.fulldummy</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fulldummy(term::CategoricalTerm)
```


Assign &quot;contrasts&quot; that include all indicator columns (dummy variables) and an intercept column.

This will result in an under-determined set of contrasts, which is not a problem in the random effects because of the regularization, or &quot;shrinkage&quot;, of the conditional modes.

The interaction of `fulldummy` with complex random effects is subtle and complex with numerous potential edge cases. As we discover these edge cases, we will document and determine their behavior. Until such time, please check the model summary to verify that the expansion is working as you expected. If it is not, please report a use case on GitHub.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/randomeffectsterm.jl#L195-L207" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
fit(MixedModel, @formula(reaction ~ 1 + days + (1 + fulldummy(days)|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -805.3991  1610.7982  1764.7982  1882.5630  2010.6559

Variance components:
            Column    Variance  Std.Dev.   Corr.
subj     (Intercept)   727.90862 26.97978
         days: 0       499.02116 22.33878 -0.22
         days: 1       381.88986 19.54200 -0.07 +0.44
         days: 2       463.96481 21.53984 -0.34 +0.05 +0.53
         days: 3       379.47587 19.48014 +0.28 -0.45 +0.18 +0.63
         days: 4       327.96770 18.10988 +0.65 -0.82 -0.38 -0.07 +0.66
         days: 5      1051.33792 32.42434 +0.42 -0.52 -0.30 -0.27 +0.20 +0.62
         days: 6      2321.59696 48.18295 +0.27 -0.48 -0.44 -0.08 +0.36 +0.58 +0.54
         days: 7      1134.97511 33.68939 +0.27 -0.12 -0.28 +0.06 -0.03 +0.18 +0.32 +0.48
         days: 8      1870.82916 43.25308 +0.35 -0.41 -0.41 -0.27 +0.04 +0.45 +0.82 +0.55 +0.55
         days: 9      2277.31514 47.72122 +0.47 -0.14 -0.15 -0.34 -0.10 +0.32 +0.67 +0.07 +0.32 +0.78
Residual                19.21115  4.38305
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652       7.36068  34.87    <1e-99
days: 1        7.84395     5.4558    1.44    0.1505
days: 2        8.71009     7.28336   1.20    0.2317
days: 3       26.3402      8.52163   3.09    0.0020
days: 4       31.9976      9.20409   3.48    0.0005
days: 5       51.8667     11.3962    4.55    <1e-05
days: 6       55.5265     14.7077    3.78    0.0002
days: 7       62.0988     10.1303    6.13    <1e-09
days: 8       79.9777     13.3249    6.00    <1e-08
days: 9       94.1994     13.1611    7.16    <1e-12
───────────────────────────────────────────────────
```


This fit produces a better fit as measured by the objective (negative twice the log-likelihood is 1610.8) but at the expense of adding many more parameters to the model. As a result, model comparison criteria such, as `AIC` and `BIC`, are inflated.

But using `zerocorr` on the individual terms does remove the correlations between the levels:

```julia
fit(MixedModel, @formula(reaction ~ 1 + days + zerocorr(1 + days|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + zerocorr(1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -882.9138  1765.8276  1807.8276  1813.6757  1874.8797

Variance components:
            Column    Variance Std.Dev.  Corr.
subj     (Intercept)   958.4877 30.9595
         days: 1         0.0000  0.0000   .  
         days: 2         0.0000  0.0000   .     .  
         days: 3         0.0000  0.0000   .     .     .  
         days: 4         0.0000  0.0000   .     .     .     .  
         days: 5       519.6167 22.7951   .     .     .     .     .  
         days: 6      1704.1051 41.2808   .     .     .     .     .     .  
         days: 7       608.7869 24.6736   .     .     .     .     .     .     .  
         days: 8      1273.0228 35.6794   .     .     .     .     .     .     .     .  
         days: 9      1753.8569 41.8791   .     .     .     .     .     .     .     .     .  
Residual               434.8673 20.8535
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652       8.79822  29.17    <1e-99
days: 1        7.84395     6.95116   1.13    0.2591
days: 2        8.71009     6.95116   1.25    0.2102
days: 3       26.3402      6.95116   3.79    0.0002
days: 4       31.9976      6.95116   4.60    <1e-05
days: 5       51.8667      8.78557   5.90    <1e-08
days: 6       55.5265     11.9579    4.64    <1e-05
days: 7       62.0988      9.06312   6.85    <1e-11
days: 8       79.9777     10.9106    7.33    <1e-12
days: 9       94.1994     12.0729    7.80    <1e-14
───────────────────────────────────────────────────
```


```julia
fit(MixedModel, @formula(reaction ~ 1 + days + (1|subj) + zerocorr(days|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + (1 | subj) + zerocorr(days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -882.9138  1765.8276  1807.8276  1813.6757  1874.8797

Variance components:
            Column    Variance Std.Dev.  Corr.
subj     (Intercept)   958.4877 30.9595
         days: 1         0.0000  0.0000   .  
         days: 2         0.0000  0.0000   .     .  
         days: 3         0.0000  0.0000   .     .     .  
         days: 4         0.0000  0.0000   .     .     .     .  
         days: 5       519.6167 22.7951   .     .     .     .     .  
         days: 6      1704.1051 41.2808   .     .     .     .     .     .  
         days: 7       608.7869 24.6736   .     .     .     .     .     .     .  
         days: 8      1273.0228 35.6794   .     .     .     .     .     .     .     .  
         days: 9      1753.8569 41.8791   .     .     .     .     .     .     .     .     .  
Residual               434.8673 20.8535
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652       8.79822  29.17    <1e-99
days: 1        7.84395     6.95116   1.13    0.2591
days: 2        8.71009     6.95116   1.25    0.2102
days: 3       26.3402      6.95116   3.79    0.0002
days: 4       31.9976      6.95116   4.60    <1e-05
days: 5       51.8667      8.78557   5.90    <1e-08
days: 6       55.5265     11.9579    4.64    <1e-05
days: 7       62.0988      9.06312   6.85    <1e-11
days: 8       79.9777     10.9106    7.33    <1e-12
days: 9       94.1994     12.0729    7.80    <1e-14
───────────────────────────────────────────────────
```


```julia
fit(MixedModel, @formula(reaction ~ 1 + days + zerocorr(1 + fulldummy(days)|subj)), sleepstudy,
    contrasts = Dict(:days => DummyCoding()))
```


```
Linear mixed model fit by maximum likelihood
 reaction ~ 1 + days + zerocorr(1 + days | subj)
   logLik   -2 logLik     AIC       AICc        BIC    
  -878.9843  1757.9686  1801.9686  1808.4145  1872.2137

Variance components:
            Column    Variance  Std.Dev.   Corr.
subj     (Intercept)  1135.22369 33.69308
         days: 0       775.53471 27.84842   .  
         days: 1       357.33070 18.90319   .     .  
         days: 2       221.08627 14.86897   .     .     .  
         days: 3         0.00000  0.00000   .     .     .     .  
         days: 4        44.25126  6.65216   .     .     .     .     .  
         days: 5       670.03272 25.88499   .     .     .     .     .     .  
         days: 6      1739.38985 41.70599   .     .     .     .     .     .     .  
         days: 7       908.48034 30.14101   .     .     .     .     .     .     .     .  
         days: 8      1457.44632 38.17652   .     .     .     .     .     .     .     .     .  
         days: 9      2027.51717 45.02796   .     .     .     .     .     .     .     .     .     .  
Residual               181.13013 13.45846
 Number of obs: 180; levels of grouping factors: 18

  Fixed-effects parameters:
───────────────────────────────────────────────────
                 Coef.  Std. Error      z  Pr(>|z|)
───────────────────────────────────────────────────
(Intercept)  256.652      10.7804   23.81    <1e-99
days: 1        7.84395     9.11387   0.86    0.3894
days: 2        8.71009     8.68869   1.00    0.3161
days: 3       26.3402      7.95052   3.31    0.0009
days: 4       31.9976      8.10366   3.95    <1e-04
days: 5       51.8666     10.0217    5.18    <1e-06
days: 6       55.5264     12.6429    4.39    <1e-04
days: 7       62.0988     10.6622    5.82    <1e-08
days: 8       79.9777     12.0075    6.66    <1e-10
days: 9       94.1994     13.2609    7.10    <1e-11
───────────────────────────────────────────────────
```


## Fitting generalized linear mixed models {#Fitting-generalized-linear-mixed-models}

To create a GLMM representation, the distribution family for the response, and possibly the link function, must be specified.

```julia
verbagg = MixedModels.dataset(:verbagg)
verbaggform = @formula(r2 ~ 1 + anger + gender + btype + situ + mode + (1|subj) + (1|item));
gm1 = fit(MixedModel, verbaggform, verbagg, Bernoulli())
```


```
Generalized Linear Mixed Model fit by maximum likelihood (nAGQ = 1)
  r2 ~ 1 + anger + gender + btype + situ + mode + (1 | subj) + (1 | item)
  Distribution: Bernoulli{Float64}
  Link: LogitLink()

   logLik    deviance     AIC       AICc        BIC    
 -4067.9164  8135.8329  8153.8329  8153.8566  8216.2370

Variance components:
        Column   Variance Std.Dev. 
subj (Intercept)  1.793488 1.339212
item (Intercept)  0.117149 0.342271

 Number of obs: 7584; levels of grouping factors: 316, 24

Fixed-effects parameters:
──────────────────────────────────────────────────────
                   Coef.  Std. Error       z  Pr(>|z|)
──────────────────────────────────────────────────────
(Intercept)   -0.152977    0.385222    -0.40    0.6913
anger          0.0573961   0.0167528    3.43    0.0006
gender: M      0.320783    0.191207     1.68    0.0934
btype: scold  -1.05988     0.184159    -5.76    <1e-08
btype: shout  -2.10389     0.186517   -11.28    <1e-28
situ: self    -1.05438     0.151194    -6.97    <1e-11
mode: want     0.706981    0.151007     4.68    <1e-05
──────────────────────────────────────────────────────
```


The canonical link, which is `LogitLink` for the `Bernoulli` distribution, is used if no explicit link is specified.

Note that, in keeping with convention in the [`GLM` package](https://github.com/JuliaStats/GLM.jl), the distribution family for a binary (i.e. 0/1) response is the `Bernoulli` distribution. The `Binomial` distribution is only used when the response is the fraction of trials returning a positive, in which case the number of trials must be specified as the case weights.

### Optional arguments to fit {#Optional-arguments-to-fit}

An alternative approach is to create the `GeneralizedLinearMixedModel` object then call `fit!` on it. The optional arguments `fast` and/or `nAGQ` can be passed to the optimization process via both `fit` and `fit!` (i.e these optimization settings are not used nor recognized when constructing the model).

As the name implies, `fast=true`, provides a faster but somewhat less accurate fit. These fits may suffice for model comparisons.

```julia
gm1a = fit(MixedModel, verbaggform, verbagg, Bernoulli(), fast = true)
deviance(gm1a) - deviance(gm1)
```


```
0.33801335722546355
```


```julia
@benchmark fit(MixedModel, $verbaggform, $verbagg, Bernoulli())
```


```
BenchmarkTools.Trial: 4 samples with 1 evaluation per sample.
 Range (min … max):  1.380 s …  1.384 s  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.382 s             ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.382 s ± 2.029 ms  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █                   █           █                      █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.38 s        Histogram: frequency by time        1.38 s <

 Memory estimate: 29.89 MiB, allocs estimate: 472418.
```


```julia
@benchmark fit(MixedModel, $verbaggform, $verbagg, Bernoulli(), fast = true)
```


```
BenchmarkTools.Trial: 53 samples with 1 evaluation per sample.
 Range (min … max):  92.660 ms … 137.961 ms  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     93.585 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   95.931 ms ±   8.396 ms  ┊ GC (mean ± σ):  0.12% ± 0.82%

   █▂                                                           
  ▅██▇▇▁▁▁▅▁▁▁▁▁▅▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▅▁▁▁▁▁▅ ▁
  92.7 ms       Histogram: log(frequency) by time       126 ms <

 Memory estimate: 9.54 MiB, allocs estimate: 82996.
```


The optional argument `nAGQ=k` causes evaluation of the deviance function to use a `k` point adaptive Gauss-Hermite quadrature rule. This method only applies to models with a single, simple, scalar random-effects term, such as

```julia
contraception = MixedModels.dataset(:contra)
contraform = @formula(use ~ 1 + age + abs2(age) + livch + urban + (1|dist));
bernoulli = Bernoulli()
deviances = Dict{Symbol,Float64}()
b = @benchmarkable deviances[:default] = deviance(fit(MixedModel, $contraform, $contraception, $bernoulli));
run(b)
b = @benchmarkable deviances[:fast] = deviance(fit(MixedModel, $contraform, $contraception, $bernoulli, fast = true));
run(b)
b = @benchmarkable deviances[:nAGQ] = deviance(fit(MixedModel, $contraform, $contraception, $bernoulli, nAGQ=9));
run(b)
b = @benchmarkable deviances[:nAGQ_fast] = deviance(fit(MixedModel, $contraform, $contraception, $bernoulli, nAGQ=9, fast=true));
run(b)
sort(deviances)
```


```
OrderedCollections.OrderedDict{Symbol, Float64} with 4 entries:
  :default   => 2372.73
  :fast      => 2372.78
  :nAGQ      => 2372.46
  :nAGQ_fast => 2372.51
```


# Extractor functions {#Extractor-functions}

`LinearMixedModel` and `GeneralizedLinearMixedModel` are subtypes of `StatsAPI.RegressionModel` which, in turn, is a subtype of `StatsBase.StatisticalModel`. Many of the generic extractors defined in the `StatsBase` package have methods for these models.

## Model-fit statistics {#Model-fit-statistics}

The statistics describing the quality of the model fit include
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.loglikelihood' href='#StatsAPI.loglikelihood'><span class="jlbinding">StatsAPI.loglikelihood</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
loglikelihood(model::StatisticalModel)
loglikelihood(model::StatisticalModel, observation)
```


Return the log-likelihood of the model.

With an `observation` argument, return the contribution of `observation` to the log-likelihood of `model`.

If `observation` is a `Colon`, return a vector of each observation&#39;s contribution to the log-likelihood of the model. In other words, this is the vector of the pointwise log-likelihood contributions.

In general, `sum(loglikehood(model, :)) == loglikelihood(model)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L68-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.aic' href='#StatsAPI.aic'><span class="jlbinding">StatsAPI.aic</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
aic(model::StatisticalModel)
```


Akaike&#39;s Information Criterion, defined as $-2 \log L + 2k$, with $L$ the likelihood of the model, and `k` its number of consumed degrees of freedom (as returned by [`dof`](/constructors#StatsAPI.dof)).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L182-L188" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.bic' href='#StatsAPI.bic'><span class="jlbinding">StatsAPI.bic</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
bic(model::StatisticalModel)
```


Bayesian Information Criterion, defined as $-2 \log L + k \log n$, with $L$ the likelihood of the model,  $k$ its number of consumed degrees of freedom (as returned by [`dof`](/constructors#StatsAPI.dof)), and $n$ the number of observations (as returned by [`nobs`](/constructors#StatsAPI.nobs)).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L205-L212" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.dof' href='#StatsAPI.dof'><span class="jlbinding">StatsAPI.dof</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
dof(model::StatisticalModel)
```


Return the number of degrees of freedom consumed in the model, including when applicable the intercept and the distribution&#39;s dispersion parameter.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L114-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.nobs' href='#StatsAPI.nobs'><span class="jlbinding">StatsAPI.nobs</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
nobs(model::StatisticalModel)
```


Return the number of independent observations on which the model was fitted. Be careful when using this information, as the definition of an independent observation may vary depending on the model, on the format used to pass the data, on the sampling plan (if specified), etc.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L104-L111" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
loglikelihood(fm1)
```


```
-163.6635299405715
```


```julia
aic(fm1)
```


```
333.327059881143
```


```julia
bic(fm1)
```


```
337.5306520261295
```


```julia
dof(fm1)   # 1 fixed effect, 2 variances
```


```
3
```


```julia
nobs(fm1)  # 30 observations
```


```
30
```


```julia
loglikelihood(gm1)
```


```
-4067.9164291743905
```


In general the [`deviance`](https://en.wikipedia.org/wiki/Deviance_(statistics)) of a statistical model fit is negative twice the log-likelihood adjusting for the saturated model.
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.deviance-Tuple{StatisticalModel}' href='#StatsAPI.deviance-Tuple{StatisticalModel}'><span class="jlbinding">StatsAPI.deviance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
deviance(m::GeneralizedLinearMixedModel{T}, nAGQ=1)::T where {T}
```


Return the deviance of `m` evaluated by the Laplace approximation (`nAGQ=1`) or `nAGQ`-point adaptive Gauss-Hermite quadrature.

If the distribution `D` does not have a scale parameter the Laplace approximation is the squared length of the conditional modes, $u$, plus the determinant of $Λ'Z'WZΛ + I$, plus the sum of the squared deviance residuals.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L75-L84" target="_blank" rel="noreferrer">source</a></Badge>

</details>


Because it is not clear what the saturated model corresponding to a particular `LinearMixedModel` should be, negative twice the log-likelihood is called the `objective`.
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.objective' href='#MixedModels.objective'><span class="jlbinding">MixedModels.objective</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
objective(m::LinearMixedModel)
```


Return negative twice the log-likelihood of model `m`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L784-L788" target="_blank" rel="noreferrer">source</a></Badge>

</details>


This value is also accessible as the `deviance` but the user should bear in mind that this doesn&#39;t have all the properties of a deviance which is corrected for the saturated model. For example, it is not necessarily non-negative.

```julia
objective(fm1)
```


```
327.327059881143
```


```julia
deviance(fm1)
```


```
327.327059881143
```


The value optimized when fitting a `GeneralizedLinearMixedModel` is the Laplace approximation to the deviance or an adaptive Gauss-Hermite evaluation.
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.deviance!' href='#MixedModels.deviance!'><span class="jlbinding">MixedModels.deviance!</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
deviance!(m::GeneralizedLinearMixedModel, nAGQ=1)
```


Update `m.η`, `m.μ`, etc., install the working response and working weights in `m.LMM`, update `m.LMM.A` and `m.LMM.R`, then evaluate the `deviance`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L135-L140" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
MixedModels.deviance!(gm1)
```


```
8135.832858348763
```


## Fixed-effects parameter estimates {#Fixed-effects-parameter-estimates}

The `coef` and `fixef` extractors both return the maximum likelihood estimates of the fixed-effects coefficients. They differ in their behavior in the rank-deficient case. The associated `coefnames` and `fixefnames` return the corresponding coefficient names.
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.coef' href='#StatsAPI.coef'><span class="jlbinding">StatsAPI.coef</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
coef(model::StatisticalModel)
```


Return the coefficients of the model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L8-L12" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.coefnames' href='#StatsAPI.coefnames'><span class="jlbinding">StatsAPI.coefnames</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
coefnames(model::StatisticalModel)
```


Return the names of the coefficients.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L15-L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fixef' href='#MixedModels.fixef'><span class="jlbinding">MixedModels.fixef</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fixef(m::MixedModel)
```


Return the fixed-effects parameter vector estimate of `m`.

In the rank-deficient case the truncated parameter vector, of length `rank(m)` is returned. This is unlike `coef` which always returns a vector whose length matches the number of columns in `X`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L551-L559" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.fixefnames' href='#MixedModels.fixefnames'><span class="jlbinding">MixedModels.fixefnames</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
fixefnames(m::MixedModel)
```


Return a (permuted and truncated in the rank-deficient case) vector of coefficient names.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L562-L566" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
coef(fm1)
coefnames(fm1)
```


```
1-element Vector{String}:
 "(Intercept)"
```


```julia
fixef(fm1)
fixefnames(fm1)
```


```
1-element Vector{String}:
 "(Intercept)"
```


An alternative extractor for the fixed-effects coefficient is the `β` property. Properties whose names are Greek letters usually have an alternative spelling, which is the name of the Greek letter.

```julia
fm1.β
```


```
1-element Vector{Float64}:
 1527.4999999999989
```


```julia
fm1.beta
```


```
1-element Vector{Float64}:
 1527.4999999999989
```


```julia
gm1.β
```


```
7-element Vector{Float64}:
 -0.15297658506693748
  0.05739606276708647
  0.3207830298410298
 -1.059879787138412
 -2.103893961040908
 -1.054376691878555
  0.706980989641357
```


A full list of property names is returned by `propertynames`

```julia
propertynames(fm1)
```


```
(:formula, :reterms, :Xymat, :feterm, :sqrtwts, :parmap, :dims, :A, :L, :optsum, :θ, :theta, :β, :beta, :βs, :betas, :λ, :lambda, :stderror, :σ, :sigma, :σs, :sigmas, :σρs, :sigmarhos, :b, :u, :lowerbd, :X, :y, :corr, :vcov, :PCA, :rePCA, :objective, :pvalues)
```


```julia
propertynames(gm1)
```


```
(:A, :L, :theta, :beta, :coef, :λ, :σ, :sigma, :X, :y, :lowerbd, :objective, :σρs, :σs, :corr, :vcov, :PCA, :rePCA, :LMM, :β, :θ, :b, :u, :resp, :wt)
```


The variance-covariance matrix of the fixed-effects coefficients is returned by
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.vcov' href='#StatsAPI.vcov'><span class="jlbinding">StatsAPI.vcov</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
vcov(model::StatisticalModel)
```


Return the variance-covariance matrix for the coefficients of the model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L151-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
vcov(fm2)
```


```
2×2 Matrix{Float64}:
 43.9868   -1.37039
 -1.37039   2.25671
```


```julia
vcov(gm1)
```


```
7×7 Matrix{Float64}:
  0.148396    -0.00561488   -0.00982342   …  -0.0112072    -0.011347
 -0.00561488   0.000280656   7.19122e-5      -1.47964e-5    1.0241e-5
 -0.00982342   7.19122e-5    0.0365602       -8.04427e-5    5.25878e-5
 -0.0167988   -1.43709e-5   -9.25629e-5       0.000265802  -0.000172097
 -0.0166229   -2.90553e-5   -0.000162391      0.000658966  -0.000520523
 -0.0112072   -1.47964e-5   -8.04427e-5   …   0.0228597    -0.000247782
 -0.011347     1.0241e-5     5.25878e-5      -0.000247782   0.0228032
```


The standard errors are the square roots of the diagonal elements of the estimated variance-covariance matrix of the fixed-effects coefficient estimators.
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.stderror' href='#StatsAPI.stderror'><span class="jlbinding">StatsAPI.stderror</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
stderror(model::StatisticalModel)
```


Return the standard errors for the coefficients of the model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L144-L148" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
stderror(fm2)
```


```
2-element Vector{Float64}:
 6.632257724308144
 1.5022355036397133
```


```julia
stderror(gm1)
```


```
7-element Vector{Float64}:
 0.38522187715517414
 0.016752777502287746
 0.19120720019964096
 0.1841585773287559
 0.18651728135544754
 0.15119425239025674
 0.15100718227007012
```


Finally, the `coeftable` generic produces a table of coefficient estimates, their standard errors, and their ratio. The _p-values_ quoted here should be regarded as approximations.
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.coeftable' href='#StatsAPI.coeftable'><span class="jlbinding">StatsAPI.coeftable</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
coeftable(model::StatisticalModel; level::Real=0.95)
```


Return a table with coefficients and related statistics of the model. `level` determines the level for confidence intervals (by default, 95%).

The returned `CoefTable` object implements the [Tables.jl](https://github.com/JuliaData/Tables.jl/) interface, and can be converted e.g. to a `DataFrame` via `using DataFrames; DataFrame(coeftable(model))`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/statisticalmodel.jl#L22-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
coeftable(fm2)
```


```
──────────────────────────────────────────────────
                Coef.  Std. Error      z  Pr(>|z|)
──────────────────────────────────────────────────
(Intercept)  251.405      6.63226  37.91    <1e-99
days          10.4673     1.50224   6.97    <1e-11
──────────────────────────────────────────────────
```


## Covariance parameter estimates {#Covariance-parameter-estimates}

The covariance parameters estimates, in the form shown in the model summary, are a `VarCorr` object

```julia
VarCorr(fm2)
```


```
Variance components:
            Column    Variance Std.Dev.   Corr.
subj     (Intercept)  565.51066 23.78047
         days          32.68212  5.71683 +0.08
Residual              654.94145 25.59182

```


```julia
VarCorr(gm1)
```


```
Variance components:
        Column   Variance Std.Dev. 
subj (Intercept)  1.793488 1.339212
item (Intercept)  0.117149 0.342271


```


Individual components are returned by other extractors
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.varest' href='#MixedModels.varest'><span class="jlbinding">MixedModels.varest</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
varest(m::LinearMixedModel)
```


Returns the estimate of σ², the variance of the conditional distribution of Y given B.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L1287-L1291" target="_blank" rel="noreferrer">source</a></Badge>



```julia
varest(m::GeneralizedLinearMixedModel)
```


Returns the estimate of ϕ², the variance of the conditional distribution of Y given B.

For models with a dispersion parameter ϕ, this is simply ϕ². For models without a dispersion parameter, this value is `missing`. This differs from `disperion`, which returns `1` for models without a dispersion parameter.

For Gaussian models, this parameter is often called σ².


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L814-L824" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.sdest' href='#MixedModels.sdest'><span class="jlbinding">MixedModels.sdest</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
sdest(m::LinearMixedModel)
```


Return the estimate of σ, the standard deviation of the per-observation noise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L987-L991" target="_blank" rel="noreferrer">source</a></Badge>



```julia
sdest(m::GeneralizedLinearMixedModel)
```


Return the estimate of the dispersion, i.e. the standard deviation of the per-observation noise.

For models with a dispersion parameter ϕ, this is simply ϕ. For models without a dispersion parameter, this value is `missing`. This differs from `disperion`, which returns `1` for models without a dispersion parameter.

For Gaussian models, this parameter is often called σ.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/generalizedlinearmixedmodel.jl#L710-L720" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
varest(fm2)
```


```
654.9414500374786
```


```julia
sdest(fm2)
```


```
25.591823890404502
```


```julia
fm2.σ
```


```
25.591823890404502
```


## Conditional modes of the random effects {#Conditional-modes-of-the-random-effects}

The `ranef` extractor
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.ranef' href='#MixedModels.ranef'><span class="jlbinding">MixedModels.ranef</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
ranef(m::LinearMixedModel; uscale=false)
```


Return, as a `Vector{Matrix{T}}`, the conditional modes of the random effects in model `m`.

If `uscale` is `true` the random effects are on the spherical (i.e. `u`) scale, otherwise on the original scale.

For a named variant, see [`raneftables`](/constructors#MixedModels.raneftables).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/linearmixedmodel.jl#L895-L904" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
ranef(fm1)
```


```
1-element Vector{Matrix{Float64}}:
 [-16.628221011733647 0.36951602248392174 … 53.57982326003439 -42.49434258554296]
```


```julia
fm1.b
```


```
1-element Vector{Matrix{Float64}}:
 [-16.628221011733647 0.36951602248392174 … 53.57982326003439 -42.49434258554296]
```


returns the _conditional modes_ of the random effects given the observed data. That is, these are the values that maximize the conditional density of the random effects given the observed data. For a `LinearMixedModel` these are also the conditional means.

These are sometimes called the _best linear unbiased predictors_ or [`BLUPs`](https://en.wikipedia.org/wiki/Best_linear_unbiased_prediction) but that name is not particularly meaningful.

At a superficial level these can be considered as the &quot;estimates&quot; of the random effects, with a bit of hand waving, but pursuing this analogy too far usually results in confusion.

To obtain tables associating the values of the conditional modes with the levels of the grouping factor, use
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.raneftables' href='#MixedModels.raneftables'><span class="jlbinding">MixedModels.raneftables</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
raneftables(m::MixedModel; uscale = false)
```


Return the conditional means of the random effects as a `NamedTuple` of Tables.jl-compliant tables.

::: tip Note

The API guarantee is only that the NamedTuple contains Tables.jl tables and not on the particular concrete type of each table.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/MixedModels.jl/blob/v4.31.0/src/mixedmodel.jl#L166-L173" target="_blank" rel="noreferrer">source</a></Badge>

</details>


as in

```julia
DataFrame(only(raneftables(fm1)))
```

<div v-html="`&lt;div&gt;&lt;div style = &quot;float: left;&quot;&gt;&lt;span&gt;6×2 DataFrame&lt;/span&gt;&lt;/div&gt;&lt;div style = &quot;clear: both;&quot;&gt;&lt;/div&gt;&lt;/div&gt;&lt;div class = &quot;data-frame&quot; style = &quot;overflow-x: scroll;&quot;&gt;&lt;table class = &quot;data-frame&quot; style = &quot;margin-bottom: 6px;&quot;&gt;&lt;thead&gt;&lt;tr class = &quot;header&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;Row&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;batch&lt;/th&gt;&lt;th style = &quot;text-align: left;&quot;&gt;(Intercept)&lt;/th&gt;&lt;/tr&gt;&lt;tr class = &quot;subheader headerLastRow&quot;&gt;&lt;th class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;&lt;/th&gt;&lt;th title = &quot;String&quot; style = &quot;text-align: left;&quot;&gt;String&lt;/th&gt;&lt;th title = &quot;Float64&quot; style = &quot;text-align: left;&quot;&gt;Float64&lt;/th&gt;&lt;/tr&gt;&lt;/thead&gt;&lt;tbody&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;1&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;A&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;-16.6282&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;2&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;B&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;0.369516&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;3&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;C&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;26.9747&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;4&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;D&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;-21.8014&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;5&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;E&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;53.5798&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td class = &quot;rowNumber&quot; style = &quot;font-weight: bold; text-align: right;&quot;&gt;6&lt;/td&gt;&lt;td style = &quot;text-align: left;&quot;&gt;F&lt;/td&gt;&lt;td style = &quot;text-align: right;&quot;&gt;-42.4943&lt;/td&gt;&lt;/tr&gt;&lt;/tbody&gt;&lt;/table&gt;&lt;/div&gt;`"></div>

The corresponding conditional variances are returned by
<details class='jldocstring custom-block' open>
<summary><a id='MixedModels.condVar' href='#MixedModels.condVar'><span class="jlbinding">MixedModels.condVar</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



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


```julia
condVar(fm1)
```


```
1-element Vector{Array{Float64, 3}}:
 [362.3104675622471;;; 362.3104675622471;;; 362.3104675622471;;; 362.3104675622471;;; 362.3104675622471;;; 362.3104675622471]
```


## Case-wise diagnostics and residual degrees of freedom {#Case-wise-diagnostics-and-residual-degrees-of-freedom}

The `leverage` values
<details class='jldocstring custom-block' open>
<summary><a id='StatsAPI.leverage' href='#StatsAPI.leverage'><span class="jlbinding">StatsAPI.leverage</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
leverage(model::RegressionModel)
```


Return the diagonal of the projection matrix of the model.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/JuliaStats/StatsAPI.jl/blob/v1.7.0/src/regressionmodel.jl#L51-L55" target="_blank" rel="noreferrer">source</a></Badge>

</details>


```julia
leverage(fm1)
```


```
30-element Vector{Float64}:
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 ⋮
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
 0.15650534082766315
```


are used in diagnostics for linear regression models to determine cases that exert a strong influence on their own predicted response.

The documentation refers to a &quot;projection&quot;. For a linear model without random effects the fitted values are obtained by orthogonal projection of the response onto the column span of the model matrix and the sum of the leverage values is the dimension of this column span. That is, the sum of the leverage values is the rank of the model matrix and `n - sum(leverage(m))` is the degrees of freedom for residuals. The sum of the leverage values is also the trace of the so-called &quot;hat&quot; matrix, `H`. (The name &quot;hat matrix&quot; reflects the fact that $\hat{\mathbf{y}} = \mathbf{H} \mathbf{y}$.  That is, `H` puts a hat on `y`.)

For a linear mixed model the sum of the leverage values will be between `p`, the rank of the fixed-effects model matrix, and `p + q` where `q` is the total number of random effects. This number does not represent a dimension (or &quot;degrees of freedom&quot;) of a linear subspace of all possible fitted values because the projection is not an orthogonal projection. Nevertheless, it is a reasonable measure of the effective degrees of freedom of the model and `n - sum(leverage(m))` can be considered the effective residual degrees of freedom.

For model `fm1` the dimensions are

```julia
n, p, q, k = size(fm1)
```


```
(30, 1, 6, 1)
```


which implies that the sum of the leverage values should be in the range [1, 7]. The actual value is

```julia
sum(leverage(fm1))
```


```
4.695160224829895
```


For model `fm2` the dimensions are

```julia
n, p, q, k = size(fm2)
```


```
(180, 2, 36, 1)
```


providing a range of [2, 38] for the effective degrees of freedom for the model. The observed value is

```julia
sum(leverage(fm2))
```


```
28.611526266274847
```


When a model converges to a singular covariance, such as

```julia
fm3 = fit(MixedModel, @formula(yield ~ 1+(1|batch)), MixedModels.dataset(:dyestuff2))
```


```
Linear mixed model fit by maximum likelihood
 yield ~ 1 + (1 | batch)
   logLik   -2 logLik     AIC       AICc        BIC    
   -81.4365   162.8730   168.8730   169.7961   173.0766

Variance components:
            Column   Variance Std.Dev.
batch    (Intercept)   0.00000 0.00000
Residual              13.34610 3.65323
 Number of obs: 30; levels of grouping factors: 6

  Fixed-effects parameters:
───────────────────────────────────────────────
              Coef.  Std. Error     z  Pr(>|z|)
───────────────────────────────────────────────
(Intercept)  5.6656    0.666986  8.49    <1e-16
───────────────────────────────────────────────
```


the effective degrees of freedom is the lower bound.

```julia
sum(leverage(fm3))
```


```
0.9999999999999998
```


Models for which the estimates of the variances of the random effects are large relative to the residual variance have effective degrees of freedom close to the upper bound.

```julia
fm4 = fit(MixedModel, @formula(diameter ~ 1+(1|plate)+(1|sample)),
    MixedModels.dataset(:penicillin))
```


```
Linear mixed model fit by maximum likelihood
 diameter ~ 1 + (1 | plate) + (1 | sample)
   logLik   -2 logLik     AIC       AICc        BIC    
  -166.0942   332.1883   340.1883   340.4761   352.0676

Variance components:
            Column   Variance Std.Dev. 
plate    (Intercept)  0.714979 0.845565
sample   (Intercept)  3.135193 1.770648
Residual              0.302426 0.549933
 Number of obs: 144; levels of grouping factors: 24, 6

  Fixed-effects parameters:
─────────────────────────────────────────────────
               Coef.  Std. Error      z  Pr(>|z|)
─────────────────────────────────────────────────
(Intercept)  22.9722    0.744596  30.85    <1e-99
─────────────────────────────────────────────────
```


```julia
sum(leverage(fm4))
```


```
27.46531792078765
```


Also, a model fit by the REML criterion generally has larger estimates of the variance components and hence a larger effective degrees of freedom.

```julia
fm4r = fit(MixedModel, @formula(diameter ~ 1+(1|plate)+(1|sample)),
    MixedModels.dataset(:penicillin), REML=true)
```


```
Linear mixed model fit by REML
 diameter ~ 1 + (1 | plate) + (1 | sample)
 REML criterion at convergence: 330.8605889912561

Variance components:
            Column   Variance Std.Dev. 
plate    (Intercept)  0.716908 0.846704
sample   (Intercept)  3.730907 1.931555
Residual              0.302416 0.549923
 Number of obs: 144; levels of grouping factors: 24, 6

  Fixed-effects parameters:
─────────────────────────────────────────────────
               Coef.  Std. Error      z  Pr(>|z|)
─────────────────────────────────────────────────
(Intercept)  22.9722    0.808572  28.41    <1e-99
─────────────────────────────────────────────────
```


```julia
sum(leverage(fm4r))
```


```
27.472361698411703
```

