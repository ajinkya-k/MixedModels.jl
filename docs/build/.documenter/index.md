
# MixedModels.jl Documentation {#MixedModels.jl-Documentation}



_MixedModels.jl_ is a Julia package providing capabilities for fitting and examining linear and generalized linear mixed-effect models. It is similar in scope to the [_lme4_](https://github.com/lme4/lme4) package for `R`.
- [Model constructors](constructors#Model-constructors)
    - [Examples of linear mixed effects model fits](constructors#Examples-of-linear-mixed-effects-model-fits)
    - [Fitting generalized linear mixed models](constructors#Fitting-generalized-linear-mixed-models)
- [Extractor functions](constructors#Extractor-functions)
    - [Model fit statistics](constructors#Model-fit-statistics)
    - [Fixed effects parameter estimates](constructors#Fixed-effects-parameter-estimates)
    - [Covariance parameter estimates](constructors#Covariance-parameter-estimates)
    - [Conditional modes of the random effects](constructors#Conditional-modes-of-the-random-effects)
    - [Case wise diagnostics and residual degrees of freedom](constructors#Case-wise-diagnostics-and-residual-degrees-of-freedom)
- [Details of the parameter estimation](optimization#Details-of-the-parameter-estimation)
    - [The probability model](optimization#The-probability-model)
    - [Linear Mixed Effects Models](optimization#Linear-Mixed-Effects-Models)
    - [Internal structure of \Lambda_\theta and \bf Z](optimization#Internal-structure-of-\Lambda_\theta-and-\bf-Z)
    - [A blocked Cholesky factor](optimization#A-blocked-Cholesky-factor)
    - [Generalized Linear Mixed Effects Models](optimization#Generalized-Linear-Mixed-Effects-Models)
- [Normalized Gauss Hermite Quadrature](GaussHermite#Normalized-Gauss-Hermite-Quadrature)
    - [Evaluating the weights and abscissae](GaussHermite#Evaluating-the-weights-and-abscissae)
    - [Application to a model for contraception use](GaussHermite#Application-to-a-model-for-contraception-use)
- [Parametric bootstrap for mixed effects models](bootstrap#Parametric-bootstrap-for-mixed-effects-models)
    - [The parametric bootstrap](bootstrap#The-parametric-bootstrap)
    - [Reduced Precision Bootstrap](bootstrap#Reduced-Precision-Bootstrap)
    - [Distributed Computing and the Bootstrap](bootstrap#Distributed-Computing-and-the-Bootstrap)
- [Rank deficiency in mixed effects models](rankdeficiency#Rank-deficiency-in-mixed-effects-models)
    - [Fixed effects](rankdeficiency#Fixed-effects)
    - [Random effects](rankdeficiency#Random-effects)
- [Alternative display and output formats](mime#Alternative-display-and-output-formats)

