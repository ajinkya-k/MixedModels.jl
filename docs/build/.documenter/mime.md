
# Alternative display and output formats {#Alternative-display-and-output-formats}

In the documentation, we have presented the output from MixedModels.jl in the same format you will see when working in the REPL. You may have noticed, however, that output from other packages received pretty printing. For example, DataFrames are converted into nice HTML tables. In MixedModels, we recently (v3.2.0) introduced limited support for such pretty printing. (For more details on how the print and display system in Julia works, check out [this NextJournal post](https://nextjournal.com/sdanisch/julias-display-system).)

In particular, we have defined Markdown, HTML and LaTeX output, i.e. `show` methods, for our types. Note that the Markdown output can also be easily and more flexibly translated into HTML, LaTeX (e.g. with `booktabs`) or even a MS Word Document using tools such as [pandoc](https://pandoc.org/). Packages like `IJulia` and `Documenter` can often detect the presence of these display options and use them automatically.

```julia
using MixedModels
form = @formula(rt_trunc ~ 1 + spkr * prec * load +
                          (1 + load | item) +
                          (1 + spkr + prec + load | subj))
contr = Dict(:spkr => EffectsCoding(),
             :prec => EffectsCoding(),
             :load => EffectsCoding(),
             :item => Grouping(),
             :subj => Grouping())
kbm = fit(MixedModel, form, MixedModels.dataset(:kb07); contrasts=contr)
```

<div v-html="`&lt;table&gt;&lt;tr&gt;&lt;th align=&quot;left&quot;&gt;&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;Est.&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;SE&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;z&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;p&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;σ_subj&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;σ_item&lt;/th&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;&amp;#40;Intercept&amp;#41;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;2182.0690&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;78.1911&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;27.91&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;lt;1e-99&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;318.9314&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;358.3855&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;spkr: old&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;67.9659&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;19.0808&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;3.56&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.0004&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;67.0805&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;prec: maintain&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-333.7038&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;18.5949&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-17.95&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;lt;1e-71&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;58.9474&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;78.3741&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;19.1623&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;4.09&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;lt;1e-04&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;63.2114&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;19.6931&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;spkr: old &amp;amp; prec: maintain&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-21.5694&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;16.8440&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-1.28&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.2004&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;spkr: old &amp;amp; load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;18.1670&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;16.8441&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;1.08&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.2808&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;prec: maintain &amp;amp; load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;4.3165&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;16.8441&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.26&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.7977&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;spkr: old &amp;amp; prec: maintain &amp;amp; load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;23.2112&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;16.8440&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;1.38&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;0.1682&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Residual&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;712.4110&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;&#10;`"></div>

Note that the display here is more succinct than the standard REPL display:

```julia
using DisplayAs
kbm |> DisplayAs.Text
```


```
Linear mixed model fit by maximum likelihood
 rt_trunc ~ 1 + spkr + prec + load + spkr & prec + spkr & load + prec & load + spkr & prec & load + (1 + load | item) + (1 + spkr + prec + load | subj)
    logLik   -2 logLik      AIC         AICc        BIC     
 -14400.6879  28801.3758  28845.3758  28845.9489  28966.1429

Variance components:
             Column       Variance  Std.Dev.   Corr.
subj     (Intercept)     101717.2454 318.9314
         spkr: old         4499.7962  67.0805 +0.70
         prec: maintain    3474.7975  58.9474 -0.71 -0.00
         load: yes         3995.6822  63.2114 +0.28 +0.88 +0.47
item     (Intercept)     128440.1457 358.3855
         load: yes          387.8169  19.6931 +0.80
Residual                 507529.4719 712.4110
 Number of obs: 1789; levels of grouping factors: 56, 32

  Fixed-effects parameters:
────────────────────────────────────────────────────────────────────────────────
                                             Coef.  Std. Error       z  Pr(>|z|)
────────────────────────────────────────────────────────────────────────────────
(Intercept)                             2182.07        78.1911   27.91    <1e-99
spkr: old                                 67.9659      19.0808    3.56    0.0004
prec: maintain                          -333.704       18.5949  -17.95    <1e-71
load: yes                                 78.3741      19.1623    4.09    <1e-04
spkr: old & prec: maintain               -21.5694      16.844    -1.28    0.2004
spkr: old & load: yes                     18.167       16.8441    1.08    0.2808
prec: maintain & load: yes                 4.31651     16.8441    0.26    0.7977
spkr: old & prec: maintain & load: yes    23.2112      16.844     1.38    0.1682
────────────────────────────────────────────────────────────────────────────────
```


This brevity is intentional: we wanted these types to work well with traditional academic publishing constraints on tables. The summary for a model fit presented in the REPL does not mesh well with being treated as a single table (with columns shared between the random and fixed effects). In our experience, this leads to difficulties in typesetting the resulting tables. We nonetheless encourage users to report fit statistics such as the log likelihood or AIC as part of the caption of their table. If the correlation parameters in the random effects are of interest, then [`VarCorr`](/api#MixedModels.VarCorr) can also be pretty printed:

```julia
VarCorr(kbm)
```

<div v-html="`&lt;table&gt;&lt;tr&gt;&lt;th align=&quot;left&quot;&gt; &lt;/th&gt;&lt;th align=&quot;left&quot;&gt;Column&lt;/th&gt;&lt;th align=&quot;right&quot;&gt; Variance&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;Std.Dev&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;Corr.&lt;/th&gt;&lt;th align=&quot;right&quot;&gt; &lt;/th&gt;&lt;th align=&quot;right&quot;&gt; &lt;/th&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;subj&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#40;Intercept&amp;#41;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;101717.2454&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;318.9314&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt;spkr: old&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;4499.7962&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;67.0805&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;#43;0.70&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt;prec: maintain&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;3474.7975&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;58.9474&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-0.71&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;-0.00&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt;load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;3995.6822&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;63.2114&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;#43;0.28&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;#43;0.88&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;#43;0.47&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;item&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#40;Intercept&amp;#41;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;128440.1457&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;358.3855&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt;load: yes&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;387.8169&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;19.6931&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;&amp;#43;0.80&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Residual&lt;/td&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt;507529.4719&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;712.4110&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;&#10;`"></div>

Similarly for [`BlockDescription`](/api#MixedModels.BlockDescription), `OptSummary` and `MixedModels.likelihoodratiotest`:

```julia
BlockDescription(kbm)
```

<div v-html="`&lt;table&gt;&lt;tr&gt;&lt;th align=&quot;left&quot;&gt;rows&lt;/th&gt;&lt;th align=&quot;left&quot;&gt;subj&lt;/th&gt;&lt;th align=&quot;left&quot;&gt;item&lt;/th&gt;&lt;th align=&quot;left&quot;&gt;fixed&lt;/th&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;224&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;BlkDiag&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;64&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;Dense&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;BlkDiag/Dense&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;9&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;Dense&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;Dense&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;Dense&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;&#10;`"></div>

```julia
kbm.optsum
```

<div v-html="`&lt;table&gt;&lt;tr&gt;&lt;th align=&quot;left&quot;&gt;&lt;/th&gt;&lt;th align=&quot;left&quot;&gt;&lt;/th&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;&lt;b&gt;Initialization&lt;/b&gt;&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Initial parameter vector&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#91;1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0&amp;#93;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Initial objective value&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;29340.042234597688&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;&lt;b&gt;Optimizer settings&lt;/b&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Optimizer&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;code&gt;LN_BOBYQA&lt;/code&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Backend&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;code&gt;nlopt&lt;/code&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Lower bounds&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#91;0.0, -Inf, -Inf, -Inf, 0.0, -Inf, -Inf, 0.0, -Inf, 0.0, 0.0, -Inf, 0.0&amp;#93;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;ftol_rel&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;1.0e-12&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;ftol_abs&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;1.0e-8&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;xtol_rel&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;0.0&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;xtol_abs&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#91;1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10&amp;#93;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;initial_step&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#91;0.75, 1.0, 1.0, 1.0, 0.75, 1.0, 1.0, 0.75, 1.0, 0.75, 0.75, 1.0, 0.75&amp;#93;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;maxfeval&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;-1&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;maxtime&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;-1.0&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;xtol_zero_abs&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;0.001&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;ftol_zero_abs&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;1.0e-5&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;&lt;b&gt;Result&lt;/b&gt;&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Function evaluations&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;317&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Final parameter vector&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;#91;0.4477, 0.066, -0.0592, 0.0249, 0.0672, 0.0578, 0.0852, 0.0, 0.0002, 0.0, 0.5031, 0.022, 0.0167&amp;#93;&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Final objective value&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;28801.3758&lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;Return code&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&lt;code&gt;FTOL_REACHED&lt;/code&gt;&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;&#10;`"></div>

```julia
m0 = fit(MixedModel, @formula(reaction ~ 1 + (1|subj)), MixedModels.dataset(:sleepstudy))
m1 = fit(MixedModel, @formula(reaction ~ 1 + days + (1+days|subj)), MixedModels.dataset(:sleepstudy))
MixedModels.likelihoodratiotest(m0,m1)
```

<div v-html="`&lt;table&gt;&lt;tr&gt;&lt;th align=&quot;left&quot;&gt;&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;model-dof&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;deviance&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;χ²&lt;/th&gt;&lt;th align=&quot;right&quot;&gt;χ²-dof&lt;/th&gt;&lt;th align=&quot;left&quot;&gt;P&amp;#40;&amp;gt;χ²&amp;#41;&lt;/th&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;reaction ~ 1 &amp;#43; &amp;#40;1 | subj&amp;#41;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;3&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;1911&lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;right&quot;&gt; &lt;/td&gt;&lt;td align=&quot;left&quot;&gt; &lt;/td&gt;&lt;/tr&gt;&lt;tr&gt;&lt;td align=&quot;left&quot;&gt;reaction ~ 1 &amp;#43; days &amp;#43; &amp;#40;1 &amp;#43; days | subj&amp;#41;&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;6&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;1752&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;159&lt;/td&gt;&lt;td align=&quot;right&quot;&gt;3&lt;/td&gt;&lt;td align=&quot;left&quot;&gt;&amp;lt;1e-33&lt;/td&gt;&lt;/tr&gt;&lt;/table&gt;&#10;`"></div>

To explicitly invoke this behavior, we must specify the right `show` method. (The raw and not rendered output is intentionally shown here.)

```julia
show(MIME("text/markdown"), m1)
```


```
|                                        |      Est. |      SE |      z |      p |   σ_subj |   σ_item |
|:-------------------------------------- | ---------:| -------:| ------:| ------:| --------:| --------:|
| (Intercept)                            | 2182.0690 | 78.1911 |  27.91 | <1e-99 | 318.9314 | 358.3855 |
| spkr: old                              |   67.9659 | 19.0808 |   3.56 | 0.0004 |  67.0805 |          |
| prec: maintain                         | -333.7038 | 18.5949 | -17.95 | <1e-71 |  58.9474 |          |
| load: yes                              |   78.3741 | 19.1623 |   4.09 | <1e-04 |  63.2114 |  19.6931 |
| spkr: old & prec: maintain             |  -21.5694 | 16.8440 |  -1.28 | 0.2004 |          |          |
| spkr: old & load: yes                  |   18.1670 | 16.8441 |   1.08 | 0.2808 |          |          |
| prec: maintain & load: yes             |    4.3165 | 16.8441 |   0.26 | 0.7977 |          |          |
| spkr: old & prec: maintain & load: yes |   23.2112 | 16.8440 |   1.38 | 0.1682 |          |          |
| Residual                               |  712.4110 |         |        |        |          |          |
```


```julia
show(MIME("text/html"), m1)
```


```
<table><tr><th align="left"></th><th align="right">Est.</th><th align="right">SE</th><th align="right">z</th><th align="right">p</th><th align="right">σ_subj</th><th align="right">σ_item</th></tr><tr><td align="left">&#40;Intercept&#41;</td><td align="right">2182.0690</td><td align="right">78.1911</td><td align="right">27.91</td><td align="right">&lt;1e-99</td><td align="right">318.9314</td><td align="right">358.3855</td></tr><tr><td align="left">spkr: old</td><td align="right">67.9659</td><td align="right">19.0808</td><td align="right">3.56</td><td align="right">0.0004</td><td align="right">67.0805</td><td align="right"> </td></tr><tr><td align="left">prec: maintain</td><td align="right">-333.7038</td><td align="right">18.5949</td><td align="right">-17.95</td><td align="right">&lt;1e-71</td><td align="right">58.9474</td><td align="right"> </td></tr><tr><td align="left">load: yes</td><td align="right">78.3741</td><td align="right">19.1623</td><td align="right">4.09</td><td align="right">&lt;1e-04</td><td align="right">63.2114</td><td align="right">19.6931</td></tr><tr><td align="left">spkr: old &amp; prec: maintain</td><td align="right">-21.5694</td><td align="right">16.8440</td><td align="right">-1.28</td><td align="right">0.2004</td><td align="right"> </td><td align="right"> </td></tr><tr><td align="left">spkr: old &amp; load: yes</td><td align="right">18.1670</td><td align="right">16.8441</td><td align="right">1.08</td><td align="right">0.2808</td><td align="right"> </td><td align="right"> </td></tr><tr><td align="left">prec: maintain &amp; load: yes</td><td align="right">4.3165</td><td align="right">16.8441</td><td align="right">0.26</td><td align="right">0.7977</td><td align="right"> </td><td align="right"> </td></tr><tr><td align="left">spkr: old &amp; prec: maintain &amp; load: yes</td><td align="right">23.2112</td><td align="right">16.8440</td><td align="right">1.38</td><td align="right">0.1682</td><td align="right"> </td><td align="right"> </td></tr><tr><td align="left">Residual</td><td align="right">712.4110</td><td align="right"></td><td align="right"></td><td align="right"></td><td align="right"></td><td align="right"></td></tr></table>
```


Note for that LaTeX, the column labels for the random effects are slightly changed: σ is placed into math mode and escaped and the grouping variable is turned into a subscript. Similarly for the likelihood ratio test, the χ² is escaped into math mode. This transformation improves pdfLaTeX and journal compatibility, but also means that XeLaTeX and LuaTeX may use a different font at this point.

```julia
show(MIME("text/latex"), m1)
```


```
\begin{tabular}
{l | r | r | r | r | r | r}
 & Est. & SE & z & p & $\sigma_\text{subj}$ & $\sigma_\text{item}$ \\
\hline
(Intercept) & 2182.0690 & 78.1911 & 27.91 & <1e-99 & 318.9314 & 358.3855 \\
spkr: old & 67.9659 & 19.0808 & 3.56 & 0.0004 & 67.0805 &   \\
prec: maintain & -333.7038 & 18.5949 & -17.95 & <1e-71 & 58.9474 &   \\
load: yes & 78.3741 & 19.1623 & 4.09 & <1e-04 & 63.2114 & 19.6931 \\
spkr: old \& prec: maintain & -21.5694 & 16.8440 & -1.28 & 0.2004 &   &   \\
spkr: old \& load: yes & 18.1670 & 16.8441 & 1.08 & 0.2808 &   &   \\
prec: maintain \& load: yes & 4.3165 & 16.8441 & 0.26 & 0.7977 &   &   \\
spkr: old \& prec: maintain \& load: yes & 23.2112 & 16.8440 & 1.38 & 0.1682 &   &   \\
Residual & 712.4110 &  &  &  &  &  \\
\end{tabular}
```


This escaping behavior can be disabled by specifying `"text/xelatex"` as the MIME type. (Note that other symbols may still be escaped, as the internal conversion uses the `Markdown` module from the standard library, which performs some escaping on its own.)

```julia
show(MIME("text/xelatex"), m1)
```


```
\begin{tabular}
{l | r | r | r | r | r | r}
 & Est. & SE & z & p & σ\_subj & σ\_item \\
\hline
(Intercept) & 2182.0690 & 78.1911 & 27.91 & <1e-99 & 318.9314 & 358.3855 \\
spkr: old & 67.9659 & 19.0808 & 3.56 & 0.0004 & 67.0805 &   \\
prec: maintain & -333.7038 & 18.5949 & -17.95 & <1e-71 & 58.9474 &   \\
load: yes & 78.3741 & 19.1623 & 4.09 & <1e-04 & 63.2114 & 19.6931 \\
spkr: old \& prec: maintain & -21.5694 & 16.8440 & -1.28 & 0.2004 &   &   \\
spkr: old \& load: yes & 18.1670 & 16.8441 & 1.08 & 0.2808 &   &   \\
prec: maintain \& load: yes & 4.3165 & 16.8441 & 0.26 & 0.7977 &   &   \\
spkr: old \& prec: maintain \& load: yes & 23.2112 & 16.8440 & 1.38 & 0.1682 &   &   \\
Residual & 712.4110 &  &  &  &  &  \\
\end{tabular}
```


This output can also be written directly to file:

```julia
open("model.md", "w") do io
    show(io, MIME("text/markdown"), kbm)
end
```

