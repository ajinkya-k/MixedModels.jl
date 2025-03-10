import{_ as l,c as i,o as n,aA as t,j as s,a}from"./chunks/framework.BeL8ovVr.js";const Q=JSON.parse('{"title":"Rank deficiency in mixed-effects models","description":"","frontmatter":{},"headers":[],"relativePath":"rankdeficiency.md","filePath":"rankdeficiency.md","lastUpdated":null}'),o={name:"rankdeficiency.md"},r={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},p={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.052ex",height:"1.025ex",role:"img",focusable:"false",viewBox:"0 -442 465 453","aria-hidden":"true"},d={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},h={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.439ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.138ex",height:"1.439ex",role:"img",focusable:"false",viewBox:"0 -442 503 636","aria-hidden":"true"},c={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},k={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.357ex",height:"1.025ex",role:"img",focusable:"false",viewBox:"0 -442 600 453","aria-hidden":"true"},m={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},f={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.025ex"},xmlns:"http://www.w3.org/2000/svg",width:"1.357ex",height:"1.025ex",role:"img",focusable:"false",viewBox:"0 -442 600 453","aria-hidden":"true"},u={class:"MathJax",jax:"SVG",style:{direction:"ltr",position:"relative"}},g={style:{overflow:"visible","min-height":"1px","min-width":"1px","vertical-align":"-0.186ex"},xmlns:"http://www.w3.org/2000/svg",width:"5.254ex",height:"1.692ex",role:"img",focusable:"false",viewBox:"0 -666 2322.4 748","aria-hidden":"true"};function y(x,e,b,w,T,E){return n(),i("div",null,[e[25]||(e[25]=t('<h1 id="Rank-deficiency-in-mixed-effects-models" tabindex="-1">Rank deficiency in mixed-effects models <a class="header-anchor" href="#Rank-deficiency-in-mixed-effects-models" aria-label="Permalink to &quot;Rank deficiency in mixed-effects models {#Rank-deficiency-in-mixed-effects-models}&quot;">​</a></h1><p>The <em>(column) rank</em> of a matrix refers to the number of linearly independent columns in the matrix. Clearly, the rank can never be more than the number of columns; however, the rank can be less than the number of columns. In a regression context, this corresponds to a (linear) dependency in the predictors. The simplest case of rank deficiency is a duplicated predictor or a predictor that is exactly a multiple of another predictor. However, rank deficiency can also arise in more subtle ways, such as from missing cells in a two-factor experimental design. Rank deficiency can also arise as an extreme case of multicollinearity. In all cases, it is important to remember that we can only assess the numerical rank of a matrix, which may be less than its theoretical rank, and that evaluation of this numerical rank requires setting some numerical tolerance levels. These choices are not always well defined. In other words, the rank of a matrix is well-defined in theory but in practice can be difficult to evaluate.</p><p>Rank deficiency can occur in two ways in mixed-effects models: in the fixed effects and in the random effects. The implications of rank deficiency and thus the handling of it differ between these.</p><h2 id="Fixed-effects" tabindex="-1">Fixed effects <a class="header-anchor" href="#Fixed-effects" aria-label="Permalink to &quot;Fixed effects {#Fixed-effects}&quot;">​</a></h2><p>The consequences of rank deficiency in the fixed effects are similar to those in classical ordinary least squares (OLS) regression. If one or more predictors can be expressed as a linear combination of the other columns, then this column is redundant and the model matrix is rank deficient. Note however, that the redundant column is not defined uniquely. For example, in the case that of two columns <code>a</code> and <code>b</code> where <code>b = 2a</code>, then the rank deficiency can be handled by eliminating either <code>a</code> or <code>b</code>. While we defined <code>b</code> here in terms of <code>a</code>, it may be that <code>b</code> is actually the more &#39;fundamental&#39; predictor and hence we may define <code>a</code> in terms of <code>b</code> as <code>a = 0.5b</code>. The user may of course possess this information, but the choice is not apparent to the modelling software. As such, the handling of rank deficiency in <code>MixedModels.jl</code> should not be taken as a replacement for thinking about the nature of the predictors in a given model.</p>',5)),s("p",null,[e[4]||(e[4]=a("There is a widely accepted convention for how to make the coefficient estimates for these redundant columns well-defined: we set their value to zero and their standard errors to ")),e[5]||(e[5]=s("code",null,"NaN",-1)),e[6]||(e[6]=a(" (and thus also their ")),s("mjx-container",r,[(n(),i("svg",p,e[0]||(e[0]=[s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D467",d:"M347 338Q337 338 294 349T231 360Q211 360 197 356T174 346T162 335T155 324L153 320Q150 317 138 317Q117 317 117 325Q117 330 120 339Q133 378 163 406T229 440Q241 442 246 442Q271 442 291 425T329 392T367 375Q389 375 411 408T434 441Q435 442 449 442H462Q468 436 468 434Q468 430 463 420T449 399T432 377T418 358L411 349Q368 298 275 214T160 106L148 94L163 93Q185 93 227 82T290 71Q328 71 360 90T402 140Q406 149 409 151T424 153Q443 153 443 143Q443 138 442 134Q425 72 376 31T278 -11Q252 -11 232 6T193 40T155 57Q111 57 76 -3Q70 -11 59 -11H54H41Q35 -5 35 -2Q35 13 93 84Q132 129 225 214T340 322Q352 338 347 338Z",style:{"stroke-width":"3"}})])])],-1)]))),e[1]||(e[1]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"z")])],-1))]),e[7]||(e[7]=a(" and ")),s("mjx-container",d,[(n(),i("svg",h,e[2]||(e[2]=[s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D45D",d:"M23 287Q24 290 25 295T30 317T40 348T55 381T75 411T101 433T134 442Q209 442 230 378L240 387Q302 442 358 442Q423 442 460 395T497 281Q497 173 421 82T249 -10Q227 -10 210 -4Q199 1 187 11T168 28L161 36Q160 35 139 -51T118 -138Q118 -144 126 -145T163 -148H188Q194 -155 194 -157T191 -175Q188 -187 185 -190T172 -194Q170 -194 161 -194T127 -193T65 -192Q-5 -192 -24 -194H-32Q-39 -187 -39 -183Q-37 -156 -26 -148H-6Q28 -147 33 -136Q36 -130 94 103T155 350Q156 355 156 364Q156 405 131 405Q109 405 94 377T71 316T59 280Q57 278 43 278H29Q23 284 23 287ZM178 102Q200 26 252 26Q282 26 310 49T356 107Q374 141 392 215T411 325V331Q411 405 350 405Q339 405 328 402T306 393T286 380T269 365T254 350T243 336T235 326L232 322Q232 321 229 308T218 264T204 212Q178 106 178 102Z",style:{"stroke-width":"3"}})])])],-1)]))),e[3]||(e[3]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"p")])],-1))]),e[8]||(e[8]=a("-values). The values that have been defined to be zero, as opposed to evaluating to zero, are displayed as ")),e[9]||(e[9]=s("code",null,"-0.0",-1)),e[10]||(e[10]=a(" as an additional visual aid to distinguish them from the other coefficients. In practice the determination of rank and the redundant coefficients is done via a 'pivoting' scheme during a decomposition to move the surplus columns to the right side of the model matrix. In subsequent calculations, these columns are effectively ignored (as their estimates are zero and thus won't contribute to any other computations). For display purposes, this pivoting is unwound when the ")),e[11]||(e[11]=s("code",null,"coef",-1)),e[12]||(e[12]=a(" values are displayed."))]),e[26]||(e[26]=t('<p>Both the pivoted and unpivoted coefficients are available in MixedModels. The <a href="/MixedModels.jl/constructors#MixedModels.fixef"><code>fixef</code></a> extractor returns the pivoted, truncated estimates (i.e. the non redundant terms), while the <a href="/MixedModels.jl/constructors#StatsAPI.coef"><code>coef</code></a> extractor returns the unpivoted estimates (i.e. all terms, included the redundant ones). The same holds for the associated <a href="/MixedModels.jl/constructors#MixedModels.fixefnames"><code>fixefnames</code></a> and <a href="/MixedModels.jl/constructors#StatsAPI.coefnames"><code>coefnames</code></a>.</p><h3 id="Pivoting-is-platform-dependent" tabindex="-1">Pivoting is platform dependent <a class="header-anchor" href="#Pivoting-is-platform-dependent" aria-label="Permalink to &quot;Pivoting is platform dependent {#Pivoting-is-platform-dependent}&quot;">​</a></h3><p>In MixedModels.jl, we use standard numerical techniques to detect rank deficiency. We currently offer no guarantees as to which exactly of the standard techniques (pivoted QR decomposition, pivoted Cholesky decomposition, etc.) will be used. This choice should be viewed as an implementation detail. Similarly, we offer no guarantees as to which of columns will be treated as redundant. This choice may vary between releases and even between platforms (both in broad strokes of &quot;Linux&quot; vs. &quot;Windows&quot; and at the level of which BLAS options are loaded on a given processor architecture) for the same release. In other words, <em>you should not rely on the order of the pivoted columns being consistent!</em> when you switch to a different computer or a different operating system. If consistency in the pivoted columns is important to you, then you should instead determine your rank ahead of time and remove extraneous columns / predictors from your model specification.</p><p>This lack of consistency guarantees arises from a more fundamental issue: numeric linear algebra is challenging and sensitive to the underlying floating point operations. Due to rounding error, floating point arithmetic is not associative:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> +</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> +</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> ==</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> +</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> +</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">0.1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 0.3</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>false</span></span></code></pre></div><p>This means that &quot;nearly&quot; / numerically rank deficient matrices may or may not be detected as rank deficient, depending on details of the platform. Determining the rank of a matrix is the type of problem that is well-defined in theory but not in practice.</p><p>Currently, a coarse heuristic is applied to reduce the chance that the intercept column will be pivoted, but even this behavior is not guaranteed.</p><h3 id="Undetected-Rank-Deficiency" tabindex="-1">Undetected Rank Deficiency <a class="header-anchor" href="#Undetected-Rank-Deficiency" aria-label="Permalink to &quot;Undetected Rank Deficiency {#Undetected-Rank-Deficiency}&quot;">​</a></h3><p>Undetected rank deficiency in the fixed effects will lead to numerical issues, such as nonsensical estimates. A <code>PosDefException</code> may indicate rank deficiency because the covariance matrix will only be positive semidefinite and not positive definite (see <a href="/MixedModels.jl/optimization#Details-of-the-parameter-estimation">Details of the parameter estimation</a>). In other words, checking that the fixed effects are full rank is a great first step in debugging a <code>PosDefException</code>.</p><p>Note that <code>PosDefException</code> is not specific to rank deficiency and may arise in other ill-conditioned models. In any case, examining the model specification and the data to verify that they work together is the first step. For generalized linear mixed-effects models, it may also be worthwhile to try out <code>fast=true</code> instead of the default <code>fast=false</code>. See this <a href="https://github.com/JuliaStats/MixedModels.jl/issues/349" target="_blank" rel="noreferrer">GitHub issue</a> and linked Discourse discussion for more information.</p><h2 id="Random-effects" tabindex="-1">Random effects <a class="header-anchor" href="#Random-effects" aria-label="Permalink to &quot;Random effects {#Random-effects}&quot;">​</a></h2><p>Rank deficiency presents less of a problem in the random effects than in the fixed effects because the &quot;estimates&quot; (more formally, the conditional modes of the random effects given the observed data) are determined as the solution to a penalized least squares problem. The <em>shrinkage</em> effect which moves the conditional modes (group-level predictions) towards the grand mean is a form of <em>regularization</em>, which provides well-defined &quot;estimates&quot; for overparameterized models. (For more reading on this general idea, see also this <a href="https://jakevdp.github.io/blog/2015/07/06/model-complexity-myth/" target="_blank" rel="noreferrer">blog post</a> on the model complexity myth.)</p><p>The nature of the penalty in the penalized least squares solution is such that the &quot;estimates&quot; are well-defined even when the covariance matrix of the random effects converges to a &quot;singular&quot; or &quot;boundary&quot; value. In other words, singularity of the covariance matrix for the random effects, which means that there are one or more directions in which there is no variability in the random effects, is different from singularity of the model matrix for the random effects, which would affect the ability to define uniquely these coefficients. The penalty term always provides a unique solution for the random-effects coefficients.</p>',14)),s("p",null,[e[19]||(e[19]=a("In addition to handling naturally occurring rank deficiency in the random effects, the regularization allows us to fit explicitly overparameterized random effects. For example, we can use ")),e[20]||(e[20]=s("code",null,"fulldummy",-1)),e[21]||(e[21]=a(" to fit both an intercept term and ")),s("mjx-container",c,[(n(),i("svg",k,e[13]||(e[13]=[s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D45B",d:"M21 287Q22 293 24 303T36 341T56 388T89 425T135 442Q171 442 195 424T225 390T231 369Q231 367 232 367L243 378Q304 442 382 442Q436 442 469 415T503 336T465 179T427 52Q427 26 444 26Q450 26 453 27Q482 32 505 65T540 145Q542 153 560 153Q580 153 580 145Q580 144 576 130Q568 101 554 73T508 17T439 -10Q392 -10 371 17T350 73Q350 92 386 193T423 345Q423 404 379 404H374Q288 404 229 303L222 291L189 157Q156 26 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 112 180T152 343Q153 348 153 366Q153 405 129 405Q91 405 66 305Q60 285 60 284Q58 278 41 278H27Q21 284 21 287Z",style:{"stroke-width":"3"}})])])],-1)]))),e[14]||(e[14]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"n")])],-1))]),e[22]||(e[22]=a(" indicator variables in the random effects for a categorical variable with ")),s("mjx-container",m,[(n(),i("svg",f,e[15]||(e[15]=[s("g",{stroke:"currentColor",fill:"currentColor","stroke-width":"0",transform:"scale(1,-1)"},[s("g",{"data-mml-node":"math"},[s("g",{"data-mml-node":"mi"},[s("path",{"data-c":"1D45B",d:"M21 287Q22 293 24 303T36 341T56 388T89 425T135 442Q171 442 195 424T225 390T231 369Q231 367 232 367L243 378Q304 442 382 442Q436 442 469 415T503 336T465 179T427 52Q427 26 444 26Q450 26 453 27Q482 32 505 65T540 145Q542 153 560 153Q580 153 580 145Q580 144 576 130Q568 101 554 73T508 17T439 -10Q392 -10 371 17T350 73Q350 92 386 193T423 345Q423 404 379 404H374Q288 404 229 303L222 291L189 157Q156 26 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 112 180T152 343Q153 348 153 366Q153 405 129 405Q91 405 66 305Q60 285 60 284Q58 278 41 278H27Q21 284 21 287Z",style:{"stroke-width":"3"}})])])],-1)]))),e[16]||(e[16]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"n")])],-1))]),e[23]||(e[23]=a(" levels instead of the usual ")),s("mjx-container",u,[(n(),i("svg",g,e[17]||(e[17]=[t('<g stroke="currentColor" fill="currentColor" stroke-width="0" transform="scale(1,-1)"><g data-mml-node="math"><g data-mml-node="mi"><path data-c="1D45B" d="M21 287Q22 293 24 303T36 341T56 388T89 425T135 442Q171 442 195 424T225 390T231 369Q231 367 232 367L243 378Q304 442 382 442Q436 442 469 415T503 336T465 179T427 52Q427 26 444 26Q450 26 453 27Q482 32 505 65T540 145Q542 153 560 153Q580 153 580 145Q580 144 576 130Q568 101 554 73T508 17T439 -10Q392 -10 371 17T350 73Q350 92 386 193T423 345Q423 404 379 404H374Q288 404 229 303L222 291L189 157Q156 26 151 16Q138 -11 108 -11Q95 -11 87 -5T76 7T74 17Q74 30 112 180T152 343Q153 348 153 366Q153 405 129 405Q91 405 66 305Q60 285 60 284Q58 278 41 278H27Q21 284 21 287Z" style="stroke-width:3;"></path></g><g data-mml-node="mo" transform="translate(822.2,0)"><path data-c="2212" d="M84 237T84 250T98 270H679Q694 262 694 250T679 230H98Q84 237 84 250Z" style="stroke-width:3;"></path></g><g data-mml-node="mn" transform="translate(1822.4,0)"><path data-c="31" d="M213 578L200 573Q186 568 160 563T102 556H83V602H102Q149 604 189 617T245 641T273 663Q275 666 285 666Q294 666 302 660V361L303 61Q310 54 315 52T339 48T401 46H427V0H416Q395 3 257 3Q121 3 100 0H88V46H114Q136 46 152 46T177 47T193 50T201 52T207 57T213 61V578Z" style="stroke-width:3;"></path></g></g></g>',1)]))),e[18]||(e[18]=s("mjx-assistive-mml",{unselectable:"on",display:"inline",style:{top:"0px",left:"0px",clip:"rect(1px, 1px, 1px, 1px)","-webkit-touch-callout":"none","-webkit-user-select":"none","-khtml-user-select":"none","-moz-user-select":"none","-ms-user-select":"none","user-select":"none",position:"absolute",padding:"1px 0px 0px 0px",border:"0px",display:"block",width:"auto",overflow:"hidden"}},[s("math",{xmlns:"http://www.w3.org/1998/Math/MathML"},[s("mi",null,"n"),s("mo",null,"−"),s("mn",null,"1")])],-1))]),e[24]||(e[24]=a(" contrasts."))]),e[27]||(e[27]=t(`<div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">kb07 </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> MixedModels</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">dataset</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:kb07</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">contrasts </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> Dict</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(var </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=&gt;</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> HelmertCoding</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">() </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">for</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> var </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">in</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:spkr</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:prec</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:load</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">))</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">fit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(MixedModel, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@formula</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rt_raw </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> spkr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> prec </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> load </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">subj) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">prec</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">item)), kb07; contrasts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">contrasts)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Linear mixed model fit by maximum likelihood</span></span>
<span class="line"><span> rt_raw ~ 1 + spkr + prec + load + spkr &amp; prec + spkr &amp; load + prec &amp; load + spkr &amp; prec &amp; load + (1 | subj) + (1 + prec | item)</span></span>
<span class="line"><span>    logLik   -2 logLik      AIC         AICc        BIC     </span></span>
<span class="line"><span> -14846.7248  29693.4496  29719.4496  29719.6547  29790.8120</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Variance components:</span></span>
<span class="line"><span>             Column      Variance  Std.Dev.  Corr.</span></span>
<span class="line"><span>item     (Intercept)     158695.303 398.366</span></span>
<span class="line"><span>         prec: maintain   93637.345 306.002 -0.81</span></span>
<span class="line"><span>subj     (Intercept)     105909.621 325.438</span></span>
<span class="line"><span>Residual                 842822.907 918.054</span></span>
<span class="line"><span> Number of obs: 1789; levels of grouping factors: 32, 56</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Fixed-effects parameters:</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span>                                            Coef.  Std. Error      z  Pr(&gt;|z|)</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span>(Intercept)                             2225.71       85.5665  26.01    &lt;1e-99</span></span>
<span class="line"><span>spkr: old                                 74.1916     21.7062   3.42    0.0006</span></span>
<span class="line"><span>prec: maintain                          -377.212      58.2866  -6.47    &lt;1e-10</span></span>
<span class="line"><span>load: yes                                101.958      21.7062   4.70    &lt;1e-05</span></span>
<span class="line"><span>spkr: old &amp; prec: maintain               -28.0275     21.7062  -1.29    0.1966</span></span>
<span class="line"><span>spkr: old &amp; load: yes                     26.8642     21.7062   1.24    0.2159</span></span>
<span class="line"><span>prec: maintain &amp; load: yes               -18.6514     21.7062  -0.86    0.3902</span></span>
<span class="line"><span>spkr: old &amp; prec: maintain &amp; load: yes    15.4985     21.7062   0.71    0.4752</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">fit</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(MixedModel, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">@formula</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(rt_raw </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">~</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> spkr </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> prec </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">*</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> load </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">subj) </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> (</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">+</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">fulldummy</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(prec)</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">|</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">item)), kb07; contrasts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">contrasts)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>Linear mixed model fit by maximum likelihood</span></span>
<span class="line"><span> rt_raw ~ 1 + spkr + prec + load + spkr &amp; prec + spkr &amp; load + prec &amp; load + spkr &amp; prec &amp; load + (1 | subj) + (1 + prec | item)</span></span>
<span class="line"><span>    logLik   -2 logLik      AIC         AICc        BIC     </span></span>
<span class="line"><span> -14846.7248  29693.4496  29725.4496  29725.7566  29813.2802</span></span>
<span class="line"><span></span></span>
<span class="line"><span>Variance components:</span></span>
<span class="line"><span>             Column       Variance  Std.Dev.   Corr.</span></span>
<span class="line"><span>item     (Intercept)     1225798.517 1107.158</span></span>
<span class="line"><span>         prec: break      493552.854  702.533 -0.82</span></span>
<span class="line"><span>         prec: maintain  1006118.560 1003.055 -0.98 +0.80</span></span>
<span class="line"><span>subj     (Intercept)      105911.863  325.441</span></span>
<span class="line"><span>Residual                  842822.320  918.054</span></span>
<span class="line"><span> Number of obs: 1789; levels of grouping factors: 32, 56</span></span>
<span class="line"><span></span></span>
<span class="line"><span>  Fixed-effects parameters:</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span>                                            Coef.  Std. Error      z  Pr(&gt;|z|)</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span>
<span class="line"><span>(Intercept)                             2225.71       85.567   26.01    &lt;1e-99</span></span>
<span class="line"><span>spkr: old                                 74.1916     21.7062   3.42    0.0006</span></span>
<span class="line"><span>prec: maintain                          -377.212      58.287   -6.47    &lt;1e-10</span></span>
<span class="line"><span>load: yes                                101.958      21.7062   4.70    &lt;1e-05</span></span>
<span class="line"><span>spkr: old &amp; prec: maintain               -28.0275     21.7062  -1.29    0.1966</span></span>
<span class="line"><span>spkr: old &amp; load: yes                     26.8642     21.7062   1.24    0.2159</span></span>
<span class="line"><span>prec: maintain &amp; load: yes               -18.6514     21.7062  -0.86    0.3902</span></span>
<span class="line"><span>spkr: old &amp; prec: maintain &amp; load: yes    15.4985     21.7062   0.71    0.4752</span></span>
<span class="line"><span>──────────────────────────────────────────────────────────────────────────────</span></span></code></pre></div><p>This may be useful when the <code>PCA</code> property suggests a random effects structure larger than only main effects but smaller than all interaction terms. This is also similar to the functionality provided by <code>dummy</code> in <code>lme4</code>, but as in the difference between <code>zerocorr</code> in Julia and <code>||</code> in R, there are subtle differences in how this expansion interacts with other terms in the random effects.</p>`,5))])}const C=l(o,[["render",y]]);export{Q as __pageData,C as default};
