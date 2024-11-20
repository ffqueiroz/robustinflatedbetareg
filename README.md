
# robustinflatedbetareg

<!-- badges: start -->
<!-- badges: end -->

The function fit two robust estimators for inflated beta regression models (M-LSE and M-LME).
The tuning constants are selected via the data-driven algorithm. The function used to fit is <tt>fit.betainflated()</tt>.

```r
fit.betainflated(y, S, X, Z, linkmu="logit", linkphi="log", linkvartheta = "logit")
```

The script auxiliaryfunctions.R contains auxiliary functions used in <tt>fit.betainflated()</tt> 
function. In addition, the functions <tt>residuals_ZIBE()</tt> and <tt>envelope()</tt> provides
residuals and envelope plots, respectively, for fits of the robust inflated beta regression
models. 