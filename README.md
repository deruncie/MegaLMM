**MegaLMM** is a package for fitting multi-trait linear mixed models (MvLMMs) with multiple random effects based on the paper available here:

Please treat this as a Beta version and let me know of issues running the functions.

### Installation:

```{r}
devtools::install_github('deruncie/MegaLMM')
```

Note: the package requires OpenMP which causes problems on a Mac. I was able to get everything to run by following these instructions: <https://mac.r-project.org/openmp/> including the addition to `~/.R/Makevars`. Installing on Unix-like machines and maybe Windows should cause fewer problems.

### Introduction:

Please see the vignette [`MultiEnvironmentTrial.Rmd`](https://github.com/deruncie/MegaLMM/blob/master/vignettes/MultiEnvironmentTrial.Rmd) for an introduction to using **MegaLMM**. This will apply **MegaLMM** to the analysis of a very incomplete multi-environment trial, and go through many of the key functions in the package

You can simply download the .Rmd file above and run it yourself. If you would like to build the vignette, do:

```{r}
devtools::install_github('deruncie/MegaLMM', build_opts = c("--no-resave-data", "--no-manual"),force = TRUE,build_vignettes = TRUE)
```

```{r}
vignette(topic = 'MultiEnvironmentTrial',package='MegaLMM')
```

### Background

**MegaLMM** implements multivariate linear mixed models of the form:

    Y = XB + Z_1 U_1 + Z_2 U_2 + ... + E

where `Y` is a `n x t` matrix of observations for `n` individuals and `t` traits, `X` is a design matrix for `b` fixed effects (including an intercept), `Z_1` through `Z_m` are design matrices for `m` random effects, and `E` is a `n x t` matrix of residual errors. The random effects are `U_1` through `U_m` are mutually independent of eachother and the residuals, but columns of each `U_1` matrix can be correlated, and each column vector marginally follows a multivariate normal distribution with a known covariance matrix `K_m`.

MvLMMs like this are notoriously difficult to fit. We address this by re-paramterizing the MvLMM as a mixed effect factor model:

    Y = F Lambda + X B_r + Z_1 U_r1 + Z_2 U_r2 + ... + E_r
    F = X B_f + Z_1 U_f1 + Z_2 U_f2 + ... + E_f

where `F` is a `n x k` matrix of latent factor traits and `Lambda` is a `k x t` matrix of factor loadings. This is the model actually fit by **MegaLMM**. For full details, please see the accompanying paper.

![Model outline](misc/Conceptual_model_v4.png)

**Figure 1** from Runcie et al 2021.

The unique aspects of **MegaLMM** relative to other factor models are:

1.  The residuals of `Y` after accounting for the factors are not assumed to be iid, but are modeled with independent (across traits) LMMs accounting for both fixed and random effects.
2.  The factors themseleves are also not assumed to be iid, but are modeled with the same LMMs. This highlights the parallel belief that these latent factors represent traits that we just didn't measure directly.
3.  Each factor is shared by all modeled sources of variation (fixed effects, random effects and residuals), rather than being unique to a particular source.
4.  The factor loadings are strongly regularized so ensure that estimation is efficient. We accomplish this by ordering the factors from most-to-least important using a prior similar to that proposed by [Bhattarchya and Dunson (2011)](https://pubmed.ncbi.nlm.nih.gov/23049129/)
