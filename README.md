**MegaLMM** is a package for fitting multi-trait linear mixed models (MvLMMs) with multiple random effects based on the paper available here: 
    
Please treat this as a Beta version and let me know of issues running the functions.    
    
### Installation:
```{r}
devtools::install_github('deruncie/MegaLMM')
```

### Introduction:
Please see the vignette `Running_MegaLMM.Rmd` for an introduction to using **MegaLMM**

If you would like to build the vignette (see below), do:
```{r}
devtools::install_github('deruncie/MegaLMM', build_opts = c("--no-resave-data", "--no-manual"),force = TRUE,build_vignettes = TRUE)
```

```{r}
vignette(topic = 'Running_MegaLMM',package='MegaLMM')
```

