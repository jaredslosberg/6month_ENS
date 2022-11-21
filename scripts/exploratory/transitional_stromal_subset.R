---
title: "Transitional_stromal_subset"
author: "Jared Slosberg"
date: "1/4/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Subset and process cells from GCA labeled as transitional stromal. Characterize their similarity to MENs, and if this is driven by a EMT-like signature

```{r load}
library(tidyverse)
library(monocle3)
library(here)

source(here("scripts/accessory_functions/monocle_mods"))
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
