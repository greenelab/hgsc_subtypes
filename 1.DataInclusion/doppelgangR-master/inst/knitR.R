library(knitr)
opts_chunk$set(cache = TRUE)
knit("../vignette//doppelgangR.Rnw")
system("pdflatex doppelgangR.tex")
