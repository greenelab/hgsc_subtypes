#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("doppelgangR", excludePattern="AllClasses|nonexports|sn")
RSCRIPT

## Note: After doing ./runinline.sh, please remove the line \alias{doppelgangR} from man/doppelgangR-package.Rd.
grep -v '\\alias{doppelgangR}' doppelgangR/man/doppelgangR-package.Rd > doppelgangR/man/doppelgangR-package.Rd_tmp
mv doppelgangR/man/doppelgangR-package.Rd_tmp doppelgangR/man/doppelgangR-package.Rd

R CMD build doppelgangR

