# HRM 1.3.0

First release after the package was archived on CRAN on 2024-04-20. The
issues that led to the archival are addressed below.

## Breaking changes

* The graphical user interface has been removed: `hrm_GUI()` and the
  deprecated `hrm.GUI()` are gone. It was built on 'RGtk2', which was
  orphaned and archived from CRAN on 2021-12-15, and on 'tcltk'. Both
  dependencies have been dropped.

* Multivariate repeated measures have been removed. The `variable` argument
  of `hrm_test()` no longer exists, and passing it is an error. The feature
  shipped without any active test coverage; the code is retained outside the
  package for future work.

## Bug fixes

* `confint()` no longer fails for designs with a whole-plot factor. It
  previously aborted with `missing value where TRUE/FALSE needed`, which
  affected every such design, parametric and nonparametric alike, including
  the example given in the README. The cause is an upstream failure of
  `mvtnorm::qmvnorm()` on larger correlation matrices.

  `confint()` now tries progressively weaker forms of multiplicity control
  and reports which one was used in `attr(, "status")`:

  1. family-wise error rate over all factor level combinations,
  2. family-wise error rate within each whole-plot group,
  3. a Sidak correction within each group,
  4. none, in which case the interval columns are `NA`.

  Results for designs without a whole-plot factor are unchanged.

* Stale `@param` entries were removed from 16 help pages, where documented
  arguments did not exist in the corresponding functions.

* `plot()` on an unsupported design now raises an error instead of printing
  a message and returning invisibly.

## Improvements

* `confint()` is considerably faster. The correlation matrix was built with
  a double loop performing three matrix products per cell; `stats::cov2cor()`
  gives an identical result in a fraction of the time.

* Fewer dependencies. `xtable`, `tcltk`, `reshape2` and `RGtk2` are no
  longer used, and `MASS` moved to `Suggests`. `MASS`, `matrixcalc`, `plyr`
  and `ggplot2` moved from `Depends` to `Imports`, so loading HRM no longer
  attaches them to the search path. Code that relied on this must now
  attach them itself, for example `library(ggplot2)` before passing theme
  objects to `plot()`.

* R (>= 4.2.0) is now required, which reflects what the dependency tree can
  actually be installed on.

* Tests no longer assert the exact class vector of 'ggplot2' objects, which
  changed in 'ggplot2' 4.0.

* Test coverage added for `confint()` on whole-plot designs, which was
  previously untested entirely.

# HRM 1.2.1

* See the CRAN archive for earlier versions.
