# Submission of HRM 1.3.0

This is a resubmission of a package that was archived on CRAN on 2024-04-20
("issues were not corrected despite reminders"). The last version on CRAN was
1.2.1 (2020-02-06).

## What caused the archival, and how it is fixed

The archived version depended on packages and settings that are no longer
available or accepted:

* `Suggests: RGtk2, cairoDevice` — 'RGtk2' was orphaned and archived from
  CRAN on 2021-12-15. The graphical user interface that used it has been
  removed from the package entirely, along with `Imports: tcltk`.
* `SystemRequirements: C++11` — removed; the specification is no longer
  needed.

While preparing this release I also removed an untested multivariate branch,
fixed `confint()` (it aborted for every design with a whole-plot factor),
corrected 16 help pages with documented arguments that did not exist, and
reduced the dependencies. See NEWS.md for the full list, including the
breaking changes.

## Test environments

* local macOS 26.5 (aarch64), R 4.4.2
* GitHub Actions, ubuntu-latest, R oldrel-1

<!-- TODO before submitting: run these and record the results here.
     Neither has been run yet, so do not submit with this list as it stands.
* win-builder, R-devel and R-release
* R-hub: linux, windows, macos
-->

## R CMD check results

0 errors | 0 warnings | 1 note

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Martin Happ <statistics@happ.co.at>'

    New submission

    Package was archived on CRAN

    CRAN repository db overrides:
      X-CRAN-Comment: Archived on 2024-04-20 as issues were not corrected
        despite reminders.

This note is expected for a package returning from the archive; the issues
that led to the archival are described above.

## Downstream dependencies

There are no reverse dependencies on CRAN.
