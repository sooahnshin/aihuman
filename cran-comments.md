## Current Submission
This is the submission of version 1.0.1, which includes the following updates:

* Updated a function (`compute_bounds_aipw()`, to allow for a special case) implementing the methods proposed in the accompanying paper.

## Resubmission
This is the submission of version 1.0.0, which includes the following significant updates:

* Added new functions implementing the methods proposed in the accompanying paper (listed as the last reference in the CITATION file).
* Added a new vignette to demonstrate the functionality of the new methods and reflect recent changes in the package.

## Resubmission
This is a third submission. In this version I have

* Made ensure that I do not comment out some code lines in examples.
* Changed all the information messages using `cat()` to `message()`.

## Resubmission
This is a second submission. In this version I have

* Used `\donttest{}` (instead of `\dontrun{}`) for some tests listed under examples because the algorithm is (by design) computationally expensive. 
* Made ensure that I do not use more than 2 cores in the examples, vignettes, etc.

## Test environments
* local R installation (macOS), R 4.3.1
* macOS-latest (on GitHub Actions), (release)
* windows-latest (on GitHub Actions), (release)
* ubuntu-latest (on GitHub Actions), (release)
* ubuntu-latest (on GitHub Actions), (old release)
* ubuntu-latest (on GitHub Actions), (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse Dependencies
There are no reverse dependencies to check.