# R package development documentation



Go to for R package dev best practice : https://r-pkgs.org/ by [Hadley Wickham](http://hadley.nz/) and [Jenny Bryan](http://jennybryan.org/)



### Package creation

two ways to create the package.

- Call `usethis::create_package()`.

- In RStudio, do *File > New Project > New Directory > R Package*. This ultimately calls `usethis::create_package()`, so really there’s just one way.

  

### Devtools R package development functions

- [use_r()](https://usethis.r-lib.org/reference/use_r.html) creates and/or opens a script in R/.
- [load_all()](https://devtools.r-lib.org/reference/load_all.html) this is like a test drive. “lather, rinse, repeat” cycle of package development. It simulates the process of building, installing during interactive development. As your package accumulates more functions, some exported, some not, some of which call each other, some of which call functions from packages you depend on, [load_all()](https://devtools.r-lib.org/reference/load_all.html) gives you a much more accurate sense of how the package is developing than test driving functions defined in the global environment
- [check()](https://devtools.r-lib.org/reference/check.html) :incremental development of .R and .Rmd file, check to see if all the moving parts are working
- [use_mit_license()](https://usethis.r-lib.org/reference/licenses.html) also puts a copy of the full license in LICENSE.md and adds this file to .Rbuildignore. It’s considered a best practice to include a full license in your package’s source, such as on GitHub, but CRAN disallows the inclusion of this file in a package tarball.
- [document()](https://devtools.r-lib.org/reference/document.html)  add R documentation file, man/genePanel.Rd. call to [document()](https://devtools.r-lib.org/reference/document.html) updates the NAMESPACE file, based on @export tags found in roxygen comments.
- [install()](https://devtools.r-lib.org/reference/install.html) when no error from check() re-install
- [use_testthat()](https://usethis.r-lib.org/reference/use_testthat.html) : formalize install of package by as a unit test.This initializes the unit testing machinery for your package. adds Suggests: testthat to DESCRIPTION, creates the directory tests/testthat/, and adds the script tests/testthat.R. You’ll notice that testthat is probably added with a minimum version of 3.0.0 and a second DESCRIPTION field, Config/testthat/edition
- [use_package](https://usethis.r-lib.org/reference/use_package.html) Adding 'other package to Imports field in DESCRIPTION
- [use_github()](https://usethis.r-lib.org/reference/use_github.html) connect to git repository
- [use_readme_rmd()](https://usethis.r-lib.org/reference/use_readme_rmd.html) function initializes a basic, executable README.Rmd ready for you to edit
- [use_devtools()](https://usethis.r-lib.org/reference/rprofile-helper.html) creates .Rprofile adding below in this repeatedly attach devtools in every R session

```
if (interactive()) {
  suppressMessages(require(devtools))
}
```

- devtools::dev_sitrep() :If this reveals that certain tools or packages are missing or out-of-date, you are encouraged to update them.
- devtools::build() to make bundled package "source tarballs" for easy export between platforms
- to exclude a specific file or directory is to use `usethis::use_build_ignore("notes")`. this is quivalent to adding it in `.Rbuildignore`





### Name of the package requirements

Formal requirements:

There are three formal requirements:

1. The name can only consist of letters, numbers, and periods, i.e., `.`.
2. It must start with a letter.
3. It cannot end with a period.

Unfortunately, this means you can’t use either hyphens or underscores, i.e., `-` or `_`, in your package name. We recommend against using periods in package  names, due to confusing associations with file extensions and S3  methods.



Pragmatic advice:

- Pick a unique name that’s easy to Google. This makes it easy for  potential users to find your package (and associated resources) and for  you to see who’s using it.

- Don’t pick a name that’s already in use on CRAN or Bioconductor. You  may also want to consider some other types of name collision:

  - Is there an in-development package maturing on, say, GitHub that  already has some history and seems to be heading towards release?
  - Is this name already used for another piece of software or for a  library or framework in, e.g., the Python or JavaScript ecosystem?

- Avoid using both upper and lower case letters: doing so makes the package name hard to type and even harder to remember. For example,  it’s hard to remember if it’s Rgtk2 or RGTK2 or RGtk2.

- Give preference to names that are pronounceable, so people are  comfortable talking about your package and have a way to hear it inside  their head.

- Find a word that evokes the problem and modify it so that it’s unique:

  - lubridate makes dates and times easier.
  - rvest “harvests” the content from web pages.
  - r2d3 provides utilities for working with D3 visualisations.
  - forcats is an anagram of factors, which we use **for** **cat**egorical data.

- Use abbreviations:

  - Rcpp = R + C++ (plus plus)
  - brms = Bayesian Regression Models using Stan

- Add an extra R:

  - stringr provides string tools.
  - beepr plays notification sounds.
  - callr calls R, from R.

- Don’t get sued.

  - If you’re creating a package that talks to a commercial service,  check the branding guidelines. For example, rDrop isn’t called rDropbox  because Dropbox prohibits any applications from using the full  trademarked name.

  Nick Tierney presents a fun typology of package names in his [Naming Things](https://www.njtierney.com/post/2018/06/20/naming-things/) blog post; see that for more inspiring examples. He also has some experience with renaming packages, so the post [So, you’ve decided to change your r package name](https://www.njtierney.com/post/2017/10/27/change-pkg-name/) is a good resource if you don’t get this right the first time.





The [available package](https://cran.r-project.org/package=available) has a function called `available()` that helps you evaluate a potential package name from many angles:

- Checks for validity.
- Checks availability on CRAN, Bioconductor, and beyond.
- Searches various websites to help you discover any unintended  meanings. In an interactive session, the URLs you see above are opened  in browser tabs.
- Attempts to report whether name has positive or negative sentiment.

```
library(available)

available("doofus")
#> Urban Dictionary can contain potentially offensive results,
#>   should they be included? [Y]es / [N]o:
#> 1: 1
#> ── doofus ──────────────────────────────────────────────────────────────────
#> Name valid: ✔
#> Available on CRAN: ✔ 
#> Available on Bioconductor: ✔
#> Available on GitHub:  ✔ 
#> Abbreviations: http://www.abbreviations.com/doofus
#> Wikipedia: https://en.wikipedia.org/wiki/doofus
#> Wiktionary: https://en.wiktionary.org/wiki/doofus
#> Sentiment:???
```





