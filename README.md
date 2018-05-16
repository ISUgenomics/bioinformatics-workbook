# Bioinformatics Workbook

## Preface

The best way to learn bioinformatics is through examples of real world problems.  The Bioinformatics Workbook provides the reader with an in depth understanding of experimental design, data acquisition, data wrangling, data analysis and visualization.  This is accomplished through worked out example problem in each of these sections along with one or more advanced problem sets and corresponding solutions.  This books assumes that the reader has some knowledge of biology and basic understanding of the Unix command line.  However, for the beginner, the appendix contains introductory material and tips/tricks for common bioinformatic problems, that is referred to for more information throughout the book.

Please start your exploration at the [index](https://isugenomics.github.io/bioinformatics-workbook/).

## For Developer

To run the repo locally, go to the root directory of this repo

- install [bundler](https://bundler.io/) and [jekyll](https://jekyllrb.com/) by [RubyGems](http://rubygems.org/):

  ``` shell
    gem install bundler jekyll
  ```

- install the dependencies:

  ``` shell
      bundle install
  ```

- run the site:

  ``` shell
    bundle exec jekyll serve
  ```

- when developing, use the watch mode to automatically regenerate the site:

  ``` shell
    bundle exec jekyll serve --watch
  ```
