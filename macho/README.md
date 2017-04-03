# White Paper

The white paper for the MACHO dark matter portion of the hack is [here](WhitePaper.md).

# Summary of Proposed Hack

The ultimate objective of this portion of the hack it to create Metric Analysis Framework metrics to determine the value of a given opsim run for the intermediate mass MACHO science. 

It will build off the many existing maf examples contained in the [sims_maf_contrib repo](https://github.com/LSST-nonproject/sims_maf_contrib/blob/master/tutorials/Index.ipynb).

# Dependencies and Suggested Installations

## Anaconda

This is optional but will make various LSST software installations easier.

## opsim run
Download an opsim catalog from:

https://www.lsst.org/scientists/simulations/opsim/opsim-survey-data

For this hack we will use the latest: [minion_1016](http://ops2.lsst.org/runs/reference_run/minion_1016/minion_1016_sqlite.db.gz)

This sqlite database is ~2 GB and took me approximately 30 minutes to download.

## sims_maf_contrib

The bulk of the hack will likely build off the many tutorials in the sims_maf_contrib repo. I created a fork of this repo to contain the work for this hack. It is recommended that anyone working on this portion of the hack clone this fork of the repo:

https://github.com/wadawson/sims_maf_contrib

## MAF

To install MAF follow the directions at:

https://github.com/wadawson/sims_maf_contrib/blob/master/tutorials/Index.ipynb

## Install sqlite

With homebrew this can be done with:

> brew install sqlite3

## Install DB Browser for sqlite

http://sqlitebrowser.org

This can be useful if you want to hack the opsim database more than is convenient with with MAF. For example, if you want to make a triangle plot.
