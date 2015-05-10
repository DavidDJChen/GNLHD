# GNLHD
The package GNLHD. This package provides constructions for LHD, NLHD and GNLHD in R. The package  mainly includes
construction for GNLHD, which is a kind of new design originated from LHD, and two types of algrorithms to optimizing
GNLHDs. Currently, the R package can be only run on Windows platform. More details will be found in a coming paper,
GNLHD And Its Optimization written by Daijun Chen and Shifeng Xiong.
#Installation
You can also install the development version from Github, which provides daily build of
GNLHD:
install.packages("devtools")

library(devtools)

install.packages("numbers")

library(numbers)

install_github("DavidDJChen/GNLHD")

if you know GIT and R CMD build, it's handly way:

git clone http://github.com/DavidDJChen/GNLHD.git

R CMD BUILD GNLHD

R CMD INSTALL GNLHD_*.tar.gz
#Usage
