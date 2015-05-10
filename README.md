# GNLHD
The package **GNLHD**. This package provides constructions for LHD, NLHD and GNLHD in R. The package  mainly includes
construction for GNLHD, which is a kind of new design originated from LHD, and two types of algrorithms to optimizing
GNLHDs. Currently, the R package can be only run on **Windows** platform. More details will be found in a coming paper,
GNLHD And Its Optimization written by Daijun Chen and Shifeng Xiong.
##Installation
You can also install the development version from Github, which provides daily build of
GNLHD:

```s
install.packages("devtools")
library(devtools)
install.packages("numbers")
library(numbers)
install_github("DavidDJChen/GNLHD")
```

If you know GIT and R CMD build, it's handly way:

```bash
git clone http://github.com/DavidDJChen/GNLHD.git
R CMD BUILD GNLHD
R CMD INSTALL GNLHD_*.tar.gz
```

#Usage

```s
library(GNLHD)

# Construct a GNLHD class u1 with stucture s=(3,5) and q=2.
u1<-GNLHD$new(s=c(3,5),q=2)

# Show the GNLHD class u1.
u1
<GNLHD>
  Public:
    diagnose: function
    GNLH: function
    GNLH_Full: function
    GNLH_illegal_set: function
    GNLH_permutation: function
    initialize: function
    Lcm: 15
    q: 2
    s: 3 5
    StandGNLHD: function
    Swap: function
    t: 5 3

# Generate two generalized nested permutations with structure s=(3,5)
u1$GNLH_Full


```



