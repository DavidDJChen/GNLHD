# GNLHD
The package **GNLHD**. This package provides constructions for LHD, NLHD and GNLHD in R. The package  mainly includes
construction for GNLHD, which is a kind of new design originated from LHD, and two types of algrorithms to optimizing
GNLHDs. Currently, the R package can be only run on **Windows** platform. More details will be found in a coming paper,
Generalized Nested Latin Hypercube Design And Its Optimization written by Daijun Chen and Shifeng Xiong.
##Installation
You can also install the development version from Github, which provides daily build of
GNLHD:

```s
>install.packages("devtools")
>library(devtools)
>install.packages("numbers")
>library(numbers)
>install_github("DavidDJChen/GNLHD")
```

If you know GIT and R CMD build, it's handly way:

```bash
git clone http://github.com/DavidDJChen/GNLHD.git
R CMD BUILD GNLHD
R CMD INSTALL GNLHD_*.tar.gz
```

#Usage
Part 1: Construction for GNLHD and some useful functions
```s
>library(GNLHD)

# Construct a GNLHD class u1 with stucture s=(3,5) and q=2.
>u1<-GNLHD$new(s=c(3,5),q=2)

# Show the GNLHD class u1.
>u1
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

# Generate a generalzied nested permutation with structure s=(3,5).
> u1$GNLH_permutation()
 [1] 10 14  2  9  4  7 15  6  1 13  3  5  8 11 12

# Generate q generalized nested permutations with structure s=(3,5).
>u1$GNLH_Full
      [,1] [,2]
 [1,]    1    4
 [2,]   14    9
 [3,]    8   13
 [4,]   10    1
 [5,]    6   11
 [6,]   15   15
 [7,]    2   12
 [8,]    7   10
 [9,]    5    8
[10,]   12    3
[11,]    9    6
[12,]   13    2
[13,]   11    5
[14,]    4    7
[15,]    3   14

# Generate a GNLH structure s=(3,5) and q=2.
> u1$GNLH()
     [,1] [,2]
[1,]    4    4
[2,]   13    8
[3,]    9   15
[4,]   10    1
[5,]    3   10

# Generate a Standard GNLHD with s=(3,5) and q=2.
>u1$StandGNLHD()
          [,1]      [,2]
[1,] 0.2333333 0.4333333
[2,] 0.5666667 0.3000000
[3,] 0.7666667 0.7666667
[4,] 0.9000000 0.1000000
[5,] 0.1000000 0.8333333

# Do layer-in swap operation on a in the 2nd layer and  the 1st column.
> a<-u1$GNLH_Full()
> a
      [,1] [,2]
 [1,]   12   12
 [2,]    2    3
 [3,]    9    9
 [4,]   14    6
 [5,]    5   14
 [6,]    1    1
 [7,]   10    7
 [8,]   15   11
 [9,]    7    4
[10,]    6   15
[11,]    8    5
[12,]    4    8
[13,]    3   13
[14,]   13    2
[15,]   11   10
> u1$Swap(a, structure=c(3,5),column=1, Swap_type="in", Swap_layer=2,lcm=15) 
      [,1] [,2]
 [1,]   12   12
 [2,]    2    3
 [3,]    9    9
 [4,]    5    6
 [5,]   14   14
 [6,]    1    1
 [7,]   10    7
 [8,]   15   11
 [9,]    7    4
[10,]    6   15
[11,]    8    5
[12,]    4    8
[13,]    3   13
[14,]   13    2
[15,]   11   10

# Do  layer-between swap operation on a in the 2nd layer and  the 1st column.
> u1$Swap(a, structure=c(3,5),column=1, Swap_type="between", Swap_layer=2,lcm=15) 
      [,1] [,2]
 [1,]   12   12
 [2,]    2    3
 [3,]    9    9
 [4,]   14    6
 [5,]    6   14
 [6,]    1    1
 [7,]   10    7
 [8,]   15   11
 [9,]    7    4
[10,]    5   15
[11,]    8    5
[12,]    4    8
[13,]    3   13
[14,]   13    2
[15,]   11   10
```
Part 2: Calucate the Phi_p value with p=50 and t=2 for a.

```s
> Phi_p(Design=a,t=2,p=50)
[1] 0.7169776
```

Part 3: Utilizing an efficient Sequential algorithm to get a optimal GNLHD

```s
# Utilizing an efficient sequential algorithm with iteration=1000 to get a optimal GNLHD
from a.
> Optimal_GNLHD_SequentialAlg(GNLHD=u1, GNLH_Full=a, iteration=1000,T_h_initial=0.1,M=100,J=6,t=2,p=50,
                                     tolerance=0.1,alpha=c(0.8,0.9,0.7))
      [,1] [,2]
 [1,]   12   12
 [2,]    2    3
 [3,]    9    9
 [4,]   15    4
 [5,]    4   15
 [6,]    1    1
 [7,]   10    7
 [8,]   14   11
 [9,]    7    5
[10,]    6   14
[11,]    8    6
[12,]    5    8
[13,]    3   13
[14,]   13    2
[15,]   11   10

```
## License

This package is free and open source software, licensed under GPL-3.

## References

```s
Daijun Chen and Shifeng Xiong, Generalized Nested Latin Hypercube Design and Its Optimization.
```

## TODO

*1. Deploy on CRAN:

You can install the advanced version on [CRAN](http://cran.r-project.org/package=GNLHD):

```s
install.packages('GNLHD')
```

*2. Deploy on RForge:

You can also install from [RForge](http://rforge.net/GNLHD/),which provides daily build of **GNLHD**;

```s

# update all existing packages
update.packages(ask=FALSE, repos='http://cran.rstudio.org')
install.packages('GNLHD', repos=c('http://rforge.net','http://cran.rstudio.org'),type='source')









