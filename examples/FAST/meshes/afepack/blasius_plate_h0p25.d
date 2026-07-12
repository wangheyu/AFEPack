#------------------------------------------------------------#
# Thin-plate proxy for Navier-Stokes Blasius Example 7.1.4   #
# Domain: [-1,5]x[-1,1] minus [0,5]x[-0.04,0.04]            #
#------------------------------------------------------------#

#=========
| POINTS |
=========#
8 # number of points #

# Polygon boundary around a thin plate of half-thickness 0.04 #
0:  -1  -1     .750   1
1:   5  -1     .750   2
2:   5  -0.04  .250   3
3:   0  -0.04  .250   6
4:   0   0.04  .250   6
5:   5   0.04  .250   3
6:   5   1     .750   4
7:  -1   1     .750   1

#===========
| SEGMENTS |
===========#
8 # Number of segments #

# One simple closed polygon boundary #
0:   0   1   2
1:   1   2   3
2:   2   3   6
3:   3   4   6
4:   4   5   6
5:   5   6   3
6:   6   7   4
7:   7   0   1
