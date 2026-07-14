#-----------------------------------------------------------------#
# Original IFISS quad_domain half bow-tie for Exercises 5.7 / 7.6 #
# Domain: 0 <= x <= 4, 0 <= y <= 1 - 3x/16                       #
# Mesh size target: h = 0.125                                     #
#-----------------------------------------------------------------#

#=========
| POINTS |
=========#
4 # number of points #

0:   0   0.000   .125   1
1:   4   0.000   .125   2
2:   4   0.250   .125   3
3:   0   1.000   .125   4

#===========
| SEGMENTS |
===========#
4 # Number of segments #

# Boundary segments: inflow=1, walls=2, natural outflow=5 #
0:   0   1   2
1:   1   2   5
2:   2   3   2
3:   3   0   1
