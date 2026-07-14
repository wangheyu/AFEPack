#-------------------------------------------------------------#
# Contracting-channel proxy for book Exercise 1.5 / 5.7       #
# Domain: 0 <= x <= 4, |y| <= 1 - x/8                         #
# Mesh size target: h = 0.25                                  #
#-------------------------------------------------------------#

#=========
| POINTS |
=========#
4 # number of points #

# Nodes which define the boundary #
0:   0  -1.0   .250   1
1:   4  -0.5   .250   2
2:   4   0.5   .250   3
3:   0   1.0   .250   4

#===========
| SEGMENTS |
===========#
4 # Number of segments #

# Boundary segments: inflow=1, walls=2, natural outflow=5 #
0:   0   1   2
1:   1   2   5
2:   2   3   2
3:   3   0   1
