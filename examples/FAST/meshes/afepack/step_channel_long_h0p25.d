#------------------------------------------------------#
# Long L-shaped step channel for Navier-Stokes Ex. 7.4 #
# Domain: [-1,10]x[-1,1] without [-1,0]x[-1,0]        #
# Mesh size target: h = 0.25                           #
#------------------------------------------------------#

#=========
| POINTS |
=========#
6 # number of points #

# Nodes which define the boundary #
0:  -1   0   .250   1
1:   0   0   .250   2
2:   0  -1   .250   3
3:  10  -1   .250   4
4:  10   1   .250   5
5:  -1   1   .250   6

#===========
| SEGMENTS |
===========#
6 # Number of segments #

# Boundary segments #
0:   0   1   2
1:   1   2   2
2:   2   3   2
3:   3   4   5
4:   4   5   2
5:   5   0   1
