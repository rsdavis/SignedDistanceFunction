# SignedDistanceFunction
C++ class for constructing the signed distance function from a level set function.  

When working with level set functions, it is sometimes necessary to know the minimal distance to the zero level set at any point withing the domain. This is the purpose of the signed distance function. One can take any level set function (with a zero level set) and convert it into a signed distance function with the same zero set, but with distance information built into the contours.

This repository consist of a single C++ header file that contains the necessary routines for calculating the signed distance function from a given level set function. It uses an iterative scheme to converge to the level set function to within a specified tolerance.

The signed distance function can be helpful in image analysis, data modifications, and level set methods that are used in material science. More information can be found at my [personal website](http://rsdavis.mycpanel.princeton.edu/wp/?p=24).
