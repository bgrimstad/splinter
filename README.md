##Multivariate Splines
Multivariate Splines is a function approximation library implementing various multivariate splines in C++. It contains the following implementations:

1. a speedy implementation of the tensor product [B-spline](http://en.wikipedia.org/wiki/B-spline), and 
2. a simple implementation of [radial basis function splines](http://en.wikipedia.org/wiki/Radial_basis_function), including the [thin plate spline](http://en.wikipedia.org/wiki/Thin_plate_spline).

The B-spline may approximate any multivariate function sampled on a grid. The user may construct a linear (degree 1) or cubic (degree 3) spline that interpolates the data. The B-spline is constructed from the samples by solving a linear system. A modern desktop computer limits the number of samples to about 100 000 when constructing a B-spline - evaluation of the spline is very fast due to the local support property of B-splines. 

The user may create a penalized B-spline (P-spline) that smooths the data instead of interpolating it. The construction of a P-spline is more computationally demanding than the B-spline - a large least-square problem must be solved - bringing the limit on the number of samples down to about 10 000.

When sampling is expensive and/or scattered (not on a grid) the radial basis function splines may be utilized for function approximation. The user should expect a high computational cost for constructing and evaluating a radial basis function spline, even with a modest number of samples (up to about 1 000 samples). 

The library is based on the C++ template linear algebra library [Eigen](http://eigen.tuxfamily.org); its sparse matrix support is particularly important for the speed of the tensor product B-spline implementation.

###Author's note:
Multivariate Splines implements various splines for function approximation with the purpose of utilizing the splines/approximations in mathematical programming (nonlinear optimization). Thus, special attention has been given to functionality that may support a nonlinear optimization solver. For example, the B-spline implementation includes evaluation of the Jacobian and Hessian.

NOTE: focus has not been on curve fitting, NURBS, etc.

NOTE: general implementation which is readily extended with new functionality. The author welcomes new contributions and improvements to the code. 

NOTE: the goal is to create an open, general, and fast library for multivariate splines.

###Requirements for use:
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)


###Guides
* [Installation](docs/install.md)
* [Basic usage](docs/basic_usage.md)
