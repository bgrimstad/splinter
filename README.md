##SPLINTER
SPLINTER (SPLine INTERpolation) is a library for multivariate function approximation implemented in C++. The library can be used for function approximation, regression, data smoothing, and much more. Currently, the library contains the following implementations:

1. [tensor product B-splines](http://en.wikipedia.org/wiki/B-spline), 
2. [radial basis functions](http://en.wikipedia.org/wiki/Radial_basis_function), including the [thin plate spline](http://en.wikipedia.org/wiki/Thin_plate_spline), and
3. [polynomial regression](http://en.wikipedia.org/wiki/Polynomial_regression).

The coefficients in these models are computed using ordinary least squares (OLS). The name of the library, SPLINTER, originates from the tensor product B-spline implementation, which was the first of the methods to be implemented.

The B-spline may approximate any sampled multivariate function. The user may construct a linear (degree 1), quadratic (degree 2), cubic (degree 3) or higher degree B-spline that smoothes or interpolates the data. The B-spline is constructed from the samples by solving a linear system. On a modern desktop computer the practical limit on the number of samples is about 100 000 when constructing a B-spline. This translates to a practical limit of 6-7 variables. Evaluation time, however, is independent of the number of samples due to the local support property of B-splines. That is, only samples neighbouring the evaluation point affect the B-spline value. Evaluation do however scale with the degree and number of variables of the B-spline.

The user may create a penalized B-spline (P-spline) that smooths the data instead of interpolating it. The construction of a P-spline is more computationally demanding than the B-spline - a large least-square problem must be solved - bringing the practical limit on the number of samples down to about 10 000.

When sampling is expensive and/or scattered (not on a grid) a radial basis function may be utilized for function approximation. The user should expect a high computational cost for constructing and evaluating a radial basis function spline, even with a modest number of samples (up to about 1 000 samples). 

![Illustration of a B-spline](assets/bspline.png)
Figure: Illustration of a cubic B-spline generated with the SPLINTER library.

###Sharing
SPLINTER is the result of several years of development towards a fast and general library for multivariate function approximation. The initial intention with the library was to build splines for use in mathematical programming (nonlinear optimization). Thus, some effort has been put into functionality that supports this, e.g. Jacobian and Hessian computations for the B-spline. 

By making SPLINTER publicly available we hope to help anyone looking for a multivariate function approximation library. In return, we expect nothing but your suggestions, improvements, and feature requests. If you use SPLINTER in a scientific work we kindly ask you to cite it. You can cite it as shown in the bibtex entry below (remember to update the date accessed).
```
@misc{SPLINTER,
  title={{SPLINTER: a library for multivariate function approximation}},
  author={Bjarne Grimstad and others},
  howpublished={\url{http://github.com/bgrimstad/splinter}},
  year={2015},
  note={Accessed: 2015-05-16}
}
```
###Contributing
Everyone is welcome to use and contribute to SPLINTER. We believe that collective effort over time is the only way to create a great library: one that makes multivariate function approximation more accessible to the programming community.

The current goals with the library are:

1. To make the library more accessible by improving the interfaces and documentation
2. To improve the current code via testing
3. To implement new function approximation methods such as Kriging and k-nearest neighbors regression

The simplest way to contribute to SPLINTER is to use it and give us feedback on the experience. If you would like to contribute by coding, you can get started by picking a suitable issue from the [list of issues](https://github.com/bgrimstad/splinter/issues). The issues are labeled with the type of work (`Bug`, `Docs`, `Enhancement`, `New feature`, `Refactoring`, `Tests`) and level of difficulty (`Beginner`, `Intermediate`, `Advanced`). Some issues are also labeled as `Critical`, which means that they deserve our attention and prioritization.

###Requirements for use
A standards compliant C++11 compiler.

###Guides
* [Basic usage](docs/basic_usage.md)
* [Compilation](docs/compile.md)
* [MatLab interface](docs/MATLAB.md)
* [Python interface](docs/python.md)
* [C interface](docs/C.md)
