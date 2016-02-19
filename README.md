##SPLINTER
SPLINTER (SPLine INTERpolation) is a library for *multivariate function approximation* implemented in C++. The library can be used for function approximation, regression, data smoothing, data reduction, and much more. The library contains the following approximation methods:

- [tensor product B-splines](docs/bspline.md)
- [radial basis function networks](docs/rbfnetwork.md)
- [polynomial regression](docs/polynomial_regression.md)

A shared feature of these methods is that they are based on models that are linear in the coefficients. The models differ in which type of basis functions that is used. For example, the B-spline uses piecewise polynomial basis functions while the radial basis function network may use various radial basis functions. The coefficients in these models are computed using ordinary least squares (OLS), possibly with Tikhonov regularization (regularization is partially supported today). The name of the library, SPLINTER, originates from the tensor product B-spline implementation, which was the first of the methods to be implemented.

![Illustration of a B-spline](assets/bspline.png)
Figure: Illustration of a bicubic B-spline generated with the SPLINTER library.

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
3. To implement new B-spline features

The simplest way to contribute to SPLINTER is to use it and give us feedback on the experience. If you would like to contribute by coding, you can get started by picking a suitable issue from the [list of issues](https://github.com/bgrimstad/splinter/issues). The issues are labeled with the type of work (`Bug`, `Docs`, `Enhancement`, `New feature`, `Refactoring`, `Tests`) and level of difficulty (`Beginner`, `Intermediate`, `Advanced`). Some issues are also labeled as `Critical`, which means that they deserve our attention and prioritization.

###Requirements for use
A standards compliant C++11 compiler.

###Guides
* [Basic usage](docs/basic_usage.md)
* [Compilation](docs/compile.md)
* [MatLab interface](docs/MATLAB.md)
* [Python interface](docs/python.md)
* [C interface](docs/C.md)
