##Radial basis function network
The [radial basis function network](http://en.wikipedia.org/wiki/Radial_basis_function_network) is a widely used multivariate function approximation technique that inherently supports construction from scattered (unstructured) data points. This method is also popular when the function to be approximated is expensive to evaluate.

Arguably, a disadvantage with radial basis function networks is that they have global support (evaluation of the network requires an evaluation of all radial basis functions). Furthermore, the network is trained by solving a linear system that often becomes ill-conditioned. Thus, the user should expect a high computational cost for constructing and evaluating a radial basis function network, even with a modest number of samples (up to about 1 000 samples). 

SPLINTER implements the following [radial basis functions](http://en.wikipedia.org/wiki/Radial_basis_function):
- Gaussian
- Multiquadric
- Inverse quadratic
- Inverse multiquadric
- Thin plate spline

The following Python example shows how to build a normalized radial basis function network:
```
Show example here
```

Mention here that a regularization term can be added to the OLS objective function.


