### TODO list

#### Important
- Add facilities for assessing model accuracy (MSE, cross-validation, etc.)
- Finish implementation of operator== overloads in operator_overloads.cpp
- Extract symbolic function functionality into its own project

#### Normal
- Implement Hessian for radial basis function splines
- Remove code redundancy by wrapping affine operations on the control points (bspline) and tensor products (bsplinebasis)
- Improve constructors of BSpline
- Implement support for scattered data interpolation with B-Splines (look into MBA algorithm)
- Implement NURBS
- Implement natural end conditions for B-spline
- Investigate integrated B-splines and tension B-splines: http://www.cas.mcmaster.ca/~modersit/Pubs/2009-SPIE-BFM.pdf
- Consider creating a class for interpolating B-splines
- Replace assertions with exceptions (assertions halts execution when testing)
- Test executable runtime speed on different compilers
- Add tests that provoke errors
- Automatically download Catch? See https://github.com/philsquared/Catch/blob/master/docs/build-systems.md#cmake
- Provide precompiled debug binaries?
