### TODO list

#### Important
- Add facilities for assessing model accuracy (MSE, cross-validation, etc.)
- Finish implementation of operator== overloads in operator_overloads.cpp
- Add Python interface (2 and 3?)  (Investigate if a single C interface may be used by both Python and MatLab)

#### Normal
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
- Default compilation target architecture should be equal to the architecture of the host computer
