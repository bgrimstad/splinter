##Python interface
We have made an interface so you can use SPLINTER through Python. If you follow the instructions provided below you should be good to go in a matter of minutes (unless you choose to compile the library yourself).
The interface has been tested and should work with both Python 2.7 and 3.x.

###Setup
First, head to the [releases tab](https://github.com/bgrimstad/splinter/releases) and download a copy of splinter-python. Or, you can compile the library for yourself (as described in the [documentation](../docs/compile.md)). After doing the step `make install`, you should have a directory structure that looks like this:
- main
  - lib
    - windows
      - x86
        - splinter-x-y.dll
      - x86-64
        - splinter-x-y.dll
    - linux
      - x86
        - libsplinter-x-y.so
      - x86-64
        - libsplinter-x-y.so
    - osx
      - x86
        - libsplinter-x-y.so
      - x86-64
        - libsplinter-x-y.so
  - include
      - cinterface.h
  - python
      - splinter
          - All files from the matlab directory in the repository
    
- The numbers in the file name (x-y) corresponds to the SPLINTER version, where x is the major and y is the minor version.

Make sure the folder called python is in your path, or that that folder is your current working directory. Then you can do
`import splinter`
and it should automatically load all classes along with the binary for you.

### Basic usage
You can then start using the library by making a DataTable and populate it with samples using addSample().

```
import splinter

def f(x, y):
    return x**2*y + y**3

d = splinter.DataTable()

for i in range(10):
    for j in range(10):
        d.addSample([i,j], f(i,j))

bspline = splinter.BSpline(d, 3) # Make a BSpline of degree 3 in all dimensions

# Evaluate one point:
approxVal = bspline.eval([0.5,0.5])

# Evaluate three points:
approxVal = bspline.eval([0.3,0.2, 0.4,0.2, 7.8,7.9])
# or, equivalent:
approxVal = bspline.eval([[0.3,0.2], [0.4,0.2], [7.8,7.9]])

# Save BSpline to file for loading it later:
bspline.save("myfile.myextension")

# Load BSpline from file:
loadedBSpline = splinter.BSpline("myfile.myextension")

```
Notice that if you are going to evaluate the approximant in more than one point, it is preferred to call eval once, instead of n times. This is because you then only make a call to the binary one time, instead of n times.

When you are done populating the DataTable you can create an Approximant of your choosing:
- BSpline
- PSpline
- RadialBasisFunction
- PolynomialRegression

All these derive from the Approximant base class, and their usage mainly differ in the signature of the constructor.

Please consult the documentation for the C++ version of the library if you still have unanswered questions after reading this document.
