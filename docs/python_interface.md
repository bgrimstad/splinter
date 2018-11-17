## Python interface
NOTE: If you intend to use SPLINTER with Python on Linux x64, you can install SPLINTER through PyPI: `pip install splinterpy`

We have made an interface so you can use SPLINTER through Python. If you follow the instructions provided below you should be good to go in a matter of minutes (unless you choose to compile the library yourself).
The interface has been tested and should work with both Python 2.7 and 3.x.

### Setup
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
          - All files from the python directory in the repository
    
- The numbers in the filename (x-y) corresponds to the SPLINTER version, where x is the major and y is the minor version.

Make sure the folder called python is in your path, or that that folder is your current working directory. Then you can do
`import splinter`
and it should automatically load all classes along with the binary for you.

### Basic usage

```python
import splinter

# Example with one variable
def f1(x):
    return -1. + 2*x + 0.1*(x**2) + 10*np.random.rand(1)[0]

x = np.arange(0, 11, 1)
y = np.array([f1(xi) for xi in x])

# Build cubic B-spline that interpolates the data
bspline = splinter.bspline_interpolator(x, y, degree=3)

# Evaluate the B-spline
xd = np.arange(0, 10, .01)
yd = bspline.eval(xd)

# Save BSpline to file for loading it later:
filename = "bspline.json"
bspline.to_json(filename)

# Load BSpline from file:
loaded_bspline = splinter.BSpline.from_json(filename)
```
Notice that if you are going to evaluate the BSpline in more than one point, it is preferred to call eval once, instead of n times.

More examples can be found in the directory `/python/examples`.

Please consult the documentation for the C++ version of the library if you still have unanswered questions after reading this document.
