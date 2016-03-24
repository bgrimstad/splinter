##MatLab interface
We have made an interface so you can use SPLINTER through MatLab. If you follow the instructions provided below you should be good to go in a matter of minutes (unless you choose to compile the library yourself).

###Setup
First, head to the [releases tab](https://github.com/bgrimstad/splinter/releases) and download a copy of splinter-matlab. Or, you can compile the library for yourself (as described in the [documentation](../docs/compile.md)). After doing the step `make install`, you should have a directory structure that looks like this:
- main
  - lib
    - windows
      - x86
        - splinter-x-y.dll
      - x86-64
        - splinter-x-y.dll
    - linux
      - x86-64
        - libsplinter-x-y.so
    - osx
      - x86-64
        - libsplinter-x-y.so
  - include
      - cinterface.h
  - matlab
    - All files from the matlab directory in the repository
    
- Note that you only need the folders that correspond to your platform.
- 32 bit MatLab only exist on Windows.
- The numbers in the file name (x-y) corresponds to the SPLINTER version, where x is the major and y is the minor version.

Run setup() (defined in setup.m), which should load SPLINTER for you automatically. You must run this once for each time you start MatLab.

### Basic usage
You can then start using the library by making a DataTable and populate it with samples using add_sample(). To speed up the process of adding samples you should call preallocate() before adding samples, if you have any idea how many samples you are going to add. Doing so will avoid having to resize the matrix keeping track of the samples.

When you are done populating the DataTable you can create an Approximant of your choosing:
- BSpline
- PSpline

All these derive from the Approximant base class, and their usage mainly differ in the signature of the constructor.

Please consult the documentation for the C++ version of the library if you still have unanswered questions after reading this document.
