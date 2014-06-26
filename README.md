
Steps:
#1: git clone
#2: cd bspline
#3: cmake .
#4: make
#5: sudo make install

EIGEN_DIRECTORY (default: "/usr/local/include/eigen3")
EIGEN_DIRECTORY (default: "include/bspline")
EIGEN_DIRECTORY (default: "lib/bspline")

Troubleshooting:
If you get "fatal error: Eigen/Dense: No such file or directory", that is because the compiler could not find Eigen. You need to run step # again with -DEIGEN_DIRECTORY:STRING="/path/to/eigen"

