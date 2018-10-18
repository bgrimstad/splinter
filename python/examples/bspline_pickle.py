###################################################################################
# Pickle/unpickle functionality for the BSpline class.                            #
# Works by creating a named temporary file which the BSpline is saved to.         #
# The temp file is then read and the data is returned to Pickle.                  #
# When unpickling the data we get from Pickle is written to a temporary file      #
# which is then read by the SPLINTER backend.                                     #
#                                                                                 #
# This is a hack and tested for 5 minutes before I went to bed, so use with care. #
###################################################################################

import pickle
import copyreg
import tempfile

import splinter

# TODO: Change me
splinter.load("/home/awenhaug/projects/splinter/build/splinter-python/lib/linux/x86-64/libsplinter-3-0.so")


def constructor(serialized_data) -> splinter.BSpline:
    print(f"Unpickling BSpline")
    with tempfile.NamedTemporaryFile() as temp:
        with open(temp.name, "wb") as f:
            f.write(serialized_data)

        return splinter.BSpline(temp.name)


def reducer(bspline: splinter.BSpline):
    print(f"Pickling BSpline")
    with tempfile.NamedTemporaryFile() as temp:
        bspline.save(temp.name)
        with open(temp.name, "rb") as f:
            data = f.read()

    return constructor, (data, )


# Register reducer as the reducer for objects of type BSpline
copyreg.pickle(splinter.BSpline, reducer)


xs = list(range(10))
ys = [4.1*x**3 - 1.3*x**2 for x in xs]

filename = "bspline_pickle_example.p"

try:
    bspline = pickle.load(open(filename, "rb"))
    print(f"Loaded BSpline from {filename}")
except FileNotFoundError:
    bspline = splinter.BSplineBuilder(xs, ys).build()
    pickle.dump(bspline, open(filename, "wb"))
    print(f"Saved BSpline to {filename}")


print(bspline.eval([(x1+x0)/2 for x0, x1 in zip(xs, xs[1:])]))
