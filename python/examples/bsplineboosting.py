# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
import matplotlib.pyplot as plt
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinter

# Only for dev purposes
splinter.load("/home/bjarne/Code/C++/splinter4/splinter/bin/Release/libsplinter-3-0.so")
# splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-3-0.so")


class BSplineBoosting:
    def __init__(self, loss='ls', learning_rate=0.1, n_estimators=100, subsample=1.0, alpha=1.0):
        """
        Class for stochastic gradient boosting with B-spline learners
        :param loss: loss function, 'ls' for least squares loss function
        :param learning_rate: Learning rate applied to new estimators. Advised to be in the range (0, 0.1).
        :param n_estimators: Number of estimators. Inversely proportional to learning_rate, e.g. a high number of estimators is required if the learning rate is low.
        :param subsample: Amount of samples (drawn randomly) to use when fitting a new estimator. Must be a number in (0, 1].
        :param alpha: Regularization parameter for P-spline learners. Must be a strictly positive number.
        TODO: add init and warm_start as in sklearn.ensemble.GradientBoostingRegressor
        """
        # super(Function, self).__init__()
        self._loss = ls  # Function pointer
        self._learning_rate = learning_rate
        self._n_estimators = n_estimators
        self._subsample = subsample
        self._alpha = alpha
        self._learner = 'pspline'
        self._estimators = [None] * self._n_estimators
        self._oob_improvement = np.array((self._n_estimators,))
        self._train_score = np.array((self._n_estimators,))

    def fit(self, x: np.ndarray, y: np.ndarray):
        """
        Fit to data (x, y)
        """
        n = x.shape[0]

        # Initialize first base learner
        self._estimators[0] = Const(np.mean(y) / self._learning_rate)

        for i in range(1, self._n_estimators):
            u_hat = y - self.eval(x)

            learners = [None] * n
            goodness = np.zeros((n,))
            ss_tot = np.sum(np.apply_along_axis(np.square, 0, u_hat - np.mean(u_hat)))

            for j in range(n):
                learners[j] = splinter.BSplineBuilder(x, u_hat, smoothing=splinter.BSplineBuilder.Smoothing.PSPLINE, alpha=self._alpha, knot_spacing=splinter.BSplineBuilder.KnotSpacing.EQUIDISTANT, num_basis_functions=20).build()
                u_hat_est = learners[j].eval(x)
                ss_res = np.sum(np.apply_along_axis(np.square, 0, u_hat - u_hat_est))
                goodness[j] = 1 - (ss_res / ss_tot)

            best_learner = 0
            for j in range(1, n):
                if goodness[j] > goodness[best_learner]:
                    best_learner = j

            self._estimators[i] = learners[best_learner]

    def eval(self, x):
        """
        Evaluate ensemble at x
        """
        y = np.zeros(x.shape)
        for i in range(self._n_estimators):
            if self._estimators[i] is not None:
                y += self._learning_rate * np.array(self._estimators[i].eval(x))
            else:
                break
        return y

    def predict(self, x):
        """
        Predict value at x
        """
        return self.eval(x)


class Const:
    def __init__(self, constant):
        self._const = constant

    def eval(self, x):
        return self._const * np.ones(x.shape)


def ls(x):
    return sse(x, np.zeros(x.shape))


def sse(x: np.ndarray, y: np.ndarray):
    assert(x.shape == y.shape)
    assert(x.ndim == 1)
    return np.sum(np.apply_along_axis(np.square, 0, x - y))


def mse(x: np.ndarray, y: np.ndarray):
    assert(x.shape == y.shape)
    assert(x.ndim == 1)
    return sse(x, y) / x.shape[0]


def underlying_func(x: np.array):
    return 0.5*x + 0.1*np.apply_along_axis(np.square, 0, x)


def noisy_func(x: np.array):
    return underlying_func(x) + np.random.randn(x.shape[0])


if __name__ == '__main__':

    print("Testing stochastic gradient boosting with B-spline learners")
    fh = mse
    a = np.array([3, 4, 5])
    b = np.array([2, 3, 3])
    print(fh(a, b))
    print(ls(b))

    # Sampling
    x = np.arange(0, 10, 0.1)
    xd = np.arange(0, 10, 0.01)
    y = noisy_func(x)
    yd = underlying_func(xd)

    # Just one P-spline
    pspline = splinter.BSplineBuilder(x, y, smoothing=splinter.BSplineBuilder.Smoothing.PSPLINE, alpha=0.1,
                                      knot_spacing=splinter.BSplineBuilder.KnotSpacing.EQUIDISTANT,
                                      num_basis_functions=20).build()

    yd_pspline = pspline.eval(xd)

    # Boosting
    bb = BSplineBoosting(learning_rate=0.1, n_estimators=100, alpha=1.0)
    bb.fit(x, y)

    # Prediction
    yd_boost = bb.eval(xd)

    # Plotting
    plt.plot(x, y, '*', label='Data points')
    plt.plot(xd, yd, label='Unknown function')
    plt.plot(xd, yd_pspline, label='P-spline')
    plt.plot(xd, yd_boost, label='Boosted spline')
    plt.legend(loc='upper left')
    plt.show()
