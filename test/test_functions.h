/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TEST_FUNCTIONS_H
#define SPLINTER_TEST_FUNCTIONS_H

#include "generaldefinitions.h"
#include "function.h"
#include <string.h>
#include <iostream>

using namespace SPLINTER;

class Term;

class TestFunction : public Function
{
public:
    TestFunction(unsigned int n, Term *func);
    TestFunction(unsigned int n, Term &func);

    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override;

    unsigned int getNumVariables() const override;

private:
    unsigned int numVariables;

    Term *f;
    std::vector<Term *> df;
    std::vector<std::vector<Term *>> ddf;

    void calculateJacobian();
    void calculateHessian();
};


class Plus;
class Mul;
class Exp;

class Term
{
public:
    virtual double eval(const std::vector<double> &x) const = 0;

    // Returns the derivative with respect to the variable x
    virtual Term *derivative(int x) const = 0;

    virtual Term *clone() const = 0;

    virtual Term *simplify() = 0;

    virtual void pretty_text(std::ostream &out) const = 0;


    friend Plus operator+(const Term &lhs, const Term &rhs);
    friend Plus operator+(const Term &lhs, double rhs);
    friend Plus operator+(double lhs, const Term &rhs);

    friend Plus operator-(const Term &lhs, const Term &rhs);
    friend Plus operator-(const Term &lhs, double rhs);
    friend Plus operator-(double lhs, const Term &rhs);

    friend Mul operator*(const Term &lhs, const Term &rhs);
    friend Mul operator*(const Term &lhs, double rhs);
    friend Mul operator*(double lhs, const Term &rhs);

    friend Mul operator/(const Term &lhs, const Term &rhs);
    friend Mul operator/(const Term &lhs, double rhs);
    friend Mul operator/(double lhs, const Term &rhs);

    friend Exp operator^(const Term &lhs, const Term &rhs);
    friend Exp operator^(const Term &base, double exp);
    friend Exp operator^(double base, const Term &exp);

    virtual ~Term() {};
};

class Plus : public Term
{
public:
    Plus(const Term *left, const Term *right);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *left;
    const Term *right;
};


class Mul : public Term
{
public:
    Mul(const Term *left, const Term *right);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *left;
    const Term *right;
};

class Var : public Term
{
public:
    Var(int varNum);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

    inline int getVarNum() const { return varNum; }

private:
    int varNum;
};


class Const : public Term
{
public:
    Const(double val);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

    inline double getVal() const { return val; }

private:
    double val;
};

class Exp : public Term
{
public:
    Exp(const Term *base, const Term *exponent);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *base;
    const Term *exponent;
};

class Log : public Term
{
public:
    Log(double base, const Term *arg);
    Log(double base, const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    double base;
    const Term *arg;
};

class Sin : public Term
{
public:
    Sin(const Term *arg);
    Sin(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *arg;
};

class Cos : public Term
{
public:
    Cos(const Term *arg);
    Cos(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *arg;
};

class Tan : public Term
{
public:
    Tan(const Term *arg);
    Tan(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *arg;
};

class E : public Term
{
public:
    E(const Term *arg);
    E(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(int x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    const Term *arg;
};


#endif //SPLINTER_TEST_FUNCTIONS_H