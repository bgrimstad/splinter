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
#include <testfunction.h>

using namespace SPLINTER;

extern std::vector<TestFunction *> testFunctions;

class Var;
class Term
{
public:
    virtual double eval(const std::vector<double> &x) const = 0;

    // Returns the derivative with respect to the variable x
    virtual Term *derivative(Var x) const = 0;

    virtual Term *clone() const = 0;

    virtual Term *simplify() = 0;

    bool isConstant() const;

    virtual void pretty_text(std::ostream &out) const = 0;

    virtual ~Term();
};

std::ostream &operator<<(std::ostream &out, const Term &term);

class Plus : public Term
{
public:
    Plus();
    Plus(const Term &lhs, const Term &rhs);
    Plus(const Term &lhs, Term *rhs);
    Plus(Term *lhs, const Term &rhs);
    Plus(Term *lhs, Term *rhs);

    virtual ~Plus();

    void add(const Term &term);
    void add(Term *term);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

    inline const std::vector<Term *> &getTerms() const { return terms; }

private:
    std::vector<Term *> terms;

    void steal_children();
};


class Mul : public Term
{
public:
    Mul();
    Mul(const Term &lhs, const Term &rhs);
    Mul(const Term &lhs, Term *rhs);
    Mul(Term *lhs, const Term &rhs);
    Mul(Term *lhs, Term *rhs);

    virtual ~Mul();

    void add(const Term &term);
    void add(Term *term);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    std::vector<Term *> terms;

    void steal_children();
    Term *multiply_by_plus();
};

class Var : public Term
{
public:
    Var(int varNum);
    Var(int varNum, const char *name);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

    inline int getVarNum() const { return varNum; }

private:
    int varNum;
    const char *name;
};


class Const : public Term
{
public:
    Const(double val);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

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
    Exp(const Term &base, const Term &exponent);
    Exp(const Term &base, Term *exponent);
    Exp(Term *base, const Term &exponent);
    Exp(Term *base, Term *exponent);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

    inline Term *getBase() { return base; }
    inline Term *getExponent() { return exponent; }

private:
    Term *base;
    Term *exponent;
};

class Log : public Term
{
public:
    Log(double base, Term *arg);
    Log(double base, const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    double base;
    Term *arg;
};

class Sin : public Term
{
public:
    Sin(Term *arg);
    Sin(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    Term *arg;
};

class Cos : public Term
{
public:
    Cos(Term *arg);
    Cos(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    Term *arg;
};

class Tan : public Term
{
public:
    Tan(Term *arg);
    Tan(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    Term *arg;
};

class E : public Term
{
public:
    E(Term *arg);
    E(const Term &arg);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *simplify() override;

    void pretty_text(std::ostream &out) const override;

private:
    Term *arg;
};

Plus operator+(const Term &lhs, const Term &rhs);
Plus operator+(const Term &lhs, double rhs);
Plus operator+(double lhs, const Term &rhs);

Plus operator-(const Term &lhs, const Term &rhs);
Plus operator-(const Term &lhs, double rhs);
Plus operator-(double lhs, const Term &rhs);

Mul operator*(const Term &lhs, const Term &rhs);
Mul operator*(const Term &lhs, double rhs);
Mul operator*(double lhs, const Term &rhs);

Mul operator/(const Term &lhs, const Term &rhs);
Mul operator/(const Term &lhs, double rhs);
Mul operator/(double lhs, const Term &rhs);

Exp operator^(const Term &lhs, const Term &rhs);
Exp operator^(const Term &base, double exp);
Exp operator^(double base, const Term &exp);

#endif //SPLINTER_TEST_FUNCTIONS_H
