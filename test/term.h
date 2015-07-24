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
#include "termfunction.h"
#include <set>

using namespace SPLINTER;


class Var;
class Term
{
public:
    virtual double eval(const std::vector<double> &x) const = 0;

    // Returns the derivative with respect to the variable x
    virtual Term *derivative(Var x) const = 0;

    virtual Term *clone() const = 0;

    Term *simplify();

    virtual Term *concatenate() = 0;

    virtual Term *flatten() = 0;

    virtual void getVariables(std::set<Var> &vars) const = 0;

    virtual inline bool isConstant() const { return false; }

    virtual inline double getConstValue() const { return 0.0; }

    virtual inline bool isConstDegree() const { return true; }

    virtual inline double getConstDegree() const { return 1.0; }

    virtual inline void sort() {}

    virtual void pretty_text(std::ostream &out) const = 0;

    virtual ~Term();
};

class MultiTerm : public Term
{
public:
    virtual ~MultiTerm();

    //void add(Term *term);
    void getVariables(std::set<Var> &vars) const override;
    Term *concatenate() override;
    Term *flatten() override;
    bool isConstDegree() const override;
    inline const std::vector<Term *> &getTerms() const { return terms; }

protected:
    std::vector<Term *> terms;
};

class Plus : public MultiTerm
{
public:
    Plus();
    Plus(Term *lhs, Term *rhs);

    virtual ~Plus();

    void add(Term *term);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void sort() override;

    Term *concatenate() override;

    Term *flatten() override;

    double getConstDegree() const override;

    void pretty_text(std::ostream &out) const override;

private:
    void steal_children();
};


class Mul : public MultiTerm
{
public:
    Mul();
    Mul(double coefficient);
    Mul(Term *lhs, Term *rhs);

    virtual ~Mul();

    void add(Term *term);

    inline void setCoefficient(double newCoefficient) { coefficient = newCoefficient; }
    inline double getCoefficient() { return coefficient; }

    void multiply(double coefficient);

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    void sort() override;

    Term *clone() const override;

    Term *concatenate() override;

    Term *flatten() override;

    bool isConstant() const override;

    inline double getConstValue() const override { return coefficient; }

    double getConstDegree() const override;

    void pretty_text(std::ostream &out) const override;

private:
    double coefficient;

    void steal_children();
    Term *multiply_by_plus();
};


class Var : public Term
{
public:
    Var(int varNum);
    Var(int varNum, const char *name);

    Var(const Var &var);

    // Used by std::sort
    Var &operator=(Var &&other);

    virtual ~Var();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    inline Term *concatenate() { return this; }

    inline Term *flatten() { return this; }

    void getVariables(std::set<Var> &vars) const override;

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

    inline Term *concatenate() { return this; }

    inline Term *flatten() { return this; }

    inline bool isConstant() const override { return true; }

    inline double getConstValue() const override { return val; }

    inline double getConstDegree() const override { return 0.0; }

    void getVariables(std::set<Var> &vars) const override {}

    void pretty_text(std::ostream &out) const override;

    inline double getVal() const { return val; }

private:
    double val;
};

class Exp : public Term
{
public:
    Exp(Term *base, Term *exponent);

    virtual ~Exp();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    Term *concatenate() override;

    Term *flatten() override;

    void getVariables(std::set<Var> &vars) const override;

    bool isConstDegree() const { return exponent->isConstant(); }

    double getConstDegree() const { return exponent->getConstValue(); }

    void sort() override;

    void pretty_text(std::ostream &out) const override;

    inline Term *getBase() { return base; }
    inline Term *getExponent() { return exponent; }
    inline void setExponent(Term *newExponent) { exponent = newExponent; }
    inline void setBase(Term *newBase) { base = newBase; }

private:
    Term *base;
    Term *exponent;
};


class FuncTerm : public Term
{
public:
    FuncTerm(Term *arg);
    virtual ~FuncTerm();

    Term *concatenate() override;
    inline Term *flatten() override { arg->flatten(); return this; }
    void getVariables(std::set<Var> &vars) const override;
    inline void sort() override { arg->sort(); }
    inline Term *getArg() { return arg; }

protected:
    Term *arg;
};


class Log : public FuncTerm
{
public:
    Log(Term *arg);
    Log(const Term &arg);
    Log(double base, Term *arg);
    Log(double base, const Term &arg);

    virtual ~Log();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void pretty_text(std::ostream &out) const override;

    inline double getBase() const { return base; }

private:
    double base;
};


class Sin : public FuncTerm
{
public:
    Sin(Term *arg);
    Sin(const Term &arg);

    virtual ~Sin();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void pretty_text(std::ostream &out) const override;
};


class Cos : public FuncTerm
{
public:
    Cos(Term *arg);
    Cos(const Term &arg);

    virtual ~Cos();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void pretty_text(std::ostream &out) const override;
};


class Tan : public FuncTerm
{
public:
    Tan(Term *arg);
    Tan(const Term &arg);

    virtual ~Tan();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void pretty_text(std::ostream &out) const override;
};


class E : public FuncTerm
{
public:
    E(Term *arg);
    E(const Term &arg);

    virtual ~E();

    double eval(const std::vector<double> &x) const override;

    Term *derivative(Var x) const override;

    Term *clone() const override;

    void pretty_text(std::ostream &out) const override;
};

bool equal(Term *lhs, Term *rhs);

bool plusCompare(Term *lhs, Term *rhs);
bool mulCompare(Term *lhs, Term *rhs);
bool compare(Term *lhs, Term *rhs);

bool operator<(const Var &lhs, const Var &rhs);
std::ostream &operator<<(std::ostream &out, const Term &term);

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
