/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <test_functions.h>
#include <iostream>
#include <testingutilities.h>

using namespace SPLINTER;

TestFunction::TestFunction(unsigned int n, Term *func)
        : numVariables(n),
          f(func)
{
    calculateJacobian();
    calculateHessian();
}
TestFunction::TestFunction(unsigned int n, Term &func)
        : TestFunction(n, &func)
{
}

double TestFunction::eval(DenseVector x) const
{
    std::vector<double> xvec(numVariables);
    for(size_t i = 0; i < numVariables; i++) {
        xvec.at(i) = x(i);
    }

    return f->eval(xvec);
}

DenseMatrix TestFunction::evalJacobian(DenseVector x) const
{
    DenseMatrix jacobian(1, numVariables);

    std::vector<double> xvec(numVariables);
    for(size_t i = 0; i < numVariables; i++) {
        xvec.at(i) = x(i);
    }

    for(size_t i = 0; i < numVariables; i++) {
        jacobian(0, i) = df.at(i)->eval(xvec);
    }

    return jacobian;
}

DenseMatrix TestFunction::evalHessian(DenseVector x) const
{
    DenseMatrix hessian(numVariables, numVariables);

    std::vector<double> xvec(numVariables);
    for(size_t i = 0; i < numVariables; i++) {
        xvec.at(i) = x(i);
    }

    for(size_t i = 0; i < numVariables; i++) {
        for(size_t j = 0; j < numVariables; j++) {
            hessian(i, j) = ddf.at(i).at(j)->eval(xvec);
        }
    }

    return hessian;
}

unsigned int TestFunction::getNumVariables() const
{
    return numVariables;
}

void TestFunction::calculateJacobian()
{
    df.clear();
    for(size_t i = 0; i < numVariables; i++) {
        df.push_back(f->derivative(i));
    }
}

void TestFunction::calculateHessian()
{
    ddf.clear();

    for(size_t i = 0; i < numVariables; i++) {
        ddf.push_back(std::vector<Term *>());

        for(size_t j = 0; j < numVariables; j++) {
            ddf.at(i).push_back(df.at(i)->derivative(j));
        }
    }
}


Plus::Plus(const Term *left, const Term *right)
    : left(left->clone()),
      right(right->clone())
{
}

double Plus::eval(const std::vector<double> &x) const
{
    return left->eval(x) + right->eval(x);
}

Term *Plus::derivative(int x) const
{
    auto d_left = left->derivative(x);
    auto d_right = right->derivative(x);

    auto result = new Plus(d_left, d_right);

    delete d_left;
    delete d_right;

    return result;
}

Term *Plus::clone() const
{
    return new Plus(left, right);
}

Term *Plus::simplify()
{
    auto l = const_cast<Term *>(left);
    auto r = const_cast<Term *>(right);

    left = l->simplify();

    right = r->simplify();

    auto l_const = dynamic_cast<const Const *>(left);
    auto r_const = dynamic_cast<const Const *>(right);

    if(l_const != nullptr) {
        // a + b where a and b are constants
        double l_const_val = l_const->getVal();
        if(r_const != nullptr) {
            double tot_const_val = l_const_val + r_const->getVal();
            delete left;
            delete right;
            delete this;
            return new Const(tot_const_val);
        }
        // 0 + x
        else {
            if(l_const_val == 0.0) {
                auto toReturn = const_cast<Term *>(right);
                delete left;
                delete this;
                return toReturn;
            }
        }
    }
    // x + 0
    if(r_const != nullptr) {
        double r_const_val = r_const->getVal();
        if(r_const_val == 0.0) {
            auto toReturn = const_cast<Term *>(left);
            delete right;
            delete this;
            return toReturn;
        }
    }

    return this;
}

void Plus::pretty_text(std::ostream &out) const
{
    out << "(";
    left->pretty_text(out);
    out << " + ";
    right->pretty_text(out);
    out << ")";
}

Mul::Mul(const Term *left, const Term *right)
    : left(left->clone()),
      right(right->clone())
{
}

double Mul::eval(const std::vector<double> &x) const
{
    return left->eval(x) * right->eval(x);
}

Term *Mul::derivative(int x) const
{
    auto d_left = left->derivative(x);
    auto d_right = right->derivative(x);

    // if f(x,y) = x*y then
    // df/dx = dx/dx * dy/dx + x * dy/dx
    auto d_left_times_right = Mul(d_left, right);
    auto left_times_d_right = Mul(left, d_right);

    auto result = new Plus(&d_left_times_right, &left_times_d_right);

    delete d_left;
    delete d_right;

    return result;
}

Term *Mul::clone() const
{
    return new Mul(left, right);
}

Term *Mul::simplify()
{
    auto l = const_cast<Term *>(left);
    auto r = const_cast<Term *>(right);

    left = l->simplify();

    right = r->simplify();

    auto l_const = dynamic_cast<const Const *>(left);
    auto r_const = dynamic_cast<const Const *>(right);

    if(l_const != nullptr) {
        double const_val = l_const->getVal();

        // a*b where a and b are constants
        if(r_const != nullptr) {
            double tot_const = const_val * r_const->getVal();
            delete left;
            delete right;
            delete this;

            return new Const(tot_const);
        }

        // 0 * x
        if(const_val == 0.0) {
            delete left;
            delete right;
            delete this;
            return new Const(0.0);
        }
        // 1 * x
        else if(const_val == 1.0) {
            auto toReturn = const_cast<Term *>(right);
            delete left;
            delete this; // Commenting out this line avoids segfault with f = Log(10, x);
            return toReturn;
        }
    }
    if(r_const != nullptr) {
        double const_val = r_const->getVal();
        // x * 0
        if (const_val == 0.0) {
            delete left;
            delete right;
            delete this;
            return new Const(0.0);
        }
        // x * 1
        else if(const_val == 1.0) {
            auto toReturn = const_cast<Term *>(left);
            delete right;
            delete this;
            return toReturn;
        }
    }

    return this;
}

void Mul::pretty_text(std::ostream &out) const
{
    left->pretty_text(out);
    out << "*";
    right->pretty_text(out);
}

Var::Var(int varNum)
    : varNum(varNum)
{
}

double Var::eval(const std::vector<double> &x) const
{
    return x.at(varNum);
}

Term *Var::derivative(int x) const
{
    if(x == varNum) {
        return new Const(1);
    }
    else
    {
        return new Const(0);
    }
}

Term *Var::clone() const
{
    return new Var(varNum);
}

Term *Var::simplify()
{
    return this;
}

void Var::pretty_text(std::ostream &out) const
{
    out << "x" << varNum;
}

Const::Const(double val)
        : val(val)
{
}

double Const::eval(const std::vector<double> &x) const
{
    return val;
}

Term *Const::derivative(int x) const
{
    return new Const(0);
}

Term *Const::clone() const
{
    return new Const(val);
}

Term *Const::simplify()
{
    return this;
}

void Const::pretty_text(std::ostream &out) const
{
    out << val;
}

Exp::Exp(const Term *base, const Term *exponent)
    : base(base->clone()),
      exponent(exponent->clone())
{
}

double Exp::eval(const std::vector<double> &x) const
{
    return std::pow(base->eval(x), exponent->eval(x));
}

Term *Exp::derivative(int x) const
{
    auto f = base;
    auto df = f->derivative(x);
    auto g = exponent;
    auto dg = g->derivative(x);

    // f^g
    auto fg = this;

    auto minus_one = Const(-1);

    // 1/f
    auto one_over_f = Exp(f, &minus_one);

    auto g_over_f = Mul(g, &one_over_f);

    // df * g/f
    auto df_mul_g_over_f = Mul(df, &g_over_f);

    // ln(f)
    auto lnf = Log(std::exp(1.0), f);

    // dg * ln(f)
    auto dg_ln_f = Mul(dg, &lnf);

    auto temp = Plus(&df_mul_g_over_f, &dg_ln_f);

    auto result = new Mul(fg, &temp);

    delete df;
    delete dg;

    return result;
}

Term *Exp::clone() const
{
    return new Exp(base, exponent);
}

Term *Exp::simplify()
{
    auto l = const_cast<Term *>(base);
    auto r = const_cast<Term *>(exponent);

    base = l->simplify();

    exponent = r->simplify();

    return this;
}

void Exp::pretty_text(std::ostream &out) const
{
    base->pretty_text(out);
    out << "^(";
    exponent->pretty_text(out);
    out << ")";
}

Log::Log(double base, const Term *arg)
        : base(base),
          arg(arg->clone())
{
}

Log::Log(double base, const Term &arg)
        : Log(base, &arg)
{
}

double Log::eval(const std::vector<double> &x) const
{
    return log(base, arg->eval(x));
}

Term *Log::derivative(int x) const
{
    // d(log(f(x)))/dx = df(x) / (f(x) * ln(base))
    auto d_arg = arg->derivative(x);
    auto minus_one = Const(-1);

    auto log_base_of_e = Const(log(base, std::exp(1.0)));

    auto one_over_arg = Exp(arg, &minus_one);

    auto temp = Mul(d_arg, &one_over_arg);

    auto result = new Mul(&log_base_of_e, &temp);

    delete d_arg;

    return result;
}

Term *Log::clone() const
{
    return new Log(base, arg);
}

Term *Log::simplify()
{
    auto r = const_cast<Term *>(arg);

    arg = r->simplify();

    return this;
}

void Log::pretty_text(std::ostream &out) const
{
    out << "log_" << base << "(";
    arg->pretty_text(out);
    out << ")";
}


Sin::Sin(const Term *arg)
        : arg(arg->clone())
{
}

Sin::Sin(const Term &arg)
        : Sin(&arg)
{
}

double Sin::eval(const std::vector<double> &x) const
{
    return std::sin(arg->eval(x));
}

Term *Sin::derivative(int x) const
{
    // dsin(f(x))/dx = df(x)/dx * cos(f(x))
    auto d_arg = arg->derivative(x);

    auto cos = Cos(arg);

    auto result = new Mul(d_arg, &cos);

    delete d_arg;

    return result;
}

Term *Sin::clone() const
{
    return new Sin(arg);
}

Term *Sin::simplify()
{
    auto r = const_cast<Term *>(arg);

    arg = r->simplify();

    return this;
}

void Sin::pretty_text(std::ostream &out) const
{
    out << "sin(";
    arg->pretty_text(out);
    out << ")";
}


Cos::Cos(const Term *arg)
        : arg(arg->clone())
{
}

Cos::Cos(const Term &arg)
        : Cos(&arg)
{
}

double Cos::eval(const std::vector<double> &x) const
{
    return std::cos(arg->eval(x));
}

Term *Cos::derivative(int x) const
{
    // dcos(f(x))/dx = -1 * df(x)/dx * sin(f(x))
    auto d_arg = arg->derivative(x);

    auto sin = Sin(arg);

    auto minus_one = Const(-1);

    auto temp = Mul(&minus_one, &sin);

    auto result = new Mul(d_arg, &temp);

    delete d_arg;

    return result;
}

Term *Cos::clone() const
{
    return new Cos(arg);
}

Term *Cos::simplify()
{
    auto r = const_cast<Term *>(arg);

    arg = r->simplify();

    return this;
}

void Cos::pretty_text(std::ostream &out) const
{
    out << "cos(";
    arg->pretty_text(out);
    out << ")";
}


Tan::Tan(const Term *arg)
        : arg(arg->clone())
{
}

Tan::Tan(const Term &arg)
        : Tan(&arg)
{
}

double Tan::eval(const std::vector<double> &x) const
{
    return std::tan(arg->eval(x));
}

Term *Tan::derivative(int x) const
{
    // dtan(f(x))/dx = df(x) / cos^2(x)
    auto d_arg = arg->derivative(x);

    auto cos = Cos(arg);

    auto denominator = Mul(&cos, &cos);

    auto minus_one = Const(-1);

    auto temp = Exp(&denominator, &minus_one);

    auto result = new Mul(d_arg, &temp);

    delete d_arg;

    return result;
}

Term *Tan::clone() const
{
    return new Tan(arg);
}

Term *Tan::simplify()
{
    auto r = const_cast<Term *>(arg);

    arg = r->simplify();

    return this;
}

void Tan::pretty_text(std::ostream &out) const
{
    out << "tan(";
    arg->pretty_text(out);
    out << ")";
}


E::E(const Term *arg)
        : arg(arg->clone())
{
}

E::E(const Term &arg)
        : E(&arg)
{
}

double E::eval(const std::vector<double> &x) const
{
    return std::exp(arg->eval(x));
}

Term *E::derivative(int x) const
{
    // de(f(x))/dx = df(x) * e(f(x))
    auto d_arg = arg->derivative(x);

    auto e = E(arg);

    auto result = new Mul(d_arg, &e);

    delete d_arg;

    return result;
}

Term *E::clone() const
{
    return new E(arg);
}

Term *E::simplify()
{
    auto r = const_cast<Term *>(arg);

    arg = r->simplify();

    return this;
}

void E::pretty_text(std::ostream &out) const
{
    out << "exp(";
    arg->pretty_text(out);
    out << ")";
}





Plus operator+(const Term &lhs, const Term &rhs)
{
    return Plus(&lhs, &rhs);
}

Plus operator+(const Term &lhs, double rhs)
{
    Const temp(rhs);
    return Plus(&lhs, &temp);
}

Plus operator+(double lhs, const Term &rhs)
{
    Const temp(lhs);
    return Plus(&temp, &rhs);
}


Plus operator-(const Term &lhs, const Term &rhs)
{
    Const temp(-1);
    return Plus(&lhs, new Mul(&temp, &rhs));
}

Plus operator-(const Term &lhs, double rhs)
{
    Const temp(-rhs);
    return Plus(&lhs, &temp);
}

Plus operator-(double lhs, const Term &rhs)
{
    Const temp(lhs);
    Const temp2(-1);
    return Plus(&temp, new Mul(&temp2, &rhs));
}


Mul operator*(const Term &lhs, const Term &rhs)
{
    return Mul(&lhs, &rhs);
}

Mul operator*(double lhs, const Term &rhs)
{
    Const temp(lhs);
    return Mul(&temp, &rhs);
}

Mul operator*(const Term &lhs, double rhs)
{
    Const temp(rhs);
    return Mul(&lhs, &temp);
}


Mul operator/(const Term &lhs, const Term &rhs)
{
    Const temp(-1);
    return Mul(&lhs, new Exp(&rhs, &temp));
}

Mul operator/(double lhs, const Term &rhs)
{
    Const temp(lhs);
    Const temp2(-1);
    return Mul(&temp, new Exp(&rhs, &temp2));
}

Mul operator/(const Term &lhs, double rhs)
{
    Const temp(rhs);
    Const temp2(-1);
    return Mul(&lhs, new Exp(&temp, &temp2));
}


Exp operator^(const Term &lhs, const Term &rhs)
{
    return Exp(&lhs, &rhs);
}

Exp operator^(const Term &base, double exp)
{
    Const temp(exp);
    return Exp(&base, &temp);
}

Exp operator^(double base, const Term &exp)
{
    Const temp(base);
    return Exp(&temp, &exp);
}