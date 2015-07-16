/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*
 * Regarding the suiciding the objects in this file does:
 * https://isocpp.org/wiki/faq/freestore-mgmt#delete-this
 */

#include <test_functions.h>
#include <iostream>
#include <testingutilities.h>

bool operator<(const Var &lhs, const Var &rhs) {
    return lhs.getVarNum() < rhs.getVarNum();
}

std::ostream &operator<<(std::ostream &out, const Term &term)
{
    term.pretty_text(out);
    return out;
}

bool Term::isConstant() const
{
    auto temp = const_cast<Term *>(this);
    return dynamic_cast<Const *>(temp) != nullptr;
}

Term::~Term()
{
}

Plus::Plus()
{
}

Plus::Plus(const Term &lhs, const Term &rhs)
    : Plus(lhs.clone(), rhs.clone())
{
}

Plus::Plus(const Term &lhs, Term *rhs)
    : Plus(lhs.clone(), rhs)
{
}

Plus::Plus(Term *lhs, const Term &rhs)
    : Plus(lhs, rhs.clone())
{
}

Plus::Plus(Term *lhs, Term *rhs)
{
    add(lhs);
    add(rhs);
}

Plus::~Plus()
{
//    for(auto &term : terms) {
//        delete term;
//    }
}

void Plus::add(const Term &term)
{
    terms.push_back(term.clone());
}

void Plus::add(Term *term)
{
    terms.push_back(term);
}

/* TODO: Maybe better to evaluate in a tree-like fashion so
 * we avoid round off errors? F.ex. when we have 4 terms:
 * temp0 = term0 + term1
 * temp1 = term2 + term3
 * return temp0 + temp1
 */
double Plus::eval(const std::vector<double> &x) const
{
    auto tot = 0.0;
    for(auto &term : terms) {
        tot += term->eval(x);
    }
    return tot;
}

Term *Plus::derivative(Var x) const
{
    auto result = new Plus();

    for(auto &term : terms) {
        result->add(term->derivative(x));
    }

    return result;
}

Term *Plus::clone() const
{
    auto result = new Plus();
    result->terms.reserve(terms.size());
    for(auto &term : terms) {
        result->terms.push_back(term->clone());
    }
    return result;
}

// Find all terms that are Plus *, and "steal" their children.
// This can be done because x0 + (x1 + x2) = x0 + x1 + x2
void Plus::steal_children() {
    // Store the current size. It will increase if we need to steal terms,
    // and there's no need to see if those need stealing, as that has been
    // done in the simplify loop above.
    size_t size = terms.size();
    for(int i = size - 1; i >= 0; --i) {
        // See if the term is a Plus, and if so, we "steal" all its children
        auto plus = dynamic_cast<Plus *>(terms.at(i));
        if(plus != nullptr) {
            // Resize now for efficiency
            terms.erase(terms.begin() + i); // Erase the Plus entry
            terms.reserve(terms.size() + plus->terms.size());
            for(int j = plus->terms.size() - 1; j >= 0; --j) {
                terms.push_back(plus->terms.at(j));
                plus->terms.erase(plus->terms.begin() + j);
            }
            delete plus;
        }
    }
}

Term *Plus::simplify()
{
    for(auto &term : terms) {
        term = term->simplify();
    }

    steal_children();

    size_t size = terms.size();
    double totConstVal = 0.0;
    for(int i = size - 1; i >= 0; --i) {
        auto constant = dynamic_cast<Const *>(terms.at(i));
        if(constant != nullptr) {
            totConstVal += constant->getVal();

            auto term = terms.at(i);
            terms.erase(terms.begin() + i);
            delete term;
        }
    }

    if(totConstVal != 0.0) {
        add(new Const(totConstVal));
    }

    // TODO: Sort terms, then concatenate as possible

    // If this is empty it means all terms summed to 0.0
    if(terms.size() == 0) {
        delete this;
        return new Const(0.0);

    } else if(terms.size() == 1) {
        auto temp = terms.at(0)->clone();
        delete this;
        return temp;
    }

    return this;
}

std::set<Var> Plus::getVariables() const {
    auto vars = std::set<Var>();

    for(auto &term : terms) {
        auto termVars = term->getVariables();
        for(auto &termVar : termVars) {
            vars.insert(termVar);
        }
    }

    return vars;
}

void Plus::pretty_text(std::ostream &out) const
{
    out << "(";
    for(auto it = terms.cbegin(); it != terms.cend(); ++it) {
        (*it)->pretty_text(out);
        if(it + 1 != terms.cend()) {
            out << " + ";
        }
    }
    out << ")";
}


Mul::Mul()
{
}

Mul::Mul(const Term &lhs, const Term &rhs)
    : Mul(lhs.clone(), rhs.clone())
{
}

Mul::Mul(const Term &lhs, Term *rhs)
    : Mul(lhs.clone(), rhs)
{
}

Mul::Mul(Term *lhs, const Term &rhs)
    : Mul(lhs, rhs.clone())
{
}

Mul::Mul(Term *lhs, Term *rhs)
{
    add(lhs);
    add(rhs);
}

Mul::~Mul()
{
//    for(auto &term : terms) {
//        delete term;
//    }
}

void Mul::add(const Term &term)
{
    terms.push_back(term.clone());
}

void Mul::add(Term *term)
{
    terms.push_back(term);
}

/* TODO: Maybe better to evaluate in a tree-like fashion so
 * we avoid round off errors? F.ex. when we have 4 terms:
 * temp0 = term0 * term1
 * temp1 = term2 * term3
 * return temp0 * temp1
 */
double Mul::eval(const std::vector<double> &x) const
{
    auto tot = 1.0;
    for(auto &term : terms) {
        tot *= term->eval(x);
    }
    return tot;
}

Term *Mul::derivative(Var x) const
{
    Plus *result = new Plus();
    for(int i = 0; i < terms.size(); ++i) {

        Mul *temp = new Mul();
        for(int j = 0; j < terms.size(); ++j) {
            if(i == j) {
                temp->add(terms.at(j)->derivative(x));
            } else {
                temp->add(terms.at(j)->clone());
            }
        }

        result->add(temp);
    }

    return result;
}

Term *Mul::clone() const
{
    auto result = new Mul();
    result->terms.reserve(terms.size());
    for(auto &term : terms) {
        result->terms.push_back(term->clone());
    }
    return result;
}

// Find all terms that are Mul *, and "steal" their children.
// This can be done because x0 * (x1 * x2) = x0 * x1 * x2
void Mul::steal_children() {
    // Store the current size. It will increase if we need to steal terms,
    // and there's no need to see if those need stealing, as that has been
    // done in the simplify loop above.
    size_t size = terms.size();
    for(int i = size - 1; i >= 0; --i) {
        // See if the term is a Mul, and if so, we "steal" all its children
        auto mul = dynamic_cast<Mul *>(terms.at(i));
        if(mul != nullptr) {
            // Resize now for efficiency
            terms.erase(terms.begin() + i); // Erase the Mul entry
            terms.reserve(terms.size() + mul->terms.size());
            for(int j = mul->terms.size() - 1; j >= 0; --j) {
                terms.push_back(mul->terms.at(j));
                mul->terms.erase(mul->terms.begin() + j);
            }
            delete mul;
        }
    }
}

/* See if there is a plus that is a child term, and if so,
 * multiply by it so we get a resulting Plus with the same number of terms
 * as the previous Plus
 */
Term *Mul::multiply_by_plus() {
    for(auto it = terms.end() - 1; it >= terms.begin(); --it) {
        auto term = *it;
        auto plus = dynamic_cast<Plus *>(term);

        if(plus != nullptr) {
            auto result = new Plus();

            auto plus_terms = plus->getTerms();

            // Important to erase this _before_ the call to
            // this->clone below, else the plus will get
            // cloned too and start an infinite loop
            terms.erase(it);

            for(auto &plus_term : plus_terms) {
                auto temp = new Mul();
                temp->add(plus_term->clone());
                temp->add(this->clone());
                result->add(temp);
            }

            //delete this;
            return result->simplify();
        }
    }

    return this;
}

Term *Mul::simplify()
{
    /* TODO:
     * look for Plus terms to multiply by,
     * concatenate identical variables
     */
    for(auto &term : terms) {
        term = term->simplify();
    }

    steal_children();

    /* TODO:
     * Sort and concatenate */

    /* Multiply all the constants and see what we end up with */
    size_t size = terms.size();
    double totConstVal = 1.0;
    for(int i = size - 1; i >= 0; --i) {
        auto constant = dynamic_cast<Const *>(terms.at(i));
        if(constant != nullptr) {
            totConstVal *= constant->getVal();

            auto term = terms.at(i);
            terms.erase(terms.begin() + i);
            delete term;
        }
    }

    if(totConstVal == 0.0) {
        delete this;
        return new Const(0.0);
    }
    if(totConstVal != 1.0) {
        terms.insert(terms.begin(), new Const(totConstVal));
    }


    for(int i = 0; i < terms.size(); ++i) {
        Var *var0 = dynamic_cast<Var *>(terms.at(i));
        Term *exp0 = nullptr;

        if(var0 == nullptr) {
            auto exp = dynamic_cast<Exp *>(terms.at(i));
            if(exp != nullptr) {
                var0 = dynamic_cast<Var *>(exp->getBase());
                exp0 = exp->getExponent();
            }
        }

        if(var0 != nullptr) {
            auto replacementExponent = new Plus();
            if(exp0 != nullptr) {
                replacementExponent->add(exp0->clone());
            } else {
                replacementExponent->add(new Const(1.0));
            }

            for(int j = i+1; j < terms.size(); ++j) {
                Var *var1 = dynamic_cast<Var *>(terms.at(j));
                Term *exp1 = nullptr;

                if(var1 == nullptr) {
                    auto exp = dynamic_cast<Exp *>(terms.at(j));
                    if(exp != nullptr) {
                        var1 = dynamic_cast<Var *>(exp->getBase());
                        exp1 = exp->getExponent();
                    }
                }

                if(var1 != nullptr) {
                    if(var0->getVarNum() == var1->getVarNum()) {
                        if(exp1 != nullptr) {
                            replacementExponent->add(exp1->clone());
                        } else {
                            replacementExponent->add(new Const(1.0));
                        }
                        delete terms.at(j); // free memory
                        terms.erase(terms.begin() + j);
                    }
                }
            }

            auto replacement = new Exp(var0->clone(), replacementExponent); // MUST be done before the next call
            //delete terms.at(i); // free memory (This causes segfault!)

            terms.at(i) = replacement->simplify();
        }
    }

    // If this is empty it means all terms was a constant with value 1.0
    if(terms.size() == 0) {
        delete this;
        return new Const(1.0);

    } else if(terms.size() == 1) {
        auto temp = terms.at(0)->clone();
        delete this;
        return temp;
    }

    auto temp = multiply_by_plus();
    if(temp != this) {
        // Don't delete this, it has been deleted already by multiply_by_plus()
        return temp;
    }

    return this;
}

std::set<Var> Mul::getVariables() const {
    auto vars = std::set<Var>();

    for(auto &term : terms) {
        auto termVars = term->getVariables();
        for(auto &termVar : termVars) {
            vars.insert(termVar);
        }
    }

    return vars;
}

void Mul::pretty_text(std::ostream &out) const
{
    for(auto it = terms.cbegin(); it != terms.cend(); ++it) {
        (*it)->pretty_text(out);
        if(it + 1 != terms.cend()) {
            out << "*";
        }
    }
}

Var::Var(int varNum)
    : Var(varNum, nullptr)
{
}

Var::Var(int varNum, const char *name)
    : varNum(varNum),
      name(name)
{
}

double Var::eval(const std::vector<double> &x) const
{
    return x.at(varNum);
}

Term *Var::derivative(Var x) const
{
    if(x.getVarNum() == varNum) {
        return new Const(1);
    }
    else
    {
        return new Const(0);
    }
}

Term *Var::clone() const
{
    return new Var(varNum, name); // TODO: Copy name?
}

Term *Var::simplify()
{
    return this;
}

std::set<Var> Var::getVariables() const {
    auto var = std::set<Var>();
    var.insert(*this);
    return var;
}

void Var::pretty_text(std::ostream &out) const
{
    if(name != nullptr) {
        out << name;
    } else {
        out << "x" << varNum;
    }
}

Const::Const(double val)
    : val(val)
{
}

double Const::eval(const std::vector<double> &x) const
{
    return val;
}

Term *Const::derivative(Var x) const
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

std::set<Var> Const::getVariables() const {
    return std::set<Var>();
}

void Const::pretty_text(std::ostream &out) const
{
    out << val;
}

Exp::Exp(const Term &base, const Term &exponent)
    : Exp(base.clone(), exponent.clone())
{
}

Exp::Exp(const Term &base, Term *exponent)
    : Exp(base.clone(), exponent)
{
}

Exp::Exp(Term *base, const Term &exponent)
    : Exp(base, exponent.clone())
{
}

Exp::Exp(Term *base, Term *exponent)
    : base(base),
      exponent(exponent)
{
}

double Exp::eval(const std::vector<double> &x) const
{
    auto baseVal = base->eval(x);
    auto expVal = exponent->eval(x);

    return std::pow(baseVal, expVal);
}

Term *Exp::derivative(Var x) const
{
    auto f = base;
    auto g = exponent;

    // f^g
    auto fg = this;

    // 1/f
    auto one_over_f = new Exp(f->clone(), new Const(-1));

    auto g_over_f = new Mul(g->clone(), one_over_f);

    // df * g/f
    auto df_mul_g_over_f = new Mul(f->derivative(x), g_over_f);

    // ln(f)
    auto lnf = new Log(std::exp(1.0), f->clone());

    // dg * ln(f)
    auto dg_ln_f = new Mul(g->derivative(x), lnf);

    auto temp = new Plus(df_mul_g_over_f, dg_ln_f);

    return new Mul(fg->clone(), temp);
}

Term *Exp::clone() const
{
    return new Exp(base->clone(), exponent->clone());
}

Term *Exp::simplify()
{
    base = base->simplify();

    exponent = exponent->simplify();

    auto constantBase = dynamic_cast<Const *>(base);
    auto constantExponent = dynamic_cast<Const *>(exponent);
    if(constantExponent != nullptr) {
        // Both are constants, eval and return constant
        if(constantBase != nullptr) {
            std::vector<double> dummy;
            auto constVal = eval(dummy);
            delete this;
            return new Const(constVal);
        }
        auto exponentVal = constantExponent->getVal();
        if(exponentVal == 0.0) {
            delete this;
            return new Const(1.0);
        }
        // Don't do anything is exponentVal is 1, it will get
        // converted from x^1 ->x -> x^1 forever
    }
    // Handle case where base = 0, so:
    // 0^x
    else if(constantBase != nullptr) {
        if(constantBase->getVal() == 0.0) {
            delete this;
            return new Const(0.0);
        }
    }

    return this;
}

std::set<Var> Exp::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = base->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }
    termVars = exponent->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void Exp::pretty_text(std::ostream &out) const
{
    if(exponent->isConstant() && dynamic_cast<Const *>(exponent)->getVal() == 1.0) {
        base->pretty_text(out);
    } else {
        out << "(";
        base->pretty_text(out);
        out << "^";
        exponent->pretty_text(out);
        out << ")";
    }
}

Log::Log(double base, Term *arg)
        : base(base),
          arg(arg)
{
}

Log::Log(double base, const Term &arg)
        : Log(base, arg.clone())
{
}

double Log::eval(const std::vector<double> &x) const
{
    return log(base, arg->eval(x));
}

Term *Log::derivative(Var x) const
{
    // d(log(f(x)))/dx = df(x) / (f(x) * ln(base))
    auto log_base_of_e = new Const(log(base, std::exp(1.0)));

    auto one_over_arg = new Exp(arg->clone(), new Const(-1));

    auto darg_over_arg = new Mul(arg->derivative(x), one_over_arg);

    return new Mul(log_base_of_e, darg_over_arg);
}

Term *Log::clone() const
{
    return new Log(base, arg->clone());
}

Term *Log::simplify()
{
    arg = arg->simplify();

    if(arg->isConstant()) {
        std::vector<double> dummy;
        auto constantVal = eval(dummy);
        delete this;
        return new Const(constantVal);
    }

    return this;
}

std::set<Var> Log::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = arg->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void Log::pretty_text(std::ostream &out) const
{
    out << "log_" << base << "(";
    arg->pretty_text(out);
    out << ")";
}


Sin::Sin(Term *arg)
        : arg(arg)
{
}

Sin::Sin(const Term &arg)
        : Sin(arg.clone())
{
}

double Sin::eval(const std::vector<double> &x) const
{
    return std::sin(arg->eval(x));
}

Term *Sin::derivative(Var x) const
{
    // dsin(f(x))/dx = df(x)/dx * cos(f(x))
    return new Mul(arg->derivative(x), new Cos(arg->clone()));
}

Term *Sin::clone() const
{
    return new Sin(arg->clone());
}

Term *Sin::simplify()
{
    arg = arg->simplify();

    return this;
}

std::set<Var> Sin::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = arg->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void Sin::pretty_text(std::ostream &out) const
{
    out << "sin(";
    arg->pretty_text(out);
    out << ")";
}


Cos::Cos(Term *arg)
        : arg(arg)
{
}

Cos::Cos(const Term &arg)
        : Cos(arg.clone())
{
}

double Cos::eval(const std::vector<double> &x) const
{
    return std::cos(arg->eval(x));
}

Term *Cos::derivative(Var x) const
{
    // dcos(f(x))/dx = -1 * df(x)/dx * sin(f(x))
    auto temp = new Mul(new Const(-1), new Sin(arg->clone()));

    return new Mul(arg->derivative(x), temp);
}

Term *Cos::clone() const
{
    return new Cos(arg->clone());
}

Term *Cos::simplify()
{
    arg = arg->simplify();

    return this;
}

std::set<Var> Cos::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = arg->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void Cos::pretty_text(std::ostream &out) const
{
    out << "cos(";
    arg->pretty_text(out);
    out << ")";
}


Tan::Tan(Term *arg)
        : arg(arg)
{
}

Tan::Tan(const Term &arg)
        : Tan(arg.clone())
{
}

double Tan::eval(const std::vector<double> &x) const
{
    return std::tan(arg->eval(x));
}

Term *Tan::derivative(Var x) const
{
    // dtan(f(x))/dx = df(x) / cos^2(x)
    auto denominator = new Mul(new Cos(arg->clone()), new Cos(arg->clone()));

    auto temp = new Exp(denominator, new Const(-1));

    return new Mul(arg->derivative(x), temp);
}

Term *Tan::clone() const
{
    return new Tan(arg->clone());
}

Term *Tan::simplify()
{
    arg = arg->simplify();

    return this;
}

std::set<Var> Tan::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = arg->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void Tan::pretty_text(std::ostream &out) const
{
    out << "tan(";
    arg->pretty_text(out);
    out << ")";
}


E::E(Term *arg)
        : arg(arg)
{
}

E::E(const Term &arg)
        : E(arg.clone())
{
}

double E::eval(const std::vector<double> &x) const
{
    return std::exp(arg->eval(x));
}

Term *E::derivative(Var x) const
{
    // de(f(x))/dx = df(x) * e(f(x))
    return new Mul(arg->derivative(x), new E(arg->clone()));
}

Term *E::clone() const
{
    return new E(arg->clone());
}

Term *E::simplify()
{
    arg = arg->simplify();

    return this;
}

std::set<Var> E::getVariables() const {
    auto vars = std::set<Var>();

    auto termVars = arg->getVariables();
    for(auto &termVar : termVars) {
        vars.insert(termVar);
    }

    return vars;
}

void E::pretty_text(std::ostream &out) const
{
    out << "exp(";
    arg->pretty_text(out);
    out << ")";
}





Plus operator+(const Term &lhs, const Term &rhs)
{
    return Plus(lhs, rhs);
}

Plus operator+(const Term &lhs, double rhs)
{
    Const temp(rhs);
    return Plus(lhs, temp);
}

Plus operator+(double lhs, const Term &rhs)
{
    Const temp(lhs);
    return Plus(temp, rhs);
}


Plus operator-(const Term &lhs, const Term &rhs)
{
    Const temp(-1);
    return Plus(lhs, new Mul(temp, rhs));
}

Plus operator-(const Term &lhs, double rhs)
{
    Const temp(-rhs);
    return Plus(lhs, temp);
}

Plus operator-(double lhs, const Term &rhs)
{
    Const temp(lhs);
    Const temp2(-1);
    return Plus(temp, new Mul(temp2, rhs));
}


Mul operator*(const Term &lhs, const Term &rhs)
{
    return Mul(lhs, rhs);
}

Mul operator*(double lhs, const Term &rhs)
{
    Const temp(lhs);
    return Mul(temp, rhs);
}

Mul operator*(const Term &lhs, double rhs)
{
    Const temp(rhs);
    return Mul(lhs, temp);
}


Mul operator/(const Term &lhs, const Term &rhs)
{
    Const temp(-1);
    return Mul(lhs, new Exp(rhs, temp));
}

Mul operator/(double lhs, const Term &rhs)
{
    Const temp(lhs);
    Const temp2(-1);
    return Mul(temp, new Exp(rhs, temp2));
}

Mul operator/(const Term &lhs, double rhs)
{
    Const temp(rhs);
    Const temp2(-1);
    return Mul(lhs, new Exp(temp, temp2));
}


Exp operator^(const Term &lhs, const Term &rhs)
{
    return Exp(lhs, rhs);
}

Exp operator^(const Term &base, double exp)
{
    Const temp(exp);
    return Exp(base, temp);
}

Exp operator^(double base, const Term &exp)
{
    Const temp(base);
    return Exp(temp, exp);
}
