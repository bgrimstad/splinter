/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "term.h"
#include <iostream>
#include <testingutilities.h>


Term::~Term()
{
}

Term *Term::simplify()
{
    auto flattened = flatten();
//    std::cout << "Flattened: ";
//    flattened->pretty_text(std::cout);
//    std::cout << std::endl;

    auto concatenated = flattened->concatenate();
//    std::cout << "Concatenated: ";
//    concatenated->pretty_text(std::cout);
//    std::cout << std::endl;

    if(concatenated != flattened && flattened != this) {
        delete flattened;
    }

    return concatenated;
}

MultiTerm::~MultiTerm()
{
    for(auto &term : terms) {
        delete term;
    }
}

//void MultiTerm::add(Term *term)
//{
//    // Make sure all Var's are represented as Exp(Var, 1)
//    auto var = dynamic_cast<Var *>(term);
//    if(var != nullptr) {
//        auto temp = new Mul();
//        temp->add(new Exp(var, new Mul(1.0)));
//        terms.push_back(temp);
//    }
//    else
//    {
//        auto constant = dynamic_cast<Const *>(term);
//        if(constant != nullptr) {
//            terms.push_back(new Mul(constant->getConstValue()));
//        }
//        else
//        {
//            auto exp = dynamic_cast<Exp *>(term);
//            if(exp != nullptr) {
//                auto temp = new Mul();
//                temp->terms.push_back(exp);
//                terms.push_back(temp);
//            }
//            else {
//                terms.push_back(term);
//            }
//        }
//    }
//}

void MultiTerm::getVariables(std::set<Var> &vars) const
{
    for(auto &term : terms) {
        term->getVariables(vars);
    }
}

Term *MultiTerm::concatenate()
{
    for(auto &term : terms) {
        auto temp = term->concatenate();
        if(temp != term) {
            delete term;
        }
        term = temp;
    }

    return this;
}

Term *MultiTerm::flatten()
{
    for(auto &term : terms) {
        auto temp = term->flatten();
        if(temp != term) {
            delete term;
        }
        term = temp;
    }

    return this;
}

bool MultiTerm::isConstDegree() const
{
    for(auto &term : terms) {
        if(!term->isConstDegree()) {
            return false;
        }
    }

    return true;
}

/*
 * Addition: term0 + term1 + ... + termn
 */
Plus::Plus()
{
}

Plus::Plus(Term *lhs, Term *rhs)
{
    add(lhs);
    add(rhs);
}

Plus::Plus(const Plus &other)
{
	terms.reserve(other.terms.size());
	for (auto term : other.terms) {
		terms.push_back(term->clone());
	}
}

Plus::~Plus()
{
}

void Plus::add(Term *term)
{
    auto mul = dynamic_cast<Mul *>(term);
    auto plus = dynamic_cast<Plus *>(term);
    if(mul == nullptr && plus == nullptr) {
        auto temp = new Mul();
        temp->add(term);
        term = temp;
    }
    terms.push_back(term);
}

/* TODO: Maybe better to evaluate in a tree-like fashion so
 * we avoid round off errors? F.ex. when we have 4 terms:
 * temp0 = term0 + term1
 * temp1 = term2 + term3
 * return temp0 + temp1
 * instead of
 * temp = 0.0
 * temp += termi for all i
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

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
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

void Plus::sort()
{
    for(auto &term : terms) {
        term->sort();
    }

    std::sort(terms.begin(), terms.end(), plusCompare);
}

// Find all terms that are Plus *, and "steal" their children.
// This can be done because x0 + (x1 + x2) = x0 + x1 + x2
void Plus::steal_children() {
    // Store the current size. It will increase if we need to steal terms,
    // and there's no need to see if those need stealing, as that has been
    // done in the concatenate loop above.
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

Term *Plus::concatenate()
{
    MultiTerm::concatenate();

    sort();

//    std::cout << "Before Plus::concatenate: ";
//    pretty_text(std::cout);
//    std::cout << std::endl;
//    std::cout << std::endl;

    int term_size = terms.size();
    for(int i = 0, new_i = 0; i < term_size; ++i) {
        new_i = i;
        // All children of a Plus are of type Mul if the Plus is in proper form
        auto cur = dynamic_cast<Mul *>(terms.at(i));

        // See if cur and next can be summed:
        // cur == 1.1*y^1*x^2
        // next == -0.3*y^1*x^2
        // new cur = 0.8*y^1*x^2
        if(i+1 < term_size) {
            auto next = dynamic_cast<Mul *>(terms.at(i+1));

            if(equal(cur, next)) {
                auto newCoefficient = cur->getCoefficient() + next->getCoefficient();
                cur->setCoefficient(newCoefficient);

                delete next;
                terms.erase(terms.begin() + i+1);
                --term_size;
                new_i = i-1;
            }
        }

        // Remove terms that have a coefficient of 0
        if(cur->getCoefficient() == 0.0) {
            delete cur;
            terms.erase(terms.begin() + i);
            --term_size;
            new_i = i-1;
        }

        i = new_i;
    }

//    std::cout << "After Plus::concatenate: ";
//    pretty_text(std::cout);
//    std::cout << std::endl;
//    std::cout << std::endl;
//    std::cout << std::endl;

    if(terms.size() == 1) {
        auto temp = terms.at(0);
        terms.at(0) = nullptr;
        return temp;
    }

    if(terms.size() == 0) {
        return new Mul(0.0);
    }

    return this;
}

Term *Plus::flatten()
{
    MultiTerm::flatten();

    steal_children();

    return this;
}

double Plus::getConstDegree() const
{
    double highest = 0.0;
    for(auto &term : terms) {
        if(term->isConstDegree()) {
            if(term->getConstDegree() > highest) {
                highest = term->getConstDegree();
            }
        }
    }

    return highest;
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


/*
 * Multiplication: coefficient * term0 * term1 * ... * termn
 */
Mul::Mul()
    : Mul(1.0)
{
}

Mul::Mul(double coefficient)
    : coefficient(coefficient)
{
}

Mul::Mul(Term *lhs, Term *rhs)
    : coefficient(1.0)
{
    add(lhs);
    add(rhs);
}

Mul::Mul(const Mul &other)
	: coefficient(other.coefficient)
{
	terms.reserve(other.terms.size());
	for (auto term : other.terms) {
		terms.push_back(term->clone());
	}
}

Mul::~Mul()
{
}

void Mul::add(Term *term)
{
    auto exp = dynamic_cast<Exp *>(term);
    auto mul = dynamic_cast<Mul *>(term);
    auto plus = dynamic_cast<Plus *>(term);
    if(exp == nullptr && mul == nullptr && plus == nullptr) {
        term = new Exp(term, new Mul(1.0));
    }
    terms.push_back(term);
}

void Mul::multiply(double coefficient)
{
    this->coefficient *= coefficient;
}

/* TODO: Maybe better to evaluate in a tree-like fashion so
 * we avoid round off errors? F.ex. when we have 4 terms:
 * temp0 = term0 * term1
 * temp1 = term2 * term3
 * return temp0 * temp1
 * instead of
 * temp = 0.0
 * temp *= termi for all i
 */
double Mul::eval(const std::vector<double> &x) const
{
    auto tot = coefficient;
    for(auto &term : terms) {
        tot *= term->eval(x);
    }
    return tot;
}

Term *Mul::derivative(Var x) const
{
    if(terms.size() == 0) {
        return new Mul(0.0);
    }

    Plus *result = new Plus();
    for(int i = 0; i < terms.size(); ++i) {

        Mul *temp = new Mul(coefficient);
        for(int j = 0; j < terms.size(); ++j) {
            if(i == j) {
                temp->add(terms.at(j)->derivative(x));
            } else {
                temp->add(terms.at(j)->clone());
            }
        }

        result->add(temp);
    }

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Mul::clone() const
{
    auto result = new Mul(coefficient);
    result->terms.reserve(terms.size());
    for(auto &term : terms) {
        result->terms.push_back(term->clone());
    }
    return result;
}

void Mul::sort()
{
    for(auto &term : terms) {
        term->sort();
    }

    std::sort(terms.begin(), terms.end(), mulCompare);
}

// Find all terms that are Mul *, and "steal" their children.
// This can be done because x0 * (x1 * x2) = x0 * x1 * x2
void Mul::steal_children() {
    // Store the current size. It will increase if we need to steal terms,
    // and there's no need to see if those need stealing, as that has been
    // done in the concatenate loop above.
    size_t size = terms.size();
    for(int i = size - 1; i >= 0; --i) {
        // See if the term is a Mul, and if so, we "steal" all its children
        auto mul = dynamic_cast<Mul *>(terms.at(i));
        if(mul != nullptr) {
            this->coefficient *= mul->coefficient;
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
    if(terms.size() <= 0) {
        return this;
    }
	for (int i = terms.size() - 1; i >= 0; --i) {
		auto term = terms.at(i);
        auto plus = dynamic_cast<Plus *>(term);

        if(plus != nullptr) {
            auto result = new Plus();

            auto plus_terms = plus->getTerms();

            // Important to erase this _before_ the call to
            // this->clone below, else the plus will get
            // cloned too and start an infinite loop
            terms.erase(terms.begin() + i);

            for(auto &plus_term : plus_terms) {
                auto temp = new Mul();
                temp->add(plus_term->clone());
                temp->add(this->clone());
                result->add(temp);
            }

            // Delete the Plus that we copied from
            delete plus;

            auto temp = result->flatten();
            if(temp != result) {
                delete result;
            }
            return temp;
        }
    }

    return this;
}

Term *Mul::concatenate()
{
    MultiTerm::concatenate();

    // This must be done because when we concatenate the children,
    // some children may be constant and replace themselves with
    // a Mul(constantValue)
    steal_children();

    sort();

//    std::cout << "Before Mul::concatenate: ";
//    pretty_text(std::cout);
//    std::cout << std::endl;

    /*
     * Concatenate
     */
    int term_size = terms.size() - 1;
    for(int i = 0; i < term_size; ++i) {
        // All children of a Mul should be of type Exp if the Mul is in proper form
        auto cur = dynamic_cast<Exp *>(terms.at(i));
        auto next = dynamic_cast<Exp *>(terms.at(i+1));
        if(equal(cur->getBase(), next->getBase())) {
            auto newExponent = new Plus();
            newExponent->add(cur->getExponent());
            newExponent->add(next->getExponent());
            auto temp = newExponent->simplify();
            if(temp != newExponent) {
                delete newExponent;
            }
            cur->setExponent(temp);


            // Set to nullptr to avoid deleting
            next->setExponent(nullptr);
            delete next;
            terms.erase(terms.begin() + i+1);

            if(cur->isConstDegree() && cur->getConstDegree() == 0.0) {
                delete cur;
                terms.erase(terms.begin() + i);
                --term_size;
            }
            --term_size;
            --i;
        }
    }

//    std::cout << "After Mul::concatenate: ";
//    pretty_text(std::cout);
//    std::cout << std::endl;

    return this;
}

// (x+y)((x+y)*x) = (x+y)(x+y)*x
// = x*(x+y)*x + y*(x+y)*x
Term *Mul::flatten()
{
    MultiTerm::flatten();

    steal_children();

    return multiply_by_plus();
}

bool Mul::isConstant() const
{
    return terms.size() == 0;
}

double Mul::getConstDegree() const
{
    double degree = 0.0;
    for(auto &term : terms) {
        degree += term->getConstDegree();
    }
    return degree;
}

void Mul::pretty_text(std::ostream &out) const
{
    if(terms.size() == 0) {
        out << coefficient;
    }
    else {
        if(coefficient != 1.0) {
            if(coefficient == -1.0) {
                out << "-";
            } else {
                out << coefficient << "*";
            }
        }
    }

    for(auto it = terms.cbegin(); it != terms.cend(); ++it) {
        (*it)->pretty_text(out);
        if(it + 1 < terms.cend()) {
            out << "*";
        }
    }
}


/*
 * Variable
 */
Var::Var(int varNum)
    : Var(varNum, nullptr)
{
}

Var::Var(int varNum, const char *name)
        : varNum(varNum),
          name(strdup(name))
{
}

Var::Var(const Var &var)
        : Var(var.varNum, var.name)
{
}

Var &Var::operator=(Var &&other)
{
    name = std::move(other.name);
    other.name = nullptr;
    return *this;
}

Var::~Var()
{
    free((void *) name);
}

double Var::eval(const std::vector<double> &x) const
{
    return x.at(varNum);
}

Term *Var::derivative(Var x) const
{
    if(x.getVarNum() == varNum) {
        return new Mul(1);
    }
    else
    {
        return new Mul(0);
    }
}

Term *Var::clone() const
{
    return new Var(varNum, name);
}

void Var::getVariables(std::set<Var> &vars) const {
    vars.insert(*this);
}

void Var::pretty_text(std::ostream &out) const
{
    if(name != nullptr) {
        out << name;
    } else {
        out << "x" << varNum;
    }
}


/*
 * Constant value
 */
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
    return new Mul(0);
}

Term *Const::clone() const
{
    return new Const(val);
}

void Const::pretty_text(std::ostream &out) const
{
    out << val;
}


/*
 * Power: x^y
 * TODO: Rename to Pow?
 */
Exp::Exp(Term *base, Term *exponent)
    : base(base),
      exponent(exponent)
{
}

Exp::~Exp()
{
    delete base;
    delete exponent;
}

double Exp::eval(const std::vector<double> &x) const
{
    auto baseVal = base->eval(x);
    auto expVal = exponent->eval(x);

    return std::pow(baseVal, expVal);
}

// https://en.wikipedia.org/wiki/Differentiation_rules#Generalized_power_rule
Term *Exp::derivative(Var x) const
{
    auto f = base;
    auto g = exponent;

    // f^g
    auto fg = this;

    // 1/f
    auto one_over_f = new Exp(f->clone(), new Mul(-1));

    // df * g/f
    auto df_mul_g_over_f = new Mul();
    df_mul_g_over_f->add(f->derivative(x));
    df_mul_g_over_f->add(g->clone());
    df_mul_g_over_f->add(one_over_f);

    // ln(f)
    auto lnf = new Log(f->clone());

    // dg * ln(f)
    auto dg_ln_f = new Mul(g->derivative(x), lnf);

    auto temp = new Plus(df_mul_g_over_f, dg_ln_f);

    auto result = new Mul(fg->clone(), temp);

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Exp::clone() const
{
    return new Exp(base->clone(), exponent->clone());
}

Term *Exp::concatenate()
{
    auto temp = base->concatenate();
    if(temp != base) {
        delete base;
        base = temp;
    }

    temp = exponent->concatenate();
    if(temp != exponent) {
        delete exponent;
        exponent = temp;
    }

    if(base->isConstant()) {
        // a^b
        if(exponent->isConstant()) {
            std::vector<double> dummy;
            auto constVal = eval(dummy);
            return new Mul(constVal);
        }

        // 0^x
        // Ignore the case of x == 0
        if(base->getConstValue() == 0.0) {
            return new Mul(0.0);
        }
    }
    // x^0
    else if(exponent->isConstant() && exponent->getConstValue() == 0.0) {
        return new Mul(1.0);
    }

    return this;
}

Term *Exp::flatten()
{
    auto temp = base->flatten();
    if(temp != base) {
        delete base;
        base = temp;
    }

    temp = exponent->flatten();
    if(temp != exponent) {
        delete exponent;
        exponent = temp;
    }

    return this;
}

void Exp::getVariables(std::set<Var> &vars) const {
    base->getVariables(vars);
    exponent->getVariables(vars);
}

void Exp::sort()
{
    base->sort();
    exponent->sort();
}

void Exp::pretty_text(std::ostream &out) const
{
    if(isConstDegree() && getConstDegree() == 1.0) {
        base->pretty_text(out);
    } else {
        out << "(";
        base->pretty_text(out);
        out << "^";
        exponent->pretty_text(out);
        out << ")";
    }
}

FuncTerm::FuncTerm(Term *arg)
    : arg(arg)
{
}

FuncTerm::~FuncTerm()
{
    delete arg;
}

Term *FuncTerm::concatenate()
{
    auto temp = arg->concatenate();
    if(temp != arg) {
        delete arg;
        arg = temp;
    }

    if(arg->isConstant()) {
        std::vector<double> dummy;
        auto constantVal = eval(dummy);
        return new Mul(constantVal);
    }

    return this;
}

void FuncTerm::getVariables(std::set<Var> &vars) const
{
    arg->getVariables(vars);
}


/*
 * Logarithm with base n: log_n(x)
 */
Log::Log(Term *arg)
    : Log(std::exp(1.0), arg)
{
}

Log::Log(const Term &arg)
    : Log(arg.clone())
{
}

Log::Log(double base, Term *arg)
        : base(base),
          FuncTerm(arg)
{
}

Log::Log(double base, const Term &arg)
    : Log(base, arg.clone())
{
}

Log::~Log()
{
}

double Log::eval(const std::vector<double> &x) const
{
    return log(base, arg->eval(x));
}

Term *Log::derivative(Var x) const
{
    // d(log(f(x)))/dx = log_base(e) * (df(x) / f(x))
    // Actually it is 1/ln(base) * (df(x) / f(x)), but because
    // log_base(e) = ln(e) / ln(base) = 1/ln(base) we avoid the extra division

    auto one_over_arg = new Exp(arg->clone(), new Mul(-1));

    auto result = new Mul(log(base, std::exp(1.0)));
    result->add(arg->derivative(x));
    result->add(one_over_arg);

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Log::clone() const
{
    return new Log(base, arg->clone());
}

void Log::pretty_text(std::ostream &out) const
{
    if(base == std::exp(1.0)) {
        out << "ln(";
    } else {
        out << "log_" << base << "(";
    }
    arg->pretty_text(out);
    out << ")";
}


/*
 * Sinus: sin(x)
 */
Sin::Sin(Term *arg)
    : FuncTerm(arg)
{
}

Sin::Sin(const Term &arg)
    : Sin(arg.clone())
{
}

Sin::~Sin()
{
}

double Sin::eval(const std::vector<double> &x) const
{
    return std::sin(arg->eval(x));
}

Term *Sin::derivative(Var x) const
{
    // dsin(f(x))/dx = df(x)/dx * cos(f(x))
    auto result = new Mul(arg->derivative(x), new Cos(arg->clone()));

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Sin::clone() const
{
    return new Sin(arg->clone());
}

void Sin::pretty_text(std::ostream &out) const
{
    out << "sin(";
    arg->pretty_text(out);
    out << ")";
}


/*
 * Cosinus: cos(x)
 */
Cos::Cos(Term *arg)
    : FuncTerm(arg)
{
}

Cos::Cos(const Term &arg)
    : Cos(arg.clone())
{
}

Cos::~Cos()
{
}

double Cos::eval(const std::vector<double> &x) const
{
    return std::cos(arg->eval(x));
}

Term *Cos::derivative(Var x) const
{
    // dcos(f(x))/dx = -1 * df(x)/dx * sin(f(x))
    auto result = new Mul(arg->derivative(x), new Sin(arg->clone()));
    result->multiply(-1.0);

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Cos::clone() const
{
    return new Cos(arg->clone());
}

void Cos::pretty_text(std::ostream &out) const
{
    out << "cos(";
    arg->pretty_text(out);
    out << ")";
}


/*
 * Tangential: tan(x)
 */
Tan::Tan(Term *arg)
    : FuncTerm(arg)
{
}

Tan::Tan(const Term &arg)
    : Tan(arg.clone())
{
}

Tan::~Tan()
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

    auto temp = new Exp(denominator, new Mul(-1));

    auto result = new Mul(arg->derivative(x), temp);

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *Tan::clone() const
{
    return new Tan(arg->clone());
}

void Tan::pretty_text(std::ostream &out) const
{
    out << "tan(";
    arg->pretty_text(out);
    out << ")";
}


/*
 * Exponential function e^x
 */
E::E(Term *arg)
        : FuncTerm(arg)
{
}

E::E(const Term &arg)
    : E(arg.clone())
{
}

E::~E()
{
}

double E::eval(const std::vector<double> &x) const
{
    return std::exp(arg->eval(x));
}

Term *E::derivative(Var x) const
{
    // de(f(x))/dx = df(x) * e(f(x))
    auto result = new Mul(arg->derivative(x), new E(arg->clone()));

    auto simplified = result->simplify();
    if(simplified != result) {
        delete result;
    }
    return simplified;
}

Term *E::clone() const
{
    return new E(arg->clone());
}

void E::pretty_text(std::ostream &out) const
{
    out << "exp(";
    arg->pretty_text(out);
    out << ")";
}


// Returns true if the terms are of equal type, and their children are of equal type
bool equal(Term *lhs, Term *rhs) {
//    std::cout << "Comparing ";
//    lhs->pretty_text(std::cout);
//    std::cout << " to ";
//    rhs->pretty_text(std::cout);

    bool termsEqual = true;

    auto lconst = dynamic_cast<Const *>(lhs);
    auto lvar = dynamic_cast<Var *>(lhs);
    auto lexp = dynamic_cast<Exp *>(lhs);
    auto lmul = dynamic_cast<Mul *>(lhs);
    auto lplus = dynamic_cast<Plus *>(lhs);
    auto le = dynamic_cast<E *>(lhs);
    auto llog = dynamic_cast<Log *>(lhs);
    auto lsin = dynamic_cast<Sin *>(lhs);
    auto lcos = dynamic_cast<Cos *>(lhs);
    auto ltan = dynamic_cast<Tan *>(lhs);

    auto rconst = dynamic_cast<Const *>(rhs);
    auto rvar = dynamic_cast<Var *>(rhs);
    auto rexp = dynamic_cast<Exp *>(rhs);
    auto rmul = dynamic_cast<Mul *>(rhs);
    auto rplus = dynamic_cast<Plus *>(rhs);
    auto re = dynamic_cast<E *>(rhs);
    auto rlog = dynamic_cast<Log *>(rhs);
    auto rsin = dynamic_cast<Sin *>(rhs);
    auto rcos = dynamic_cast<Cos *>(rhs);
    auto rtan = dynamic_cast<Tan *>(rhs);

    // Sort by degree increasing
    double ldegree = lhs->getConstDegree();
    double rdegree = rhs->getConstDegree();

    // Constant degrees before non constant degrees
    if(lhs->isConstDegree() && rhs->isConstDegree()) {
//        std::cout << "both are const degrees!";
//        std::cout << "ldegree: " << ldegree << ", rdegree: " << rdegree;
        termsEqual = ldegree == rdegree;
    }

    // If we have come this far, it means the degree is equal

    // lpri < rpri means lhs should be before rhs
    int lpri =
            (lconst != nullptr) * 1
            + (lvar != nullptr) * 2
            + (lexp != nullptr) * 3
            + (lmul != nullptr) * 4
            + (lplus != nullptr) * 5
            + (le != nullptr) * 6
            + (llog != nullptr) * 7
            + (lsin != nullptr) * 8
            + (lcos != nullptr) * 9
            + (ltan != nullptr) * 10
            + (lplus != nullptr) * 11;
    int rpri =
            (rconst != nullptr) * 1
            + (rvar != nullptr) * 2
            + (rexp != nullptr) * 3
            + (rmul != nullptr) * 4
            + (rplus != nullptr) * 5
            + (re != nullptr) * 6
            + (rlog != nullptr) * 7
            + (rsin != nullptr) * 8
            + (rcos != nullptr) * 9
            + (rtan != nullptr) * 10
            + (rplus != nullptr) * 11;

    if(termsEqual) {
        if(lpri != rpri) {
            termsEqual = false;
        }

        else if(lvar != nullptr) {
            termsEqual = lvar->getVarNum() == rvar->getVarNum();
        }

        else if(lexp != nullptr) {
            termsEqual = equal(lexp->getBase(), rexp->getBase())
                    && equal(lexp->getExponent(), rexp->getExponent());
        }

        else if(lmul != nullptr) {
            auto lterms = lmul->getTerms();
            auto rterms = rmul->getTerms();
            if(lterms.size() != rterms.size()) {
                termsEqual = false;
            }
            else
            {
                int n_terms = lterms.size();
                for(int i = n_terms-1; i >= 0; --i) {
                    if(!equal(lterms.at(i), rterms.at(i))) {
                        termsEqual = false;
                        break;
                    }
                }
            }
        }

        else if(le != nullptr || lsin != nullptr || lcos != nullptr || ltan != nullptr) {
            auto lFunc = dynamic_cast<FuncTerm *>(lhs);
            auto rFunc = dynamic_cast<FuncTerm *>(rhs);
            return equal(lFunc->getArg(), rFunc->getArg());
        }

        else if(llog != nullptr) {
            if(llog->getBase() != rlog->getBase()) {
                termsEqual = false;
            }
            else {
                termsEqual = equal(llog->getArg(), rlog->getArg());
            }
        }

        else {
            assert(lplus != nullptr && rplus != nullptr);

            if(lplus->getTerms().size() != rplus->getTerms().size()) {
                termsEqual = false;
            }

            for(int i = lplus->getTerms().size() - 1; i >= 0; --i) {
                auto lterm = lplus->getTerms().at(i);
                auto rterm = rplus->getTerms().at(i);
                if(!equal(lterm, rterm)) {
                    termsEqual = false;
                    break;
                }
            }
        }
    }


//    if(termsEqual) {
//        std::cout << ", Equal!";
//    } else {
//        std::cout << ", Not equal!";
//    }
//    std::cout << std::endl;

    return termsEqual;
}

// Returns true if lhs should be before rhs in a sorted range
bool plusCompare(Term *lhs, Term *rhs) {
    auto lmul = dynamic_cast<Mul *>(lhs);
    auto rmul = dynamic_cast<Mul *>(rhs);

    assert(lmul != nullptr && rmul != nullptr);

    double ldegree = lhs->getConstDegree();
    double rdegree = rhs->getConstDegree();

    // Constant degrees before non constant degrees
    if(lhs->isConstDegree() != rhs->isConstDegree()) {
        return lhs->isConstDegree();
    }

    // Sort by degree descending
    if(ldegree != rdegree) {
        return ldegree > rdegree;
    }

    // If we have come this far, it means the degree is equal
    auto lterms = lmul->getTerms();
    auto rterms = rmul->getTerms();
    if(lterms.size() != rterms.size()) {
        return lterms.size() < rterms.size();
    }
    else
    {
        int n_terms = lterms.size();
        for(int i = n_terms-1; i >= 0; --i) {
            if(compare(lterms.at(i), rterms.at(i))) {
                return true;
            }
            if(compare(rterms.at(i), lterms.at(i))) {
                return false;
            }
        }

        // Equal, default to false to avoid stopping further comparison in the caller
        return false;
    }
}

// Returns true if lhs should be before rhs in a sorted range
bool mulCompare(Term *lhs, Term *rhs) {
    auto lexp = dynamic_cast<Exp *>(lhs);
    auto rexp = dynamic_cast<Exp *>(rhs);

    assert(lexp != nullptr && rexp != nullptr);

    if(compare(lexp->getBase(), rexp->getBase()))
    {
        return true;
    }
    else if(compare(rexp->getBase(), lexp->getBase())) {
        return false;
    }
    else {
        if(compare(lexp->getExponent(), rexp->getExponent()))
        {
            return true;
        }
        else if(compare(rexp->getExponent(), lexp->getExponent())) {
            return false;
        }
        else
        {
            // Equal, default to false to avoid stopping further comparison in the caller
            return false;
        }
    }
}

// Returns true if lhs should be before rhs in a sorted range
bool compare(Term *lhs, Term *rhs) {
    auto lconst = dynamic_cast<Const *>(lhs);
    auto lvar = dynamic_cast<Var *>(lhs);
    auto lexp = dynamic_cast<Exp *>(lhs);
    auto lmul = dynamic_cast<Mul *>(lhs);
    auto lplus = dynamic_cast<Plus *>(lhs);
    auto le = dynamic_cast<E *>(lhs);
    auto llog = dynamic_cast<Log *>(lhs);
    auto lsin = dynamic_cast<Sin *>(lhs);
    auto lcos = dynamic_cast<Cos *>(lhs);
    auto ltan = dynamic_cast<Tan *>(lhs);

    auto rconst = dynamic_cast<Const *>(rhs);
    auto rvar = dynamic_cast<Var *>(rhs);
    auto rexp = dynamic_cast<Exp *>(rhs);
    auto rmul = dynamic_cast<Mul *>(rhs);
    auto rplus = dynamic_cast<Plus *>(rhs);
    auto re = dynamic_cast<E *>(rhs);
    auto rlog = dynamic_cast<Log *>(rhs);
    auto rsin = dynamic_cast<Sin *>(rhs);
    auto rcos = dynamic_cast<Cos *>(rhs);
    auto rtan = dynamic_cast<Tan *>(rhs);

    double ldegree = lhs->getConstDegree();
    double rdegree = rhs->getConstDegree();

    // Constant degrees before non constant degrees
    if(lhs->isConstDegree() != rhs->isConstDegree()) {
        return lhs->isConstDegree();
    }

    if(ldegree != rdegree) {
        return ldegree < rdegree;
    }

    // If we have come this far, it means the degree is equal

    // lpri < rpri means lhs should be before rhs
    int lpri =
            (lconst != nullptr) * 1
            + (lvar != nullptr) * 2
            + (lexp != nullptr) * 3
            + (lmul != nullptr) * 4
            + (lplus != nullptr) * 5
            + (le != nullptr) * 6
            + (llog != nullptr) * 7
            + (lsin != nullptr) * 8
            + (lcos != nullptr) * 9
            + (ltan != nullptr) * 10
            + (lplus != nullptr) * 11;
    int rpri =
            (rconst != nullptr) * 1
            + (rvar != nullptr) * 2
            + (rexp != nullptr) * 3
            + (rmul != nullptr) * 4
            + (rplus != nullptr) * 5
            + (re != nullptr) * 6
            + (rlog != nullptr) * 7
            + (rsin != nullptr) * 8
            + (rcos != nullptr) * 9
            + (rtan != nullptr) * 10
            + (rplus != nullptr) * 11;

    if(lpri != rpri) {
        return lpri < rpri;
    }

    else if(lconst != nullptr) {
        return lconst->getVal() < rconst->getVal();
    }

    else if(lvar != nullptr) {
        return lvar->getVarNum() < rvar->getVarNum();
    }

    else if(lexp != nullptr) {
        if(compare(lexp->getBase(), rexp->getBase()))
        {
            return true;
        }
        else if(compare(rexp->getBase(), lexp->getBase())) {
            return false;
        }
        else {
            if(compare(lexp->getExponent(), rexp->getExponent()))
            {
                return true;
            }
            else if(compare(rexp->getExponent(), lexp->getExponent())) {
                return false;
            }
            else {

                // Equal, default to false to avoid stopping further comparison
                return false;
            }
        }
    }

    else if(lmul != nullptr) {
        auto lterms = lmul->getTerms();
        auto rterms = rmul->getTerms();
        if(lterms.size() != rterms.size()) {
            return lterms.size() < rterms.size();
        }
        else
        {
            int n_terms = lterms.size();
            for(int i = n_terms-1; i >= 0; --i) {
                if(compare(lterms.at(i), rterms.at(i))) {
                    return true;
                }
                if(compare(rterms.at(i), lterms.at(i))) {
                    return false;
                }
            }

            // Equal, default to false to avoid stopping further comparison
            return false;
        }
    }

    else if(le != nullptr || lsin != nullptr || lcos != nullptr || ltan != nullptr) {
        auto lFunc = dynamic_cast<FuncTerm *>(lhs);
        auto rFunc = dynamic_cast<FuncTerm *>(rhs);
        return compare(lFunc->getArg(), rFunc->getArg());
    }

    else if(llog != nullptr) {
        if(llog->getBase() < rlog->getBase()) {
            return true;
        }
        else {
            return compare(llog->getArg(), rlog->getArg());
        }
    }

    else {
        assert(lplus != nullptr && rplus != nullptr);

        if(lplus->getTerms().size() != rplus->getTerms().size()) {
            return lplus->getTerms().size() < rplus->getTerms().size();
        }

        for(int i = lplus->getTerms().size() - 1; i >= 0; --i) {
            auto lterm = lplus->getTerms().at(i);
            auto rterm = rplus->getTerms().at(i);
            if(compare(lterm, rterm)) {
                return true;
            }
            if(compare(rterm, lterm)) {
                return false;
            }
        }

        // Equal, default to false to avoid stopping further comparison
        return false;
    }
}

bool operator<(const Var &lhs, const Var &rhs) {
    return lhs.getVarNum() < rhs.getVarNum();
}

std::ostream &operator<<(std::ostream &out, const Term &term)
{
    term.pretty_text(out);
    return out;
}

Plus operator+(const Term &lhs, const Term &rhs)
{
    return Plus(lhs.clone(), rhs.clone());
}

Plus operator+(const Term &lhs, double rhs)
{
    return Plus(lhs.clone(), new Mul(rhs));
}

Plus operator+(double lhs, const Term &rhs)
{
    return Plus(new Mul(lhs), rhs.clone());
}


Plus operator-(const Term &lhs, const Term &rhs)
{
    auto temp = new Mul(-1.0);
    temp->add(rhs.clone());
    return Plus(lhs.clone(), temp);
}

Plus operator-(const Term &lhs, double rhs)
{
    return Plus(lhs.clone(), new Mul(-rhs));
}

Plus operator-(double lhs, const Term &rhs)
{
    auto temp = new Mul(lhs);
    auto temp2 = new Mul(-1.0);
    temp2->add(rhs.clone());
    return Plus(temp, temp2);
}


Mul operator*(const Term &lhs, const Term &rhs)
{
    return Mul(lhs.clone(), rhs.clone());
}

Mul operator*(double lhs, const Term &rhs)
{
    auto temp = Mul(lhs);
    temp.add(rhs.clone());
    return temp;
}

Mul operator*(const Term &lhs, double rhs)
{
    auto temp = Mul();
    temp.multiply(rhs);
    temp.add(lhs.clone());
    return temp;
}


Mul operator/(const Term &lhs, const Term &rhs)
{
    return Mul(lhs.clone(), new Exp(rhs.clone(), new Mul(-1.0)));
}

Mul operator/(double lhs, const Term &rhs)
{
    auto temp = new Mul(lhs);
    auto temp2 = new Mul(-1.0);
    return Mul(temp, new Exp(rhs.clone(), temp2));
}

Mul operator/(const Term &lhs, double rhs)
{
    return Mul(lhs.clone(), new Exp(new Mul(rhs), new Mul(-1.0)));
}


Exp operator^(const Term &lhs, const Term &rhs)
{
    return Exp(lhs.clone(), rhs.clone());
}

Exp operator^(const Term &base, double exp)
{
    return Exp(base.clone(), new Mul(exp));
}

Exp operator^(double base, const Term &exp)
{
    return Exp(new Mul(base), exp.clone());
}
