#pragma once
#include "../sample_src/Shekel/ShekelProblem.hpp"
#include <functional>

extern double kShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double aShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double cShekel[NUM_SHEKEL_PROBLEMS][NUM_SHEKEL_COEFF];
extern double minShekel[NUM_SHEKEL_PROBLEMS][2];
extern double maxShekel[NUM_SHEKEL_PROBLEMS][2];
extern double lConstantShekel[NUM_SHEKEL_PROBLEMS];

struct ConsTrial {
    double x, z;
    int index;
    ConsTrial(double _x = 0, double _z = 0, int _index = 1) {
        x = _x;
        z = _z;
        index = _index;
    }
};

class MyOptFunction {
public:
    virtual double operator() (double x) = 0;
    virtual double Compute(double x) = 0;
};

class MyConstrainedProblem
{
protected:
    double OptimalPoint;
    double OptimalValue;
    uint m; // Number of constraints
    std::vector<MyOptFunction* > g_vec;
    MyOptFunction* func;
    double LoBound;
    double UpBound;
public:
    MyConstrainedProblem(
        double _OptimalPoint, double _OptimalValue,
        double _LoBound, double _UpBound,
        uint _m, // Quantity of constraints
        const std::vector<MyOptFunction*>& _g_vec,
        MyOptFunction* _func
    ) : OptimalPoint(_OptimalPoint), OptimalValue(_OptimalValue), LoBound(_LoBound),
        UpBound(_UpBound), m(_m), g_vec(_g_vec), func(_func) {} // Yeah, I know...

    double ComputeFunction(double x) const;
    double ComputeGjConstr(uint i, double x) const; // i==m return ComputeFunction

    /// Get global minimizer
    double GetOptimumPoint() const;
    /// Get global minimum value
    double GetOptimumValue() const;
    uint GetNumberofConstr() const;
};

enum MyConstrPrType{SheckelOnly, HillOnly, RandomSheckelOnly,RandomHillOnly, RandomHillSheckel};



class MySheckelFunction : public MyOptFunction {
private:
    uint index;
    double delta;
public:
    MySheckelFunction(uint _index, double _delta) {
        index = _index; delta = _delta;
    }
    double operator() (double x) {
        double res = 0.0;

        for (int j = 0; j < NUM_SHEKEL_COEFF; j++)
        {
            res = res - 1 / (kShekel[index][j] * pow(x - aShekel[index][j], 2.0) +
                cShekel[index][j]);
        }
        return res - delta;
    }
    double Compute(double x) { return this->operator()(x); }
};

class MyConstrainedProblemGenerator {
protected:
    std::pair<double, double> FindOptimalPoint();
public:
    MyConstrainedProblem* Generate(MyConstrPrType type, uint m, double delta = 0.1, uint seed = 0);
    vector<MyConstrainedProblem*> GenerateNProblems(uint n, MyConstrPrType type, uint m, double delta = 0.01, uint seed = 0);
};