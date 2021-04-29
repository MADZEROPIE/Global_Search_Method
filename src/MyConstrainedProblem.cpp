#include "MyConstrainedProblem.h"
#include <iostream>

double MyConstrainedProblem::ComputeFunction(double x) const {
    return func->Compute(x);
}

double MyConstrainedProblem::ComputeGjConstr(uint i, double x) const {
    double res;
    if (i >= 0 && i < m) {
        res = this->g_vec[i]->Compute(x);
    }
    else if (i == m) {
        res = this->func->Compute(x);
    }
    else throw("Incorrect index");
    return res;
}

double MyConstrainedProblem::GetOptimumPoint() const {
    return OptimalPoint;
}

double MyConstrainedProblem::GetOptimumValue() const {
    return OptimalValue;
}

uint MyConstrainedProblem::GetNumberofConstr() const {
    return m;
}


double MyConstrainedProblemGenerator::findMinMax(const vector<MyOptFunction*>& g_vec, double LoBound, double UpBound, uint64_t steps) {
    int m = g_vec.size();
    if (m == 0) return 0;

    double minMax = INFINITY;
    double h = (UpBound - LoBound) / steps;
    for (double x = LoBound; x <= UpBound; x += h) {
        double z_max = g_vec[0]->Compute(x);
        for (int i = 1; i < m; ++i) {
            double z_tmp= g_vec[i]->Compute(x);
            if (z_tmp > z_max) { z_max = z_tmp; }
        }
        if (z_max < minMax) {
            minMax = z_max;
        }
    }
    return minMax;
}

MyConstrainedProblem* MyConstrainedProblemGenerator::Generate(MyConstrPrType type, uint m, double delta, uint seed) {
    MyConstrainedProblem* problem = nullptr;

    double OptimalPoint;
    double OptimalValue;
    std::vector<MyOptFunction*> g_vec;
    MyOptFunction* func;
    double LoBound;
    double UpBound;
    std::mt19937 gen;
    gen.seed(seed);


    if (type == SheckelOnly) {
        LoBound = 0.0;
        UpBound = 10.0;
        std::vector<int> SheckelIndex(m + 1u);
        for (uint i = 0; i <= m; ++i) {
            SheckelIndex[i] = gen() % NUM_SHEKEL_PROBLEMS;
            for (uint j = 0; j < i; ++j) {
                if (SheckelIndex[i] == SheckelIndex[j]) {
                    --i; 
                    break;
                }
            }
        }

        for (uint i = 0; i < m; ++i) {
            uint ind = SheckelIndex[i];
            MySheckelFunction* f_ptr = new MySheckelFunction(SheckelIndex[i], 0.0);
            g_vec.push_back(f_ptr);
        }
        
        func = new MySheckelFunction(SheckelIndex[m], 0.0); 
    }
    else if (HillOnly) {
        LoBound = 0.0;
        UpBound = 1.0;
        std::vector<int> HillIndex(m + 1u);
        for (uint i = 0; i <= m; ++i) {
            HillIndex[i] = gen() % NUM_HILL_PROBLEMS;
            for (uint j = 0; j < i; ++j) {
                if (HillIndex[i] == HillIndex[j]) {
                    --i;
                    j = i;
                }
            }
        }

        for (uint i = 0; i < m; ++i) {
            MyOptFunction* f_ptr = new MyHillFunction(HillIndex[i], 0.0);
            g_vec.push_back(f_ptr);
        }
        int ind = HillIndex[m];
        func = new MyHillFunction(HillIndex[m], 0.0);
    }
    const int steps = 100000;
    double minMax = this->findMinMax(g_vec, LoBound, UpBound, steps);
    for (int i = 0; i < m; ++i) {
        if (type == HillOnly) {
            auto fptr = static_cast<MyHillFunction*>(g_vec[i]);
            fptr->delta = delta + minMax;
        }
        else if (type == SheckelOnly) {
            auto fptr = static_cast<MySheckelFunction*>(g_vec[i]);
            fptr->delta = delta + minMax;
        }
    }
    //--Findin' minimum
    int k = 0;
    const double h = (UpBound - LoBound) / steps;
    bool flag = false;
    double x = LoBound;
    for (; x <= UpBound && !flag; x += h) {
        double z;
        int i = 0;
        for (; i < m && (g_vec[i]->Compute(x) <= 0.0); ++i);
        if (i == m) {
            flag = true;
            OptimalPoint = x;
            OptimalValue = func->Compute(x);
            ++k;
        }
    }
    for (; x < UpBound; x += h) {
        double z;
        int i = 0;
        for (; i < m && (z = g_vec[i]->Compute(x)) <= 0.0; ++i);
        if (i == m) {
            ++k;
            z = func->Compute(x);
            if (z < OptimalValue) {
                OptimalPoint = x;
                OptimalValue = z;
            }
        }
    }
    std::cout << k * 100.0 / (steps + 1.0) << "%\n";
    if (flag == false) {
        throw "Something bad happened while generation";
    }
    problem = new MyConstrainedProblem(OptimalPoint, OptimalValue, LoBound, UpBound, m, g_vec, func);
    
    return problem;
}

vector<MyConstrainedProblem*> MyConstrainedProblemGenerator::GenerateNProblems(uint n, MyConstrPrType type, uint m, double delta, uint seed)
{
    vector<MyConstrainedProblem*> problemVec(n);
    std::mt19937 gen;
    gen.seed(seed);
    for (int i = 0; i < n; ++i) {
        std::cout << i << " ";
        problemVec[i] = MyConstrainedProblemGenerator::Generate(type, m, delta, gen());
    }
    return problemVec;
}
