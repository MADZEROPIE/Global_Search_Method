#include "MyConstrainedProblem.h"
#include <iostream>

double MyConstrainedProblem::ComputeFunction(double x) const
{
    return func->Compute(x);
}

double MyConstrainedProblem::ComputeGjConstr(uint i, double x) const
{
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

double MyConstrainedProblem::GetOptimumPoint() const
{
    return OptimalPoint;
}

double MyConstrainedProblem::GetOptimumValue() const
{
    return OptimalValue;
}

uint MyConstrainedProblem::GetNumberofConstr() const
{
    return m;
}


MyConstrainedProblem* MyConstrainedProblemGenerator::Generate(MyConstrPrType type, uint m, double delta, uint seed)
{
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
        try
        {
        LoBound = 0.0;
        UpBound = 10.0;
        std::vector<int> SheckelIndex(m + 1u);
        for (uint i = 0; i <= m; ++i) {
            SheckelIndex[i] = gen() % NUM_SHEKEL_PROBLEMS;
            for (uint j = 0; j < i; ++j) {
                if (SheckelIndex[i] == SheckelIndex[j]) {
                    --i; 
                    j = i;
                }
            }
        }
        double minmaxSh = maxShekel[SheckelIndex[0]][0];
        for (uint i = 1u; i < m; ++i) {
            if (minmaxSh > maxShekel[SheckelIndex[i]][0])
                minmaxSh = maxShekel[SheckelIndex[i]][0];
        }
        for (uint i = 0; i < m; ++i) {
            uint ind = SheckelIndex[i];
            MyOptFunction* f_ptr = new MySheckelFunction(SheckelIndex[i], delta + minmaxSh);
            g_vec.push_back(f_ptr);
        }
        
        func = new MySheckelFunction(SheckelIndex[m], 0.0); 

        //--Findin' minimum
        const int steps = 1000000;
        const double h = (UpBound - LoBound) / steps;
        bool flag = false;
        double x = LoBound;
        for (; x < UpBound && !flag; x += h) {
            double z;
            int i = 0;
            for (; i < m; ++i) {
                if (g_vec[i]->Compute(x) > 0.0) break;
            }
            if (i == m) {
                flag = true;
                OptimalPoint = x;
                OptimalValue = func->Compute(x);
            }
        }
        for (; x < UpBound; x += h) {
            double z;
            int i = 0;
            for (; i < m && (z = g_vec[i]->Compute(x)) <= 0.0; ++i);
            if (i == m) {
                z = func->Compute(x);
                if (z < OptimalValue) {
                    OptimalPoint = x;
                    OptimalValue = z;
                }
            }
        }
        if (flag == false) {
            for (auto el : SheckelIndex) {
                std::cout << el << "\n";
            }
            std::cout << "\n\n\n\n";
            throw SheckelIndex;
        }

        problem = new MyConstrainedProblem(OptimalPoint, OptimalValue, LoBound, UpBound, m, g_vec, func);
        //std::cout << OptimalPoint << " " << OptimalValue << '\n' << minShekel[SheckelIndex[m]][1] << " " << minShekel[SheckelIndex[m]][0] << '\n';
        }
        catch (std::vector<int>& SheckelIndex)
        {
            
        }
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
        double minmaxHl = maxHill[HillIndex[0]][0];
        for (uint i = 1u; i < m; ++i) {
            if (minmaxHl > maxHill[HillIndex[i]][0])
                minmaxHl = maxHill[HillIndex[i]][0];
        }
        for (uint i = 0; i < m; ++i) {
            uint ind = HillIndex[i];
            MyOptFunction* f_ptr = new MyHillFunction(HillIndex[i], delta + minmaxHl);
            g_vec.push_back(f_ptr);
        }
        int ind = HillIndex[m];
        func = new MyHillFunction(HillIndex[m], 0.0);

        //--Findin' minimum
        const int steps = 100000;
        const double h = (UpBound - LoBound) / steps;
        bool flag = false;
        double x = LoBound;
        for (; x < UpBound && !flag; x += h) {
            double z;
            int i = 0;
            for (; i < m && (g_vec[i]->Compute(x) <= 0.0); ++i);
            if (i == m) {
                flag = true;
                OptimalPoint = x;
                OptimalValue = func->Compute(x);
            }
        }
        for (; x < UpBound; x += h) {
            double z;
            int i = 0;
            for (; i < m && (z = g_vec[i]->Compute(x)) <= 0.0; ++i);
            if (i == m) {
                z = func->Compute(x);
                if (z < OptimalValue) {
                    OptimalPoint = x;
                    OptimalValue = z;
                }
            }
        }
        if (flag == false) {
            throw "Something bad happened while generation";
        }
        problem = new MyConstrainedProblem(OptimalPoint, OptimalValue, LoBound, UpBound, m, g_vec, func);
        //std::cout << OptimalPoint << " " << OptimalValue << '\n' << minHill[HillIndex[m]][1] << " " << minHill[HillIndex[m]][0] << '\n';
    }
    return problem;
}

vector<MyConstrainedProblem*> MyConstrainedProblemGenerator::GenerateNProblems(uint n, MyConstrPrType type, uint m, double delta, uint seed)
{
    vector<MyConstrainedProblem*> problemVec(n);
    for (int i = 0; i < n; ++i) {
        problemVec[i] = MyConstrainedProblemGenerator::Generate(type, m, delta, seed+i);
    }
    return problemVec;
}