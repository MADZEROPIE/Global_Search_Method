#pragma once
#include "IOptProblem.hpp"

#include <iostream>
#include <functional> //for std::function
#include <cmath>  //for math functions e.g. sin() or cos()
#include <vector> 
#include <algorithm>
#include <fstream>


#include <omp.h>

#define NUMTH 6
//About style of naming. There is no style.

using std::vector;
using std::pair;

//typedef pair<double, double> Trial;
struct TrialD {
    vector<double> x;
    double z;
    TrialD(const vector<double>& _x = { 0.0 }, double _z = 0.0) {
        x = _x;
        z = _z;
    }
};

class MinimazerD { //Don't ask me why... But only because using a sledge-hammer to crack a nut sounds fun.

protected:
    //std::function<double(double)> func; //Function that needs findin' minimum
    vector<double> a; //Beginning of the segment
    vector<double> b; //End of the segment

    double r; //Coefficient of the method

    double eps;
    IOptProblem* IOPPtr; //Problem that needs findin' minimum
    TrialD sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number
    int dim;

public:
    MinimazerD(IOptProblem* _IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500) {
        IOPPtr = _IOPPtr;
        IOPPtr->GetBounds(a, b);
        r = (_r > 1.0) ? _r : 2.0;
        eps = _eps;
        NMax = _NMax;
        dim = IOPPtr->GetDimension();
    }

    TrialD find_glob_min_fixed_index(std::vector<double>& x, int index, bool stop_crit = false) {
        if (index >= dim) return TrialD(x, IOPPtr->ComputeFunction(x));
        vector<TrialD> vec;
        auto lb = x;
        auto rb = x;
        lb[index] = a[index];
        rb[index] = b[index];
        vec.push_back(TrialD(lb, IOPPtr->ComputeFunction( lb )));
        vec.push_back(TrialD(rb, IOPPtr->ComputeFunction( rb )));
        double M = 0;
        size_t k = 2;
        size_t t = 0;
        M = abs((vec[1].z - vec[0].z) / (vec[1].x[index] - vec[0].x[index]));
        for (; ((stop_crit && abs(vec[t + 1].x[index] - IOPPtr->GetOptimumPoint()[index]) > eps)
            || (!stop_crit && (vec[t + 1].x[index] - vec[t].x[index] >= eps))); ++k) {
            find_glob_min_fixed_index(x, index + 1, stop_crit);
            for (size_t i = 0; i < (k - 1u); ++i) {
                double M_tmp = abs((vec[i + 1].z - vec[i].z) / (vec[i + 1].x[index] - vec[i].x[index]));
                if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R = m * (vec[1].x[index] - vec[0].x[index]) + (pow((vec[1].z - vec[0].z), 2) / (m * (vec[1].x[index] - vec[0].x[index]))) - 2 * (vec[1].z + vec[0].z);
            for (size_t i = 1; i < (k - 1u); ++i) {
                double R_tmp = m * (vec[i + 1].x[index] - vec[i].x[index]) +
                    (pow((vec[i + 1].z - vec[i].z), 2) / (m * (vec[i + 1].x[index] - vec[i].x[index]))) - 2 * (vec[i + 1].z + vec[i].z);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            x[index] = (vec[t].x[index] + vec[t + 1].x[index]) / 2 - (vec[t + 1].z - vec[t].z) / (2 * m);
            TrialD t1_pair(x, IOPPtr->ComputeFunction( x ));
            ++count;
            //vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const Trial& a, const Trial& b) {return a.x <= b.x; }), t1_pair); //No need for sorting, only to insert

            vec.insert(vec.begin() + t + 1, t1_pair);
            /*double M_tmp = abs((vec[t + 1].z - vec[t].z) / (vec[t + 1].x - vec[t].x));
            if (M_tmp > M) M = M_tmp;
            M_tmp = abs((vec[t + 2].z - vec[t+1].z) / (vec[t + 2].x - vec[t+1].x));
            if (M_tmp > M) M = M_tmp;*/
        }
        auto min = vec[t + 1];
        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i].z < min.z) {
                min = vec[i];
            }
        }
        sol = min; solved = true;
        return min;
    }

    TrialD find_glob_min(bool stop_crit = false) {
        vector<TrialD> vec;
        vec.push_back(TrialD(a, IOPPtr->ComputeFunction( a )));
        vec.push_back(TrialD(b, IOPPtr->ComputeFunction( b )));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;
        vector<double> x = a;
        M = abs((vec[1].z - vec[0].z) / (vec[1].x[0] - vec[0].x[0]));
        for (; ((stop_crit && abs(vec[t + 1].x[0] - IOPPtr->GetOptimumPoint()[0]) > eps)
            || (!stop_crit && (vec[t + 1].x[0] - vec[t].x[0] >= eps)))
            && k < NMax; ++k) {
            find_glob_min_fixed_index(x, 1, stop_crit);
            for (size_t i = 0; i < (k - 1u); ++i) {
                double M_tmp = abs((vec[i + 1].z - vec[i].z) / (vec[i + 1].x[0] - vec[i].x[0]));
                if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R = m * (vec[1].x[0] - vec[0].x[0]) + (pow((vec[1].z - vec[0].z), 2) / (m * (vec[1].x[0] - vec[0].x[0]))) - 2 * (vec[1].z + vec[0].z);
            for (size_t i = 1; i < (k - 1u); ++i) {
                double R_tmp = m * (vec[i + 1].x[0] - vec[i].x[0]) +
                    (pow((vec[i + 1].z - vec[i].z), 2) / (m * (vec[i + 1].x[0] - vec[i].x[0]))) - 2 * (vec[i + 1].z + vec[i].z);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            x[0] = (vec[t].x[0] + vec[t + 1].x[0]) / 2 - (vec[t + 1].z - vec[t].z) / (2 * m);
            TrialD t1_pair(x, IOPPtr->ComputeFunction(x));
            ++count;
            //vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const Trial& a, const Trial& b) {return a.x <= b.x; }), t1_pair); //No need for sorting, only to insert

            vec.insert(vec.begin() + t + 1, t1_pair);
            /*double M_tmp = abs((vec[t + 1].z - vec[t].z) / (vec[t + 1].x - vec[t].x));
            if (M_tmp > M) M = M_tmp;
            M_tmp = abs((vec[t + 2].z - vec[t+1].z) / (vec[t + 2].x - vec[t+1].x));
            if (M_tmp > M) M = M_tmp;*/
        }
        auto min = vec[t + 1];
        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i].z < min.z) {
                min = vec[i];
            }
        }
        sol = min; solved = true;
        return min;
    }

    void Show_info() {
        if (!solved) this->find_glob_min();
        //std::cout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
        std::cout << "�������� ���������� eps = " << eps << std::endl;
        std::cout << "�������� ������ r = " << r << std::endl;
        std::cout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            //fout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
            fout << "�������� ���������� eps = " << eps << std::endl;
            fout << "�������� ������ r = " << std::fixed << r << std::endl;
            fout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
        }
    }

    TrialD GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }

    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};