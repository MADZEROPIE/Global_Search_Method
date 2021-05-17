#pragma once
#include "IOptProblem.hpp"
#include "Map.h"

#include <iostream>
#include <functional>  //for std::function
#include <cmath>   //for math functions e.g. sin() or cos()
#include <vector> 
#include <algorithm>
#include <fstream>

#include <omp.h>

//#define NUMTH 6
//About style of naming. There is no style.

using std::vector;
using std::pair;

//typedef pair<double, double> Trial;
struct TrialD {
    vector<double> x;
    double z;
    TrialD(const vector<double>& _x = {}, double _z = 0.0) {
        x = _x;
        z = _z;
    }
};

struct TrialP {  // Trial for Peano. Same as Trial, so ...
    double x;
    double z;
    TrialP(double _x = 0, double _z = 0.0) {
        x = _x;
        z = _z;
    }
    TrialD toTrialD(const vector <double>& lb, const vector <double>& rb, int m=10) {
        TrialD res;
        res.x.resize(lb.size());
        res.z = z;
        vector<double> xd(lb.size());
        mapd(x, m, res.x.data(), lb.size(), 1);
        return res;
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
    vector<TrialD> all_trials;

public:
    MinimazerD(IOptProblem* _IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500) {
        IOPPtr = _IOPPtr;
        IOPPtr->GetBounds(a, b);
        r = (_r > 1.0) ? _r : 2.0;
        eps = _eps;
        NMax = _NMax;
        dim = IOPPtr->GetDimension();
    }

     void find_glob_min_fixed_index(std::vector<double>& x, int index, bool save_trials = false) {
        if (index >= dim) return;
        vector<TrialD> vec;
        auto lb = x;
        auto rb = x;
        lb[index] = a[index];
        rb[index] = b[index];
        find_glob_min_fixed_index(lb, index + 1, save_trials);
        find_glob_min_fixed_index(rb, index + 1, save_trials);
        vec.push_back(TrialD(lb, IOPPtr->ComputeFunction( lb )));
        vec.push_back(TrialD(rb, IOPPtr->ComputeFunction( rb )));
        double M = 0;
        size_t k = 2;
        size_t t = 0;
        M = abs((vec[1].z - vec[0].z) / (vec[1].x[index] - vec[0].x[index]));
        for (; ((vec[t + 1].x[index] - vec[t].x[index]) >= eps) && count < NMax; ++k) {
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
            find_glob_min_fixed_index(x, index + 1, save_trials);
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
                x = vec[i].x;
            }
            if (save_trials)
                all_trials.push_back(vec[i]);
        }
        sol = min;
    }

    TrialD find_glob_min(bool save_trials = false) {
        vector<TrialD> vec;
        auto lb = a, rb = b;
        find_glob_min_fixed_index(lb, 1 , save_trials);
        find_glob_min_fixed_index(rb, 1, save_trials);
        vec.push_back(TrialD(lb, IOPPtr->ComputeFunction( lb )));
        vec.push_back(TrialD(rb, IOPPtr->ComputeFunction( rb )));
        count = 2;
        double M = 0;
        int k = 2;
        int t = 0;
        vector<double> x = a;
        for (; (vec[t + 1].x[0] - vec[t].x[0]) >= eps && count < NMax; ++k) {  
            for (int i = 0; i < (k - 1); ++i) {
                double M_tmp = abs((vec[i + 1].z - vec[i].z) / (vec[i + 1].x[0] - vec[i].x[0]));
                if (M_tmp > M) M = M_tmp;
            }
            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R = m * (vec[1].x[0] - vec[0].x[0]) + (pow((vec[1].z - vec[0].z), 2) / (m * (vec[1].x[0] - vec[0].x[0]))) - 2 * (vec[1].z + vec[0].z);
            for (int i = 1; i < (k - 1); ++i) {
                double R_tmp = m * (vec[i + 1].x[0] - vec[i].x[0]) +
                    (pow((vec[i + 1].z - vec[i].z), 2) / (m * (vec[i + 1].x[0] - vec[i].x[0]))) - 2 * (vec[i + 1].z + vec[i].z);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            x[0] = (vec[t].x[0] + vec[t + 1].x[0]) / 2 - (vec[t + 1].z - vec[t].z) / (2 * m);
            find_glob_min_fixed_index(x, 1, save_trials);
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
            if (save_trials)
                all_trials.push_back(vec[i]);
        }
        sol = min; solved = true;
        return sol;
    }

    TrialD find_glob_min_Peano(bool save_trials = false) {  // SS
        vector<TrialP> vec;
        count = 2;
        double M = 0;
        int k = 2;
        int t = 0;
        double curr_diff = 1;
        vector<double> x1(dim);
        double x = 0;
        
        mapd(0, 10, x1.data(), dim, 1);
        //
        for (int i = 0; i < dim; ++i) {
            x1[i] = (b[i] - a[i]) * x1[i] + (a[i] + b[i]) / 2;  // b[i]+-a[i] are const, so...
            std::cout << x1[i] << " ";
        }
        std::cout << '\n';
        vec.push_back(TrialP(0, IOPPtr->ComputeFunction(x1)));

        mapd(1, 10, x1.data(), dim, 1);
        //
        for (int i = 0; i < dim; ++i) {
            x1[i] = (b[i] - a[i]) * x1[i] + (a[i] + b[i]) / 2;  // b[i]+-a[i] are const, so...
            std::cout << x1[i] << " ";
        }
        std::cout << '\n';
        vec.push_back(TrialP(1, IOPPtr->ComputeFunction(x1)));

        for (; curr_diff >= eps && k < NMax; ++k) {
            for (int i = 0; i < (k - 1); ++i) {
                double M_tmp = abs((vec[i + 1].z - vec[i].z) / pow((vec[i + 1].x - vec[i].x), 1.0/dim));
                if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = M;
            t = 0;
            double dist = pow(vec[1].x - vec[0].x, 1.0 / dim);
            double R = m * dist + (pow((vec[1].z - vec[0].z), 2) / (m * dist)) - 2 * (vec[1].z + vec[0].z);
            for (size_t i = 1; i < (k - 1u); ++i) {
                double dist = pow(vec[i + 1].x - vec[i].x, 1.0 / dim);
                double R_tmp = r*m * dist + (pow((vec[i + 1].z - vec[i].z), 2) / (r * m * dist)) - 2 * (vec[i + 1].z + vec[i].z);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            double x_t1 = (vec[t].x + vec[t + 1].x) / 2;
            if((vec[t + 1].z - vec[t].z) < 0)
                x_t1 += pow(abs(vec[t + 1].z - vec[t].z) /  m, dim) / (2*r);
            else 
                x_t1 -= pow(abs(vec[t + 1].z - vec[t].z) / m, dim) / (2*r);
            mapd(x, 10, x1.data(), dim, 1);
            //
            for (int i = 0; i < dim; ++i) {
                x1[i] = (b[i] - a[i]) * x1[i] + (a[i]+b[i]) / 2;  // b[i]+-a[i] are const, so...
            }
            ++count;
            curr_diff = std::min(pow(vec[t + 1].x - x_t1, 1.0/dim), pow((x_t1 - vec[t].x), 1.0 / dim));
            TrialP t1_pair(x_t1, IOPPtr->ComputeFunction(x1));
            vec.insert(vec.begin() + t + 1, t1_pair);
            
        }

        // Findin' min in vec
        auto min = vec[t + 1];
        for (int i = 0; i < vec.size(); ++i) {
            if (vec[i].z < min.z) {
                min = vec[i];
            }
            if (save_trials)
                all_trials.push_back(vec[i].toTrialD(a,b));
        }
        count = k;
        sol = min.toTrialD(a,b); solved = true;
        return sol;
    }

    void Show_info() {
        if (!solved) this->find_glob_min();
        //std::cout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
        std::cout << "Точность вычислений eps = " << eps << std::endl;
        std::cout << "Параметр метода r = " << r << std::endl;
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            //fout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
            fout << "Точность вычислений eps = " << eps << std::endl;
            fout << "Параметр метода r = " << std::fixed << r << std::endl;
            fout << "Функция была посчитана " << count << " раз(а)" << std::endl;
        }
    }

    TrialD GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }

    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
    vector<TrialD>& GetTrials() { return all_trials; }
    void saveTrialsInFile(std::ofstream& fout) {
        auto num_trials = all_trials.size();
        for (auto &tr : all_trials) {
            //auto& tr = all_trials[i];
            //std::cout << tr.x.size() << " ";
            for (int k = 0; k < dim - 1; ++k) {
                fout << tr.x[k] << ',';
            }
            fout << tr.x[dim - 1] << std::endl;
            
        }
        int a = 1;
        ++a;
    }

};
