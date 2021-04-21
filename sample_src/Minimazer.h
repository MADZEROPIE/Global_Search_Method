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
struct Trial {
    double x, z;
    Trial(double _x = 0.0, double _z = 0.0) {
        x = _x;
        z = _z;
    }
};

class Minimazer { //Don't ask me why... But only because using a sledge-hammer to crack a nut sounds fun.

protected:
    //std::function<double(double)> func; //Function that needs findin' minimum
    double a; //Beginning of the segment
    double b; //End of the segment

    double r; //Coefficient of the method

    double eps;
    IOptProblem* IOPPtr; //Problem that needs findin' minimum
    Trial sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number

public:
    Minimazer(IOptProblem* _IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500) {
        IOPPtr = _IOPPtr;
        vector<double> tmp_lb, tmp_rb;
        IOPPtr->GetBounds(tmp_lb, tmp_rb);
        a = tmp_lb[0]; b = tmp_rb[0];
        r = (_r > 1.0) ? _r : 2.0;
        eps = _eps;
        NMax = _NMax;
    }

    Trial find_glob_min(bool stop_crit = false) {
        vector<Trial> vec;
        vec.push_back(Trial(a, IOPPtr->ComputeFunction({ a })));
        vec.push_back(Trial(b, IOPPtr->ComputeFunction({ b })));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;
        M = abs((vec[1].z - vec[0].z) / (vec[1].x - vec[0].x));
        for (; ((stop_crit && abs(vec[t + 1].x - IOPPtr->GetOptimumPoint()[0]) > eps)
            || (!stop_crit && (vec[t + 1].x - vec[t].x >= eps)))
            && k < NMax; ++k) {
            for (size_t i = 0; i < (k - 1u); ++i) {
                double M_tmp = abs((vec[i + 1].z - vec[i].z) / (vec[i + 1].x - vec[i].x));
                if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R = m * (vec[1].x - vec[0].x) + (pow((vec[1].z - vec[0].z), 2) / (m * (vec[1].x - vec[0].x))) - 2 * (vec[1].z + vec[0].z);
            for (size_t i = 1; i < (k - 1u); ++i) {
                double R_tmp = m * (vec[i + 1].x - vec[i].x) + (pow((vec[i + 1].z - vec[i].z), 2) / (m * (vec[i + 1].x - vec[i].x))) - 2 * (vec[i + 1].z + vec[i].z);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            double x_t1 = (vec[t].x + vec[t + 1].x) / 2 - (vec[t + 1].z - vec[t].z) / (2 * m);
            Trial t1_pair(x_t1, IOPPtr->ComputeFunction({ (x_t1) }));
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

    Trial find_glob_min_threadver() {
        vector<Trial> vec;
        count = NUMTH + 1;
        double M = 0;
        int k = NUMTH + 1;
        int iter_count = 0;
        size_t t = 0;
        int tj_size = NUMTH;
        std::vector<Trial> tj_vec;
        double h = (b - a) / NUMTH;
        omp_set_num_threads(NUMTH);
        for (int i = 0; i <= NUMTH; ++i) {
            vec.push_back(Trial(a + i * h, IOPPtr->ComputeFunction({ a + i * h })));
        }
        for (; iter_count < NMax; iter_count = k = vec.size()) {
            //int tj_size = (tj_vec.NUMTH() < NUMTH) ? tj_vec.NUMTH() : NUMTH;
            for (int i = 0; i < (k - 1); ++i) {
                if ((vec[i + 1].x - vec[i].x) < eps) {  // Можно искать только по интервалам tj, но ...
                    Trial min = vec[0];
                    for (int j = 1; j < k; ++j) {
                        if (vec[j].z < min.z)
                            min = vec[j];
                    }
                    count = k;
                    solved = true;
                    sol = min;
                    return min;
                }
            }
            //#pragma omp parallel shared(vec) num_threads(NUMTH)
            {
                //#pragma omp for
                for (int i = 0; i < (k - 1); ++i) {
                    double M_tmp = abs((vec[i + 1].z - vec[i].z) / (vec[i + 1].x - vec[i].x));
                    if (M_tmp > M)
                        M = M_tmp;
                }
            }
            double m = 1.0;
            if (M != 0.0)
                m = r * M;
            tj_vec.resize(k - 1);
            double R;
            for (int i = 0; i < (k - 1); ++i) {
                R = m * (vec[i + 1].x - vec[i].x) + (pow((vec[i + 1].z - vec[i].z), 2) /
                    (m * (vec[i + 1].x - vec[i].x))) - 2 * (vec[i + 1].z + vec[i].z);
                tj_vec[i] = (Trial(i, R));
            }

            for (int j = 0; j < tj_size; ++j) {  // Ставим на первые tj_size мест максимальные R
                for (int l = j + 1; l < (k - 1); ++l) {
                    if (tj_vec[l].z > tj_vec[j].z) {
                        std::swap(tj_vec[l], tj_vec[j]);
                    }
                }
            }

            std::vector<Trial> tmp_vec(tj_size);
            //std::vector <std::thread> th_vec;
            for (int i = 0; i < tj_size; ++i) {
                tmp_vec[i].x = (vec[tj_vec[i].x + 1].x + vec[tj_vec[i].x].x) / 2 -
                    (vec[tj_vec[i].x + 1].z - vec[tj_vec[i].x].z) / (2 * m);
            }
#pragma omp parallel shared(tmp_vec) num_threads(NUMTH)
            {
#pragma omp for
                for (int i = 0; i < tj_size; ++i) {
                    tmp_vec[i].z = IOPPtr->ComputeFunction({ tmp_vec[i].x });
                }
            }

            for (auto& t_pair : tmp_vec) {
                vec.insert(std::lower_bound(vec.begin(), vec.end(), t_pair,
                    [](const Trial& a, const Trial& b) {
                        return a.x <= b.x;
                    }), t_pair);  // No need for sorting, only to insert
            }
            ++iter_count;
        }

        Trial min = vec[0];
        for (int j = 1; j < k; ++j) {
            if (vec[j].z < min.z)
                min = vec[j];
        }
        count = k;
        solved = true;
        sol = min;
        return min;
    }

    void Show_info() {
        if (!solved) this->find_glob_min();
        std::cout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
        std::cout << "Точность вычислений eps = " << eps << std::endl;
        std::cout << "Параметр метода r = " << r << std::endl;
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
            fout << "Точность вычислений eps = " << eps << std::endl;
            fout << "Параметр метода r = " << std::fixed << r << std::endl;
            fout << "Функция была посчитана " << count << " раз(а)" << std::endl;
        }
    }
    Trial GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};