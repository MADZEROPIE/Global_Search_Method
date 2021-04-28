#pragma once
#include "MyConstrainedProblem.h"

#include <functional> //for std::function
#include <cmath>  //for math functions e.g. sin() or cos()
#include <vector> 
#include <algorithm>
#include <fstream>
#include <iostream>


#include <omp.h>

class MyMinimazer { // TODO: Ancestor (parent) class. Or no.
protected:
    //std::function<double(double)> func; // Function that needs findin' minimum
    double a; // Beginning of the segment
    double b; // End of the segment
    //std::vector<std::function<double(double)> > g_vec; // g functions vector

    MyConstrainedProblem* MCPtr; //Problem that needs findin' minimum
    uint m;
    std::vector<double> r; // Coefficient of the method. VECTOR OF COEFFs ?
    std::vector<double> e_vec; // Reserve vector (?)
    double eps;
    ConsTrial sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number
    double q = 0.01;

    ConsTrial make_trial(double x) {
        ConsTrial tr;
        tr.x = x;
        for (uint i = 0; i < m; ++i) {
            tr.z = MCPtr->ComputeGjConstr(i, x);
            if (tr.z > 0.0) {
                tr.index = i;
                return tr;
            }
        }
        tr.z = MCPtr->ComputeGjConstr(m, x);
        tr.index = m;
        return tr;
    }
public:
    MyMinimazer(MyConstrainedProblem* _MCPtr, std::vector<double> _r, double _eps = 0.01, uint64_t _NMax = 500) {
        MCPtr = _MCPtr;
        m = MCPtr->GetNumberofConstr();
        MCPtr->GetBounds(a, b);
        r = _r;
        eps = _eps;
        NMax = _NMax;
    }

    ConsTrial find_glob_min(bool stop_crit = false) {
        vector<ConsTrial> vec;
        vector<vector<ConsTrial>> Vvec(m + 1);
        vector<double> mv(m + 1), zv(m + 1), ev(m + 1);
        vec.push_back(this->make_trial(a));
        vec.push_back(this->make_trial(b));
        Vvec[vec[0].index].push_back(vec[0]);
        Vvec[vec[1].index].push_back(vec[1]);
        count = 2;
        for (int i = 0; i <= m; ++i) {
            mv[i] = 0;
        }
        double max_eps = b - a;
        int t = 0;
        while (max_eps >= eps) {
            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                mv[i] = 0;
                int s = Vvec[i].size();
                for (int j = 0; j < s - 1; ++j) {
                    double M_tmp = abs(Vvec[i][j + 1].z - Vvec[i][j].z) / (Vvec[i][j + 1].x - Vvec[i][j].x);
                    if (M_tmp > mv[i])
                        mv[i] = M_tmp;

                }
                if (mv[i] == 0) mv[i] = 1; //???
                //else mv[i] = r[i] * mv[i];
                ev[i] = mv[i] * q * eps;
            }


            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                if (Vvec[i].size() != 0) {
                    double z_tmp = Vvec[i][0].z;
                    int s = Vvec[i].size();
                    for (int j = 1; j < s; ++j)
                        if (Vvec[i][j].z < z_tmp)
                            z_tmp = Vvec[i][j].z;

                    if (z_tmp <= 0)
                        zv[i] = -ev[i];
                    else
                        zv[i] = z_tmp;
                }
            }

            // Finally R. But hey, SR is better
            int s = vec.size();
            double R = -INFINITY;

            for (int i = 0; i < s - 1; ++i) {
                double R_tmp;
                if (vec[i].index == vec[i + 1].index) {
                    auto v = vec[i].index;
                    R_tmp = vec[i + 1].x - vec[i].x + (vec[i + 1].z - vec[i].z) * (vec[i + 1].z - vec[i].z) /  // I know, I know...
                        (r[v] * r[v] * mv[v] * mv[v] * (vec[i + 1].x - vec[i].x)) - 2 * (vec[i + 1].z + vec[i].z - 2 * zv[v]) / (r[v] * mv[v]);
                }
                else if (vec[i].index < vec[i + 1].index) {
                    auto v = vec[i + 1].index;
                    R_tmp = 2 * (vec[i + 1].x - vec[i].x) - 4 * (vec[i + 1].z - zv[v]) / (r[v] * mv[v]);
                }
                else {
                    auto v = vec[i].index;
                    R_tmp = 2 * (vec[i + 1].x - vec[i].x) - 4 * (vec[i].z - zv[v]) / (r[v] * mv[v]);
                }
                if (R_tmp > R) {
                    t = i;
                    R = R_tmp;
                }
            }
            double x_new;
            if (vec[t].index != vec[t + 1].index) {
                x_new = (vec[t].x + vec[t + 1].x) / 2;
            }
            else {
                x_new = (vec[t].x + vec[t + 1].x) / 2 - (vec[t + 1].z + vec[t].z) / (2 * r[vec[t].index] * mv[vec[t].index]);
            }

            auto tr = make_trial(x_new);
            max_eps = std::min(vec[t + 1].x - x_new, x_new - vec[t].x);
            vec.insert(vec.begin() + t + 1, tr);
            Vvec[tr.index].push_back(tr);  // TODO: INSERT WITH BINARY SEARCH
            std::sort(Vvec[tr.index].begin(), Vvec[tr.index].end(), [](const ConsTrial& a, const ConsTrial& b) { return a.x < b.x; });
            ++count;
            //std::cout << tr.x << " " << tr.z << '\n';

        }
        if (Vvec[m].size() > 0) {
            sol = Vvec[m][0];
            for (int i = 1; i < Vvec[m].size(); ++i) {
                if (sol.z > Vvec[m][i].z)
                    sol = Vvec[m][i];
            }
            solved = true;
        }
        //sol = vec[t + 1];
        return sol;
    }

    ConsTrial find_glob_min_BF(uint64_t steps = 100000) {
        double h = (b - a) / steps;
        double x = a;
        bool flag = false;
        for (; x <= b && !flag; x += h) {
            auto tr = make_trial(x);
            if (tr.index == m) {
                sol = tr;
                flag = true;
                solved = true;
            }
        }
        for (; x <= b; x += h) {
            auto tr = make_trial(x);
            if (tr.index == m && tr.z < sol.z) {
                sol = tr; 
            }
        }
        count = steps + 1;
        return sol;
    }

    void Show_info() {
        //if (!solved) this->find_glob_min();
        std::cout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
        std::cout << "Точность вычислений eps = " << eps << std::endl;
        /*std::cout << "Параметр метода r = " << r << std::endl;*/
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        //if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "Минимум функции = " << sol.z << " в точке x = " << sol.x << std::endl;
            fout << "Точность вычислений eps = " << eps << std::endl;
            /*fout << "Параметр метода r = " << std::fixed << r << std::endl;*/
            fout << "Функция была посчитана " << count << " раз(а)" << std::endl;
        }
    }
    ConsTrial GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};
