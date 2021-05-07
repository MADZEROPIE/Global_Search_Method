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
    double q = 10;

    ConsTrial make_trial(double x) {
        ConsTrial tr;
        tr.x = x;
        tr.z.resize(m+1);
        for (uint i = 0; i < m; ++i) {
            tr.z[i]=(MCPtr->ComputeGjConstr(i, x));
            if (tr.z[i] > 0.0) {
                tr.index = i;
                return tr;
            }
        }
        tr.z[m] = (MCPtr->ComputeFunction(x));
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
        vector<vector<ConsTrial>> Vvec(m + 1);
        vector<ConsTrial>& vec = Vvec[0];
        vector<double> mv(m + 1), zv(m + 1), ev(m + 1);
        vec.push_back(this->make_trial(a));
        vec.push_back(this->make_trial(b));
        count = 2;
        for (int i = 1; i <= m; ++i) {
            if (vec[0].index >= i)
                Vvec[i].push_back(vec[0]);
            if (vec[1].index >= i)
                Vvec[i].push_back(vec[1]);
        }

        double max_eps = b - a;
        int t = 0;
        while (max_eps >= eps) {
            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                mv[i] = 0;
                int s = Vvec[i].size();
                for (int j = 0; j < s - 1; ++j) {
                    double M_tmp = abs(Vvec[i][j + 1].z[i] - Vvec[i][j].z[i]) / (Vvec[i][j + 1].x - Vvec[i][j].x);
                    if (M_tmp > mv[i])
                        mv[i] = M_tmp;
                }
                if (mv[i] == 0) mv[i] = 1;
                //else mv[i] = r[i] * mv[i];
                ev[i] = mv[i] * q * eps;
            }

            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                if (Vvec[i].size() != 0) {
                    zv[i] = INFINITY;
                    double z_tmp = Vvec[i][0].z[i];
                    int s = Vvec[i].size();
                    for (int j = 1; j < s; ++j) {
                        if (Vvec[i][j].z[i] < z_tmp)
                            z_tmp = Vvec[i][j].z[i];
                    }

                    if (z_tmp > 0 || i == m) {
                        zv[i] = z_tmp;
                    }
                    else {
                        std::cout << i << " ";
                        zv[i] = -ev[i];
                    }
                }
            }

            // Finally R. But hey, SR is better
            int s = vec.size();
            double R = -INFINITY;

            for (int i = 0; i < s - 1; ++i) {
                double R_tmp;
                double del = vec[i + 1].x - vec[i].x;
                if (vec[i].index == vec[i + 1].index) {
                    auto v = vec[i].index;
                    std::cout << v << " " << zv[v] << std::endl;
                    R_tmp = del + (vec[i + 1].z[v] - vec[i].z[v]) * (vec[i + 1].z[v] - vec[i].z[v]) /
                        (r[v] * r[v] * mv[v] * mv[v] * del) - 2 * (vec[i + 1].z[v] + vec[i].z[v] - 2 * zv[v]) / (r[v] * mv[v]);
                }
                else if (vec[i].index < vec[i + 1].index) {
                    auto v = vec[i + 1].index;
                    std::cout << v << " " << zv[v] << std::endl;
                    R_tmp = 2 * del - 4 * (vec[i + 1].z[v] - zv[v]) / (r[v] * mv[v]); // (vec[i + 1].z[v]
                }
                else {
                    auto v = vec[i].index;
                    std::cout << v << " " << zv[v] << std::endl;
                    R_tmp = 2 * del - 4 * (vec[i].z[v] - zv[v]) / (r[v] * mv[v]);
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
                int v = vec[t].index;
                x_new = (vec[t].x + vec[t + 1].x) / 2 - (vec[t + 1].z[v] - vec[t].z[v]) / (2 * r[v] * mv[v]);
            }

            auto tr = make_trial(x_new);
            max_eps = std::min(vec[t + 1].x - x_new, x_new - vec[t].x);
            vec.insert(vec.begin() + t + 1, tr);
            for (int i = 1; i <= m; ++i) {
                if (tr.index >= i) {
                    Vvec[i].push_back(tr);  // TODO: INSERT WITH BINARY SEARCH
                    std::sort(Vvec[i].begin(), Vvec[i].end(), [](const ConsTrial& a, const ConsTrial& b) { return a.x < b.x; });
                }
                else break;
            }
            ++count;
            //std::cout << tr.x << " " << tr.z << '\n';

        }
        if (Vvec[m].size() > 0) {
            sol = Vvec[m][0];
            for (int i = 1; i < Vvec[m].size(); ++i) {
                if (sol.z[m] > Vvec[m][i].z[m])
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
            if (tr.index == m && tr.z[m] < sol.z[m]) {
                sol = tr; 
            }
        }
        count = steps + 1;
        return sol;
    }

    void Show_info() {
        if (solved) {
            std::cout << "������� ������� = " << sol.z[m] << " � ����� x = " << sol.x << std::endl;
            std::cout << "�������� ���������� eps = " << eps << std::endl;
            /*std::cout << "�������� ������ r = " << r << std::endl;*/
            std::cout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
        }
    }
    void Show_info_in_file(std::ofstream& fout) {
        //if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "������� ������� = " << sol.z[m] << " � ����� x = " << sol.x << std::endl;
            fout << "�������� ���������� eps = " << eps << std::endl;
            /*fout << "�������� ������ r = " << std::fixed << r << std::endl;*/
            fout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
        }
    }
    ConsTrial GetMin() {
        //if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};
