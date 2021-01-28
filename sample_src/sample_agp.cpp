#include <iostream>
//#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "Hill/HillProblemFamily.hpp"
#include "Shekel/ShekelProblemFamily.hpp"


//Code from https://github.com/MADZEROPIE/Global-Search-Method
#include <functional> //for std::function
#include <cmath>  //for math functions e.g. sin() or cos()
#include <vector> 
#include <algorithm>
#include <fstream>

#include <thread>
#include <ctime>
#include <omp.h>

#define NUMTH 6
//About style of naming. There is no style.

using std::vector;
using std::pair;
using std::thread;

typedef pair<double, double> dpair;


class Minimazer { //Don't ask me why... But only because using a sledge-hammer to crack a nut sounds fun.

protected:
    //std::function<double(double)> func; //Function that needs findin' minimum
    double a; //Beginning of the segment
    double b; //End of the segment

    double r; //Coefficient of the method

    double eps;
    IOptProblem* IOPPtr; //Problem that needs findin' minimum
    dpair sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number

public:
    Minimazer(IOptProblem*  _IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax=500) {
        IOPPtr = _IOPPtr;
        vector<double> tmp_lb, tmp_rb;
        IOPPtr->GetBounds(tmp_lb, tmp_rb);
        a = tmp_lb[0]; b = tmp_rb[0];
        r = (_r > 1.0) ? _r : 2.0;
        eps = _eps;
        NMax = _NMax;
    }
   
    dpair find_glob_min(bool stop_crit = false) {
        vector<pair<double, double> > vec;
        vec.push_back(dpair(a, IOPPtr->ComputeFunction({ a })));
        vec.push_back(dpair(b, IOPPtr->ComputeFunction({ b })));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;

        for (; ((stop_crit && abs(vec[t + 1].first - IOPPtr->GetOptimumPoint()[0]) > eps) 
            || (!stop_crit && (vec[t+1].first-vec[t].first>=eps)))
            && k<NMax; ++k) {
            for (size_t i = 0; i < (k - 1u); ++i) {
                double M_tmp = abs((vec[i + 1].second - vec[i].second) / (vec[i + 1].first - vec[i].first));
                if (M_tmp > M) M = M_tmp;
            }

            double m = 1;
            if (M != 0) m = r * M;
            t = 0;
            double R = m * (vec[1].first - vec[0].first) + (pow((vec[1].second - vec[0].second), 2) / (m * (vec[1].first - vec[0].first))) - 2 * (vec[1].second + vec[0].second);
            for (size_t i = 1; i < (k - 1u); ++i) {
                double R_tmp = m * (vec[i + 1].first - vec[i].first) + (pow((vec[i + 1].second - vec[i].second), 2) / (m * (vec[i + 1].first - vec[i].first))) - 2 * (vec[i + 1].second + vec[i].second);
                if (R_tmp > R) { t = i; R = R_tmp; }
            }

            double x_t1 = (vec[t].first + vec[t + 1].first) / 2 - (vec[t + 1].second - vec[t].second) / (2 * m);
            dpair t1_pair(x_t1, IOPPtr->ComputeFunction({ (x_t1) }));
            ++count;
            vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const dpair& a, const dpair& b) {return a.first <= b.first; }), t1_pair); //No need for sorting, only to insert
        }
        auto min = vec[t+1];
        /*for (int i = 1; i < vec.size(); ++i) {
            if (vec[i].second < min.second) {
                min = vec[i];
            }
        }*/
        sol = min; solved = true;
        return min;
    }

    dpair find_glob_min_threadver() {
        vector<pair<double, double> > vec;
        //vec.push_back(dpair(a, IOPPtr->ComputeFunction({ a })));
        //vec.push_back(dpair(b, IOPPtr->ComputeFunction({ b })));
        count = NUMTH + 1;
        double M = 0;
        int k = NUMTH+1;
        int iter_count = 0;
        size_t t = 0;
        int tj_size = NUMTH;
        std::vector<std::pair<int, double> > tj_vec;
        double h = (b - a) / NUMTH;
        for (int i = 0; i <= NUMTH; ++i) {
            vec.push_back(dpair(a + i*h, IOPPtr->ComputeFunction({ a + i * h })));
        }
        for (; iter_count < NMax; k = vec.size()) {
            //int tj_size = (tj_vec.NUMTH() < NUMTH) ? tj_vec.NUMTH() : NUMTH;
            for (int i = 0; i < (k - 1); ++i) {
                if ((vec[i + 1].first - vec[i].first) < eps) {  // Можно искать только по интервалам tj, но ...
                    dpair min = vec[0];
                    for (int j = 1; j < k; ++j) {
                        if (vec[j].second < min.second)
                            min = vec[j];
                    }
                    count = k;
                    solved = true;
                    sol = min;
                    return min;
                }
            }
            for (int i = 0; i < (k - 1); ++i) {
                double M_tmp = abs((vec[i + 1].second - vec[i].second) / (vec[i + 1].first - vec[i].first));
                if (M_tmp > M)
                    M = M_tmp;
            }
            double m = 1.0;
            if (M != 0.0)
                m = r * M;
            tj_vec.resize(k - 1);
            double R;
            for (int i = 0; i < (k - 1); ++i) {
                R = m * (vec[i + 1].first - vec[i].first) + (pow((vec[i + 1].second - vec[i].second), 2) /
                    (m * (vec[i + 1].first - vec[i].first))) - 2 * (vec[i + 1].second + vec[i].second);
                tj_vec[i] = (dpair(i, R));
            }
            
            for (int j = 0; j < tj_size; ++j) {  // Ставим на первые tj_size мест максимальные R
                for (int l = j + 1; l < (k - 1); ++l) {
                    if (tj_vec[l].second > tj_vec[j].second) {
                        std::swap(tj_vec[l], tj_vec[j]);
                    }
                }
            }
            
            std::vector<dpair> tmp_vec(tj_size);
            //std::vector <std::thread> th_vec;
            for (int i = 0; i < tj_size; ++i) {
                tmp_vec[i].first = (vec[tj_vec[i].first + 1].first + vec[tj_vec[i].first].first) / 2 -
                    (vec[tj_vec[i].first + 1].second - vec[tj_vec[i].first].second) / (2 * m);
            }
            #pragma omp parallel shared(tmp_vec) num_threads(NUMTH)
            {
                #pragma omp for
                for (int i = 0; i < tj_size; ++i) {
                    tmp_vec[i].second = IOPPtr->ComputeFunction({ tmp_vec[i].first });
                }
            }
            //for (int i = 0; i < tj_size; ++i) {
            //    th_vec.push_back((thread([](dpair* pair, IOptProblem* IOPPtr, double x_t) {
            //        // //   std::cout << "1";
            //        *pair = dpair(x_t, IOPPtr->ComputeFunction({ x_t }));
            //        }, &tmp_vec[i], IOPPtr, tmp_vec[i].first)));
            //    //tmp_vec[i] = dpair(x_t, IOPPtr->ComputeFunction({ x_t }));
            //}
            //for (auto& th:th_vec) {
            //    th.join();
            //}
            for (auto& t_pair : tmp_vec) {
                vec.insert(std::lower_bound(vec.begin(), vec.end(), t_pair,
                    [](const dpair& a, const dpair& b) {
                        return a.first <= b.first;
                    }), t_pair);  // No need for sorting, only to insert
            }
            ++iter_count;
        }

        dpair min = vec[0];
        for (int j = 1; j < k; ++j) {
            if (vec[j].second < min.second)
                min = vec[j];
        }
        count = k;
        solved = true;
        sol = min;
        return min;
    }

    void Show_info() {
        if (!solved) this->find_glob_min();
        std::cout << "Минимум функции = " << sol.second << " в точке x = " << sol.first << std::endl;
        std::cout << "Точность вычислений eps = " << eps << std::endl;
        std::cout << "Параметр метода r = " << r << std::endl;
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "Минимум функции = " << sol.second << " в точке x = " << sol.first << std::endl;
            fout << "Точность вычислений eps = " << eps << std::endl;
            fout << "Параметр метода r = " << std::fixed <<r << std::endl;
            fout << "Функция была посчитана " << count << " раз(а)" << std::endl;
        }
    }
    dpair GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};
//End of code from https://github.com/MADZEROPIE/Global-Search-Method


class Minimazer2 { // TODO: Ancestor (parent) class
protected:
    std::function<double(double)> func; // Function that needs findin' minimum
    double a; // Beginning of the segment
    double b; // End of the segment
    std::vector<std::function<double(double)> > g_vec; // g functions vector
    std::vector<double> r; // Coefficient of the method. VECTOR OF COEFFs ?
    std::vector<double> e_vec; // Reserve vector (?)
    double eps;
    //IOptProblem* IOPPtr; // Problem that needs findin' minimum
    dpair sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number

public:
    Minimazer2(const std::function<double(double)>& _func, const std::vector<std::function<double(double)> >& _g_vec,
        const std::vector<double>& _e_vec, std::vector<double> _r,
        double _a = 0.0, double _b = 1.0, double _eps = 0.01, uint64_t _NMax = 500) :
        func(_func), g_vec(_g_vec), e_vec(_e_vec), r(_r)
    {
        a = _a; b = _b;
        eps = _eps;
        NMax = _NMax;
    }

    dpair find_glob_min(bool stop_crit = false) {
        vector<dpair> vec;
        vec.push_back(dpair(a, func(a)));
        vec.push_back(dpair(b, func(b)));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;


        solved = true;
        return sol;
    }

    void Show_info() {
        //if (!solved) this->find_glob_min();
        std::cout << "Минимум функции = " << sol.second << " в точке x = " << sol.first << std::endl;
        std::cout << "Точность вычислений eps = " << eps << std::endl;
        /*std::cout << "Параметр метода r = " << r << std::endl;*/
        std::cout << "Функция была посчитана " << count << " раз(а)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        //if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "Минимум функции = " << sol.second << " в точке x = " << sol.first << std::endl;
            fout << "Точность вычислений eps = " << eps << std::endl;
            /*fout << "Параметр метода r = " << std::fixed << r << std::endl;*/
            fout << "Функция была посчитана " << count << " раз(а)" << std::endl;
        }
    }
    dpair GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};


class Tester {
private:
    Minimazer Min;
    dpair expected;
    dpair deviation;
    double eps;

public:
    Tester(IOptProblem* IOPPtr, double _eps = 0.01,double _r = 2.0, uint64_t _NMax = 500): Min(IOPPtr, _eps,_r,_NMax) {
        eps = _eps;
        auto exp_tmp = IOPPtr->GetOptimumPoint();
        if (abs(IOPPtr->GetOptimumValue() - IOPPtr->ComputeFunction(exp_tmp)) > 0.001) {
            std::cout.precision(10);
            std::cout << IOPPtr->GetOptimumValue() << " " << IOPPtr->ComputeFunction(exp_tmp) << std::endl;
            std::cout << "INCORRECT OPTIMUM VALUE OR POINT" << std::endl;
        }
        expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
    }

    bool Test(bool stop_crit) {
        dpair res = Min.find_glob_min(stop_crit);
        
        double dev = (res.first - expected.first);
        //std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < 2.0 * eps);    
    }

    bool Test_par(size_t size = 5 ) {
        dpair res = Min.find_glob_min_threadver();

        double dev = (res.first - expected.first);
        //std::cout << std::this_thread::get_id()<< ((abs(dev) < eps*2)? " YEEEEEEEEEEEEEEEEEES": " NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < (eps * 2));
    }

    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res = Min.GetMin();
        std::cout << "Ожидаемый результат: y = " << expected.second << " в точке x = " << expected.first << std::endl;
        std::cout << "Отклонение по x = " << res.first - expected.first << std::endl;
        std::cout << "Отклонение по y = " << res.second -expected.second << std::endl;
        std::cout << std::endl;
    }

    void Show_info_in_file(std::ofstream& fout) {
        //if (!Min.IsSolved()) Test();
        if (fout.is_open()) {
            Min.Show_info_in_file(fout);
            auto res = Min.GetMin();
            fout << "Ожидаемый результат: y = " << expected.second << " в точке x = " << expected.first << std::endl;
            fout << "Отклонение по x = " << res.first - expected.first << std::endl;
            fout << "Отклонение по y = " << res.second - expected.second << std::endl;
            fout << std::endl;
        }
    }
    unsigned long long GetCount() { return Min.GetCount(); }


};


void func(IOptProblemFamily* IOPFPtr, std::string filepath , double r, double eps, uint64_t NMax, const std::string& family_name, bool stop_crit) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    vector<int> CountVec1(NMax+1);
    for (int i = 0; i <= NMax; ++i) {
        CountVec1[i] = 0;
    }
    for (size_t i = 0; i < IOPFPtr->GetFamilySize(); ++i) {
        //std::cout << "Тестируется " << family_name << " Problem" << i << std::endl;
        Tester Tes(IOPFPtr->operator[](i), eps, r, NMax);
        //bool tmp = Tes.Test_par(10);
        //Tes.Show_info();
        bool tmp = Tes.Test(stop_crit);
        if (tmp) {
            ++CorrectCount;
            //std::cout << "YEP\n";

            ++CountVec1[Tes.GetCount()];
        }
        else {
            std::cout << i<<" NOPE\n";
            Tes.Show_info();
        }
    }
    std::cout << "Правильно решено " << CorrectCount << " из " << IOPFPtr->GetFamilySize() <<" " << family_name 
        << " family." << std::endl << std::endl;
    //file << "Правильно решено " << CorrectCount << " из " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;

    file << "sep=,\n";
    file << "0,0\n";
    double tmp;
    int i;
    for (i = 1; i < CountVec1.size(); ++i) {
        CountVec1[i] += CountVec1[i - 1];
    }
    for (i = 0; i < CountVec1.size(); ++i) {
        tmp = double(CountVec1[i]) / double(IOPFPtr->GetFamilySize());
        file << i << ',' << tmp << "\n";
        for (; (i + 1) < CountVec1.size() && CountVec1[i + 1] == CountVec1[i]; ++i);
    }
    file << '\n';
}



int main(int argc,char* argv[]) {
    
    std::string filepath = "results";
    
    if (argc > 1){ filepath = argv[1];}
    double r = 5;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.001;
    if (argc > 3) eps = std::stod(argv[3]);
    bool stop_crit = false;
    if (argc > 4) {
        stop_crit = (argv[4] == "ON");
    }

    setlocale(LC_ALL, "Russian");

    //file.precision(6);

    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t NMax = 500000;
    THansenProblemFamily HFam;
    THillProblemFamily HillFam;
    TShekelProblemFamily ShekFam;

    vector<IOptProblemFamily*> vec = { &HFam, & HillFam, & ShekFam };

    vector<std::string> names_vec = { "Hansen" ,"Hill","Shekel" };
    auto t1 = omp_get_wtime();
    for (size_t i = 0; i < vec.size();++i) {
        func(vec[i], filepath + names_vec[i] + ".csv", r, eps, NMax, names_vec[i], stop_crit);
    }
    auto t2 = omp_get_wtime();
    std::cout <<"Time: "<< t2 - t1 << std::endl;
    return 0;
}
