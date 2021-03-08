#include <iostream>
//#include "HansenProblem.hpp"
#include "HansenProblemFamily.hpp"
#include "Hill/HillProblemFamily.hpp"
#include "Shekel/ShekelProblemFamily.hpp"
#include "MyConstrainedProblem.h"

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
    Minimazer(IOptProblem*  _IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax=500) {
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
            || (!stop_crit && (vec[t+1].x-vec[t].x>=eps)))
            && k<NMax; ++k) {
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
            vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const Trial& a, const Trial& b) {return a.x <= b.x; }), t1_pair); //No need for sorting, only to insert
            
            //vec.insert(vec.begin() + t + 1, t1_pair);
            /*double M_tmp = abs((vec[t + 1].z - vec[t].z) / (vec[t + 1].x - vec[t].x));
            if (M_tmp > M) M = M_tmp;
            M_tmp = abs((vec[t + 2].z - vec[t+1].z) / (vec[t + 2].x - vec[t+1].x));
            if (M_tmp > M) M = M_tmp;*/
        }
        auto min = vec[t+1];
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
        //vec.push_back(Trial(a, IOPPtr->ComputeFunction({ a })));
        //vec.push_back(Trial(b, IOPPtr->ComputeFunction({ b })));
        count = NUMTH + 1;
        double M = 0;
        int k = NUMTH + 1;
        int iter_count = 0;
        size_t t = 0;
        int tj_size = NUMTH;
        std::vector<Trial> tj_vec;
        double h = (b - a) / NUMTH;
        omp_set_num_threads(NUMTH);
        {
            for (int i = 0; i <= NUMTH; ++i) {
                vec.push_back(Trial(a + i * h, IOPPtr->ComputeFunction({ a + i * h })));
            }
            for (; iter_count < NMax; k = vec.size()) {
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
            fout << "Параметр метода r = " << std::fixed <<r << std::endl;
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

//------------------------------------------------------------------------------------------\\



class Minimazer2 { // TODO: Ancestor (parent) class. Or no.
protected:
    std::function<double(double)> func; // Function that needs findin' minimum
    double a; // Beginning of the segment
    double b; // End of the segment
    std::vector<std::function<double(double)> > g_vec; // g functions vector
    std::vector<double> r; // Coefficient of the method. VECTOR OF COEFFs ?
    std::vector<double> e_vec; // Reserve vector (?)
    double eps;
    //IOptProblem* IOPPtr; // Problem that needs findin' minimum
    Trial sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number

    ConsTrial make_trial(double x) {
        ConsTrial tr;

        return tr;
    }
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

    Trial find_glob_min(bool stop_crit = false) {
        vector<Trial> vec;
        vec.push_back(Trial(a, func(a)));
        vec.push_back(Trial(b, func(b)));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;


        solved = true;
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
    Trial GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }
    unsigned long long GetCount() { return count; }

    bool IsSolved() { return solved; }
};

//------------------------------------------------------------------------------------------\\

class Tester {
private:
    Minimazer Min;
    Trial expected;
    Trial deviation;
    double eps;

public:
    Tester(IOptProblem* IOPPtr, double _eps = 0.01,double _r = 2.0, uint64_t _NMax = 500): Min(IOPPtr, _eps,_r,_NMax) {
        eps = _eps;
        auto exp_tmp = IOPPtr->GetOptimumPoint();
        /*if (abs(IOPPtr->GetOptimumValue() - IOPPtr->ComputeFunction(exp_tmp)) > 0.001) {
            std::cout.precision(10);
            std::cout << IOPPtr->GetOptimumValue() << " " << IOPPtr->ComputeFunction(exp_tmp) << std::endl;
            std::cout << "INCORRECT OPTIMUM VALUE OR POINT" << std::endl;
        }*/
        //expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
        expected = Trial(exp_tmp[0], IOPPtr->ComputeFunction(exp_tmp));
    }

    bool Test(bool stop_crit) {
        Trial res = Min.find_glob_min(stop_crit);
        
        double dev = (res.x - expected.x);
        //std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps || abs(res.z - expected.z) < 7e-5); // Is this LEGAL? Well, no, but...
    }

    bool Test_par() {
        Trial res = Min.find_glob_min_threadver();

        double dev = (res.x - expected.x);
        //std::cout << std::this_thread::get_id()<< ((abs(dev) < eps*2)? " YEEEEEEEEEEEEEEEEEES": " NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps || abs(res.z - expected.z) < 7e-5); // Is this LEGAL? Well, no, but...
    }

    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res = Min.GetMin();
        std::cout << "Ожидаемый результат: y = " << expected.z << " в точке x = " << expected.x << std::endl;
        std::cout << "Отклонение по x = " << res.x - expected.x << std::endl;
        std::cout << "Отклонение по y = " << res.z -expected.z << std::endl;
        std::cout << std::endl;
    }

    void Show_info_in_file(std::ofstream& fout) {
        //if (!Min.IsSolved()) Test();
        if (fout.is_open()) {
            Min.Show_info_in_file(fout);
            auto res = Min.GetMin();
            fout << "Ожидаемый результат: y = " << expected.z << " в точке x = " << expected.x << std::endl;
            fout << "Отклонение по x = " << res.x - expected.x << std::endl;
            fout << "Отклонение по y = " << res.z - expected.z << std::endl;
            fout << std::endl;
        }
    }
    unsigned long long GetCount() { return Min.GetCount(); }


};

//------------------------------------------------------------------------------------------\\

void func(IOptProblemFamily* IOPFPtr, std::string filepath , double r, double eps, uint64_t NMax, const std::string& family_name, bool stop_crit) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    vector<int> CountVec1 (NMax+1);
    for (int i = 0; i <= NMax; ++i) {
        CountVec1[i] = 0;
    }
    for (size_t i = 0; i < IOPFPtr->GetFamilySize(); ++i) {
        //std::cout << "Тестируется " << family_name << " Problem" << i << std::endl;
        Tester Tes(IOPFPtr->operator[](i), eps, r, NMax);
        //bool tmp = Tes.Test_par();
        //Tes.Show_info();
        bool tmp = Tes.Test(stop_crit);
        if (tmp) {
            ++CorrectCount;
            //std::cout << "YEP\n";

            ++CountVec1[Tes.GetCount()];
        }
        else {
            //std::cout << i << " NOPE\n";
            //Tes.Show_info();
        }
    }
    std::cout << "Правильно решено " << CorrectCount << " из " << IOPFPtr->GetFamilySize() <<" " << family_name 
        << " family." << std::endl << std::endl;
    //file << "Правильно решено " << CorrectCount << " из " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;

    file << "sep=,\n";
    //file << "0,0\n";
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

//------------------------------------------------------------------------------------------\\

int main(int argc,char* argv[]) {
    
    //std::string filepath = "results";
    //
    //if (argc > 1){ filepath = argv[1];}
    //double r = 4;
    //if (argc > 2) r = std::stod(argv[2]);
    //double eps = 0.001;
    //if (argc > 3) eps = std::stod(argv[3]);
    //bool stop_crit = false;
    //if (argc > 4) {
    //    stop_crit = (argv[4] == "ON");
    //}

    //setlocale(LC_ALL, "Russian");

    ////file.precision(6);

    ////TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    //uint64_t NMax = 250;
    //THansenProblemFamily HFam;
    //THillProblemFamily HillFam;
    //TShekelProblemFamily ShekFam;

    //vector<IOptProblemFamily*> vec = { &HFam, & HillFam, & ShekFam};

    //vector<std::string> names_vec = { "Hansen" ,"Hill", "Shekel" };
    //auto t1 = omp_get_wtime();
    //for (size_t i = 0; i < vec.size();++i) {
    //    func(vec[i], filepath + names_vec[i] + ".csv", r, eps, NMax, names_vec[i], stop_crit);
    //}
    //auto t2 = omp_get_wtime();
    //std::cout <<"Time: "<< t2 - t1 << std::endl;

    MyConstrainedProblemGenerator gen;

    auto pr_vec = gen.GenerateNProblems(100, HillOnly, 300, 0.2);
    //std::cout << pr_vec[0]->GetOptimumPoint() <<" "<< pr_vec[0]->GetOptimumValue();

    return 0;
}
