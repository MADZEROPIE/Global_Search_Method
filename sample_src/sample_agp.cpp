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

//About style of naming. There is no style.

using std::vector;
using std::pair;

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
   
    dpair find_glob_min() {
        vector<pair<double, double> > vec;
        vec.push_back(dpair(a, IOPPtr->ComputeFunction({ a })));
        vec.push_back(dpair(b, IOPPtr->ComputeFunction({ b })));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;

        for (; abs(vec[t + 1].first - IOPPtr->GetOptimumPoint()[0]) > eps && k<NMax; ++k) {
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
        sol = vec[t+1]; solved = true;
        return vec[t+1];
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



class Tester {
private:
    Minimazer Min;
    dpair expected;
    dpair deviation;
    double eps;
public:
    Tester(IOptProblem* IOPPtr, double _eps = 0.01,double _r=2.0, uint64_t _NMax=500): Min(IOPPtr, _eps,_r,_NMax) {
        eps = _eps;
        auto exp_tmp = IOPPtr->GetOptimumPoint();
        expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
    }

    bool Test() {
        dpair res=Min.find_glob_min();
        
        double dev = (res.first - expected.first);
        //std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps);
        
        
    }
    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res=Min.GetMin();
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


void func(IOptProblemFamily* IOPFPtr,std::string filepath , double r, double eps, uint64_t NMax) {
	std::ofstream file;
	file.open(filepath);
    uint64_t CorrectCount = 0;
    vector<unsigned long long> CountVec1(IOPFPtr->GetFamilySize());
    for (size_t i = 0; i < IOPFPtr->GetFamilySize(); ++i) {
        //std::cout << "Тестируется THansenProblem" << i << "..." << std::endl;
        Tester Tes(IOPFPtr->operator[](i), eps, r);
        if (Tes.Test()) ++CorrectCount;
        Tes.Show_info();
        //if (file.is_open()) {
        //    file << "THansenProblem" << i << std::endl;
        //    Tes.Show_info_in_file(file);
       // }
        CountVec1[i] = Tes.GetCount();
    }
    std::cout << "Правильно решено " << CorrectCount << " из " << IOPFPtr->GetFamilySize() << "." << std::endl << std::endl;
    //file << "Правильно решено " << CorrectCount << " из " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;
    std::sort(CountVec1.begin(), CountVec1.end());

    file << "sep=,\n";
    for (size_t i = 0; i < CountVec1.size() && CountVec1[i] != NMax; ++i) {
        while (i + 1 < CountVec1.size() && CountVec1[i] == CountVec1[i + 1]) { ++i; }
        file << CountVec1[i] << ',' << double(i + 1) / double(CountVec1.size()) << "\n";
    }
    file << '\n';
}



int main(int argc,char* argv[]) {
   
    
    
    std::string filepath = "results";
    
    if (argc > 1){ filepath = argv[1];}
    double r = 2.0;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.01;
    if(argc > 3) eps = std::stod(argv[3]);

    setlocale(LC_ALL, "Russian");

    
   
    //file.precision(6);
    
    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t CorrectCount=0;
    uint64_t NMax=500;
    THansenProblemFamily HFam;
    THillProblemFamily HillFam;
    TShekelProblemFamily ShekFam;
    
    vector<IOptProblemFamily*> vec = { &HFam,&HillFam,&ShekFam };

    vector<std::string> names_vec = { "Hansen","Hill","Shekel" };

    for (size_t i=0;i<vec.size();++i) {
		func(vec[i], filepath + names_vec[i] + ".csv", r, eps, NMax);
	}


   
    return 0;
}