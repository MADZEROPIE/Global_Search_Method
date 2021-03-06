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

public:
    Minimazer(IOptProblem*  _IOPPtr, double _eps = 0.01, double _r = 2.0) {
        IOPPtr = _IOPPtr;
        vector<double> tmp_lb, tmp_rb;
        IOPPtr->GetBounds(tmp_lb, tmp_rb);
        a = tmp_lb[0]; b = tmp_rb[0];
        r = (_r > 1.0) ? _r : 2.0;
        eps = _eps;
    }
   
    dpair find_glob_min() {
        vector<pair<double, double> > vec;
        vec.push_back(dpair(a, IOPPtr->ComputeFunction({ a })));
        vec.push_back(dpair(b, IOPPtr->ComputeFunction({ b })));
        count = 2;
        double M = 0;
        size_t k = 2;
        size_t t = 0;

        for (; ((vec[t + 1].first) - (vec[t].first)) > eps; ++k) {
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
        std::cout << "������� ������� = " << sol.second << " � ����� x = " << sol.first << std::endl;
        std::cout << "�������� ���������� eps = " << eps << std::endl;
        std::cout << "�������� ������ r = " << r << std::endl;
        std::cout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "������� ������� = " << sol.second << " � ����� x = " << sol.first << std::endl;
            fout << "�������� ���������� eps = " << eps << std::endl;
            fout << "�������� ������ r = " << std::fixed <<r << std::endl;
            fout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
        }
    }
    dpair GetMin() {
        if (!solved) this->find_glob_min();
        return sol;
    }

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
    Tester(IOptProblem* IOPPtr, double _eps = 0.01,double _r=2.0): Min(IOPPtr, _eps,_r) {
        eps = _eps;
        auto exp_tmp = IOPPtr->GetOptimumPoint();
        expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
    }

    bool Test() {
        dpair res=Min.find_glob_min();
        
        double dev = (res.first - expected.first);
       // std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps);
        
        
    }
    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res=Min.GetMin();
        std::cout << "��������� ���������: y = " << expected.second << " � ����� x = " << expected.first << std::endl;
        std::cout << "���������� �� x = " << res.first - expected.first << std::endl;
        std::cout << "���������� �� y = " << res.second -expected.second << std::endl;
        std::cout << std::endl;
    }

    void Show_info_in_file(std::ofstream& fout) {
        //if (!Min.IsSolved()) Test();
        if (fout.is_open()) {
            Min.Show_info_in_file(fout);
            auto res = Min.GetMin();
            fout << "��������� ���������: y = " << expected.second << " � ����� x = " << expected.first << std::endl;
            fout << "���������� �� x = " << res.first - expected.first << std::endl;
            fout << "���������� �� y = " << res.second - expected.second << std::endl;
            fout << std::endl;
        }
    }

};



int main(int argc,char* argv[]) {
   
    
    std::ofstream file;
    std::string filepath = "results.csv";
    
    if (argc > 1) filepath = argv[1];
    double r = 3.0;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.01;
    if(argc>3) eps = std::stod(argv[3]);

    setlocale(LC_ALL, "Russian");

    
    file.open(filepath, std::ofstream::app);
    //file.precision(6);

    uint64_t CorrectCount=0;
    THansenProblemFamily HFam;
    for (size_t i = 0; i < HFam.GetFamilySize();++i) {
        //std::cout << "����������� THansenProblem" << i << "..." << std::endl;
        Tester Tes(HFam[i], eps, r);
        if (Tes.Test()) ++CorrectCount;
        if (file.is_open()) {
            file << "THansenProblem" << i << std::endl;
            Tes.Show_info_in_file(file);
        }
    }
    std::cout << "��������� ������ " << CorrectCount << " �� " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;


    CorrectCount = 0;
    THillProblemFamily HillFam;
    for (size_t i = 0; i < HillFam.GetFamilySize(); ++i) {
        //std::cout << "����������� THillProblem" << i << "..." << std::endl;
        Tester Tes(HillFam[i], eps, r);
        if (Tes.Test()) ++CorrectCount;
        if (file.is_open()) {
            file << "THillProblem" << i << std::endl;
            Tes.Show_info_in_file(file);
        }
    }
    std::cout << "��������� ������ " << CorrectCount << " �� " << HillFam.GetFamilySize() << " THillProblem." << std::endl << std::endl;

    CorrectCount = 0;
    TShekelProblemFamily ShekFam;
    for (size_t i = 0; i < ShekFam.GetFamilySize(); ++i) {
        //std::cout << "����������� TShekelProblem" << i << "..." << std::endl;
        Tester Tes(ShekFam[i], eps, r);
        if (Tes.Test()) ++CorrectCount;
        if (file.is_open()) {
            file << "TShekelProblem" << i << std::endl;
            Tes.Show_info_in_file(file);
        }
    }
    std::cout << "��������� ������ " << CorrectCount << " �� " << ShekFam.GetFamilySize() << " TShekelProblem." << std::endl << std::endl;



    file.close();
    return 0;
}