#pragma once
#include "Minimazer2.h"

#include "MyConstrainedProblem.h"
#include "MyConstrainedProblemFamily.h"


class MyTester {
private:
    MyMinimazer Min;
    double exp_x, exp_z;
    ConsTrial deviation;
    double eps;
    int m;

public:
    MyTester(MyConstrainedProblem* MCPtr, std::vector<double> _r, double _eps = 0.01, uint64_t _NMax = 500) : Min(MCPtr, _r, _eps, _NMax) {
        eps = _eps;
        auto exp_tmp = MCPtr->GetOptimumPoint();
        m = MCPtr->GetNumberofConstr();
        if (abs(MCPtr->GetOptimumValue() - MCPtr->ComputeFunction(exp_tmp)) > 0.0001) {
            std::cout.precision(10);
            std::cout << MCPtr->GetOptimumValue() << " " << MCPtr->ComputeFunction(exp_tmp) << std::endl;
            std::cout << "INCORRECT OPTIMUM VALUE OR POINT" << std::endl;
        }
        //expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
        exp_x = exp_tmp;
        exp_z = MCPtr->ComputeFunction(exp_tmp);
    }

    bool Test(bool stop_crit) {
        ConsTrial res = Min.find_glob_min(stop_crit);

        double dev = abs(res.x - exp_x);
        bool r1 = (abs(dev) < eps);
        std::cout<<((r1)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return r1; // || abs(res.z - expected.z) < eps); // Is this LEGAL? Well, no, but...
    }

    bool Test_BF(uint64_t steps = 1000000) {
        ConsTrial res = Min.find_glob_min_BF(steps);

        double dev = abs(res.x - exp_x);
        bool r1 = (abs(dev) < eps);
        std::cout << ((r1) ? "YEEEEEEEEEEEEEEEEEES" : "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return r1; // || abs(res.z - expected.z) < eps); // Is this LEGAL? Well, no, but...
    }


    void Show_info() {
        if (Min.IsSolved()) {
            Min.Show_info();
            auto res = Min.GetMin();
            std::cout << "Ожидаемый результат: y = " << exp_z << " в точке x = " << exp_x << std::endl;
            std::cout << "Отклонение по x = " << res.x - exp_x << std::endl;
            std::cout << "Отклонение по y = " << res.z[m] - exp_z << std::endl;
            std::cout << std::endl;
        }
    }

    void Show_info_in_file(std::ofstream& fout) {
        //if (!Min.IsSolved()) Test();
        if (fout.is_open()) {
            Min.Show_info_in_file(fout);
            auto res = Min.GetMin();
            fout << "Ожидаемый результат: y = " << exp_z << " в точке x = " << exp_x << std::endl;
            fout << "Отклонение по x = " << res.x - exp_x << std::endl;
            fout << "Отклонение по y = " << res.z[m] - exp_z << std::endl;
            fout << std::endl;
        }
    }
    unsigned long long GetCount() { return Min.GetCount(); }

};

void Myfunc(MyConstrainedProblemFamily* IOPFPtr, std::string filepath, std::vector<double> r, double eps, uint64_t NMax, const std::string& family_name, bool stop_crit) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    //vector<int> CountVec1(NMax + 1);
    /*for (int i = 0; i <= NMax; ++i) {
        CountVec1[i] = 0;
    }*/
    for (size_t i = 0; i < IOPFPtr->GetFamilySize(); ++i) {
        //std::cout << "Тестируется " << family_name << " Problem" << i << std::endl;
        MyTester Tes(IOPFPtr->operator[](i), r, eps, NMax);
        bool tmp = Tes.Test(stop_crit);
        //bool tmp = Tes.Test_BF();

        if (tmp) {
            ++CorrectCount;
            std::cout << i << " YEP\n";
            Tes.Show_info();
            //++CountVec1[Tes.GetCount()];
        }
        else {
            std::cout << i << " NOPE\n";
            Tes.Show_info();
        }
    }
    std::cout << "Правильно решено " << CorrectCount << " из " << IOPFPtr->GetFamilySize() << " " << family_name
        << " family." << std::endl << std::endl;
    //file << "Правильно решено " << CorrectCount << " из " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;

    //file << "sep=,\n";
    //file << "0,0\n";
    double tmp;
    int i;
    /*for (i = 1; i < CountVec1.size(); ++i) {
        CountVec1[i] += CountVec1[i - 1];
    }
    for (i = 0; i < CountVec1.size(); ++i) {
        tmp = double(CountVec1[i]) / double(IOPFPtr->GetFamilySize());
        file << i << ',' << tmp << "\n";
        for (; (i + 1) < CountVec1.size() && CountVec1[i + 1] == CountVec1[i]; ++i);
    }*/
    //file << '\n';
}
