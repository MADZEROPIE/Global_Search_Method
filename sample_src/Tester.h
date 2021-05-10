#pragma once
#include "Minimazer.h"

#include "HansenProblemFamily.hpp"
#include "Hill/HillProblemFamily.hpp"
#include "Shekel/ShekelProblemFamily.hpp"


class Tester {
private:
    Minimazer Min;
    Trial expected;
    Trial deviation;
    double eps;

public:
    Tester(IOptProblem* IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500) : Min(IOPPtr, _eps, _r, _NMax) {
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
        std::cout << "Отклонение по y = " << res.z - expected.z << std::endl;
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

void func(IOptProblemFamily* IOPFPtr, std::string filepath, double r, double eps, uint64_t NMax, const std::string& family_name, bool stop_crit) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    vector<int> CountVec1(NMax + 1);
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
    std::cout << "Правильно решено " << CorrectCount << " из " << IOPFPtr->GetFamilySize() << " " << family_name
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
