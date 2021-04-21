#pragma once
#include "Minimazer2.h"

#include "MyConstrainedProblem.h"
#include "MyConstrainedProblemFamily.h"


class Tester2 {
private:
    Minimazer2 Min;
    ConsTrial expected;
    ConsTrial deviation;
    double eps;

public:
    Tester2(MyConstrainedProblem* MCPtr, std::vector<double> _r, double _eps = 0.01, uint64_t _NMax = 500) : Min(MCPtr, _r, _eps, _NMax) {
        eps = _eps;
        auto exp_tmp = MCPtr->GetOptimumPoint();
        /*if (abs(IOPPtr->GetOptimumValue() - IOPPtr->ComputeFunction(exp_tmp)) > 0.001) {
            std::cout.precision(10);
            std::cout << IOPPtr->GetOptimumValue() << " " << IOPPtr->ComputeFunction(exp_tmp) << std::endl;
            std::cout << "INCORRECT OPTIMUM VALUE OR POINT" << std::endl;
        }*/
        //expected = std::make_pair(exp_tmp[0], IOPPtr->GetOptimumValue());
        expected = ConsTrial(exp_tmp, MCPtr->ComputeFunction(exp_tmp));
    }

    bool Test(bool stop_crit) {
        ConsTrial res = Min.find_glob_min(stop_crit);

        double dev = (res.x - expected.x);
        //std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps || abs(res.z - expected.z) < 7e-5); // Is this LEGAL? Well, no, but...
    }


    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res = Min.GetMin();
        std::cout << "��������� ���������: y = " << expected.z << " � ����� x = " << expected.x << std::endl;
        std::cout << "���������� �� x = " << res.x - expected.x << std::endl;
        std::cout << "���������� �� y = " << res.z - expected.z << std::endl;
        std::cout << std::endl;
    }

    void Show_info_in_file(std::ofstream& fout) {
        //if (!Min.IsSolved()) Test();
        if (fout.is_open()) {
            Min.Show_info_in_file(fout);
            auto res = Min.GetMin();
            fout << "��������� ���������: y = " << expected.z << " � ����� x = " << expected.x << std::endl;
            fout << "���������� �� x = " << res.x - expected.x << std::endl;
            fout << "���������� �� y = " << res.z - expected.z << std::endl;
            fout << std::endl;
        }
    }
    unsigned long long GetCount() { return Min.GetCount(); }

};

void func2(MyConstrainedProblemFamily* IOPFPtr, std::string filepath, std::vector<double> r, double eps, uint64_t NMax, const std::string& family_name, bool stop_crit) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    vector<int> CountVec1(NMax + 1);
    for (int i = 0; i <= NMax; ++i) {
        CountVec1[i] = 0;
    }
    for (size_t i = 0; i < IOPFPtr->GetFamilySize(); ++i) {
        //std::cout << "����������� " << family_name << " Problem" << i << std::endl;
        Tester2 Tes(IOPFPtr->operator[](i), r, eps, NMax);
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
    std::cout << "��������� ������ " << CorrectCount << " �� " << IOPFPtr->GetFamilySize() << " " << family_name
        << " family." << std::endl << std::endl;
    //file << "��������� ������ " << CorrectCount << " �� " << HFam.GetFamilySize() << " THansenProblem." << std::endl << std::endl;

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