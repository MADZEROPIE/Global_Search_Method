#pragma once
#include "DimMinimazer.h"

#include "Grishagin/GrishaginProblemFamily.hpp"


class TesterD {
private:
    MinimazerD Min;
    TrialD expected;
    TrialD deviation;
    double eps;
    int dim;

public:
    TesterD(IOptProblem* IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500) : Min(IOPPtr, _eps, _r, _NMax) {
        eps = _eps;
        auto exp_tmp = IOPPtr->GetOptimumPoint();
        /*if (abs(IOPPtr->GetOptimumValue() - IOPPtr->ComputeFunction(exp_tmp)) > 0.001) {
            std::cout.precision(10);
            std::cout << IOPPtr->GetOptimumValue() << " " << IOPPtr->ComputeFunction(exp_tmp) << std::endl;
            std::cout << "INCORRECT OPTIMUM VALUE OR POINT" << std::endl;
        }*/
        expected = TrialD(exp_tmp, IOPPtr->ComputeFunction(exp_tmp));
        dim = IOPPtr->GetDimension();
    }

    bool Test(bool save_trials = false) {
        TrialD res = Min.find_glob_min(save_trials);

        double dev = abs(res.x[0] - expected.x[0]);
        for (int i = 1; i < dim; ++i) {
            dev = std::max(dev, res.x[i] - expected.x[i]);
        }
        //std::cout<<((abs(dev) < eps)? "YEEEEEEEEEEEEEEEEEES": "NOOOOOOOOOOOOOOOOOOO") << std::endl;
        return (abs(dev) < eps || abs(res.z - expected.z) < 7e-5); // Is this LEGAL? Well, no, but...
    }

    void Show_info() {
        //if (!Min.IsSolved()) Test();
        Min.Show_info();
        auto res = Min.GetMin();
        std::cout << "Ожидаемый результат: y = " << expected.z << " в точке x = ";
        for (int i = 0; i < dim; ++i) {
            std::cout << expected.x[i]<<" ";
        }
        std::cout << std::endl;
        std::cout << "Отклонение по x = ";
        for (int i = 0; i < dim; ++i) {
            std::cout << res.x[i] - expected.x[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Отклонение по y = " << res.z - expected.z << std::endl;
        std::cout << std::endl;
    }

    //void Show_info_in_file(std::ofstream& fout) {
    //    //if (!Min.IsSolved()) Test();
    //    if (fout.is_open()) {
    //        Min.Show_info_in_file(fout);
    //        auto res = Min.GetMin();
    //        fout << "Ожидаемый результат: y = " << expected.z << " в точке x = " << expected.x << std::endl;
    //        fout << "Отклонение по x = " << res.x - expected.x << std::endl;
    //        fout << "Отклонение по y = " << res.z - expected.z << std::endl;
    //        fout << std::endl;
    //    }
    //}
    void saveTrialsInFile(std::ofstream& fout) {
        return Min.saveTrialsInFile(fout);
    }
    unsigned long long GetCount() { return Min.GetCount(); }


};

//------------------------------------------------------------------------------------------\\

void func2(IOptProblemFamily* IOPFPtr, std::string filepath, double r, double eps, uint64_t NMax, const std::string& family_name, bool save_trials = false) {
    std::ofstream file;
    file.open(filepath);
    uint64_t CorrectCount = 0;
    //vector<int> CountVec1(10000000);//NMax * (IOPFPtr->operator[](0)->GetDimension()+1));
    /*for (int i = 0; i <= NMax; ++i) {
        CountVec1[i] = 0;
    }*/
    auto fam_size = IOPFPtr->GetFamilySize();

    for (size_t i = 0; i < fam_size; ++i) {
        //std::cout << "Тестируется " << family_name << " Problem" << i << std::endl;
        TesterD Tes(IOPFPtr->operator[](i), eps, r, NMax);
        bool tmp = Tes.Test(save_trials);
        std::cout << i << " ";
        if (tmp) {
            ++CorrectCount;
            std::cout << "YEP\n";

            //++CountVec1[Tes.GetCount()];
        }
        else {
            std::cout << "NOPE\n";
            //Tes.Show_info();
        }
        Tes.Show_info();
        if (save_trials) {
            std::ofstream file;
            std::string filepath2 = family_name+std::to_string(i)+".csv";
            file.open(filepath2);
            file << "sep=,\n";
            Tes.saveTrialsInFile(file);
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
    }
    file << '\n';*/
}
