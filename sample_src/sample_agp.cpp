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

#include "Tester.h"

#include "MyConstrainedProblemFamily.h"
#include "Tester2.h"

#include "TesterD.h"
#include "GKLS/GKLSProblemFamily.hpp"


//------------------------------------------------------------------------------------------\\

int main(int argc,char* argv[]) {
    
    std::string filepath = "results";
    if (argc > 1){ filepath = argv[1];}
    double r = 4.5;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.0001;
    if (argc > 3) eps = std::stod(argv[3]);
    bool stop_crit = false;
    if (argc > 4) {
        stop_crit = (argv[4] == "ON");
    }

    setlocale(LC_ALL, "Russian");

    ////file.precision(6);

    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t NMax = 500;
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


    //--TESTS--

    int n = 5;
    int m = 4;
    std::cout << "Генерируется семейство функций...\n";
    MyConstrainedProblemFamily MCPFam1(n, HillOnly, m, 5.5, 11);
    std::cout << "Генерация завершена...\n";
    vector<double> cr(m + 1);
    for (auto& r1 : cr) {
        r1 = r;
    }
    auto t1 = omp_get_wtime();
    Myfunc(&MCPFam1, "HillConstr", cr, eps, NMax, "HillConstr", false);
    auto t2 = omp_get_wtime();
    std::cout << "Time: " << t2 - t1 << std::endl;

    //TGrishaginProblemFamily fam;
    //func2(&fam, "Grishagin.txt", r, eps, 5000000, "Grishagin");

    //int a;
    //std::cin >> a;

    //TGKLSProblemFamily fam2(3);
    //func2(&fam2, "GKLS.txt", r, eps, 50000000000, "GKLS", false);

    //THansenProblemFamily HFam;
    //THillProblemFamily HillFam;
    //TShekelProblemFamily ShekFam;

    //vector<IOptProblemFamily*> vec = { &HFam, & HillFam, & ShekFam};

    //vector<std::string> names_vecD = { "HansenD" ,"HillD", "ShekelD" };
    //auto t1 = omp_get_wtime();
    //for (size_t i = 0; i < vec.size();++i) {
    //    func2(vec[i], filepath + names_vecD[i] + ".csv", r, eps, NMax, names_vecD[i]);
    //}
    //auto t2 = omp_get_wtime();
    //std::cout <<"Time: "<< t2 - t1 << std::endl;



    // ----- OLD THINGS. DELETE BEFORE RELEASE -----

    /*//int maxind = 0;
    //double maxdiff = 0.0;
    //for (int i = 0; i < ShekFam.GetFamilySize();++i) {
    //    auto IOPPtr = ShekFam.operator[](i);
    //    double a, b;
    //    vector<double> tmp_lb, tmp_rb;
    //    IOPPtr->GetBounds(tmp_lb, tmp_rb);
    //    a = tmp_lb[0]; b = tmp_rb[0];
    //    double h = (b - a) / 100000;
    //    double mz = IOPPtr->ComputeFunction({ a }), mx=a;
    //    for (double x = a + h; x <= b; x += h) {
    //        double z = IOPPtr->ComputeFunction({ x });
    //        if (z > mz) {
    //            mz = z;
    //            mx = x;  
    //        }
    //    }
    //    if (IOPPtr->ComputeFunction({ b }) > mz) {
    //        mz = IOPPtr->ComputeFunction({ b });
    //        mx = b;
    //    }
    //    double pmx = IOPPtr->GetMaxPoint()[0];
    //    double pmz = IOPPtr->GetMaxValue();
    //    if (abs(pmz - mz) > maxdiff) {
    //        maxdiff = abs(pmz - mz);
    //        maxind = i;
    //    }
    //}
    //std::cout << maxdiff << " " << maxind;*/

    //auto pr_vec = gen.GenerateNProblems(1000000, SheckelOnly, 11, 0.35);
    //std::cout << pr_vec[0]->GetOptimumPoint() <<" "<< pr_vec[0]->GetOptimumValue();

    //std::vector<int> broken = { 133, 452, 800, 603, 688, 398, 875, 569, 304, 990, 843, 557 };
/*
    double minmaxSh = maxShekel[broken[0]][0], maxmaxSh = maxShekel[broken[0]][0];
    uint j = 0;
    for (uint i = 0; i < 11u; ++i) {
        if (minmaxSh > maxShekel[broken[i]][0]) {
            minmaxSh = maxShekel[broken[i]][0];
            j = i;
        }
        if (maxmaxSh < maxShekel[broken[i]][0]) 
            maxmaxSh = maxShekel[broken[i]][0];
    }
    std::cout << j<<" "<<minmaxSh << "\n";
    //minmaxSh = -0.2464;
    double delta = 0.2;

    int m = 11;
    double LoBound = 0.0;
    double UpBound = 10.0;
    double h = 1e-5;
    double x = LoBound;
    delta += minmaxSh;
    std::cout << delta << '\n';
    bool flag = false;
    double mindiff = 10000000000.0, maxdiff=0.0;
    uint64_t incl=0, all=0;
    for (; x <= UpBound; x += h) {
        double z;
        int i = 0;
        for (; i < m; ++i) {
            double tmp = ShekFam[broken[i]]->ComputeFunction({ x }) - delta;
            if (tmp > 0.0) { 
                if (tmp < mindiff) mindiff = tmp;
                if (tmp > maxdiff) maxdiff = tmp;
                //std::cout <<i<<" "<<x<<" "<< ShekFam[broken[i]]->ComputeFunction({ x }) - delta <<"\n";
                break;
            }
        }
        if (i == m) {
            ++incl;
            flag = true;
        }
        ++all;
    }
    std::cout << incl << " " << all << " " << double(incl) / all * 100.0 << "%\n";
    std::cout << maxdiff << " " << maxmaxSh << " " << maxmaxSh - delta << "\n";
    std::cout << mindiff << " " << flag;

    std::cout << "\n\n\n\n\n\n";
    */

    //double minmaxH = maxHill[broken[0]][0], maxmaxH = maxHill[broken[0]][0];
    //uint j = 0;
    //for (uint i = 0; i < 11u; ++i) {
    //    if (minmaxH > maxHill[broken[i]][0]) {
    //        minmaxH = maxHill[broken[i]][0];
    //        j = i;
    //    }
    //    if (maxmaxH < maxHill[broken[i]][0])
    //        maxmaxH = maxHill[broken[i]][0];
    //}
    //std::cout << j << " " << minmaxH << "\n";
    ////minmaxSh = -0.2464;
    //double delta = 0.353;

    //int m = 11;
    //double LoBound = 0.0;
    //double UpBound = 1.0;
    //double h = 1e-5;
    //double x = LoBound;
    //delta += minmaxH;
    //std::cout << delta << '\n';
    //bool flag = false;
    //double mindiff = 10000000000.0, maxdiff = 0.0;
    //uint64_t incl = 0, all = 0;
    //for (; x <= UpBound; x += h) {
    //    double z;
    //    int i = 0;
    //    for (; i < m; ++i) {
    //        double tmp = HillFam[broken[i]]->ComputeFunction({ x }) - delta;
    //        if (tmp > 0.0) {
    //            if (tmp < mindiff) mindiff = tmp;
    //            if (tmp > maxdiff) maxdiff = tmp;
    //            //std::cout <<i<<" "<<x<<" "<< ShekFam[broken[i]]->ComputeFunction({ x }) - delta <<"\n";
    //            break;
    //        }
    //    }
    //    if (i == m) {
    //        ++incl;
    //        flag = true;
    //    }
    //    ++all;
    //}
    //std::cout << incl << " " << all << " " << double(incl) / all * 100.0 << "%\n";
    //std::cout << maxdiff << " " << maxmaxH << " " << maxmaxH - delta << "\n";
    //std::cout << mindiff << " " << flag;
    return 0;
}
