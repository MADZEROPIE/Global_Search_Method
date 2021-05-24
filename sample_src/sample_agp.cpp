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

int main(int argc, char* argv[]) {

    std::string filepath = "../Results/";
    //if (argc > 1) { filepath = argv[1]; }
    //
    //if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.001;
    //if (argc > 3) eps = std::stod(argv[3]);
    bool stop_crit = false;
    //if (argc > 4) {
    //    stop_crit = (argv[4] == "ON");
    //}

    setlocale(LC_ALL, "Russian");

    ////file.precision(6);

    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t NMax = 50000;
    THansenProblemFamily HFam;
    THillProblemFamily HillFam;
    TShekelProblemFamily ShekFam;

    vector<IOptProblemFamily*> vec = { &HFam, &HillFam, &ShekFam };

    vector <int> r_vec = { 2,3,4,5,6 };
    
    for(int r=2;r<6;++r)
    {  // Casual(?)
    vector<std::string> names_vec = { "Hansen" ,"Hill", "Shekel" };
    auto t1 = omp_get_wtime();
    for (size_t i = 0; i < vec.size(); ++i) {
        func(vec[i], filepath + names_vec[i]+"_r"+std::to_string(r) + ".csv", r, eps, NMax, names_vec[i], stop_crit);
    }
    auto t2 = omp_get_wtime();
    std::cout << "Time: " << t2 - t1 << std::endl;
    }

    int n = 30;
    int m = 4;
    vector<std::string> names_vecC = { "SheckConstr", "HillConstr" };
    std::cout << "Генерируется семейство функций...\n";
    auto type = SheckelOnly;
    MyConstrainedProblemFamily MCPFam1(n, type, m, 0.5, 11);
    std::cout << "Генерация завершена...\n";
    double r = 3.5;
    std::vector<int> q_vec = { 1, 10, 100 };

    for (auto q : q_vec) { // Constrained
        vector<double> cr(m + 1);
        for (int i = 0; i <= m; ++i) {
            cr[i] = r;
        }
        auto t_1 = omp_get_wtime();
        Myfunc(&MCPFam1, filepath + names_vecC[type] + "_q" + std::to_string(q) + ".csv", cr, eps, NMax, names_vecC[type], false, q);
        auto t_2 = omp_get_wtime();
        std::cout << "Time: " << t_2 - t_1 << std::endl;
    }

    for (auto r : r_vec)
    { // D
    TGrishaginProblemFamily famGr;
    func2(&famGr, filepath + "Grishagin" + "_r" + std::to_string(r) +".csv", r, eps*10, 5000000, "Grishagin", true);

    TGKLSProblemFamily famGKLS(2);
    func2(&famGKLS, filepath + "GKLS" + "_r" + std::to_string(r) + ".csv", r, eps, 5000000, "GKLS", true);

    }
    //{  // D = 1
    //    vector<std::string> names_vecD = { "HansenD" ,"HillD", "ShekelD" };
    //    auto t3 = omp_get_wtime();
    //    for (size_t i = 0; i < vec.size(); ++i) {
    //        //func2(vec[i], filepath + names_vecD[i] + ".csv", r, eps, NMax, names_vecD[i]);
    //    }
    //    auto t4 = omp_get_wtime();
    //    std::cout << "Time: " << t4 - t3 << std::endl;
    //}
    return 0;
}
