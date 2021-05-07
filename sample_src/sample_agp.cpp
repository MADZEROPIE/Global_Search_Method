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
    double r = 2.5;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.01;
    if (argc > 3) eps = std::stod(argv[3]);
    bool stop_crit = false;
    if (argc > 4) {
        stop_crit = (argv[4] == "ON");
    }

    setlocale(LC_ALL, "Russian");

    ////file.precision(6);

    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t NMax = 5000;
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

    int n = 1;
    int m = 4;
    r = 4.5;
    eps = 0.01;
    std::cout << "Генерируется семейство функций...\n";
    MyConstrainedProblemFamily MCPFam1(n, SheckelOnly, m, 0.5, 11);
    std::cout << "Генерация завершена...\n";
    vector<double> cr(m + 1);
    for (int i = 0; i <= m;++i) {
        cr[i] = r;
    }
    auto t1 = omp_get_wtime();
    Myfunc(&MCPFam1, "HillConstr", cr, eps, NMax, "HillConstr", false);
    auto t2 = omp_get_wtime();
    std::cout << "Time: " << t2 - t1 << std::endl;

    /*TGrishaginProblemFamily fam;
    func2(&fam, "Grishagin.txt", r, eps, 5000000, "Grishagin", true);*/

    //int a;
    //std::cin >> a;

    //TGKLSProblemFamily fam2(2);
    //func2(&fam2, "GKLS.csv", r, eps, 50000000000, "GKLS", true);

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

    return 0;
}
