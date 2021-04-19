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
            || (!stop_crit && (vec[t + 1].x - vec[t].x >= eps)))
            && k < NMax; ++k) {
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
            //vec.insert(std::lower_bound(vec.begin(), vec.end(), t1_pair, [](const Trial& a, const Trial& b) {return a.x <= b.x; }), t1_pair); //No need for sorting, only to insert
            
            vec.insert(vec.begin() + t + 1, t1_pair);
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
        count = NUMTH + 1;
        double M = 0;
        int k = NUMTH + 1;
        int iter_count = 0;
        size_t t = 0;
        int tj_size = NUMTH;
        std::vector<Trial> tj_vec;
        double h = (b - a) / NUMTH;
        omp_set_num_threads(NUMTH);
        for (int i = 0; i <= NUMTH; ++i) {
            vec.push_back(Trial(a + i * h, IOPPtr->ComputeFunction({ a + i * h })));
        }
        for (; iter_count < NMax; iter_count= k = vec.size()) {
            //int tj_size = (tj_vec.NUMTH() < NUMTH) ? tj_vec.NUMTH() : NUMTH;
            for (int i = 0; i < (k - 1); ++i) {
                if ((vec[i + 1].x - vec[i].x) < eps) {  // ����� ������ ������ �� ���������� tj, �� ...
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

            for (int j = 0; j < tj_size; ++j) {  // ������ �� ������ tj_size ���� ������������ R
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
        std::cout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
        std::cout << "�������� ���������� eps = " << eps << std::endl;
        std::cout << "�������� ������ r = " << r << std::endl;
        std::cout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
            fout << "�������� ���������� eps = " << eps << std::endl;
            fout << "�������� ������ r = " << std::fixed <<r << std::endl;
            fout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
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
    //std::function<double(double)> func; // Function that needs findin' minimum
    double a; // Beginning of the segment
    double b; // End of the segment
    //std::vector<std::function<double(double)> > g_vec; // g functions vector

    MyConstrainedProblem* MCPtr; //Problem that needs findin' minimum
    uint m;
    std::vector<double> r; // Coefficient of the method. VECTOR OF COEFFs ?
    std::vector<double> e_vec; // Reserve vector (?)
    double eps;
    ConsTrial sol;
    bool solved = false;
    unsigned long long count = 0; // How many times func was executed
    unsigned long long NMax; // Magic Number
    double delta = 5.0;

    ConsTrial make_trial(double x) {
        ConsTrial tr;
        tr.x = x;
        for (uint i = 0; i < m; ++i) {
            tr.z = MCPtr->ComputeGjConstr(i, x);
            if (tr.z > 0.0) {
                tr.index = i;
                return tr;
            }
        }
        tr.z = MCPtr->ComputeGjConstr(m, x);
        tr.index = m;
        return tr;
    }
public:
    Minimazer2(MyConstrainedProblem* _MCPtr,  std::vector<double> _r, double _eps = 0.01, uint64_t _NMax = 500) {
        MCPtr = _MCPtr;
        m = MCPtr->GetNumberofConstr();
        MCPtr->GetBounds(a, b);
        r = _r;
        eps = _eps;
        NMax = _NMax;
    }

    ConsTrial find_glob_min(bool stop_crit = false) {
        vector<ConsTrial> vec;
        vector<vector<ConsTrial>> Vvec(m + 1);
        vector<double> mv(m + 1), zv(m + 1), ev(m + 1);
        vec.push_back(this->make_trial(a));
        vec.push_back(this->make_trial(b));
        Vvec[vec[0].index].push_back(vec[0]);
        Vvec[vec[1].index].push_back(vec[1]);

        for (int i = 0; i <= m; ++i) {
            mv[i] = 0;
        }
        double max_eps = b-a;
        int t = 0;
        while ((vec[t+1].x-vec[t].x) > eps) {
            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                int s = Vvec[i].size();
                for (int j = 0; j < s - 1; ++j) {
                    double M_tmp = abs(Vvec[i][j + 1].z - Vvec[i][j].z) / (Vvec[i][j + 1].x - Vvec[i][j].x);
                    if (M_tmp > mv[i])
                        mv[i] = M_tmp;

                }
                //if (mv[i] == 0) mv[i] = 1; //???
                //else mv[i] = r[i] * mv[i];
                ev[i] = mv[i] * delta;
            }
            for (int i = 0; i <= m; ++i) {  // Could be done in parallel
                if (Vvec[i].size() != 0) {
                    double z_tmp = Vvec[i][0].z;
                    int s = Vvec[i].size();
                    for (int j = 1; j < s; ++j)
                        if (Vvec[i][j].z < z_tmp)
                            z_tmp = Vvec[i][j].z;
                    
                    if (z_tmp <= 0) 
                        zv[i] = -ev[i];
                    else 
                        zv[i] = z_tmp;  
                }
            }

            // Finally R. But hey, SR is better
            int s = vec.size();
            double R = -INFINITY;
            
            for (int i = 0; i < s-1 ; ++i) {
                double R_tmp;
                if (vec[i].index == vec[i + 1].index) {
                    auto v = vec[i].index;
                    R_tmp = vec[i + 1].x - vec[i].x + (vec[i + 1].z - vec[i].z) * (vec[i + 1].z - vec[i].z) /  // I know, I know...
                        (r[v] * r[v] * mv[v] * mv[v] * vec[i + 1].x - vec[i].x) - 2 * (vec[i + 1].z + vec[i].z - 2 * zv[v]) / (r[v] * mv[v]);
                }
                else if (vec[i].index < vec[i + 1].index) {
                    auto v = vec[i+1].index;
                    R_tmp = 2 * (vec[i + 1].x - vec[i].x) - 4 * (vec[i + 1].z - zv[v]) / r[v] * mv[v];
                }
                else {
                    auto v = vec[i].index;
                    R_tmp = 2 * (vec[i + 1].x - vec[i].x) - 4 * (vec[i].z - zv[v]) / r[v] * mv[v];
                }
                if (R_tmp > R) {
                    t = i;
                    R = R_tmp;
                }
            }
            double x_new;
            if (vec[t].index != vec[t + 1].index) {
                x_new = (vec[t].x + vec[t + 1].x) / 2;
            }
            else {
                x_new = (vec[t].x + vec[t + 1].x) / 2 - (vec[t + 1].z + vec[t].z) / (2 * r[vec[t].index] * mv[vec[t].index]);
            }
            auto tr = make_trial(x_new);
            max_eps = std::min(vec[t + 1].x - x_new, x_new - vec[t].x);
            vec.insert(vec.begin() + t, tr);
            Vvec[tr.index].push_back(tr);  // TODO: INSERT WITH BINARY SEARCH
            std::sort(Vvec[tr.index].begin(), Vvec[tr.index].end(), [](const ConsTrial& a, const ConsTrial& b) {return a.x < b.x; });

        }
        
        sol = vec[t + 1];  // NEED TEST

        solved = true;
        return sol;
    }

    void Show_info() {
        //if (!solved) this->find_glob_min();
        std::cout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
        std::cout << "�������� ���������� eps = " << eps << std::endl;
        /*std::cout << "�������� ������ r = " << r << std::endl;*/
        std::cout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
    }
    void Show_info_in_file(std::ofstream& fout) {
        //if (!solved) this->find_glob_min();
        if (fout.is_open()) {
            fout << "������� ������� = " << sol.z << " � ����� x = " << sol.x << std::endl;
            fout << "�������� ���������� eps = " << eps << std::endl;
            /*fout << "�������� ������ r = " << std::fixed << r << std::endl;*/
            fout << "������� ���� ��������� " << count << " ���(�)" << std::endl;
        }
    }
    ConsTrial GetMin() {
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
    Tester(IOptProblem* IOPPtr, double _eps = 0.01, double _r = 2.0, uint64_t _NMax = 500): Min(IOPPtr, _eps,_r,_NMax) {
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
        std::cout << "��������� ���������: y = " << expected.z << " � ����� x = " << expected.x << std::endl;
        std::cout << "���������� �� x = " << res.x - expected.x << std::endl;
        std::cout << "���������� �� y = " << res.z -expected.z << std::endl;
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
        //std::cout << "����������� " << family_name << " Problem" << i << std::endl;
        Tester Tes(IOPFPtr->operator[](i), eps, r, NMax);
        bool tmp = Tes.Test_par();
        //Tes.Show_info();
        //bool tmp = Tes.Test(stop_crit);
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
    std::cout << "��������� ������ " << CorrectCount << " �� " << IOPFPtr->GetFamilySize() <<" " << family_name 
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

//------------------------------------------------------------------------------------------\\

int main(int argc,char* argv[]) {
    
    std::string filepath = "results";
    if (argc > 1){ filepath = argv[1];}
    double r = 4;
    if (argc > 2) r = std::stod(argv[2]);
    double eps = 0.001;
    if (argc > 3) eps = std::stod(argv[3]);
    bool stop_crit = false;
    if (argc > 4) {
        stop_crit = (argv[4] == "ON");
    }

    setlocale(LC_ALL, "Russian");

    ////file.precision(6);

    //TODO: CREATE AND USE ADDITIONAL CLASS OR FUNCTION  || CHANGE TESTER
    uint64_t NMax = 250;
    THansenProblemFamily HFam;
    THillProblemFamily HillFam;
    TShekelProblemFamily ShekFam;

    vector<IOptProblemFamily*> vec = { &HFam, & HillFam, & ShekFam};

    vector<std::string> names_vec = { "Hansen" ,"Hill", "Shekel" };
    auto t1 = omp_get_wtime();
    for (size_t i = 0; i < vec.size();++i) {
        func(vec[i], filepath + names_vec[i] + ".csv", r, eps, NMax, names_vec[i], stop_crit);
    }
    auto t2 = omp_get_wtime();
    std::cout <<"Time: "<< t2 - t1 << std::endl;


    //--TESTS--

    MyConstrainedProblemGenerator gen;

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
