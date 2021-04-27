#pragma once
#include "MyConstrainedProblem.h"

class MyConstrainedProblemFamily {
	MyConstrainedProblemGenerator MCPGen;
	vector<MyConstrainedProblem*> fam;
public:
	MyConstrainedProblemFamily(uint n, MyConstrPrType type, uint m, double delta = 0.01, uint seed = 0) {
		fam = MCPGen.GenerateNProblems(n, type, m, delta, seed);
	}

	MyConstrainedProblem* operator[](int ind) const{
		return fam[ind];
	}

	int GetFamilySize() const {
		return fam.size();
	}

};
