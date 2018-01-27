#pragma once
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

vector<vector<int>> parallel_prefix(vector<vector<vector<int>>> Matrices, int rank, int p, int mod_by);
