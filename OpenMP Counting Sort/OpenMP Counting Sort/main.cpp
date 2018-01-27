#include <stdio.h>
#include <time.h>
#include <ctime>
#include <omp.h>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace std;


void CountingSort(int, int, int, int);

int main(int argc, char *argv[])
{
	cout.precision(20);
	omp_set_dynamic(0);
	int max, n, p, seed;

	// Aruguments are, max, how many to sort, how many threads, and the seed.
	if (argc > 1)
	{
		max = atoll(argv[1]);
		n = atoi(argv[2]);
		p = atoi(argv[3]);
		seed = atoi(argv[4]);
	}
	else
	{
		max = 1000;
		n = 100;
		p = 4;
		seed = 262;
	}
	std::cout << "The range is (0, " << max << ").  " << n << " numbers to be sorted, with " << p << " threads." << std::endl;
	srand(time(NULL));

	double time = omp_get_wtime();

	CountingSort(max, n, p, seed);

	time = omp_get_wtime() - time;
	printf("Total time = %f seconds \n ", time);

	return 0;
}


void CountingSort(int max, int n, int p, int seed) {

	omp_set_num_threads(p);
	std::vector<int> toSort, counts;
	vector<vector<int>> sorted;
	int perProcess = n / p;
	cout << perProcess << endl;


	// Fill up the array with random numbers bounded by the maximum
	toSort.resize(n);
	toSort[0] = ((137 * seed) + 94732) % max;
	for (int i = 1; i < n; i++)
	{
		toSort[i] = (137 * toSort[i - 1] + 94732) % max;
		cout << toSort[i] << " ";
	}

	// Initialize the counting array to size of max. make all values zero.
	counts.resize(max);
	sorted.resize(p);

	double timeStatic = omp_get_wtime();
#pragma omp parallel for schedule(static)  num_threads(p)
	for (int i = 0; i < toSort.size(); i++)
	{
		counts[toSort[i]]++;
	}
	cout << "Counts done." << endl;
#pragma omp parallel for schedule(static) num_threads(p) shared(sorted)
	for (int i = 0; i < counts.size(); i++)
	{
		int me = omp_get_thread_num();
		int number = i;
		int times = counts[i];
		for (int j = 0; j < times; j++)
		{
			sorted[me].push_back(number);
		}
	}
	// Now we have to concatenate the vectors into a single sorted vector.
	vector<int> final;
	if (p > 1)
	{
		final.reserve(n);
		for (int i = 0; i < sorted.size(); i++)
		{
			final.insert(final.end(), sorted[i].begin(), sorted[i].end());
		}
	}
	timeStatic = omp_get_wtime() - timeStatic;


	cout << "Time Taken to sort: " << timeStatic << " seconds." << endl;
	for (int i = 0; i < final.size(); i++)
	{
		cout << final[i] << " ";
	}

	return;

}