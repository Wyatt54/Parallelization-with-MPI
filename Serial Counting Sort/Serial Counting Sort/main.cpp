#include <string>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <vector>

int main(int argc, char *argv[])
{
	// Simple counting sort. First input is the Maximal value generated, and the second is
	// how many numbers to generate.
	int n, max, seed;
	std::vector<int> toSort, counts, sorted;
	clock_t totalTime, sortTime;
	if (argc > 1)
	{
		max = std::atoi(argv[1]);
		n = std::atoi(argv[2]);
		seed = std::atoi(argv[3]);
	}
	else
	{
		max = 100;
		n = 16777216;
		seed = 262;
	}
	std::cout << "The range is (0, " << max << ").  " << n << " numbers to be sorted." << std::endl;
	srand(time(NULL));
	
	totalTime = clock();
	// Fill up the array with random numbers bounded by the maximum
	toSort.resize(n);
	toSort[0] = ((137 * seed) + 865) % max;
	for (int i = 1; i < n; i++)
	{
		toSort[i] = ((34 * toSort[i - 1]) + 865) % max;
	}

	// Initialize the counting array to size of max. make all values zero.
	counts.resize(max);

	// Here is the sorting part. Go through the array toSort and increment
	// the count for each encounter.
	sortTime = clock();
	for (int i = 0; i < toSort.size(); i++)
	{
		counts[toSort[i]]++;
	}

	// Now that we have all the counts we output them to a sorted array.
	for (int i = 0; i < counts.size(); i++)
	{
		for (int j = 0; j < counts[i]; j++)
		{
			sorted.push_back(i);
		}
	}
	sortTime = clock() - sortTime;
	totalTime = clock() - totalTime;
	std::cout << "Total time taken: " << (float)totalTime / 1000 << " seconds." << std::endl;
	std::cout << "Time to sort: " << (float)sortTime / 1000 << " seconds." << std::endl;


	return 0;
}