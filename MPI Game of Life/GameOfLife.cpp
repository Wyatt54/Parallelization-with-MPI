// GameOfLife.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "functions.cpp"


int main(int argc, char *argv[])
{
	// Establishes what rank it is, and how many processes are running.
	int rank, p, n, G;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	
	// The first argument is the size of the matrix, and the second argument is how many generations you want to simulate.
	n = std::atoi(argv[1]);
	G = std::atoi(argv[2]);
	
	

	// These actions will only be done by the first process.
	// This includes printing the input variables, constructing 
	// the matrix, and sending the information to the other processes.
	if (rank == 0)
	{
		cout << "An " << n << " by " << n << " square matrix will be used." << endl;
		cout << G << " generations will be simulated." << endl;
		printf("Rank=%d: number of processes =%d\n", rank, p);
	}

	// Generate the board in its initial state
	vector<vector<int>> myBoard = GenerateInitialGoL(n);

	clock_t SimulationTime;
	// Simulates G generations
	SimulationTime = clock();
	Simulate(myBoard, n, G);
	SimulationTime = clock() - SimulationTime;
	
	cout << "Total Simulation + Display time: ";
	printf("%f\n", (float)SimulationTime/1000);

	MPI_Finalize();

	return 0;
}