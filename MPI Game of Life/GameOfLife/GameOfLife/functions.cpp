#include "stdafx.h"
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

// This function displays the board for a single process.
// Won't be used in the final project, just for my own testing.
static void PrintMyBoard(vector<vector<int>> board, int columns, int rows)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			cout << board[j][i] << " ";
			if (j == columns - 1)
			{
				cout << endl;
			}
		}
	}
}

// This function generates p random numbers and distributes
// them across all available ranks.  Then each rank uses
// the random number as a seed to generate their portion
// of the board.
static vector<vector<int>> GenerateInitialGoL(int n)
{
	int rank, p;
	int mySeed = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	clock_t RandomSendTime, RandomReceiveTime;
	// First generate p random numbers to use as seeds
	// for generating the board

	// These actions are only completed by rank 0
	if (rank == 0)
	{
		// Generate the numbers
		vector<int> r_numbers;
		r_numbers.resize(p);
		srand(time(NULL));

		for (int i = 0; i < p; i++)
		{
			r_numbers[i] = rand();
		}
		mySeed = r_numbers[0];

		// Now Distribute them
		RandomSendTime = clock();
		for (int i = 1; i < p; i++)
		{
			MPI_Send(&r_numbers[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		RandomSendTime = clock() - RandomSendTime;
	}
	// All the other ranks recieve the random seed
	else if (rank > 0)
	{
		MPI_Status status;
		RandomReceiveTime = clock();
		MPI_Recv(&mySeed, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		RandomReceiveTime = clock() - RandomReceiveTime;
	}

	// Now all ranks have a different random seed.
	// Now each rank needs to generate their portion
	// of the board.

	// First implant the seed and determine how many
	// columns each process gets.
	srand(mySeed);
	int columns = n / p;

	// Create a 2D vector with x = columns, and y = n
	vector<vector<int>> myBoard;
	myBoard.resize(columns);
	for (int i = 0; i < columns; i++)
	{
		myBoard[i].resize(n);
	}
	

	// Now loop through the board and give each cell
	// a value of 0 for dead, or 1 for alive,
	// generated randomly
	for (int i = 0; i < columns; i++)
	{
		for (int j = 0; j < n; j++)
		{
			myBoard[i][j] = rand() % 2;
		}
	}
	// PrintMyBoard(myBoard, columns, n); Testing
	if (rank == 0)
	{
		float sTime = (float)RandomSendTime / 1000;
		cout << "Random number send time: " << sTime << endl;
	}
	else if (rank > 0)
	{
		float rTime = (float)RandomReceiveTime / 1000;
		cout << "Random number receive time: " << rTime << endl;
	}
	

	return myBoard;
}

// This function determines the state of the cell based
// off of the previous generation of the board.
static int DetermineState(vector<vector<int>> oldBoard, int x, int y, vector<int> prevColumn, vector<int> nextColumn, int columns)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// I'm splitting the neighbor check into 3 categories.  Check left
	// is NW, W, SW, Check middle is N, S, and check right is NE, E, and SE.
	int livingNeighbors = 0, state = 0;
	state = oldBoard[x][y];
	// Check left section.

	// if x is zero then all the left neighbors are in the previous column.
	if (x == 0)
	{
		// Check W
		livingNeighbors += prevColumn[y];
		
		// Check NW
		// Check if the current cell is at the top
		if (y == 0)
		{
			livingNeighbors += prevColumn[prevColumn.size() - 1];
		}
		else
		{
			livingNeighbors += prevColumn[y - 1];
		}

		// Check SW
		// Check to see if the current cell is at the bottom
		if (y == prevColumn.size() - 1)
		{
			livingNeighbors += prevColumn[0];
		}
		else
		{
			livingNeighbors += prevColumn[y + 1];
		}
	}
	else
	{
		// Here the previous column is on the board belonging
		// to the current process
		// Check NW
		if (y == 0)
		{
			livingNeighbors += oldBoard[x-1][prevColumn.size() - 1];
		}
		else
		{
			livingNeighbors += oldBoard[x-1][y - 1];
		}

		// Check W
		livingNeighbors += oldBoard[x - 1][y];

		// Check SW
		if (y == prevColumn.size() - 1)
		{
			livingNeighbors += oldBoard[x-1][0];
		}
		else
		{
			livingNeighbors += oldBoard[x-1][y + 1];
		}
	}

	// Check the middle section
	// Check N
	if (y == 0)
	{
		livingNeighbors += oldBoard[x][prevColumn.size() - 1];
	}
	else
	{
		livingNeighbors += oldBoard[x][y - 1];
	}
	// Check S
	if (y == prevColumn.size() - 1)
	{
		livingNeighbors += oldBoard[x][0];
	}
	else
	{
		livingNeighbors += oldBoard[x][y + 1];
	}

	// Check the right section
	// if x is n/p then all the right neighbors are in the next column.
	if (x == columns - 1)
	{
		// Check E
		livingNeighbors += nextColumn[y];

		// Check NE
		// Check if the current cell is at the top
		if (y == 0)
		{
			livingNeighbors += nextColumn[nextColumn.size() - 1];
		}
		else
		{
			livingNeighbors += nextColumn[y - 1];
		}

		// Check SE
		// Check to see if the current cell is at the bottom
		if (y == nextColumn.size() - 1)
		{
			livingNeighbors += nextColumn[0];
		}
		else
		{
			livingNeighbors += nextColumn[y + 1];
		}
	}
	else
	{
		// Here the next column is on the board belonging
		// to the current process
		// Check NE
		if (y == 0)
		{
			livingNeighbors += oldBoard[x + 1][prevColumn.size() - 1];
		}
		else
		{
			livingNeighbors += oldBoard[x + 1][y - 1];
		}

		// Check E
		livingNeighbors += oldBoard[x + 1][y];

		// Check SE
		if (y == prevColumn.size() - 1)
		{
			livingNeighbors += oldBoard[x + 1][0];
		}
		else
		{
			livingNeighbors += oldBoard[x + 1][y + 1];
		}
	}
	// Prints out living neighbors.  Just for testing.
	// cout << "Living Neighbors " << livingNeighbors << endl;

	// Now we determine the new state of the cell based off of
	// the old state, and the number of living neighbors.
	if (state == 1)
	{
		// Cell dies with either less than 2 or greater than 3
		// living neighbors.
		if (livingNeighbors < 2 || livingNeighbors > 3)
		{
			//cout << "Rank: " << rank << " Location " << x << y << " Old state: " << state << " Living Neighbors: " << livingNeighbors << " New state: 0" << endl;
			state = 0;
		}
		else
		{
			//cout << "Rank: " << rank << " Location " << x << y << " Old state: " << state << " Living Neighbors: " << livingNeighbors << " New state: " << state << endl;
		}
	}
	else if (state == 0)
	{
		if (livingNeighbors == 3)
		{
			//cout << "Rank: " << rank << " Location " << x << y << " Old state: " << state << " Living Neighbors: " << livingNeighbors << " New state: 1" << endl;
			// Comes back to life with 3 living neighbors
			state = 1;
		}
		else
		{
			//cout << "Rank: " << rank << " Location " << x << y << " Old state: " << state << " Living Neighbors: " << livingNeighbors << " New state: " << state << endl;
		}
	}

	return state;
}

// This function gathers the entire board into rank 0
// and then displays the board.
static void DisplayGoL(vector<vector<int>> myBoard, int n)
{
	int rank, p;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	int totalCells = n * n;
	int columns = n / p;
	clock_t GatherTime;

	// First initialize the entire board
	vector<int> wholeBoard;
	wholeBoard.resize(totalCells);

	// We need to copy myBoard from a 2d vector
	// to a 1d vector to make it MPI_Send friendly
	vector<int> SendableBoard;
	int CellsToSend = n * columns;
	SendableBoard.resize(CellsToSend);

	// Loop through myBoard and add each element
	// to SendableBoard.
	int k = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			SendableBoard[k] = myBoard[j][i];
			k++;
		}
	}

	
	GatherTime = clock();
	// Do an MPI_Gather to assimilate the board
	MPI_Gather(SendableBoard.data(), CellsToSend, MPI_INT, wholeBoard.data(), CellsToSend, MPI_INT, 0, MPI_COMM_WORLD);
	GatherTime = clock() - GatherTime;
	if (rank == 0)
	{
		cout << "Gather Time: " << (float)GatherTime / 1000 << endl;
	}

	if (rank == 0)
	{
		for (int i = 0; i < totalCells; i++)
		{
			if (i > 0 && i % n == 0)
			{
				cout << endl;
			}
			cout << wholeBoard[i] << " ";
		}
		cout << endl << endl;
	}
	
	return;
}

// This function simulates the game for G generations
static void Simulate(vector<vector<int>> myBoard, int n, int G)
{
	int rank, p;
	int mySeed = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	int columns = n / p;
	clock_t NeighborColumnSendTime, NeighborColumnReceiveTime;
	float TotalCSTime = 0.0, TotalCRTime = 0.0;

	// Initialize timers
	clock_t DisplayTime, GenerationTime;
	float TotalDisplayTime = 0.0, TotalGenerationTime = 0.0;

	// Set up containers for the previous and next columns
	vector<int> prevColumn, nextColumn;
	prevColumn.resize(n);
	nextColumn.resize(n);

	// The first thing to do is to make sure that each cell 
	// has access to all of their neighbors.
	// First figure out where the previous column and next
	// column will come from.
	int previous = rank - 1;
	int next = rank + 1;
	if (previous == -1)
	{
		previous = p - 1;
	}
	if (next == p)
	{
		next = 0;
	}
	// Print out the original state of the board
	if (rank == 0)
	{
		cout << "Starting board." << endl;
	}
	DisplayGoL(myBoard, n);
	
	// Here is where we simulate the generations.
	for (int g = 1; g < G + 1; g++)
	{
		// Start the generation timer.
		GenerationTime = clock();

		if (p > 1)
		{
			// Send the first column to previous, and the last
			// column to next.
			NeighborColumnSendTime = clock();
			MPI_Send(myBoard[0].data(), n, MPI_INT, previous, 0, MPI_COMM_WORLD);
			MPI_Send(myBoard[columns - 1].data(), n, MPI_INT, next, 0, MPI_COMM_WORLD);
			NeighborColumnSendTime = clock() - NeighborColumnSendTime;
			TotalCSTime += (float)NeighborColumnSendTime / 1000;

			// Receive the necessary columns from the previous and
			// next ranks.
			MPI_Status status;
			NeighborColumnReceiveTime = clock();
			MPI_Recv(prevColumn.data(), n, MPI_INT, previous, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(nextColumn.data(), n, MPI_INT, next, 0, MPI_COMM_WORLD, &status);
			NeighborColumnReceiveTime = clock() - NeighborColumnReceiveTime;
			TotalCRTime += (float)NeighborColumnReceiveTime / 1000;
		}
		else if (p == 1)
		{
			prevColumn = myBoard[myBoard.size() - 1];
			nextColumn = myBoard[0];
		}
		

		// Print out the columns received.  Just for testing
		/*for (int i = 0; i < n; i++)
		{
		cout << prevColumn[i] << " ";
		} */

		vector<vector<int>> oldBoard = myBoard;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < columns; j++)
			{
				myBoard[j][i] = DetermineState(oldBoard, j, i, prevColumn, nextColumn, columns);
			}
		}
		// Reset Old Board to the current board, and move to the next generation
		oldBoard = myBoard;

		if (g % 10 == 0)
		{
			if (rank == 0)
			{
				cout << "Generation " << g << " completed." << endl;
			}
			DisplayTime = clock();
			DisplayGoL(myBoard, n);
			DisplayTime = clock() - DisplayTime;
			TotalDisplayTime += (float)DisplayTime / 1000;
		}
		// Do a MPI_Barrier to ensure all processes are on
		// the same generation.
		MPI_Barrier(MPI_COMM_WORLD);

		// End the generation timer
		GenerationTime = clock() - GenerationTime;
		TotalGenerationTime += (float)GenerationTime / 1000;
	}
	float AverageDisplayTime = TotalDisplayTime / (G / 10);
	cout << "Average Display Time: " << AverageDisplayTime << endl;

	float AverageGenerationTime = TotalGenerationTime / G;
	cout << "Average Generation Time: " << AverageGenerationTime << endl;

	if (p > 1)
	{
		cout << "Total Neighbor Column Send Time: " << TotalCSTime << endl;

		cout << "Total Neighbor Column Receive Time: " << TotalCRTime << endl;
	}
	


	return;
}