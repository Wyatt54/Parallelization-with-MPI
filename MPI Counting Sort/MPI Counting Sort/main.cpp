#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

vector<vector<int>> parallel_prefix(vector<vector<vector<int> > >, int, int, int);
vector<int> MPI_RNG(int);

int main(int argc, char *argv[])
{
	// Establishes what rank it is, and how many processes are running.
	int rank, p, n, max, perProcess, size;
	std::vector<int> toSort, bigSort, counts, allCounts, sorted, finalSorted, sizes;
	clock_t totalTime, sortTime, reduceTime, countTime;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL));

	totalTime = clock();
	if (argc > 1)
	{
		max = std::atoi(argv[1]);
		n = std::atoi(argv[2]);
	}
	else
	{
		max = 1000;
		n = 8388608;
	}
	if (rank == 0)
	{
		std::cout << "The range is (0, " << max << ").  " << n << " numbers to be sorted." << std::endl;
	}
	perProcess = n / p;

	toSort.resize(perProcess);
	toSort = MPI_RNG(n);
	MPI_Barrier(MPI_COMM_WORLD);

	// Now we distribute, and begin the timer for sorting.
	sortTime = clock();
	
	// Now every process has a portion of the array to sort. Next step
	// is to initialize the counting array.
	// Initialize the counting arrays to size of max.
	counts.resize(max);
	allCounts.resize(max);
	finalSorted.resize(n);
	sizes.resize(p);

	
	sortTime = clock();
	countTime = clock();
	// Here is the sorting part. Go through the array toSort and increment
	// the count for each encounter.
	for (int i = 0; i < toSort.size(); i++)
	{
		counts[toSort[i]]++;
	}
	countTime = clock() - countTime;

	// Now we do an all reduce so all processes have all the counts.
	reduceTime = clock();
	MPI_Allreduce(counts.data(), allCounts.data(), max, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	reduceTime = clock() - reduceTime;

	perProcess = max / p;
	int start = rank * perProcess;
	int end = start + perProcess;
	for (start = rank * perProcess; start < end; start++)
	{
		for (int j = 0; j < allCounts[start]; j++)
		{
			sorted.push_back(start);
		}
	}

	// Now we have to gather all the arrays in rank 0.  The problem here is that the 
	// messages are of variable size.  Therefore we must first send the sizes, and then
	// recieve the messages.
	size = sorted.size();
	MPI_Gather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Now we can send the sorted arrays
	if (rank > 0)
	{
		MPI_Send(sorted.data(), size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	
	// Now rank 0 can receive the messages.
	if (rank == 0)
	{
		MPI_Status status;
		finalSorted = sorted;
		finalSorted.resize(n);
		int sizeIncrement = sizes[0];
		for (int i = 1; i < p; i++)
		{
			MPI_Recv(finalSorted.data() + sizeIncrement, sizes[i], MPI_INT, i, 0, MPI_COMM_WORLD, &status);
			sizeIncrement += sizes[i];
		}
		cout << "Array is sorted." << endl;
		totalTime = clock() - totalTime;
		sortTime = clock() - sortTime;
		std::cout << "Total time taken: " << (float)totalTime / 1000 << " seconds." << std::endl;
		std::cout << "Time to sort: " << (float)sortTime / 1000 << " seconds." << std::endl;
	}
	

	MPI_Finalize();
	return 0;
}

// This is parallel prefix with the operator being matrix multiplication
vector<vector<int>> parallel_prefix(vector<vector<vector<int>>> Matrices, int rank, int p, int mod_by)
{
	int xor_by = 1;
	// The first step is a local multiplication of all M values.
	// In a matrix represented by:
	// [ a b ]
	// [ c d ]
	// The new matrix will be this:
	// [ a^2+bc ab+bd ]
	// [ ca+dc cb+d^2 ]
	// So the first step will be to complete this operation once for every matrix M in M_values

	static vector<vector<int>> global_sum;
	global_sum = Matrices[0];
	for (static int i = 1; i < Matrices.size(); i++)
	{
		vector<vector<int>> temp_vector;
		temp_vector.resize(2);
		temp_vector[0].resize(2);
		temp_vector[1].resize(2);
		temp_vector = global_sum;
		temp_vector[0][0] = (global_sum[0][0] * Matrices[i][0][0] + global_sum[1][0] * Matrices[i][0][1]) % 1000;
		temp_vector[0][1] = (global_sum[0][1] * Matrices[i][0][0] + global_sum[1][1] * Matrices[i][0][1]) % 1000;
		temp_vector[1][0] = (global_sum[0][0] * Matrices[i][1][0] + global_sum[0][1] * Matrices[i][1][1]) % 1000;
		temp_vector[1][1] = (global_sum[0][1] * Matrices[i][1][0] + global_sum[1][1] * Matrices[i][1][1]) % 1000;

		global_sum = temp_vector;
	}
	// Initialize the local sum as the identity matrix
	static vector<vector<int>> local_sum;
	local_sum = global_sum;
	local_sum[0][0] = 1;
	local_sum[0][1] = 0;
	local_sum[1][0] = 0;
	local_sum[1][1] = 1;
	//cout << "Send Vector: [ " << local_sum[0][0] << " " << local_sum[0][1] << " " << local_sum[1][0] << " " << local_sum[1][1] << " " << endl;
	// Now that all the local sums have been computed we can start step 2: communication.

	// Determine how many steps it will take
	int steps = 0;
	int j = 1;
	while (j < p)
	{
		j *= 2;
		steps++;
	}
	int k = 0;
	while (k < steps)
	{
		// First determine the rank's mate.
		static int mate;
		mate = rank ^ (xor_by << k);
		//cout << "Rank " << rank << " mate is: " << mate << endl;

		// Now we send the global sum to mate, and receive our mate's global sum.
		// First modify the local sum vector to a vector that can be sent.
		// Send vector syntax is [ a c b d ]
		static vector<int> send_vector, recv_vector;
		send_vector.resize(4);
		recv_vector.resize(4);
		send_vector[0] = global_sum[0][0];
		send_vector[1] = global_sum[0][1];
		send_vector[2] = global_sum[1][0];
		send_vector[3] = global_sum[1][1];
		//cout << "Rank " << rank << " Send Vector: [ " << send_vector[0] << " " << send_vector[1] << " " << send_vector[2] << " " << send_vector[3] << " " << endl;

		// Send the vector to your mate, and recieve a vector from your mate.
		static MPI_Status status;

		// Update the local sum if your mate rank is lower than your rank.
		if (mate < rank)
		{
			// If mate is less than rank we first receive, then we send.
			MPI_Recv(recv_vector.data(), 4, MPI_INT, mate, 0, MPI_COMM_WORLD, &status);
			MPI_Send(send_vector.data(), 4, MPI_INT, mate, 1, MPI_COMM_WORLD);
			static vector<vector<int>> temp_vector;
			temp_vector = local_sum;
			temp_vector[0][0] = (local_sum[0][0] * recv_vector[0] + local_sum[1][0] * recv_vector[1]) % 1000;
			temp_vector[0][1] = (local_sum[0][1] * recv_vector[0] + local_sum[1][1] * recv_vector[1]) % 1000;
			temp_vector[1][0] = (local_sum[0][0] * recv_vector[2] + local_sum[0][1] * recv_vector[3]) % 1000;
			temp_vector[1][1] = (local_sum[0][1] * recv_vector[2] + local_sum[1][1] * recv_vector[3]) % 1000;

			local_sum = temp_vector;
			//cout << "Rank " << rank << " New Local Sum [ " << local_sum[0][0] << " " << local_sum[0][1] << " " << local_sum[1][0] << " " << local_sum[1][1] << " " << endl;
		}
		else
		{
			// If mate > rank then we first send then receive.
			MPI_Send(send_vector.data(), 4, MPI_INT, mate, 0, MPI_COMM_WORLD);
			MPI_Recv(recv_vector.data(), 4, MPI_INT, mate, 1, MPI_COMM_WORLD, &status);
		}

		// No matter what each rank updates its global sum
		static vector<vector<int>> temp_vector;
		temp_vector = global_sum;
		temp_vector[0][0] = (global_sum[0][0] * recv_vector[0] + global_sum[1][0] * recv_vector[1]) % 1000;
		temp_vector[0][1] = (global_sum[0][1] * recv_vector[0] + global_sum[1][1] * recv_vector[1]) % 1000;
		temp_vector[1][0] = (global_sum[0][0] * recv_vector[2] + global_sum[0][1] * recv_vector[3]) % 1000;
		temp_vector[1][1] = (global_sum[0][1] * recv_vector[2] + global_sum[1][1] * recv_vector[3]) % 1000;

		global_sum = temp_vector;
		//cout << "Rank " << rank << " New Local Sum [ " << local_sum[0][0] << " " << local_sum[0][1] << " " << local_sum[1][0] << " " << local_sum[1][1] << " " << endl;
		//cout << "Rank " << rank << " New global Sum [ " << global_sum[0][0] << " " << global_sum[0][1] << " " << global_sum[1][0] << " " << global_sum[1][1] << " " << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		k++;
		// After completion of this loop the local sum is the parallel prefix output for each process.
	}
	//cout << "Rank: " << rank << " Local_Sum matrix: " << local_sum[0][0] << " " << local_sum[0][1] << " " << local_sum[1][0] << " " << local_sum[1][1] << endl;
	//cout << "Rank: " << rank << " Global_Sum matrix: " << global_sum[0][0] << " " << global_sum[0][1] << " " << global_sum[1][0] << " " << global_sum[1][1] << " " << endl;


	return local_sum;
}

vector<int> MPI_RNG(int size)
{
	// Establishes what rank it is, and how many processes are running.
	static int rank, p, per_Process;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	static vector<int> Broadcast_data;
	per_Process = size / p;


	// The first and second arguments are constants for number generation, the third is a large prime to mod by, and the fourth is a random seed. x1 is calculated based off x0.
	// All provided by the user except x1.
	// Rank 0 broadcasts the data to all processes.
	if (rank == 0)
	{
		Broadcast_data.push_back(34);
		Broadcast_data.push_back(865);
		Broadcast_data.push_back(1000);
		Broadcast_data.push_back(262);
		Broadcast_data.push_back((34 *262) % 1000);

		// NOTE: THIS PUSH BACK IS HOW MANY RANDOM NUMBERS WILL BE GENERATED
		Broadcast_data.push_back(size);
		cout << "Rank " << rank << " Broadcast Data: ";
		for (static int i = 0; i < 6; i++)
		{
			cout << Broadcast_data[i] << " ";
		}
		cout << endl;
	}
	else
	{
		Broadcast_data.resize(6);
	}
	MPI_Bcast(Broadcast_data.data(), 6, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	

	// Initialize an array of n/p values at every process.  Each of the n/p values is the matrix M.
	// M is this 2 dimmensional array:
	// [ a 1 ]
	// [ b 0 ]
	static vector<vector<int>> M;
	M.resize(2);
	M[0].resize(2);
	M[1].resize(2);
	M[0][0] = Broadcast_data[0];
	M[0][1] = Broadcast_data[1];
	M[1][0] = 1;
	M[1][1] = 0;

	// Now we must initialize the array of these M values.  Notation might get complex here
	// as we are dealing with 3D arrays.
	static vector<vector<vector<int>>> M_values;
	for (static int i = 0; i < per_Process; i++)
	{
		M_values.push_back(M);
	}
	// Now we are ready for the parallel prefix operation.  Note that the operator here
	// is matrix multiplication.
	static vector<vector<int>> prefix;
	prefix = parallel_prefix(M_values, rank, p, Broadcast_data[3]);

	// Now that we have the parallel prefix value, we need to multiply it by every matrix
	// in M values to get the C matrix from which we can get our random numbers.
	// But first we have to do a local parallel prefix on M_values so that each value
	// in M_values are different.

	static vector<vector<int>> local_sum;
	local_sum = M_values[0];
	for (static int i = 1; i < M_values.size(); i++)
	{
		vector<vector<int>> temp_vector;
		temp_vector.resize(2);
		temp_vector[0].resize(2);
		temp_vector[1].resize(2);
		temp_vector = local_sum;
		temp_vector[0][0] = (local_sum[0][0] * M_values[i][0][0] + local_sum[1][0] * M_values[i][0][1]) % 1000;
		temp_vector[0][1] = (local_sum[0][1] * M_values[i][0][0] + local_sum[1][1] * M_values[i][0][1]) % 1000;
		temp_vector[1][0] = (local_sum[0][0] * M_values[i][1][0] + local_sum[0][1] * M_values[i][1][1]) % 1000;
		temp_vector[1][1] = (local_sum[0][1] * M_values[i][1][0] + local_sum[1][1] * M_values[i][1][1]) % 1000;

		local_sum = temp_vector;
		M_values[i] = temp_vector;
	}

	for (int i = 0; i < M_values.size(); i++)
	{
		vector<vector<int>> temp_vector;
		temp_vector.resize(2);
		temp_vector[0].resize(2);
		temp_vector[1].resize(2);
		temp_vector = M_values[i];
		temp_vector[0][0] = (prefix[0][0] * M_values[i][0][0] + prefix[1][0] * M_values[i][0][1]) % 1000;
		temp_vector[0][1] = (prefix[0][1] * M_values[i][0][0] + prefix[1][1] * M_values[i][0][1]) % 1000;
		temp_vector[1][0] = (prefix[0][0] * M_values[i][1][0] + prefix[0][1] * M_values[i][1][1]) % 1000;
		temp_vector[1][1] = (prefix[0][1] * M_values[i][1][0] + prefix[1][1] * M_values[i][1][1]) % 1000;

		M_values[i] = temp_vector;
	}

	// Now all we have to do are local operations to generate the random numbers
	// [xi xi-1] = [x1 x0]M^i-1
	// So xi = ax1 + cx0
	int x0 = 262;
	int x1 = (34 *x0 )% 1000;
	vector<int> random_numbers;
	for (int i = 0; i < per_Process; i++)
	{
		int xi;
		xi = x1 * M_values[i][0][0] + x0 * M_values[i][0][1];
		xi = xi % 1000;
		random_numbers.push_back(xi);
	}



	return random_numbers;
}