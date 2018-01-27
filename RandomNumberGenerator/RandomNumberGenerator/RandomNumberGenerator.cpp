#include "functions.h"

int main(int argc, char *argv[])
{
	// Establishes what rank it is, and how many processes are running.
	static int rank, p, n, per_Process;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	static vector<int> Broadcast_data;
	n = 1048576;
	per_Process = n / p;

	// Initialize timers
	clock_t TotalTime, BroadcastTime, ParallelPrefixTime;
	float t_time = 0.0, b_time = 0.0, p_time = 0.0;
	TotalTime = clock();

	// The first and second arguments are constants for number generation, the third is a large prime to mod by, and the fourth is a random seed. x1 is calculated based off x0.
	// All provided by the user except x1.
	// Rank 0 broadcasts the data to all processes.
	if (rank == 0)
	{
		for (static int i = 1; i < 5; i++)
		{
			Broadcast_data.push_back(std::atoi(argv[i]));
		}
		Broadcast_data.push_back(std::atoi(argv[1]) *std::atoi(argv[4]) % std::atoi(argv[3]));

		// NOTE: THIS PUSH BACK IS HOW MANY RANDOM NUMBERS WILL BE GENERATED
		Broadcast_data.push_back(n);
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
	BroadcastTime = clock();
	MPI_Bcast(Broadcast_data.data(), 6, MPI_INT, 0, MPI_COMM_WORLD);
	BroadcastTime = clock() - BroadcastTime;
	b_time += (float)BroadcastTime / 1000;
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
	ParallelPrefixTime = clock();
	prefix = parallel_prefix(M_values, rank, p, Broadcast_data[3]);
	ParallelPrefixTime = clock() - ParallelPrefixTime;
	p_time += (float)ParallelPrefixTime / 1000;


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
		temp_vector[0][0] = local_sum[0][0] * M_values[i][0][0] + local_sum[1][0] * M_values[i][0][1];
		temp_vector[0][1] = local_sum[0][1] * M_values[i][0][0] + local_sum[1][1] * M_values[i][0][1];
		temp_vector[1][0] = local_sum[0][0] * M_values[i][1][0] + local_sum[0][1] * M_values[i][1][1];
		temp_vector[1][1] = local_sum[0][1] * M_values[i][1][0] + local_sum[1][1] * M_values[i][1][1];

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
		temp_vector[0][0] = prefix[0][0] * M_values[i][0][0] + prefix[1][0] * M_values[i][0][1];
		temp_vector[0][1] = prefix[0][1] * M_values[i][0][0] + prefix[1][1] * M_values[i][0][1];
		temp_vector[1][0] = prefix[0][0] * M_values[i][1][0] + prefix[0][1] * M_values[i][1][1];
		temp_vector[1][1] = prefix[0][1] * M_values[i][1][0] + prefix[1][1] * M_values[i][1][1];

		M_values[i] = temp_vector;
	}
	
	// Now all we have to do are local operations to generate the random numbers
	// [xi xi-1] = [x1 x0]M^i-1
	// So xi = ax1 + cx0
	int x0 = std::atoi(argv[4]);
	int x1 = std::atoi(argv[1]) *x0 % std::atoi(argv[3]);
	vector<int> random_numbers;
	for (int i = 0; i < per_Process; i++)
	{
		int xi;
		xi = x1 * M_values[i][0][0] + x0 * M_values[i][0][1];
		xi = xi % std::atoi(argv[3]);
		random_numbers.push_back(xi);
	}
	/*cout << "Rank " << rank << " Random numbers: ";
	for (int i = 0; i < per_Process; i++)
	{
		cout << random_numbers[i] << " ";
	}
	cout << endl;*/

	TotalTime = clock() - TotalTime;
	t_time += (float)TotalTime / 1000;

	cout << "Rank " << rank << " Broadcast time: " << b_time << endl;
	cout << "Rank " << rank << " Total time: " << t_time << endl;
	cout << "Rank " << rank << " Total parallel prefix time: " << p_time << endl;

	MPI_Finalize();

	return 0;
}