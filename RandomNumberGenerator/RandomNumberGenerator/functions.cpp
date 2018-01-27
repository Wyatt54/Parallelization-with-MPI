#include "functions.h"

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
		temp_vector[0][0] = global_sum[0][0] * Matrices[i][0][0] + global_sum[1][0] * Matrices[i][0][1];
		temp_vector[0][1] = global_sum[0][1] * Matrices[i][0][0] + global_sum[1][1] * Matrices[i][0][1];
		temp_vector[1][0] = global_sum[0][0] * Matrices[i][1][0] + global_sum[0][1] * Matrices[i][1][1];
		temp_vector[1][1] = global_sum[0][1] * Matrices[i][1][0] + global_sum[1][1] * Matrices[i][1][1];

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
			temp_vector[0][0] = local_sum[0][0] * recv_vector[0] + local_sum[1][0] * recv_vector[1];
			temp_vector[0][1] = local_sum[0][1] * recv_vector[0] + local_sum[1][1] * recv_vector[1];
			temp_vector[1][0] = local_sum[0][0] * recv_vector[2] + local_sum[0][1] * recv_vector[3];
			temp_vector[1][1] = local_sum[0][1] * recv_vector[2] + local_sum[1][1] * recv_vector[3];

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
		temp_vector[0][0] = global_sum[0][0] * recv_vector[0] + global_sum[1][0] * recv_vector[1];
		temp_vector[0][1] = global_sum[0][1] * recv_vector[0] + global_sum[1][1] * recv_vector[1];
		temp_vector[1][0] = global_sum[0][0] * recv_vector[2] + global_sum[0][1] * recv_vector[3];
		temp_vector[1][1] = global_sum[0][1] * recv_vector[2] + global_sum[1][1] * recv_vector[3];

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