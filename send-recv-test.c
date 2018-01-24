/*
This program measures how long it takes MPI_Send and MPI_Recv to work.
This program loops through message sizes from 1 byte to 4096 bytes, sends each message 1000 times, averages them, and outputs it to the screen.
*/


#include <stdio.h>
#include <mpi.h>
#include <assert.h>
#include <sys/time.h>

int main(int argc,char *argv[])
{

   int rank,p;
   int i=1, j=0;
   struct timeval t1,t2;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   MPI_Comm_size(MPI_COMM_WORLD,&p);

   assert(p>=2);

   while(i<4097)
   {
        int tSendTotal=0, tRecvTotal=0, tSendAvg=0, tRecvAvg=0;
        while(j<1000)
        {
                 if(rank==7) {
                                char x[4097] = { 'A' };
                                int dest = 0, k = 0;
                                while(k<i)
                                {
                                        x[k] = 'A';
                                        k++;
                                }
                                gettimeofday(&t1,NULL);
                                MPI_Send(&x[0],i,MPI_CHAR,dest,0,MPI_COMM_WORLD);
                                gettimeofday(&t2,NULL);
                                tSendTotal += (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                              } else
                 if (rank==0) {
                                char y[4097]={ 0 };
                                MPI_Status status;
                                gettimeofday(&t1,NULL);
                                MPI_Recv(&y[0],i,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
                                gettimeofday(&t2,NULL);
                                tRecvTotal += (t2.tv_sec-t1.tv_sec)*1000000 + (t2.tv_usec-t1.tv_usec);
                              }
           j++;
         }
   j=0;
   tSendAvg = tSendTotal/1000;
   tRecvAvg = tRecvTotal/1000;
   if (tSendAvg > 0)
   {
   printf("Rank=7: send message %d to rank 0; Send time %d microsec\n", i, tSendAvg);
   }
   if (tRecvAvg > 0)
   {
   printf("Rank=0: received message %d from rank 0; Recv time %d microsec\n", i, tRecvAvg);
   }
   i=i*2;
   sleep(1);
   }

   MPI_Finalize();
}
