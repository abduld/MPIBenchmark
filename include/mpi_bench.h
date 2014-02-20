
#ifndef __MPI_BENCH_H__
#define __MPI_BENCH_H__

#ifdef MPI_PROFILE_DEBUG
#define MPI_DEBUG(...)    printf("::DEBUG:: " __VA_ARGS__)
#else
#define MPI_DEBUG(...)
#endif

#ifdef MPI_PROFILE

#include "stdint.h"

#ifdef __cplusplus
#define EXTERN_C extern "C"
#define BEGIN_EXTERN_C EXTERN_C {
#define END_EXTERN_C }
#else
#define EXTERN_C
#define BEGIN_EXTERN_C
#define END_EXTERN_C
#endif /* __cplusplus */


#define IMPI_hrtime PMPI_hrtime
#define PMPI_RequestWait()
#define PMPI_Init(...)      pMPI_Init(__LINE__, __VA_ARGS__)
#define PMPI_Init_thread(...) pMPI_Init_thread(__LINE__, __VA_ARGS__)
#define PMPI_Bcast(...)     pMPI_Bcast(__LINE__, __VA_ARGS__)
#define PMPI_Send(...)      pMPI_Send(__LINE__, __VA_ARGS__)
#define PMPI_Recv(...)      pMPI_Recv(__LINE__, __VA_ARGS__)
#define PMPI_Scatter(...)   pMPI_Scatter(__LINE__, __VA_ARGS__)
#define PMPI_Gather(...)    pMPI_Gather(__LINE__, __VA_ARGS__)
#define PMPI_Reduce(...)    pMPI_Reduce(__LINE__, __VA_ARGS__)
#define PMPI_Isend(...)     pMPI_Isend(__LINE__, __VA_ARGS__)
#define PMPI_Irecv(...)     pMPI_Irecv(__LINE__, __VA_ARGS__)
#define PMPI_Waitall(...)   pMPI_Waitall(__LINE__, __VA_ARGS__)
#define PMPI_Finalize()     pMPI_Finalize(__LINE__)


#define IMPI_RequestWait() iMPI_RequestWait(__LINE__)
#define IMPI_Init(...)      iMPI_Init(__LINE__, __VA_ARGS__)
#define IMPI_Init_thread(...) iMPI_Init_thread(__LINE__, __VA_ARGS__)
#define IMPI_Bcast(...)     iMPI_Bcast(__LINE__, __VA_ARGS__)
#define IMPI_Send(...)      iMPI_Send(__LINE__, __VA_ARGS__)
#define IMPI_Recv(...)      iMPI_Recv(__LINE__, __VA_ARGS__)
#define IMPI_Scatter(...)   iMPI_Scatter(__LINE__, __VA_ARGS__)
#define IMPI_Gather(...)    iMPI_Gather(__LINE__, __VA_ARGS__)
#define IMPI_Reduce(...)    iMPI_Reduce(__LINE__, __VA_ARGS__)
#define IMPI_Isend(...)     iMPI_Isend(__LINE__, __VA_ARGS__)
#define IMPI_Irecv(...)     iMPI_Irecv(__LINE__, __VA_ARGS__)
#define IMPI_Waitall(...)   iMPI_Waitall(__LINE__, __VA_ARGS__)
#define IMPI_Finalize()     iMPI_Finalize(__LINE__)

BEGIN_EXTERN_C
    uint64_t PMPI_hrtime(void);
    int PMPI_getRank(void);
    void PMPI_addTime(const char * callName, int line, size_t byteCount, uint64_t start, uint64_t end, int from, int to);
    int pMPI_Init(int line, int *argc, char ***argv);
    int pMPI_Init_thread(int line, int *argc, char ***argv, int required, int *provided);
    int pMPI_Finalize(int line);
    int pMPI_Bcast(int line, void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
    int pMPI_Send(int line, void *buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    int pMPI_Recv(int line, void *buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
    int pMPI_Irecv(int line, void *buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
    int pMPI_Isend(int line, void *buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
    int pMPI_Waitall(int line, int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
    int pMPI_Scatter(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int pMPI_Gather(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int pMPI_Reduce(int line, void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

    int iMPI_RequestWait(int line);
    int iMPI_Init(int line, int *argc, char ***argv);
    int iMPI_Init_thread(int line, int *argc, char ***argv, int required, int *provided);
    int iMPI_Finalize(int line);
    int iMPI_Bcast(int line, void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
    int iMPI_Send(int line, void *buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
    int iMPI_Recv(int line, void *buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
#define iMPI_Irecv pMPI_Irecv
#define iMPI_Isend pMPI_Isend
#define iMPI_Waitall pMPI_WaitAll
    int iMPI_Scatter(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int iMPI_Gather(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm);
    int iMPI_Reduce(int line, void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
END_EXTERN_C


#define startTimeBlock(tic, blk)  do {                              \
                                     tic = PMPI_hrtime();           \
                                     blk;                           \
                                  } while(0)

#define endTimeBlock(tic, toc, msg, line, blk) do {                         \
                                            blk;                            \
                                            toc = PMPI_hrtime();            \
                                            PMPI_addTime(msg, line, 0, tic, toc, -1, -1);\
                                         } while(0)

#define detailedTimeBlock(msg, line, byteCount, from, to, blk) do {           \
                                    uint64_t tic, toc;              \
                                    tic = PMPI_hrtime();		    \
				    MPI_DEBUG("MASTER(%d) ::: (%d) before processing %s\n", line, from, msg);	\
                                    blk;		  	                \
                                    toc = PMPI_hrtime();		    \
				    MPI_DEBUG("MASTER(%d) ::: (%d) after processing %s\n", line, from, msg);	\
                                    PMPI_addTime(msg, line, byteCount, tic, toc, from, to);\
                                } while(0)

#define timeBlock(msg, blk)     do {                                \
				    detailedTimeBlock(msg, __LINE__, 0, -1, -1, blk); \
                                } while(0)


#define recordTime(msg)         PMPI_addTime(msg, __LINE__, 0, PMPI_hrtime(), 0, -1, -1)

#else /* MPI_PROFILE */

#define PMPI_RequestWait()
#define PMPI_Init      MPI_Init
#define PMPI_Init_thread MPI_Init_thread
#define PMPI_Bcast     MPI_Bcast
#define PMPI_Send      MPI_Send
#define PMPI_Recv      MPI_Recv
#define PMPI_Scatter   MPI_Scatter
#define PMPI_Reduce    MPI_Reduce
#define PMPI_Gather    MPI_Gather
#define PMPI_Isend     MPI_Isend
#define PMPI_Irecv     MPI_Irecv
#define PMPI_WaitAll   MPI_WaitAll
#define PMPI_Finalize  MPI_Finalize

#define IMPI_RequestWait()
#define IMPI_Init      MPI_Init
#define IMPI_Init_thread MPI_Init_thread
#define IMPI_Bcast     MPI_Bcast
#define IMPI_Send      MPI_Send
#define IMPI_Recv      MPI_Recv
#define IMPI_Scatter   MPI_Scatter
#define IMPI_Reduce    MPI_Reduce
#define IMPI_Gather    MPI_Gather
#define IMPI_Isend     MPI_Isend
#define IMPI_Irecv     MPI_Irecv
#define IMPI_WaitAll   MPI_WaitAll
#define IMPI_Finalize  MPI_Finalize


#define startTimeBlock(tic, blk)  blk
#define endTimeBlock(tic, toc, msg, line, blk) blk
#define detailedTimeBlock(msg, line, byteCount, from, to, blk) blk
#define timeBlock(msg, blk) blk
#define recordTime(msg)


#endif /* MPI_PROFILE */

#endif /* __MPI_BENCH_H__ */



