

#ifdef MPI_PROFILE

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <vector>	// std::vector
#include <mpi.h>
#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "time.h"

#include "mpi_bench.h"

static int rank;
static std::stringstream ss;

static uint64_t mpiStartTime, mpiEndTime;

#define ALL_RANKS	-1

#define isMasterQ       ((rank) == 0)

#ifdef __APPLE__
#include    <mach/mach_time.h>
static double o_timebase = 0;
static uint64_t o_timestart = 0;
#endif /* __APPLE__ */

EXTERN_C uint64_t PMPI_hrtime(void) {
#define NANOSEC ((uint64_t) 1e9)
    struct timespec ts;
#ifdef __APPLE__
#define O_NANOSEC   (+1.0E-9)
#define O_GIGA      UINT64_C(1000000000)
    if (!o_timestart) {
        mach_timebase_info_data_t tb = { 0 };
        mach_timebase_info(&tb);
        o_timebase = tb.numer;
        o_timebase /= tb.denom;
        o_timestart = mach_absolute_time();
    }
    double diff = (mach_absolute_time() - o_timestart) * o_timebase;
    ts.tv_sec = diff * O_NANOSEC;
    ts.tv_nsec = diff - (ts.tv_sec * O_GIGA);
#undef O_NANOSEC
#undef O_GIGA
#else /* __APPLE__ */
    clock_gettime(CLOCK_MONOTONIC, &ts);
#endif /* __APPLE__ */
    return (((uint64_t) ts.tv_sec) * NANOSEC + ts.tv_nsec);
#undef NANOSEC
}

EXTERN_C int PMPI_getRank() {
	return rank;
}

static inline int rankCount() {
    int nRanks;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    return nRanks;
}

static inline size_t calcByteCount(int count, MPI_Datatype dt) {
    MPI_Aint sz;
    MPI_Type_extent(dt, &sz);
    return count*sz;
}

EXTERN_C void PMPI_addTime(const char * callName, int line, size_t byteCount, uint64_t start, uint64_t end, int from, int to) {
    ss << callName << "," << rank << "," << line << "," << byteCount << "," << start << "," << end << "," << from << "," << to << std::endl;
    return ;
}

static inline void doPrint() {
    const std::string tmp = ss.str();
    const char* cstr = tmp.c_str();
    if (isMasterQ) {
        printf("%s", cstr);
    } else {
        // we are going to send the string to the master
        int len = (int) tmp.length() + 1;
        MPI_Send(&len, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send((void *) cstr, len, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
    }

    return ;
}

static inline void doInit() {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (isMasterQ) {
        ss << "XXX MPI_TIMER_START XXX" << std::endl;
        ss << "call,rank,line,bytes,start,end,from,to" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    mpiStartTime = PMPI_hrtime();
    return ;
}

EXTERN_C int pMPI_Init(int line, int *argc, char ***argv) {
    int err = MPI_Init(argc, argv);
    doInit();
    return err;
}

EXTERN_C int pMPI_Init_thread(int line, int *argc, char ***argv, int required, int *provided) {
    int err = MPI_Init_thread(argc, argv, required, provided);
    doInit();
    return err;
}

EXTERN_C int pMPI_Finalize(int line) {
    mpiEndTime = PMPI_hrtime();
    PMPI_addTime("MPITime", 0, 0, mpiStartTime, mpiEndTime, -1, -1);
    doPrint();

    if (isMasterQ) {
        int ii, nranks;
        nranks = rankCount();
        for (int ii = 1; ii < nranks; ii++) {
            char * buf;
            int bufSize;
            MPI_Recv(&bufSize, 1, MPI_INT, ii, ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            buf = (char *) calloc(bufSize, sizeof(char));
            MPI_Recv(buf, bufSize, MPI_CHAR, ii, ii, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s", buf);
            free(buf);
        }
        printf("XXX MPI_TIMER_END XXX\n");
    }
    int err = MPI_Finalize();
    return err;
}

EXTERN_C int pMPI_Bcast(int line, void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    int err;
    size_t byteCount = calcByteCount(count, datatype);
    detailedTimeBlock("MPI_Bcast", line, byteCount, rank, ALL_RANKS, {
        err = MPI_Bcast(buffer, count, datatype, root, comm);
    });
    return err;
}

EXTERN_C int pMPI_Send(int line, void *buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    int err;
    size_t byteCount = calcByteCount(count, datatype);
    detailedTimeBlock("MPI_Send", line, byteCount, rank, dest, {
        err = MPI_Send(buffer, count, datatype, dest, tag, comm);
    });
    return err;
}

EXTERN_C int pMPI_Recv(int line, void *buffer, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Status *status) {
    int err;
    size_t byteCount = calcByteCount(count, datatype);
    detailedTimeBlock("MPI_Recv", line, byteCount, source, rank, {
        err = MPI_Recv(buffer, count, datatype, source, tag, comm, status);
    });
    return err;
}

EXTERN_C int pMPI_Irecv(int line, void *buffer, int count, MPI_Datatype datatype,
                        int source, int tag, MPI_Comm comm, MPI_Request *request) {
    int err;
    size_t byteCount = calcByteCount(count, datatype);
    detailedTimeBlock("MPI_Irecv", line, byteCount, source, rank, {
        err = MPI_Irecv(buffer, count, datatype, source, tag, comm, request);
    });
    return err;
}


EXTERN_C int pMPI_Isend(int line, void *buffer, int count, MPI_Datatype datatype,
                        int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    int err;
    size_t byteCount = calcByteCount(count, datatype);
    detailedTimeBlock("MPI_Isend", line, byteCount, rank, dest, {
        err = MPI_Isend(buffer, count, datatype, dest, tag, comm, request);
    });
    return err;
}

EXTERN_C int pMPI_Waitall(int line, int count, MPI_Request array_of_requests[],
                          MPI_Status array_of_statuses[]) {
    int err;
    detailedTimeBlock("MPI_Waitall", line, 0, -1, -1, {
        err = MPI_Waitall(count, array_of_requests, array_of_statuses);
    });
    return err;
}

EXTERN_C int pMPI_Scatter(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                          void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
                          MPI_Comm comm) {
    int to;
    int err;
    size_t byteCount;
    if (isMasterQ) {
        int nranks = rankCount();
        to = ALL_RANKS;
        byteCount = nranks * calcByteCount(sendcnt, sendtype);
    } else {
        to = root;
        byteCount = calcByteCount(recvcnt, recvtype);
    }
    detailedTimeBlock("MPI_Scatter", line, byteCount, rank, to, {
        err = MPI_Scatter(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
    });
    return err;
}

EXTERN_C int pMPI_Gather(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                         void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                         int root, MPI_Comm comm) {
    int to;
    int err;
    size_t byteCount;
    if (isMasterQ) {
        int nranks = rankCount();
        to = ALL_RANKS;
        byteCount = nranks*calcByteCount(recvcnt, recvtype);
    } else {
        to = root;
        byteCount = calcByteCount(sendcnt, sendtype);
    }
    detailedTimeBlock("MPI_Gather", line, byteCount, rank, to, {
        err = MPI_Gather(sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm);
    });
    return err;
}

EXTERN_C int pMPI_Reduce(int line, void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                         MPI_Op op, int root, MPI_Comm comm) {
    int to;
    int err;
    int byteCount = calcByteCount(count, datatype);
    if (isMasterQ) {
        to = ALL_RANKS;
    } else {
        to = root;
    }
    detailedTimeBlock("MPI_Reduce", line, byteCount, rank, to, {
        err = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    });
    return err;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


static std::vector<MPI_Request> requests;

static inline void addRequest(MPI_Request req) {
    requests.push_back(req);
    return ;
}

static inline int doWaitForRequests(int line) {
    int err = MPI_SUCCESS;
    int numRequests = requests.size();
    if (numRequests > 0) {
        MPI_Request * reqs = requests.data();
        err = pMPI_Waitall(line, numRequests, reqs, MPI_STATUSES_IGNORE);
    }
    requests.clear();
    return err;
}

static inline int doWaitForRequests() {
    return doWaitForRequests(-1);
}

EXTERN_C int iMPI_RequestWait(int line) {
    return doWaitForRequests(line);
}

EXTERN_C int iMPI_Init(int line, int *argc, char ***argv) {
    requests.clear();
    return pMPI_Init(line, argc, argv);
}

EXTERN_C int iMPI_Init_thread(int line, int *argc, char ***argv, int required, int *provided) {
    requests.clear();
    return pMPI_Init_thread(line, argc, argv, required, provided);
}

EXTERN_C int iMPI_Finalize(int line) {
    doWaitForRequests();
    return pMPI_Finalize(line);
}

EXTERN_C int iMPI_Send(int line, void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    MPI_Request req;
    int err = pMPI_Isend(line, buf, count, datatype, dest, tag, comm, &req);
    addRequest(req);
    return err;
}

EXTERN_C int iMPI_Recv(int line, void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    MPI_Request req;
    int err = pMPI_Irecv(line, buf, count, datatype, source, tag, comm, &req);
    addRequest(req);
    return err;
}

EXTERN_C int iMPI_Bcast(int line, void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
    if (isMasterQ) {
        int ii;
        int err;
        int nRanks = rankCount();
        for (ii = 1; ii < nRanks; ii++) {
            err = iMPI_Send(line, buffer, count, datatype, ii, ii, comm);
            if (err != MPI_SUCCESS) {
                break ;
            }
        }
        return err;
    } else {
        MPI_Status st;
        return iMPI_Recv(line, buffer, count, datatype, root, rank, comm, &st);
    }
}

template <typename T>
static inline T * doReduce(MPI_Op op, T * out, const T * in, int count) {
    int ii, jj;
    int nRanks;
    const T * iter;

    assert(isMasterQ);
    assert(out != NULL);

    nRanks = rankCount();
    if (op == MPI_MAX) { // cannot use a switch here...
        for (jj = 1; jj < nRanks; jj++) {
            iter = &in[jj*count];
            for (ii = 1; ii < count; ii++) {
                if (out[ii] < iter[ii]) {
                    out[ii] = iter[ii];
                }
            }
        }
    } else if(op == MPI_MIN) {
        for (jj = 1; jj < nRanks; jj++) {
            iter = &in[jj*count];
            for (ii = 1; ii < count; ii++) {
                if (out[ii] > iter[ii]) {
                    out[ii] = iter[ii];
                }
            }
        }
    } else if(op == MPI_SUM) {
        for (jj = 1; jj < nRanks; jj++) {
            iter = &in[jj*count];
            for (ii = 1; ii < count; ii++) {
                out[ii] += iter[ii];
            }
        }
    } else if(op == MPI_PROD) {
        for (jj = 1; jj < nRanks; jj++) {
            iter = &in[jj*count];
            for (ii = 1; ii < count; ii++) {
                out[ii] *= iter[ii];
            }
        }
    } else {
        printf("::ERROR:: Unsupported operation...\n");
        return NULL;
    }
    return out;
}

static inline void * doReduce(MPI_Datatype datatype, MPI_Op op, void * out, const void * in, int count) {
    if (datatype == MPI_INT) { // cannot use a switch here...
        return doReduce<int>(op, (int *) out, (const int *) in, count);
    } else if (datatype == MPI_FLOAT) {
        return doReduce<float>(op, (float *) out, (const float *) in, count);
    } else if (datatype == MPI_UNSIGNED_LONG_LONG) {
        return doReduce<unsigned long long>(op, (unsigned long long *) out, (const unsigned long long *) in, count);
    } else {
        printf("::ERROR:: Unsupported type...\n");
    }
    return NULL;
}

EXTERN_C int iMPI_Reduce(int line, void *sendbuf, void *recvbuf,
                         int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    if (recvbuf == NULL) {
        return iMPI_Send(line, sendbuf, count, datatype, root, rank, comm);
    } else {
        int ii;
        int err;
        void * red;
        int nRanks = rankCount();
        size_t stride = calcByteCount(count, datatype);
        size_t byteCount = nRanks*stride;
        char * buf = (char *) malloc(byteCount);
        for (ii = 1; ii < nRanks; ii++) {
            MPI_Status st;
            err = iMPI_Recv(line, buf + ii * stride, count, datatype, ii, ii, comm, &st);
            if (err != MPI_SUCCESS) {
                return err;
            }
        }
        memcpy(recvbuf, sendbuf, stride);
        doWaitForRequests(line);
        red = doReduce(datatype, op, recvbuf, buf, count);
        assert(red != NULL);
        free(buf);
        return err;
    }
}

EXTERN_C int iMPI_Gather(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                         void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                         int root, MPI_Comm comm) {
    if (recvbuf == NULL) {
        MPI_Status st;
        assert(!isMasterQ);
        return iMPI_Send(line, sendbuf, sendcnt, sendtype, root, rank, comm);
    } else {
        int ii;
        int err = MPI_SUCCESS;

        int nRanks = rankCount();
        size_t sendStride = calcByteCount(sendcnt, sendtype);
        size_t recvStride = calcByteCount(sendcnt, sendtype);
        if (sendStride > recvStride) {
            memcpy(recvbuf, sendbuf, recvStride);
        } else {
            memcpy(recvbuf, sendbuf, sendStride);
        }
        for (ii = 1; ii < nRanks; ii++) {
	    MPI_Status st;
            void *recv_subarray = (char *)recvbuf + ii*recvStride;
       	    err = iMPI_Recv(line, recv_subarray, recvcnt, recvtype, ii, ii, comm, &st);
            if (err != MPI_SUCCESS) {
                return err;
            }
        }
        return err;
    }
}

EXTERN_C int iMPI_Scatter(int line, void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                          void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root,
                          MPI_Comm comm) {
    if (sendbuf == NULL) {
        MPI_Status st;
        assert(!isMasterQ);
        return iMPI_Recv(line, recvbuf, recvcnt, recvtype, root, rank, comm, &st);
    } else {
        int ii;
        int err = MPI_SUCCESS;

        int nRanks = rankCount();
        size_t sendStride = calcByteCount(sendcnt, sendtype);
        size_t recvStride = calcByteCount(sendcnt, sendtype);
        if (sendStride > recvStride) {
            memcpy(recvbuf, sendbuf, recvStride);
        } else {
            memcpy(recvbuf, sendbuf, sendStride);
        }
        for (ii = 1; ii < nRanks; ii++) {
            void *send_subarray = (char *)sendbuf + ii*sendStride;
            err = iMPI_Send(line, send_subarray, sendcnt, sendtype, ii, ii, comm);
            if (err != MPI_SUCCESS) {
                return err;
            }
        }
        return err;
    }
}


#endif /* MPI_PROFILE */

