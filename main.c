#include <stdbool.h>
#include <time.h>
#include "mpi.h"
#include "matrix.h"

#define ERR_ARGV -1
#define ERR_MEMORY_ALLOCATION -2
#define ERR_FILE_OPEN -3
#define ERR_MPI_INIT -4

#define MAIN_PROC 0

typedef enum lu_print_mode {
    ALU_PRINT_MODE_L,
    ALU_PRINT_MODE_U,
    ALU_PRINT_MODE_LU
} lu_print_mode_te;

void print_alu_matrix(int dim, float **alu, lu_print_mode_te mode) {
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            float el;

            if (mode == ALU_PRINT_MODE_L) {
                el = i == j
                     ? 1
                     : j > i
                       ? 0
                       : alu[i][j];
            } else if (mode == ALU_PRINT_MODE_U) {
                el = j < i
                     ? 0
                     : alu[i][j];
            } else { // (mode == ALU_PRINT_MODE_LU)
                el = alu[i][j];
            }

            printf("%*.*f", 10, 2, el);
        }

        printf("\n");
    }
}

int fprint(char *filename, int nprocs, int dim, double time, char *name, int name_len) {
    FILE *file = fopen(filename, "a");

    if (file == NULL) {
        return ERR_FILE_OPEN;
    }

    printf("%15s -> filename %15s[%d] -> name %15d -> nprocs %15d -> dim %15.2f -> time\n", filename, name, name_len, nprocs, dim, time);
    fprintf(file, "%15s[%d] -> name %15d -> nprocs %15d -> dim %15.2f -> time\n", name, name_len, nprocs, dim, time);
    fclose(file);

    return 0;
}

/*
 * A = |  2    4   -4 |                          L = |   1   0   0 |
 *     |{ 1}  -4    3 | // R_1 - (1 / 2) * R_0       | 0.5   1   0 |
 *     | -6   -9    5 |                              |   X   X   1 |
 *
 * A = |  2    4   -4 |                           L = |   1   0   0 |
 *     |  0   -6    5 |                               | 0.5   1   0 |
 *     |{-6}  -9    5 | // R_2 - (-6 / 2) * R_0       |  -3   X   1 |
 *
 * A = |  2    4   -4 |                           L = |    1    0    0 |
 *     |  0   -6    5 |                               |  0.5    1    0 |
 *     |  0  { 3}  -7 | // R_2 - (3 / -6) * R_1       |   -3 -0.5    1 |
 *
 * U = |    2      4     -4 |   LU = |    2      4     -4 |
 *     |    0     -6      5 |        |  0.5     -6      5 |
 *     |    0      0   -4.5 |        |   -3   -0.5   -4.5 |
 */

double lu_decomposition_seq(int dim, float **alu) {
    clock_t exec_start = clock();

    for (int j = 0; j < dim - 1; j++) {
        for (int i = j + 1; i < dim; i++) {
            alu[i][j] = alu[i][j] / alu[j][j]; // => L

            for (int k = j + 1; k < dim; k++) {
                alu[i][k] = alu[i][k] - alu[i][j] * alu[j][k]; // => U
            }
        }
    }

    return (double) (clock() - exec_start) / CLOCKS_PER_SEC;
}

int mpi_init(int *rank, int *nprocs, char *name, int *name_len) {
    int status;

    if ((status = MPI_Init(NULL, NULL)) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, status);
    } else if ((status = MPI_Comm_rank(MPI_COMM_WORLD, rank)) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, status);
    } else if ((status = MPI_Comm_size(MPI_COMM_WORLD, nprocs)) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, status);
    } else if ((status = MPI_Get_processor_name(name, name_len)) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, status);
    }

    return status;
}

// A lot of time is spent on sending matrices when running on multiple machines.
double lu_decomposition_mpi_original(int dim, float **alu, int rank, int nprocs) {
    double exec_start, exec_time = 0;

    if (rank == MAIN_PROC) {
        exec_start = MPI_Wtime();
    }

    for (int j = 0; j < dim - 1; j++) {
        for (int i = j + 1; i < dim; i++) {
            if ((i % nprocs) == rank) {
                alu[i][j] = alu[i][j] / alu[j][j]; // => L

                for (int k = j + 1; k < dim; k++) {
                    alu[i][k] = alu[i][k] - alu[i][j] * alu[j][k]; // => U
                }
            }
        }

        for (int i = j + 1; i < dim; i++) {
            MPI_Bcast(&alu[i][j], dim - j, MPI_FLOAT, i % nprocs, MPI_COMM_WORLD);
        }
    }

    if (rank == MAIN_PROC) {
        exec_time = MPI_Wtime() - exec_start;
    }

    return exec_time;
}

double lu_decomposition_mpi_chunk(int dim, float **alu, int rank, int nprocs) {
    double exec_start, exec_time = 0;

    if (rank == MAIN_PROC) {
        exec_start = MPI_Wtime();
    }

    for (int j = 0; j < dim - 1; j++) {
        int row_count = dim - (j + 1);
        int chunk_size = row_count / nprocs;
        chunk_size = chunk_size == 0 ? 1 : chunk_size;
        int chunk_count = (row_count / chunk_size) + (row_count % chunk_size == 0 ? 0 : 1);
        chunk_count = chunk_count > nprocs ? nprocs : chunk_count;
        int last_chunk_size = chunk_size + (row_count - chunk_size * chunk_count);

        for (int i = j + 1; i < dim; i++) {
            int chunk_index = (i - (j + 1)) / chunk_size;

            if (chunk_index >= (nprocs - 1)) {
                chunk_index = nprocs - 1;
            }

            //if (rank == MAIN_PROC) {
            //    printf("[%d][%d]   row_count %d   chunk_count %d   chunk_size %d   last_chunk_size %d   chunk_index %d\n",
            //           i,
            //           j,
            //           row_count,
            //           chunk_count,
            //           chunk_size,
            //           last_chunk_size,
            //           chunk_index);
            //}

            if (chunk_index == rank) {
                alu[i][j] = alu[i][j] / alu[j][j]; // => L

                for (int k = j + 1; k < dim; k++) {
                    alu[i][k] = alu[i][k] - alu[i][j] * alu[j][k]; // => U
                }
            }
        }

        for (int i = j + 1; i < dim;) {
            int chunk_index = (i - (j + 1)) / chunk_size;

            if (chunk_index >= (nprocs - 1)) {
                chunk_index = nprocs - 1;
                chunk_size = last_chunk_size;
            }

            MPI_Bcast(&alu[i][j], chunk_size * dim - j, MPI_FLOAT, chunk_index, MPI_COMM_WORLD);

            i += chunk_size;
        }
    }

    if (rank == MAIN_PROC) {
        exec_time = MPI_Wtime() - exec_start;
    }

    return exec_time;
}

// Uniform loading of processors is not ensured.
// The solution is to periodically recalculate the ranges.
double lu_decomposition_mpi_optimise(int dim, float **alu, int rank, int nprocs) {
    double exec_start, exec_time = 0;

    if (rank == MAIN_PROC) {
        exec_start = MPI_Wtime();
    }

    int row_count = dim - 1;
    int range_size = row_count / nprocs;
    range_size = row_count % nprocs == 0
            ? range_size
            : (range_size + 1);

    for (int j = 0; j < dim - 1; j++) {
        for (int i = j + 1; i < dim; i++) {
            if (((i - 1) / range_size) == rank) {
                alu[i][j] = alu[i][j] / alu[j][j]; // => L

                for (int k = j + 1; k < dim; k++) {
                    alu[i][k] = alu[i][k] - alu[i][j] * alu[j][k]; // => U
                }
            }
        }

        MPI_Bcast(&alu[j + 1][0], dim, MPI_FLOAT, j / range_size, MPI_COMM_WORLD);
    }

    if (rank == MAIN_PROC) {
        exec_time = MPI_Wtime() - exec_start;
    }

    return exec_time;
}

int main(int argc, char **argv) {
    char *char_buffer = NULL;
    int dim = 5;
    bool print_matrix = false;

    if (argc == 2) {
        if ((dim = (int) strtol(argv[1], &char_buffer, 10)) < 2 || *char_buffer != '\0') {
            return ERR_ARGV;
        }
    } else if (argc == 3) {
        if ((dim = (int) strtol(argv[1], &char_buffer, 10)) < 2 || *char_buffer != '\0') {
            return ERR_ARGV;
        }
        if (strcmp("true", argv[2]) == 0) {
            print_matrix = true;
        } else if (strcmp("false", argv[2]) == 0) {
            print_matrix = false;
        } else {
            return ERR_ARGV;
        }
    }

    float **a1 = alloc_matrix_f(dim, dim);

    if (a1 == NULL) {
        return ERR_MEMORY_ALLOCATION;
    }

    srand(time(NULL));
    a1 = random_matrix_init_f(dim, dim, a1);
    //a1[0][0] = 2;
    //a1[0][1] = 4;
    //a1[0][2] = -4;
    //a1[1][0] = 1;
    //a1[1][1] = -4;
    //a1[1][2] = 3;
    //a1[2][0] = -6;
    //a1[2][1] = -9;
    //a1[2][2] = 5;

    float **a2 = copy_matrix_f(dim, dim, a1);

    if (a2 == NULL) {
        free(a1);

        return ERR_MEMORY_ALLOCATION;
    }

    float **a3 = copy_matrix_f(dim, dim, a1);

    if (a3 == NULL) {
        free(a1);
        free(a2);

        return ERR_MEMORY_ALLOCATION;
    }

    int rank, nprocs, name_len;
    char name[MPI_MAX_PROCESSOR_NAME];

    if (mpi_init(&rank, &nprocs, name, &name_len) != MPI_SUCCESS) {
        free(a1);
        free(a2);
        free(a3);

        return ERR_MPI_INIT;
    }

    if (rank == MAIN_PROC) {
        if (print_matrix) {
            printf("A:\n");
            print_matrix_f(dim, dim, a1);
        }

        double seq_time = lu_decomposition_seq(dim, a1);

        if (fprint("seq.log", 1, dim, seq_time, name, name_len) != 0) {
            MPI_Abort(MPI_COMM_WORLD, ERR_FILE_OPEN);
            free(a1);
            free(a2);
            free(a3);

            return ERR_FILE_OPEN;
        }

        if (print_matrix) {
            printf("L (seq):\n");
            print_alu_matrix(dim, a1, ALU_PRINT_MODE_L);

            printf("U (seq):\n");
            print_alu_matrix(dim, a1, ALU_PRINT_MODE_U);
        }
    }

    double mpi_time = lu_decomposition_mpi_optimise(dim, a2, rank, nprocs);

    if (rank == MAIN_PROC) {
        if (fprint("mpi_optimise.log", nprocs, dim, mpi_time, name, name_len) != 0) {
            MPI_Abort(MPI_COMM_WORLD, ERR_FILE_OPEN);
            free(a1);
            free(a2);
            free(a3);

            return ERR_FILE_OPEN;
        }

        if (print_matrix) {
            printf("L (optimise):\n");
            print_alu_matrix(dim, a2, ALU_PRINT_MODE_L);

            printf("U (optimise):\n");
            print_alu_matrix(dim, a2, ALU_PRINT_MODE_U);
        }
    }

    double mpi_chunk_time = lu_decomposition_mpi_chunk(dim, a3, rank, nprocs);

    if (rank == MAIN_PROC) {
        if (fprint("mpi_chunk.log", nprocs, dim, mpi_chunk_time, name, name_len) != 0) {
            MPI_Abort(MPI_COMM_WORLD, ERR_FILE_OPEN);
            free(a1);
            free(a2);
            free(a3);

            return ERR_FILE_OPEN;
        }

        if (print_matrix) {
            printf("L (mpi chunk):\n");
            print_alu_matrix(dim, a3, ALU_PRINT_MODE_L);

            printf("U (mpi chunk):\n");
            print_alu_matrix(dim, a3, ALU_PRINT_MODE_U);
        }
    }

    MPI_Finalize();
    free(a1);
    free(a2);
    free(a3);

    return 0;
}
