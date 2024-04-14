#include <iostream>
#include <cstdlib>
#include <fstream>
#include <omp.h>

using namespace std;

// Global variables
int N, M;
char **current_iteration;
char **new_iteration;

// Function prototypes
void alloc_memory(char *filename);
void read_initial_state(char *filename);
void add_border();
char analyze_cell(int i, int j);
void apply_algorithm();
void copy_matrix();
void write_final_state(char *filename);

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <iterations> <output_file>" << endl;
        return 1;
    }

    double one_thread_time;

    int iterations = atoi(argv[2]);

    alloc_memory(argv[1]);

    read_initial_state(argv[1]);

    // Start time measurement

    for (int j = 1; j <= 8; j++)
    {
        // set no of threads
        double start_time = omp_get_wtime();
        omp_set_num_threads(j);
        for (int i = 0; i < iterations; i++)
        {
            // Update boundary cells
            add_border();

            // Apply Conway's Game of Life rules to update grid
            apply_algorithm();

            // Copy new state to current state
            copy_matrix();
        }

        // Write final state to output file
        write_final_state(argv[3]);

        int num_threads;
#pragma omp parallel
        {
#pragma omp single
            {
                num_threads = omp_get_num_threads();
            }
        }
        cout << "Number of processing elements: " << num_threads << endl;

        // End time measurement
        double end_time = omp_get_wtime();

        // Calculate the total time taken
        double total_time = end_time - start_time;

        if (j == 1)
            one_thread_time = total_time;

        cout << "Total time taken: " << total_time << " seconds" << endl;
        cout << "Speed Up: " << one_thread_time / total_time << " seconds" << endl;
        cout << "Parallel Efficiency: " << ((one_thread_time / total_time) / j) * 100 << endl;
    }

    return 0;
}

void alloc_memory(char *filename)
{
    ifstream f(filename);

    f >> N >> M;

    f.close();

    N = N + 2;
    M = M + 2;

    current_iteration = new char *[N];
    new_iteration = new char *[N];

    for (int i = 0; i < N; i++)
    {
        current_iteration[i] = new char[M];
    }

    for (int i = 0; i < N; i++)
    {
        new_iteration[i] = new char[M];
    }
}

void read_initial_state(char *filename)
{
    ifstream f(filename);

    int dummy;

    f >> dummy >> dummy;

    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            f >> current_iteration[i][j];
        }
    }

    f.close();
}

void add_border()
{
    int i;
// top line
#pragma omp parallel for private(i)
    for (int i = 1; i < M - 1; i++)
    {
        current_iteration[0][i] = current_iteration[N - 2][i];
    }

// bottom line
#pragma omp parallel for private(i)
    for (int i = 1; i < M - 1; i++)
    {
        current_iteration[N - 1][i] = current_iteration[1][i];
    }

// left column
#pragma omp parallel for private(i)
    for (int i = 1; i < N - 1; i++)
    {
        current_iteration[i][0] = current_iteration[i][M - 2];
    }

// right column
#pragma omp parallel for private(i)
    for (int i = 1; i < N - 1; i++)
    {
        current_iteration[i][M - 1] = current_iteration[i][1];
    }

    // corners
    current_iteration[0][0] = current_iteration[N - 2][M - 2]; // upper left
    current_iteration[0][M - 1] = current_iteration[N - 2][1]; // upper right
    current_iteration[N - 1][0] = current_iteration[1][M - 2]; // lower left
    current_iteration[N - 1][M - 1] = current_iteration[1][1]; // lower right
}

char analyze_cell(int i, int j)
{
    int alive_neighbours = 0, dead_neighbours = 0;

    current_iteration[i - 1][j - 1] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i - 1][j] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i - 1][j + 1] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i][j + 1] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i + 1][j + 1] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i + 1][j] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i + 1][j - 1] == '.' ? dead_neighbours++ : alive_neighbours++;
    current_iteration[i][j - 1] == '.' ? dead_neighbours++ : alive_neighbours++;

    if (alive_neighbours < 2)
        return '.';

    else if (alive_neighbours > 3)
        return '.';

    else if (current_iteration[i][j] == 'X' && (alive_neighbours == 2 || alive_neighbours == 3))
        return 'X';

    else if (current_iteration[i][j] == '.' && alive_neighbours == 3)
        return 'X';

    return '.';
}

void apply_algorithm()
{
    int i, j = 0;

#pragma omp parallel for private(i, j)
    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            new_iteration[i][j] = analyze_cell(i, j);
        }
    }
}

void copy_matrix()
{
    int i, j = 0;

#pragma omp parallel for private(i, j)
    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            current_iteration[i][j] = new_iteration[i][j];
        }
    }
}

void write_final_state(char *filename)
{
    ofstream f(filename);

    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            f << current_iteration[i][j] << " ";
        }
        f << endl;
    }

    f.close();
}
