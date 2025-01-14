#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int N, M;

char **current_iteration;
char **new_iteration;

// g++ -o gol gol.cpp
// gol input_file.txt 20 output_file.txt

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

    for (int i = 1; i < M - 1; i++)
    {
        current_iteration[0][i] = current_iteration[N - 2][i];
    }

    // bottom line

    for (int i = 1; i < M - 1; i++)
    {
        current_iteration[N - 1][i] = current_iteration[1][i];
    }

    // left column

    for (int i = 1; i < N - 1; i++)
    {
        current_iteration[i][0] = current_iteration[i][M - 2];
    }

    // right column

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

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <iterations> <output_file>" << endl;
        return 1;
    }

    // Get the number of iterations from command-line arguments
    int iterations = atoi(argv[2]);

    // Start time measurement
    auto start = high_resolution_clock::now();

    // Allocate memory for the grid
    alloc_memory(argv[1]);

    // Read initial state from file
    read_initial_state(argv[1]);

    // Run the simulation for the specified number of iterations
    for (int i = 0; i < iterations; i++)
    {
        // Update boundary cells
        add_border();

        // Apply Conway's Game of Life rules to update grid
        apply_algorithm();

        // Copy new state to current state
        copy_matrix();
    }

    // End time measurement
    auto stop = high_resolution_clock::now();

    // Calculate the total time taken
    auto duration = duration_cast<seconds>(stop - start);

    cout << "Total time taken: " << duration.count() << " seconds" << endl;

    // Write final state to output file
    write_final_state(argv[3]);

    return 0;
}
