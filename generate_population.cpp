#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>

int main()
{
    int N = 100;
    std::srand(std::time(0));
    std::vector<std::vector<char>> matrix(N, std::vector<char>(N));

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            matrix[i][j] = (std::rand() % 2 == 0) ? '.' : 'X';
        }
    }

    std::ofstream file("input_file.txt");
    if (file.is_open())
    {
        file << N << " " << N << '\n';
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                file << matrix[i][j];
            }
            file << '\n';
        }
        file.close();
    }
    else
    {
        std::cout << "Unable to open file";
    }

    return 0;
}
