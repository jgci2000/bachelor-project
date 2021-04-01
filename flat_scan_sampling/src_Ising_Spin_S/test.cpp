#include <iostream>

int main()
{
    int **a = new int*[3];
    for (int i = 0; i < 3; i++)
        a[i] = new int[i + 1];
    
    int c = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < i + 1; j++)
            a[i][j] = c++;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < i + 1; j++)
            std::cout << a[i][j] << " ";
        std::cout << std::endl;;
    }
    
    for (int i = 0; i < 3; i++)
        delete[] a[i];

    delete[] a;
    return 0;
}

