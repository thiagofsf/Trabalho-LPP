#include <omp.h>
#include <bits/stdc++.h>

using namespace std;

//Função que imprime uma matriz
//Entradas: tamanho n da matriz, ponteiro para a matriz
void imprime_matriz(int n, int** mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

//Função que aloca uma matriz na memória
//Entradas: tamanho n da matriz
//saída: ponteiro para matriz alocada
int** alocar_matriz(int n)
{
    //Aloca vetor com tamanho total de elementos, aloca vator principal da matriz
    int* data = (int*)malloc(n * n * sizeof(int));
    int** array = (int**)malloc(n * sizeof(int*));
    //associando
    for (int i = 0; i < n; i++)
    {
        array[i] = &(data[n * i]);
    }
    return array;
}

//Preenche a matriz com o valor desejado
//entradas: tamanho n da matriz, ponteiro para matriz a ser preenchida, valor a preencher
void preenche_matriz(int n, int**& mat, int val)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            mat[i][j] = val;
        }
    }
}

//libera matriz
//entradas: tamanho n da matriz, ponteiro para matriz
void libera_matriz(int n, int** mat)
{
    free(mat[0]);
    free(mat);
}

//Multiplicação direta (usada para fim da recursão no strassen) - Aqui implementada usando paralelismo também
//entradas: tamanho n da matriz, ponteiros para matriz 1 e matriz 2
//saida: ponteiro para matriz produto
int** multiplica_matrizes(int n, int** mat1, int** mat2)
{
    int** prod = alocar_matriz(n);

    int i, j;

#pragma omp parallel for collapse(2)
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            prod[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                prod[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return prod;
}

//obtém submatriz quadrante de uma matriz maior
//entradas: tamanho n da matriz, ponteiro para matriz, valores de começo de linha e coluna (offset)
//saída: ponteiro para submatriz resultante
int** obter_submatriz(int n, int** mat, int lin, int col)
{
    int m = n / 2;
    int** slice = alocar_matriz(m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            slice[i][j] = mat[lin + i][col + j];
        }
    }
    return slice;
}

//Soma duas matrizes e retorna a matriz resultante
//Entradas: tamanho n da matriz, ponteiros para matriz 1 e matriz 2, flag (true para soma, false para subtração)
//Saida: ponteiro para matriz resultado
int** soma_matrizes(int n, int** mat1, int** mat2, bool add)
{
    int** result = alocar_matriz(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (add)
                result[i][j] = mat1[i][j] + mat2[i][j];
            else
                result[i][j] = mat1[i][j] - mat2[i][j];
        }
    }

    return result;
}

//Combina 4 submatrizes tamanho m em 1 matriz tamanho n=m*2
//entradas: tamanho m das submatrizes, ponteiros para submatrizes quadrantes: 11, 12, 21, 22
//saída: matriz resultado da combinaçao das submatrizes
int** combina_matrizes(int m, int** c11, int** c12, int** c21, int** c22)
{
    int n = 2 * m;
    int** result = alocar_matriz(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < m && j < m)
                result[i][j] = c11[i][j];
            else if (i < m)
                result[i][j] = c12[i][j - m];
            else if (j < m)
                result[i][j] = c21[i - m][j];
            else
                result[i][j] = c22[i - m][j - m];
        }
    }

    return result;
}

//Implementação paralela de Strassen usando OpenMP
//Entradas: tamanho n da matriz, ponteiros para matriz 1 e matriz 2
int** strassen(int n, int** mat1, int** mat2)
{
    //Condição final de recursão
    if (n <= 64)
    {
        return multiplica_matrizes(n, mat1, mat2);
    }

    int m = n / 2;

    //Obtendo submatrizes quadrantes
    int** a = obter_submatriz(n, mat1, 0, 0);
    int** b = obter_submatriz(n, mat1, 0, m);
    int** c = obter_submatriz(n, mat1, m, 0);
    int** d = obter_submatriz(n, mat1, m, m);
    int** e = obter_submatriz(n, mat2, 0, 0);
    int** f = obter_submatriz(n, mat2, 0, m);
    int** g = obter_submatriz(n, mat2, m, 0);
    int** h = obter_submatriz(n, mat2, m, m);

    //Obtendo P1
    int** p1;
    #pragma omp task shared(p1)
    {
        int** bds = soma_matrizes(m, b, d, false);
        int** gha = soma_matrizes(m, g, h, true);
        p1 = strassen(m, bds, gha);
        libera_matriz(m, bds);
        libera_matriz(m, gha);
    }

    //Obtendo P2
    int** p2;
    #pragma omp task shared(p2)
    {
        int** ada = soma_matrizes(m, a, d, true);
        int** eha = soma_matrizes(m, e, h, true);
        p2 = strassen(m, ada, eha);
        libera_matriz(m, ada);
        libera_matriz(m, eha);
    }

    //Obtendo P3
    int** p3;
    #pragma omp task shared(p3)
    {
        int** acs = soma_matrizes(m, a, c, false);
        int** efa = soma_matrizes(m, e, f, true);
        p3 = strassen(m, acs, efa);
        libera_matriz(m, acs);
        libera_matriz(m, efa);
    }

    //Obtendo P4
    int** p4;
    #pragma omp task shared(p4)
    {
        int** aba = soma_matrizes(m, a, b, true);
        p4 = strassen(m, aba, h);
        libera_matriz(m, aba);
    }

    //Obtendo P5
    int** p5;
    #pragma omp task shared(p5)
    {
        int** fhs = soma_matrizes(m, f, h, false);
        p5 = strassen(m, a, fhs);
        libera_matriz(m, fhs);
    }

    //Obtendo P6
    int** p6;
    #pragma omp task shared(p6)
    {
        int** ges = soma_matrizes(m, g, e, false);
        p6 = strassen(m, d, ges);
        libera_matriz(m, ges);
    }

    //Obtendo P7
    int** p7;
    #pragma omp task shared(p7)
    {
        int** cda = soma_matrizes(m, c, d, true);
        p7 = strassen(m, cda, e);
        libera_matriz(m, cda);
    }

    //Aguardar threads processarem tarefas
    #pragma omp taskwait 

    libera_matriz(m, a);
    libera_matriz(m, b);
    libera_matriz(m, c);
    libera_matriz(m, d);
    libera_matriz(m, e);
    libera_matriz(m, f);
    libera_matriz(m, g);
    libera_matriz(m, h);

    //Matrizes P1-P7 obtidas, agora vamos gerar as submatrizes C11, C12, C21, C22

    //Gerando C11
    int** c11;
    #pragma omp task shared(c11)
    {
        int** s1s2a = soma_matrizes(m, p1, p2, true);
        int** s6s4s = soma_matrizes(m, p6, p4, false);
        c11 = soma_matrizes(m, s1s2a, s6s4s, true);
        libera_matriz(m, s1s2a);
        libera_matriz(m, s6s4s);
    }

    //Gerando C12
    int** c12;
    #pragma omp task shared(c12)
    {
        c12 = soma_matrizes(m, p4, p5, true);
    }

    //Gerando C21
    int** c21;
    #pragma omp task shared(c21)
    {
        c21 = soma_matrizes(m, p6, p7, true);
    }

    //Gerando C22
    int** c22;
    #pragma omp task shared(c22)
    {
        int** s2s3s = soma_matrizes(m, p2, p3, false);
        int** s5s7s = soma_matrizes(m, p5, p7, false);
        c22 = soma_matrizes(m, s2s3s, s5s7s, true);
        libera_matriz(m, s2s3s);
        libera_matriz(m, s5s7s);
    }

    //Aguardar threads processarem tarefas
    #pragma omp taskwait

    libera_matriz(m, p1);
    libera_matriz(m, p2);
    libera_matriz(m, p3);
    libera_matriz(m, p4);
    libera_matriz(m, p5);
    libera_matriz(m, p6);
    libera_matriz(m, p7);

    //Agora vamos obter a matriz produto pela combinação das submatrizes C
    int** prod = combina_matrizes(m, c11, c12, c21, c22);

    libera_matriz(m, c11);
    libera_matriz(m, c12);
    libera_matriz(m, c21);
    libera_matriz(m, c22);

    return prod;
}

//Função que checa se a multiplicação está correta, usada para testes
bool check(int n, int** prod1, int** prod2)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (prod1[i][j] != prod2[i][j])
                return false;
        }
    }
    return true;
}

//Procedimento Principal
int main()
{
    int p;
    cout << "\nInsira o número de threads: ";
    cin >> p;

    int n;
    cout << "\nInsira a dimensão n da matriz: ";
    cin >> n;

    int** mat1 = alocar_matriz(n);
    preenche_matriz(n, mat1, 2);

    int** mat2 = alocar_matriz(n);
    preenche_matriz(n, mat2, 4);

    //Imprimindo matrizes pequenas para conferir
    if(n<32){
        cout << "\nImprimindo matriz A:\n";
        imprime_matriz(n, mat1);
        cout << "\nImprimindo matriz B:\n";
        imprime_matriz(n, mat2);
    }

    //Obtendo tempo inicial
    double startParStrassen = omp_get_wtime();
    
    int** prod;
    //setando numero de threads
    omp_set_num_threads(p);

    //iniciando a multiplicação em paralelo
    #pragma omp parallel
    {
    #pragma omp single
        {
            prod = strassen(n, mat1, mat2);
        }
    }

    //obtendo tempo final
    double endParStrassen = omp_get_wtime();

    //Imprimindo matrizes pequenas para conferir
    if(n<32){
        cout << "\nImprimindo matriz C:\n";
        imprime_matriz(n, prod);
    }

    cout << "\nTempo de execução do Strassen Paralelo (OMP): ";
    cout << setprecision(5) << endParStrassen - startParStrassen << endl;

    cout << endl;

    return 0;
}