#include <bits/stdc++.h>
#include <time.h>
#include <sys/resource.h>

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

//Multiplicação serial direta (usada para fim da recursão no strassen)
//entradas: tamanho n da matriz, ponteiros para matriz 1 e matriz 2
//saida: ponteiro para matriz produto
int** multiplica_matrizes(int n, int** mat1, int** mat2)
{
    int** prod = alocar_matriz(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
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

//Função strassen Serial
//entradas: tamanho n da matriz, ponteiro para matriz 1 e ponteiro para matriz 2
//saida: ponteiro para matriz produto da multiplicação
int** strassen(int n, int** mat1, int** mat2)
{

    //Condição de parada de recursão
    if (n <= 64)
    {
        return multiplica_matrizes(n, mat1, mat2);
    }

    int m = n / 2;

    //Obtendo submatrizes quadrantes:
    int** a = obter_submatriz(n, mat1, 0, 0);
    int** b = obter_submatriz(n, mat1, 0, m);
    int** c = obter_submatriz(n, mat1, m, 0);
    int** d = obter_submatriz(n, mat1, m, m);
    int** e = obter_submatriz(n, mat2, 0, 0);
    int** f = obter_submatriz(n, mat2, 0, m);
    int** g = obter_submatriz(n, mat2, m, 0);
    int** h = obter_submatriz(n, mat2, m, m);

    //Obtendo as submatrizes P1 a P7
    
    int** bds = soma_matrizes(m, b, d, false);
    int** gha = soma_matrizes(m, g, h, true);
    int** p1 = strassen(m, bds, gha);
    libera_matriz(m, bds);
    libera_matriz(m, gha);

    int** ada = soma_matrizes(m, a, d, true);
    int** eha = soma_matrizes(m, e, h, true);
    int** p2 = strassen(m, ada, eha);
    libera_matriz(m, ada);
    libera_matriz(m, eha);

    int** acs = soma_matrizes(m, a, c, false);
    int** efa = soma_matrizes(m, e, f, true);
    int** p3 = strassen(m, acs, efa);
    libera_matriz(m, acs);
    libera_matriz(m, efa);

    int** aba = soma_matrizes(m, a, b, true);
    int** p4 = strassen(m, aba, h);
    libera_matriz(m, aba);
    libera_matriz(m, b);

    int** fhs = soma_matrizes(m, f, h, false);
    int** p5 = strassen(m, a, fhs);
    libera_matriz(m, fhs);
    libera_matriz(m, a);
    libera_matriz(m, f);
    libera_matriz(m, h);

    int** ges = soma_matrizes(m, g, e, false);
    int** p6 = strassen(m, d, ges);
    libera_matriz(m, ges);
    libera_matriz(m, g);

    int** cda = soma_matrizes(m, c, d, true);
    int** p7 = strassen(m, cda, e);
    libera_matriz(m, cda);
    libera_matriz(m, c);
    libera_matriz(m, d);
    libera_matriz(m, e);

    //Matrizes P1 - P7 calculadas, agora vamos obter as quatro submatrizes C11, C12, C21, C22

    int** s1s2a = soma_matrizes(m, p1, p2, true);
    int** s6s4s = soma_matrizes(m, p6, p4, false);
    int** c11 = soma_matrizes(m, s1s2a, s6s4s, true);
    libera_matriz(m, s1s2a);
    libera_matriz(m, s6s4s);
    libera_matriz(m, p1);

    int** c12 = soma_matrizes(m, p4, p5, true);
    libera_matriz(m, p4);

    int** c21 = soma_matrizes(m, p6, p7, true);
    libera_matriz(m, p6);

    int** s2s3s = soma_matrizes(m, p2, p3, false);
    int** s5s7s = soma_matrizes(m, p5, p7, false);
    int** c22 = soma_matrizes(m, s2s3s, s5s7s, true);
    libera_matriz(m, s2s3s);
    libera_matriz(m, s5s7s);
    libera_matriz(m, p2);
    libera_matriz(m, p3);
    libera_matriz(m, p5);
    libera_matriz(m, p7);

    //Após obtidas as submatrizes serão combinadas para gerar a matriz produto
    int** prod = combina_matrizes(m, c11, c12, c21, c22);

    libera_matriz(m, c11);
    libera_matriz(m, c12);
    libera_matriz(m, c21);
    libera_matriz(m, c22);

    return prod;
}

void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total)
{
    long seg_CPU, seg_sistema, mseg_CPU, mseg_sistema;
    struct rusage ptempo;

    getrusage(0,&ptempo);

    seg_CPU = ptempo.ru_utime.tv_sec;
    mseg_CPU = ptempo.ru_utime.tv_usec;
    seg_sistema = ptempo.ru_stime.tv_sec;
    mseg_sistema = ptempo.ru_stime.tv_usec;

    *seg_CPU_total     = (seg_CPU + 0.000001 * mseg_CPU);
    *seg_sistema_total = (seg_sistema + 0.000001 * mseg_sistema);
}

//Procedimento principal
int main(int argc, char* argv[])
{
    double s_CPU_inicial=0, s_total_inicial=0, s_CPU_final=0, s_total_final=0;

    int n;
    //definindo n
    cout << endl;
    cout << "Inserir dimensão n da matriz: ";
    cin >> n;

    //alocar matrizes
    int** mat1 = alocar_matriz(n);
    int** mat2 = alocar_matriz(n);

    //preencher matrizes
    preenche_matriz(n, mat1, 2);
    preenche_matriz(n, mat2, 4);

    if(n<32){
        cout << "\nImprimindo matriz A:\n";
        imprime_matriz(n, mat1);
        cout << "\nImprimindo matriz B:\n";
        imprime_matriz(n, mat2);
    }

    //obter tempo inicial
    Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);

    //realizar multiplicação
    int** prod;

    prod = strassen(n, mat1, mat2);

    //obter tempo final
    Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);

    cout << "\nTempo de execução do Strassen Sequencial: ";
    cout << setprecision(5) << s_CPU_final - s_CPU_inicial << endl;
    cout << endl;

    if(n<32){
        cout << "\nImprimindo matriz C:\n";
        imprime_matriz(n, prod);
    }

    return 0;
}