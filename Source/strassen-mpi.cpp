#include <mpi.h>
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

//Função strassen Serial, usada nos processos
//entradas: tamanho n da matriz, ponteiro para matriz 1 e ponteiro para matriz 2
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

    //Obtenso as submatrizes P1 a P7
    
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

//Função Strassen paralela implementada com logica MPI
//Entradas: tamanho n da matriz, ponteiros para matriz 1, matriz 2, matriz produto, o número do processo e a quantidade de processadores disponiveis
void strassen(int n, int** mat1, int** mat2, int**& prod, int rank, int num_process)
{

    //Condição base de strassen
    if (n == 1)
    {
        prod = alocar_matriz(1);
        prod[0][0] = mat1[0][0] * mat2[0][0];
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

    //Alocando P1 - P7
    int** p1 = alocar_matriz(m);
    int** p2 = alocar_matriz(m);
    int** p3 = alocar_matriz(m);
    int** p4 = alocar_matriz(m);
    int** p5 = alocar_matriz(m);
    int** p6 = alocar_matriz(m);
    int** p7 = alocar_matriz(m);

    //Definir de acordo com a quantidade de processadores como paralelizar:
    if(num_process == 1){
        if (rank == 0){
            //Calculando P1
            int** bds = soma_matrizes(m, b, d, false);
            int** gha = soma_matrizes(m, g, h, true);
            p1 = strassen(m, bds, gha);
            //Calculando P2
            int** ada = soma_matrizes(m, a, d, true);
            int** eha = soma_matrizes(m, e, h, true);
            p2 = strassen(m, ada, eha);
            //Calculando P3
            int** acs = soma_matrizes(m, a, c, false);
            int** efa = soma_matrizes(m, e, f, true);
            p3 = strassen(m, acs, efa);
             //Calculando P4
            int** aba = soma_matrizes(m, a, b, true);
            p4 = strassen(m, aba, h);
             //Calculando P5
            int** fhs = soma_matrizes(m, f, h, false);
            p5 = strassen(m, a, fhs);
            //Calculando P6
            int** ges = soma_matrizes(m, g, e, false);
            p6 = strassen(m, d, ges);
            //Calculando P7
            int** cda = soma_matrizes(m, c, d, true);
            p7 = strassen(m, cda, e);

            libera_matriz(m, bds);
            libera_matriz(m, gha);
            libera_matriz(m, ada);
            libera_matriz(m, eha);
            libera_matriz(m, acs);
            libera_matriz(m, efa);
            libera_matriz(m, aba);
            libera_matriz(m, fhs);
            libera_matriz(m, ges);
            libera_matriz(m, cda);
        }
        libera_matriz(m, b);
        libera_matriz(m, a);
        libera_matriz(m, f);
        libera_matriz(m, h);
        libera_matriz(m, g);
        libera_matriz(m, c);
        libera_matriz(m, d);
        libera_matriz(m, e);
    }
    else if(num_process == 2){
        //Processo 0 verifica e recebe resultados dos outros processos
        if (rank == 0)
        {
            cout << "Numero de processos: " << num_process << endl;
            MPI_Recv(&(p1[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p2[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p3[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p4[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p5[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p6[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p7[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //processo 1 calcula P1 - P7
        if (rank == 1)
        {
            //Calculando P1
            int** bds = soma_matrizes(m, b, d, false);
            int** gha = soma_matrizes(m, g, h, true);
            p1 = strassen(m, bds, gha);
            //Calculando P2
            int** ada = soma_matrizes(m, a, d, true);
            int** eha = soma_matrizes(m, e, h, true);
            p2 = strassen(m, ada, eha);
            //Calculando P3
            int** acs = soma_matrizes(m, a, c, false);
            int** efa = soma_matrizes(m, e, f, true);
            p3 = strassen(m, acs, efa);
             //Calculando P4
            int** aba = soma_matrizes(m, a, b, true);
            p4 = strassen(m, aba, h);
             //Calculando P5
            int** fhs = soma_matrizes(m, f, h, false);
            p5 = strassen(m, a, fhs);
            //Calculando P6
            int** ges = soma_matrizes(m, g, e, false);
            p6 = strassen(m, d, ges);
            //Calculando P7
            int** cda = soma_matrizes(m, c, d, true);
            p7 = strassen(m, cda, e);

            libera_matriz(m, bds);
            libera_matriz(m, gha);
            libera_matriz(m, ada);
            libera_matriz(m, eha);
            libera_matriz(m, acs);
            libera_matriz(m, efa);
            libera_matriz(m, aba);
            libera_matriz(m, fhs);
            libera_matriz(m, ges);
            libera_matriz(m, cda);
            
            MPI_Send(&(p1[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p2[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p3[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p4[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p5[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p6[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p7[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        libera_matriz(m, b);
        libera_matriz(m, a);
        libera_matriz(m, f);
        libera_matriz(m, h);
        libera_matriz(m, g);
        libera_matriz(m, c);
        libera_matriz(m, d);
        libera_matriz(m, e);
    }
    else if(num_process == 4){
        //Processo 0 verifica e recebe resultados dos outros processos
        if (rank == 0)
        {
            cout << "Numero de processos: " << num_process << endl;
            MPI_Recv(&(p1[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p2[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p3[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p4[0][0]), m * m, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p5[0][0]), m * m, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p6[0][0]), m * m, MPI_INT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p7[0][0]), m * m, MPI_INT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //processo 1 calcula P1, P2 e P3
        if (rank == 1)
        {
            //Calculando P1
            int** bds = soma_matrizes(m, b, d, false);
            int** gha = soma_matrizes(m, g, h, true);
            p1 = strassen(m, bds, gha);
            //Calculando P2
            int** ada = soma_matrizes(m, a, d, true);
            int** eha = soma_matrizes(m, e, h, true);
            p2 = strassen(m, ada, eha);
            //Calculando P3
            int** acs = soma_matrizes(m, a, c, false);
            int** efa = soma_matrizes(m, e, f, true);
            p3 = strassen(m, acs, efa);

            libera_matriz(m, bds);
            libera_matriz(m, gha);
            libera_matriz(m, ada);
            libera_matriz(m, eha);
            libera_matriz(m, acs);
            libera_matriz(m, efa);
            
            MPI_Send(&(p1[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p2[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p3[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 2 calcula P4 e P5
        if (rank == 2)
        {
            //Calculando P4
            int** aba = soma_matrizes(m, a, b, true);
            p4 = strassen(m, aba, h);
            //Calculando P5
            int** fhs = soma_matrizes(m, f, h, false);
            p5 = strassen(m, a, fhs);

            libera_matriz(m, aba);
            libera_matriz(m, fhs);

            MPI_Send(&(p4[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p5[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, b);
        libera_matriz(m, a);
        libera_matriz(m, f);
        libera_matriz(m, h);

        //Processo 3 calcula P6 e P7
        if (rank == 3)
        {   
            //Calculando P6
            int** ges = soma_matrizes(m, g, e, false);
            p6 = strassen(m, d, ges);
            //Calculando P7
            int** cda = soma_matrizes(m, c, d, true);
            p7 = strassen(m, cda, e);

            libera_matriz(m, ges);
            libera_matriz(m, cda);

            MPI_Send(&(p6[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p7[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, g);
        libera_matriz(m, c);
        libera_matriz(m, d);
        libera_matriz(m, e);
    }
    else if(num_process == 6){
        //Processo 0 verifica e recebe resultados dos outros processos
        if (rank == 0)
        {
            cout << "Numero de processos: " << num_process << endl;
            MPI_Recv(&(p1[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p2[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p3[0][0]), m * m, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p4[0][0]), m * m, MPI_INT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p5[0][0]), m * m, MPI_INT, 4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p6[0][0]), m * m, MPI_INT, 5, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p7[0][0]), m * m, MPI_INT, 5, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //processo 1 calcula P1 e P2
        if (rank == 1)
        {
            //Calculando P1
            int** bds = soma_matrizes(m, b, d, false);
            int** gha = soma_matrizes(m, g, h, true);
            p1 = strassen(m, bds, gha);
            //Calculando P2
            int** ada = soma_matrizes(m, a, d, true);
            int** eha = soma_matrizes(m, e, h, true);
            p2 = strassen(m, ada, eha);

            libera_matriz(m, bds);
            libera_matriz(m, gha);
            libera_matriz(m, ada);
            libera_matriz(m, eha);
            
            MPI_Send(&(p1[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p2[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 2 calcula P3
        if (rank == 2)
        {
            int** acs = soma_matrizes(m, a, c, false);
            int** efa = soma_matrizes(m, e, f, true);
            p3 = strassen(m, acs, efa);
            libera_matriz(m, acs);
            libera_matriz(m, efa);
            MPI_Send(&(p3[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 3 calcula P4
        if (rank == 3)
        {
            int** aba = soma_matrizes(m, a, b, true);
            p4 = strassen(m, aba, h);
            libera_matriz(m, aba);
            MPI_Send(&(p4[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, b);

        //Processo 4 calcula P5
        if (rank == 4)
        {
            int** fhs = soma_matrizes(m, f, h, false);
            p5 = strassen(m, a, fhs);
            libera_matriz(m, fhs);
            MPI_Send(&(p5[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, a);
        libera_matriz(m, f);
        libera_matriz(m, h);

        //Processo 5 calcula P6 e P7
        if (rank == 5)
        {   
            //Calculando P6
            int** ges = soma_matrizes(m, g, e, false);
            p6 = strassen(m, d, ges);
            //Calculando P7
            int** cda = soma_matrizes(m, c, d, true);
            p7 = strassen(m, cda, e);

            libera_matriz(m, ges);
            libera_matriz(m, cda);

            MPI_Send(&(p6[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&(p7[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, g);
        libera_matriz(m, c);
        libera_matriz(m, d);
        libera_matriz(m, e);
    }
    else if(num_process >=8){
        //Processo 0 verifica e recebe resultados dos outros processos
        if (rank == 0)
        {
            cout << "Numero de processos: " << num_process << endl;
            MPI_Recv(&(p1[0][0]), m * m, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p2[0][0]), m * m, MPI_INT, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p3[0][0]), m * m, MPI_INT, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p4[0][0]), m * m, MPI_INT, 4, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p5[0][0]), m * m, MPI_INT, 5, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p6[0][0]), m * m, MPI_INT, 6, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&(p7[0][0]), m * m, MPI_INT, 7, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //processo 1 calcula P1
        if (rank == 1)
        {
            int** bds = soma_matrizes(m, b, d, false);
            int** gha = soma_matrizes(m, g, h, true);
            p1 = strassen(m, bds, gha);
            libera_matriz(m, bds);
            libera_matriz(m, gha);
            MPI_Send(&(p1[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 2 calcula P2
        if (rank == 2)
        {
            int** ada = soma_matrizes(m, a, d, true);
            int** eha = soma_matrizes(m, e, h, true);
            p2 = strassen(m, ada, eha);
            libera_matriz(m, ada);
            libera_matriz(m, eha);
            MPI_Send(&(p2[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 3 calcula P3
        if (rank == 3)
        {
            int** acs = soma_matrizes(m, a, c, false);
            int** efa = soma_matrizes(m, e, f, true);
            p3 = strassen(m, acs, efa);
            libera_matriz(m, acs);
            libera_matriz(m, efa);
            MPI_Send(&(p3[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }

        //processo 4 calcula P4
        if (rank == 4)
        {
            int** aba = soma_matrizes(m, a, b, true);
            p4 = strassen(m, aba, h);
            libera_matriz(m, aba);
            MPI_Send(&(p4[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, b);

        //Processo 5 calcula P5
        if (rank == 5)
        {
            int** fhs = soma_matrizes(m, f, h, false);
            p5 = strassen(m, a, fhs);
            libera_matriz(m, fhs);
            MPI_Send(&(p5[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, a);
        libera_matriz(m, f);
        libera_matriz(m, h);

        //Processo 6 calcula P6
        if (rank == 6)
        {
            int** ges = soma_matrizes(m, g, e, false);
            p6 = strassen(m, d, ges);
            libera_matriz(m, ges);
            MPI_Send(&(p6[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, g);

        //Processo 7 calcula P7
        if (rank == 7)
        {
            int** cda = soma_matrizes(m, c, d, true);
            p7 = strassen(m, cda, e);
            libera_matriz(m, cda);
            MPI_Send(&(p7[0][0]), m * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        libera_matriz(m, c);
        libera_matriz(m, d);
        libera_matriz(m, e);

    }
    else{
        cout << "Quantidade de processadores incompatível" << endl;
        abort;
    }

    
    MPI_Barrier(MPI_COMM_WORLD); //BARREIRA DE SINCRONIA

    //Após obter resultados, P0 finaliza
    if (rank == 0)
    {
        //Criando as submatrizes C11, C12, C21, C22
        int** s1s2a = soma_matrizes(m, p1, p2, true);
        int** s6s4s = soma_matrizes(m, p6, p4, false);
        int** c11 = soma_matrizes(m, s1s2a, s6s4s, true);
        libera_matriz(m, s1s2a);
        libera_matriz(m, s6s4s);

        int** c12 = soma_matrizes(m, p4, p5, true);

        int** c21 = soma_matrizes(m, p6, p7, true);

        int** s2s3s = soma_matrizes(m, p2, p3, false);
        int** s5s7s = soma_matrizes(m, p5, p7, false);
        int** c22 = soma_matrizes(m, s2s3s, s5s7s, true);
        libera_matriz(m, s2s3s);
        libera_matriz(m, s5s7s);
        
        // Submatrizes C criadas agora serão combinadas para gerar a matriz produto resultante da multiplicação
        prod = combina_matrizes(m, c11, c12, c21, c22);

        libera_matriz(m, c11);
        libera_matriz(m, c12);
        libera_matriz(m, c21);
        libera_matriz(m, c22);
    }

    libera_matriz(m, p1);
    libera_matriz(m, p2);
    libera_matriz(m, p3);
    libera_matriz(m, p4);
    libera_matriz(m, p5);
    libera_matriz(m, p6);
    libera_matriz(m, p7);
}

//Função que checa se a multiplicação está correta
//Usada para testes
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

//Procedimento principal
int main(int argc, char* argv[])
{
    int p_rank;
    int num_process;

    //Inicializando ambiente MPI
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        printf("MPI-INIT Failed\n");
        return 0;
    }

    //Definindo rank e numero de processos
    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_process);

    int n;
    //definindo n
    if (p_rank == 0)
    {
        cout << endl;
        cout << "Inserir dimensão n da matriz: ";
        cin >> n;
    }

    //aguardar definição de n, enviar n para todos os processos
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //alocar matrizes
    int** mat1 = alocar_matriz(n);
    int** mat2 = alocar_matriz(n);

    if (p_rank == 0)
    {   
        //preencher matrizes
        preenche_matriz(n, mat1, 2);
        preenche_matriz(n, mat2, 4);

        if(n<32){
            cout << "\nImprimindo matriz A:\n";
            imprime_matriz(n, mat1);
            cout << "\nImprimindo matriz B:\n";
            imprime_matriz(n, mat2);
        }

    }

    //enviar matrizes para todos os processos
    MPI_Bcast(&(mat1[0][0]), n * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(mat2[0][0]), n * n, MPI_INT, 0, MPI_COMM_WORLD);

    //obter tempo inicial
    double startTime = MPI_Wtime();

    //realizar multiplicação
    int** prod;
    strassen(n, mat1, mat2, prod, p_rank, num_process);

    //obter tempo final
    double endTime = MPI_Wtime();

    if (p_rank == 0)
    {
        cout << "\nTempo de execução do Strassen Paralelo (MPI): ";
        cout << setprecision(5) << endTime - startTime << endl;
        cout << endl;

        if(n<32){
            cout << "\nImprimindo matriz C:\n";
            imprime_matriz(n, prod);
        }
        
    }
    

    MPI_Finalize();

    return 0;
}