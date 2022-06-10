#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/resource.h>
#include <malloc.h>

int** cria_matriz(int n, int valor);
void popula_matriz(int n, int valor, int **matriz);
void imprime_matriz(int n, int **matriz);
int** multiplica_matrizes(int n, int **A, int **B);
int** dividir_matriz(int n, int **mat, int lin, int col);
int** multiplica_strassen(int n, int **A, int **B);
int** multiplica_strassen_rec(int n, int **A, int** BB);
int** multiplica_divconq(int n, int **A, int **B);
void compor_matriz(int n, int** A, int** B,int lin,int col);
void Tempo_CPU_Sistema(double *seg_CPU_total, double *seg_sistema_total);

int main() {
    printf ("Comparação de Algoritmos para Multiplicação de Matrizes\n\n");

    //Inicialização de Variaveis
    int n, k, va, vb, i, j, escreve;
    int **a, **b;
    double s_CPU_inicial=0, s_total_inicial=0, s_CPU_final=0, s_total_final=0;

    FILE *arq;
    arq = fopen("resultado.txt", "a+");
    if(arq == NULL){
        printf("Problemas na CRIACAO/ABERTURA do arquivo\n");
        return 0;
    }
    char resultado[100];
    char temp[100];

    escreve = fprintf(arq, "\n\n\n[ ['Tamanho da Matriz', 'Algoritmo Direto', 'Strassen'], \n");
    if(escreve == EOF){
        printf("Erro na Gravacao\n");
    }

    //LOOP PRINCIPAL
    int controle = 1;
    while(controle!=0){
        //Interacao com o usuario: Definindo tamanho e valores presentes na matriz A e B
        printf ("Defina o k para tamanho 2^k x 2^k da matriz: ");
        scanf ("%d", &k);
        printf ("\nDefina o Valor para o popular a matriz A: ");
        scanf ("%d", &va);
        printf ("\nDefina o Valor para o popular a matriz B: ");
        scanf ("%d", &vb);

        //Define n = 2^k
        n = pow(2, k);

        escreve = fprintf(arq, "[' %d x %d ', ",n,n);
        if(escreve == EOF){
            printf("Erro na Gravacao\n");
        }

        //criando e populando matriz A:
        a = cria_matriz(n, va);
        //criando e populando matriz B:
        b = cria_matriz(n, vb);

        //printf("\nImprimindo Matriz A:\n");
        //imprime_matriz(n, a);
        //printf("\nImprimindo Matriz B:\n");
        //imprime_matriz(n, b);

        printf("\nCriando Matriz C = A x B, Calculada por algoritmo Direto de multiplicação\n");
        Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
        int **c = multiplica_matrizes(n, a, b);
        Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
        double temponormal = s_CPU_final - s_CPU_inicial;
        //imprime_matriz(n, c);

        /*
        printf("\nCriando Matriz D = A x B, Calculada por algoritmo Normal de Divisão e Conquista\n");
        Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
        int **d = multiplica_divconq(n, a, b);
        Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
        double tempondivconq = s_CPU_final - s_CPU_inicial;
        //imprime_matriz(n, d);
        */

        printf("\nCriando Matriz E = A x B, Calculada por algoritmo Strassen de multiplicação\n");
        Tempo_CPU_Sistema(&s_CPU_inicial, &s_total_inicial);
        int **e = multiplica_strassen(n, a, b);
        Tempo_CPU_Sistema(&s_CPU_final, &s_total_final);
        double tempostrassen = s_CPU_final - s_CPU_inicial;
        //imprime_matriz(n, e);

        printf("\nTempo de execução do algoritmo Direto: %f", temponormal);
        //printf("\nTempo de execução do algoritmo Divisão e Conquista Comum: %f", tempondivconq);
        printf("\nTempo de execução do algoritmo strassen: %f", tempostrassen);

        escreve = fprintf(arq, " %f , ", temponormal);
        if(escreve == EOF){
            printf("Erro na Gravacao\n");
        }
        /*
        escreve = fprintf(arq, "' %f ', ", tempondivconq);
        if(escreve == EOF){
            printf("Erro na Gravacao\n");
        }
        */
        escreve = fprintf(arq, " %f ", tempostrassen);
        if(escreve == EOF){
            printf("Erro na Gravacao\n");
        }

        escreve = fprintf(arq, "], \n");
        if(escreve == EOF){
            printf("Erro na Gravacao\n");
        }

        fflush(arq);

        printf("\n\nReiniciar procedimento: 1");
        printf("\nEncerrar programa: 0:");
        printf("\nEscolha uma opção:");
        scanf ("%d", &controle);
    }
    escreve = fprintf(arq, "]");
    if(escreve == EOF){
        printf("Erro na Gravacao\n");
    }
    fclose(arq);
    return 0;
}

//Função que cria uma matriz e a preenche com zeros;
//Entradas: o tamanho da matriz nxn a ser criada, o valor para preencher a matriz
//Saída: Ponteiro para nova Matriz gerada
int **cria_matriz(int n, int valor){
    int i;
    int ** novaMatriz;
    // aloca um vetor de n ponteiros para linhas
    novaMatriz = malloc (n * sizeof (int*)) ;
    // aloca cada uma das linhas (vetores de COL inteiros)
    for (i=0; i < n; i++){
        novaMatriz[i] = malloc (n * sizeof (int)) ;
    }
    popula_matriz(n, valor, novaMatriz);
    return novaMatriz;
}

//Função que popula uma matriz com um valor especifico
//Entradas: tamanho da matriz nxn, valor a colocar em todos os elementos, ponteiro para a saida
//Esta função é void, portanto não retorna e não deve ser chamada em modo atribuição
void popula_matriz(int n, int valor, int **matriz){
    int i,j;
    for (i=0; i<n; i++){
        for(j=0; j<n; j++){
            matriz[i][j] = valor;
        }
    }
}

//Função que imprime uma matriz na tela;
//Recebe como parametros o n tamanho da matriz nxn e o ponteiro para a matriz
//Essa função não retorna
void imprime_matriz(int n, int **matriz){
    int i, j;
    for(i=0; i<n; i++){
        printf("|");
        for(j=0; j<n; j++){
            printf("  ");
            printf("%d",matriz[i][j]);
            printf("  ");
            printf("|");
        }
        printf("\n");
    }
}

//Função que multiplica uma matriz por outra matriz
//entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a Matriz B
//Saída: Retorna o ponteiro para uma matriz C, resultado da multiplicação
int** multiplica_matrizes(int n, int **A, int **B){
    int i, j, x, aux;
    //Cria a matriz C resposta e a preenche com zeros
    int **C = cria_matriz(n, 0);
    //inicia a auxiliar usada na multiplicacao
    aux = 0;
    //loop da multiplicacao
    for(i = 0; i < n; i++) {
        for(j = 0; j < n; j++) {
            for(x = 0; x < n; x++) {
                aux +=  A[i][x] * B[x][j];
            }
            C[i][j] = aux;
            aux = 0;
        }
    }
    //retorna matriz resposta C
    return C;
}

//Função que Divide uma matriz
//entradas: tamanho da matriz nxn, ponteiro para a matriz a ser dividida, valor da linha e coluna para obter o quadrante especifico
//Saída: Retorna o ponteiro para uma matriz matRes do quadrante obtido por lin e col, resultado da fragmentacao
int** dividir_matriz(int n, int **mat, int lin, int col) {
    int n2 = n/2;
    //Cria matriz resultado e a povoa com zeros
    int ** matRes = cria_matriz(n2, 0);
    int i,j,l=lin,c=col;
    for(i = 0;i < n2; i++) {
        c=col;
        for(j = 0; j < n2; j++) {
            matRes[i][j] = mat[l][c];
            c++;
        }
        l++;
    }
    return matRes;
}

//Função que Soma matrizes
//entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a Matriz B
//Saída: Retorna o ponteiro para uma matriz C, resultado da soma
int** soma_matriz(int n, int** A, int** B){
    //cria uma matriz C e a preenche com zeros
    int ** C = cria_matriz(n, 0);
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            C[i][j]=A[i][j]+B[i][j];
        }
    }
    return C;
}

//Função que Subtrai matrizes
//entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a Matriz B
//Saída: Retorna o ponteiro para uma matriz C, resultado da subtração
int** subtrai_matriz(int n, int** A, int** B){
    int ** C = cria_matriz(n, 0);
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            C[i][j]=A[i][j]-B[i][j];
        }
    }
    return C;
}

int** multiplica_divconq(int n, int **A, int **B){
    //Cria Matriz Resposta C e a popula com zeros
    int **C = cria_matriz(n, 0);

    //Se n for maior que 1, divide a matriz
    if(n>64) {
        int ** a11 = dividir_matriz(n, A, 0, 0);
        int ** a12 = dividir_matriz(n, A, 0, (n/2));
        int ** a21 = dividir_matriz(n, A, (n/2), 0);
        int ** a22 = dividir_matriz(n, A, (n/2), (n/2));
        int ** b11 = dividir_matriz(n, B, 0, 0);
        int ** b12 = dividir_matriz(n, B, 0, n/2);
        int ** b21 = dividir_matriz(n, B, n/2, 0);
        int ** b22 = dividir_matriz(n, B, n/2, n/2);

        //Chamada Recursiva para dividir e conquistar
        int ** c11 = soma_matriz(n/2, multiplica_divconq(n/2, a11, b11), multiplica_divconq(n/2, a12, b21));
        int ** c12 = soma_matriz(n/2, multiplica_divconq(n/2, a11, b12), multiplica_divconq(n/2, a12, b22));
        int ** c21 = soma_matriz(n/2, multiplica_divconq(n/2, a21, b11), multiplica_divconq(n/2, a22, b21));
        int ** c22 = soma_matriz(n/2, multiplica_divconq(n/2, a21, b12), multiplica_divconq(n/2, a22, b22));

        //Compor (juntar) as Matrizes
        compor_matriz(n/2, c11, C, 0, 0);
        compor_matriz(n/2, c12, C,0,n/2);
        compor_matriz(n/2, c21, C, n/2, 0);
        compor_matriz(n/2, c22, C, n/2, n/2);
    }
    else {
        multiplica_matrizes(n, A, B);
    }
    return C;
}

//Funcão que multiplica pelo algoritmo de Strassen
//Essa função encapsula a função recursiva do algoritmo
//Entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a matriz B
//Saida: Matriz C, resultado da operação A x B
int ** multiplica_strassen(int n, int **A, int **B){
    int **C = multiplica_strassen_rec(n, A, B);
    return C;
}

//Função que implementa a solução divisão e conquista do algoritmo de Strassen
//Entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a matriz B
//Saida: Matriz C, resultado da operação A x B
int** multiplica_strassen_rec(int n, int **A, int** B){

    //Cria Matriz Resposta C e a popula com zeros
    int **C = cria_matriz(n, 0);

    //Se n for maior que 64, divide a matriz
    if(n>64) {
        int ** a11 = dividir_matriz(n, A, 0, 0);
        int ** a12 = dividir_matriz(n, A, 0, (n/2));
        int ** a21 = dividir_matriz(n, A, (n/2), 0);
        int ** a22 = dividir_matriz(n, A, (n/2), (n/2));
        int ** b11 = dividir_matriz(n, B, 0, 0);
        int ** b12 = dividir_matriz(n, B, 0, n/2);
        int ** b21 = dividir_matriz(n, B, n/2, 0);
        int ** b22 = dividir_matriz(n, B, n/2, n/2);

        //Chamada Recursiva para dividir e Conquistar
        int** P= multiplica_strassen_rec(n/2, soma_matriz(n/2, a11, a22),soma_matriz(n/2, b11, b22));
        int** Q= multiplica_strassen_rec(n/2, soma_matriz(n/2, a21, a22), b11);
        int** R= multiplica_strassen_rec(n/2, a11,subtrai_matriz(n/2, b12, b22));
        int** S= multiplica_strassen_rec(n/2, a22,subtrai_matriz(n/2, b21, b11));
        int** T= multiplica_strassen_rec(n/2, soma_matriz(n/2, a11, a12),b22);
        int** U= multiplica_strassen_rec(n/2,subtrai_matriz(n/2, a21, a11),soma_matriz(n/2, b11, b12));
        int** V= multiplica_strassen_rec(n/2, subtrai_matriz(n/2, a12, a22), soma_matriz(n/2, b21, b22));

        int** c11 = soma_matriz(n/2, subtrai_matriz(n/2, soma_matriz(n/2, P, S), T), V);
        int** c12 = soma_matriz(n/2, R,T);
        int** c21 = soma_matriz(n/2, Q,S);
        int** c22 = soma_matriz(n/2, subtrai_matriz(n/2, soma_matriz(n/2, P, R), Q), U);
        //Compor (juntar) as Matrizes
        compor_matriz(n/2, c11, C, 0, 0);
        compor_matriz(n/2, c12, C,0,n/2);
        compor_matriz(n/2, c21, C, n/2, 0);
        compor_matriz(n/2, c22, C, n/2, n/2);
    }
    else {
        //Condição para fim da recursão.
        multiplica_matrizes(n, A, B);
    }
    return C;
}

//Função realiza a composição de uma sub-matriz a uma matriz maior (inverso da divisão)
//Entradas: tamanho n da matriz nxn, ponteiro para a submatriz A a ser inserida, ponteiro para a matriz B que recebe a submatriz A, posição lin e col do inicio da inserção
//Essa função é void e não retorna
void compor_matriz(int n, int** A, int** B,int lin,int col){
    int i,j,l=lin,c=col;
    for(i = 0; i < n; i++){
        c=col;
        for(j = 0; j < n; j++){
            B[l][c]=A[i][j];
            c++;
        }
        l++;
    }
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