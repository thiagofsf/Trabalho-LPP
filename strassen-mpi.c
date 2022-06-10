#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define Nconst 256


// Função Criada para Alocação de memória das matrizes
// Entrada: o tamaho da matriz nxn a ser criada
// Saida: ponteiro para matriz alocada
int **mem_alloc(int size){
	int i;
	int **array = (int **) malloc (size * sizeof (int *));
    for (i=0; i< size; i++) array[i]=malloc(size*sizeof(int));
	return array;
}

//Função que inicializa (aloca) e preenche uma matriz;
//Entradas: o tamanho da matriz nxn a ser criada, o valor para preencher a matriz
//Saída: Ponteiro para a matriz criada
int ** cria_matriz(int n, int valor){
	int ** array = mem_alloc(n);
	int i,j;
	for (i = 0; i < n; i++){
		for (j =0; j<n; j++) array[i][j] = valor;
	}
	return array;
}

//Função que Divide uma matriz
//entradas: tamanho da matriz nxn, ponteiro para a matriz a ser dividida, valor da linha e coluna para obter o quadrante especifico
//Saída: Retorna o ponteiro para uma matriz matRes do quadrante obtido por lin e col, resultado da fragmentacao
int** dividir_matriz(int n, int **mat, int lin, int col){
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
//entradas: ponteiro para a matriz A, ponteiro para a Matriz B, tamanho da matriz nxn
//Saída: Retorna o ponteiro para uma matriz C, resultado da soma
int** soma_matriz(int** A, int** B, int n){
    //Aloca uma matriz C
    int ** C = mem_alloc(n);
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            C[i][j]=A[i][j]+B[i][j];
        }
    }
    return C;
}

//Função que Subtrai matrizes
//entradas: ponteiro para a matriz A, ponteiro para a Matriz B, tamanho da matriz nxn
//Saída: Retorna o ponteiro para uma matriz C, resultado da subtração
int** subtrai_matriz(int** A, int** B, int n){
    int ** C = mem_alloc(n);
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            C[i][j]=A[i][j]-B[i][j];
        }
    }
    return C;
}

//Função serial direta que multiplica uma matriz por outra matriz
//entradas: tamanho da matriz nxn, ponteiro para a matriz A, ponteiro para a Matriz B
//Saída: Retorna o ponteiro para uma matriz C, resultado da multiplicação
int** multiplica_matrizes(int **A, int **B, int n){
    int i, j, x, aux;
    //Cria a matriz C resposta e a preenche com zeros
    int **C = mem_alloc(n);
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

//Função que implementa o algoritmo de Strassen
//Entradas: ponteiro para a matriz A, ponteiro para a matriz B, tamanho da matriz nxn,
//Saida: Matriz C, resultado da operação A x B
int **Strassen(int **matrixA, int **matrixB, int n){
		int **P1,**P2,**P3,**P4,**P5,**P6,**P7;
		int **C11,**C12,**C21,**C22;
		int **S1,**S2,**S3,**S4,**S5,**S6,**S7,**S8,**S9,**S10;
		int **A11,**A12,**A21,**A22,**B11,**B12,**B21,**B22;
		int i,j;
		int ** res = mem_alloc(n);
		int new_n = n/2;
		if(n>64) {
				A11 = mem_alloc(new_n);
				A12 = mem_alloc(new_n);
				A21 = mem_alloc(new_n);
				A22 = mem_alloc(new_n);
				B11 = mem_alloc(new_n);
				B12 = mem_alloc(new_n);
				B21 = mem_alloc(new_n);
				B22 = mem_alloc(new_n);
				for (i = 0; i < new_n; i++) {
					for(j = 0; j<new_n; j++){
						A11[i][j] = matrixA[i][j];
						A12[i][j] = matrixA[i][j + new_n];
						A21[i][j] = matrixA[i + new_n][j];
						A22[i][j] = matrixA[i + new_n][j + new_n];
						B11[i][j] = matrixB[i][j];
						B12[i][j] = matrixB[i][j + new_n];
						B21[i][j] = matrixB[i + new_n][j];
						B22[i][j] = matrixB[i + new_n][j + new_n];
					}
				}
				S1 = mem_alloc(new_n);
				S2 = mem_alloc(new_n);
				S3 = mem_alloc(new_n);
				S4 = mem_alloc(new_n);
				S5 = mem_alloc(new_n);
				S6 = mem_alloc(new_n);
				S7 = mem_alloc(new_n);
				S8 = mem_alloc(new_n);
				S9 = mem_alloc(new_n);
				S10 = mem_alloc(new_n);
				for (i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
						S1[i][j] = B12[i][j] - B22[i][j];
						S2[i][j] = A11[i][j] + A12[i][j];
						S3[i][j] = A21[i][j] + A22[i][j];
						S4[i][j] = B21[i][j] - B11[i][j];
						S5[i][j] = A11[i][j] + A22[i][j];
						S6[i][j] = B11[i][j] + B22[i][j];
						S7[i][j] = A12[i][j] - A22[i][j];
						S8[i][j] = B21[i][j] + B22[i][j];
						S9[i][j] = A11[i][j] - A21[i][j];
						S10[i][j] = B11[i][j] + B12[i][j];
					}
				}
				P1 = Strassen(A11,S1,new_n);
				P2 = Strassen(S2,B22,new_n);
				P3 = Strassen(S3,B11,new_n);
				P4 = Strassen(A22,S4,new_n);
				P5 = Strassen(S5,S6,new_n);
				P6 = Strassen(S7,S8,new_n);
				P7 = Strassen(S9,S10,new_n);
				
                //Liberar espaço de memória não necessário
				free(A11);free(A12);free(A21);free(A22);free(B11);free(B12);free(B21);free(B22);
				free(S1);free(S2);free(S3);free(S4);free(S5);free(S6);free(S7);free(S8);free(S9);free(S10);

				C11 = mem_alloc(new_n);
				C12 = mem_alloc(new_n);
				C21 = mem_alloc(new_n);
				C22 = mem_alloc(new_n);
				for(i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
						C11[i][j] = P5[i][j] + P4[i][j] - P2[i][j] + P6[i][j];
						C12[i][j] = P1[i][j] + P2[i][j];
						C21[i][j] = P3[i][j] + P4[i][j];
						C22[i][j] = P5[i][j] + P1[i][j] - P3[i][j] - P7[i][j];
					}
				}
                
                // Liberando espaço de memória não necessário
				free(P1);free(P2);free(P3);free(P4);free(P5);free(P6);free(P7);

				for (i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
					res[i][j] = C11[i][j];
					res[i][j+new_n] = C12[i][j];
					res[new_n+i][j] = C21[i][j];
					res[new_n+i][new_n+j] = C22[i][j];
					}
				}
		    }
		else {
            multiplica_matrizes(matrixA, matrixB, n);
        }
		return res;
}

int main(int argc, char *argv[]) {
            
            //Definição de variáveis da interação
            int n, va, vb;
            n = 128;
            va= 2;
            vb= 4;
            int new_n = n/2,rank,comm_size;
			int i,j;
			srand(0);

			//Inicializando as Matrizes A e B
			int **A = cria_matriz(n, va);
			int **B = cria_matriz(n, vb);

			//Alocando Memória para matrizes do processo
			int **P1 = mem_alloc(new_n);
			int **P2 = mem_alloc(new_n);
			int **P3 = mem_alloc(new_n);
			int **P4 = mem_alloc(new_n);
			int **P5 = mem_alloc(new_n);
			int **P6 = mem_alloc(new_n);
			int **P7 = mem_alloc(new_n);
			int **C11 = mem_alloc(new_n);
			int **C12 = mem_alloc(new_n);
			int **C21 = mem_alloc(new_n);
			int **C22 = mem_alloc(new_n);
			int **C_parallel = mem_alloc(n);
			float parallel_start;
			
            //MPI
			MPI_Init(&argc,&argv);
			MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);

		  //Iniciando Matrizes
		if (rank == 0){
			parallel_start = MPI_Wtime();
			for (i=1;i<comm_size;i++){
				// Enviando as matrizes A,B para todos os nós
				MPI_Send(&(A[0][0]),(n*n),MPI_INT,i,1,MPI_COMM_WORLD);
				MPI_Send(&(B[0][0]),(n*n),MPI_INT,i,2,MPI_COMM_WORLD);
			}

			//Criar submatrizes para raizes
			int **A11 = mem_alloc(new_n);
			int **A12 = mem_alloc(new_n);
			int **A21 = mem_alloc(new_n);
			int **A22 = mem_alloc(new_n);
			int **B11 = mem_alloc(new_n);
			int **B12 = mem_alloc(new_n);
			int **B21 = mem_alloc(new_n);
			int **B22 = mem_alloc(new_n);

			for (i = 0; i < new_n; i++) {
				for(j = 0; j < new_n; j++){
					A11[i][j] = A[i][j];
					A12[i][j] = A[i][j + new_n];
					A21[i][j] = A[i + new_n][j];
					A22[i][j] = A[i + new_n][j + new_n];
					B11[i][j] = B[i][j];
					B12[i][j] = B[i][j + new_n];
					B21[i][j] = B[i + new_n][j];
					B22[i][j] = B[i + new_n][j + new_n];
				}
			}

			free(A);free(B);

			// Dependendo do número de raizes, calcular adequadamente
			// Si, Pi
			if (comm_size == 2){
				int **S1 = mem_alloc(new_n);
				int **S2 = mem_alloc(new_n);
				int **S3 = mem_alloc(new_n);
				int **S4 = mem_alloc(new_n);
				for (i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
						S1[i][j] = B12[i][j] - B22[i][j];
						S2[i][j] = A11[i][j] + A12[i][j];
						S3[i][j] = A21[i][j] + A22[i][j];
						S4[i][j] = B21[i][j] - B11[i][j];
					}
				}
				P1 = Strassen(A11,S1,new_n);
				P2 = Strassen(S2,B22,new_n);
				P3 = Strassen(S3,B11,new_n);
				P4 = Strassen(A22,S4,new_n);
				free(A11);free(A12);free(A21);free(A22);
				free(B11);free(B12);free(B21);free(B22);
				free(S1);free(S2);free(S3);free(S4);
			}
			else if (comm_size == 4){
				int **S1 = mem_alloc(new_n);
				int **S2 = mem_alloc(new_n);
				for (i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
						S1[i][j] = B12[i][j] - B22[i][j];
						S2[i][j] = A11[i][j] + A12[i][j];
					}
				}
			P1 = Strassen(A11,S1,new_n);
			P2 = Strassen(S2,B22,new_n);
			free(A11);free(A12);free(A21);free(A22);
			free(B11);free(B12);free(B21);free(B22);
			free(S1);free(S2);
			}
			else if (comm_size > 4){
				int **S1 = mem_alloc(new_n);
				free(A12);free(A21);free(A22);
				free(B11);free(B21);
				for (i=0;i<new_n;i++){
					for(j=0;j<new_n;j++){
						S1[i][j] = B12[i][j] - B22[i][j];
					}
				}
				P1 = Strassen(A11,S1,new_n);
				free(A11);free(S1);
				free(B12);free(B22);
			}
		}
		else{	
                // Alocar memoria para A,B para todos os nós
				int **local_A = mem_alloc(n);
				int **local_B = mem_alloc(n);
				MPI_Recv(&(local_A[0][0]),(n*n),MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(local_B[0][0]),(n*n),MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				
                // Criar cópias locais das submatrizes de A e B
				int **local_A11 = mem_alloc(new_n);
				int **local_A12 = mem_alloc(new_n);
				int **local_A21 = mem_alloc(new_n);
				int **local_A22 = mem_alloc(new_n);
				int **local_B11 = mem_alloc(new_n);
				int **local_B12 = mem_alloc(new_n);
				int **local_B21 = mem_alloc(new_n);
				int **local_B22 = mem_alloc(new_n);
				for (i = 0; i < new_n; i++) {
					for(j = 0; j<new_n; j++){
						local_A11[i][j] = local_A[i][j];
						local_A12[i][j] = local_A[i][j + new_n];
						local_A21[i][j] = local_A[i + new_n][j];
						local_A22[i][j] = local_A[i + new_n][j + new_n];
						local_B11[i][j] = local_B[i][j];
						local_B12[i][j] = local_B[i][j + new_n];
						local_B21[i][j] = local_B[i + new_n][j];
						local_B22[i][j] = local_B[i + new_n][j + new_n];
					}
				}
				free(local_A);free(local_B);
				// Novamente, de acordo com o número de nós, calcular adequadamente
				if (comm_size == 2){
						int **local_S5 = soma_matriz(local_A11,local_A22,new_n);
						int **local_S6 = soma_matriz(local_B11,local_B22,new_n);
						int **local_S7 = subtrai_matriz(local_B21,local_B22,new_n);
						int **local_S8 = soma_matriz(local_B21,local_B22,new_n);
						int **local_S9 = subtrai_matriz(local_A11,local_A21,new_n);
						int **local_S10 = soma_matriz(local_B11,local_B12,new_n);
						free(local_A11);free(local_A12);free(local_A21);free(local_A22);
						free(local_B11);free(local_B12);free(local_B21);free(local_B22);
						int **local_P5 = Strassen(local_S5,local_S6,new_n);
						int **local_P6 = Strassen(local_S7,local_S8,new_n);
						int **local_P7 = Strassen(local_S9,local_S10,new_n);
						// Mandar os coeficientes Pi calculados a raiz
						MPI_Send(&(local_P5[0][0]),(new_n*new_n),MPI_INT,0,5,MPI_COMM_WORLD);
						MPI_Send(&(local_P6[0][0]),(new_n*new_n),MPI_INT,0,6,MPI_COMM_WORLD);
						MPI_Send(&(local_P7[0][0]),(new_n*new_n),MPI_INT,0,7,MPI_COMM_WORLD);
						free(local_S5);free(local_S6);free(local_S7);free(local_S8);
						free(local_P5);free(local_P6);free(local_P7);
				}
				else if (comm_size == 4){
					if (rank == 1){
							int **local_S3 = soma_matriz(local_A21,local_A22,new_n);
							int **local_S4 = subtrai_matriz(local_B21,local_B11,new_n);
							int **local_P3 = Strassen(local_S3,local_B11,new_n);
							int **local_P4 = Strassen(local_A22,local_S4,new_n);
							free(local_A11);free(local_A12);free(local_A21);free(local_A22);
							free(local_B11);free(local_B12);free(local_B21);free(local_B22);
							free(local_S3);free(local_S4);
							// Mandar os coeficientes Pi calculados a raiz
							MPI_Send(&(local_P3[0][0]),(new_n*new_n),MPI_INT,0,3,MPI_COMM_WORLD);
							MPI_Send(&(local_P4[0][0]),(new_n*new_n),MPI_INT,0,4,MPI_COMM_WORLD);
							free(local_P3);free(local_P4);
					}
					else if(rank == 2){
							int **local_S5 = soma_matriz(local_A11,local_A22,new_n);
							int **local_S6 = soma_matriz(local_B11,local_B22,new_n);
							int **local_S7 = subtrai_matriz(local_A12,local_A22,new_n);
							int **local_S8 = soma_matriz(local_B21,local_B22,new_n);
							free(local_A11);free(local_A12);free(local_A21);free(local_A22);
							free(local_B11);free(local_B12);free(local_B21);free(local_B22);
							int **local_P5 = Strassen(local_S5,local_S6,new_n);
							int **local_P6 = Strassen(local_S7,local_S8,new_n);
							free(local_S5);free(local_S6);free(local_S7);free(local_S8);
							// Mandar os coeficientes Pi calculados a raiz
							MPI_Send(&(local_P5[0][0]),(new_n*new_n),MPI_INT,0,5,MPI_COMM_WORLD);
							MPI_Send(&(local_P6[0][0]),(new_n*new_n),MPI_INT,0,6,MPI_COMM_WORLD);
							free(local_P5);free(local_P6);
					}
					else if(rank == 3){
							int **local_S9 = subtrai_matriz(local_A11,local_A21,new_n);
							int **local_S10 = soma_matriz(local_B11,local_B12,new_n);
							int **local_P7 = Strassen(local_S9,local_S10,new_n);
							// Mandar os coeficientes Pi calculados a raiz
							MPI_Send(&(local_P7[0][0]),(new_n*new_n),MPI_INT,0,7,MPI_COMM_WORLD);
							free(local_S9);free(local_S10);free(local_P7);
					}
				}
				else if (comm_size > 4){ 
                    // Cada nó calcula um Pi e envia a raiz
					if (rank == 1){
							int **local_S2 = soma_matriz(local_A11,local_A12,new_n);
							int **local_P2 = Strassen(local_S2,local_B22,new_n);
							MPI_Send(&(local_P2[0][0]),(new_n*new_n),MPI_INT,0,2,MPI_COMM_WORLD);
					}
					else if (rank == 2){
							int **local_S3 = soma_matriz(local_A21,local_A22,new_n);
							int **local_P3 = Strassen(local_S3,local_B11,new_n);
							MPI_Send(&(local_P3[0][0]),(new_n*new_n),MPI_INT,0,3,MPI_COMM_WORLD);
					}
					else if (rank == 3){
							int **local_S4 = subtrai_matriz(local_B21,local_B11,new_n);
							int **local_P4 = Strassen(local_A22,local_S4,new_n);
							MPI_Send(&(local_P4[0][0]),(new_n*new_n),MPI_INT,0,4,MPI_COMM_WORLD);
					}
					else if(rank == 4){
							int **local_S5 = soma_matriz(local_A11,local_A22,new_n);
							int **local_S6 = soma_matriz(local_B11,local_B22,new_n);
							int **local_P5 = Strassen(local_S5,local_S6,new_n);
							MPI_Send(&(local_P5[0][0]),(new_n*new_n),MPI_INT,0,5,MPI_COMM_WORLD);
					}
					else if(rank == 5){
							int **local_S7 = subtrai_matriz(local_A12,local_A22,new_n);
							int **local_S8 = soma_matriz(local_B21,local_B22,new_n);
							int **local_P6 = Strassen(local_S7,local_S8,new_n);
							MPI_Send(&(local_P6[0][0]),(new_n*new_n),MPI_INT,0,6,MPI_COMM_WORLD);
					}
					else if(rank == 6){
							int **local_S9 = subtrai_matriz(local_A11,local_A21,new_n);
							int **local_S10 = soma_matriz(local_B11,local_B12,new_n);
							int **local_P7 = Strassen(local_S9,local_S10,new_n);
							MPI_Send(&(local_P7[0][0]),(new_n*new_n),MPI_INT,0,7,MPI_COMM_WORLD);
					}
				}
			}
		if (rank == 0){
			// Dependendo do número de nós, receber adequadamente os coeficientes Pi
			if (comm_size == 2){
				MPI_Recv(&(P5[0][0]),(new_n*new_n),MPI_INT,1,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P6[0][0]),(new_n*new_n),MPI_INT,1,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P7[0][0]),(new_n*new_n),MPI_INT,1,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			else if (comm_size == 4){
				MPI_Recv(&(P3[0][0]),(new_n*new_n),MPI_INT,1,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P4[0][0]),(new_n*new_n),MPI_INT,1,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P5[0][0]),(new_n*new_n),MPI_INT,2,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P6[0][0]),(new_n*new_n),MPI_INT,2,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P7[0][0]),(new_n*new_n),MPI_INT,3,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			else if (comm_size > 4){
				MPI_Recv(&(P2[0][0]),(new_n*new_n),MPI_INT,1,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P3[0][0]),(new_n*new_n),MPI_INT,2,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P4[0][0]),(new_n*new_n),MPI_INT,3,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P5[0][0]),(new_n*new_n),MPI_INT,4,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P6[0][0]),(new_n*new_n),MPI_INT,5,6,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				MPI_Recv(&(P7[0][0]),(new_n*new_n),MPI_INT,6,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			// Calcular as submatrizes C(i,j)
			for (i = 0; i<new_n; i++){
				for (j = 0;j<new_n; j++){
				C11[i][j] = P5[i][j] + P4[i][j] - P2[i][j] + P6[i][j];
				C12[i][j] = P1[i][j] + P2[i][j];
				C21[i][j] = P3[i][j] + P4[i][j];
				C22[i][j] = P5[i][j] + P1[i][j] - P3[i][j] - P7[i][j];
				}
			}
			// Criar uma matriz C, formata pelas submatrizes C(i,j)
			for (i=0;i<new_n;i++){
				for(j=0;j<new_n;j++){
					C_parallel[i][j] = C11[i][j];
					C_parallel[i][j+new_n] = C12[i][j];
					C_parallel[new_n+i][j] = C21[i][j];
					C_parallel[new_n+i][new_n+j] = C22[i][j];
				}
			}
			float parallel_end = MPI_Wtime();
			printf("With %d core(s), parallel calcuation time is %f \n",comm_size,parallel_end - parallel_start);
		}
		MPI_Finalize();
}