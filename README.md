# Trabalho-LPP

Trabalho prático usado para disciplina Laboratório de Programação Paralela

# Strassen

Sejam A e B matrizes quadradas de ordem 2n x 2n, e seja C o produto dessas matrizes, para calcular esse produto particiona-se A, B e C em quatro submatrizes de mesmo tamanho.
Então, Definem-se então as matrizes:

P1 = (A11 + A22) x (B11 + B22)

P2 = (A21 + A22) x B11

P3 = A11 x (B12 - B22)

P4 = A22 x (B21 - B11)

P5 = (A11 + A12) X B22

P6 = (A21 - A11) x (B11 + B12)

P7 = (A12 - A22) x (B21 + B22)
 
Que serão usadas para expressar Ci,j em termos dos Pk.
Devido a definição das matrizes P, pode-se eliminar uma multiplicação de matrizes e reduzir para 7 a sua quantidade (Uma multiplicação para cada Pk),
expressando assim os Ci,j como:

C11 = P1 + P4 - P5 + P7

C12 = P3 + P5

C21 = P2 + P4

C22 = P1 - P2 + P3 + P6

Ao todo, o Algoritmo de Strassen realiza 7 operações de multiplicação e 18 operações de soma/subtração, essas menos custosas.
Em implementações práticas do método de Strassen, a multiplicação de submatrizes de tamanho suficientemente pequenas é feita pelo método usual,
pois nesses casos, o método usual se mostra mais eficiente.

