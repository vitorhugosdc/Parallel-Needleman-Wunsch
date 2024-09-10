#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>


#define A 0 // representa uma base Adenina
#define T 1 // representa uma base Timina
#define G 2 // representa uma base Guanina
#define C 3 // representa uma base Citosina
#define X 4 // representa um gap

#define sair 11

#define maxSeq 10000

int rank, size;

char mapaBases[5] = {'A', 'T', 'G', 'C', '-'};

int seqMaior[maxSeq] = {1, 2, 3, 1, 2, 3, 0, 1, 3, 2, 0, 2};  // TGCTGCATCGAG
int seqMenor[maxSeq] = {1, 2, 3, 0, 1, 2, 0, 2};  // TGCATGAG


int alinhaGMaior[maxSeq],
    alinhaGMenor[maxSeq];

int matrizEscores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior = 12,
    tamSeqMenor = 8,
    tamAlinha,
    penalGap = 0,
    grauMuta = 0,
    diagEscore,
    linEscore,
    colEscore;

int matrizPesos[4][4] = {
    {2, -1, -1, -1},
    {-1, 2, -1, -1},
    {-1, -1, 2, -1},
    {-1, -1, -1, 2}
};

/* leitura do tamanho da sequencia maior */
void leTamMaior(void) {
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do {
        printf("\nDigite 0 < valor < %d = ", maxSeq);
        scanf("%d", &tamSeqMaior);
    } while (tamSeqMaior < 1 || tamSeqMaior > maxSeq);
}

/* leitura do tamanho da sequencia menor */
void leTamMenor(void) {
    printf("\nLeitura do Tamanho da Sequencia Menor:");
    do {
        printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
        scanf("%d", &tamSeqMenor);
    } while (tamSeqMenor < 1 || tamSeqMenor > tamSeqMaior);
}

/* leitura do valor da penalidade de gap */
void lePenalidade(void) {
    printf("\nLeitura da Penalidade de Gap:");
    do {
        printf("\nDigite valor >= 0 = ");
        scanf("%d", &penalGap);
    } while (penalGap < 0);
}

/* leitura da matriz de pesos */
void leMatrizPesos(void) {
    int i, j;
    printf("\nLeitura da Matriz de Pesos:\n");
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
            printf("Digite valor %c x %c = ", mapaBases[i], mapaBases[j]);
            scanf("%d", &(matrizPesos[i][j]));
        }
        printf("\n");
    }
}

/* mostra da matriz de pesos */
void mostraMatrizPesos(void) {
    int i, j;
    printf("\nMatriz de Pesos Atual:");
    printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
    for (i = 0; i < 4; i++) {
        printf("%4c", mapaBases[i]);
        for (j = 0; j < 4; j++) {
            printf("%4d", matrizPesos[i][j]);
        }
        printf("\n");
    }
}

void leGrauMutacao(void) {
    int prob;
    printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
    do {
        printf("\nDigite 0 <= valor <= 100 = ");
        scanf("%d", &prob);
    } while (prob < 0 || prob > 100);
    grauMuta = prob;
}

void leSequencias(void) {
    int i, erro;
    char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];

    printf("\nLeitura das Sequencias:\n");

    do {
        printf("\nPara a Sequencia Maior,");
        printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
        do {
            printf("\n> ");
            fgets(seqMaiorAux, maxSeq, stdin);
            tamSeqMaior = strlen(seqMaiorAux) - 1;
        } while (tamSeqMaior < 1);
        i = 0;
        erro = 0;
        do {
            switch (seqMaiorAux[i]) {
                case 'A': seqMaior[i] = A; break;
                case 'T': seqMaior[i] = T; break;
                case 'G': seqMaior[i] = G; break;
                case 'C': seqMaior[i] = C; break;
                default: erro = 1; break;
            }
            i++;
        } while (erro == 0 && i < tamSeqMaior);
    } while (erro == 1);

    do {
        printf("\nPara a Sequencia Menor, ");
        printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
        do {
            printf("\n> ");
            fgets(seqMenorAux, maxSeq, stdin);
            tamSeqMenor = strlen(seqMenorAux) - 1;
        } while (tamSeqMenor < 1 || tamSeqMenor > tamSeqMaior);
        i = 0;
        erro = 0;
        do {
            switch (seqMenorAux[i]) {
                case 'A': seqMenor[i] = A; break;
                case 'T': seqMenor[i] = T; break;
                case 'G': seqMenor[i] = G; break;
                case 'C': seqMenor[i] = C; break;
                default: erro = 1; break;
            }
            i++;
        } while (erro == 0 && i < tamSeqMenor);
    } while (erro == 1);
}

void geraSequencias(void) {
    int i, dif, probAux, ind, nTrocas;

    printf("\nGeracao Aleatoria das Sequencias:\n");

    for (i = 0; i < tamSeqMaior; i++) {
        seqMaior[i] = rand() % 4;
    }

    dif = tamSeqMaior - tamSeqMenor;

    ind = 0;
    if (dif > 0) {
        ind = rand() % dif;
    }

    for (i = 0; i < tamSeqMenor; i++) {
        seqMenor[i] = seqMaior[ind + i];
    }

    i = 0;
    nTrocas = 0;
    while (i < tamSeqMenor && nTrocas < ((grauMuta * tamSeqMenor) / 100)) {
        probAux = rand() % 100 + 1;
        if (probAux <= grauMuta) {
            seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
            nTrocas++;
        }
        i++;
    }

    printf("\nSequencias Geradas: Dif = %d, Ind = %d, NTrocas = %d\n", dif, ind, nTrocas);
}

void leSequenciasArquivo(void) {
    FILE *file;
    char filename[100];
    char seqMaiorAux[maxSeq];
    char seqMenorAux[maxSeq];

    printf("\nDigite o nome do arquivo para leitura das sequencias: ");
    scanf("%s", filename);

    file = fopen(filename, "r");
    if (file == NULL) {
        printf("\nErro ao abrir o arquivo.\n");
        return;
    }

    if (fgets(seqMaiorAux, maxSeq, file) != NULL) {
        seqMaiorAux[strcspn(seqMaiorAux, "\n")] = '\0';
        tamSeqMaior = strlen(seqMaiorAux);

        for (int i = 0; i < tamSeqMaior; i++) {
            switch (seqMaiorAux[i]) {
                case 'A': seqMaior[i] = A; break;
                case 'T': seqMaior[i] = T; break;
                case 'G': seqMaior[i] = G; break;
                case 'C': seqMaior[i] = C; break;
                default: 
                    printf("\nCaractere inválido encontrado na sequência maior: %c\n", seqMaiorAux[i]);
                    fclose(file);
                    return;
            }
        }
    } else {
        printf("\nErro ao ler a sequência maior.\n");
        fclose(file);
        return;
    }

    if (fgets(seqMenorAux, maxSeq, file) != NULL) {
        seqMenorAux[strcspn(seqMenorAux, "\n")] = '\0';
        tamSeqMenor = strlen(seqMenorAux);

        for (int i = 0; i < tamSeqMenor; i++) {
            switch (seqMenorAux[i]) {
                case 'A': seqMenor[i] = A; break;
                case 'T': seqMenor[i] = T; break;
                case 'G': seqMenor[i] = G; break;
                case 'C': seqMenor[i] = C; break;
                default: 
                    printf("\nCaractere inválido encontrado na sequência menor: %c\n", seqMenorAux[i]);
                    fclose(file);
                    return;
            }
        }
    } else {
        printf("\nErro ao ler a sequência menor.\n");
        fclose(file);
        return;
    }

    fclose(file);

    printf("\nSequências lidas com sucesso:\n");
    printf("Sequência Maior (%d): %s\n", tamSeqMaior, seqMaiorAux);
    printf("Sequência Menor (%d): %s\n", tamSeqMenor, seqMenorAux);
}

void mostraSequencias(void) {
    int i;

    printf("\nSequencias Atuais:\n");
    printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
    printf("%c", mapaBases[seqMaior[0]]);
    for (i = 1; i < tamSeqMaior; i++) {
        printf("%c", mapaBases[seqMaior[i]]);
    }
    printf("\n");

    printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
    printf("%c", mapaBases[seqMenor[0]]);
    for (i = 1; i < tamSeqMenor; i++) {
        printf("%c", mapaBases[seqMenor[i]]);
    }
    printf("\n");
}

void geraMatrizEscores()
{
    int tamBloco;
    int lin, col, peso;
    MPI_Status status;
    MPI_Request request;

    if (rank == 0) {
        printf("Insira o tamanho do bloco: ");
        scanf("%d", &tamBloco);
        for (int i = 1; i < size; i++) {
            MPI_Send(&tamBloco, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&tamBloco, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    if (rank == 0){
        for (col = 0; col <= tamSeqMaior; col++){
            matrizEscores[0][col] = -col * penalGap;
        }

        for (lin = 1; lin <= tamSeqMenor; lin++){
            matrizEscores[lin][0] = -lin * penalGap;
        }

        for (int i = 1; i < size; i++){
            MPI_Send(matrizEscores[0], tamSeqMaior + 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            for (lin = 1; lin <= tamSeqMenor; lin++){
                MPI_Send(&matrizEscores[lin][0], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        }
    }
    else{
        MPI_Recv(matrizEscores[0], tamSeqMaior + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        for (lin = 1; lin <= tamSeqMenor; lin++){
            MPI_Recv(&matrizEscores[lin][0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
    }

    if (rank > 0){
        for (lin = rank; lin <= tamSeqMenor; lin += size - 1){
            for (col = 1; col <= tamSeqMaior; col++){
                if (col % tamBloco == 1 && lin != 1){
                    int bloco;
                    if(col + tamBloco <= tamSeqMaior){
                        bloco = tamBloco;
                    }else{
                        bloco = tamSeqMaior - col + 1;
                    }

                    int proc;
                    if(rank == 1){
                        proc = size - 1;
                    }else{
                        proc = rank -1;
                    }

                    MPI_Recv(&matrizEscores[lin - 1][col], bloco, MPI_INT, proc, 0, MPI_COMM_WORLD, &status);
                }

                peso = matrizPesos[seqMenor[lin - 1]][seqMaior[col - 1]];
                int diagEscore = matrizEscores[lin - 1][col - 1] + peso;
                int linEscore = matrizEscores[lin][col - 1] - penalGap;
                int colEscore = matrizEscores[lin - 1][col] - penalGap;

                if (diagEscore >= linEscore && diagEscore >= colEscore) {
                    matrizEscores[lin][col] = diagEscore;
                } else if (linEscore > colEscore) {
                    matrizEscores[lin][col] = linEscore;
                } else {
                    matrizEscores[lin][col] = colEscore;
                }

                if (col % tamBloco == 0 || col == tamSeqMaior){
                    int bloco;
                    if(col % tamBloco == 0){
                        bloco = tamBloco;
                    }else{
                        bloco = (tamSeqMaior - col + tamBloco - 1) % tamBloco;
                    }
                    MPI_Isend(&matrizEscores[lin][col - bloco + 1], bloco, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

                    if (rank < size - 1){
                        MPI_Isend(&matrizEscores[lin][col - bloco + 1], bloco, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &request);
                    }
                    else if (rank == size - 1){
                        MPI_Isend(&matrizEscores[lin][col - bloco + 1], bloco, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
                    }
                }
            }
        }
    }

    if (rank == 0){
        for (int proc = 1; proc < size; proc++){
            for (lin = proc; lin <= tamSeqMenor; lin += size - 1){
                for (int b = 0; b < tamSeqMaior; b += tamBloco){

                    int bloco;
                    if(b + tamBloco <= tamSeqMaior){
                        bloco = tamBloco;
                    }else{
                        bloco = tamSeqMaior - b;
                    }
                    MPI_Recv(&matrizEscores[lin][b + 1], bloco, MPI_INT, proc, 0, MPI_COMM_WORLD, &status);
                }
            }
        }
        int linMaior = 1, colMaior = 1, maior = matrizEscores[1][1];
        for (lin = 1; lin <= tamSeqMenor; lin++) {
            for (col = 1; col <= tamSeqMaior; col++) {
                if (maior <= matrizEscores[lin][col]) {
                    linMaior = lin;
                    colMaior = col;
                    maior = matrizEscores[lin][col];
                }
            }
        }
        printf("\nMatriz de escores Gerada.");
        printf("\nUltimo Maior escore = %d na celula [%d,%d]",  maior, linMaior, colMaior);
    }
}

/* salva a matriz de escores em um arquivo */
void salvaMatrizEscores(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("\nErro ao abrir o arquivo para escrita.\n");
        return;
    }

    fprintf(file, "%4c%4c", ' ', ' ');
    for (int i = 0; i < tamSeqMaior; i++) {
        fprintf(file, "%4c", mapaBases[seqMaior[i]]);
    }
    fprintf(file, "\n");

    for (int lin = 0; lin <= tamSeqMenor; lin++) {
        if (lin == 0)
            fprintf(file, "%4c", '-');
        else
            fprintf(file, "%4c", mapaBases[seqMenor[lin - 1]]);
        for (int col = 0; col <= tamSeqMaior; col++) {
            fprintf(file, "%4d", matrizEscores[lin][col]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    printf("\nMatriz de escores salva em %s.\n", filename);
}

void mostraMatrizEscores(void) {
    int i, lin, col;

    printf("\nMatriz de escores Atual:\n");

    printf("%4c%4c", ' ', ' ');
    for (i = 0; i <= tamSeqMaior; i++) {
        printf("%4d", i);
    }
    printf("\n");

    printf("%4c%4c%4c", ' ', ' ', '-');
    for (i = 0; i < tamSeqMaior; i++) {
        printf("%4c", mapaBases[seqMaior[i]]);
    }
    printf("\n");

    printf("%4c%4c", '0', '-');
    for (col = 0; col <= tamSeqMaior; col++) {
        printf("%4d", matrizEscores[0][col]);
    }
    printf("\n");

    for (lin = 1; lin <= tamSeqMenor; lin++) {
        printf("%4d%4c", lin, mapaBases[seqMenor[lin - 1]]);
        for (col = 0; col <= tamSeqMaior; col++) {
            printf("%4d", matrizEscores[lin][col]);
        }
        printf("\n");
    }
}

void traceBack()
{ int tbLin, tbCol, peso, pos, posAux, aux, i;

  printf("\nGeracao do Alinhamento Global:\n");
  tbLin=tamSeqMenor;
  tbCol=tamSeqMaior;
  pos=0;
  do
  {
    peso=matrizPesos[(seqMenor[tbLin-1])][(seqMaior[tbCol-1])];
    diagEscore = matrizEscores[tbLin-1][tbCol-1]+peso;
    linEscore = matrizEscores[tbLin][tbCol-1]-penalGap;
    colEscore = matrizEscores[tbLin-1][tbCol]-penalGap;

      if ((diagEscore>=linEscore)&&(diagEscore>=colEscore))
      {
        alinhaGMenor[pos]=seqMenor[tbLin-1];
        alinhaGMaior[pos]=seqMaior[tbCol-1];
        tbLin--;
        tbCol--;
        pos++;
      }
      else if (linEscore>colEscore)
           {
              alinhaGMenor[pos]=X;
              alinhaGMaior[pos]=seqMaior[tbCol-1];
              tbCol--;
              pos++;
           }
           else
           {
              alinhaGMenor[pos]=seqMenor[tbLin-1];
              alinhaGMaior[pos]=X;
              tbLin--;
              pos++;
           }

  } while ((tbLin!=0)&&(tbCol!=0));

  /* descarrega o restante de gaps da linha 0, se for o caso */
  while (tbLin>0)
  {
    alinhaGMenor[pos]=seqMenor[tbLin-1];
    alinhaGMaior[pos]=X;
    tbLin--;
    pos++;
  }

  /* descarrega o restante de gaps da coluna 0, se for o caso */
  while (tbCol>0)
  {
    alinhaGMenor[pos]=X;
    alinhaGMaior[pos]=seqMaior[tbCol-1];
    tbCol--;
    pos++;
  }

  tamAlinha=pos;

  /* o alinhamento foi feito de tras para frente e deve ser
     invertido, conforme segue */
  for (i=0;i<(tamAlinha/2);i++)
  {
    aux=alinhaGMenor[i];
    alinhaGMenor[i]=alinhaGMenor[tamAlinha-i-1];
    alinhaGMenor[tamAlinha-i-1]=aux;

    aux=alinhaGMaior[i];
    alinhaGMaior[i]=alinhaGMaior[tamAlinha-i-1];
    alinhaGMaior[tamAlinha-i-1]=aux;
  }

  printf("\nAlinhamento Global Gerado.");
}

void mostraAlinhamentoGlobal(void)
{   int i;

  printf("\nAlinhamentos Atuais - Tamanho = %d:\n", tamAlinha);

  printf("%c",mapaBases[alinhaGMaior[0]]);
  for (i=1; i<tamAlinha; i++)
    printf("%c",mapaBases[alinhaGMaior[i]]);
  printf("\n");

  printf("%c",mapaBases[alinhaGMenor[0]]);
  for (i=1; i<tamAlinha; i++)
    printf("%c",mapaBases[alinhaGMenor[i]]);
  printf("\n");
}

int menuOpcao(void) {
    int op;
    char enter;

    do {
        printf("\nMenu de Opcao:");
        printf("\n<01> Ler Matriz de Pesos");
        printf("\n<02> Mostrar Matriz de Pesos");
        printf("\n<03> Ler Penalidade de Gap");
        printf("\n<04> Mostrar Penalidade");
        printf("\n<05> Definir Sequencias Genomicas");
        printf("\n<06> Mostrar Sequencias");
        printf("\n<07> Gerar Matriz de Escores");
        printf("\n<08> Mostrar Matriz de Escores");
        printf("\n<09> Gerar Alinhamento Global");
        printf("\n<10> Mostrar Alinhamento Global");
        printf("\n<11> Salvar Matriz de Escores em Arquivo");
        printf("\n<12> Sair");
        printf("\nDigite a opcao => ");
        scanf("%d", &op);
        scanf("%c", &enter);
    } while (op < 1 || op > sair);

    return (op);
}

void trataOpcao(int op) {
    int resp;
    char enter;

    switch (op) {
        case 1:
            leMatrizPesos();
            break;
        case 2:
            mostraMatrizPesos();
            break;
        case 3:
            lePenalidade();
            break;
        case 4:
            printf("\nPenalidade = %d", penalGap);
            break;
        case 5:
            printf("\nDeseja Definicao: <1>MANUAL, <2>ALEATORIA ou <3>ARQUIVO? = ");
            scanf("%d", &resp);
            scanf("%c", &enter);
            if (resp == 1) {
                leSequencias();
            } else if (resp == 2) {
                leTamMaior();
                leTamMenor();
                leGrauMutacao();
                geraSequencias();
            } else if (resp == 3) {
                leSequenciasArquivo();
            }
            break;
        case 6:
            mostraSequencias();
            break;
        case 7:
            geraMatrizEscores();
            salvaMatrizEscores("matriz_escores.txt");
            break;
        case 8:
            mostraMatrizEscores();
            break;
        case 9:
            traceBack();
            break;
        case 10:
            mostraAlinhamentoGlobal();
            break;
    }
}

int main(int argc, char **argv) {
    int opcao;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    srand(time(NULL));

    if (rank == 0) {
        do {
            printf("\n\nPrograma Needleman-Wunsch Paralelizado com MPI\n");
            opcao = menuOpcao();

            MPI_Bcast(&opcao, 1, MPI_INT, 0, MPI_COMM_WORLD);

            trataOpcao(opcao);

        } while (opcao != sair);
    } else {
        do {
            MPI_Bcast(&opcao, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (opcao == 7) {
                geraMatrizEscores();
            }

        } while (opcao != sair);
    }

    MPI_Finalize();

    return 0;
}
