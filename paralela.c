#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#define A 0 // representa uma base Adenina
#define T 1 // representa uma base Timina
#define G 2 // representa uma base Guanina
#define C 3 // representa uma base Citosina
#define X 4 // representa um gap

#define sair 13

#define maxSeq 10000

char mapaBases[5] = {'A', 'T', 'G', 'C', '-'};

int seqMaior[maxSeq],
    seqMenor[maxSeq];

int alinhaGMaior[10][maxSeq],
    alinhaGMenor[10][maxSeq];

int matrizEscores[maxSeq + 1][maxSeq + 1];

int tamSeqMaior,
    tamSeqMenor,
    tamAlinha[10], // Armazena o tamanho de cada alinhamento
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

int numAlignments = 1;
int numThreads = 1;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int thread_id;
    int alignment_index;
    int start_lin;
    int start_col;
    int end_lin;
} thread_data_t;

int multiple_paths_possible(int lin, int col);

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
    printf("\nDigite o nome do arquivo para leitura das sequencias: ");
    scanf("%s", filename);

    file = fopen(filename, "r");
    if (file == NULL) {
        printf("\nErro ao abrir o arquivo.\n");
        return;
    }

    fscanf(file, "%d %d", &tamSeqMaior, &tamSeqMenor);
    for (int i = 0; i < tamSeqMaior; i++) {
        fscanf(file, " %c", (char *)&seqMaior[i]);
        switch (seqMaior[i]) {
            case 'A': seqMaior[i] = A; break;
            case 'T': seqMaior[i] = T; break;
            case 'G': seqMaior[i] = G; break;
            case 'C': seqMaior[i] = C; break;
        }
    }
    for (int i = 0; i < tamSeqMenor; i++) {
        fscanf(file, " %c", (char *)&seqMenor[i]);
        switch (seqMenor[i]) {
            case 'A': seqMenor[i] = A; break;
            case 'T': seqMenor[i] = T; break;
            case 'G': seqMenor[i] = G; break;
            case 'C': seqMenor[i] = C; break;
        }
    }

    fclose(file);
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

void *preencheMatrizEscores(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    int lin, col, peso;

    if (data->start_lin > tamSeqMenor) {
        pthread_exit(NULL);
    }

    for (lin = data->start_lin; lin <= data->end_lin; lin++) {
        for (col = 1; col <= tamSeqMaior; col++) {
            peso = matrizPesos[seqMenor[lin - 1]][seqMaior[col - 1]];
            diagEscore = matrizEscores[lin - 1][col - 1] + peso;
            linEscore = matrizEscores[lin][col - 1] - penalGap;
            colEscore = matrizEscores[lin - 1][col] - penalGap;

            if (diagEscore >= linEscore && diagEscore >= colEscore) {
                matrizEscores[lin][col] = diagEscore;
            } else if (linEscore > colEscore) {
                matrizEscores[lin][col] = linEscore;
            } else {
                matrizEscores[lin][col] = colEscore;
            }
        }
    }
    pthread_exit(NULL);
}

void geraMatrizEscores(void) {
    int col;
    pthread_t threads[numThreads];
    thread_data_t thread_data[numThreads];

    printf("\nGeracao da Matriz de escores:\n");

    for (col = 0; col <= tamSeqMaior; col++) {
        matrizEscores[0][col] = -1 * (col * penalGap);
    }

    for (int lin = 0; lin <= tamSeqMenor; lin++) {
        matrizEscores[lin][0] = -1 * (lin * penalGap);
    }

    int linhas_por_thread = (tamSeqMenor + numThreads - 1) / numThreads;
    for (int i = 0; i < numThreads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].start_lin = i * linhas_por_thread + 1;
        thread_data[i].end_lin = (i + 1) * linhas_por_thread;
        if (thread_data[i].end_lin > tamSeqMenor) {
            thread_data[i].end_lin = tamSeqMenor;
        }
        pthread_create(&threads[i], NULL, preencheMatrizEscores, (void *)&thread_data[i]);
    }

    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    int linMaior = 1;
    int colMaior = 1;
    int maior = matrizEscores[1][1];
    for (int lin = 1; lin <= tamSeqMenor; lin++) {
        for (int col = 1; col <= tamSeqMaior; col++) {
            if (maior <= matrizEscores[lin][col]) {
                linMaior = lin;
                colMaior = col;
                maior = matrizEscores[lin][col];
            }
        }
    }
    printf("\nMatriz de escores Gerada.");
    printf("\nUltimo Maior escore = %d na celula [%d,%d]", maior, linMaior, colMaior);
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

void *traceBackMultiThread(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    int tbLin = data->start_lin;
    int tbCol = data->start_col;
    int pos = 0;
    int localAlinhaGMaior[maxSeq], localAlinhaGMenor[maxSeq];
    int localTamAlinha = 0;
    
    printf("\nGeracao do Alinhamento Global (Thread %d):\n", data->alignment_index);

    while (tbLin > 0 && tbCol > 0) {
        int peso = matrizPesos[seqMenor[tbLin - 1]][seqMaior[tbCol - 1]];
        int diagEscore = matrizEscores[tbLin - 1][tbCol - 1] + peso;
        int linEscore = matrizEscores[tbLin][tbCol - 1] - penalGap;
        int colEscore = matrizEscores[tbLin - 1][tbCol] - penalGap;

        if (diagEscore >= linEscore && diagEscore >= colEscore) {
            localAlinhaGMaior[pos] = seqMaior[tbCol - 1];
            localAlinhaGMenor[pos] = seqMenor[tbLin - 1];
            tbLin--;
            tbCol--;
        } else if (linEscore >= colEscore) {
            localAlinhaGMaior[pos] = seqMaior[tbCol - 1];
            localAlinhaGMenor[pos] = X;
            tbCol--;
        } else {
            localAlinhaGMaior[pos] = X;
            localAlinhaGMenor[pos] = seqMenor[tbLin - 1];
            tbLin--;
        }
        pos++;
    }

    while (tbLin > 0) {
        localAlinhaGMaior[pos] = X;
        localAlinhaGMenor[pos] = seqMenor[tbLin - 1];
        tbLin--;
        pos++;
    }

    while (tbCol > 0) {
        localAlinhaGMaior[pos] = seqMaior[tbCol - 1];
        localAlinhaGMenor[pos] = X;
        tbCol--;
        pos++;
    }

    localTamAlinha = pos;

    for (int i = 0; i < localTamAlinha / 2; i++) {
        int aux = localAlinhaGMenor[i];
        localAlinhaGMenor[i] = localAlinhaGMenor[localTamAlinha - i - 1];
        localAlinhaGMenor[localTamAlinha - i - 1] = aux;

        aux = localAlinhaGMaior[i];
        localAlinhaGMaior[i] = localAlinhaGMaior[localTamAlinha - i - 1];
        localAlinhaGMaior[localTamAlinha - i - 1] = aux;
    }

    pthread_mutex_lock(&mutex);
    if (data->alignment_index < numAlignments) {
        // Store results only if it's one of the top-k paths
        memcpy(alinhaGMaior[data->alignment_index], localAlinhaGMaior, sizeof(int) * localTamAlinha);
        memcpy(alinhaGMenor[data->alignment_index], localAlinhaGMenor, sizeof(int) * localTamAlinha);
        tamAlinha[data->alignment_index] = localTamAlinha;
    }
    pthread_mutex_unlock(&mutex);

    printf("\nAlinhamento Global (Thread %d) Gerado. Tamanho = %d:\n", data->alignment_index, localTamAlinha);
    for (int i = 0; i < localTamAlinha; i++) {
        printf("%c", mapaBases[localAlinhaGMaior[i]]);
    }
    printf("\n");
    for (int i = 0; i < localTamAlinha; i++) {
        printf("%c", mapaBases[localAlinhaGMenor[i]]);
    }
    printf("\n");

    pthread_exit(NULL);
}

void traceBackParallel() {
    pthread_t threads[numAlignments];
    thread_data_t thread_data[numAlignments];
    int alignment_count = 0;

    for (int i = tamSeqMenor; i > 0 && alignment_count < numAlignments; i--) {
        for (int j = tamSeqMaior; j > 0 && alignment_count < numAlignments; j--) {
            if (multiple_paths_possible(i, j)) {
                thread_data[alignment_count].thread_id = alignment_count;
                thread_data[alignment_count].alignment_index = alignment_count;
                thread_data[alignment_count].start_lin = i;
                thread_data[alignment_count].start_col = j;
                pthread_create(&threads[alignment_count], NULL, traceBackMultiThread, (void *)&thread_data[alignment_count]);
                alignment_count++;
            }
        }
    }

    for (int i = 0; i < alignment_count; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("\nTodos os alinhamentos foram gerados.\n");
}

int multiple_paths_possible(int lin, int col) {
    if (lin == 0 || col == 0) return 0;

    int peso = matrizPesos[seqMenor[lin - 1]][seqMaior[col - 1]];
    int diagEscore = matrizEscores[lin - 1][col - 1] + peso;
    int linEscore = matrizEscores[lin][col - 1] - penalGap;
    int colEscore = matrizEscores[lin - 1][col] - penalGap;

    return (diagEscore == linEscore || diagEscore == colEscore || linEscore == colEscore);
}

void mostraAlinhamentoGlobal(void) {
    for (int k = 0; k < numAlignments; k++) {
        printf("\nAlinhamento %d - Tamanho = %d:\n", k + 1, tamAlinha[k]);

        printf("%c", mapaBases[alinhaGMaior[k][0]]);
        for (int i = 1; i < tamAlinha[k]; i++) {
            printf("%c", mapaBases[alinhaGMaior[k][i]]);
        }
        printf("\n");

        printf("%c", mapaBases[alinhaGMenor[k][0]]);
        for (int i = 1; i < tamAlinha[k]; i++) {
            printf("%c", mapaBases[alinhaGMenor[k][i]]);
        }
        printf("\n");
    }
}

void leNumAlignments(void) {
    printf("\nDigite o numero de alinhamentos a serem exibidos: ");
    do {
        scanf("%d", &numAlignments);
        if (numAlignments < 1) {
            printf("Numero invalido. Digite novamente: ");
        }
    } while (numAlignments < 1);
}

void leNumThreads(void) {
    printf("\nDigite o numero de threads a serem usados: ");
    do {
        scanf("%d", &numThreads);
        if (numThreads < 1) {
            printf("Numero invalido. Digite novamente: ");
        }
    } while (numThreads < 1);
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
        printf("\n<12> Definir Numero de Threads");
        printf("\n<13> Sair");
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
            break;
        case 8:
            mostraMatrizEscores();
            break;
        case 9:
            leNumAlignments();
            traceBackParallel();
            break;
        case 10:
            mostraAlinhamentoGlobal();
            break;
        case 11:
            salvaMatrizEscores("matriz_escores.txt");
            break;
        case 12:
            leNumThreads();
            break;
    }
}

int main(void) {
    int opcao;

    srand(time(NULL));

    do {
        printf("\n\nPrograma Needleman-Wunsch Sequencial\n");
        opcao = menuOpcao();
        trataOpcao(opcao);

    } while (opcao != sair);

    pthread_mutex_destroy(&mutex);

    return 0;
}
