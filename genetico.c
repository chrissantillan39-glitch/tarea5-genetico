#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

// ---------- Parametros ----------
#define N_IND          50     // tamano de la poblacion
#define P_CRUZA        0.9    // probabilidad de cruza
#define P_MUTACION     0.1    // probabilidad de que un hijo mute (se voltea 1 bit)
#define GENERACIONES   100    // cuantas veces se repite el proceso
#define T              16     // bits por variable (t en la formula)
#define LONG_GEN       (2 * T) // genes totales: T bits para x + T bits para y
#define A_LIM          (-5.12) // limite inferior del dominio
#define B_LIM          (5.12)  // limite superior del dominio

// ---------- Individuo ----------
typedef struct {
    int genes[LONG_GEN]; // cromosoma binario: [0..T-1] = x, [T..2T-1] = y
    double x, y;         // valores decodificados
    double f;            // f(x,y) evaluada
} Individuo;

// ---------- Utilidades ----------
double aleatorio01(void) {
    return (double)rand() / RAND_MAX;
}

// X = A + (valor_decimal / (2^T - 1)) * (B - A)
double decodificar(const int *bits) {
    unsigned long decimal = 0;
    for (int i = 0; i < T; i++) {
        decimal = (decimal << 1) | bits[i];
    }
    unsigned long maximo = (1UL << T) - 1;
    return A_LIM + ((double)decimal / maximo) * (B_LIM - A_LIM);
}

// Rastrigin de dos variables
double rastrigin(double x, double y) {
    return 20.0 + (x * x - 10.0 * cos(2.0 * PI * x))
                + (y * y - 10.0 * cos(2.0 * PI * y));
}

// Decodifica los genes y evalua la funcion
void evaluar(Individuo *ind) {
    ind->x = decodificar(&ind->genes[0]);
    ind->y = decodificar(&ind->genes[T]);
    ind->f = rastrigin(ind->x, ind->y);
}

// Comparador para qsort: menor f(x,y) primero (minimizamos)
int comparar(const void *a, const void *b) {
    double fa = ((const Individuo *)a)->f;
    double fb = ((const Individuo *)b)->f;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

// ---------- Operadores geneticos ----------
void poblacion_inicial(Individuo *pob) {
    for (int i = 0; i < N_IND; i++) {
        for (int j = 0; j < LONG_GEN; j++) {
            pob[i].genes[j] = rand() % 2;
        }
        evaluar(&pob[i]);
    }
}

// Cruza de un punto: con probabilidad P_CRUZA intercambia las colas
void cruzar(const Individuo *p1, const Individuo *p2, Individuo *h1, Individuo *h2) {
    *h1 = *p1;
    *h2 = *p2;
    if (aleatorio01() < P_CRUZA) {
        int punto = 1 + rand() % (LONG_GEN - 1); // punto entre 1 y LONG_GEN-1
        for (int j = punto; j < LONG_GEN; j++) {
            h1->genes[j] = p2->genes[j];
            h2->genes[j] = p1->genes[j];
        }
    }
}

// Mutacion: con probabilidad P_MUTACION voltea un bit al azar
void mutar(Individuo *ind) {
    if (aleatorio01() < P_MUTACION) {
        int pos = rand() % LONG_GEN;
        ind->genes[pos] = 1 - ind->genes[pos];
    }
}

int main(void) {
    srand((unsigned)time(NULL));

    // Arreglo de 2*N: primeros N = poblacion, siguientes N = hijos
    Individuo *pob = (Individuo *)malloc(2 * N_IND * sizeof(Individuo));
    if (pob == NULL) {
        fprintf(stderr, "Error: sin memoria\n");
        return 1;
    }

    poblacion_inicial(pob);
    qsort(pob, N_IND, sizeof(Individuo), comparar);

    for (int gen = 1; gen <= GENERACIONES; gen++) {
        // 3 y 4. Cruza (padres aleatorios) y mutacion de los hijos
        for (int k = 0; k < N_IND; k += 2) {
            const Individuo *p1 = &pob[rand() % N_IND];
            const Individuo *p2 = &pob[rand() % N_IND];

            Individuo *h1 = &pob[N_IND + k];
            Individuo *h2 = &pob[N_IND + k + 1];

            cruzar(p1, p2, h1, h2);
            mutar(h1);
            mutar(h2);
        }

        // 5. Evaluar a los hijos (los padres ya estan evaluados)
        for (int i = N_IND; i < 2 * N_IND; i++) {
            evaluar(&pob[i]);
        }

        // 6. Ordenar toda la poblacion (padres + hijos) y quedarnos con los N mejores
        qsort(pob, 2 * N_IND, sizeof(Individuo), comparar);
        // Los primeros N_IND son los mejores; el resto se descarta en la siguiente vuelta

        if (gen % 10 == 0 || gen == 1) {
            printf("Gen %3d | mejor f = %.6f en (x=%.4f, y=%.4f)\n",
                   gen, pob[0].f, pob[0].x, pob[0].y);
        }
    }

    printf("\nMejor individuo final:\n");
    printf("x = %.6f, y = %.6f, f(x,y) = %.6f\n", pob[0].x, pob[0].y, pob[0].f);
    printf("(El minimo global de Rastrigin es 0 en (0,0))\n");

    free(pob);
    return 0;
}
