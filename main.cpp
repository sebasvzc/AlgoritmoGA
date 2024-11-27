#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

using namespace std;

#include "Generacion.h"

int main() {
    srand(static_cast<unsigned int>(time(0)));

    // Parámetros del Algoritmo Genético
    int tamanoPoblacion = 100;  // Tamaño de la población
    int maximoGeneraciones = 100;  // Número máximo de generaciones

    // Variables específicas del problema
    int numeroCultivos = 5;                   // Número de cultivos
    int meses = 8;                            // Número de meses
    int dimension = numeroCultivos * meses;   // Dimensión total del cromosoma

    Cultivacion cultivacion(meses, numeroCultivos);

    // Inicializar la población y la estructura cultivoPlantado
    Generacion poblacion(tamanoPoblacion, dimension);

    poblacion.inicializarCromosomas(numeroCultivos, meses, cultivacion);
    poblacion.inicializarValoresObjetivo(numeroCultivos, meses, cultivacion);

    Cromosoma mejorCromosoma = poblacion.encontrarMejorCromosoma();
    double mejorAptitud = mejorCromosoma.valorObjetivo;

    // Bucle externo: iterar a través de las generaciones
    for (int generacion = 0; generacion < maximoGeneraciones; ++generacion) {
        poblacion.obtenerNuevaGeneracion(tamanoPoblacion, numeroCultivos, meses, cultivacion);
        poblacion.inicializarValoresObjetivo(numeroCultivos, meses, cultivacion);

        if (poblacion.encontrarMejorCromosoma().valorObjetivo < mejorAptitud) {
            mejorCromosoma = poblacion.encontrarMejorCromosoma();
            mejorAptitud = mejorCromosoma.valorObjetivo;
        }
    }

    mejorCromosoma.imprimirDetallesCromosoma(numeroCultivos, meses,
                                             cultivacion.areaTotalDisponible,
                                             cultivacion.requerimientoAgua,
                                             cultivacion.mesesCultivo,
                                             cultivacion.maxCosechaPorArea);
    return 0;
}
