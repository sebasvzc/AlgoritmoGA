#ifndef GENERACION_H
#define GENERACION_H

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

#include "Cromosoma.h"
#include "Cultivacion.h"

class Generacion {
   public:
    int tamanoPoblacion = 100;       // Tamano de la poblacion
    double tasaMutacion = 0.05;      // Parametros del algoritmo genetico
    double tasaCruce = 0.8;          // Parametros del algoritmo genetico
    vector<Cromosoma> poblacion;     // Vector de cromosomas
    vector<double> valoresObjetivo;  // Valores de la funcion objetivo

    Generacion() : tamanoPoblacion(100) {
        poblacion.reserve(tamanoPoblacion);
        valoresObjetivo.reserve(tamanoPoblacion);
    }

    Generacion(int tamanoPoblacion, int dimension) : tamanoPoblacion(tamanoPoblacion) {
        for (int i = 0; i < tamanoPoblacion; ++i) {
            poblacion.emplace_back(dimension);
        }
        valoresObjetivo.resize(tamanoPoblacion, 0.0);
    }

    void inicializarCromosomas(int numeroCultivos, int meses, Cultivacion& cultivacion) {
        int dimension = numeroCultivos * meses;

        for (int k = 0; k < tamanoPoblacion; ++k) {
            Cromosoma nuevoCromosoma = Cromosoma::inicializar(dimension, numeroCultivos, meses,
                                                              cultivacion.mesesCultivo,
                                                              cultivacion.requerimientoAgua,
                                                              cultivacion.cultivable,
                                                              cultivacion.aguaInicialDisponible,
                                                              cultivacion.areaTotalDisponible);
            poblacion.push_back(nuevoCromosoma);
        }
    }

    void inicializarHijoValidado(Cromosoma& hijoValidado, const Cromosoma& hijo) {
        hijoValidado.genes.assign(hijo.genes.size(), 0.0);
        hijoValidado.cultivoPlantado.assign(hijo.cultivoPlantado.size(), 0.0);
    }

    vector<int> generarSecuenciaAleatoriaCultivos(int numeroCultivos, mt19937& gen) {
        vector<int> secuenciaCultivos(numeroCultivos);
        iota(secuenciaCultivos.begin(), secuenciaCultivos.end(), 0);  // Llenar con 0, 1, 2, ..., numeroCultivos - 1
        shuffle(secuenciaCultivos.begin(), secuenciaCultivos.end(), gen);
        return secuenciaCultivos;
    }

    double ajustarAreaAsignada(double areaAsignada, double areaDisponible, mt19937& gen) {
        if (areaAsignada > areaDisponible) {
            chi_squared_distribution<> dist(5);
            double prcAreaUsada;
            do {
                prcAreaUsada = 8 * dist(gen) / 100.0;
            } while (prcAreaUsada > 1);
            areaAsignada = prcAreaUsada * areaDisponible;
        }
        return areaAsignada;
    }

    void actualizarHijoValidado(Cromosoma& hijoValidado, vector<double>& areaDisponible, double areaAsignada,
                                int cultivo, int mes, int numeroCultivos, int meses, const Cultivacion& cultivacion) {
        for (int m = 0; m < cultivacion.mesesCultivo[cultivo] && (mes + m) < meses; ++m) {
            int indiceSubsecuente = cultivo + numeroCultivos * (mes + m);

            // Actualizar cromosoma para los meses actuales y subsecuentes
            hijoValidado.genes[indiceSubsecuente] += areaAsignada;  // Sumar area asignada a los meses subsecuentes
            if (m == 0) {
                hijoValidado.cultivoPlantado[indiceSubsecuente] = areaAsignada;  // Establecer el area porcentual
            }

            // Restar recursos para los meses actuales y subsecuentes
            areaDisponible[mes + m] -= areaAsignada;
        }
    }

    void validarHijo(const Cromosoma& hijo, int numeroCultivos, int meses, Cultivacion& cultivacion, Cromosoma& hijoValidado) {
        // Inicializar hijo validado
        inicializarHijoValidado(hijoValidado, hijo);

        random_device rd;
        mt19937 gen(rd());

        vector<double> areaDisponible(meses, 1.0);

        for (int mes = 0; mes < meses; ++mes) {
            vector<int> secuenciaCultivos = generarSecuenciaAleatoriaCultivos(numeroCultivos, gen);

            for (int cultivo : secuenciaCultivos) {
                int indice = cultivo + numeroCultivos * mes;
                if (hijo.cultivoPlantado[indice] > 0.0) {
                    double areaAsignada = ajustarAreaAsignada(hijo.cultivoPlantado[indice], areaDisponible[mes], gen);

                    actualizarHijoValidado(hijoValidado, areaDisponible, areaAsignada, cultivo, mes, numeroCultivos, meses, cultivacion);
                }
            }
        }
    }

    int seleccionarGenNoCeroAleatorio(const Cromosoma& cromosoma, int numeroCultivos, mt19937& gen) {
        vector<int> indicesNoCero;
        for (int i = 0; i < cromosoma.cultivoPlantado.size(); ++i) {
            if (cromosoma.cultivoPlantado[i] > 0.0) {
                indicesNoCero.push_back(i);
            }
        }
        if (indicesNoCero.empty()) return -1;  // No hay genes no cero para mutar

        uniform_int_distribution<> indiceAleatorio(0, indicesNoCero.size() - 1);
        return indicesNoCero[indiceAleatorio(gen)];
    }

    double reducirAreaGen(Cromosoma& cromosoma, int indiceSeleccionado, mt19937& gen) {
        uniform_real_distribution<> reduccionAleatoria(0.01, 0.1);  // Porcentaje de reduccion aleatoria
        double porcentajeReduccion = reduccionAleatoria(gen);
        double areaAReducir = cromosoma.cultivoPlantado[indiceSeleccionado] * porcentajeReduccion;
        cromosoma.cultivoPlantado[indiceSeleccionado] -= areaAReducir;
        return areaAReducir;
    }

    void actualizarGenTrasMutacion(Cromosoma& cromosoma, int indiceSeleccionado, double areaAReducir,
                                   const Cultivacion& cultivacion, int numeroCultivos, int meses) {
        int mesSeleccionado = indiceSeleccionado / numeroCultivos;
        int cultivoSeleccionado = indiceSeleccionado % numeroCultivos;

        for (int m = 0; m < cultivacion.mesesCultivo[cultivoSeleccionado] && (mesSeleccionado + m) < meses; ++m) {
            int indice = cultivoSeleccionado + numeroCultivos * (mesSeleccionado + m);
            cromosoma.genes[indice] -= areaAReducir;
        }
    }

    void reinicializarCromosoma(Cromosoma& cromosoma, int numeroCultivos, int meses, Cultivacion& cultivacion, mt19937& gen) {
        vector<double> areaDisponible(meses, 1.0);                          // Reiniciar area disponible
        vector<double> aguaDisponible = cultivacion.aguaInicialDisponible;  // Reiniciar disponibilidad de agua

        for (int mes = 0; mes < meses; ++mes) {
            // Recalcular area disponible para el mes actual y todos los meses siguientes
            for (int m = mes; m < meses; ++m) {
                areaDisponible[m] = 1.0;  // Reiniciar area disponible al maximo (100%)
                for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
                    int indice = cultivo + numeroCultivos * m;
                    if (cromosoma.genes[indice] > 0.0)
                        areaDisponible[m] -= cromosoma.genes[indice];
                }
            }
            // Agregar nueva area aleatoria si es necesario
            if (Cromosoma::debeEntrarAlBucleDeInicializacion(areaDisponible[mes])) {
                // Elegir un cultivo aleatorio
                int cultivo = rand() % numeroCultivos;
                int periodoCrecimiento = cultivacion.mesesCultivo[cultivo];

                // Validar si es cultivable y hay suficiente agua
                if (!Cromosoma::esCultivable(cultivacion.cultivable, cultivo, mes, periodoCrecimiento, numeroCultivos))
                    continue;

                // Determinar el area minima disponible durante el periodo de crecimiento
                double areaMinimaDisponible = areaDisponible[mes];
                for (int m = 1; m < periodoCrecimiento && (mes + m) < meses; ++m)
                    areaMinimaDisponible = min(areaMinimaDisponible, areaDisponible[mes + m]);

                // Determinar el area a usar
                chi_squared_distribution<> dist(5);
                double prcAreaUsada = 8 * dist(gen) / 100.0;
                double areaUsada = (prcAreaUsada > 1 ? 0.0 : prcAreaUsada) * areaMinimaDisponible;

                // Validar suficiencia de agua
                if (!Cromosoma::esAguaSuficiente(aguaDisponible, cultivacion.requerimientoAgua, cultivo, mes, periodoCrecimiento, areaUsada, cultivacion.areaTotalDisponible))
                    continue;

                // Asignar area a genes y cultivoPlantado
                for (int m = 0; m < periodoCrecimiento && (mes + m) < meses; ++m) {
                    int indice = cultivo + numeroCultivos * (mes + m);
                    cromosoma.genes[indice] += areaUsada;
                    areaDisponible[mes + m] -= areaUsada;

                    double areaEnHectareas = areaUsada * cultivacion.areaTotalDisponible;
                    double aguaADeducir = cultivacion.requerimientoAgua[cultivo] * areaEnHectareas;
                    if (aguaDisponible[mes + m] < aguaADeducir)
                        aguaADeducir = aguaDisponible[mes + m];
                    aguaDisponible[mes + m] -= aguaADeducir;
                }

                int indicePlantacion = cultivo + numeroCultivos * mes;
                cromosoma.cultivoPlantado[indicePlantacion] += areaUsada;
            }
        }
    }

    void mutarCromosoma(Cromosoma& cromosoma, int numeroCultivos, int meses, Cultivacion& cultivacion, mt19937& gen) {
        // Seleccionar un gen aleatorio para mutar
        int indiceSeleccionado = seleccionarGenNoCeroAleatorio(cromosoma, numeroCultivos, gen);
        if (indiceSeleccionado == -1) return;  // No es posible mutar

        // Reducir el area del gen seleccionado
        double areaAReducir = reducirAreaGen(cromosoma, indiceSeleccionado, gen);

        // Actualizar el cromosoma y cultivoPlantado
        actualizarGenTrasMutacion(cromosoma, indiceSeleccionado, areaAReducir, cultivacion, numeroCultivos, meses);

        // Reinicializar el cromosoma mutado
        reinicializarCromosoma(cromosoma, numeroCultivos, meses, cultivacion, gen);
    }

    pair<Cromosoma, Cromosoma> seleccionarPadres() {
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dist(0, tamanoPoblacion - 1);

        Cromosoma padre1 = poblacion[dist(gen)];
        Cromosoma padre2 = poblacion[dist(gen)];

        return std::make_pair(padre1, padre2);
    }

    pair<Cromosoma, Cromosoma> realizarCruce(const Cromosoma& padre1, const Cromosoma& padre2, int numeroCultivos, int meses) {
        random_device rd;
        mt19937 gen(rd());

        if (uniform_real_distribution<>(0.0, 1.0)(gen) >= tasaCruce) {
            return std::make_pair(padre1, padre2);  // No se realiza cruce, devolver los padres como hijos
        }

        uniform_int_distribution<> distMes(0, meses - 1);
        int puntoCruce = distMes(gen);

        Cromosoma hijo1 = padre1;
        Cromosoma hijo2 = padre2;

        for (int mes = puntoCruce; mes < meses; ++mes) {
            for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
                int indice = cultivo + numeroCultivos * mes;
                swap(hijo1.genes[indice], hijo2.genes[indice]);
                swap(hijo1.cultivoPlantado[indice], hijo2.cultivoPlantado[indice]);
            }
        }

        return std::make_pair(hijo1, hijo2);
    }

    Cromosoma validarYMutar(Cromosoma hijo, int numeroCultivos, int meses, Cultivacion& cultivacion) {
        Cromosoma hijoValidado;
        validarHijo(hijo, numeroCultivos, meses, cultivacion, hijoValidado);

        random_device rd;
        mt19937 gen(rd());
        if (uniform_real_distribution<>(0, 1)(gen) < tasaMutacion) {
            mutarCromosoma(hijoValidado, numeroCultivos, meses, cultivacion, gen);
        }

        return hijoValidado;
    }

    void agregarAPoblacion(const Cromosoma& cromosoma) {
        poblacion.push_back(cromosoma);
    }

    void combinarGeneraciones(Generacion& siguienteGeneracion, int numeroCultivos, int meses, Cultivacion& cultivacion) {
        Generacion generacionCombinada;
        generacionCombinada.poblacion.insert(generacionCombinada.poblacion.end(), poblacion.begin(), poblacion.end());
        generacionCombinada.poblacion.insert(generacionCombinada.poblacion.end(), siguienteGeneracion.poblacion.begin(), siguienteGeneracion.poblacion.end());

        generacionCombinada.inicializarValoresObjetivo(numeroCultivos, meses, cultivacion);

        sort(generacionCombinada.poblacion.begin(), generacionCombinada.poblacion.end(), [](const Cromosoma& a, const Cromosoma& b) {
            return a.valorObjetivo < b.valorObjetivo;
        });

        // Seleccionar los mejores individuos para la siguiente generación
        poblacion.clear();
        poblacion.insert(poblacion.end(), generacionCombinada.poblacion.begin(), generacionCombinada.poblacion.begin() + tamanoPoblacion);
    }

    void obtenerNuevaGeneracion(int tamanoPoblacion, int numeroCultivos, int meses, Cultivacion& cultivacion) {
        Generacion siguienteGeneracion;

        for (int i = 0; i < tamanoPoblacion / 2; ++i) {
            // Seleccionar padres
            pair<Cromosoma, Cromosoma> padres = seleccionarPadres();
            Cromosoma padre1 = padres.first;
            Cromosoma padre2 = padres.second;

            // Realizar cruce
            pair<Cromosoma, Cromosoma> hijos = realizarCruce(padre1, padre2, numeroCultivos, meses);
            Cromosoma hijo1 = hijos.first;
            Cromosoma hijo2 = hijos.second;

            // Validar y mutar hijos
            hijo1 = validarYMutar(hijo1, numeroCultivos, meses, cultivacion);
            hijo2 = validarYMutar(hijo2, numeroCultivos, meses, cultivacion);

            // Agregar hijos a la siguiente generación
            siguienteGeneracion.agregarAPoblacion(hijo1);
            siguienteGeneracion.agregarAPoblacion(hijo2);
        }

        // Combinar generaciones actual y siguiente
        combinarGeneraciones(siguienteGeneracion, numeroCultivos, meses, cultivacion);
    }

    // Calcular el agua total requerida para todos los cultivos en un mes
    double calcularAguaTotalRequerida(const Cromosoma& cromosoma, int numeroCultivos, int mes, double areaTotalDisponible,
                                      const vector<double>& requerimientoAgua) const {
        double aguaTotalRequerida = 0.0;
        for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
            int indice = cultivo + numeroCultivos * mes;
            double areaAsignada = cromosoma.genes[indice];
            if (areaAsignada > 0) {
                aguaTotalRequerida += requerimientoAgua[cultivo] * areaAsignada * areaTotalDisponible;
            }
        }
        return aguaTotalRequerida;
    }

    // Calcular el coeficiente de agua para un mes
    double calcularCoeficienteAgua(double aguaTotalRequerida, double aguaTotalDisponible) const {
        return aguaTotalRequerida > 0 ? min(1.0, max(0.0, aguaTotalDisponible / aguaTotalRequerida)) : 1.0;
    }

    // Calcular cosecha esperada y real para cada cultivo en un mes
    double calcularCosechaCultivo(const Cromosoma& cromosoma, int numeroCultivos, int mes, double areaTotalDisponible,
                                  double coeficienteAgua, double conductividadElectrica,
                                  const vector<int>& mesesCultivo, const vector<double>& maxCosechaPorArea,
                                  const vector<double>& susceptibilidadAgua, const vector<double>& reduccionRendimiento,
                                  const vector<double>& salinidadCritica) const {
        double cosechaMensual = 0.0;
        for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
            int indice = cultivo + numeroCultivos * mes;
            double areaAsignada = cromosoma.genes[indice];

            if (areaAsignada <= 0) continue;

            // Calcular cosecha esperada
            double cosechaEsperada = (maxCosechaPorArea[cultivo] * areaAsignada) / mesesCultivo[cultivo];

            // Calcular efecto del agua
            double factorExponente = (coeficienteAgua * susceptibilidadAgua[cultivo]) / areaAsignada;
            double efectoAgua = 1 - exp(-factorExponente);

            // Calcular efecto de la salinidad
            double impactoSalinidad = reduccionRendimiento[cultivo] * (conductividadElectrica - salinidadCritica[cultivo]);
            double efectoSalinidad = min(1.0, max(0.0, 1.0 - impactoSalinidad / 100.0));

            // Calcular cosecha real
            double cosechaReal = cosechaEsperada * efectoAgua * efectoSalinidad;
            cosechaMensual += cosechaReal;
        }
        return cosechaMensual;
    }

    // Actualizar la salinidad para el siguiente mes
    double actualizarSalinidad(const Cromosoma& cromosoma, int numeroCultivos, int mes, double areaTotalDisponible,
                                const vector<double>& cambioSalinidadPorArea) const {
        double cambioTotalSalinidad = 0.0;
        for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
            int indice = cultivo + numeroCultivos * mes;
            double areaAsignada = cromosoma.genes[indice];
            cambioTotalSalinidad += cambioSalinidadPorArea[cultivo] * (areaAsignada * areaTotalDisponible);
        }
        return cambioTotalSalinidad;
    }

    // Transferir agua no utilizada al siguiente mes
    void transferirAguaSobrante(vector<double>& aguaDisponible, int mes, double aguaTotalRequerida) const {
        if (mes < aguaDisponible.size() - 1) {
            aguaDisponible[mes + 1] += max(0.0, aguaDisponible[mes] - aguaTotalRequerida);
        }
    }

    double funcionObjetivo(const Cromosoma& cromosoma, int numeroCultivos, int meses, Cultivacion& cultivacion) {
        double cosechaTotal = 0.0;
        double conductividadElectrica = cultivacion.conductividadElectrica;
        vector<double> aguaDisponible = cultivacion.aguaInicialDisponible;

        for (int mes = 0; mes < meses; ++mes) {
            // Calcular el agua total requerida
            double aguaTotalRequerida = calcularAguaTotalRequerida(cromosoma, numeroCultivos, mes, cultivacion.areaTotalDisponible, cultivacion.requerimientoAgua);

            // Calcular el coeficiente de agua
            double coeficienteAgua = calcularCoeficienteAgua(aguaTotalRequerida, aguaDisponible[mes]);

            // Calcular la cosecha del cultivo para el mes
            double cosechaMensual = calcularCosechaCultivo(cromosoma, numeroCultivos, mes, cultivacion.areaTotalDisponible,
                                                           coeficienteAgua, conductividadElectrica, cultivacion.mesesCultivo,
                                                           cultivacion.maxCosechaPorArea, cultivacion.susceptibilidadAgua, cultivacion.reduccionRendimiento, cultivacion.salinidadCritica);

            // Actualizar la salinidad para el siguiente mes
            if (mes < meses - 1) {
                double cambioSalinidad = actualizarSalinidad(cromosoma, numeroCultivos, mes, cultivacion.areaTotalDisponible, cultivacion.cambioSalinidadPorArea);
                conductividadElectrica += cambioSalinidad;
            }

            // Transferir agua no utilizada
            transferirAguaSobrante(aguaDisponible, mes, aguaTotalRequerida);

            // Acumular la cosecha mensual
            cosechaTotal += cosechaMensual;
        }

        // Devolver el negativo de la cosecha total porque buscamos minimizar la función objetivo
        return -cosechaTotal;
    }

    void actualizarValorObjetivo(size_t indice, int numeroCultivos, int meses, Cultivacion& cultivacion) {
        double valor = funcionObjetivo(poblacion[indice], numeroCultivos, meses, cultivacion);
        poblacion[indice].valorObjetivo = valor;  // Actualizar el valor objetivo del cromosoma
        valoresObjetivo[indice] = valor;          // Actualizar el vector de valores objetivo
    }

    void inicializarValoresObjetivo(int numeroCultivos, int meses, Cultivacion& cultivacion) {
        valoresObjetivo.resize(poblacion.size());
        for (size_t i = 0; i < poblacion.size(); ++i) {
            actualizarValorObjetivo(i, numeroCultivos, meses, cultivacion);
        }
    }

    void imprimirPoblacion(int numeroCultivos, double areaTotalDisponible) const {
        for (const Cromosoma& cromosoma : poblacion) {
            cromosoma.imprimirCromosoma(numeroCultivos);
        }
    }

    Cromosoma encontrarMejorCromosoma() const {
        auto mejorCromosoma = poblacion[0];
        for (const auto& cromosoma : poblacion) {
            if (cromosoma.valorObjetivo < mejorCromosoma.valorObjetivo) {
                mejorCromosoma = cromosoma;
            }
        }
        return mejorCromosoma;
    }
};

#endif /* GENERACION_H */
