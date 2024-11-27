#ifndef CROMOSOMA_H
#define CROMOSOMA_H

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

class Cromosoma {
   public:
    vector<double> genes;            // Representa el cromosoma
    vector<double> cultivoPlantado;  // Arreglo de cultivoPlantado para el cromosoma
    double valorObjetivo;            // Almacena el valor de la funcion objetivo

    Cromosoma() : genes(0), cultivoPlantado(0), valorObjetivo(0.0) {}

    Cromosoma(int dimension)
        : genes(dimension, 0.0), cultivoPlantado(dimension, 0.0), valorObjetivo(0.0) {}

    // Validar si el cultivo puede crecer en los meses actuales y subsiguientes
    static bool esCultivable(const vector<int>& cultivable, int cultivo, int mes, int periodoCrecimiento, int numeroCultivos) {
        for (int m = 0; m < periodoCrecimiento && (mes + m) < cultivable.size() / numeroCultivos; ++m) {
            if (cultivable[(cultivo) + numeroCultivos * (mes + m)] == 0) {
                return false;  // El cultivo no puede ser cultivado en este mes
            }
        }
        return true;  // El cultivo puede ser cultivado en todos los meses requeridos
    }

    // Calcular si se debe entrar al bucle de inicializacion basado en la condicion exponencial
    static bool debeEntrarAlBucleDeInicializacion(double areaDisponible) {
        double resultado = -0.7 * exp(-6 * areaDisponible + 5.25) + 107;
        return resultado > (rand() % 100);
    }

    // Validar si el agua disponible es suficiente para el crecimiento del cultivo durante el periodo de crecimiento con una probabilidad de continuar basada en la escasez
    static bool esAguaSuficiente(const vector<double>& aguaDisponible, const vector<double>& requerimientoAgua, int cultivo, int mes, int periodoCrecimiento, double areaUsada, double areaTotalDisponible) {
        double areaEnHectareas = areaUsada * areaTotalDisponible;  // Convertir porcentaje de area usada a hectareas
        random_device rd;
        mt19937 gen(rd());

        for (int m = 0; m < periodoCrecimiento && (mes + m) < aguaDisponible.size(); ++m) {
            double aguaRequerida = requerimientoAgua[cultivo] * areaEnHectareas;  // Agua requerida para este cultivo en el area en hectareas
            double disponible = aguaDisponible[mes + m];

            if (disponible < aguaRequerida) {
                double porcentajeEscasez = (aguaRequerida - disponible) / aguaRequerida;  // Calcular cuanto el agua requerida excede el agua disponible
                double probabilidadContinuar = 1.0 - porcentajeEscasez;                   // Probabilidad de continuar el bucle

                // Verificar si debemos continuar basado en la probabilidad calculada
                uniform_real_distribution<> dis(0.0, 1.0);
                if (dis(gen) > probabilidadContinuar) {
                    return false;  // No hay suficiente agua y la probabilidad de continuar fallo
                }
            }
        }
        return true;  // Suficiente agua para todo el periodo de crecimiento o continuado basado en la probabilidad
    }

    // Metodo estatico para inicializar una luciernaga y devolverla
    static Cromosoma inicializar(int dimension, int numeroCultivos, int meses, const vector<int>& mesesCultivo,
                                 const vector<double>& requerimientoAgua, const vector<int>& cultivable,
                                 const vector<double>& aguaInicialDisponible, double areaTotalDisponible) {
        Cromosoma nuevoCromosoma(dimension);                    // Crear un nuevo objeto Cromosoma
        vector<double> areaDisponible(meses, 1.0);              // Inicializar area disponible al 100% para cada mes
        vector<double> aguaDisponible = aguaInicialDisponible;  // Copiar disponibilidad inicial de agua

        random_device rd;
        mt19937 gen(rd());
        chi_squared_distribution<> dist(5);

        // Inicializar los arreglos genes y cultivoPlantado
        for (int mes = 0; mes < meses; ++mes) {
            while (debeEntrarAlBucleDeInicializacion(areaDisponible[mes])) {
                // Seleccionar un cultivo aleatorio
                int cultivo = rand() % numeroCultivos;
                int periodoCrecimiento = mesesCultivo[cultivo];

                // Validar si el cultivo puede ser cultivado
                if (!esCultivable(cultivable, cultivo, mes, periodoCrecimiento, numeroCultivos)) {
                    continue;
                }

                // Determinar el porcentaje de area a usar
                double prcAreaUsada = 8 * dist(gen) / 100.0;
                double areaUsada = (prcAreaUsada > 1 ? 0.0 : prcAreaUsada) * areaDisponible[mes];

                // Verificar suficiencia de agua
                if (!esAguaSuficiente(aguaDisponible, requerimientoAgua, cultivo, mes, periodoCrecimiento, areaUsada, areaTotalDisponible)) {
                    continue;
                }

                // Asignar area a los arreglos genes y cultivoPlantado para los meses actuales y subsiguientes
                for (int m = 0; m < periodoCrecimiento && (mes + m) < meses; ++m) {
                    int indice = cultivo + numeroCultivos * (mes + m);
                    nuevoCromosoma.genes[indice] += areaUsada;
                    areaDisponible[mes + m] -= areaUsada;

                    double areaEnHectareas = areaUsada * areaTotalDisponible;
                    double aguaADeducir = requerimientoAgua[cultivo] * areaEnHectareas;

                    if (aguaDisponible[mes + m] < aguaADeducir) {
                        aguaADeducir = aguaDisponible[mes + m];
                    }
                    aguaDisponible[mes + m] -= aguaADeducir;
                }

                // Actualizar el arreglo cultivoPlantado
                int indicePlantacion = cultivo + numeroCultivos * mes;
                nuevoCromosoma.cultivoPlantado[indicePlantacion] = areaUsada;
            }

            // Transferir agua no usada al siguiente mes
            if (mes < meses - 1) {
                aguaDisponible[mes + 1] += aguaDisponible[mes];  // Agregar agua restante al siguiente mes
            }
        }

        return nuevoCromosoma;
    }

    void imprimirCromosoma(int numeroCultivos) const {
        cout << fixed << setprecision(2);
        cout << "(";
        for (size_t j = 0; j < genes.size(); ++j) {
            cout << genes[j];
            if (j != genes.size() - 1) {
                cout << ", ";
            }
            if ((j + 1) % numeroCultivos == 0 && j != genes.size() - 1) {
                cout << "| ";
            }
        }
        cout << ")" << endl;
    }

    void imprimirDetallesCromosoma(int numeroCultivos, int meses, double areaTotalDisponible,
                                    const vector<double>& requerimientoAgua, const vector<int>& mesesCultivo,
                                    const vector<double>& maxCosechaPorArea) const {
        cout << "Informe detallado de la luciernaga:" << endl;

        double totalAguaUsadaMes = 0.0;
        double cosechaTotal = 0.0;
        vector<double> cosechaPorCultivo(numeroCultivos, 0.0);

        for (int mes = 0; mes < meses; ++mes) {
            cout << "Mes " << mes + 1 << ":" << endl;
            totalAguaUsadaMes = 0.0;

            for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
                int indice = cultivo + numeroCultivos * mes;
                double areaAsignadaPorcentaje = genes[indice];
                double areaRealAsignada = areaAsignadaPorcentaje * areaTotalDisponible;


                if (areaAsignadaPorcentaje <= 0) continue;


                double aguaUsada = requerimientoAgua[cultivo] * areaRealAsignada;
                totalAguaUsadaMes += aguaUsada;


                double cosechaEsperada = (maxCosechaPorArea[cultivo] * areaRealAsignada) / mesesCultivo[cultivo];
                cosechaTotal += cosechaEsperada;
                cosechaPorCultivo[cultivo] += cosechaEsperada;

                cout << "  Cultivo " << cultivo + 1 << ":" << endl;
                cout << "    Area Asignada (porcentaje): " << areaAsignadaPorcentaje * 100 << "%" << endl;
                cout << "    Area Real Asignada (hectareas): " << areaRealAsignada << " hectareas" << endl;
                cout << "    Agua Usada (metros cubicos): " << aguaUsada << " metros cubicos" << endl;
            }

            cout << "Agua total usada en el mes " << mes + 1 << ": " << totalAguaUsadaMes << " metros cubicos" << endl;
            cout << endl;
        }

        // Imprimir la cosecha total por cultivo
        cout << "Cosecha total por cultivo:" << endl;
        for (int cultivo = 0; cultivo < numeroCultivos; ++cultivo) {
            cout << "  Cultivo " << cultivo + 1 << ": " << cosechaPorCultivo[cultivo] << " toneladas" << endl;
        }

        //Imprimir el rendimiento
        cout << "Cosecha total de todos los cultivos: " << cosechaTotal << " toneladas" << endl;
    }
};

#endif /* CROMOSOMA_H */

