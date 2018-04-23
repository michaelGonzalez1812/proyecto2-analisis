//
// Created by mike on 12/04/18.
//

#ifndef TAREA03_ROBOT_H
#define TAREA03_ROBOT_H

#include "LUDoolittle.hpp"

namespace anpi {


    /**
     * Mapea los indices de los dos nodos terminales de una resistencia a un indice
     * lineal para acceder al vector x.
     * @tparam T
     * @param m1 [in] fila del primer nodo.
     * @param n1 [in] columna del primer nodo.
     * @param m2 [in] fila del segundo nodo.
     * @param M [in] Cantidad de filas de la matriz.
     * @param N [in] Cantidad de columnas de la matriz.
     * @param x [out] posicion en el vector;
     */
    void mapeoMatrizVector(const unsigned int& m1, const unsigned int& n1,
                           const unsigned int& m2,
                           const unsigned int& M, const unsigned int& N,
                           unsigned int& x);

    /**
     * Mapea un indice en el vector x en las coordenadas en la matriz de ambos nodos.
     * @tparam T
     * @param m1 [out] fila del primer nodo.
     * @param n1 [out] columna del primer nodo.
     * @param m2 [out] fila del segundo nodo.
     * @param n2 [out] columna del segundo nodo.
     * @param M [in] Cantidad de filas de la matriz.
     * @param N [in] Cantidad de columnas en la matriz.
     * @param x [in] posixion en el arreglo x.
     */
    void mapeoVectorMatriz(unsigned int& m1, unsigned int& n1,
                           unsigned int& m2, unsigned int& n2,
                           const unsigned int& M, const unsigned int& N,
                           const unsigned int& x);

    /**
     *
     * @tparam T
     * @param M [out] matriz de factores de las ecuaciones
     * @param r [in] Arreglo con los valores de resistencias
     * @param Rm [in] Filas de la matriz de nodos
     * @param Rn [in] Columnas de la matriz de nodos
     */
    template <typename T>
    void ecuacionesMallas(anpi::Matrix<T>& M, const std::vector<T>& r,
                          const unsigned int Rm, const unsigned int Rn) {
        unsigned int posEcuacion = (Rm * Rn-1);
        unsigned int x = 0;

        for (unsigned int m = 0; m < Rm-1; ++m) {
            for (unsigned int n = 0; n < Rn-1; ++n) {
                anpi::mapeoMatrizVector(m, n, m, Rm, Rn, x);
                M[posEcuacion][x] = r[x];
                anpi::mapeoMatrizVector(m, n+1, m+1, Rm, Rn, x);
                M[posEcuacion][x] = r[x];
                anpi::mapeoMatrizVector(m, n, m+1, Rm, Rn, x);
                M[posEcuacion][x] = -r[x];
                anpi::mapeoMatrizVector(m+1, n, m+1, Rm, Rn, x);
                M[posEcuacion][x] = -r[x];
                ++posEcuacion;
            }
        }
    }

    /**
     * Genera las ecuaciones para los nodos
     * @tparam T
     * @param M [out] matriz de factores de las ecuaciones
     * @param b [out] lado derecho de las ecuaciones
     * @param Rm [in] Filas de la matriz de nodos
     * @param Rn [in] Columnas de la matriz de nodos
     */
    template <typename T>
    void ecuacionesNodos(anpi::Matrix<T>& M, const std::vector<T>& b,
                         const unsigned int Rm, const unsigned int Rn){

        unsigned int posIr = 0; //posicion en el arreglo x
        unsigned int posEcuacion = 0;
        bool excluida = false;

        for(unsigned int m = 0; m < Rm; ++m) {
            for(unsigned int n = 0; n < Rn; ++n) {
                if (m == 0) {
                    if (n == 0) { // Esquina superior izquierda
                        if (!excluida && b[posEcuacion] == 0) {
                            excluida = true;
                            continue;
                        }
                        anpi::mapeoMatrizVector(m, n, 1, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    } else if (n == Rn - 1) { // Esquina superior derecha
                        if (!excluida && b[posEcuacion] == 0) {
                            excluida = true;
                            continue;
                        }
                        anpi::mapeoMatrizVector(m, n, 1, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                        anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                    } else { // Borde superior
                        anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m, n, m+1, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    }
                } else if (n == 0) {
                    if (m == Rm - 1) { // Esquina inferior izquierda
                        if (!excluida && b[posEcuacion] == 0) {
                            excluida = true;
                            continue;
                        }
                        anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    } else { // Borde izquierdo
                        anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                        anpi::mapeoMatrizVector(m, n, m+1, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    }
                } else if (m == Rm - 1) {
                    if (n == Rn - 1) { // Esquina inferior derecha
                        if (!excluida && b[posEcuacion] == 0) {
                            excluida = true;
                            continue;
                        }
                        anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                    } else { // Borde inferior
                        anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = 1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    }
                } else if (n == Rn - 1) { // Borde derecho
                    anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = 1;
                    anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = 1;
                    anpi::mapeoMatrizVector(m, n, m+1, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = -1;
                } else { // Internos
                    anpi::mapeoMatrizVector(m-1, n, m, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = 1;
                    anpi::mapeoMatrizVector(m, n-1, m, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = 1;
                    anpi::mapeoMatrizVector(m, n, m+1, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = -1;
                    anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                    M[posEcuacion][posIr] = -1;
                }
                ++posEcuacion;
            }
        }
    }
}


#endif //TAREA03_ROBOT_H
