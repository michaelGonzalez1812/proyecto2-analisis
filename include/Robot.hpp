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
}


#endif //TAREA03_ROBOT_H
