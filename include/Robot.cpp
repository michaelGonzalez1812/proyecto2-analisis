//
// Created by mike on 12/04/18.
//

#include "Robot.hpp"




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
    void anpi::mapeoMatrizVector(const unsigned int& m1, const unsigned int& n1,
                           const unsigned int& m2,
                           const unsigned int& Mm, const unsigned int& Nn,
                           unsigned int& x){
        unsigned int M = Mm - 1;
        unsigned int N = Nn - 1;
        x = 0;
        if (m1 == m2)
            x = m1 * N + n1;
        else{
            x = M * N + N;
            x += n1 * M + m1;
        }
    }

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
    void anpi::mapeoVectorMatriz(unsigned int& m1, unsigned int& n1,
                           unsigned int& m2, unsigned int& n2,
                           const unsigned int& Mm, const unsigned int& Nn,
                           const unsigned int& x){
        unsigned int M = Mm - 1;
        unsigned int N = Nn - 1;
        unsigned int inicioCol = M * N + N;

        if (x < inicioCol){
            n1 = x % N;
            n2 = n1 + 1;
            m1 = (x - n1) / N;
            m2 = m1;
        } else {
            unsigned int pos = x - inicioCol;
            m1 = pos % M;
            m2 = m1;
            n1 = (pos - m1) / M;
            n2 = n1;
        }
    }
