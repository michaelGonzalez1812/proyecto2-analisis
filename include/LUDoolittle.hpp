/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: 
 * @Date  : 03.03.2018
 */

#include <cmath>
#include <limits>
#include <functional>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"

#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {


    /**
     * Auxiliary method used to debug LU decomposition.
     *
     * It separates a packed LU matrix into the lower triangular matrix
     * L and the upper triangular matrix U, such that the diagonal of L
     * is composed by 1's.
     */
    template<typename T>
    void unpackDoolittle(const Matrix<T> &LU,
                         Matrix<T> &L,
                         Matrix<T> &U) {
        if (LU.cols() == LU.rows()) {
            L = LU;
            U = LU;
            for (unsigned int i = 0; i < LU.rows(); i++) {
                for (unsigned int j = 0; j < LU.cols(); j++) {
                    if (i == j) { // estamos en la diagonal
                        L[i][j] = 1;
                        U[i][j] = LU[i][j];
                    } else if (i < j) { // matriz triangular superior
                        L[i][j] = 0;
                        U[i][j] = LU[i][j];
                    } else { // matriz triangular inferior
                        L[i][j] = LU[i][j];
                        U[i][j] = 0;
                    }
                }
            }
        } else throw anpi::Exception("La matriz no es cuadrada");
    }

    /**
     *-o
     * a single matrix LU.
     *
     * The L matrix will have in the Doolittle's LU decomposition a
     * diagonal of 1's
     *
     * @param[in] A a square matrix
     * @param[out] LU matrix encoding the L and U matrices
     * @param[out] permut permutation vector, holding the indices
     * of the
     *             original matrix falling into the corresponding element.
     *             For example if permut[5]==3 holds, then the fifth row
     *             of the LU decomposition in fact is dealing with the third
     *             row of the original matrix.
     *
     * @throws anpi::Exception if matrix cannot be decomposed, or input
     *         matrix is not square.
     */
    template<typename T>
    void luDoolittle(const Matrix<T> &A,
                     Matrix<T> &LU,
                     std::vector<size_t> &permut) {
        if (A.rows() == A.cols()) {
            LU = A;
            T cte = T(0);
            for (unsigned int k = 0; k < A.rows() - 1; k++) { // iteracion para crear la matriz triangular superior
                for (unsigned int i = k + 1; i < A.rows(); i++) {
                    for (unsigned int j = k; j < A.rows(); j++) {
                        if (LU[k][k] != 0) { // evitar division entre cero
                            if (j == k) {
                                cte = LU[i][k] / LU[k][k]; // cte para generar ceros en la matriz
                                LU[i][j] = cte; // creacion de L
                            } else {
                                LU[i][j] = LU[i][j] - cte * LU[k][j]; // creacion de U
                            }

                        }

                    }

                }

            }
        } else throw anpi::Exception("La matriz no es cuadrada");
    }

    template<typename T>
    void pivot(const Matrix<T> &A, std::vector<T>& permut, std::vector<T>& mayores, const int& k){
        unsigned int pos = k;
        T aux;
        T big = fabs(A[permut[k]][k] / mayores[permut[k]]); //se escala el valor -> s[permut[k]] = el mayor de la fila

        //busca el mayor de la columna
        for (unsigned int i = k+1; i < A.rows(); ++i) {
            aux = fabs(A[permut[i]][k] / mayores[permut[i]]);
            if (aux > big) {
                big = aux;
                pos = i;
            }
        }

        //switch
        aux = permut[pos];
        permut[pos] = permut[k];
        permut[k] = aux;
    }
}

#endif

