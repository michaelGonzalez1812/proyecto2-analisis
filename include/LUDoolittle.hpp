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
#include <vector>

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
                     std::vector<unsigned int> &permut) {

        const T tol = 0.0001; //tolerancia para division entre 0

        if (A.rows() == A.cols()) {
            std::vector<T> bigxrow(A.rows());
            LU = A;
            T cte = T(0);

            //encuentra los valores mas grandes de cada fila
            for (unsigned int i = 0; i < A.rows(); ++i) {
                permut[i] = i;
                bigxrow[i] = fabs(A[i][0]);
                for (unsigned int j = 1; j < A.rows(); ++j) {
                    if(fabs(A[i][j]) > bigxrow[i]) bigxrow[i] = fabs(A[i][j]);
                }
            }

            for (unsigned int k = 0; k < A.rows() - 1; ++k) { // iteracion para crear la matriz triangular superior
                pivot(LU, permut, bigxrow, k);
                for (unsigned int i = k + 1; i < A.rows(); ++i) {

                    if (fabs(LU[permut[k]][k]) > fabs(tol))             // evitar division entre cero
                        cte = LU[permut[i]][k] / LU[permut[k]][k]; // cte para generar ceros en la matriz
                    else
                        throw anpi::Exception("Division entre cero");
                    LU[permut[i]][k] = cte; // creacion de L
                    for (unsigned int j = k+1; j < A.rows(); ++j)
                        LU[permut[i]][j] = LU[permut[i]][j] - cte * LU[permut[k]][j]; // creacion de U
                }
            }
        } else throw anpi::Exception("La matriz no es cuadrada");
    }

    /**
     *
     * @param A [in]
     * @param permut [out] order of the rows
     * @param mayores [in] Bigger abslutes values by each row.
     * @param k [in] Column.
     */
    template<typename T>
    void pivot(const Matrix<T> &A, std::vector<unsigned int>& permut, std::vector<T>& mayores, const int& k){
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

    /**
     * This method solve a the system using the LU decomposition of the matrix.
     * @param A [in]
     * @param x [out] solution
     * @param b [in] elements on right side of the equal
     * @param permut [in] order of the rows
     */
    template <typename T>
    void solveLU(const anpi::Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b,
                 const std::vector<unsigned int> &permut) {
        T sum;
        std::vector<T> d = b;
        for (unsigned int i = 1; i < A.rows(); ++i) {
            sum = d[permut[i]];
            for (unsigned int j = 0; j <= i-1; ++j) {
                sum = sum - A[permut[i]][j] * d[permut[j]];
            }
            d[permut[i]] = sum;
        }

        x[A.rows()-1] = d[permut[A.rows()-1]] / A[permut[A.rows()-1]][A.rows()-1];
        for (int i = A.rows()-2; i >= 0; --i) {
            sum = 0;
            for (unsigned int j = i+1; j < A.rows(); ++j) {
                sum = sum + A[permut[i]][j] * x[j];
            }
            x[i] = (d[permut[i]] - sum) / A[permut[i]][i];
        }
    }
}

#endif

