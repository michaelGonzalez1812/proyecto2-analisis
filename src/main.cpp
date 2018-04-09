/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <LUCrout.hpp>
#include "LUDoolittle.hpp"

/**
 * Herramienta para visualizar matrices
 * @tparam T
 * @param A
 */
template<typename T>
void printMatriz(anpi::Matrix<T> &A) {
    for (int unsigned i = 0; i < A.rows(); i++) {
        std::cout << "{ ";
        for (unsigned int j = 0; j < A.cols(); j++) {
            if (j == A.cols() - 1) std::cout << A[i][j];
            else std::cout << A[i][j] << ", ";
        }
        std::cout << " }" << std::endl;
    }
}

/**
 * Herramienta para visualizar vectores
 * @tparam T
 * @param A
 */
template<typename T>
void printVector(std::vector<T> &A) {
    std::cout << "{";
    for (int unsigned i = 0; i < A.size(); i++) {
        if (i == A.size() - 1) std::cout << A[i];
        else std::cout << A[i] << ", ";
    }
    std::cout << "}" << std::endl;
}

int main() {
//***********************************************************
//    Some example code
//***********************************************************

//LUDoolittle con pivoteo------------------------------------
    anpi::Matrix<float> A = {{2.5, -4, 5.3},
                             {6.8,  4, 3.6},
                             {3,    3, 3  }};
    anpi::Matrix<float> LU;
    std::vector<unsigned int> permut(3);

    anpi::luDoolittle(A, LU, permut);
    std::cout << "orden:   " << std::endl;

    printVector(permut);
    std::cout << std::endl;
    printMatriz(LU);
//-----------------------------------------------------------

    return EXIT_FAILURE;
}

//    anpi::Matrix<float> A = {{-1, -2, 1, 2},
//                             {2,  0,  1, 2},
//                             {-1, -1, 0, 1},
//                             {1,  1,  1, 1}};
//    std::cout << "------------TEST LU------------" << std::endl;
//    std::cout << "-Matriz A original-" << std::endl;
//    printMatriz(A);
//    anpi::Matrix<float> LU;
//    std::vector<size_t> p;
//    std::cout << "- LU Dolittle de A -" << std::endl;
//    anpi::luDoolittle(A, LU, p);
//    printMatriz(LU);
//    std::cout << "- Unpack Dolittle de LU -" << std::endl;
//    anpi::Matrix<float> L, U, LUA;
//    anpi::unpackDoolittle(LU, L, U);
//    std::cout << "L->" << std::endl;
//    printMatriz(L);
//    std::cout << "U->" << std::endl;
//    printMatriz(U);
//    std::cout << "- Comprobacion de A=L*U -" << std::endl;
//    LUA = L * U;
//    printMatriz(LUA);
//    std::cout << "-----------------------------" << std::endl;
//    std::cout << "- LU Crout de A -" << std::endl;
//    anpi::luCrout(A, LU, p);
//    printMatriz(LU);
//    std::cout << "- Unpack Crout de LU -" << std::endl;
//    anpi::unpackCrout(LU, L, U);
//    std::cout << "L->" << std::endl;
//    printMatriz(L);
//    std::cout << "U->" << std::endl;
//    printMatriz(U);
//    std::cout << "- Comprobacion de A=L*U -" << std::endl;
//    LUA = L * U;
//    printMatriz(LUA);
//    std::cout << "-------TEST Producto Matricial--------" << std::endl;
//    std::cout << "-Matriz A -" << std::endl;
//    printMatriz(A);
//    std::cout << "Vector O -" << std::endl;
//    std::vector<float> O = {1, 2, 3, 4};
//    printVector(O);
//    std::cout << "-A*O -" << std::endl;
//    std::vector<float> R = A * O;
//    printVector(R);
//    std::cout << "-----------------------------" << std::endl;
//    std::cout << "A*A" << std::endl;
//    LUA=A*A;
//    printMatriz(LUA);
