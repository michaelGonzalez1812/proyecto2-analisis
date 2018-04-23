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
#include "Robot.hpp"

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
//    anpi::Matrix<float> A = {{1,  2, -1},
//                             {3,  1,  1},
//                             {1, -1,  2}};
//    anpi::Matrix<float> LU;
//    std::vector<unsigned int> permut(3);
//    std::vector<float> x(3);
//    std::vector<float> b{-3, 4, 6};
//
//    anpi::luDoolittle(A, LU, permut);
//
//    std::cout << std::endl;
//    printMatriz(LU);
//
//    anpi::solveLU(LU, x, b, permut);
//
//    std::cout << "orden:   " << std::endl;
//    printVector(permut);
//    std::cout << std::endl << "resultado" << std::endl;
//    printVector(x);
//-----------------------------------------------------------

//mapeo matriz vector----------------------------------------
//    const unsigned int m1 = 1;
//    const unsigned int n1 = 0;
//    const unsigned int m2 = 1;
//    const unsigned int M = 5;
//    const unsigned int N = 4;
//    unsigned int x = 0;
//    anpi::mapeoMatrizVector(m1, n1, m2, M, N, x);
//
//    std::cout << x;
//------------------------------------------------------------


//mapeo vector matriz-----------------------------------------
//    unsigned int m1 = 0;
//    unsigned int n1 = 0;
//    unsigned int m2 = 0;
//    unsigned int n2 = 0;
//    const unsigned int M = 4;
//    const unsigned int N = 3;
//    unsigned int x = 7;
//    anpi::mapeoVectorMatriz(m1, n1, m2, n2, M, N, x);
//
//    std::cout << "m1= " << m1 << "   n1: " <<  n1 << std::endl;
//    std::cout << "m2= " << m2 << "   n2: " <<  n2 << std::endl;
//--------------------------------------------------------------

//Montar ecuaciones---------------------------------------------
    anpi::Matrix<float> M(7, 7);
    std::vector<float> x(7);
    std::vector<float> b {-1, 0, 0, 1, 0, 0, 0};
    std::vector<unsigned int> permut(7);
    anpi::Matrix<float> LU;
    std::vector<float> r {10, 10, 100, 10, 100, 100, 10};

    anpi::ecuacionesNodos(M, b, 2, 3);
    anpi::ecuacionesMallas(M, r, 2, 3);

    anpi::luDoolittle(M, LU, permut);

    anpi::solveLU(LU, x, b, permut);

    printMatriz(M);
    std::cout << std::endl << "x:   ";
    printVector(x);
//---------------------------------------------------------------
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
