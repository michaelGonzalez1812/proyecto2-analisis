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
#include "PlotPy.hpp"

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
    int unsigned i;
    std::cout << "{";
    for (i = 0; i < A.size(); i++) {
        if (i == A.size() - 1) std::cout << "i[" << i << "]: " << A[i];
        else std::cout << "i[" << i << "]: " << A[i] << ", ";
        if (i == 19) std::cout << std::endl;
    }
    std::cout << "}" << std::endl << i << std::endl;
}

void castVector(const std::vector<unsigned int>& x, std::vector<int>& y) {
    for (unsigned int i = 0; i < x.size(); ++i) {
        y.push_back(x[i]);
    }
}


int main() {
    anpi::Matrix<float> M(38, 38); //Matriz para los coeficientes de las ecuacines
    std::vector<float> i(38); //vector de corrientes
    std::vector<float> b(38); //parte derecha de las ecuaciones
    const unsigned int minicial = 0;
    const unsigned int ninicial = 1;
    const unsigned int mfinal = 3;
    const unsigned int nfinal = 4;
    const unsigned int filasMatrizNodos = 4; //Tamaño de la matriz de nodos
    const unsigned int columnasMatrizNodos = 6; //Tamaño de la matriz de nodos

//Descomentar para trabajar con estrategia 1--------------------------------------------
    std::vector<unsigned int> x;
    std::vector<unsigned int> y;
//--------------------------------------------------------------------------------------

//Descomentar para trabajar con estrategia 2--------------------------------------------
//    const float alpha = 1; //tamaño del paso
//    std::vector<float> x;
//    std::vector<float> y;
//    anpi::Matrix<float> dx(4, 6);
//    anpi::Matrix<float> dy(4, 6);
//--------------------------------------------------------------------------------------



    std::vector<unsigned int> permut(38);
    anpi::Matrix<float> LU;
    std::vector<float> r { 1000, 1000, 1000, 1000, 1000,
                           1000, 1000, 1000, 1000, 1000,
                           1000, 1, 1, 1, 1000,
                           1000, 1000, 1000, 1000, 1000,
                           1000, 1000, 1000,
                           1, 1, 1000,
                           1000, 1000, 1000,
                           1000, 1000, 1000,
                           1000, 1000, 1,
                           1000, 1000, 1000};

    anpi::ecuacionesNodos(M, b, 4, 6, 0, 1, 3, 4);
    anpi::ecuacionesMallas(M, r, 4, 6);

    anpi::luDoolittle(M, LU, permut);
    anpi::solveLU(LU, i, b, permut);

//Descomentar para trabajar con estrategia 1--------------------------------------------
    anpi::estrategia1(i, minicial, ninicial, mfinal, nfinal, filasMatrizNodos, columnasMatrizNodos, x, y);
//--------------------------------------------------------------------------------------

//Descomentar para trabajar con estrategia 2--------------------------------------------
//    anpi::desplazamientos(i, dx, dy, 4, 6);
//    anpi::estrategia2(minicial, ninicial, mfinal, nfinal, filasMatrizNodos, columnasMatrizNodos, dx, dy, alpha, x, y);
//--------------------------------------------------------------------------------------

    std::cout<< std::endl << "ruta:    ";
    for (unsigned int j = 0; j < x.size(); ++j) {
        std::cout << "(" << x[j] << ", " << y[j] << ") ";
    }

    // Codigo de ejemplo de graficacion
    std::vector<int> c;
    std::vector<int> d;

    anpi::Plot2d <int> p;
    p.initialize(001);
    p.setTitle("Trayectoria del Robot");
    p.setXLabel("Desplazamiento en X");
    p.setYLabel("Desplazamiento en Y");
    castVector(x,c);
    castVector(y,d);
    p.plot(c,d,"Trayectoria","green");
    p.show();

    return EXIT_FAILURE;
}