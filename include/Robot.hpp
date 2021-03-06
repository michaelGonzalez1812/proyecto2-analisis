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
    void ecuacionesNodos(anpi::Matrix<T>& M, std::vector<T>& b,
                         const unsigned int Rm, const unsigned int Rn,
                         const unsigned int minicial, const unsigned int ninicial,
                         const unsigned int mfinal, const unsigned int nfinal){

        unsigned int posIr = 0; //posicion en el arreglo x
        unsigned int posEcuacion = 0;
        bool excluida = false;

        for(unsigned int m = 0; m < Rm; ++m) {
            for(unsigned int n = 0; n < Rn; ++n) {
                if (m == minicial && n == ninicial) b[posEcuacion] = -1;
                else if (m == mfinal && n == nfinal) b[posEcuacion] = 1;
                else b[posEcuacion] = 0;
                if (m == 0) {
                    if (n == 0) { // Esquina superior izquierda
                        if (!excluida && (m != minicial || n != ninicial) && (m != mfinal && n != nfinal)) {
                            excluida = true;
                            continue;
                        }
                        anpi::mapeoMatrizVector(m, n, 1, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                        anpi::mapeoMatrizVector(m, n, m, Rm, Rn, posIr);
                        M[posEcuacion][posIr] = -1;
                    } else if (n == Rn - 1) { // Esquina superior derecha
                        if (!excluida && (m != minicial || n != ninicial) && (m != mfinal && n != nfinal)) {
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
                        if (!excluida && (m != minicial || n != ninicial) && (m != mfinal && n != nfinal)) {
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
                        if (!excluida && (m != minicial || n != ninicial) && (m != mfinal && n != nfinal)) {
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

    /**
     * Determina un path entre el nodo inicial y final tomando como criterio la direccion
     * de mayor corriente en cada nodo.
     * @tparam T
     * @param I [in] corrientes
     * @param permut [in] orden del arreglo de corrientes
     * @param initialM [in] Posicion inicial en "m"
     * @param initialN [in] Posicion inicial en "n"
     * @param finalM [in] Posicion final en "n"
     * @param finalN [in] Posicion final en "n"
     * @param M [in] Cantidad de filas en la matriz de nodos
     * @param N [in] Cantidad de columnas en la matriz de nodos
     * @param x [out] arreglo de posiciones x del path
     * @param y [out] Arreglo de posiciones y del path
     */
    template <typename T>
    void estrategia1(const std::vector<T>& I,
                     const int& initialM, const int& initialN,
                     const unsigned int& finalM, const unsigned int& finalN,
                     const unsigned int& M, const unsigned int& N,
                     std::vector<unsigned int>& x,
                     std::vector<unsigned int>& y) {

        x = {(unsigned int)initialM};
        y = {(unsigned int)initialN};
        unsigned int arriba;
        unsigned int derecha;
        unsigned int abajo;
        unsigned int izquierda;

        while (x.back() != finalM || y.back() != finalN) {
            if (x.back() == 0) {
                if (y.back() == 0) { // Esquina superior izquierda
                    if (I[0] < 0) {
                        x.push_back(1);
                        y.push_back(0);
                    } else {
                        x.push_back(0);
                        y.push_back(1);
                    }
                } else if (y.back() == N - 1) { // Esquina superior derecha
                    anpi::mapeoMatrizVector(x.back(), y.back()-1, x.back(), M, N, izquierda);
                    if (I[izquierda] < 0) {
                        x.push_back(0);
                        y.push_back(y.back()-1);
                    } else {
                        x.push_back(1);
                        y.push_back(y.back());
                    }
                } else { // Borde superior
                    anpi::mapeoMatrizVector(0, y.back()-1, 0, M, N, izquierda);
                    anpi::mapeoMatrizVector(0, y.back(), 1, M, N, abajo);
                    anpi::mapeoMatrizVector(0, y.back(), 0, M, N, derecha);

                    if (-I[izquierda] > I[abajo] && -I[izquierda] > I[derecha]) {
                        x.push_back(0);
                        y.push_back(y.back()-1);
                    } else if (I[abajo] > I[derecha]) {
                        x.push_back(1);
                        y.push_back(y.back());
                    } else {
                        x.push_back(0);
                        y.push_back(y.back()+1);
                    }
                }
            } else if (y.back() == 0) {
                if (x.back() == M - 1) { // Esquina inferior izquierda
                    mapeoMatrizVector(x.back(), 0, x.back(), M, N, derecha);
                    if (I[derecha] > 0){
                        x.push_back(x.back());
                        y.push_back(1);
                    } else {
                        x.push_back(x.back()-1);
                        y.push_back(0);
                    }
                } else { // Borde izquierdo
                    mapeoMatrizVector(x.back()-1, 0, x.back(), M, N, arriba);
                    mapeoMatrizVector(x.back(), 0, x.back(), M, N, derecha);
                    mapeoMatrizVector(x.back(), 0, x.back()+1, M, N, abajo);

                    if (-I[arriba] > I[derecha] && -I[arriba] > I[abajo]) {
                        x.push_back(x.back()-1);
                        y.push_back(0);
                    } else if (I[derecha] > I[abajo]) {
                        x.push_back(x.back());
                        y.push_back(1);
                    } else {
                        x.push_back(x.back()+1);
                        y.push_back(0);
                    }
                }
            } else if (x.back() == M - 1) {
                if (y.back() == N - 1) { // Esquina inferior derecha
                    anpi::mapeoMatrizVector(x.back()-1, y.back(), x.back(), M, N, arriba);
                    if (I[arriba] > 0) {
                        x.push_back(x.back());
                        y.push_back(y.back()-1);
                    } else {
                        x.push_back(x.back()-1);
                        y.push_back(y.back());
                    }
                } else { // Borde inferior
                    anpi::mapeoMatrizVector(x.back(), y.back()-1, x.back(), M, N, izquierda);
                    anpi::mapeoMatrizVector(x.back()-1, y.back(), x.back(), M, N, arriba);
                    anpi::mapeoMatrizVector(x.back(), y.back(), x.back(), M, N, derecha);
                    if (-I[izquierda] > -I[arriba] && -I[izquierda] > I[derecha]) {
                        x.push_back(x.back());
                        y.push_back(y.back()-1);
                    } else if (-I[arriba] > I[derecha]) {
                        x.push_back(x.back()-1);
                        y.push_back(y.back());
                    } else {
                        x.push_back(x.back());
                        y.push_back(y.back()+1);
                    }
                }
            } else if (y.back() == N - 1) { // Borde derecho
                anpi::mapeoMatrizVector(x.back()-1, y.back(), x.back(), M, N, arriba);
                anpi::mapeoMatrizVector(x.back(), y.back()-1, x.back(), M, N, izquierda);
                anpi::mapeoMatrizVector(x.back(), y.back(), x.back()+1, M, N, abajo);
                if (-I[arriba] > -I[izquierda] && -I[arriba] > I[abajo]) {
                    x.push_back(x.back()-1);
                    y.push_back(y.back());
                } else if (-I[izquierda] > I[abajo]) {
                    x.push_back(x.back());
                    y.push_back(y.back()-1);
                } else {
                    x.push_back(x.back()+1);
                    y.push_back(x.back());
                }
            } else { // Internos
                anpi::mapeoMatrizVector(x.back()-1, y.back(), x.back(), M, N, arriba);
                anpi::mapeoMatrizVector(x.back(), y.back(), x.back(), M, N, derecha);
                anpi::mapeoMatrizVector(x.back(), y.back(), x.back()+1, M, N, abajo);
                anpi::mapeoMatrizVector(x.back(), y.back()-1, x.back(), M, N, izquierda);
                if(-I[arriba] > I[derecha] && -I[arriba] > I[abajo] && -I[arriba] > -I[izquierda]) {
                    x.push_back(x.back()-1);
                    y.push_back(y.back());
                } else if (I[derecha] > I[abajo] && I[derecha] > -I[izquierda]) {
                    x.push_back(x.back());
                    y.push_back(y.back()+1);
                } else if (I[abajo] > -I[izquierda]) {
                    x.push_back(x.back()+1);
                    y.push_back(y.back());
                } else {
                    x.push_back(x.back());
                    y.push_back(y.back()-1);
                }
            }
        }
    }

    template <typename T>
    void desplazamientos(const std::vector<T>& I,
                         anpi::Matrix<T>& dx,
                         anpi::Matrix<T>& dy,
                         const unsigned int& M, const unsigned int& N) {

        unsigned int arriba;
        unsigned int derecha;
        unsigned int abajo;
        unsigned int izquierda;

        for (unsigned int m = 0; m < M; ++m) {
            for (unsigned int n = 0; n < N; ++n) {
                if (m == 0) {
                    if (n == 0) { // Esquina superior izquierda
                        anpi::mapeoMatrizVector(0, 0, 1, M, N, abajo);
                        dx[m][n] = I[0];
                        dy[m][n] = I[abajo];
                    } else if (n == N - 1) { // Esquina superior derecha
                        anpi::mapeoMatrizVector(m, n-1, m, M, N, izquierda);
                        anpi::mapeoMatrizVector(m, n, 1, M, N, abajo);
                        dx[m][n] = I[izquierda];
                        dy[m][n] = I[abajo];
                    } else { // Borde superior
                        anpi::mapeoMatrizVector(0, n-1, 0, M, N, izquierda);
                        anpi::mapeoMatrizVector(0, n, 1, M, N, abajo);
                        anpi::mapeoMatrizVector(0, n, 0, M, N, derecha);
                        dx[m][n] = I[izquierda] + I[derecha];
                        dy[m][n] = I[abajo];
                    }
                } else if (n == 0) {
                    if (m == M - 1) { // Esquina inferior izquierda
                        mapeoMatrizVector(m, 0, m, M, N, derecha);
                        mapeoMatrizVector(m-1, 0, m, M, N, arriba);
                        dx[m][n] = I[derecha];
                        dy[m][n] = I[arriba];
                    } else { // Borde izquierdo
                        mapeoMatrizVector(m-1, 0, m, M, N, arriba);
                        mapeoMatrizVector(m, 0, m, M, N, derecha);
                        mapeoMatrizVector(m, 0, m+1, M, N, abajo);
                        dx[m][n] = I[derecha];
                        dy[m][n] = I[arriba] + I[abajo];
                    }
                } else if (m == M - 1) {
                    if (n == N - 1) { // Esquina inferior derecha
                        anpi::mapeoMatrizVector(m-1, n, m, M, N, arriba);
                        anpi::mapeoMatrizVector(m, n-1, m, M, N, izquierda);
                        dx[m][n] = I[izquierda];
                        dy[m][n] = I[arriba];
                    } else { // Borde inferior
                        anpi::mapeoMatrizVector(m, n-1, m, M, N, izquierda);
                        anpi::mapeoMatrizVector(m-1, n, m, M, N, arriba);
                        anpi::mapeoMatrizVector(m, n, m, M, N, derecha);
                        dx[m][n] = I[izquierda] + I[derecha];
                        dy[m][n] = I[arriba];
                    }
                } else if (n == N - 1) { // Borde derecho
                    anpi::mapeoMatrizVector(m-1, n, m, M, N, arriba);
                    anpi::mapeoMatrizVector(m, n-1, m, M, N, izquierda);
                    anpi::mapeoMatrizVector(m, n, m+1, M, N, abajo);
                    dx[m][n] = I[izquierda];
                    dy[m][n] = I[arriba] + I[abajo];
                } else { // Internos
                    anpi::mapeoMatrizVector(m-1, n, m, M, N, arriba);
                    anpi::mapeoMatrizVector(m, n, m, M, N, derecha);
                    anpi::mapeoMatrizVector(m, n, m+1, M, N, abajo);
                    anpi::mapeoMatrizVector(m, n-1, m, M, N, izquierda);
                    dx[m][n] = I[izquierda] + I[derecha];
                    dy[m][n] = I[arriba] + I[abajo];
                }
            }
        }
    }

    template <typename T>
    void estrategia2(const int& initialM, const int& initialN,
                     const int& finalM, const int& finalN,
                     const unsigned int& M, const unsigned int& N,
                     const Matrix<T>& dx, const Matrix<T>& dy,
                     const T& alpha,
                     std::vector<T, std::allocator<T> >& x, std::vector<T, std::allocator<T> >& y) {
        x = {initialM};
        y = {initialN};
        int x_0;
        int y_0;
        T x_c; //componente x de la interpolacion en y
        T y_c; //componente y de la interpolacion en y

        //puede dar problemas por los tipos de las variables
        while (std::abs(x.back() - finalM) > 0.5 || std::abs(y.back() - finalN) > 0.5) {
            x_0 = x.back();
            y_0 = y.back();
            if (x.back() - x_0 != 0 && y.back() - y_0 != 0) { //interpolar
                x_c = (dx[x_0][y_0] + dx[x_0+1][y_0] + dx[x_0][y_0+1] + dx[x_0+1][y_0+1])/4;
                y_c = (dy[x_0][y_0] + dy[x_0+1][y_0] + dy[x_0][y_0+1] + dy[x_0+1][y_0+1])/4;
            } else {
                x_c = dx[x_0][y_0];
                y_c = dy[x_0][y_0];
            }
            x.push_back(x.back()+x_c*alpha);
            y.push_back(y.back()+y_c*alpha);
        }
    }

}
#endif //TAREA03_ROBOT_H