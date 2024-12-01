#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX 100

#define MAX_DIM 100

void mostrarMenu() {
    printf("--------- Metodos Numericos ----------\n");
    printf("------------ Integrantes -------------\n");
    printf("Mariana Alejandra Lopez Ramirez\n");
    printf("Jessica Esmeralda Alcantar Hernandez\n");
      printf("Andrea Jaimes Molina\n");
}


double f(double x, int funcIndex);
void metodoDeBiseccion(double (*func)(double, int), double a, double b, int iteraciones, double tolerancia);
void secante(double a, double b, double epsilon, int imax, int funcIndex);


// Esta función implementa el método de bisección
	void metodoDeBiseccion(double (*func)(double, int), double a, double b, int iteraciones, double tolerancia) {
	    double p, fp, fa, fb;
	    int iter = 0;
	    double Er = 1.0;
	
	    printf("Iter\t a\t\t b\t\t p\t\t f(p)\t\t Error relativo\n");
	    while (iter < iteraciones && Er > tolerancia) {
	        p = (a + b) / 2.0;
	        fa = func(a, 1);
	        fb = func(b, 1);
	        fp = func(p, 1);
	
	        printf("%d\t %.6lf\t %.6lf\t %.6lf\t %.6lf\t %.6lf\n", iter + 1, a, b, p, fp, Er);
	
	        if (iter > 0) {
	            Er = fabs((b - a) / p);
	        }
	
	        if (fa * fp < 0) {
	            b = p;
	        } else {
	            a = p;
	        }
	        iter++;
	    }
	
	    printf("\nLa raiz obtenida es: %.6lf\n", p);
	    printf("Se alcanzo en la iteracion: %d\n", iter);
	    printf("Con una tolerancia de: %.6lf\n", Er);
	}
	
void secante(double a, double b, double epsilon, int imax, int funcIndex) {
    double x0 = a;
    double x1 = b;
    double f0 = f(x0, funcIndex);
    double f1 = f(x1, funcIndex);
    double x2, error;

    // Encabezado de la tabla
    printf("\nIter\t x0\t\t f(x0)\t\t x1\t\t f(x1)\t\t x2\t\t f(x2)\t\t Error\n");
	int i;
    for (i = 1; i <= imax; i++) {
        // Comprobar que no hay división por cero
        if (f1 - f0 == 0) {
            printf("Error: Division por cero en la iteracion %d.\n", i);
            return;
        }

        // Calcular la siguiente aproximación
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        double f2 = f(x2, funcIndex);
        error = fabs(f2);

        // Imprimir la fila de la tabla
        printf("%d\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\t %.10f\n", 
               i, x0, f0, x1, f1, x2, f2, error);

        // Verificar si se ha alcanzado la convergencia
        if (error < epsilon) {
            printf("Se encontro una raiz en x = %.10f\n", x2);
            printf("En %d iteraciones\n", i);
            printf("El valor de la funcion es %.10f\n", f2);
            return;
        }

        // Actualizar para la siguiente iteración
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }

    printf("Despues de %d iteraciones, no se encontro ninguna raiz dentro del criterio de convergencia\n", imax);
}

// Función para evaluar f(x) según el índice de función
double f(double x, int funcIndex) {
    switch (funcIndex) {
        case 1:
            return (x * x * cos(x) - 2 * x);
        case 2:
            if (x != 0)
                return (6 - 2 / (x * x)) * (exp(2 + x) / 4) + 1;
            else
                return INFINITY; // Para x=0, devolver infinito
        case 3:
            return (x * x * x - 3 * sin(x * x) + 1);
        case 4:
            return (x * x * x + 6 * x * x + 9.4 * x + 2.5);
        default:
            printf("indice de funcion no valido\n");
            return 0;
    }
}

///////////////////////////PROGRAMA 2//////////////////
void leermatriz(int n, double matriz[][n], double vector[]) {
    printf("Introduce los coeficientes de la matriz:\n");
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("Elemento [%d][%d]: ", i + 1, j + 1);
            if (scanf("%lf", &matriz[i][j]) != 1) {
                printf("Entrada invalida.\n");
                exit(1);
            }
        }
        printf("Introduce el valor del vector independiente en [%d]: ", i + 1);
        if (scanf("%lf", &vector[i]) != 1) {
            printf("Entrada invalida.\n");
            exit(1);
        }
    }
}

void mostrarmatriz(int n, double matriz[][n], double vector[]) {
    printf("\nMatriz:\n");
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%8.3lf ", matriz[i][j]);
        }
        printf("| %8.3lf\n", vector[i]);
    }
}

void corregircoeficiente(int n, double matriz[][n], double vector[]) {
    int fila, columna;
    printf("Fila a corregir [i]: ");
    scanf("%d", &fila);
    printf("Columna a corregir [j] (o 0 para el vector independiente): ");
    scanf("%d", &columna);

    if (columna == 0) {
        printf("Introduce el nuevo valor del vector independiente en [%d]: ", fila);
        scanf("%lf", &vector[fila - 1]);
    } else {
        printf("Ingresar el nuevo valor para [%d][%d]: ", fila, columna);
        scanf("%lf", &matriz[fila - 1][columna - 1]);
    }
}

int DD(int n, double matriz[][n]) {
    int i, j;
    for (i = 0; i < n; i++) {
        double suma = 0;
        for (j = 0; j < n; j++) {
            if (i != j) {
                suma += fabs(matriz[i][j]);
            }
        }
        if (fabs(matriz[i][i]) < suma) {
            return 0;
        }
    }
    return 1;
}

double determinante(int n, double matriz[][n]) {
    double copiamatriz[n][n];
    int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            copiamatriz[i][j] = matriz[i][j];
        }
    }

    double det = 1;
    for (i = 0; i < n; i++) {
        for (k = i + 1; k < n; k++) {
            if (copiamatriz[i][i] == 0) return 0;
            double factor = copiamatriz[k][i] / copiamatriz[i][i];
            for (j = i; j < n; j++) {
                copiamatriz[k][j] -= factor * copiamatriz[i][j];
            }
        }
        det *= copiamatriz[i][i];
    }
    return det;
}

void resolver(int n, double matriz[][n], double vector[]) {
    double soluciones[n];
    int i, j;

    for (i = n - 1; i >= 0; i--) {
        soluciones[i] = vector[i];
        for (j = i + 1; j < n; j++) {
            soluciones[i] -= matriz[i][j] * soluciones[j];
        }
        soluciones[i] /= matriz[i][i];
    }

    printf("\nSolución del sistema:\n");
    for (i = 0; i < n; i++) {
        printf("x[%d] = %lf\n", i + 1, soluciones[i]);
    }
}

///////////////////////////PROGRAMA 3//////////////////

void jacobi(int n, double matriz[n][n], double vector[n], double x_inicial[n], int max_iter, double tol) {
    double x[n], x_new[n], error;
    int i,j,k;
    for ( i = 0; i < n; i++) {
        x[i] = x_inicial[i];
    }

    printf("\nIteracion\tVector solucion\t\t\t\tError\n");

    for (k = 1; k <= max_iter; k++) {
        for (i = 0; i < n; i++) {
            double suma = 0;
            for (j = 0; j < n; j++) {
                if (i != j) {
                    suma += matriz[i][j] * x[j];
                }
            }
            x_new[i] = (vector[i] - suma) / matriz[i][i];
        }

        error = 0;
        for ( i = 0; i < n; i++) {
            error = fmax(error, fabs(x_new[i] - x[i]));
            x[i] = x_new[i];
        }

        printf("%d\t\t", k);
        for ( i = 0; i < n; i++) {
            printf("%.6lf ", x[i]);
        }
        printf("\t%.6lf\n", error);

        if (error < tol) {
            printf("\nConvergencia alcanzada en %d iteraciones.\n", k);
            printf("Solucion aproximada:\n");
            for ( i = 0; i < n; i++) {
                printf("x[%d] = %.6lf\n", i + 1, x[i]);
            }
            return;
        }
    }

    printf("\nNo se alcanzo la convergencia en el numero maximo de iteraciones\n");
}


///////////////////////////PROGRAMA POTENCIAS//////////////////

#define MAX_DIM 100

// Función para capturar una matriz cuadrada
void capturarMatrizUnica(int dimension, double matrizUnica[MAX_DIM][MAX_DIM]) {
    printf("Introduce los elementos de la matriz cuadrada:\n");
    for (int fila = 0; fila < dimension; fila++) {
        for (int columna = 0; columna < dimension; columna++) {
            printf("Elemento [%d, %d]: ", fila + 1, columna + 1);
            if (scanf("%lf", &matrizUnica[fila][columna]) != 1) {
                printf("Entrada inválida. Finalizando.\n");
                exit(1);
            }
        }
    }
}

// Función para mostrar una matriz cuadrada
void mostrarMatrizUnica(int dimension, double matrizUnica[MAX_DIM][MAX_DIM]) {
    printf("\nMatriz capturada:\n");
    for (int fila = 0; fila < dimension; fila++) {
        for (int columna = 0; columna < dimension; columna++) {
            printf("%.2lf\t", matrizUnica[fila][columna]);
        }
        printf("\n");
    }
}

// Función para corregir un elemento de la matriz
void corregirElementoMatriz(int dimension, double matrizUnica[MAX_DIM][MAX_DIM]) {
    int filaCorregir, columnaCorregir;
    double nuevoValor;
    printf("Introduce la fila del coeficiente incorrecto (1-%d): ", dimension);
    if (scanf("%d", &filaCorregir) != 1 || filaCorregir < 1 || filaCorregir > dimension) {
        printf("Entrada inválida para la fila. Finalizando.\n");
        exit(1);
    }
    printf("Introduce la columna del coeficiente incorrecto (1-%d): ", dimension);
    if (scanf("%d", &columnaCorregir) != 1 || columnaCorregir < 1 || columnaCorregir > dimension) {
        printf("Entrada inválida para la columna. Finalizando.\n");
        exit(1);
    }
    printf("Introduce el nuevo valor para el coeficiente [%d, %d]: ", filaCorregir, columnaCorregir);
    if (scanf("%lf", &nuevoValor) != 1) {
        printf("Entrada inválida. Finalizando.\n");
        exit(1);
    }
    matrizUnica[filaCorregir - 1][columnaCorregir - 1] = nuevoValor;
    printf("Coeficiente actualizado correctamente.\n");
}

// Función para normalizar un vector
void normalizarVector(double vector[MAX_DIM], int n) {
    double norma = 0.0;
    for (int i = 0; i < n; i++) {
        norma += vector[i] * vector[i];
    }
    norma = sqrt(norma);
    if (fabs(norma - 1.0) > 1e-6) {
        printf("\tEl vector inicial no tiene norma 1. Se normalizará automáticamente.\n");
        for (int i = 0; i < n; i++) {
            vector[i] /= norma;
        }
    }
}

// Método de potencias para calcular valores propios
void metodoPotencias(double matriz[MAX_DIM][MAX_DIM], double vector[MAX_DIM], int n, double tol, int max_iter) {
    double lambda_anterior = 0.0, lambda, nuevo_vector[MAX_DIM];
    printf("\n\tIteración\tValor Propio\tVector Propio\n");
    printf("\t---------------------------------------------------\n");

    for (int iter = 1; iter <= max_iter; iter++) {
        for (int i = 0; i < n; i++) {
            nuevo_vector[i] = 0.0;
            for (int j = 0; j < n; j++) {
                nuevo_vector[i] += matriz[i][j] * vector[j];
            }
        }

        double norma = 0.0;
        for (int i = 0; i < n; i++) {
            norma += nuevo_vector[i] * nuevo_vector[i];
        }
        norma = sqrt(norma);
        if (fabs(norma) < 1e-10) {
            printf("\tError: La norma del vector es 0. No se puede normalizar.\n");
            return;
        }
        for (int i = 0; i < n; i++) {
            nuevo_vector[i] /= norma;
        }

        lambda = 0.0;
        for (int i = 0; i < n; i++) {
            lambda += nuevo_vector[i] * vector[i];
        }

        printf("\t%d\t\t%.6lf\t[", iter, lambda);
        for (int i = 0; i < n; i++) {
            printf(" %.6lf", nuevo_vector[i]);
            if (i < n - 1) printf(",");
        }
        printf(" ]\n");

        if (fabs(lambda - lambda_anterior) < tol) {
            printf("\n\tConvergencia alcanzada en iteración: %d\n", iter);
            return;
        }

        lambda_anterior = lambda;
        for (int i = 0; i < n; i++) {
            vector[i] = nuevo_vector[i];
        }
    }
    printf("\n\tNo se alcanzó convergencia en el número máximo de iteraciones.\n");
}



int main(){
	
	int op;
    int programa,n, corregir, capturarOtra, n3, corregir3, n1, max_iter, fila, columna, i, j,capturarOtraMatriz;
    double matriz[MAX][MAX], vector[MAX], tol;
	do{

   	mostrarMenu();
    printf("\n\tPaquete de programas\n");
    printf("1-Solucion de ecuaciones\n");
    printf("2-Solucion de sistemas de ecuaciones\n");
    printf("3-Metodo de potencias\n");
    printf("4-Salir\n");
	printf("\nQue programa desea abrir?\n");
    scanf("%d", &programa);
    system("pause");
	system("cls");
    
    
    switch(programa){
    	case 1:
	    	do{

		        printf("\n\tMenu de metodos\n");
		        printf("1 - Metodo de Biseccion\n");
		        printf("2 - Metodo de la Secante\n");
		        printf("Que metodo desea realizar?: ");
		        scanf("%d", &op);
		        system("pause");
		        system("cls");
		
		        if (op == 1) {
		            int funcion, iteraciones;
		            double a, b, tolerancia;
		
		            printf("\n\tMenu de Funciones\n");
		            printf("1 - f(x) = x^2 * cos(x) - 2x\n");
		            printf("2 - f(x) = (6 - 2/x^2) * (exp(2+x)/4) + 1\n");
		            printf("3 - f(x) = x^3 - 3 * sin(x^2) + 1\n");
		            printf("4 - f(x) = x^3 + 6x^2 + 9.4x + 2.5\n");
		            printf("Que funcion desea calcular?\n");
		            scanf("%d", &funcion);
		
		            system("pause");
		            system("cls");
		
		            // Parámetros para el método de bisección
		            printf("Ingrese los valores del intervalo:\n");
		            printf("a: ");
		            scanf("%lf", &a);
		
		            printf("b: ");
		            scanf("%lf", &b);
		
		            printf("Maximo de iteraciones: ");
		            scanf("%d", &iteraciones);
		
		            printf("Tolerancia: ");
		            scanf("%lf", &tolerancia);
		
		            system("pause");
		            system("cls");
		
		            metodoDeBiseccion(f, a, b, iteraciones, tolerancia);
		        } else if (op == 2) {
		            double a, b, epsilon;
		            int imax, funcIndex;
		
		            printf("Seleccione la funcion a evaluar:\n");
		            printf("1. f(x) = x^2 * cos(x) - 2x\n");
		            printf("2. f(x) = (6 - 2 / (x * x)) * (exp(2 + x) / 4) + 1\n");
		            printf("3. f(x) = x^3 - 3 * sin(x^2) + 1\n");
		            printf("4. f(x) = x^3 + 6 * x^2 + 9.4 * x + 2.5\n");
		            printf("Que funcion desea calcular?\n");
		            scanf("%d", &funcIndex);
		
		            printf("Ingrese los valores del intervalo:\n");
		            printf("a: ");
		            scanf("%lf", &a);
		
		            printf("b: ");
		            scanf("%lf", &b);
		            
		            
		            printf("Introduzca el criterio de convergencia: ");
		            scanf("%lf", &epsilon);
		            printf("Introduzca el número maximo de iteraciones permitidas: ");
		            scanf("%d", &imax);
		
		            secante(a, b, epsilon, imax, funcIndex);
		        } else {
		            printf("Opcion incorrecta\n");
		        }
		
		        printf("\nDesea:\n1. Obtener otra raiz\n0 Salir\n");
		        scanf("%d", &op);
		
		        system("pause");
		        system("cls");
		
			    }while(op == 1);
			    printf("\nGracias por utilizar el programa!\n");
			    system("pause");
				system("cls");
	    		
				break;//termina programa 1
				
///////////////////////////PROGRAMA 2//////////////////
		case 2:
			
    		do{
	         printf("Dimension de la matriz cuadrada: ");
			    if (scanf("%d", &n) != 1 || n <= 0) {
			        printf("Entrada invalida para la dimension. Finalizando.\n");
			        return 1;
			    }
			
			    double matriz[n][n], vector[n], x_inicial[n];
			    leermatriz(n, matriz, vector);
			
			    do {
			        mostrarmatriz(n, matriz, vector);
			        printf("\nEs correcta la matriz? ([SI=1], [NO=0]): ");
			        if (scanf("%d", &corregir) != 1 || (corregir != 0 && corregir != 1)) {
			            printf("Entrada invalida. Finalizando.\n");
			            return 1;
			        }
			        if (!corregir) {
			            corregircoeficiente(n, matriz, vector);
			        }
			    } while (!corregir);
			
			    int dominante = DD(n, matriz);
			    double det = determinante(n, matriz);
			
			    if (det != 0) {
			        printf("\nEl determinante es: %lf\n", det);
			
			        if (!dominante) {
			            printf("\nLa matriz NO es EDD. La convergencia no se garantiza.\n");
			        }
			
			        printf("Introduce el vector inicial:\n");
			        int i;
			        for (i = 0; i < n; i++) {
			            printf("\nx_inicial[%d]: ", i + 1);
			            if (scanf("%lf", &x_inicial[i]) != 1) {
			                printf("Entrada invalida. Finalizando.\n");
			                return 1;
			            }
			        }
			
			        int max_iter;
			        double tol;
			        printf("Maximo de iteraciones: ");
			        if (scanf("%d", &max_iter) != 1 || max_iter <= 0) {
			            printf("Entrada invalida para el numero de iteraciones. Finalizando.\n");
			            return 1;
			        }
			        printf("Tolerancia: ");
			        if (scanf("%lf", &tol) != 1 || tol <= 0) {
			            printf("Entrada invalida para la tolerancia. Finalizando.\n");
			            return 1;
			        }
			
			        jacobi(n, matriz, vector, x_inicial, max_iter, tol);
			    } else {
			        printf("El determinante es 0. El sistema no tiene solucion unica.\n");
			    }
			    
			    printf("¿Desea capturar otra matriz? [SI-1], [NO-0]: ");
		        scanf("%d", &capturarOtra);
		        system("pause");
				system("cls");
		
			
						
			}while(capturarOtra);
						
				printf("\nGracias por utilizar el programa!\n");
				system("pause");
				system("cls");

	    	break;

///////////////////////////PROGRAMA 3//////////////////	    
	    case 3:

    do {
        int dimensionMatriz, corregirMatriz;
        printf("Dimensión de la matriz cuadrada: ");
        if (scanf("%d", &dimensionMatriz) != 1 || dimensionMatriz <= 0 || dimensionMatriz > MAX_DIM) {
            printf("Entrada inválida para la dimensión. Finalizando.\n");
            return 1;
        }

        double matrizProcesada[MAX_DIM][MAX_DIM], vectorProcesado[MAX_DIM];

        capturarMatrizUnica(dimensionMatriz, matrizProcesada);

        do {
            mostrarMatrizUnica(dimensionMatriz, matrizProcesada);
            printf("\n¿Es correcta la matriz? ([SI=1], [NO=0]): ");
            if (scanf("%d", &corregirMatriz) != 1 || (corregirMatriz != 0 && corregirMatriz != 1)) {
                printf("Entrada inválida. Finalizando.\n");
                return 1;
            }
            if (!corregirMatriz) {
                corregirElementoMatriz(dimensionMatriz, matrizProcesada);
            }
        } while (!corregirMatriz);

        printf("Introduce el vector inicial:\n");
        for (int i = 0; i < dimensionMatriz; i++) {
            printf("vectorInicial[%d]: ", i + 1);
            if (scanf("%lf", &vectorProcesado[i]) != 1) {
                printf("Entrada inválida. Finalizando.\n");
                return 1;
            }
        }

        normalizarVector(vectorProcesado, dimensionMatriz);

        int maximoIteraciones;
        double toleranciaEspecificada;
        printf("Máximo de iteraciones: ");
        if (scanf("%d", &maximoIteraciones) != 1 || maximoIteraciones <= 0) {
            printf("Entrada inválida. Finalizando.\n");
            return 1;
        }
        printf("Tolerancia: ");
        if (scanf("%lf", &toleranciaEspecificada) != 1 || toleranciaEspecificada <= 0) {
            printf("Entrada inválida. Finalizando.\n");
            return 1;
        }

        metodoPotencias(matrizProcesada, vectorProcesado, dimensionMatriz, toleranciaEspecificada, maximoIteraciones);

        printf("\n¿Desea capturar otra matriz? [SI=1], [NO=0]: ");
        if (scanf("%d", &capturarOtraMatriz) != 1) {
            printf("Entrada inválida. Finalizando.\n");
            return 1;
        }
    } while (capturarOtraMatriz);

		    
	    			system("pause");
					system("cls");
			    	break;
		
		
		case 4:
			printf("\nGracias por utilizar el programa!\n"); 
			
			break;
			    
			    default: 
			    	printf("Opcion incorrecta");
			    	system("pause");
					system("cls");
			    	break;
			    	
		
				
			    
			    
			}//llave switch
		
		
	}while (programa != 0);
	
		printf("\nGracias por utilizar el programa!\n"); 

return 0;

}






    


   
