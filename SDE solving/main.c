//
//  main.c
//  SDE solving
//
//  Created by Лабутин Антон Александрович on 05.05.2021.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ITER_CNT = 0; // количество итераций
int EQUATIONS_CNT; // количество уравнений
int TESTS_CNT = 6;

double step; // шаг сетки
double x_0;
double *y_0;

double
(*function[2]) (double, double *); // массив указателей на функцию

double
(*solution) (double);


double
f1_1(double x, double *y) // уравнение из таблицы 1 - 1
{
    return (3 - y[0] - x);
}


void
init_f1_1(void) // начальные условия для уравнения из таблицы 1 - 1
{
    EQUATIONS_CNT = 1;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 0.0;
    y_0[0] = 0.0;
}


double
sol_f1_1(double x)
{
    return (4 - x - 4 * exp(-x));
}


double
f1_2(double x, double *y) // уравнение из таблицы 1 - 2
{
    return (sin(x) - y[0]);
}


void
init_f1_2(void) // начальные условия для уравнения из таблицы 1 - 2
{
    EQUATIONS_CNT = 1;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 0.0;
    y_0[0] = 10.0;
}


double
sol_f1_2(double x)
{
    return (0.5 * (-cos(x) + sin(x) + 21 * exp(-x)));
}


double
f1_4(double x, double *y) // уравнение из таблицы 1 - 4
{
    return (y[0] - y[0] * x);
}


void
init_f1_4(void) // начальные условия для уравнения из таблицы 1 - 4
{
    EQUATIONS_CNT = 1;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 0.0;
    y_0[0] = 5.0;
}


double
sol_f1_4(double x)
{
    return (5 * exp(-0.5 * x * (-2 + x)));
}


double
f2_1_1(double x, double *y) // первое уравнение системы из таблицы 2 - 1
{
    return ((y[0] - y[1]) / x);
}


double
f2_1_2(double x, double *y) // второе уравнение системы из таблицы 2 - 1
{
    return ((y[0] + y[1]) / x);
}


void
init_f2_1(void) // начальные условия для уравнения из таблицы 2 - 1
{
    EQUATIONS_CNT = 2;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 1.0;
    y_0[0] = 1.0;
    y_0[1] = 1.0;
}


double
f2_2_1(double x, double *y) // первое уравнение системы из таблицы 2 - 2
{
    return (x * y[0] + y[1]);
}


double
f2_2_2(double x, double *y) // второе уравнение системы из таблицы 2 - 2
{
    return (y[0] - y[1]);
}


void
init_f2_2(void) // начальные условия для уравнения из таблицы 2 - 2
{
    EQUATIONS_CNT = 2;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 0.0;
    y_0[0] = 0.0;
    y_0[1] = 1.0;
}


double
f2_14_1(double x, double *y) // первое уравнение системы из таблицы 2 - 14
{
    return (exp(-( pow(y[0], 2) + pow(y[1], 2) )) + 2 * x);
}


double
f2_14_2(double x, double *y) // второе уравнение системы из таблицы 2 - 14
{
    return (2 * pow(y[0], 2) + y[1]);
}


void
init_f2_14(void) // начальные условия для уравнения из таблицы 2 - 14
{
    EQUATIONS_CNT = 2;
    y_0 = malloc(EQUATIONS_CNT * sizeof(double));
    x_0 = 0.0;
    y_0[0] = 0.5;
    y_0[1] = 1.0;
}


void
init_test(int test_num) // выбор номера теста и соответствующих начальных условий
{
    switch (test_num) {
        case 1:
            init_f1_1();
            function[0] = f1_1;
            solution = sol_f1_1;
            break;
        case 2:
            init_f1_2();
            function[0] = f1_2;
            solution = sol_f1_2;
            break;
        case 3:
            init_f1_4();
            function[0] = f1_4;
            solution = sol_f1_4;
            break;
        case 4:
            init_f2_1();
            function[0] = f2_1_1;
            function[1] = f2_1_2;
            solution = NULL;
            break;
        case 5:
            init_f2_2();
            function[0] = f2_2_1;
            function[1] = f2_2_2;
            solution = NULL;
            break;
        case 6:
            init_f2_14();
            function[0] = f2_14_1;
            function[1] = f2_14_2;
            solution = NULL;
            break;
        default:
            break;
    }
}


void
print_vector(double x, double *y) // печать вектора (x, *y)
{
    printf("(%.10g", x);
    for (int i = 0; i < EQUATIONS_CNT; ++i) {
        printf("; %.10g", y[i]);
    }
    printf(")\n");
}


void
RungeKutta_2(double x_0, double *y_0) // метод Рунге−Кутта второго порядка точности
{
    double *k = malloc(EQUATIONS_CNT * sizeof(double));
    double *y = malloc(EQUATIONS_CNT * sizeof(double));
    print_vector(x_0, y_0);

    for (int i = 0; i < ITER_CNT; ++i) {
        for (int j = 0; j < EQUATIONS_CNT; ++j) {
            k[j] = y_0[j] + function[j](x_0, y_0) * step;
        }

        for (int j = 0; j < EQUATIONS_CNT; ++j) {
            y[j] = y_0[j] + (function[j](x_0, y_0) + function[j](x_0 + step, k)) * step / 2;
        }

        x_0 += step;

        print_vector(x_0, y);

        for (int j = 0; j < EQUATIONS_CNT; ++j) {
            y_0[j] = y[j];
        }
    }

    free(y);
    free(k);
}


void
RungeKutta_4 (double x_0, double *y_0) // метод Рунге−Кутта четвертого порядка точности
{
    double *y =  malloc(EQUATIONS_CNT * sizeof(double));
    double *k1 = malloc(EQUATIONS_CNT * sizeof(double));
    double *k2 = malloc(EQUATIONS_CNT * sizeof(double));
    double *k3 = malloc(EQUATIONS_CNT * sizeof(double));
    double *k4 = malloc(EQUATIONS_CNT * sizeof(double));

    print_vector(x_0, y_0);

    for(int i = 0; i < ITER_CNT; ++i){
        for(int j = 0; j < EQUATIONS_CNT; ++j){
            k1[j] = function[j](x_0, y_0);
            y[j] = y_0[j] + step * k1[j] / 2;
        }

        for(int j = 0; j < EQUATIONS_CNT; ++j) {
            k2[j] = function[j](x_0 + step / 2, y);
            y[j] = y_0[j] + step * k2[j] / 2;
        }

        for(int j = 0; j < EQUATIONS_CNT; ++j) {
            k3[j] = function[j](x_0 + step / 2, y);
            y[j] = y_0[j] + step * k3[j];
        }

        for(int j = 0; j < EQUATIONS_CNT; ++j) {
            k4[j] = function[j](x_0 + step, y);
            y[j] = y_0[j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) * step/6;
        }

        x_0 += step;

        print_vector(x_0, y);

        for(int j = 0; j < EQUATIONS_CNT; ++j) {
            y_0[j] = y[j];
        }
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(y);
}


void
exact_solution(double x_0, double *y_0)
{
    double *y = malloc(EQUATIONS_CNT * sizeof(double));
    print_vector(x_0, y_0);

    double x = x_0;
    for (int i = 0; i < ITER_CNT; ++i) {
        x += step;
        y[0] = solution(x);
        print_vector(x, y);
    }

    free(y);
}


int
main(void)
{
    int test_num;
    do {
        printf("Input test number:\n");
        printf("1: example 1-1\n");
        printf("2: example 1-2\n");
        printf("3: example 1-4\n");
        printf("4: example 2-1\n");
        printf("5: example 2-2\n");
        printf("6: example 2-14\n");
        scanf("%d", &test_num);
    } while (test_num < 1 && test_num > TESTS_CNT);

    printf("Input count of iterations:\n");
    scanf("%d", &ITER_CNT);

    printf("Input step:\n");
    scanf("%lf", &step);

    printf("\n");
    init_test(test_num);
    printf("Runge Kutta method of the second order of accuracy:\n");
    RungeKutta_2(x_0, y_0);
    printf("\n");

    init_test(test_num);
    printf("Runge Kutta method of the fourth order of accuracy:\n");
    RungeKutta_4(x_0, y_0);
    printf("\n");

    if (solution != NULL) {
        init_test(test_num);
        printf("Exact solution:\n");
        exact_solution(x_0, y_0);
        printf("\n");
    }

    return 0;
}
