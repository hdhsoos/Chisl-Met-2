#include <iostream> 
#include <stdio.h> 
#include <locale.h> 
#include <math.h> 

const double E = 2.718281828;
double A = 0.5;
double B = 1;
double C = -1;
double Eps_1 = 0.001;
double Del_1 = 0.001;
double Eps_2 = 0.00001;
double Del_2 = 0.00001;

double f(double x) { return x * log10(x) - 1; }
double fd_1(double x) { return log10(x) + 1 / (log(10)); }
double fd_2(double x) { return 1 / (log(10) * x); }
void MPI(double x0, double Epsilon, double Delta)
{
	int n = 1;
	double x, y, z;
	x = x0;
	printf(" \t n+1 | \t X_n   |\t X_n+1 |  |X_n+1 - X_n|  | |f(x_n + 1)|  | \n"

		"------------------------------------------------------------------------ | \n");
	do {
		y = x + C * f(x);
		z = x;
		printf("%10d | %12.8lf|%13.8lf|%17.8lf|%15.8lf|\n", n++, x, y, fabs(y - x), fabs(f(y)));
		x = y;
	} while (fabs(z - x) > Epsilon || fabs(f(x)) > Delta);
	printf("\n");
}
void MN(double x0, double Epsilon, double Delta) {
	int n = 1;
	double x, y, z;
	x = x0;
	printf(" \t n+1|   \t X_n      |\t  X_n+1   |   |X_n+1 - X_n|   |   |f(x_n + 1)|  | \n"

		"------------------------------------------------------------------------------- | \n");

	if ((f(x) * fd_2(x) > 0)) {
		do {
			if (fd_1(x) != 0) {
				y = x - (f(x) / fd_1(x));
				z = x;
				printf("%10d | %12.8lf | %13.8lf | %17.8lf | %15.8lf | \n", n++, x, y, fabs(y - x), fabs(f(y)));
				x = y;
			}
			else { printf("Значение производной в данной точке равно нулю\n"); }
		} while (fabs(z - x) > Epsilon || fabs(f(x)) > Delta);
	}
	else { printf("Не выполняется достаточное условие сходимости метода Ньютона\n"); }
	printf("\n");
}
void MMN(double x0, double Epsilon, double Delta) {
	int n = 1;
	double x, y, z = 0;
	x = x0;
	double f0 = fd_1(x0);
	printf(" \t n+1|   \t X_n      |\t  X_n+1   |   |X_n+1 - X_n|    |   |f(x_n + 1)|  | \n"

		"------------------------------------------------------------------------------- | \n");

	if ((f(x) * fd_2(x) > 0)) {
		do {
			if (f0 != 0) {
				y = x - f(x) / f0;
				z = x;
				printf("%10d | %12.8lf | %13.8lf | %17.8lf | %15.8lf | \n", n++, x, y, fabs(y - x), fabs(f(y)));
				x = y;
			}
			else { printf("Значение производной в данной точке равно нулю\n"); }
		} while (fabs(z - x) > Epsilon || fabs(f(x)) > Delta);
	}
	else {
		printf("Не выполняется достаточное условие сходимости модифицированного метода Ньютона\n");
	}printf("\n");
}
int main() {
	setlocale(LC_CTYPE, "Russian");
	double x0;
	do {
		while (true) {
			printf("Введите начальное приближение X0 из отрезка [%2.4lf, %2.4lf] : ", A, B);
			scanf_s("%lf", &x0);
			if ((x0 >= A) & (x0 <= B)) { break; }
			else {
				system("cls");
				printf("Вы вышли за границы допустимых значений\n\n\n");
			}
		}
		printf("\n Начальное приближение X0: %lf Точность: Эпсилон = %lf, Дельта = %lf\n\n", x0, Eps_1, Del_1);
		printf("\n\tМетод простых итераций\n\n");
		MPI(x0, Eps_1, Del_1);
		printf("\n\tМетод Ньютона\n\n");
		MN(x0, Eps_1, Del_1);
		printf("\n\tМодифицированный метод Ньютона\n\n");
		MMN(x0, Eps_1, Del_1);
		printf("\n Начальное приближение X0: %lf Точность: Эпсилон = %lf, Дельта = %lf\n\n", x0, Eps_2, Del_2);
		printf("\n\tМетод простых итераций\n\n");
		MPI(x0, Eps_2, Del_2);
		printf("\n\tМетод Ньютона\n\n");
		MN(x0, Eps_2, Del_2);
		printf("\n\tМодифицированный метод Ньютона\n\n");
		MMN(x0, Eps_2, Del_2);
		int s = getchar();
		int c = getchar();
		if (s == 27)
			break;
		system("cls");
	} while (true);
	return 0;
}
