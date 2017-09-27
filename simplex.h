#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#include <vector>
#include <string>
void show_task(std::vector<double> c, std::vector<double> b, double** a);
void set_data(std::vector<double> &c, std::vector<double> &b, double** a, unsigned int &n, unsigned int &m);
double * simplex_method(std::vector<double> c, std::vector<double> b, double ** a);
bool not_reference(double ** arr, unsigned int m);
bool not_optimal(double ** arr, unsigned int m, unsigned int n);
std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n);
static unsigned int ri = {0}; // Разрешающая строка
static unsigned int rj = {0}; // Разрешающий столбец
static std::string min_max; // Ищем min/max
#endif