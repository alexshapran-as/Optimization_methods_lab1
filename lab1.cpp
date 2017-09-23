#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>

double * simplex_method(std::vector<double> c, std::vector<double> b, double ** a);
bool not_reference(double ** arr, unsigned int m);
bool not_optimal(double ** arr, unsigned int m, unsigned int n);
std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n);
unsigned int ri = {0};
unsigned int rj = {0};

int main()
{
	using std::cout;
	using std::cin;
	using std::endl;
	using std::vector;
	system("chcp 65001");

	vector<double> c;
	double temp = {0.0};
	vector<double> b;
	double ** a;
	unsigned int n = {0};
	unsigned int m = {0};

	cout << endl << "Заполните вектор коэффициентов c = [c1,c2,...,cn] целевой функции" 
	" (ЦФ) F (q - для прекращения ввода):" << endl << "c[" << 1 << "] = ";
	cin >> temp;
	c.push_back(temp);
	cout << "c[" << 2 << "] = ";
	unsigned int i = {1};
	for (; cin >> temp; ++i)
	{
		c.push_back(temp);
		if (i == 1)
			cout << "c[" << i + 2 << "] = ";
		else
			cout << "c[" << i + 1 << "] = ";
	}
	cin.clear();
	while(cin.get() != '\n');
	n = i;

	cout << endl << "Заполните вектор правой части системы ограничений" 
	" b = [b1,b2,...,bm] (q - для прекращения ввода):" << endl << "b[" << 1 << "] = ";
	cin >> temp;
	b.push_back(temp);
	cout << "b[" << 2 << "] = ";
	unsigned int j = {1};
	for (; cin >> temp; ++j)
	{
		b.push_back(temp);
		if (j == 1)
			cout << "b[" << j + 2 << "] = ";
		else
			cout << "b[" << j + 1 << "] = ";
	}
	cin.clear();
	while(cin.get() != '\n');
	m = j;

	if (m > 2*n)
	{
		cout << endl << "Число линейно независимых уравнений больше числа переменных."
		" Такая СЛАУ несовместна.";
		exit(EXIT_FAILURE);
	}

	a = new double * [m];
	for (unsigned int i = {0}; i < m; ++i)
		a[i] = new double [n];
	cout << endl << "Заполните матрицу системы ограничений"
	" A = [[a11,...,a1n]...[am1,...,amn]] (q - для прекращения ввода):" << endl;
	for (unsigned int i = {0}; i < m; ++i)
	{
		for (unsigned int j = {0}; j < n; ++j)
		{
			cout << "a" << i+1 << j+1 << " = ";
			if (!(cin >> a[i][j]))
			{
				cout << endl << "Ввод не удался!";
				exit(EXIT_FAILURE);
			}
		}
	}

	simplex_method(c,b,a);

	for (unsigned int i = 0; i < m; ++i) 
		delete a[i]; 
	delete [] a;	
	return 0;
}

double * simplex_method(std::vector<double> c, std::vector<double> b, double ** a)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;
	vector<string> x(b.size() + c.size() + 1);

	double ** arr = new double * [b.size()+1];
	for (unsigned int i = {0}; i < b.size()+1; ++i)
		arr[i] = new double [c.size()+1];
	for (unsigned int i = {0}; i < b.size(); ++i)
		arr[i][0] = b[i];
	arr[b.size()][0] = {0};
	for (unsigned int i = {0}; i < b.size(); ++i)
	{
		unsigned int k = {1};
		for (unsigned int j = {0}; j < c.size(); ++j)
		{
			arr[i][k] = a[i][j];
			k++;
		}
	} 
	unsigned int k = {1};
	for (unsigned  int i = {0}; i < c.size(); ++i)
	{
		arr[b.size()][k] = -c[i];
		k++;
	}

	print(arr, b.size(), c.size());

	while (not_reference(arr,b.size()))
	{
		cout << endl << "Недопустимое решение. Производим замену базиса," 
		" чтобы получить опорное решение.";
		unsigned int i = {0};
		while (arr[i][0] > 0 && i < b.size() + 1)
			i++;
		unsigned int j = {1}; // Разрешающий столбец
		unsigned int k = {0};
		for (; j < c.size() + 1; ++j)
		{
			if (arr[i][j] < 0)
			{
				k++;
				break;
			}
		}
		if (k == 0)
		{
			cout << endl << "Задача не имеет допустимых решений!";
			exit(EXIT_FAILURE);
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			temp.push_back(arr[i][0] / arr[i][j]);
		}
		std::sort(temp.begin(),temp.end());
		unsigned int m = {0};
		while (temp[m] < 0)
			m++;
		i = {0}; // Разрешающая строка
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j)
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j)
				{
					arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j)
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j];
		ri = i;
		rj = j;
		print(arr, b.size(), c.size());
	}

	cout << endl << "Найдено опорное решение. Проводим проверку на оптимальность.";
	while (not_optimal(arr,b.size(),c.size()))
	{
		cout << endl << "Неоптимальное решение. Производим замену базиса.";
		unsigned int j = {0}; // Новый разрешающий столбец
		for (; j < c.size() + 1; ++j)
		{
			if (arr[b.size()][j] > 0)
				break;
		}
		vector<double> temp;
		for (unsigned int i = {0}; i < b.size(); ++i)
		{
			temp.push_back(arr[i][0] / arr[i][j]);
		}
		std::sort(temp.begin(),temp.end());
		unsigned int m = {0};
		while (temp[m] < 0)
			m++;
		unsigned int i = {0}; // Новая разрешающая строка
		for (; temp[m] != arr[i][0] / arr[i][j]; ++i);
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n != j)
				{
					arr[m][n] = arr[m][n] - (arr[m][j] * arr[i][n]) / (arr[i][j]);
				}
			}
		}
		for (unsigned int m = {0}; m < b.size() + 1; ++m)
		{
			for (unsigned int n = {0}; n < c.size() + 1; ++n)
			{
				if (m != i && n == j)
				{
					arr[m][n] = -arr[m][n] / arr[i][j];
				}
				if (m == i && n != j)
				{
					arr[m][n] = arr[m][n] / arr[i][j];
				}
			}
		}
		arr[i][j] = 1.0 / arr[i][j];
		ri = i;
		rj = j;
		x = print(arr, b.size(), c.size());
	}

	cout << endl << "Найдено оптимальное решение!" << endl;
	for (unsigned int i = {0}; i < c.size(); ++i)
		cout << x[i] << " = 0; ";
	unsigned int l = {0};
	for (unsigned int j = {c.size()}; j < b.size()+c.size(); ++j)
	{
		cout << endl << x[j] << " = " << arr[l++][0] << ";";
	}
}

bool not_reference(double ** arr, unsigned int m)
{
	unsigned int accept = {0};
	for (unsigned int i = {0}; i < m; ++i)
	{
		if (arr[i][0] < 0)
			accept++;
	}
	if (accept > 0)
		return true;
	else
		return false;
}

bool not_optimal(double ** arr, unsigned int m, unsigned int n)
{
	unsigned int check = {0};
	for (unsigned int j = {1}; j < n + 1; ++j)
	{
		if (arr[m][j] > 0)
			check++;
	}
	if (check > 0)
		return true;
	else 
		return false;
}

std::vector<std::string> print(double ** arr, unsigned int m, unsigned int n)
{
	using std::cout;
	using std::endl;
	using std::vector;
	using std::string;

	static unsigned int count = {0};
	cout << endl << "Simplex table:" << endl << "Si0 ";
	static vector<string> x(m+n+1);
	if (count == 0)
	{
		for (unsigned int i = {0}; i < m+n; ++i)
		{
			char temp[24]; 
			itoa(i+1,temp,10);
			x[i] = "X";
			x[i] += temp;
		}
		x[m+n] = "F";
	}
	if (count > 0)
	{
		for (unsigned int i = {0}; i < n; ++i)
		{
			if (i == rj - 1)
			{
				for (unsigned int j = {0}; j < m+n; ++j)
				{
					if (j == ri)
					{
						string temp = x[i];
						x[i] = x[j+n];
						x[j+n] = temp; 
					}
				}
			}
		}
	}
	count++;
	unsigned int k = {0};
	for (; k < n; ++k)
		cout << x[k] << " ";
	cout << endl;
	for (unsigned int i = {0}; i < m+1; ++i)
	{
		for (unsigned int j = {0}; j < n+1; ++j)
		{
			if (arr[i][j] > 0)
				cout << arr[i][j] << "   ";
			else
				cout << arr[i][j] << "  ";
		}
		cout << x[k++] << endl;
	}
	return x;
}