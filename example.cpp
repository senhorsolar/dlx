#include "dlx.h"
#include <iostream>
#include <random>
#include <string>

template<class T>
using Matrix=std::vector<std::vector<T>>;

template<class T>
Matrix<T> gen_rand_matrix(size_t n)
{
    Matrix<T> mat;
    for (size_t i = 0; i < n / 2; ++i) {
	std::vector<T> row1;
	std::vector<T> row2;
	for (size_t j = 0; j < n; ++j) {
	    float r = (float)rand() / (float)RAND_MAX;
	    row1.push_back((T)((r > 0.5) ? 1 : 0));
	    row2.push_back((T)((row1[j] == 0) ? 1 : 0));
	}
	mat.push_back(std::move(row1));
	mat.push_back(std::move(row2));
    }
    return mat;
}

std::vector<std::string> gen_column_names(size_t n)
{
    std::vector<std::string> col_names;
    for (size_t i = 0; i < n; ++i)
	col_names.push_back(std::to_string(i));
    return col_names;
}

int main()
{
    using namespace std;
    
    vector<string> column_names;
    column_names = {"A", "B", "C", "D", "E", "F", "G"};
    vector<vector<bool> > test_matrix{{0, 0, 1, 0, 1, 1, 0},
				      {1, 0, 0, 1, 0, 0, 1},
				      {0, 1, 1, 0, 0, 1, 0},
				      {1, 0, 0, 1, 0, 0, 0},
				      {0, 1, 0, 0, 0, 0, 1},
				      {0, 0, 0, 1, 1, 0, 1}};

    // size_t n = 1000;
    // test_matrix = gen_rand_matrix<bool>(n);
    // column_names = gen_column_names(n);
    
    auto dlx = Dlx::DancingLinks(test_matrix);
    dlx.Search();

    std::cout << dlx.SolutionCount() << '\n';

    for (const auto& solution : dlx.GetSolutions()) {
	auto decoded_solution = dlx.DecodeSolution(column_names, solution);
	for (auto& row : decoded_solution) {
	    for (auto& s : row) {
		std::cout << s << " ";
	    }
	    std::cout << '\n';
	}
    }

    return 0;
}
	
