#include <string>
#include <vector>
#include "dlx.h"

using namespace std;

int main() {

	// Constructor 1
	vector<string> column_names;
	column_names = {"A", "B", "C", "D", "E", "F", "G"};
	vector<vector<int> > test_matrix{ {0, 0, 1, 0, 1, 1, 0},
					  {1, 0, 0, 1, 0, 0, 1},
					  {0, 1, 1, 0, 0, 1, 0},
					  {1, 0, 0, 1, 0, 0, 0},
					  {0, 1, 0, 0, 0, 0, 1},
					  {0, 0, 0, 1, 1, 0, 1}};
	
	DancingLinks dlx1 = DancingLinks(test_matrix, column_names);
	dlx1.search();

	// Constructor 2
	std::string data = "A,B,C,D,E,F,G\n"
		           "0,0,1,0,1,1,0\n"
		           "1,0,0,1,0,0,1\n"
		           "0,1,1,0,0,1,0\n"
		           "1,0,0,1,0,0,0\n"
		           "0,1,0,0,0,0,1\n"
		           "0,0,0,1,1,0,1";

	DancingLinks dlx2 = DancingLinks(data);
	dlx2.search();
		
	return 0;
}
