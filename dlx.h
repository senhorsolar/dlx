// Dancing Links solver for the exact cover problem.

// Finds subset of rows in a binary matrix such that each row is
// disjoint, and whose union covers each column.


#ifndef DLX_H
#define DLX_H

#include <algorithm>
#include <iostream>
#include <vector>

namespace Dlx {

using Index=size_t;
using BinaryMatrix=std::vector<std::vector<bool>>;

struct Node
{
    Index left;
    Index right;
    Index up;
    Index down;

    Index row;
    Index column;
    
    Index size; // column size

    // ctor for columns and root
    Node(Index self_idx) :
	left(self_idx), right(self_idx), up(self_idx), down(self_idx),
	column(self_idx), row(-1), size(0)
    {}

    // ctor for every other node
    Node(Index self_idx, Index row_idx, Index col_idx) :
	left(self_idx), right(self_idx), up(self_idx), down(self_idx),
	column(col_idx), row(row_idx), size(0)
    {}
};

class DancingLinks
{
public:
    
    using Row=size_t;
    using Solution=std::vector<Row>;
    
public:

    // ctor
    DancingLinks(const BinaryMatrix& bin_mat);

    // Find all solutions and store in memory
    void Search();

    // Get all solutions
    const std::vector<Solution>& GetSolutions() const {
	return m_solutions;
    }

    // Return number of solutions
    const size_t SolutionCount() const { return m_solutions.size(); }

    const std::vector<size_t>& RowToColumns(size_t r) const { return m_rows.at(r);}
    const Solution& GetSolution(size_t idx) const { return m_solutions.at(idx);}

    // Helper function for converting a solution (subset of rows) to
    // whatever the user defines the covered columns to be
    template<class T>
    std::vector<std::vector<T>> DecodeSolution(const std::vector<T>& columns, const Solution&) const;

    
private:
    std::vector<Node> m_nodes; // quad linked list
    std::vector<Solution> m_solutions; // soltuions
    std::vector<size_t> m_cur_solution; // for building current solution
    std::vector<std::vector<size_t>> m_rows; // map from row to covered columns

    // dlx specific functions
    void Cover(Index col_idx); // cover column and remove conflicting rows
    void Uncover(Index col_idx); // undo effects of Cover()
    Index ChooseColumn(); // heuristic for choosing column
    void Search_(size_t k); // recursive search function
};

template<class T>
std::vector<std::vector<T>>
DancingLinks::DecodeSolution(const std::vector<T>& columns, const Solution& solution) const
{
    std::vector<std::vector<T>> decoded_solution;
    for (size_t row : solution) {
	std::vector<T> decoded_row;
	const auto& cols = RowToColumns(row);
	for (size_t col : cols) {
	    decoded_row.push_back(columns[col]);
	}
	decoded_solution.push_back(std::move(decoded_row));
    }
    return decoded_solution;
}
    
} // namespace

#endif // DLX_H
