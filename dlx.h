
// Dancing Links solver for the exact cover problem.

// Finds subset of rows in a binary matrix such that each row is
// disjoint, and whose union covers each column.


#ifndef DLX_H
#define DLX_H

#include <algorithm>
#include <iostream>
#include <vector>

#include "memory_pool.h"

namespace Dlx {

using Index=size_t;
using BinaryMatrix=std::vector<std::vector<bool>>;

class Node
{
public:
    
    Node* left;
    Node* right;
    Node* up;
    Node* down;

    Node* column;
    size_t size; // column size

    Index row_idx; // useful for generating solution
    
    // ctor
    Node() :
	left(this), right(this), up(this), down(this), column(this),
	size(0), row_idx(-1)
	{}

};

class DancingLinks
{
public:
    
    using Row=size_t;
    using Solution=std::vector<Row>;
    
public:

    // ctor
    DancingLinks(const BinaryMatrix& bin_mat)
    {
	if (bin_mat.size() == 0 or bin_mat[0].size() == 0) {
	    // \todo throw error
	    return;
	}

	size_t nrows = bin_mat.size();
	size_t ncols = bin_mat[0].size();

	// Determine how many nodes (nonzero entries in bin_mat)
	// Concurrently, determine which columns each row covers
	size_t node_count = 0;
	m_rows.reserve(nrows);
	
	for (size_t i = 0; i < nrows; ++i) {

	    std::vector<size_t> row;
	    
	    for (size_t j = 0; j < ncols; ++j) {
		if (bin_mat[i][j]) {
		    ++node_count;
		    row.push_back(j);
		}
	    }
	    
	    m_rows.push_back(std::move(row));
	}

	// Initialize memory pool and root node
	m_pool.Resize(node_count + 1 /* for root */ + ncols);
	m_root = m_pool.New();
	m_root->size = (size_t)(-1);

	// For keeping track of last node in column
	std::vector<Node*> prev_row(ncols);

	// Initialize column nodes
	Node* last_node = m_root;
	for (size_t j = 0; j < ncols; ++j) {

	    Node* column_node = m_pool.New();
	    prev_row[j] = column_node;

	    last_node->right = column_node;
	    column_node->left = last_node;

	    last_node = column_node;   
	}
	last_node->right = m_root;
	m_root->left = last_node;

	// Construct a node for each nonzero element in binary matrix
	for (size_t i = 0; i < nrows; ++i) {

	    Node* first_in_row = nullptr;
	    Node* last_in_row = nullptr;
	    
	    for (size_t j : m_rows[i]) {

		Node* cur_node = m_pool.New();
		cur_node->row_idx = i;

		// Up
		Node* up_node = prev_row[j];
		prev_row[j] = cur_node;
		cur_node->up = up_node;
		up_node->down = cur_node;

		// Column
		Node* column_node = up_node->column;
		cur_node->column = column_node;
		++column_node->size;
		
		// Left (if not the first seen in row thus far)
		if (last_in_row == nullptr) {
		    first_in_row = cur_node;
		}
		else {
		    cur_node->left = last_in_row;
		    last_in_row->right = cur_node;
		}

		last_in_row = cur_node;
	    }

	    if (first_in_row != nullptr) {
		first_in_row->left = last_in_row;
		last_in_row->right = first_in_row;
	    }	    
	}

	// Connect last nodes in columns to respective columns
	for (size_t j = 0; j < ncols; ++j) {
	    Node* column_node = prev_row[j]->column;
	    prev_row[j]->down = column_node;
	    column_node->up = prev_row[j];
	}
    }

    // Find all solutions and store in memory
    void Search()
    {
	Search_(0);
    }

    // Get all solutions
    const std::vector<Solution>& GetSolutions() const { return m_solutions;}

    // Return number of solutions
    const size_t SolutionCount() const { return m_solutions.size(); }

    const std::vector<size_t>& RowToColumns(size_t r) const { return m_rows.at(r);}
    const Solution& GetSolution(size_t idx) const { return m_solutions.at(idx);}

    // Helper function for converting a solution (subset of rows) to
    // whatever the user defines the covered columns to be
    template<class T>
    std::vector<std::vector<T>> DecodeSolution(const std::vector<T>& columns,
					       const Solution& solution) const
    {
	std::vector<std::vector<T>> decoded_solution;
	for (size_t row : solution) {
	    std::vector<T> decoded_row;
	    const auto& cols = RowToColumns(row);
	    for (size_t col : cols) {
		decoded_row.push_back(columns.at(col));
	    }
	    decoded_solution.push_back(std::move(decoded_row));
	}
	return decoded_solution;
    }

    
private:
   
    /*
      dlx specific functions
    */

    // cover column and remove conflicting rows
    void Cover(Node* col_node)
    {

	col_node->right->left = col_node->left;
	col_node->left->right = col_node->right;

	for (Node* i_node = col_node->down; i_node != col_node; i_node = i_node->down) {
	    for (Node* j_node = i_node->right; j_node != i_node; j_node = j_node->right) {
		j_node->down->up = j_node->up;
		j_node->up->down = j_node->down;

		--col_node->size;
	    }
	}	
    }

    // undo effects of Cover()
    void Uncover(Node* col_node)
    {
	for (Node* i_node = col_node->up; i_node != col_node; i_node = i_node->up) {
	    for (Node* j_node = i_node->left; j_node != i_node; j_node = j_node->left) {
		++col_node->size;

		j_node->down->up = j_node;
		j_node->up->down = j_node;
	    }
	}

	col_node->right->left = col_node;
	col_node->left->right = col_node;
    }

    // heuristic for choosing column
    Node* ChooseColumn()
    {
	Node* best_node = m_root;
	for (Node* col_node = m_root->right; col_node != m_root; col_node = col_node->right) {
	    if (col_node->size < best_node->size) {
		best_node = col_node;
	    }
	}
	return best_node;
    }

    // Recursively search for solutions
    void Search_(size_t k)
    {	
	Node* col_node = ChooseColumn();

	// Potential solution found
	if (col_node == m_root) {
	    Solution solution;
	    std::transform(m_cur_solution.begin(), m_cur_solution.begin() + k,
			   std::back_inserter(solution),
			   [](Node* r) -> Index { return r->row_idx;});
	    if (solution.size()) {
		m_solutions.push_back(std::move(solution));
	    }
	    return;	    
	}
	
	Cover(col_node);
	
	for (Node* r_node = col_node->down; r_node != col_node; r_node = r_node->down) {
	    
	    if (k >= m_cur_solution.size()) {
		m_cur_solution.resize(k * 2 + 1);
	    }
	    
	    m_cur_solution[k] = r_node;
	    
	    for (Node* c_node = r_node->right; c_node != r_node; c_node = c_node->right) {
		Cover(c_node->column);
	    }
	    
	    Search_(k+1);
	    
	    // \todo These two lines necessary?
	    r_node = m_cur_solution[k];
	    col_node = r_node->column;
	    
	    for (Node* c_node = r_node->left; c_node != r_node; c_node = c_node->left) {
		Uncover(c_node->column);
	    }
	}
	
	Uncover(col_node);
    }
    
    /*
      Member variables
     */
    MemoryPool<Node> m_pool; // For managing linked list memory
    Node* m_root;
    
    std::vector<Solution> m_solutions; // soltuions
    std::vector<Node*> m_cur_solution; // for building current solution
    std::vector<std::vector<size_t>> m_rows; // map from row to covered columns
};
    
} // namespace

#endif // DLX_H
