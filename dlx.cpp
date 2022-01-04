// dlx.cpp

#include "dlx.h"

namespace Dlx {

DancingLinks::DancingLinks(const BinaryMatrix& bin_mat)
{

    if (bin_mat.size() == 0 or bin_mat[0].size() == 0) {
	// \todo throw error
	return;
    }

    size_t nrows = bin_mat.size();
    size_t ncols = bin_mat[0].size();
    
    // Root node
    m_nodes.emplace_back(0);

    // Ensure root node connected to column nodes
    if (ncols > 0) {
	m_nodes[0].right = 1;
	m_nodes[0].left = ncols;
    }

    // For keeping track of last node in column
    std::vector<Index> prev_row;
    prev_row.reserve(ncols);
				
    // Column nodes
    for (size_t j = 1; j <= ncols; ++j) {
	m_nodes.emplace_back(j);
	m_nodes[j].left = j-1;
	m_nodes[j].right = (j == ncols) ? 0 : (j+1);

	prev_row.push_back(j);
    }

    // Construct
    m_rows.reserve(nrows);
    Index cur_idx = ncols + 1;
    for (size_t i = 0; i < nrows; ++i) {
	
	Index first_in_row = cur_idx;
	Index last_in_row = cur_idx;
	bool valid_row = false;

	std::vector<size_t> row;
	
	for (size_t j = 0; j < ncols; ++j) {
	    
	    Index col_idx = j + 1;
	    
	    if (bin_mat[i][j]) {

		valid_row = true;
		row.push_back(j);
		
		// left/right values are not necessarily correct
		// Will be wrong for first and last nodes in row
		// This gets corrected afterwards
		Index left = cur_idx - 1;
		Index right = cur_idx + 1;
		Index up = prev_row[j];

		m_nodes.emplace_back(cur_idx, i, col_idx);
		m_nodes[cur_idx].left = left;
		m_nodes[cur_idx].right = right;
		m_nodes[cur_idx].up = up;
		
		m_nodes[up].down = cur_idx;

		prev_row[j] = cur_idx;
		last_in_row = cur_idx;

		++m_nodes[col_idx].size;	 
		++cur_idx;
	    }
	}
	// Correction for first and last nodes in row
	if (valid_row) {
	    m_nodes[first_in_row].left = last_in_row;
	    m_nodes[last_in_row].right = first_in_row;
	}
	m_rows.push_back(std::move(row));
    }

    // Connect last nodes in columns to respective column nodes
    for (size_t j = 0; j < ncols; ++j) {
	Index col_idx = j + 1;
	m_nodes[prev_row[j]].down = col_idx;
	m_nodes[col_idx].up = prev_row[j];
    }
}

void
DancingLinks::Search()
{
    Search_(0);
}

void
DancingLinks::Cover(Index col_idx)
{
    m_nodes[m_nodes[col_idx].right].left = m_nodes[col_idx].left;
    m_nodes[m_nodes[col_idx].left].right = m_nodes[col_idx].right;

    for (Index i = m_nodes[col_idx].down; i != col_idx; i = m_nodes[i].down) {
	for (Index j = m_nodes[i].right; j != i; j = m_nodes[j].right) {
	    m_nodes[m_nodes[j].down].up = m_nodes[j].up;
	    m_nodes[m_nodes[j].up].down = m_nodes[j].down;

	    --m_nodes[col_idx].size;
	}
    }
}

void
DancingLinks::Uncover(Index col_idx)
{
    for (Index i = m_nodes[col_idx].up; i != col_idx; i = m_nodes[i].up) {
	for (Index j = m_nodes[i].left; j != i; j = m_nodes[j].left) {
	    ++m_nodes[col_idx].size;

	    m_nodes[m_nodes[j].down].up = j;
	    m_nodes[m_nodes[j].up].down = j;
	}
    }

    m_nodes[m_nodes[col_idx].right].left = col_idx;
    m_nodes[m_nodes[col_idx].left].right = col_idx;
}

Index
DancingLinks::ChooseColumn()
{
    Index best_column = m_nodes[0].right;
    for (Index col_idx = best_column; col_idx != 0; col_idx = m_nodes[col_idx].right) {
	if (m_nodes[col_idx].size > m_nodes[best_column].size) {
	    best_column = col_idx;
	}
    }
    return best_column;
}

void
DancingLinks::Search_(size_t k)
{
    // Solution found
    if (m_nodes[0].right == 0) {
	Solution solution;
	std::transform(m_cur_solution.begin(), m_cur_solution.begin()+k,
		       std::back_inserter(solution),
		       [&](Index r) -> Index { return m_nodes[r].row;});
	m_solutions.push_back(std::move(solution));
	return;
    }

    Index col_idx = ChooseColumn();
    Cover(col_idx);

    for (Index r = m_nodes[col_idx].down; r != col_idx; r = m_nodes[r].down) {

	if (k >= m_cur_solution.size()) {
	    m_cur_solution.resize(k * 2 + 1);
	}

	// Add row to current solution
	m_cur_solution[k] = r;

	// Remove conflicting rows (those that cover the same columns)
	for (Index j = m_nodes[r].right; j != r; j = m_nodes[j].right) {
	    Cover(m_nodes[j].column);
	}

	Search_(k+1);

	// \todo These two lines necessary?
	r = m_cur_solution[k]; 
	col_idx = m_nodes[r].column;

	// Restore conflicting rows
	for (Index j = m_nodes[r].left; j != r; j = m_nodes[j].left) {
	    Uncover(m_nodes[j].column);
	}
    }
    Uncover(col_idx);
}
    
} // namespace
