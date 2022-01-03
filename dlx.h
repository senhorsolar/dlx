// dlx.h

#ifndef DLX_H
#define DLX_H

#include <cassert>
#include <climits>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Dlx {

using Index=size_t;
using BinaryMatrix=std::vector<std::vector<bool>>;
using Solution=std::vector<std::vector<size_t>>;
    
struct Node
{
    Index left;
    Index right;
    Index up;
    Index down;

    Index row;
    Index column;
    Index size;

    Node(Index self_idx, Index row_idx, Index col_idx) :
	left(self_idx), right(self_idx), up(self_idx), down(self_idx),
	column(col_idx), row(row_idx), size(0)
    {}
    Node(Index self_idx) :
	left(self_idx), right(self_idx), up(self_idx), down(self_idx),
	column(self_idx), row(-1), size(0)
    {}
};

class DancingLinks
{

public:

    DancingLinks(const BinaryMatrix& bin_mat);

    void Search();

    const std::vector<Solution>& GetSolutions() const {
	return m_solutions;
    }
    
private:
    std::vector<Node> m_nodes;
    std::vector<Solution> m_solutions;
    std::vector<size_t> m_cur_solution;

    void Cover(Index col_idx);
    void Uncover(Index col_idx);
    Index ChooseColumn();
    void Search_(size_t k);
};

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
    Index cur_idx = ncols + 1;
    for (size_t i = 0; i < nrows; ++i) {
	
	Index first_in_row = cur_idx;
	Index last_in_row = cur_idx;
	bool valid_row = false;
	
	for (size_t j = 0; j < ncols; ++j) {
	    
	    Index col_idx = j + 1;
	    
	    if (bin_mat[i][j]) {

		valid_row = true;
		
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
    }

    // Connect last nodes in columns to respective column nodes
    for (size_t j = 0; j < ncols; ++j) {
	Index col_idx = j + 1;
	m_nodes[prev_row[j]].down = col_idx;
	m_nodes[col_idx].up = prev_row[j];
    }

    // for (int i = 0; i < m_nodes.size(); ++i) {
    // 	auto& node = m_nodes[i];
    // 	std::cout << "(" << node.row << "," << node.column-1 << ")\n";
    // 	std::cout << "\tself:\t" << i << '\n';
    // 	std::cout << "\tup:\t" << node.up << '\n';
    // 	std::cout << "\tdown:\t" << node.down << '\n';
    // 	std::cout << "\tleft:\t" << node.left << '\n';
    // 	std::cout << "\tright:\t" << node.right << '\n';
    // }
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
	// std::transform(m_cur_solution.begin(), m_cur_solution.end(),
	// 	       m_cur_solution.begin(),
	// 	       [&](Index r) -> Index { return m_nodes[r].row;});
	// m_solutions.emplace_back(m_cur_solution.begin(),
	// 			 m_cur_solution.begin()+k);
	Solution s;
	for (size_t q = 0; q < k; ++q) {
	    std::vector<size_t> row;
	    Index j = m_cur_solution[q];
	    size_t col = m_nodes[j].column - 1;
	    row.push_back(col);
	    //std::cout << col << " ";
	    for (Index i = m_nodes[j].right; i != j; i = m_nodes[i].right) {
		size_t col = m_nodes[i].column - 1;
		//std::cout << col << " ";
		row.push_back(col);
	    }
	    //std::cout << '\n';
	    s.push_back(std::move(row));
	}
	m_solutions.push_back(std::move(s));
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
    
// class Node
// {
// public:	
// 	Node* left;
// 	Node* right;
// 	Node* up;
// 	Node* down;
// 	Node* column;
	
// 	std::string name;
// 	int size;

// 	Node() {
// 		this->left = this;
// 		this->right = this;
// 		this->up = this;
// 		this->down = this;
// 	}
	
// 	Node(Node* column_node) : Node() {
// 		this->column = column_node;
// 	}

// 	// For constructing column nodes
// 	Node(std::string column_name) : Node() {
// 		this->name = column_name;
// 		this->size = INT_MAX;
// 	}	
// };

// class DancingLinks
// {
// protected:
// 	std::vector<Node*> solution;
// 	std::vector<std::vector<int> > mat;
// 	std::vector<std::string> column_names;
		
// public:

// 	DancingLinks(std::string data) : DancingLinks(std::stringstream (data)) {}
	
// 	DancingLinks(std::stringstream ss) {

// 		std::string line;
		
// 		// First line
// 		std::getline(ss, line);
// 		std::stringstream s (line);

// 		// Parse column names
// 		std::string column_name;
// 		while (std::getline(s, column_name, ',')) {
// 			column_names.push_back(column_name);
// 		}
		
// 		// Parse data
// 		std::string val;
// 		while (std::getline(ss, line)) {
// 			std::stringstream s (line);
// 			std::vector<int> v;
// 			while (getline(s, val, ',')) {
// 				v.push_back(std::stoi(val));
// 			}
// 			mat.push_back(v);
// 		}
// 		createNodes();
// 	}

// 	/*
// 	  Initialize doubly linked list and set the root node.
// 	 */
// 	DancingLinks(std::vector<std::vector<int> > mat, std::vector<std::string> column_names) : mat(mat), column_names(column_names) {
// 		createNodes();
// 	}
	
//         /* 
// 	  Destroy doubly linked lists
// 	*/
// 	~DancingLinks() {

// 		if (root_node) {
// 			Node* cur_col_node = root_node->right;
			
// 			while (cur_col_node != root_node) {
				
// 				Node* cur_row_node = cur_col_node->down;
				
// 				while (cur_row_node != cur_col_node) {
// 					cur_row_node = cur_row_node->down;
// 					delete cur_row_node->up;
// 				}
				
// 				cur_col_node = cur_col_node->right;
// 				delete cur_col_node->left;
// 			}
			
// 			delete root_node;
// 		}
// 	}

// 	/*
// 	  Find solution(s)
// 	 */
// 	int search() {
// 		search_recursive(0);
// 		return nsolutions;
// 	}

// 	int solution_count() {
// 		return nsolutions;
// 	}

// private:
// 	Node* root_node;
// 	unsigned int nsolutions = 0;

// 	void createNodes() {
// 		root_node = new Node("h");		
		
// 		// Save locations of nodes
// 		std::map<std::pair <int, int>, Node*> nodemap;
		
// 		// Assuming mat is a 2d matrix
// 		int nrows = mat.size();
// 		int ncols = column_names.size();
		
// 		// Check matrix dims
// 		for (int i = 0; i < nrows; i++) {
// 			assert(mat[i].size() == ncols);
// 		}

// 		// Initialize column nodes and create up/down connections.
// 		// Saving creation of left/right connections for the next step.
// 		Node* prev_col_node = root_node;
		
// 		for (int j = 0; j < ncols; j++) {
			
// 			Node* cur_col_node = new Node(column_names[j]);

// 			// We know beforehand how column nodes are connected to each other,
// 			// since they're adjacent to each other.
// 			// This isn't as obvious for the regular nodes since they are sparse.
// 			prev_col_node->right = cur_col_node;
// 			cur_col_node->left = prev_col_node;

// 			// Create up/down connections
// 			Node* prev_row_node = cur_col_node;
			
// 			for (int i = 0; i < nrows; i++) {
				
// 				if (mat[i][j]) {
					
// 					Node* cur_row_node = new Node(cur_col_node);
					
// 					cur_row_node->up = prev_row_node;
// 					prev_row_node->down = cur_row_node;

// 					// Save locations of nodes for later
// 					std::pair<int, int> key (i, j);
// 					nodemap.insert(make_pair(key, cur_row_node));
// 					prev_row_node = cur_row_node;
					
// 					cur_col_node->size++;
// 				}
// 			}

// 			// Ensure circular connection between vertical nodes
// 			prev_row_node->down = cur_col_node;
// 			cur_col_node->up = prev_row_node;
// 			prev_col_node = cur_col_node;
			
// 		}

// 	        // Ensure circular connection between column nodes      
// 		prev_col_node->right = root_node;
// 		root_node->left = prev_col_node;

// 		// Establish left/right connections
// 		for (int i = 0; i < nrows; i++) {
			
// 			Node* prev_node = nullptr;
// 			Node* first_node;
// 			Node* cur_node;
			
// 			for (int j = 0; j < ncols; j++) {
				
// 				if (mat[i][j]) {

// 					// Access node at this location
// 					std::pair<int, int> key (i, j);
// 					cur_node = nodemap.at(key);

// 					// Get first node in row
// 					if (prev_node == nullptr) {
// 						first_node = cur_node;
// 						prev_node = cur_node;
// 					} else {
// 						prev_node->right = cur_node;
// 						cur_node->left = prev_node;
// 						prev_node = cur_node;
// 					}
					
// 				}
				
// 			}

// 			// Ensure circular connection between horizontal nodes (if that row exists)
// 			if (prev_node != nullptr) {
// 				first_node->left = cur_node;
// 				cur_node->right = first_node;
// 			}
			
// 		}
		
// 	}
	
// 	virtual void display_solution(int q) {

// 		std::cout << "solution found:" << '\n';
		
// 		for (int k = 0; k < q; k++) {
// 			Node* j_node = solution[k];

// 			std::cout << j_node->column->name << " ";

// 			for (Node* i_node = j_node->right; i_node != j_node; i_node = i_node->right) {
// 				std::cout << i_node->column->name << " ";
// 			}

// 			std::cout << '\n';
// 		}
// 	}

// 	static inline
// 	void cover(Node* col_node) {
// 		col_node->right->left = col_node->left;
// 		col_node->left->right = col_node->right;
		
// 		for (Node* i_node = col_node->down; i_node != col_node; i_node = i_node->down) {
// 			for (Node* j_node = i_node->right; j_node != i_node; j_node = j_node->right) {
// 				j_node->down->up = j_node->up;
// 				j_node->up->down = j_node->down;
				
// 				col_node->size -= 1;
// 			}
// 		}	
// 	}

// 	static inline
// 	void uncover(Node* col_node) {
// 		for (Node* i_node = col_node->up; i_node != col_node; i_node = i_node->up) {
// 			for (Node* j_node = i_node->left; j_node != i_node; j_node = j_node->left) {
// 				col_node->size += 1;
				
// 				j_node->down->up = j_node;
// 				j_node->up->down = j_node;
// 			}
// 		}
		
// 		col_node->right->left = col_node;
// 		col_node->left->right = col_node;
// 	}
	
// 	Node* choose_column() {
// 		Node* best_node = root_node->right;
// 		Node* c = best_node;
		
// 		while (c != root_node) {
// 			c = c->right;
// 			if (c->size < best_node->size)
// 				best_node = c;
// 		}
		
// 		return best_node;
// 	}

// 	void search_recursive(int k) {
		
// 		// Solution found
// 		if (root_node->right == root_node) {
// 			display_solution(k);
// 			nsolutions += 1;
// 		}
		
// 		else {
// 			Node* c = choose_column();
// 			cover(c);

// 			for (Node* r = c->down; r != c; r = r->down) {

// 				if (k >= solution.size())
// 					solution.resize(k * 2 + 1);
				
// 				solution[k] = r;
				
// 				for (Node* j = r->right; j != r; j = j->right)
// 					cover(j->column);
				
// 				search_recursive(k+1);
				
// 				r = solution[k];
// 				c = r->column;
				
// 				for (Node* j = r->left; j != r; j = j->left)
// 					uncover(j->column);
// 			}
			
// 			uncover(c);
			
// 		}
// 	}
// };

} // namespace

#endif // DLX_H
