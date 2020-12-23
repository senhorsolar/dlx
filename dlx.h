// dlx.h
#pragma once

#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>


class Node
{
public:	
	Node* left;
	Node* right;
	Node* up;
	Node* down;
	Node* column;
	
	std::string name;
	int size;

	Node() {
		this->left = this;
		this->right = this;
		this->up = this;
		this->down = this;
	}
	
	Node(Node* column_node) : Node() {
		this->column = column_node;
	}

	// For constructing column nodes
	Node(std::string column_name) : Node() {
		this->name = column_name;
		this->size = INT_MAX;
	}	
};

class DancingLinks
{
protected:
	std::vector<Node*> solution;
	std::vector<std::vector<int> > mat;
	std::vector<std::string> column_names;
		
public:

	DancingLinks(std::string data) : DancingLinks(std::stringstream (data)) {}
	
	DancingLinks(std::stringstream ss) {

		std::string line;
		
		// First line
		std::getline(ss, line);
		std::stringstream s (line);

		// Parse column names
		std::string column_name;
		while (std::getline(s, column_name, ',')) {
			column_names.push_back(column_name);
		}
		
		// Parse data
		std::string val;
		while (std::getline(ss, line)) {
			std::stringstream s (line);
			std::vector<int> v;
			while (getline(s, val, ',')) {
				v.push_back(std::stoi(val));
			}
			mat.push_back(v);
		}
		createNodes();
	}

	/*
	  Initialize doubly linked list and set the root node.
	 */
	DancingLinks(std::vector<std::vector<int> > mat, std::vector<std::string> column_names) : mat(mat), column_names(column_names) {
		createNodes();
	}
	
        /* 
	  Destroy doubly linked lists
	*/
	~DancingLinks() {

		if (root_node) {
			Node* cur_col_node = root_node->right;
			
			while (cur_col_node != root_node) {
				
				Node* cur_row_node = cur_col_node->down;
				
				while (cur_row_node != cur_col_node) {
					cur_row_node = cur_row_node->down;
					delete cur_row_node->up;
				}
				
				cur_col_node = cur_col_node->right;
				delete cur_col_node->left;
			}
			
			delete root_node;
		}
	}

	/*
	  Find solution(s)
	 */
	int search() {
		search_recursive(0);
		return nsolutions;
	}

private:
	Node* root_node;
	unsigned int nsolutions = 0;

	void createNodes() {
		root_node = new Node("h");		
		
		// Save locations of nodes
		std::map<std::pair <int, int>, Node*> nodemap;
		
		// Assuming mat is a 2d matrix
		int nrows = mat.size();
		int ncols = column_names.size();
		
		// Check matrix dims
		for (int i = 0; i < nrows; i++) {
			assert(mat[i].size() == ncols);
		}

		// Initialize column nodes and create up/down connections.
		// Saving creation of left/right connections for the next step.
		Node* prev_col_node = root_node;
		
		for (int j = 0; j < ncols; j++) {
			
			Node* cur_col_node = new Node(column_names[j]);

			// We know beforehand how column nodes are connected to each other,
			// since they're adjacent to each other.
			// This isn't as obvious for the regular nodes since they are sparse.
			prev_col_node->right = cur_col_node;
			cur_col_node->left = prev_col_node;

			// Create up/down connections
			Node* prev_row_node = cur_col_node;
			
			for (int i = 0; i < nrows; i++) {
				
				if (mat[i][j]) {
					
					Node* cur_row_node = new Node(cur_col_node);
					
					cur_row_node->up = prev_row_node;
					prev_row_node->down = cur_row_node;

					// Save locations of nodes for later
					std::pair<int, int> key (i, j);
					nodemap.insert(make_pair(key, cur_row_node));
					prev_row_node = cur_row_node;
					
					cur_col_node->size++;
				}
			}

			// Ensure circular connection between vertical nodes
			prev_row_node->down = cur_col_node;
			cur_col_node->up = prev_row_node;
			prev_col_node = cur_col_node;
			
		}

	        // Ensure circular connection between column nodes      
		prev_col_node->right = root_node;
		root_node->left = prev_col_node;

		// Establish left/right connections
		for (int i = 0; i < nrows; i++) {
			
			Node* prev_node = nullptr;
			Node* first_node;
			Node* cur_node;
			
			for (int j = 0; j < ncols; j++) {
				
				if (mat[i][j]) {

					// Access node at this location
					std::pair<int, int> key (i, j);
					cur_node = nodemap.at(key);

					// Get first node in row
					if (prev_node == nullptr) {
						first_node = cur_node;
						prev_node = cur_node;
					} else {
						prev_node->right = cur_node;
						cur_node->left = prev_node;
						prev_node = cur_node;
					}
					
				}
				
			}

			// Ensure circular connection between horizontal nodes (if that row exists)
			if (prev_node != nullptr) {
				first_node->left = cur_node;
				cur_node->right = first_node;
			}
			
		}
		
	}
	
	virtual void display_solution(int q) {

		std::cout << "solution found:" << '\n';
		
		for (int k = 0; k < q; k++) {
			Node* j_node = solution[k];

			std::cout << j_node->column->name << " ";

			for (Node* i_node = j_node->right; i_node != j_node; i_node = i_node->right) {
				std::cout << i_node->column->name << " ";
			}

			std::cout << '\n';
		}
	}
	
	void cover(Node* col_node) {
		col_node->right->left = col_node->left;
		col_node->left->right = col_node->right;
		
		for (Node* i_node = col_node->down; i_node != col_node; i_node = i_node->down) {
			for (Node* j_node = i_node->right; j_node != i_node; j_node = j_node->right) {
				j_node->down->up = j_node->up;
				j_node->up->down = j_node->down;
				
				col_node->size -= 1;
			}
		}	
	}
	
	void uncover(Node* col_node) {
		for (Node* i_node = col_node->up; i_node != col_node; i_node = i_node->up) {
			for (Node* j_node = i_node->left; j_node != i_node; j_node = j_node->left) {
				col_node->size += 1;
				
				j_node->down->up = j_node;
				j_node->up->down = j_node;
			}
		}
		
		col_node->right->left = col_node;
		col_node->left->right = col_node;
	}
	
	Node* choose_column() {
		Node* best_node = root_node->right;
		Node* c = best_node;
		
		while (c != root_node) {
			c = c->right;
			if (c->size < best_node->size)
				best_node = c;
		}
		
		return best_node;
	}

	void search_recursive(int k) {
		
		// Solution found
		if (root_node->right == root_node) {
			display_solution(k);
			nsolutions += 1;
		}
		
		else {
			Node* c = choose_column();
			cover(c);

			for (Node* r = c->down; r != c; r = r->down) {

				if (k >= solution.size())
					solution.resize(k * 2 + 1);
				
				solution[k] = r;
				
				for (Node* j = r->right; j != r; j = j->right)
					cover(j->column);
				
				search_recursive(k+1);
				
				r = solution[k];
				c = r->column;
				
				for (Node* j = r->left; j != r; j = j->left)
					uncover(j->column);
			}
			
			uncover(c);
			
		}
	}
};
