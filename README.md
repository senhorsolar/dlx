# Dancing Links (dlx)
Exact cover solver using the dancing links (DLX) algorithm proposed by Donald Knuth [1]. 

Given a list of subsets over some __U__, an exact cover is a subset of these subsets, such that their union is equal to __U__, and each subset is disjoint from the other. The exact cover is also known as the independent set cover. 

Finding the exact cover is a NP complete problem, so no known polynomial-time solver exists. Donald Knuth detailed a simple and memory efficient way to find all exact covers for a problem represented as a binary matrix. 

The items to cover, or __U__, are represented by the columns in a binary matrix; each row represents a subset, with index _(i, j)_ set to true if row (subset) _i_ contains column (element) _j_. A solution is then a list of rows which satisfy the exact cover constraints (union equal to __U__, and all disjoint).

See example.cpp for an example usage.

# Optimizations

- The circular linked lists are implemented so that the nodes all reside in a single chunk of memory (`std::vector`). This provides a bit of a speed up due to locality.
- The solutions are represented as a list of indices indicating the rows, as opposed to a subset of the rows themselves.

# References

1. https://arxiv.org/abs/cs/0011047
