// Simple memory pool for dlx nodes

#ifndef MEMORY_POOL_H
#define MEMORY_POOL_H

#include <vector>

namespace Dlx {

template<class T>
class MemoryPool
{
public:

    MemoryPool() = default;
    
    MemoryPool(size_t n)
    {
	m_data.resize(n);
    }

    void Resize(size_t n) { m_data.resize(n); }

    T* New()
    {
	if (m_ptr < m_data.size())
	    return &m_data[m_ptr++];
	else
	    return nullptr;
    }

private:
    std::vector<T> m_data;
    size_t m_ptr = 0;
};
    
} // namespace

#endif
