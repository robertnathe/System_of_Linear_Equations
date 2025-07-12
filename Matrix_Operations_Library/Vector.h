#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>

template <typename T>
class Vector {
public:
    Vector(int size = 0);

    int getSize() const;

    T& operator[](int index);
    const T& operator[](int index) const;

    void resize(int newSize);

private:
    std::vector<T> data;
};

#include "Vector.tpp"

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec);

#endif // VECTOR_H
