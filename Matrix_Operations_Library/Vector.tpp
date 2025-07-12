#include <vector>
#include "Vector.h"

template <typename T>
Vector<T>::Vector(int size) : data(size) {}

template <typename T>
int Vector<T>::getSize() const {
    return data.size();
}

template <typename T>
T& Vector<T>::operator[](int index) {
    return data[index];
}

template <typename T>
const T& Vector<T>::operator[](int index) const {
    return data[index];
}

template <typename T>
void Vector<T>::resize(int newSize) {
    data.resize(newSize);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
    for (int i = 0; i < vec.getSize(); ++i) {
        os << vec[i] << '\n';
    }
    return os;
}
