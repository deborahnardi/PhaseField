#include <iostream>

int *ptr;
int var = 7;
int foo = 21;
ptr = &var;
std::cout << "var: " << var << " is stored at: " << ptr << " *ptr(desreferenciar) = " << *ptr << std::endl;
ptr = &foo;
std::cout << "foo: " << foo << " is stored at: " << ptr << std::endl;

int &ref = var;
std::cout << "ref: " << ref << " is stored at: " << &ref << std::endl;