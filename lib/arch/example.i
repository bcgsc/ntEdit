%module example
%{
#include <iostream>
extern void hello();
%}

extern void hello();
