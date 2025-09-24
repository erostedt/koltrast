#pragma once

#include <iostream>

#define CHECK(Expr) _check(#Expr, Expr, __FILE__, __LINE__, "")
#define CHECK_Msg(Expr, Msg) _check(#Expr, Expr, __FILE__, __LINE__, Msg)

constexpr inline void _check(const char *expr_str, bool expr, const char *file, int line, const char *msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";
        std::exit(EXIT_FAILURE);
    }
}
