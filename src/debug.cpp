#include "debug.hpp"
#include <iostream>

void to_dump(const mpz_class& ib) {
	__gmpz_dump(ib.get_mpz_t());
}

void to_dump(const mpq_class& q) {
	to_dump(q.get_num());
	std::cout << "/";
	to_dump(q.get_den());
}
std::string to_string(const mpq_class& ib) {
	return ib.get_str();
}
std::string to_string(const mpz_class& ib) {
	return ib.get_str();
}
