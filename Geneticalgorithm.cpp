#include "Geneticalgorithm.h"

Geneticalgorithm::Geneticalgorithm() {
	number = 0;
	name = "Name?";
}

Geneticalgorithm::Geneticalgorithm(int numr, std::string nam) {
	number = numr;
	name = nam;
}

std::string Geneticalgorithm::show_name()
{
	return name;
}

int Geneticalgorithm::show_number()
{
	return number;
}
