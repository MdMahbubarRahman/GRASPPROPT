#include <iostream>

#include "Tabusearch.h"
#include "Geneticalgorithm.h"
#include "Initialsolution.h"
#include "Feasibilitysearch.h"
#include "Localsearch.h"
#include "Pathrelinking.h"


int main()
{
    int num = 99;
    std::string name = "Hi there";
    
    
    Geneticalgorithm ga(num, name);
    std::cout << "Your name is : " << ga.show_name() << " and your number is :" << ga.show_number() << std::endl;

    
    TabuList list;
    list.addToTabuList(1, 45);
    list.addToTabuList(2, 69);
    list.showTabuList();
    list.clearTabuList();
    list.showTabuList();

    return 0;
}
