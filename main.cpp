#include <iostream>
#include "task.h"

int main()
{
	Task T;
	T.init();
	T.setParams();
	T.solve();
	T.saveResult();
}

