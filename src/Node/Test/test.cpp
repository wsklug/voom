#include "../Node.hpp"
#include "../LoopShellNode.hpp"
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	ifstream ifs;
	ifs.open("node.info",ios::in);
	if(!ifs.good()){
		cout << "can not open input file.." << endl;
		cout << "program terminated." << endl;
		exit(0);
	}
		
	OldNode ob(ifs);
	std::cout << ob;
	ifs >> ob;
	std::cout << ob;

	OldNode ob1(ob);
	std::cout << ob1;
	
	LoopShellNode ob2;
	ifs >> ob2;
	std::cout << ob2;

	ifs >> ob2;
	std::cout << ob2;

	LoopShellNode ob3;
	ob3 = ob2;
	std::cout << ob3;

	LoopShellNode ob4(ob2);
	std::cout << ob4;

	ifs.close();	
	return 0;
}
