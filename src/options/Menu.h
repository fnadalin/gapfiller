#ifndef MENU_H_
#define MENU_H_

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <common.h>

class Menu {

public:

	Menu ();
	void add_options ();
	void parse_options (int argc, char* argv []);

	friend class Options;

private:

	string 					_head;
	po::options_description _menu;
	po::variables_map 		_vm;
};

#endif /* MENU_H_ */
