/*
 * Auto_Zip.h
 *
 *  Created on: 14/feb/2011
 *      Author: cdf
 */

#ifndef AUTO_ZIP_H_
#define AUTO_ZIP_H_

#include <fstream>
#include <iostream>
using namespace std;

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
using namespace boost;

#include "../errors/File_Not_Found.h"
using namespace errors;

namespace useful {

class Auto_Zip {
public:
	Auto_Zip();
	Auto_Zip(const char * filename) throw (File_Not_Found);
	Auto_Zip(const string & filename) throw (File_Not_Found);
	virtual ~Auto_Zip() { close(); }
	Auto_Zip(const Auto_Zip & original) { copy(original); }

	Auto_Zip & operator=(const Auto_Zip & original) {
		if (this != &original) {
			close();
			copy(original);
		}
		return *this;
	}

	void open(const char * filename) throw (File_Not_Found);
	void open(const string & filename) throw (File_Not_Found) ;
	void close();
	ostream & filtered() {
		if (filtered_stream == NULL)
			return *output_stream;
		else
			return *filtered_stream;
	}
	ofstream & file() { return *output_stream; }

	const string & get_filename() const { return filename; }

protected:
	ofstream * output_stream;
	iostreams::filtering_ostream * filtered_stream;

	string filename;

	void copy(const Auto_Zip & copy);

};

}

#endif /* AUTO_ZIP_H_ */
