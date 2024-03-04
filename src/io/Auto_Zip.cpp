/*
 * Auto_Zip.cpp
 *
 *  Created on: 14/feb/2011
 *      Author: cdf
 */

#include "Auto_Zip.h"

namespace useful {

Auto_Zip::Auto_Zip() {
	output_stream = NULL;
	filtered_stream = NULL;
}

Auto_Zip::Auto_Zip(const string & filename) throw (File_Not_Found) {
	output_stream = NULL;
	filtered_stream = NULL;
	open(filename);
}

Auto_Zip::Auto_Zip(const char * filename) throw (File_Not_Found) {
	output_stream = NULL;
	filtered_stream = NULL;
	open(filename);
}

void Auto_Zip::copy(const Auto_Zip & copy) {
	output_stream = copy.output_stream;
	filtered_stream = copy.filtered_stream;
	filename = copy.filename;
}

void Auto_Zip::close() {
	if (filtered_stream != NULL) {
		filtered_stream->reset();
		delete filtered_stream;
		filtered_stream = NULL;
	}
	if (output_stream != NULL) {
		delete output_stream;
		output_stream = NULL;
	}
}

void Auto_Zip::open(const char * filename) throw (File_Not_Found) {
	string temp(filename);
	open(temp);
}

void Auto_Zip::open(const string & filename) throw (File_Not_Found) {
	close();

	this->filename = filename;

	if (filename.compare("-") == 0) {
		streambuf * psbuf  = cout.rdbuf();
		filtered_stream = new iostreams::filtering_ostream();
		filtered_stream->rdbuf(psbuf);
		output_stream = NULL;
	} else {
		output_stream = new ofstream(filename.c_str(),std::ios_base::binary);
		if (not (*output_stream))
			throw File_Not_Found(filename);

		filtered_stream = new iostreams::filtering_ostream();
		size_t pos = filename.rfind(".gz");
		if ((pos != string::npos) and (pos == (filename.size()-3)))
			filtered_stream->push(iostreams::gzip_compressor());
		else {
			pos = filename.rfind(".bz2");
			if ((pos != string::npos) and (pos == (filename.size()-4)))
				filtered_stream->push(iostreams::bzip2_compressor());
		}
		filtered_stream->push(*output_stream);
	}
}

/*
streampos Auto_Zip::tellg() {
	if (output_stream != NULL)
		return output_stream->tellg();
	else
		return -1;
}

streampos Auto_Zip::length() {
	if (output_stream != NULL) {
		streampos pos = output_stream->tellg();
		output_stream->seekg (0, ios::end);
		streampos length = output_stream->tellg();
		output_stream->seekg (pos);
		return length;
	} else
		return 0;
}
*/

}
