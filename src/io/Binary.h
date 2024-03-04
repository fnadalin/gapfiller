#ifndef BINARY_H_
#define BINARY_H_

namespace layout { class Contig_layout; class Layout_vector; }

#include <iostream>
#include <fstream>
#include <vector>
#include <stdint.h>
using namespace std;

#include "options/Singleton.h"

// #define BUFFER_SIZE 4194304
#define BUFFER_SIZE 4 * 8 * 1024

class Binary {

public:

	Binary () { _dim_data = 0; _i = 0; _j = 0; }
	Binary (Binary & b); // copy constructor
	~Binary () { if (_dim_data > 0) delete [] _data; }

	unsigned int 		i() { return _i; }
	unsigned int 		j() { return _j; }
	unsigned int 		dim_data () { return _dim_data; }
	uint32_t * 			data () { return _data; }
	uint32_t			block_size() { return _block_size; }

	void				initialize ();
	void				initialize (uint32_t block_size);

	static void 		char_to_uint32 (char * w, uint32_t * v, unsigned int w_size, unsigned int v_size);
	static void 		uint32_to_char (uint32_t * v, char * w, unsigned int v_size, unsigned int w_size);

	void				write_uint32_array (fstream & out);
	void				write_block (fstream & out);

	bool 				read_uint32_array (fstream & in, unsigned int dim_data); // layout vector
	bool				read_block (fstream & in); // block of contig sublists

	friend class 		Output;

	friend class 		layout::Contig_layout;
	friend class		layout::Layout_vector;

private:

	uint32_t * 			_data;
	unsigned int 		_dim_data; // number of uint32_t numbers in _data []
	uint32_t			_block_size; // number of bits stored in _data
	unsigned int 		_i; // index in _data;
	unsigned int 		_j; // index in _data[i];
};

#endif /* BINARY_H_ */
