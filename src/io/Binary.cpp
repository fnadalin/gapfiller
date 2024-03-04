#include "Binary.h"

Binary::Binary (Binary & b)
{
	_dim_data = b._dim_data;
	_data = new uint32_t [_dim_data];
	_block_size = b._block_size;
	_i = 0;
	_j = 0;
}

void
Binary::initialize ()
{
	for (uint32_t k = 0; k < _dim_data; k++) _data[k] = 0;
}

void
Binary::initialize (uint32_t block_size)
{
	_block_size = block_size;
	_dim_data = (_block_size + 31) / 32;
	_data = new uint32_t [_dim_data];
	for (uint32_t k = 0; k < _dim_data; k++) _data[k] = 0;
}

void
Binary::char_to_uint32 (char * w, uint32_t * v, unsigned int w_size, unsigned int v_size)
{
	// suppose arrays are already allocated
	unsigned int j = 0;
	for (unsigned int i = 0; i < v_size; i++) {
		unsigned int k = 0;
		v[i] = 0;
		while (k < 4 and j+k < w_size) {
			v[i] |= ((uint32_t) ((uint8_t) w[j+k]) << 8*(3-k));
			k++;
		}
		j += k;
	}
}

void
Binary::uint32_to_char (uint32_t * v, char * w, unsigned int v_size, unsigned int w_size)
{
	// suppose arrays are already allocated
	unsigned int j = 0;
	for (unsigned int i = 0; i < v_size; i++) {
		unsigned int k = 0;
		while (k < 4 and j+k < w_size) {
			uint32_t mask = ((uint32_t) 255 << 8*(3-k));
			w[j+k] = (char) ((v[i] & mask) >> 8*(3-k));
			k++;
		}
		j += k;
	}
}

void
Binary::write_uint32_array (fstream & out)
{
	if (not out.good()) return;
	_dim_data = 2;
	_data = new uint32_t [_dim_data];
	Singleton * sing = Singleton::getInstance ();
	_data[0] = sing->READ_ID_MAX_SIZE();
	_data[1] = sing->READ_LEN_MAX_SIZE();

	char * char_data = new char [4 * _dim_data];
	uint32_to_char (_data, char_data, _dim_data, 4 * _dim_data);
	out.write (char_data, 4 * _dim_data);
	delete [] char_data;
}

void
Binary::write_block (fstream & out)
{
	if (not out.good()) return;
	char * char_block_size = new char [4];
	uint32_to_char (&_block_size, char_block_size, 1, 4);
	out.write (char_block_size, 4);
	delete [] char_block_size;

	char * char_data = new char [4 * _dim_data];
	uint32_to_char (_data, char_data, _dim_data, 4 * _dim_data);
	out.write (char_data, 4 * _dim_data);
	delete [] char_data;
}

// size of data [] is known
bool
Binary::read_uint32_array (fstream & in, unsigned int dim_data)
{
	if (not in.good()) return 0;
	_dim_data = dim_data;
	char * char_data = new char [4 * dim_data];
	in.read (char_data, 4 * dim_data);

	if (not in.good()) return 0; // termina se non riesce a leggere nulla
	_data = new uint32_t [dim_data];
	char_to_uint32 (char_data, _data, 4 * dim_data, dim_data);
	delete [] char_data;
	return 1;
}

// size of data [] is printed at the beginning of the block
bool
Binary::read_block (fstream & in)
{
	if (not in.good()) return 0;
	char * char_data_size = new char [4]; // block size (in bits)
	in.read (char_data_size, 4);

	if (not in.good()) return 0;
	char_to_uint32 (char_data_size, &_block_size, 4, 1);
	delete [] char_data_size;

	_dim_data = (_block_size + 31) / 32;
	char * char_data = new char [_dim_data * 4];
	in.read (char_data, _dim_data * 4);
	_data =  new uint32_t [_dim_data];
	char_to_uint32 (char_data, _data, _dim_data * 4, _dim_data);
	delete [] char_data;
	_i = 0; _j = 0;
	return 1;
}
