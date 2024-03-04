#ifndef INCLUDE_H_
#define INCLUDE_H_

#include <list>
#include <stdint.h>
#include <cmath>
using namespace std;

// info sul layout "naive" rispetto a Reads
typedef struct {
	uint32_t	 	id;
	unsigned int	coord;
	unsigned int	length;
	bool			strand;
} read_layout;

// info sul layout "compresso"
typedef struct {
	uint32_t 	start;  // indice in _pointer (inizio sotto-lista); è lo stesso di HASHvalues
	uint32_t 	end;	// indice in _pointer (fine sotto-lista); è lo stesso di HASHvalues
	uint16_t 	dist;	// distanza dalla precedente sotto-lista
} contig_elem;

// element of dynamic list
typedef struct {
	uint32_t 	index; // position in HASHvalues
	uint16_t 	dist;
} list_element;

typedef list < list_element >::iterator list_pointer;

/*
typedef struct {
	unsigned int 	id; // read id
	unsigned int	dist;
	unsigned int	length;
	bool			strand;
} read_layout_list;
*/

typedef struct {
	uint32_t id : 31, strand : 1;
} readOriented;

// info per vector_layout (stessi elementi della lista dinamica)
typedef struct {
	readOriented 	index; // read id and strand
	uint16_t		dist;
	uint16_t		length;
} read_layout_list;

inline void
add_field_to_array (uint32_t * data, unsigned int & i, unsigned int & j, unsigned int field, unsigned int field_length)
{
	if (j + field_length <= 32) {
		// rimango sull'elemento corrente
		data[i] |= (field << (32 - (j + field_length)));
		j += field_length;
	} else {
		// prefix_len = (32 - j)
		// suffix len = field_size - (32 - j)
		if (j < 32) { // metto qualcosa anche nell'elemento i, alla fine (suffisso 32 - j)
			data[i] |= (field >> (field_length - (32 - j)));
		}
		// aggiungo il resto nell'elemento successivo
		// 1) ricavo il suffisso
		// 2) lo metto all'inizio di data[i+1]
		data[++i] |= (((field << (32 - j)) >> (32 - j)) << (64 - field_length - j));
		j -= (32 - field_length);
	}
}

inline void
get_field_from_array (uint32_t * data, unsigned int & i, unsigned int & j, uint32_t & field, unsigned int field_length)
{
	if (j + field_length <= 32) {
		//uint32_t mask = (((uint32_t) pow(2, field_length) - 1) << (32 - j - field_length));
		uint32_t mask = (((uint32_t) (1 << field_length) - 1) << (32 - j - field_length));
		field = ((data[i] & mask) >> (32 - j - field_length));
		j += field_length;
	} else {
		if (j < 32) {
			uint32_t mask = (uint32_t) pow(2, 32 - j) - 1;
			field = ((data[i] & mask) << (field_length - (32 - j)));
		} else {
			field = 0;
		}
		//uint32_t mask = (((uint32_t) pow(2, field_length - (32 - j)) - 1) << (64 - field_length - j));
		uint32_t mask = (((uint32_t) (1 << (field_length - (32 - j))) - 1) << (64 - field_length - j));
		field |= ((data[++i] & mask) >> (64 - field_length - j));
		j -= (32 - field_length);
	}
}

inline void
get_field_from_array (uint32_t * data, unsigned int & i, unsigned int & j, uint16_t & field, unsigned int field_length)
{
	if (j + field_length <= 32) {
		//uint32_t mask = (((uint32_t) pow(2, field_length) - 1) << (32 - j - field_length));
		uint32_t mask = (((uint32_t) (1 << field_length) - 1) << (32 - j - field_length));
		field = ((data[i] & mask) >> (32 - j - field_length));
		j += field_length;
	} else {
		if (j < 32) {
			uint32_t mask = (uint32_t) pow(2, 32 - j) - 1;
			field = ((data[i] & mask) << (field_length - (32 - j)));
		} else {
			field = 0;
		}
		//uint32_t mask = (((uint32_t) pow(2, field_length - (32 - j)) - 1) << (64 - field_length - j));
		uint32_t mask = (((uint32_t) (1 << (field_length - (32 - j))) - 1) << (64 - field_length - j));
		field |= ((data[++i] & mask) >> (64 - field_length - j));
		j -= (32 - field_length);
	}
}

#endif /* INCLUDE_H_ */
