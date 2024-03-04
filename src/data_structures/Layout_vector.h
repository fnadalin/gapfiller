#ifndef LAYOUT_VECTOR_H_
#define LAYOUT_VECTOR_H_

#include <vector>
#include <stdint.h>
using namespace std;

#include "data_structures/include.h"
#include "data_structures/Hash.h"
#include "io/Binary.h"

namespace layout {

class Layout_list;

class Layout_vector {

public :

	Layout_vector () { }
	Layout_vector (Layout_list * L) { _layout_list = L; }

	vector <read_layout_list> const &	get_layout_vector () { return _layout_vector; }
	vector <unsigned int> const &	 	positions_in_list () { return _positions_in_list; }
	Layout_list *						get_layout_list () { return _layout_list; }

	void								compute_positions_in_list (); 	// calcola la posizione degli elementi puntati
													// da _pointer in _list (oggetti di Layout_list)
	void								compute_layout_vector (Hash* HashHead); // converte _list in un vettore
															// inoltre aggiunge info sulle reads
	// output operations
	void 								add_layout_vector_to_array (Binary & b);
	// input operations
	void								get_layout_vector_from_array (Binary & b);
	void								get_layout_vector (fstream & in);
	// ONLY FOR DEBUG
	void								print_layout_vector (ofstream & out); // ONLY FOR DEBUG

	friend class Contig_layout;

private :

	vector <read_layout_list>	_layout_vector; // corresponding elements of _layout_list._list, with the same order and size
	vector <unsigned int>		_positions_in_list;  // position of Hashvalue index in _layout_list
	Layout_list *  				_layout_list; // la distruggo iterativamente mentre calcolo _layout_vector

};

}

#endif /* LAYOUT_VECTOR_H_ */
