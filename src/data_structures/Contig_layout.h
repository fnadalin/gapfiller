#ifndef CONTIG_LAYOUT_H_
#define CONTIG_LAYOUT_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
using namespace std;

#include "io/Binary.h"
#include "data_structures/include.h"
#include "data_structures/Contig.h"
#include "data_structures/Layout_vector.h"
using namespace contig;


namespace layout {

class Contig_layout {

public:

	Contig_layout () { }
	Contig_layout (Hash * HashHead, const vector <read_layout> & v, unsigned int contigNum);
	Contig_layout (const vector <contigRead> & reads, unsigned int contigNum);
	Contig_layout (const vector <list_element> & indexes, unsigned int contigNum)
	{
		_indexes = indexes;
		_contigNum = contigNum;
	}
	Contig_layout (const Layout_vector & L);

	inline unsigned int 			get_index (Hash * HashHead, readOriented & r);

	vector <list_element> const &	get_indexes () { return _indexes; }
	vector <contig_elem> const &	get_sublists () { return _sublists; }
	vector <contig_elem> const &	get_subvectors () { return _subvectors; }
	unsigned int					contigNum () { return _contigNum; }

	void							clear() { _sublists.clear(); _subvectors.clear(); }
	void							add_element (unsigned int start, unsigned int end, unsigned int dist);

	void							print_indexes ();
	void 							print_sublists (ofstream & out);
	void							print_subvectors (ofstream & out);

	void							print_subvectors (ofstream & out, vector <read_layout> & v);
	void					 		get_reads_layout (vector <read_layout> &); // restituisce il layout in reads a partire da subvectors

	void							compute_subvectors (); // cambia i "reference" di _sublists in modo
												   // che puntino a _layout_vector anziché a _list
	unsigned int					bits_count (); // conta i bits che servono per stampare _sublists o _subvectors

	void 							add_sublists_to_array (Binary & b);
	void							add_subvectors_to_array (Binary & b);

	void							get_sublists_from_array (Binary & b);
	void							get_subvectors_from_array (Binary & b);

private:

	vector <list_element>			_indexes;	 // vector of pairs <Hashvalue index, distance>
	vector <contig_elem>			_sublists;	 // vector of indexes to pointers delimiting sub-lists
	vector <contig_elem>			_subvectors; // vector of indexes to _layout_vector
	const Layout_vector *			_layout_vector;
	unsigned int					_contigNum;
};

}

#endif /* CONTIG_LAYOUT_H_ */

