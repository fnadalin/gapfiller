#ifndef LAYOUT_LIST_H_
#define LAYOUT_LIST_H_

#include <list>
#include <vector>
using namespace std;

#include "data_structures/include.h"
#include "data_structures/Contig.h"
#include "data_structures/Contig_layout.h"
using namespace contig;

namespace layout {

class Layout_list {

public:
	
	enum list_comparison {move, uncomparable, undef};

	Layout_list () { }
	Layout_list (const unsigned int num_reads) { _pointer.assign (2*num_reads, _list.end()); }

	void initialize (const unsigned int num_reads) { _pointer.assign (2*num_reads, _list.end()); }

	list <list_element>  		& get_list () { return _list; }
	vector <list_pointer> const & get_pointer () { return _pointer; }

	static unsigned int 		get_index (Hash * HashHead, readOriented & r);

	void						add_contig (Contig_layout & C); // aggiunge l'informazione del contig alla lista
																// e restituisce l'oggetto Contig_list creato
	list_comparison				update_list (const list_pointer list_1_start, const list_pointer list_1_end,
											 const list_pointer list_2_start, const list_pointer list_2_end, unsigned int distance);

	void						print_list ();
	void						list_stats ();
	void						print_contig_list (Hash * HashHead, Contig_layout & C);
	void						print_contigRead (Hash * HashHead, ofstream & out, Contig_layout & C);
	void						check_representation (Contig_layout & C);

private:

	/*
	_list.at(i) = i-esimo elemento della lista, contenente la distanza dall'elemento precedente (INFINITY se c'è un break o i=0)
	_pointer.at(j) = puntatore all'elemento di _list che contiene la j-esima coppia (read_id,strand)
	*/
	list	<list_element>		_list;		// ID e distanze (d = INFINITY indica che c'è un break)
	vector	<list_pointer> 		_pointer;	// indici = read ordinate per ID e strand
											// valori = iteratori (puntatori) a elementi della lista
};

}

#endif /* LAYOUT_LIST_H_ */
