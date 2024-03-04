#include "data_structures/Contig_layout.h"

namespace layout {

// ricava vector <list_element> da vector <read_layout>
Contig_layout::Contig_layout (Hash* HashHead, const vector <read_layout> & reads, unsigned int contigNum)
{
	_contigNum = contigNum;
	_indexes.resize (reads.size());
	unsigned int prev_coord = 0;
	unsigned int i = 0;
	for (vector <read_layout>::const_iterator iter = reads.begin(); iter != reads.end(); iter++) {
		list_element l;
		readOriented r;
		r.id = iter->id;
		r.strand = iter->strand;
		l.index = get_index (HashHead, r);
		l.dist = (iter == reads.begin()) ? 0 : (iter->coord - prev_coord);
		_indexes.at(i) = l;
		prev_coord = iter->coord;
		i++;
	}
}

// ricava vector <list_element> da vector <contigRead>
Contig_layout::Contig_layout (const vector <contigRead> & reads, unsigned int contigNum)
{
	_contigNum = contigNum;
	_indexes.resize (reads.size());
	unsigned int prev_pos = 0;
	unsigned int i = 0;
	for (vector <contigRead>::const_iterator iter = reads.begin(); iter != reads.end(); iter++) {
		list_element l;
		l.index = iter->index;
		l.dist = (iter == reads.begin()) ? INFTY : (iter->contPosition - prev_pos);
		_indexes.at(i) = l;
		prev_pos = iter->contPosition;
		i++;
	}
}

Contig_layout::Contig_layout (const Layout_vector & L)
{
	_layout_vector = &L;
}

inline unsigned int
Contig_layout::get_index (Hash * HashHead, readOriented & r)
{
	/*
	unsigned int fing = (r.second) ? HashHead->ComputeFingerprint (r.first, true, 0) :
									 HashHead->ComputeFingerprint (r.first, false, HashHead->overlap - HashHead->blockLength);
									 */
	unsigned int fing = r.strand ? HashHead->ComputeFingerprint (r.id, true, 0) :
									 	 	  HashHead->ComputeFingerprint (r.id, false, HashHead->overlap - HashHead->blockLength);
	unsigned int l = (fing == HashHead->q - 1) ? 2 * HashHead->numReads
											   : HashHead->HASHcounter[fing + 1];
	unsigned long int index = HashHead->HASHcounter[fing];
	//while (index < l and (HashHead->HASHvalues.at(index).first != r.first)) {
	while (index < l and (HashHead->HASHvalues.at(index).id != r.id)) {
		index++;
	}
	return index;
}


// start, end = positions in _pointer (member of Layout_list) 
// which are referred to pairs (read_id, strand)
void 
Contig_layout::add_element (unsigned int start, unsigned int end, unsigned int dist)
{
	contig_elem elem;
	elem.start = start;
	elem.end = end;
	elem.dist = dist;
	_sublists.push_back (elem);
}

void
Contig_layout::print_indexes ()
{
	cout << "Contig reads: " << endl;
	for (vector <list_element>::iterator iter = _indexes.begin(); iter != _indexes.end(); iter++) {
		cout << "(" << iter->index << "," << iter->dist << ") ";
	}
	cout << endl;
}

void 					
Contig_layout::print_sublists (ofstream & out)
{
	out << _contigNum << " " << _sublists.size() << " ";
	for (vector <contig_elem>::const_iterator	c_iter = _sublists.begin();
												c_iter != _sublists.end();
												c_iter++) {
		out << c_iter->start << " " << c_iter->end << " " << c_iter->dist << " ";
	}
	out << endl;
}

/*
void
Contig_layout::print_subvectors (ofstream & out)
{
	if (_layout_vector == NULL) return;
	out << _contigNum << " ";
	unsigned int coord = 0;
	for (vector <contig_elem>::iterator r_iter = _subvectors.begin();
										r_iter != _subvectors.end();
										r_iter++) {
		coord += r_iter->dist;
		read_layout_list r = _layout_vector->_layout_vector.at(r_iter->start);
		out << "(" << r.index.id << "," << coord << "," << r.length << "," << r.index.strand << ") ";
		unsigned int i = r_iter->start;
		while (i != r_iter->end) {
			i++;
			read_layout_list r = _layout_vector->_layout_vector.at(i);
			coord += r.dist;
			out << "(" << r.index.id << "," << coord << "," << r.length << "," << r.index.strand << ") ";
		}
	}
	out << endl;
}
*/

void
Contig_layout::print_subvectors (ofstream & out, vector <read_layout> & v)
{
	if (_layout_vector == NULL) return;
	out << _contigNum << " ";
	for (vector <read_layout>::iterator v_iter = v.begin(); v_iter != v.end(); v_iter++) {
		out << "(" << v_iter->id << "," << v_iter->coord << ","
			<< v_iter->length << "," << v_iter->strand << ") ";
	}
	out << endl;
}

void
Contig_layout::get_reads_layout (vector <read_layout> & v)
{
	if (_layout_vector == NULL) return;
	unsigned int coord = 0;
	for (vector <contig_elem>::iterator r_iter = _subvectors.begin();
										r_iter != _subvectors.end();
										r_iter++) {
		coord += r_iter->dist;
		read_layout_list r = _layout_vector->_layout_vector.at(r_iter->start);
		read_layout r1;
		r1.id = r.index.id;
		r1.coord = coord;
		r1.length = r.length;
		r1.strand = r.index.strand;
		v.push_back(r1);
		unsigned int i = r_iter->start;
		while (i != r_iter->end) {
			i++;
			read_layout_list r = _layout_vector->_layout_vector.at(i);
			coord += r.dist;
			r1.id = r.index.id;
			r1.coord = coord;
			r1.length = r.length;
			r1.strand = r.index.strand;
			v.push_back(r1);
		}
	}
}

// cambia i "reference" di _sublists in modo che puntino a _layout_vector anziché a _list
void
Contig_layout::compute_subvectors ()
{
	if (_layout_vector == NULL) return;
	unsigned int i = 0;
	_subvectors.resize(_sublists.size());
	for (vector <contig_elem>::iterator s_iter = _sublists.begin(); s_iter != _sublists.end(); s_iter++) {
		// gli elementi di s_iter ora puntano a pointer
		_subvectors.at(i).start = _layout_vector->_positions_in_list.at(s_iter->start);
		_subvectors.at(i).end = _layout_vector->_positions_in_list.at(s_iter->end);
		_subvectors.at(i).dist = s_iter->dist;
		i++;
	}
}

unsigned int
Contig_layout::bits_count ()
{
	Singleton * sing = Singleton::getInstance();
	return (sing->CONTIG_ID_MAX_SIZE() + sing->NUM_SUBLISTS_MAX_SIZE() +
			_sublists.size() * (2 * sing->READ_INDEX_MAX_SIZE() + sing->READ_DIST_MAX_SIZE()));
}

// stampa le informazioni su sublist in data
// NB: i = indice in data, j = indice del bit corrente
// stampo un contig dopo l'altro senza interruzioni (no padding)
void
Contig_layout::add_sublists_to_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	// stampo contigNum
	add_field_to_array (b._data, b._i, b._j, _contigNum, sing->CONTIG_ID_MAX_SIZE());
	add_field_to_array (b._data, b._i, b._j, _sublists.size(), sing->NUM_SUBLISTS_MAX_SIZE());
	for (vector <contig_elem>::iterator s_iter = _sublists.begin();
										s_iter != _sublists.end();
										s_iter++) {
		add_field_to_array (b._data, b._i, b._j, s_iter->start, sing->READ_INDEX_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, s_iter->end, sing->READ_INDEX_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, s_iter->dist, sing->READ_DIST_MAX_SIZE());
	}
}

void
Contig_layout::add_subvectors_to_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	// stampo contigNum
	add_field_to_array (b._data, b._i, b._j, _contigNum, sing->CONTIG_ID_MAX_SIZE());
	add_field_to_array (b._data, b._i, b._j, _subvectors.size(), sing->NUM_SUBLISTS_MAX_SIZE());
	for (vector <contig_elem>::iterator s_iter = _subvectors.begin();
										s_iter != _subvectors.end();
										s_iter++) {
		add_field_to_array (b._data, b._i, b._j, s_iter->start, sing->READ_INDEX_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, s_iter->end, sing->READ_INDEX_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, s_iter->dist, sing->READ_DIST_MAX_SIZE());
	}
}

void
Contig_layout::get_sublists_from_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	get_field_from_array (b._data, b._i, b._j, _contigNum, sing->CONTIG_ID_MAX_SIZE());
	unsigned int sublists_size;
	get_field_from_array (b._data, b._i, b._j, sublists_size, sing->NUM_SUBLISTS_MAX_SIZE());
	_sublists.resize (sublists_size);
	for (unsigned int k = 0; k <  sublists_size; k++) {
		get_field_from_array (b._data, b._i, b._j, _sublists.at(k).start, sing->READ_INDEX_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _sublists.at(k).end, sing->READ_INDEX_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _sublists.at(k).dist, sing->READ_DIST_MAX_SIZE());
	}
}

void
Contig_layout::get_subvectors_from_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	get_field_from_array (b._data, b._i, b._j, _contigNum, sing->CONTIG_ID_MAX_SIZE());
	unsigned int subvectors_size;
	get_field_from_array (b._data, b._i, b._j, subvectors_size, sing->NUM_SUBLISTS_MAX_SIZE());
	_subvectors.resize (subvectors_size);
	for (unsigned int k = 0; k < subvectors_size; k++) {
		get_field_from_array (b._data, b._i, b._j, _subvectors.at(k).start, sing->READ_INDEX_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _subvectors.at(k).end, sing->READ_INDEX_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _subvectors.at(k).dist, sing->READ_DIST_MAX_SIZE());
	}
}

}
