#include "data_structures/Layout_vector.h"
#include "data_structures/Layout_list.h"

namespace layout {

// calcola la posizione degli elementi puntati da _pointer in _list (oggetti di Layout_list)
void
Layout_vector::compute_positions_in_list ()
{
	// initialize _list_positions
	unsigned int positions_in_list_size = _layout_list->get_pointer().size();
	_positions_in_list.assign(positions_in_list_size, _layout_list->get_list().size());
	// stessa lunghezza di _pointer, gli elementi sono inizializzati a _list.size()
	// (così so quali non sono presenti nella lista)
	unsigned int list_pos = 0;
	for (list <list_element>::const_iterator	l_iter = _layout_list->get_list().begin();
												l_iter != _layout_list->get_list().end();
												l_iter++) {
		_positions_in_list.at(l_iter->index) = list_pos;
		list_pos++;
	}
}

// converte _list in un vettore, inoltre aggiunge info sulle reads
/*
void
Layout_vector::compute_layout_vector (Hash* HashHead)
{
	_layout_vector.resize (_layout_list->get_list().size());
	unsigned int i = 0;
	for (list <list_element>::const_iterator	l_iter = _layout_list->get_list().begin();
												l_iter != _layout_list->get_list().end();
												l_iter++) {
		//unsigned int position_in_readsMulti = HashHead->HASHvalues.at(l_iter->index).first;
		unsigned int position_in_readsMulti = HashHead->id(HashHead->HASHvalues.at(l_iter->index));
		_layout_vector.at(i).id = position_in_readsMulti;
		_layout_vector.at(i).dist = l_iter->dist;
		_layout_vector.at(i).length = HashHead->readsMulti.at(position_in_readsMulti).length();
		//_layout_vector.at(i).strand =  HashHead->HASHvalues.at(l_iter->index).second;
		_layout_vector.at(i).strand =  HashHead->strand(HashHead->HASHvalues.at(l_iter->index));
		i++;
	}
}
*/

// converte _list in un vettore, inoltre aggiunge info sulle reads
void
Layout_vector::compute_layout_vector (Hash* HashHead)
{
	// non alloca il vettore subito ma aggiunge gli elementi uno alla volta
	// man mano che vengono memorizzati in _layout_vector, gli elementi
	// vengono eliminati da _layout_list
	while (not _layout_list->get_list().empty()) {
		read_layout_list r;
		r.index = HashHead->HASHvalues.at(_layout_list->get_list().begin()->index);
		r.dist = _layout_list->get_list().begin()->dist;
		r.length = HashHead->readsMulti.at(r.index.id).length();
		_layout_list->get_list().erase(_layout_list->get_list().begin());
		_layout_vector.push_back(r);
	}
}

/*
void
Layout_vector::add_layout_vector_to_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	for (vector <read_layout_list>::iterator	l_iter = _layout_vector.begin();
												l_iter != _layout_vector.end();
												l_iter++) {
		add_field_to_array (b._data, b._i, b._j, l_iter->index.id, sing->READ_ID_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, l_iter->dist, sing->READ_DIST_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, l_iter->length, sing->READ_LEN_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, l_iter->index.strand, 1);
	}
}
*/

void
Layout_vector::add_layout_vector_to_array (Binary & b)
{
	Singleton * sing = Singleton::getInstance();
	for (vector <read_layout_list>::iterator	l_iter = _layout_vector.begin();
												l_iter != _layout_vector.end();
												l_iter++) {
		uint32_t index = ((l_iter->index.id << 1) | l_iter->index.strand);
		add_field_to_array (b._data, b._i, b._j, index, sing->READ_INDEX_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, l_iter->dist, sing->READ_DIST_MAX_SIZE());
		add_field_to_array (b._data, b._i, b._j, l_iter->length, sing->READ_LEN_MAX_SIZE());
	}
}

/*
void
Layout_vector::get_layout_vector_from_array (Binary & b)
{
	// suppongo di avere la struttura Singleton
	Singleton * sing = Singleton::getInstance();
	unsigned int layout_element_size = sing->READ_ID_MAX_SIZE() + sing->READ_DIST_MAX_SIZE()
									 + sing->READ_LEN_MAX_SIZE() + 1;
	unsigned int layout_vector_size = (b._block_size + layout_element_size - 1) / layout_element_size;
	_layout_vector.resize(layout_vector_size);
	for (unsigned int k = 0; k < layout_vector_size; k++) {
		uint32_t id;
		get_field_from_array (b._data, b._i, b._j, id, sing->READ_ID_MAX_SIZE());
		_layout_vector.at(k).index.id = id;
		get_field_from_array (b._data, b._i, b._j, _layout_vector.at(k).dist, sing->READ_DIST_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _layout_vector.at(k).length, sing->READ_LEN_MAX_SIZE());
		bool strand;
		get_field_from_array (b._data, b._i, b._j, strand);
		_layout_vector.at(k).index.strand = strand;
	}
}
*/

void
Layout_vector::get_layout_vector_from_array (Binary & b)
{
	// suppongo di avere la struttura Singleton
	Singleton * sing = Singleton::getInstance();
	unsigned int layout_element_size = sing->READ_ID_MAX_SIZE() + sing->READ_DIST_MAX_SIZE()
									 + sing->READ_LEN_MAX_SIZE() + 1;
	unsigned int layout_vector_size = (b._block_size + layout_element_size - 1) / layout_element_size;
	_layout_vector.resize(layout_vector_size);
	for (unsigned int k = 0; k < layout_vector_size; k++) {
		uint32_t index;
		get_field_from_array (b._data, b._i, b._j, index, sing->READ_INDEX_MAX_SIZE());
		_layout_vector.at(k).index.id = (index >> 1);
		_layout_vector.at(k).index.strand = (index & 1);
		get_field_from_array (b._data, b._i, b._j, _layout_vector.at(k).dist, sing->READ_DIST_MAX_SIZE());
		get_field_from_array (b._data, b._i, b._j, _layout_vector.at(k).length, sing->READ_LEN_MAX_SIZE());
	}
}

/*
File format: 	READ_ID_MAX_SIZE	READ_LEN_MAX_SIZE	block_size	layout_vector[]
*/
void
Layout_vector::get_layout_vector (fstream & in)
{
	if (not in.good()) return;
	{
		// estrae le costanti
		Binary b;
		if (b.read_uint32_array (in, 2)) { // stores 2 uint32_t; b._data contains READ_ID_MAX_SIZE and READ_LEN_MAX_SIZE
			Singleton * sing = Singleton::getInstance ();
			sing->set_layout_variables_for_vector (b._data[0], b._data[1]); // now parameters are stored
		}
	}
	Binary b;
	if (b.read_block (in))	// now data [] is re-defined; dim_data is updated with current block size
		get_layout_vector_from_array (b); // extract infos from data []
}

void
Layout_vector::print_layout_vector (ofstream & out)
{
	for (vector <read_layout_list>::iterator	iter = _layout_vector.begin();
												iter != _layout_vector.end();
												iter++) {
		out << iter->index.id << " " << iter->dist << " "
			<< iter->length << " " << iter->index.strand << endl;
	}

}

}
