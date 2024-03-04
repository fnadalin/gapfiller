#include "data_structures/Layout_list.h"

namespace layout {

void			
Layout_list::add_contig (Contig_layout & C)
{
	// here C contains _indexes only
	if (C.get_indexes().size() == 0)
	{
		// error in C constructor
		cerr << "Error in Contig_layout constructor" << endl;
		return;
	}

	// list_1 = previous list
	// list_2 = current list
	list_pointer list_1_start = _list.end(); // pointers to the starting and ending positions of list_1
	list_pointer list_1_end = _list.end();		
	list_pointer list_2_start = _list.end(); // pointers to the starting and ending positions of list_2
	list_pointer list_2_end = _list.end();		

	unsigned int dist_list_1_list_2 = INFTY; // distance among reads pointed by list_1_end and list_2_start
	unsigned int prev_dist = INFTY;

	// index in pointer array, corresponding to list pointers
	unsigned int pointer_pos_1_start = C.get_indexes().begin()->index;
	unsigned int pointer_pos_1_end = pointer_pos_1_start;
	unsigned int pointer_pos_2_start = pointer_pos_1_start;
	unsigned int pointer_pos_2_end = pointer_pos_1_start;

	vector <list_element>::const_iterator prev_c_iter = C.get_indexes().end(); // pointer to previous element in contig layout
	list <list_element> temp_list;

	unsigned int i = 0;
	for (vector <list_element>::const_iterator c_iter = C.get_indexes().begin(); c_iter != C.get_indexes().end(); c_iter++) {
		unsigned int pointer_pos = c_iter->index;
		if (pointer_pos >= _pointer.size()) {
			cout << "Pointer pos = " << pointer_pos << ", pointer size = " << _pointer.size() << endl;
			return;
		}
		unsigned int current_distance = INFTY;
		if (prev_c_iter != C.get_indexes().end()) {
			current_distance = c_iter->dist;
			if (dist_list_1_list_2 == INFTY) {
				dist_list_1_list_2 = current_distance;
			}
		}
		if (list_2_end != _list.end()) { // non mi trovo sul primo elemento 
			// Controlla se vale il caso 1
			list_pointer temp_p (list_2_end);
			temp_p++;
			if (_pointer.at(pointer_pos) != _list.end()) { // la read è già presente nella lista
				if (temp_p == _pointer.at(pointer_pos)) { // la read è contigua alla lista corrente
					temp_p->dist = current_distance; // aggiorno la distanza
					list_2_end++;
				} else { // la read si trova da un'altra parte
					// aggiornamento list_1 e list_2
					Layout_list::list_comparison type = update_list (list_1_start, list_1_end, list_2_start, list_2_end, dist_list_1_list_2);
					// a seconda del risultato aggiorno i puntatori a list_1
					if (type == Layout_list::uncomparable) {
						// aggiungi un nuovo elemento nel contig
						unsigned int d = (C.get_sublists().size() == 0) ? 0 : prev_dist;
						C.add_element (pointer_pos_1_start, pointer_pos_1_end, d);
					}
					if (type != Layout_list::move) {
						list_1_start = list_2_start;
						pointer_pos_1_start = pointer_pos_2_start;
						// aggiorno la distanza precedente solo se non ho unito le due liste
						prev_dist = dist_list_1_list_2;
					} 
					dist_list_1_list_2 = current_distance;
					list_1_end = list_2_end;
					pointer_pos_1_end = pointer_pos_2_end;
					list_2_start = _pointer.at(pointer_pos);
					pointer_pos_2_start = pointer_pos;
					list_2_end = _pointer.at(pointer_pos);
					pointer_pos_2_end = pointer_pos;
				}
			} else { // la read non è presente nella lista 
				if (temp_p == _list.end() or // fine della lista
					temp_p->dist == INFTY) { // break nella lista
					// inserisco l'elemento dopo list_2_end
					list_element el;
					el.index = pointer_pos;
					el.dist = current_distance;
					_list.insert (temp_p, el);
					// aggiorno pointer
					_pointer.at(pointer_pos) = ++list_2_end;
				} else { // non posso inserire dopo list_2_end
					// aggiornamento list_1 e list_2
					Layout_list::list_comparison type = update_list (list_1_start, list_1_end, list_2_start, list_2_end, dist_list_1_list_2);
					// a seconda del risultato aggiorno i puntatori a _list
					if (type == Layout_list::uncomparable) {
						// aggiungi un nuovo elemento nel contig
						unsigned int d = (C.get_sublists().size() == 0) ? 0 : prev_dist;
						C.add_element (pointer_pos_1_start, pointer_pos_1_end, d);
					}	
					if (type != Layout_list::move) {
						list_1_start = list_2_start;
						pointer_pos_1_start = pointer_pos_2_start;
						prev_dist = dist_list_1_list_2;
					} 
					dist_list_1_list_2 = current_distance;
					list_1_end = list_2_end;
					pointer_pos_1_end = pointer_pos_2_end;
					// inserisco l'elemento in coda alla lista creando una lista temporanea	
					list <list_element> temp_list;
					list_element el;
					el.index = pointer_pos;
					el.dist = INFTY;
					temp_list.push_back(el);
					// aggiorna pointer e list_2
					_pointer.at(pointer_pos) = temp_list.begin();
					list_2_start = temp_list.begin();
					pointer_pos_2_start = pointer_pos;
					list_2_end = temp_list.begin();
					pointer_pos_2_end = pointer_pos;
					// inserisce l'elemento nella lista
					_list.splice (_list.end(), temp_list);
				}
			}
		} else if (_pointer.at(pointer_pos) == _list.end()) {
			// mi trovo sul primo elemento e non è ancora stato inserito
			// crea una lista temporanea in cui inserisce l'elemento
			list <list_element> temp_list;
			list_element el;
			el.index = pointer_pos;
			el.dist = INFTY;
			temp_list.push_back(el);
			// aggiorna pointer e list_2
			_pointer.at(pointer_pos) = temp_list.begin();
			list_2_start = temp_list.begin();
			pointer_pos_2_start = pointer_pos;
			list_2_end = temp_list.begin();
			pointer_pos_2_end = pointer_pos;
			// inserisce l'elemento nella lista
			_list.splice (_list.end(), temp_list);
		} else {
			// mi trovo sul primo elemento ed è già stato inserito
			list_2_start = _pointer.at(pointer_pos);
			pointer_pos_2_start = pointer_pos;
			list_2_end = _pointer.at(pointer_pos);
			pointer_pos_2_end = pointer_pos;
		}
		prev_c_iter = c_iter;
		pointer_pos_2_end = pointer_pos;
		i++;
	}
	Layout_list::list_comparison type = update_list (list_1_start, list_1_end, list_2_start, list_2_end, dist_list_1_list_2);
	if (type == Layout_list::uncomparable) {
		// aggiungi un nuovo elemento nel contig
		unsigned int d = (C.get_sublists ().size() == 0) ? 0 : prev_dist;
		C.add_element (pointer_pos_1_start, pointer_pos_1_end, d);
	} else {
		// se non ho aggiunto un nuovo elemento, la distanza resta quella precedente
		dist_list_1_list_2 = prev_dist;
	}
	if (type != Layout_list::move) {
		pointer_pos_1_start = pointer_pos_2_start;
	}
	pointer_pos_1_end = pointer_pos_2_end;
	unsigned int d = (C.get_sublists ().size() == 0) ? 0 : dist_list_1_list_2;
	C.add_element (pointer_pos_1_start, pointer_pos_1_end, d);
}

Layout_list::list_comparison
Layout_list::update_list (const list_pointer list_1_start, const list_pointer list_1_end, 
			 const list_pointer list_2_start, const list_pointer list_2_end, unsigned int dist_list_1_list_2)
{
	//cout << "Update list" << endl;
	if (list_1_start == _list.end()) {
		//cout << "Undefined list 1" << endl;
		return Layout_list::undef;
	}
	if (list_1_end->index == list_2_start->index) {
		// caso di due reads uguali e contigue
		return Layout_list::uncomparable;
	}
	list_pointer t1 = list_1_end; 
	t1++;
	list_pointer t2 = list_2_end;
	t2++;
	if ((t1 == _list.end() or t1->dist == INFTY) and list_2_start->dist == INFTY) {
		if (t2 == _list.end() or t2->dist == INFTY) {
			// sposta list_2 dopo list_1
			list_2_start->dist = dist_list_1_list_2; // aggiorna la distanza
			if (t1 != _list.end()) { // t1->dist = INFTY
				// controlla se t1 appartiene a [list_2_start,t2), in tal caso non fa nulla
				for (list_pointer i = list_2_start; i != t2; i++) {
					if (i == t1) {
						return Layout_list::uncomparable;
					}
				}
				_list.splice (t1, _list, list_2_start, t2);
			} else { // t1 = _list.end()
				// crea una lista temporanea in cui sposta l'elemento
				list <list_element> temp_list;
				temp_list.splice (temp_list.begin(), _list, list_2_start, t2);
				// inserisce temp_list alla fine della lista
				_list.splice (_list.end(), temp_list);	
			}
			//cout << "Move list 2 after list 1" << endl;
			return Layout_list::move;
		} else if (list_1_start->dist == INFTY) {
			// sposta list_1 prima di list_2
			list_2_start->dist = dist_list_1_list_2; // aggiorna la distanza
			if (t2 != _list.end()) { // t2->dist = INFTY
				// controlla se t2 appartiene a [list_1_start,t1), in tal caso non fa nulla
				for (list_pointer i = list_1_start; i != t1; i++) {
					if (i == t2) {
						return Layout_list::uncomparable;
					}
				}
				_list.splice (list_2_start, _list, list_1_start, t1);
			} else { // t2 = _list.end()
				// crea una lista temporanea in cui sposta l'elemento
				list <list_element> temp_list;
				temp_list.splice (temp_list.begin(), _list, list_1_start, t1);
				// inserisce temp_list alla fine della lista
				_list.splice (_list.end(), temp_list);	
			}
			//cout << "Move list 1 before list 2" << endl;
			return Layout_list::move;
		} else {
			//cout << "Uncomparable lists" << endl;
			return Layout_list::uncomparable;
		}
	} else {
		//cout << "Uncomparable lists" << endl;
		return Layout_list::uncomparable;
	}
}

void						
Layout_list::print_list () 
{
	unsigned int num_pieces = 0;
	unsigned int num_reads = 0;
	cout << "List: " << endl;
	for (list <list_element>::iterator l_iter = _list.begin(); l_iter != _list.end(); l_iter++) {
		cout << "(" << l_iter->index << "," << l_iter->dist << ") ";
		num_pieces += (l_iter->dist == INFTY);
		num_reads++;
	}
	cout << endl;
	cout << "Number of pieces: " << num_pieces << endl;
	cout << "Number of used reads: " << num_reads << endl;
}

void 
Layout_list::list_stats () 
{
	unsigned long int num_pieces = 0;
	cout << "Number of used reads: " << _list.size() << endl;
	for (list <list_element>::iterator l_iter = _list.begin(); l_iter != _list.end(); l_iter++) {
		num_pieces += (l_iter->dist == INFTY);
	}
	cout << "Number of sub-lists: " << num_pieces << endl;
}

void						
Layout_list::print_contig_list (Hash * HashHead, Contig_layout & C)
{
	cout << "contig" << C.contigNum() << " layout by list: " << endl;
	for (vector <contig_elem>::const_iterator	c_iter = C.get_sublists().begin();
												c_iter != C.get_sublists().end();
												c_iter++) {
		readOriented r1= HashHead->HASHvalues.at(c_iter->start);
		readOriented r2 = HashHead->HASHvalues.at(c_iter->end);
		/*
		cout << "(" << r1.first << "," << r1.second << ") "
			 << "(" << r2.first << "," << r2.second << ") "
			 << c_iter->dist << endl;
			 */
		cout << "(" << r1.id << "," << r1.strand << ") "
			 << "(" << r2.id << "," << r2.strand << ") "
			 << c_iter->dist << endl;
	}
}

void 
Layout_list::print_contigRead (Hash * HashHead, ofstream & out, Contig_layout & C)
{
	out << C.contigNum() << " ";
	unsigned int coord = 0;
	for (vector <contig_elem>::const_iterator	c_iter = C.get_sublists().begin();
												c_iter != C.get_sublists().end();
												c_iter++) {
		list_pointer l = _pointer.at(c_iter->start); // puntatore al primo elemento della lista
		coord += c_iter->dist;
		readOriented r = HashHead->HASHvalues.at(l->index);
		/*
		out << "(" << r.first << "," << coord << "," << HashHead->readsMulti.at(r.first).length()
			<< "," << r.second << ") ";
			*/
		out << "(" << r.id << "," << coord << "," << HashHead->readsMulti.at(r.id).length()
			<< "," << r.strand << ") ";
		while (l->index != c_iter->end) { // mi fermo quando arrivo all'ultimo elemento
			l++;
			coord += l->dist;
			r = HashHead->HASHvalues.at(l->index);
			/*
			out << "(" << r.first << "," << coord << "," << HashHead->readsMulti.at(r.first).length()
				<< "," << r.second << ") ";
				*/
			out << "(" << r.id << "," << coord << "," << HashHead->readsMulti.at(r.id).length()
				<< "," << r.strand << ") ";
		}
	}
	out << endl;
}

void 
Layout_list::check_representation (Contig_layout & C)
{
	vector <list_element>::const_iterator r_iter = C.get_indexes().begin();
	for (vector <contig_elem>::const_iterator	c_iter = C.get_sublists().begin();
												c_iter != C.get_sublists().end();
												c_iter++) {
		list_pointer l = _pointer.at(c_iter->start); // puntatore al primo elemento della lista
		if (l->index != r_iter->index) {
			cerr << "contig" << C.contigNum() << ": Error in index ("
				 << r_iter->index << "," << l->index << ")" << endl;
		}
		// check distance reported by C._list
		// (except for first element which is 0 in C.list and INFTY in C.layout)
		if (c_iter == C.get_sublists().begin()) {
			if (r_iter->dist != INFTY and c_iter->dist != 0) {
				cerr << "contig" << C.contigNum() << ": Error in distance ("
					 << r_iter->dist << "," << c_iter->dist << ")" << endl;
			}
		}
		if (c_iter != C.get_sublists().begin() and c_iter->dist != r_iter->dist) {
			cerr << "contig" << C.contigNum() << ": Error in distance ("
				 << r_iter->dist << "," << c_iter->dist << ")" << endl;
		}
		while (l->index != c_iter->end) { // mi fermo quando arrivo all'ultimo elemento
			l++;
			r_iter++;
			if (l->index != r_iter->index) {
				cerr << "contig" << C.contigNum() << ": Error in index ("
					 << r_iter->index << "," << l->index << ")" << endl;
			}
			if (l->dist != r_iter->dist) {
				cerr << "contig" << C.contigNum() << ": Error in distance ("
					 << r_iter->dist << "," << l->dist << ")" << endl;
			}
		}
		r_iter++;
	}
}

}
