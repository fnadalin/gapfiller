#include "modules/Assembly_store_layout.h"

Assembly_store_layout::Assembly_store_layout (const Options & Opt, Hash * HashHead) : Assembly (Opt, HashHead)
{
	_layout_list.initialize (HashHead->numReads);
	Singleton * sing = Singleton::getInstance();
	sing->set_global_variables (Opt.k, Opt.blockLength, Opt.overlap, HashHead->maxReadLength);
	sing->set_layout_variables (HashHead->numReads, HashHead->maxReadLength);
	// sing->printInstance();
	if (not sing->feasibility_check()) {
		cerr << "Input dataset is too big: cannot handle it, sorry!" << endl;
		exit(2);
	}

	_contig_bits = 0;

	_temp_layout_name = Opt.outputPrefix;
	_temp_layout_name.append(".temp_layout");
	_layout_name = Opt.outputPrefix;
	_layout_name.append(".layout");
}

void
Assembly_store_layout::open_output_files (ofstream & stats_file,
		fstream & temp_layout_file, fstream & layout_file)
{
	// fasta_file.open (_fasta_name.c_str());
	stats_file.open (_stats_name.c_str());
	temp_layout_file.open (_temp_layout_name.c_str(), fstream::out | ios::binary);
	layout_file.open (_layout_name.c_str(), fstream::out | ios::binary);
}

void
Assembly_store_layout::close_output_files (ofstream & stats_file,
		fstream & temp_layout_file, fstream & layout_file)
{
	layout_file.close();
	temp_layout_file.close();
	stats_file.close();
	// fasta_file.close();
}

void
Assembly_store_layout::update_contig_layout (Contig * contig, Hash * HashHead)
{
	// delete reads that do not match (those trimmed during extension)
	contig->checkNonRepType (HashHead);

	// LAYOUT LIST
	if (contig->returnStatus() == Contig::MatePairFound) {
		Contig_layout C(contig->reads(), contig->contigNum()); // coord -> dist
		_layout_list.add_contig(C); // dist -> sublists
		// update Layout block
		_contig_vector.push_back(C);
		_contig_bits += C.bits_count();
	}
}

void
Assembly_store_layout::update_contig_layout_DEBUG (Contig * contig, Hash * HashHead, ofstream & File1,
		ofstream & File2, ofstream & File3)
{
	// delete reads that do not match (those trimmed during extension)
	contig->checkNonRepType (HashHead);

	// LAYOUT LIST
	if (contig->returnStatus() == Contig::MatePairFound) {
		Contig_layout C(contig->reads(), contig->contigNum()); // coord -> dist
		_layout_list.add_contig(C); // dist -> sublists
		// update Layout block
		_contig_vector.push_back(C);
		_contig_bits += C.bits_count();

		_layout_list.check_representation (C);  // check if contig layout is correctly reported
		contig->printLayout (HashHead, File1); // original layout
		C.print_sublists (File2); // layout in list format
		_layout_list.print_contigRead (HashHead, File3, C); // layout obtained from list format

	}
}

void
Assembly_store_layout::printTempLayoutBlock (bool end, fstream & temp_layout_file)
{
	if (not end and _contig_bits + 32 < 8 * BUFFER_SIZE) {
		return;
	}
	// stampo tramite un vettore di uint32_t
	Binary b;
	b.initialize (_contig_bits);
	for (vector <Contig_layout>::iterator v_iter = _contig_vector.begin();
										  v_iter != _contig_vector.end();
										  v_iter++) {
		// aggiunge la _sublist relativa al contig corrente
		v_iter->add_sublists_to_array (b);
	}
	b.write_block (temp_layout_file);
	_contig_vector.clear();
	_contig_bits = 0;
}

void
Assembly_store_layout::printLayoutVector (Layout_vector & L_vector, fstream & layout_file)
{
	Singleton * sing = Singleton::getInstance();
	{
		Binary b;
		b.write_uint32_array (layout_file); // print info: READ_ID_MAX_SIZE READ_LEN_MAX_SIZE
	}

	Binary b;
	unsigned int el_size = sing->READ_ID_MAX_SIZE() + sing->READ_DIST_MAX_SIZE()
						 + sing->READ_LEN_MAX_SIZE() + 1; // size of an element of layout_vector
	uint32_t block_size = L_vector.get_layout_vector().size() * el_size;
	b.initialize (block_size); // inizializza data[]

	// inserisce in data[] le info su layout_vector
	L_vector.add_layout_vector_to_array (b);
	b.write_block (layout_file);
}

void
Assembly_store_layout::printSubvectorsFromSublists (Layout_vector & L_vector,
		fstream & temp_layout_file, fstream & layout_file)
{
	temp_layout_file.close();
	temp_layout_file.open(_temp_layout_name.c_str(), fstream::in | ios::binary);

	while (temp_layout_file.good()) { // scandisce tutti i blocchi
		Binary b;
		if (b.read_block (temp_layout_file)) { // salvo un blocco in b._data

			Binary b_new(b); // copy of b (except for data[])
			b_new.initialize (); // allocate data[] and initialize elements to 0

			int bits_read = 0;
			unsigned int i, j; // indici precedenti in b._data -> mi servono per contare i bit

			while (bits_read < (int) b.block_size()) {
				// qui scandisce uint32_t * ottenuto da char * e riempie C
				i = b.i(); j = b.j();
				Contig_layout C(L_vector);
				C.get_sublists_from_array (b);
				// salva le info su _subvectors fino a che non raggiungo un blocco
				// di dimensione BUFFER_LENGTH, che é equivalente a un blocco di _sublists
				// (i campi hanno la stessa dimensione)
				C.compute_subvectors ();
				// aggiunge a new_data
				C.add_subvectors_to_array (b_new);
				bits_read += 32 * (b.i() - i) + (b.j() - j);
			}

			// qui stampa il blocco su _layout_file
			b_new.write_block (layout_file);
		}
	}
}

void
Assembly_store_layout::compute_layout_vector (Hash * HashHead, fstream & temp_layout_file,
		fstream & layout_file)
{
	double time1 = clock();
	// HASHcounter non mi serve più
	HashHead->delete_HASHcounter ();
	Layout_vector L_vector (&(_layout_list));
	L_vector.compute_positions_in_list ();
	L_vector.compute_layout_vector (HashHead);
	printLayoutVector (L_vector, layout_file); // stampa layout_vector su layout_file
	printSubvectorsFromSublists (L_vector, temp_layout_file, layout_file); // estrae le info sui contig da temp_layout_file
														 // e stampa su layout_file
	double time2 = clock();
	cout << endl;
	cout << "Layout computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl << endl;
}

// calcola il vettore di Contig_layout a partire dal file binario
// poi bisognerà applicare get_reads_layout per ottenere il layout in reads
void
Assembly_store_layout::get_contigs_layout (ofstream & File6, ofstream & File7,
		fstream & layout_file, vector <Contig_layout> & C_vector)
{
	layout_file.close();
	layout_file.open(_layout_name.c_str(), fstream::in | ios::binary);

	Layout_vector L;
	L.get_layout_vector (layout_file); // extract Layout_vector info from binary file
	L.print_layout_vector (File6);

	while (layout_file.good()) {
		Binary b;
		if (b.read_block (layout_file)) {
			while (32 * b.i() + b.j() < b.block_size()) {
				Contig_layout C(L);
				C.get_subvectors_from_array (b);
				// ottiene il layout in reads a partire dai subvectors
				C_vector.push_back(C);
			}
		}
	}
}

// DEBUG FUNCTIONS

void
Assembly_store_layout::print_layout_sublists (ofstream & File5, fstream & temp_layout_file)
{
	temp_layout_file.close();
	temp_layout_file.open(_temp_layout_name.c_str(), fstream::in | ios::binary);
	while (temp_layout_file.good()) {
		Binary b;
		if (b.read_block (temp_layout_file)) { // se ho letto tutto il blocco dal file
			unsigned int i, j;	// i = indice in data [] al passo precedente,
								// j = indice in data [i] al passo precedente
			int bits_read = 0;
			Contig_layout C;
			while (bits_read < (int) b.block_size()) {
				// qui scandisce uint32_t * ottenuto da char * e riempie C
				i = b.i(); j = b.j();
				C.get_sublists_from_array (b);
				C.print_sublists (File5);
				bits_read += 32 * (b.i() - i) + (b.j() - j);
			}
		}
	}
}

void
Assembly_store_layout::print_layout_from_subvectors (ofstream & File6, ofstream & File7,
		fstream & layout_file)
{
	layout_file.close();
	layout_file.open(_layout_name.c_str(), fstream::in | ios::binary);

	Layout_vector L;
	L.get_layout_vector (layout_file); // extract Layout_vector info from binary file
	L.print_layout_vector (File6);

	while (layout_file.good()) {
		Binary b;
		if (b.read_block (layout_file)) {
			while (32 * b.i() + b.j() < b.block_size()) {
				Contig_layout C(L);
				C.get_subvectors_from_array (b);
				// ottiene il layout in reads a partire dai subvectors
				vector <read_layout> v;
				C.get_reads_layout (v);
				C.print_subvectors (File7, v);
				//C.print_subvectors (File7);
			}
		}
	}
}

void
Assembly_store_layout::compute_layout_vector_DEBUG (Hash * HashHead, fstream & temp_layout_file,
	fstream & layout_file, ofstream & File4, ofstream & File5, ofstream & File6, ofstream & File7)
{
	double time1 = clock();

	// cout << "Initialize Layout_vector" << endl;
	Layout_vector L_vector (&(_layout_list));
	// cout << "Compute positions in _list" << endl;
	L_vector.compute_positions_in_list ();
	// cout << "Compute layout_vector" << endl;
	L_vector.compute_layout_vector (HashHead);
	// cout << "Print layout vector" << endl;
	L_vector.print_layout_vector (File4); // ONLY FOR DEBUG
	// cout << "Print layout vector to binary" << endl;
	printLayoutVector (L_vector, layout_file); // stampa layout_vector su layout_file
	//cout << "Print layout sublists (to text)" << endl;
	print_layout_sublists (File5, temp_layout_file); // ONLY FOR DEBUG
	//cout << "Print subvectors" << endl;
	printSubvectorsFromSublists (L_vector, temp_layout_file, layout_file); // estrae le info sui contig da temp_layout_file
														 // e stampa su layout_file
	// cout << "Print layout from subvectors" << endl;
	print_layout_from_subvectors (File6, File7, layout_file); // ONLY FOR DEBUG

	double time2 = clock();
	cout << endl;
	cout << "Layout computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl << endl;
}
