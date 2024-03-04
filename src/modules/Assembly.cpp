#include "Assembly.h"

Assembly::Assembly (const Options & Opt, Hash * HashHead)
{
	_seed_fasta = Opt.seed_fasta;
	_seed_sam = Opt.seed_sam;
	_limitedExtension = Opt.limitedExtension;
	_maxReadsToExtend = Opt.maxReadsToExtend;
	_seed_fw_files = Opt.seed_fw_files;
	_seed_rv_files = Opt.seed_rv_files;
	_seed_sam_files = Opt.seed_sam_files;
	_reverse_seed = Opt.reverse_seed;
	_reverse_mate = Opt.reverse_mate;
	_verbose = Opt.verbose;

	_fasta_vector_size = 0;
	_trash_vector_size = 0;
	_gzip = Opt.gzip;
	_bzip2 = Opt.bzip2;

	Singleton * sing = Singleton::getInstance();
	sing->set_global_variables (Opt.k, Opt.blockLength, Opt.overlap, HashHead->maxReadLength);

	Contig::set_static_vars (HashHead->q, HashHead->h, Opt.blockLength, Opt.overlap, HashHead->maxReadLength, Opt.seed_ins,
			Opt.seed_var, Opt.globalMismatch, Opt.no_cycle);
	Extension::set_static_vars (Opt.extThreshold, Opt.seed_ins, Opt.seed_var);

	_numberOfExtendedReads = 0;
	_numberOfExtensions = 0;
	_stopNotPossibleExtensions = 0;
	_stopRepetitiveSequence = 0;
	_stopLengthMaximumExceed = 0;
	_stopMatePairFound = 0;
	_stopReadCycleFound = 0;

	_fasta_name = Opt.outputPrefix;
	_fasta_name.append(".fasta");
	_trash_name = Opt.outputPrefix;
	_trash_name.append("_trash.fasta");
	if (_gzip) {
		_fasta_name.append(".gz");
		_trash_name.append(".gz");
	} else if (_bzip2) {
		_fasta_name.append(".bz2");
		_trash_name.append(".bz2");
	}

	_stats_name = Opt.outputPrefix;
	_stats_name.append(".stats");

}

void
Assembly::printFastaBlock (bool END_OF_INPUT, ostream & fasta_file)
{
	if (not END_OF_INPUT and _fasta_vector_size < 8 * BUFFER_SIZE) { // 4 * 4KB alla volta
		return;
	}
	char* data = new char [_fasta_vector_size];
	unsigned int i = 0;
	for (vector <string>::iterator iter = _fasta_vector.begin(); iter != _fasta_vector.end(); iter++){
		for (unsigned int j = 0; j < iter->size(); j++) {
			data[i] = iter->at(j);
			i++;
		}
	}
	if (fasta_file.good()) {
		fasta_file.write (data, _fasta_vector_size);
		_fasta_vector.clear();
		_fasta_vector_size = 0;
	}
	delete [] data;
}

void
Assembly::printTrashBlock (bool END_OF_INPUT, ostream & trash_file)
{
	if (not END_OF_INPUT and _trash_vector_size < 8 * BUFFER_SIZE) { // 4 * 4KB alla volta
		return;
	}
	char* data = new char [_trash_vector_size];
	unsigned int i = 0;
	for (vector <string>::iterator iter = _trash_vector.begin(); iter != _trash_vector.end(); iter++){
		for (unsigned int j = 0; j < iter->size(); j++) {
			data[i] = iter->at(j);
			i++;
		}
	}
	if (trash_file.good()) {
		trash_file.write (data, _trash_vector_size);
		_trash_vector.clear();
		_trash_vector_size = 0;
	}
	delete [] data;
}

void
Assembly::printStatsHeader (ofstream & stats_file)
{
	if (stats_file.is_open()) {
		stats_file	<< "read_extended mean_extension_length StopbyNoMoreExt StopByRepSeq "
					<< "StopByLengthExceed StopByPairFound StopbyReadCycle" << endl;
	}
}

void
Assembly::printStatsLine (ofstream & stats_file)
{
	if (stats_file.is_open()) {
		double PartialMean = (double) _numberOfExtensions / _numberOfExtendedReads;
		stats_file	<< _numberOfExtendedReads		<< " " << PartialMean  				<< " "
					<< _stopNotPossibleExtensions	<< " " << _stopRepetitiveSequence	<< " "
					<< _stopLengthMaximumExceed 	<< " " << _stopMatePairFound 		<< " "
					<< _stopReadCycleFound 			<< endl;
	}
}

void
Assembly::execute (Hash * HashHead)
{
	double time1 = clock();

	// ofstream fasta_file;
	Auto_Zip A(_fasta_name);
	ostream &fasta_file = A.filtered();
	Auto_Zip T(_trash_name);
	ostream &trash_file = T.filtered();
	ofstream stats_file;
	fstream temp_layout_file;
	fstream layout_file;
	// open_output_files (fasta_file, stats_file, temp_layout_file, layout_file);
	open_output_files (stats_file, temp_layout_file, layout_file);
	printStatsHeader (stats_file);

	unsigned int num_contigs = HashHead->numReads / 2;
	unsigned int i = 0;
	bool mate_found = false;
	// itera su tutti i seed
	while (i < num_contigs and (not _limitedExtension or i < _maxReadsToExtend)) {
		i++;
		mate_found = false;
		string seed_string = HashHead->readsMulti.at(2*(i-1)).toString (not _reverse_seed);
		string seed_mate_string = HashHead->readsMulti.at(2*i-1).toString (not _reverse_mate);
		// extend only sufficiently long seeds
		if (seed_string.size() >= Contig::overlap()) {
			mate_found = new_contig (seed_string, seed_mate_string, i, HashHead, true);
			if (mate_found) { // TRUE -> mate found
				// print layout info
				printTempLayoutBlock (false, temp_layout_file); // print layout in sublists (temporary file)
				// print fasta sequences
				printFastaBlock (false, fasta_file);
			} else {
				// print fasta sequences
				printTrashBlock (false, trash_file);
			}
			// print statistics
			if ((_numberOfExtendedReads % 10000) == 0) {
				cout << "."; cout.flush();
				if ((_numberOfExtendedReads % 1000000) == 0) {
					cout << endl;
				}
				printStatsLine (stats_file);
			}
		}
	}

	printStatsLine (stats_file);
	printFastaBlock (true, fasta_file);
	printTrashBlock (true, trash_file);

	double time2 = clock();
	cout << endl;
	cout << "Assembly computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl;

	printTempLayoutBlock (true, temp_layout_file);
	compute_layout_vector (HashHead, temp_layout_file, layout_file);
	close_output_files (stats_file, temp_layout_file, layout_file);

}

void
Assembly::execute_DEBUG (Hash * HashHead)
{
	double time1 = clock();

	Auto_Zip A (_fasta_name);
	ostream & fasta_file = A.filtered();
	ofstream stats_file;
	fstream temp_layout_file;
	fstream layout_file;
	open_output_files (stats_file, temp_layout_file, layout_file);
	printStatsHeader (stats_file);

	///////////////// DEBUG
	ofstream File1;
	ofstream File2;
	ofstream File3;
	ofstream File4;
	ofstream File5;
	ofstream File6;
	ofstream File7;
	File1.open("GF_1_orig_layout");
	File2.open("GF_2_list_layout");
	File3.open("GF_3_layout_from_lists");
	File4.open("GF_4_layout_vector");
	File5.open("GF_5_lists_from_binary");
	File6.open("GF_6_layout_vector_from_binary");
	File7.open("GF_7_layout_from_subvectors");
	//////////////////////////

	unsigned int num_contigs = HashHead->numReads / 2;
	unsigned int i = 0;
	// itera su tutti i seed
	while (i < num_contigs and (not _limitedExtension or i < _maxReadsToExtend)) {
		i++;
		string seed_string = HashHead->readsMulti.at(2*(i-1)).toString (not _reverse_seed);
		string seed_mate_string = HashHead->readsMulti.at(2*i-1).toString (not _reverse_mate);
		// extend only sufficiently long seeds
		if (seed_string.size() >= Contig::overlap()) {
			if (new_contig_DEBUG (seed_string, seed_mate_string, i, HashHead, true, File1, File2, File3)) { // TRUE -> mate found
				// print layout info
				printTempLayoutBlock (false, temp_layout_file); // print layout in sublists (temporary file)
			}
			// print fasta sequences
			printFastaBlock (false, fasta_file);
			// print statistics
			if ((_numberOfExtendedReads % 10000) == 0) {
				cout << "."; cout.flush();
				if ((_numberOfExtendedReads % 1000000) == 0) {
					cout << endl;
				}
				printStatsLine (stats_file);
			}
		}
	}

	printStatsLine (stats_file);
	printFastaBlock (true, fasta_file);

	double time2 = clock();
	cout << endl;
	cout << "Assembly computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl;

	printTempLayoutBlock (true, temp_layout_file);
	compute_layout_vector_DEBUG (HashHead, temp_layout_file, layout_file, File4, File5, File6, File7);
	close_output_files (stats_file, temp_layout_file, layout_file);

	////////////////// DEBUG
	File1.close();
	File2.close();
	File3.close();
	File4.close();
	File5.close();
	File6.close();
	File7.close();
	/////////////////////////
}

// ritorna TRUE se ho trovato la mate
bool
Assembly::new_contig (const string & seed_string, const string & seed_mate_string,
		unsigned int contig_id, Hash * HashHead, bool seed_as_query)
{
	Contig * contig = new Contig (HashHead, seed_string, seed_mate_string, contig_id,
									seed_as_query, _reverse_seed, _reverse_mate, _verbose);
	_numberOfExtendedReads++; // increment number of reads that have been extended
	while (contig->returnStatus() == Contig::Continue) {
		// Mate check only in the positions where it is supposed to occur (insert size)
		contig->MateSearch (HashHead, seed_as_query, _verbose);
		// Extension time
		if (contig->returnStatus() == Contig::Continue) {
			_numberOfExtensions++; // increment extensions number
			contig->InitializeExtension(); // initialize the extension phase
			contig->ComputeTemporaryReads (HashHead, _verbose);
			if (contig->returnStatus() == Contig::Continue) {
				contig->ComputeExtensionReads (HashHead, _verbose);
			}
			contig->DeleteExtension();
		}
		// Update statistics
		switch (contig->returnStatus()) {
		case Contig::LengthExceed : 		_stopLengthMaximumExceed++; break;
		case Contig::NoMoreExtensions :		_stopNotPossibleExtensions++; break;
		case Contig::RepetitiveSequence :	_stopRepetitiveSequence++; break;
		case Contig::MatePairFound : 		_stopMatePairFound++; break;
		case Contig::ReadCycleFound : 		_stopReadCycleFound++; break;
		default : 							contig->setStatus(Contig::Continue); break;
		}
	}
	bool mate_found = (contig->returnStatus() == Contig::MatePairFound);
	if (mate_found) {
		contig->addContig (_fasta_vector, _fasta_vector_size);
	} else {
		contig->addContig (_trash_vector, _trash_vector_size);
	}
	// update contig layout
	update_contig_layout (contig, HashHead);
	delete contig;
	return mate_found;
}

// ritorna TRUE se ho trovato la mate
bool
Assembly::new_contig_DEBUG (const string & seed_string, const string & seed_mate_string,
		unsigned int contig_id, Hash * HashHead, bool seed_as_query, ofstream & File1,
		ofstream & File2, ofstream & File3)
{
	bool mate_found = false;
	Contig * contig = new Contig (HashHead, seed_string, seed_mate_string, contig_id,
									seed_as_query, _reverse_seed, _reverse_mate, _verbose);
	_numberOfExtendedReads++; // increment number of reads that have been extended
	while (contig->returnStatus() == Contig::Continue) {
		// Mate check only in the positions where it is supposed to occur (insert size)
		contig->MateSearch (HashHead, seed_as_query, _verbose);
		// Extension time
		if (contig->returnStatus() == Contig::Continue) {
			_numberOfExtensions++; // increment extensions number
			contig->InitializeExtension(); // initialize the extension phase
			contig->ComputeTemporaryReads (HashHead, _verbose);
			if (contig->returnStatus() == Contig::Continue) {
				contig->ComputeExtensionReads (HashHead, _verbose);
			}
			contig->DeleteExtension();
		}
		// Update statistics
		switch (contig->returnStatus()) {
		case Contig::LengthExceed : 		_stopLengthMaximumExceed++; break;
		case Contig::NoMoreExtensions :		_stopNotPossibleExtensions++; break;
		case Contig::RepetitiveSequence :	_stopRepetitiveSequence++; break;
		case Contig::MatePairFound : 		_stopMatePairFound++; break;
		case Contig::ReadCycleFound : 		_stopReadCycleFound++; break;
		default : 							contig->setStatus(Contig::Continue); break;
		}
	}
	contig->addContig (_fasta_vector, _fasta_vector_size);
	// update contig layout
	// update_contig_layout (contig, HashHead);
	update_contig_layout_DEBUG (contig, HashHead, File1, File2, File3);
	delete contig;
	return mate_found;
}
