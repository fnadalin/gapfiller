#include "Assembly_only_sequence.h"

Assembly_only_sequence::Assembly_only_sequence (const Options & Opt, Hash * HashHead) : Assembly (Opt, HashHead)
{
	_seed_as_query = Opt.seed_as_query;
}

void
Assembly_only_sequence::open_output_files (ofstream & stats_file,
		fstream & temp_layout_file, fstream & layout_file)
{
	// fasta_file.open (_fasta_name.c_str());
	stats_file.open (_stats_name.c_str());
}

void
Assembly_only_sequence::close_output_files (ofstream & stats_file, fstream & temp_layout_file, fstream & layout_file)
{
	stats_file.close();
	// fasta_file.close();
}

// assemble contigs from fasta seeds and use query for extension
void
Assembly_only_sequence::execute_fasta_query (Hash * HashHead)
{
	double time1 = clock();

	Auto_Zip A(_fasta_name);
	ostream & fasta_file = A.filtered();
	ofstream stats_file;
	fstream temp_layout_file;
	fstream layout_file;
	open_output_files (stats_file, temp_layout_file, layout_file);
	printStatsHeader (stats_file);

	vector<string>::const_iterator iter_seed_fw = _seed_fw_files.begin();
	vector<string>::const_iterator iter_seed_rv = _seed_rv_files.begin();
	// so già che seed_fw e seed_rv hanno la stessa dimensione
	bool file_vector_end = (iter_seed_fw == _seed_fw_files.end());

	Layout_list L; // lista vuota
	unsigned long int i = 0; // contig ID

	// itera su tutti i file
	while (not file_vector_end) {

		Auto_Unzip seed_file_fw;
		Auto_Unzip seed_file_rv;
		file_vector_end = (iter_seed_fw == _seed_fw_files.end());

		if (not file_vector_end) {
			seed_file_fw.open (*iter_seed_fw);
			seed_file_rv.open (*iter_seed_rv);
			istream & seedFileFW = seed_file_fw.filtered();
			istream & seedFileRV = seed_file_rv.filtered();
			Fastq seed_fw;
			Fastq seed_rv;
			Fastq::FASTQ_type type_fw = Fastq::check_FASTQ_type_file(iter_seed_fw->c_str());
			Fastq::FASTQ_type type_rv = Fastq::check_FASTQ_type_file(iter_seed_rv->c_str());
			seed_fw.set_input_type(type_fw);
			seed_rv.set_input_type(type_rv);

			// compute one couple of seed files
			while (not seedFileFW.eof() and not seedFileRV.eof() and
				   (not _limitedExtension or i < _maxReadsToExtend)) {
				i++;
				seedFileFW >> seed_fw;
				string seed_string = (not _reverse_seed) ? seed_fw.get_sequence() : seed_fw.reverse_complement();
				seedFileRV >> seed_rv;
				string seed_mate_string = (_reverse_mate) ? seed_rv.reverse_complement() : seed_rv.get_sequence();

				// extend only sufficiently long seeds
				if (seed_string.size() >= Contig::overlap()) {
					new_contig (seed_string, seed_mate_string, i, HashHead, false);
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
			iter_seed_fw++; iter_seed_rv++;
		}
	}

	printStatsLine (stats_file);
	printFastaBlock (true, fasta_file);
	close_output_files (stats_file, temp_layout_file, layout_file);

	double time2 = clock();
	cout << endl;
	cout << "Assembly computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl << endl;
}

// assemble contigs from sam seeds and use query for extension
void
Assembly_only_sequence::execute_sam_query (Hash * HashHead)
{
	double time1 = clock();

	Auto_Zip A (_fasta_name);
	ostream & fasta_file = A.filtered();
	ofstream stats_file;
	fstream temp_layout_file;
	fstream layout_file;
	open_output_files (stats_file, temp_layout_file, layout_file);
	printStatsHeader (stats_file);

	vector<string>::const_iterator iter_seed_sam = _seed_sam_files.begin();
	bool file_vector_end = (iter_seed_sam == _seed_sam_files.end());

	unsigned long int i = 0; // contig ID
	// itera su tutti i file
	while (not file_vector_end) {

		Sam samFile;
		file_vector_end = (iter_seed_sam == _seed_sam_files.end());
		if (not file_vector_end) {
			samFile.open (*iter_seed_sam);
			if (not samFile.is_open()) {
				cerr << "Error while opening file " << *iter_seed_sam << endl;
				exit(1);
			}
		}

		if (not file_vector_end) {

			// compute one couple of seed files
			while (samFile.read() and
				   (not _limitedExtension or i < _maxReadsToExtend)) {
				i++;
				string seed_string = (not _reverse_seed) ?
										samFile.sequence(not samFile.is_reverse()) :
										samFile.sequence(samFile.is_reverse());
				samFile.read();
				string seed_mate_string = (_reverse_mate) ?
										samFile.sequence(samFile.is_reverse()) :
										samFile.sequence(not samFile.is_reverse());

				// extend only sufficiently long seeds
				if (seed_string.size() >= Contig::overlap()) {
					new_contig (seed_string, seed_mate_string, i, HashHead, false);
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
			iter_seed_sam++;
		}
	}

	printStatsLine (stats_file);
	printFastaBlock (true, fasta_file);
	close_output_files (stats_file, temp_layout_file, layout_file);

	double time2 = clock();
	cout << endl;
	cout << "Assembly computation time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "VM: " << vm_usage << "; RSS: " << rss_usage << endl << endl;
}
