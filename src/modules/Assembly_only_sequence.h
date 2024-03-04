#ifndef ASSEMBLY_ONLY_SEQUENCE_H_
#define ASSEMBLY_ONLY_SEQUENCE_H_

#include "Assembly.h"

class Assembly_only_sequence : public Assembly {

public:

	Assembly_only_sequence (const Options & Opt, Hash * HashHead);

	// void 	open_output_files (ofstream &, ofstream &, fstream &, fstream &);
	void 	open_output_files (ofstream &, fstream &, fstream &);
	// void 	close_output_files (ofstream &, ofstream &, fstream &, fstream &);
	void 	close_output_files (ofstream &, fstream &, fstream &);
	void	execute_fasta_query (Hash * HashHead);
	void	execute_sam_query (Hash * HashHead);

	// functions used in Assembly_store_layout only
	void	printTempLayoutBlock (bool, fstream & temp_layout_file) { }
	void	update_contig_layout (Contig * contig, Hash * HashHead) { }
	void	compute_layout_vector (Hash *, fstream &, fstream &) { }

	void	update_contig_layout_DEBUG (Contig *, Hash *, ofstream &, ofstream &, ofstream &) { }
	void 	compute_layout_vector_DEBUG (Hash *, fstream &, fstream &, ofstream &,
				ofstream &, ofstream &, ofstream &) { }

private:

	bool	_seed_as_query; // mi serve solo qui -> non uso un file separato per le query
							// in Assembly_store_layout
};

#endif /* ASSEMBLY_ONLY_SEQUENCE_H_ */
