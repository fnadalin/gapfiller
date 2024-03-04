#ifndef ASSEMBLY_STORE_LAYOUT_H_
#define ASSEMBLY_STORE_LAYOUT_H_

#include "modules/Assembly.h"

class Assembly_store_layout : public Assembly {

public:

	Assembly_store_layout (const Options & Opt, Hash * HashHead);

	// void	open_output_files (ofstream & , ofstream & , fstream & , fstream & );
	void	open_output_files (ofstream & , fstream & , fstream & );
	// void 	close_output_files (ofstream & , ofstream & , fstream & , fstream & );
	void	close_output_files (ofstream & , fstream & , fstream & );
	void	update_contig_layout (Contig * contig, Hash * HashHead);
	void	update_contig_layout_DEBUG (Contig * contig, Hash * HashHead, ofstream & File1,
				ofstream & File2, ofstream & File3);
	void	compute_layout_vector (Hash * , fstream & , fstream & );
	void	printTempLayoutBlock (bool, fstream & temp_layout_file);
	void	printLayoutVector (Layout_vector & L_vector, fstream & layout_file);
	void	printSubvectorsFromSublists (Layout_vector & , fstream & , fstream & );

	// ricava il vettore dei Contig_layout a partire dal file binario
	void 	get_contigs_layout (ofstream & File6, ofstream & File7,
				fstream & layout_file, vector <Contig_layout> & C_vector);

	// Print infos for layout check
	void	compute_layout_vector_DEBUG (Hash *, fstream &, fstream &, ofstream &,
				ofstream &, ofstream &, ofstream &);
	void	print_layout_sublists (ofstream & File5, fstream & temp_layout_file);
	void 	print_layout_from_subvectors (ofstream &, ofstream &, fstream &);

private:

	vector <Contig_layout> 	_contig_vector;
	unsigned int 			_contig_bits; // bits for layout representation
	Layout_list 			_layout_list;

	string					_temp_layout_name;
	string					_layout_name;
};

#endif /* ASSEMBLY_STORE_LAYOUT_H_ */
