#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_

#include "options/Options.h"
#include "data_structures/Contig.h"
#include "data_structures/Layout_list.h"
#include "io/Auto_Zip.h"
#include "io/Binary.h"
using namespace contig;
using namespace layout;

class Assembly {

public:

	Assembly (const Options & Opt, Hash * HashHead); // here store common options
	virtual ~Assembly () { } // do nothing

	void	printFastaBlock (bool END_OF_INPUT, ostream & _fasta_file);
	void	printTrashBlock (bool END_OF_INPUT, ostream & _trash_file);
	void	printStatsHeader (ofstream & stats_file);
	void	printStatsLine (ofstream & stats_file);
	virtual void printTempLayoutBlock (bool, fstream & temp_layout_file) = 0;

	virtual void open_output_files (ofstream & , fstream & , fstream & ) = 0;
	virtual void close_output_files (ofstream & , fstream & , fstream & ) = 0;
	void	execute (Hash * HashHead);
	bool	new_contig (const string &, const string &, unsigned int, Hash *, bool);
	virtual void update_contig_layout (Contig * contig, Hash * HashHead) = 0;
	virtual void compute_layout_vector (Hash *, fstream &, fstream &) = 0;

	// Print info for layout check
	void	execute_DEBUG (Hash * HashHead);
	bool	new_contig_DEBUG (const string &, const string &, unsigned int, Hash *,
						bool, ofstream &, ofstream &, ofstream &);
	virtual void update_contig_layout_DEBUG (Contig *, Hash *, ofstream &, ofstream &, ofstream &) = 0;
	virtual void compute_layout_vector_DEBUG (Hash *, fstream &, fstream &, ofstream &,
												ofstream &, ofstream &, ofstream &) = 0;
protected:

	vector <string> 	_fasta_vector;
	unsigned int 		_fasta_vector_size;
	vector <string>		_trash_vector;
	unsigned int 		_trash_vector_size;
	bool 				_gzip;
	bool				_bzip2;

	bool				_seed_fasta;
	bool				_seed_sam;
	bool				_limitedExtension;
	unsigned int		_maxReadsToExtend;
	vector <string> 	_seed_fw_files;
	vector <string> 	_seed_rv_files; 
	vector <string> 	_seed_sam_files;
	bool				_reverse_seed;
	bool				_reverse_mate;
	bool				_verbose;

	unsigned long int	_numberOfExtendedReads;
	unsigned long int	_numberOfExtensions;
	unsigned long int	_stopNotPossibleExtensions;
	unsigned long int	_stopRepetitiveSequence;
	unsigned long int	_stopLengthMaximumExceed;
	unsigned long int	_stopMatePairFound;
	unsigned long int	_stopReadCycleFound;

	string				_fasta_name;
	string				_trash_name;
	string				_stats_name;
};

#endif /* ASSEMBLY_H_ */
