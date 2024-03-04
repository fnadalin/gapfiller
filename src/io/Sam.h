#ifndef SAM_H
#define SAM_H

#include "../samtools/sam.h"
#include "../common.h"

#include <string>
using namespace std;

class Sam {

	public:

	Sam ();
	Sam (const string & s);
	~Sam();

	void open (const string & s);
	bool is_open () { return _is_open; }
	bool is_sam () { return _is_sam; }
	bool read ();
	string id ();
	string sequence (bool strand);
	bool is_first ();
	bool is_reverse ();
	// bool is_strand_ok (bool pair_end); // reverse complement if the sequence is not compatible with a paired read

	private:

	samfile_t * _filesam;
	bam_header_t _h;
	bam1_t * _b;
	unsigned _is_open: 1, _is_sam: 1;

};

#endif 
