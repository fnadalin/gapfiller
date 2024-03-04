#ifndef CONTIG_H_
#define CONTIG_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <stdint.h>
using namespace std;

#include "data_structures/Hash.h"
#include "options/Singleton.h"

#define perc_identity 0.6

typedef pair <string, vector<uint64_t> > contig_layout;

namespace contig {

typedef struct {
	unsigned int index;
	unsigned int extPosition;
} extRead;

typedef struct {
	unsigned int index;
	unsigned int contPosition;
} contigRead;

class Extension {

public:

	enum ExtReadType {LowRep, NonRep, Discarded};

	Extension (unsigned  int size);
	~Extension ();

	static void							set_static_vars (unsigned int extThreshold, unsigned int seed_ins, unsigned int seed_var);

	inline vector <extRead> const &		reads () { return _reads; }
	inline string const &				consensus () { return _consensus; }
	inline vector <pair<unsigned int, bool> > const & low_rep() { return _low_rep; } // <position, represented>
																					 // represented = true : low rep, represented = false : non rep
	inline unsigned int					get_farthest () { return farthest; }
	inline void							set_farthest (unsigned int f) { farthest = f; }
	inline unsigned int					max_contig_length () { return _max_contig_length; }

	void								clear_consensus ();
	inline void							insertRead (unsigned int extPosition, unsigned int index);
	void								update_count (unsigned int extPosition, const string & read);

	void								computeConsensus (unsigned int numFingerprints, unsigned int& cons_start);
	ExtReadType							SelectReadsForConsensus (unsigned int extPosition, string & sequence, unsigned int cons_start, unsigned int mismatches);

protected:

	static unsigned int					_extThreshold;
	static unsigned int					_max_contig_length;

	unsigned int *						A;
	unsigned int *						C;
	unsigned int *						G;
	unsigned int *						T;
	unsigned int *						Total;
	unsigned int 						size;
	unsigned int						farthest;
	vector <extRead>					_reads;
	vector < pair<unsigned int, bool> > _low_rep; // position and status (low rep/non rep) of low represented positions
	string								_consensus;

	bool 								checkRepeatStatus (int a, int c, int g, int t, bool& represented);
	char								returnMAX (int, int, int, int);
};

class Contig {

public:

	enum ExtensionStatus {Continue, ReadCycleFound, NoMoreExtensions, RepetitiveSequence, MatePairFound, LengthExceed};

	//Contig (const string & seed, const string & seed_mate, unsigned int contigNum);
	Contig (Hash * HashHead, const string & seed, const string & seed_mate, unsigned int contigNum, bool seeds_as_query,
			bool reverse_seed, bool reverse_mate, bool verbose);
	static void							set_static_vars (unsigned int qq, unsigned int hh, unsigned int BL, unsigned int overlap,
													 	 unsigned int max_read_length, unsigned int seed_ins, unsigned int seed_var,
													 	 unsigned int globalMismatch, bool no_cycle);

	inline unsigned int					contigNum () { return _contigNum; }
	inline string const &				sequence () { return _contig; }
	inline unsigned int 				length() { return _contig.length(); }
	inline ExtensionStatus 				returnStatus() { return _ContigStatus; }
	inline unsigned int					BL() { return _BL; }
	static unsigned int 				overlap() { return _overlap; }
	inline unsigned int 				extension_length() { return _extensionLength; }
	inline unsigned int					first_extended_position () { return _firstExtendedPosition; }
	inline unsigned int					leftmost () { return _leftmost; }
	inline vector <contigRead> const &	reads () { return _reads; }

	bool 								is_read_used (unsigned int index); // index in HashValues
	inline void 						resize (unsigned int new_length) { if (_contig.size() > new_length) _contig.resize(new_length); }
	inline void 						setStatus (ExtensionStatus status) { _ContigStatus = status; }
	inline void 						set_leftmost (unsigned int l) { _leftmost = l; }

	void								MateSearch (Hash * HashHead, bool seeds_as_query, bool verbose);

	void								InitializeExtension();
	void								DeleteExtension();

	unsigned long int					ComputeFingerprint (int n);
	unsigned long int					ComputeFingerprintForCheck (unsigned int start);
	unsigned long int					ComputeFingerprintOneStep (unsigned int start, unsigned int fingerprint);

	bool 								MaybeInsertRead(const string & read, unsigned int start, unsigned int extStart, unsigned int localMismatches, unsigned int index);
	void 								InsertRead(const string & read, unsigned int start, unsigned int extStart, unsigned int index);
	void 								extension (unsigned int first_pos);

	void								ComputeTemporaryReads (Hash* HashHead, bool verbose);
	void								ComputeExtensionReads (Hash* HashHead, bool verbose);
	void 								checkNonRepType (Hash * HashHead);

	void								printConsensus (unsigned int numfing, unsigned int first_pos, unsigned int cons_start);
	void								printContig (ofstream & ostrFile);
	void								addContig (vector <string> & fasta_vector, unsigned int & size);
	void								printLayout (Hash* HashHead, ofstream & out); // ONLY FOR DEBUG

	Extension *							contigExtension;

private:

	vector <contigRead>					_reads;  // successfully inserted reads
	vector <unsigned int>				_to_be_checked; // reads cut in correspondence of a non-rep position
													   // => check whether they actually match the final sequence
	string 								_contig;
	string								_mate;

	static unsigned long int 			q;
	static unsigned long int 			h;
	static unsigned int	 				_BL;
	static unsigned int					_overlap;
	static unsigned int 				_extensionLength;
	static unsigned int					_seed_ins;
	static unsigned int					_seed_var;
	static unsigned short				_globalMismatch;
	static bool							_no_cycle;

	unsigned int						_leftmost; // position from which overlapping reads are considered
	unsigned int 						_firstExtendedPosition; // starting position of the first read found (global coordinate)
	unsigned int						_seed_hashvalue;
	unsigned int						_seed_mate_hashvalue;
	unsigned int						_firstSearchPosition; // starting position to search for the mate
	ExtensionStatus						_ContigStatus;
	unsigned int 						_contigNum;

};

}

#endif /* CONTIG_H_ */

