#include "data_structures/Contig.h"

unsigned int		contig::Extension::_extThreshold = 2;
unsigned int		contig::Extension::_max_contig_length = 600;

unsigned long int	contig::Contig::q = 1;
unsigned long int	contig::Contig::h = 1;
unsigned int 		contig::Contig::_BL = 15;
unsigned int		contig::Contig::_overlap = 30;
unsigned int		contig::Contig::_extensionLength = 200;
unsigned int		contig::Contig::_seed_ins = 400;
unsigned int		contig::Contig::_seed_var = 200;
unsigned short		contig::Contig::_globalMismatch = 5;
bool				contig::Contig::_no_cycle = false;

namespace contig {

Extension::Extension (unsigned  int size)
{
	farthest = 0;
	_consensus = "";
	this->size = size;
	A = new unsigned int[size];
	memset(A, 0, (size)*sizeof(int));
	C = new unsigned int[size];
	memset(C, 0, (size)*sizeof(int));
	G = new unsigned  int[size];
	memset(G, 0, (size)*sizeof(int));
	T = new unsigned  int[size];
	memset(T, 0, (size)*sizeof(int));
	Total = new unsigned int[size];
	memset(Total, 0, (size)*sizeof(int));
}

Extension::~Extension()
{
	if (A != NULL)
		delete [] A;
	if (C != NULL)
		delete [] C;
	if (G != NULL)
		delete [] G;
	if (T != NULL)
		delete [] T;
	if (Total != NULL)
		delete [] Total;
}

void
Extension::set_static_vars (unsigned int extThreshold, unsigned int seed_ins, unsigned int seed_var)
{
	_extThreshold = extThreshold;
	_max_contig_length = seed_ins + seed_var;
}

void
Extension::clear_consensus ()
{
	for (unsigned int i = 0; i < size; i++) {
		A[i] = 0;
		C[i] = 0;
		G[i] = 0;
		T[i] = 0;
		Total[i] = 0; // altrimenti mi sballa tutto!!! TODO: CONTROLLARE!!!!
					  // inoltre in questo modo non rischio di estendere oltre extThr
	}
	farthest = 0; // devo ricalcolarlo quando ricavo il consenso dalle reads estratte
}

inline void
Extension::insertRead (unsigned int extStart, unsigned int index)
{
	extRead r;
	r.index = index;
	r.extPosition = extStart;
	_reads.push_back(r); // (index in HashValues, position in extension)
}

void
Extension::update_count (unsigned int extStart, const string & read)
{
	unsigned int j = extStart;
	for(unsigned int i = 0; i < read.length(); i++) {
		switch(read.at(i)) {
		case 'A': case 'a': A[j]++; Total[j]++; j++; break;
		case 'C': case 'c': C[j]++; Total[j]++; j++; break;
		case 'G': case 'g': G[j]++; Total[j]++; j++; break;
		case 'T': case 't': T[j]++; Total[j]++; j++; break;
		}
	}
	if (j > farthest) {
		farthest = j;
	}
}

void
Extension::computeConsensus (unsigned int numFingerprints, unsigned int &cons_start)
{
	_consensus = "";
	_low_rep.clear();
	unsigned int i = 0;
	while (Total[i] == 0 and i < farthest) {
		i++;
	}
	cons_start = i;
	while((Total[i] >= _extThreshold or i < numFingerprints) and i < farthest) {
		_consensus += returnMAX(A[i],C[i],G[i],T[i]);
		bool represented;
		if (checkRepeatStatus(A[i],C[i],G[i],T[i], represented)){
			pair <unsigned int, bool> el(i, represented);
			_low_rep.push_back(el);
		}
		i++;
	}
}

bool
Extension::checkRepeatStatus(int a, int c, int g, int t, bool& represented)
{
	float total = a + c + g + t;
	int max = a;
	int id = 0;
	if(c > max) {max = c; id = 1;}
	if(g > max) {max = g; id = 2;}
	if(t > max) {max = t; id = 3;}

	float res = max/total;
	if (res < perc_identity){ // distinguish low represented chars from those not represented at all
		represented = false;
		if (max == 1){
			return 1;
		}
	} else {
		represented = true;
	}

	// for every percentage < 1 and >= perc_identity a char is considered to be low represented
	if(total-max > 1) { // try to distinguish sequencing errors from SNPs
	//if (represented and total-max > 1){ // at least 2 chars different from consensus, otherwise it's a sequencing error
		int second_max = 0;
		if ((a != max or id != 0) and a > second_max) {second_max = a;}
		if ((c != max or id != 1) and c > second_max) {second_max = c;}
		if ((g != max or id != 2) and g > second_max) {second_max = g;}
		if ((t != max or id != 3) and t > second_max) {second_max = t;}
		if (second_max > 1){ // try to distinguish sequencing errors from SNPs
			return 1;
		}
	}
	return 0;
}

char
Extension::returnMAX (int a, int c, int g, int t)
{
	int max = a;
	if (c > max) max = c;
	if (g > max) max = g;
	if (t > max) max = t;

	if (max == a) return 'A';
	if (max == c) return 'C';
	if (max == g) return 'G';
	return 'T';
}

Extension::ExtReadType
Extension::SelectReadsForConsensus (unsigned int Extposition, string & sequence, unsigned int cons_start, unsigned int mismatches)
{
	if (low_rep().empty()) return LowRep;

	// indicates whether the read contains a caracter != consensus in a low represented position
	bool read_ok = true;
	vector <pair <unsigned int, bool> >::const_iterator l = low_rep().begin();
	unsigned int i = (Extposition < cons_start) ? (cons_start - Extposition) : 0;
	while (l != low_rep().end() and l->first < Extposition + i){
		l++;
	}

	// non mi serve controllare tutta la sequenza perché le posizioni con i mismatch
	// sono tutte e sole quelle low_rep
	while (l != low_rep().end() and
			l->first < Extposition + sequence.size() and
			l->first < cons_start + consensus().size()) {
		read_ok = (sequence.at(l->first-Extposition) == consensus().at(l->first-cons_start));
		if (read_ok and not l->second){ // the read is ok but covers a not-represented char
			string s = sequence.substr(0,l->first-Extposition);
			sequence = s;
			return NonRep;
		}
		l++;
	}

	/*
	unsigned short errors = 0;
	while (read_ok and i < sequence.size() and Extposition + i < cons_start + consensus().size()){
		errors += (sequence.at(i) != consensus().at(Extposition+i-cons_start));
		read_ok = (errors <= mismatches); // discard the current read if too many mismatches occur
		if (read_ok and l != low_rep().end() and Extposition+i == l->first){
			read_ok = (sequence.at(l->first-Extposition) == consensus().at(l->first-cons_start));
			if (read_ok and not l->second){ // the read is ok but covers a not-represented char
				string s = sequence.substr(0,l->first-Extposition);
				sequence = s;
				return NonRep;
			}
			l++;
		}
		i++;
	}
	*/

	return (read_ok) ? LowRep : Discarded;
}

Contig::Contig (Hash * HashHead, const string & seed, const string & seed_mate, unsigned int contigNum,
		bool seeds_as_query, bool reverse_seed, bool reverse_mate, bool verbose)
{
	_contig = seed;
	_mate = seed_mate;
	_contigNum = contigNum;
	_ContigStatus = Continue;
	// first position to be searched for overlapping reads
	_leftmost = 0;
	// first position to be searched for the mate
	_firstSearchPosition = (_seed_var < _seed_ins and _seed_ins - _seed_var > _mate.size()) ?
							_seed_ins - _seed_var - _mate.size() : 0;
	// compute fingerprint and position in hashvalue for seed and seed_mate
	if (seeds_as_query) {
		// seed_mate_id = 2*(_contigNum-1);
		unsigned long int fing = (not reverse_seed) ?
				HashHead->ComputeFingerprint (2*(_contigNum-1), true, 0) :
				HashHead->ComputeFingerprint (2*(_contigNum-1), false, _overlap - _BL);
		unsigned int l = (fing == q - 1) ? 2 * HashHead->numReads : HashHead->HASHcounter[fing + 1];
		_seed_hashvalue = HashHead->HASHcounter[fing];
		// while (_seed_hashvalue < l and (HashHead->HASHvalues.at(_seed_hashvalue).first != 2*(_contigNum-1))) {
		while (_seed_hashvalue < l and
				(HashHead->HASHvalues.at(_seed_hashvalue).id != 2*(_contigNum-1))) {
			_seed_hashvalue++;
		}
		// seed_mate_id = 2*_contigNum-1;
		fing = (not reverse_mate) ?
				HashHead->ComputeFingerprint (2*_contigNum-1, true, 0) :
				HashHead->ComputeFingerprint (2*_contigNum-1, false, _overlap - _BL);
		l = (fing == q - 1) ? 2 * HashHead->numReads : HashHead->HASHcounter[fing + 1];
		_seed_mate_hashvalue = HashHead->HASHcounter[fing];
		//while (_seed_mate_hashvalue < l and (HashHead->HASHvalues.at(_seed_mate_hashvalue).first != 2*_contigNum-1)) {
		while (_seed_mate_hashvalue < l and
				(HashHead->HASHvalues.at(_seed_mate_hashvalue).id != 2*_contigNum-1)) {
			_seed_mate_hashvalue++;
		}
	}
	if (verbose) {
		cout << endl;
		cout << "seed_string: " <<  endl;
		cout << seed << endl;
		cout << "seed_mate_string: " << endl;
		cout << seed_mate << endl;
		cout << "seed_mate_hashvalue: " << _seed_mate_hashvalue << endl;
		cout << endl;
		cout << "**** CONTIG " << _contigNum << " ****" << endl;
	}
}

void
Contig::set_static_vars (unsigned int qq, unsigned int hh, unsigned int BL, unsigned int overlap,
						 unsigned int max_read_length, unsigned int seed_ins, unsigned int seed_var,
						 unsigned int globalMismatch, bool no_cycle)
{
	q = qq;
	h = hh;
	_BL = BL;
	_overlap = overlap;
	_extensionLength = 2*max_read_length;
	_seed_ins = seed_ins;
	_seed_var = seed_var;
	_globalMismatch = globalMismatch;
	_no_cycle = no_cycle;
}

bool
Contig::is_read_used (unsigned int index) // index in HashValues
{
	vector <contigRead>::iterator it = _reads.begin();
	while (it != _reads.end() and it->index != index) {
		it++;
	}
	return (it != _reads.end());
}

void
Contig::InitializeExtension()
{
	contigExtension = new Extension(_extensionLength);
	_firstExtendedPosition = 0;
}

void
Contig::DeleteExtension()
{
	if(contigExtension != NULL) {
		delete contigExtension;
	}
	_firstExtendedPosition = 0;
}

unsigned long int
Contig::ComputeFingerprint (int n)
{
	unsigned int start = _contig.length() - _overlap - n;
	unsigned long int fingerprint = 0;
	for (unsigned int j = start; j < start + _BL; j++) {
		char c = 0;
		switch (_contig.at(j)) {
		case 'a' : case 'A' : c = 0; break;
		case 'c' : case 'C' : c = 1; break;
		case 'g' : case 'G' : c = 2; break;
		case 't' : case 'T' : c = 3; break;
		}
		fingerprint = ((fingerprint << 2) + c) % q;
	}
	return fingerprint;

}

unsigned long int
Contig::ComputeFingerprintForCheck (unsigned int start)
{
	unsigned long int fingerprint = 0;
	//unsigned int end = (start + _BL < _contig.size()) ? start + _BL : _contig.size();
	unsigned int end = start + _BL;
	for (unsigned int j = start; j < end; j++) {
		char c = 0;
		switch (_contig.at(j)) {
		case 'a' : case 'A' : c = 0; break;
		case 'c' : case 'C' : c = 1; break;
		case 'g' : case 'G' : c = 2; break;
		case 't' : case 'T' : c = 3; break;
		}
		fingerprint = ((fingerprint << 2) + c) % q;
	}
	return fingerprint;
}

unsigned long int
Contig::ComputeFingerprintOneStep (unsigned int start, unsigned int fingerprint)
{
	char first_digit = 0, c = 0;
	switch (_contig.at(start)) {
	case 'a' : case 'A' : first_digit = 0; break;
	case 'c' : case 'C' : first_digit = 1; break;
	case 'g' : case 'G' : first_digit = 2; break;
	case 't' : case 'T' : first_digit = 3; break;
	}
	switch (_contig.at(start+_BL)) {
	case 'a' : case 'A' : c = 0; break;
	case 'c' : case 'C' : c = 1; break;
	case 'g' : case 'G' : c = 2; break;
	case 't' : case 'T' : c = 3; break;
	}
	return ((q + fingerprint - first_digit * h) * 4 + c) % q;
}

bool
Contig::MaybeInsertRead (const string & read, unsigned int start, unsigned int extStart,
		unsigned int localMismatches, unsigned int index)
{
	unsigned int mismatches = 0;
	unsigned int el = 0;
	// controllo di TUTTO l'overlap
	unsigned int end = (_contig.length() > read.length() + start) ? (start + read.length()) : _contig.length();
	for (unsigned int st = start; st < end; st++) {
		mismatches += (_contig.at(st) != read.at(el++));
	}
	if (mismatches <= localMismatches) {
		// Controllo readCycle solo al momento dell'inserimento definitivo
		contigExtension->insertRead (extStart, index); // update temporary reads
		contigExtension->update_count (extStart, read); // update temporary consensus
		return true;
	}
	return false;
}

void
Contig::InsertRead (const string & read, unsigned int start, unsigned int extStart, unsigned int index)
{
	if (not _no_cycle or not is_read_used(index)) { // controlla se la read fa già parte di quelle definitive inserite
		contigExtension->update_count (extStart, read);
		contigRead r;
		r.index = index;
		r.contPosition = start;
		_reads.push_back (r);
	} else {
		_ContigStatus = ReadCycleFound;
	}
	return;
}

void
Contig::MateSearch (Hash * HashHead, bool seeds_as_query, bool verbose)
{
	if (not (length() >= _firstSearchPosition + _mate.size() and
			 length() + _seed_var >= _seed_ins and
			 length() <= _seed_ins + _seed_var) ) {
		// return false; // contig is too short
		return;
	}
	bool mateFound = false; // true if I found the mate either in hash or as a sequence
	unsigned int mate_pos;
	if (seeds_as_query and not _reads.empty()) { // if some extension have been made then search hashvalue first
		vector <contigRead>::reverse_iterator r_iter = _reads.rbegin();
		while (not mateFound and r_iter != _reads.rend() and // comincia a controllare dalla fine
			   r_iter->contPosition >= _firstSearchPosition) {
			if (r_iter->index == _seed_mate_hashvalue and
				r_iter->contPosition + _mate.size() <= length()) { // la mate non deve sforare dal contig!!
				mateFound = true;
			}
			mate_pos = r_iter->contPosition;
			r_iter++;
		}
	}
	if (mateFound) { // the mate has been used to extend
		_ContigStatus = MatePairFound;
		// cut the contig sequence after the mate
		_contig.resize(mate_pos + _mate.size());
		// drop the reads ending after the mate
		int ii = _reads.size() - 1;
		while (ii >= 0) {
			unsigned int r_start = _reads.at(ii).contPosition;
			//size_t r_size = HashHead->readsMulti.at(HashHead->HASHvalues.at(_reads.at(ii).index).first).length();
			uint32_t id = HashHead->HASHvalues.at(_reads.at(ii).index).id;
			size_t r_size = HashHead->readsMulti.at(id).length();
			if (r_start + r_size > mate_pos + _mate.size()) {
				// erase read
				_reads.erase(_reads.begin() + ii);
			}
			ii--;
		}
		if (verbose) {
			cout << "MATE FOUND!!!!" << endl;
			cout << "current contig: " << endl;
			cout << _contig << endl;
			cout << "contig length: " << length() << endl;
		}
		// return true;
	} else { // Mate not found: search pos-by-pos
		mate_pos = _firstSearchPosition;
		while (not mateFound and mate_pos + _mate.size() <= length()) { // non deve sforare!!
			unsigned int mismatches = 0;
			for (unsigned int el = 0; el < _mate.size() ; el++) {
				if(_contig.at(mate_pos+el) != _mate.at(el)) {
					mismatches++;
				}
			}
			mateFound = (mismatches * 100 <= _globalMismatch * _mate.length());
			mate_pos += (not mateFound);
		}
		if (mateFound) {
			_ContigStatus = MatePairFound;
			// cut the contig sequence after the mate
			_contig.resize(mate_pos + _mate.size());
			// drop the reads ending after the mate
			int insert_position = _reads.size();
			if (not _reads.empty()) {
				int ii = _reads.size() - 1;
				while (ii >= 0) {
					unsigned int r_start = _reads.at(ii).contPosition;
					//size_t r_size = HashHead->readsMulti.at(HashHead->HASHvalues.at(_reads.at(ii).index).first).length();
					uint32_t id = HashHead->HASHvalues.at(_reads.at(ii).index).id;
					size_t r_size = HashHead->readsMulti.at(id).length();
					if (r_start + r_size > mate_pos + _mate.size()) {
						// erase read
						_reads.erase(_reads.begin() + ii);
						insert_position--; // if I have to insert after a dropped read
					} else {
						// find the position in which the seed mate should be inserted
						if (seeds_as_query and r_start >= mate_pos and
							// not HashHead->HASHvalues.at(_reads.at(ii).index).second and
							(r_start > mate_pos or _reads.at(ii).index > _seed_mate_hashvalue) ) {
							//insert_position = ii+1;
							insert_position = ii; // in questo modo inserisce PRIMA di ii
						}
					}
					ii--;
				}
			}
			// add the seed to the list of reads
			if (seeds_as_query) {
				if (_reads.empty()) {
					// add the seed because no extension has been made yet
					contigRead c;
					c.index = _seed_hashvalue;
					c.contPosition = 0;
					_reads.push_back (c); // add the seed to the reads used for extension
					insert_position++;
				}
				contigRead c;
				c.index = _seed_mate_hashvalue;
				c.contPosition = mate_pos;
				if (insert_position < (int) _reads.size()) {
					_reads.insert (_reads.begin() + insert_position, c);
				} else {
					_reads.push_back (c);
				}
			}
			// return true;
		}
		// return false;
	}
	if (mateFound and verbose) {
		for (vector <contigRead>::iterator iter = _reads.begin(); iter != _reads.end(); iter++) {
			cout << "(" << iter->index << "," << iter->contPosition << ")" << endl;
		}
	}
}

void
Contig::extension (unsigned int first_pos)
{
	if(_ContigStatus != ReadCycleFound) { // if I realize that a cycle as been found don't perform extension at all
		unsigned int cons_start;
		// vector< pair<unsigned int, bool> > low_rep;
		//string extension = contigExtension->returnConsensus (length() - first_pos, low_rep, cons_start);
		contigExtension->computeConsensus (length() - first_pos, cons_start);
		_firstExtendedPosition = first_pos + cons_start;
		// compute extension and number of errors of such extension
		if (_firstExtendedPosition + contigExtension->consensus().length() > _contig.length()) {
			_contig = _contig.substr(0, _firstExtendedPosition);
			_contig.append(contigExtension->consensus());
			// length cutoff to avoid loops
			_ContigStatus = (_contig.length() < contigExtension->max_contig_length()) ? Continue : LengthExceed;
		} else { // if number of errors is higher than threshold stop contig extension
			_ContigStatus = RepetitiveSequence;
		}
	}
}

void
Contig::ComputeTemporaryReads (Hash* HashHead, bool verbose)
{
	unsigned int Extposition = 0; // position of an overlapping read in the current extension
	unsigned long int for_fing = 0, rev_fing = 0;
	if (verbose) {
		cout << "TEMPORARY READS" << endl;
	}
	for (unsigned int pos = _leftmost; pos < length() - _overlap + 1; pos++) {
		unsigned int localMismatch = (unsigned int) _globalMismatch * (length() - pos - _BL) / 100;
		// READS FORWARD
		for_fing = (pos == _leftmost) ? ComputeFingerprintForCheck (pos)
									  : ComputeFingerprintOneStep (pos - 1, for_fing);
		unsigned int p = HashHead->HASHcounter[for_fing];
		unsigned int l = (for_fing == HashHead->q - 1) ? 2 * HashHead->numReads : HashHead->HASHcounter[for_fing + 1];
		for (unsigned int j = p; j < l; j++) {
			// NUOVO: cerca il seed ed evita di aggiungerlo prima
			// (in questo modo: 1) sono sicura di inserirlo nella posizione corretta
			// 					2) sono sicura di trovarlo perché controllo tutte le posizioni senza usare lo slack)
			//if (HashHead->HASHvalues.at(j).second) { // solo reads +
			if (HashHead->HASHvalues.at(j).strand) { // solo reads +
				//unsigned int index = HashHead->HASHvalues.at(j).first;
				unsigned int id = HashHead->HASHvalues.at(j).id;
				Reads f = HashHead->readsMulti.at(id);
				if (f.length() >= _overlap) { // exclude short reads
					//string electedSequence = f.toString(HashHead->HASHvalues.at(j).second);
					string electedSequence = f.toString(HashHead->HASHvalues.at(j).strand);
					if (MaybeInsertRead (electedSequence, pos, Extposition, localMismatch, j) and verbose) {
						for (unsigned int kk = 0; kk < pos; kk++) {
							cout << "*";
						}
						cout << electedSequence << endl;
					}
				}
			}
		}
		// READS REVERSE
		rev_fing = (pos == _leftmost) ? ComputeFingerprintForCheck (pos + _overlap - _BL)
									  : ComputeFingerprintOneStep (pos + _overlap - _BL - 1, rev_fing);
		p = HashHead->HASHcounter[rev_fing];
		l = (rev_fing == HashHead->q - 1) ? 2 * HashHead->numReads : HashHead->HASHcounter[rev_fing + 1];
		for (unsigned int j = p; j < l; j++) {
			//if (not HashHead->HASHvalues.at(j).second) { // solo reads -
			if (not HashHead->HASHvalues.at(j).strand) { // solo reads -
				//unsigned int index = HashHead->HASHvalues.at(j).first;
				unsigned int id = HashHead->HASHvalues.at(j).id;
				Reads f = HashHead->readsMulti.at(id);
				if (f.length() >= _overlap) {
					//string electedSequence = f.toString(HashHead->HASHvalues.at(j).second);
					string electedSequence = f.toString(HashHead->HASHvalues.at(j).strand);
					if (MaybeInsertRead (electedSequence, pos, Extposition, localMismatch, j) and verbose) {
						for (unsigned int kk = 0; kk < pos; kk++) {
							cout << "*";
						}
						cout << electedSequence << endl;
					}
				}
			}
		}
		Extposition++;
	}
}

void
Contig::ComputeExtensionReads (Hash* HashHead, bool verbose)
{
	// TRUE EXTENSION
	unsigned int cons_start;
	contigExtension->computeConsensus (length() - _leftmost, cons_start);
	if (verbose) {
		printConsensus (length() - _leftmost, _leftmost, cons_start);
	}
	if (_leftmost + cons_start + contigExtension->consensus().length() <= length()) {
		_ContigStatus = NoMoreExtensions;
	} else {
		if (verbose) {
			cout << "EXTRACTED READS" << endl;
		}
		contigExtension->clear_consensus();
		unsigned int first_pos = _leftmost;
		for (vector <extRead>::const_iterator	r =  contigExtension->reads().begin();
											 	r != contigExtension->reads().end();
											 	r++) {
			//Reads f = HashHead->readsMulti.at(HashHead->HASHvalues.at(r->index).first);
			Reads f = HashHead->readsMulti.at(HashHead->HASHvalues.at(r->index).id);
			//string sequence = f.toString(HashHead->HASHvalues.at(r->index).second);
			string sequence = f.toString(HashHead->HASHvalues.at(r->index).strand);
			Extension::ExtReadType type = Extension::LowRep; // read is ok
			// controllo di non scartare il seed
			if (r->index != _seed_hashvalue) {
				unsigned int mismatches = _globalMismatch * (sequence.size() - _BL) / 100;
				type = contigExtension->SelectReadsForConsensus(r->extPosition, sequence, cons_start, mismatches);
			}
			if (type != Extension::Discarded) {
				if (type == Extension::NonRep) {
					// add current index to the set of reads to be checked at the end
					_to_be_checked.push_back(_reads.size());
				}
				unsigned int contigPosition = first_pos + r->extPosition;
				InsertRead (sequence, contigPosition, r->extPosition, r->index);
				set_leftmost (contigPosition+1); // update leftmost for subsequent extension
				if (verbose) {
					vector <contigRead>::reverse_iterator r_iter = _reads.rbegin();
					cout << "index = " << r_iter->index << ", contigPosition = " << r_iter->contPosition << endl;
					for (unsigned int kk = 0; kk < r_iter->contPosition; kk++) {
						cout << "*";
					}
					cout << sequence << endl;
				}
			}
		}
		extension (first_pos);
		if (verbose) {
			printConsensus (length() - first_pos, first_pos, cons_start);
			cout << "Current contig:" << endl;
			cout << _contig << endl;
		}
	}
}

void
Contig::checkNonRepType (Hash * HashHead)
{
	if (_ContigStatus != MatePairFound or _to_be_checked.empty()) {
		// don't care about layout if the mate wasn't found
		return;
	}
	// parto dalla fine così non ho problemi negli indici di _reads
	vector <unsigned int>::reverse_iterator t_iter = _to_be_checked.rbegin();
	int i = _reads.size() - 1;
	while (t_iter != _to_be_checked.rend() and i > -1) {
		if ((unsigned int) i == *t_iter) { // ho trovato una read che è stata trimmata
			//unsigned int index = HashHead->HASHvalues.at(_reads.at(i).index).first;
			unsigned int index = HashHead->HASHvalues.at(_reads.at(i).index).id;
			//bool strand = HashHead->HASHvalues.at(_reads.at(i).index).second;
			bool strand = HashHead->HASHvalues.at(_reads.at(i).index).strand;
			string sequence = HashHead->readsMulti.at(index).toString(strand);
			unsigned short local_mismatch = (unsigned short) _globalMismatch * sequence.size() / 100;
			unsigned short errors = 0; // non serve fare il resize della read
									   // perché se è nel layout non termina dopo il contig
			unsigned int j = 0;
			while (errors < local_mismatch + 1 and j < sequence.size()) {
				errors += (sequence.at(j) != _contig.at(_reads.at(i).contPosition + j));
				j++;
			}
			if (errors > local_mismatch) {
				// elimina la read
				_reads.erase(_reads.begin() + i);
			}
			t_iter++;
		}
		i--;
	}
}

void
Contig::printConsensus (unsigned int numfing, unsigned int first_pos, unsigned int cons_start)
{
	// if perc < perc_identity, the char is not represented at all
	cout << "CONSENSUS SEQUENCE" << endl;
	if (not contigExtension->consensus().empty()) {
		for (unsigned int s = 0; s < first_pos + cons_start; s++) {
			cout << "*";
		}
		cout << contigExtension->consensus() << endl;
		unsigned int j = 0;
		unsigned int s = 0;
		while (s < first_pos + cons_start + contigExtension->consensus().size()) {
			if (s < first_pos + cons_start) {
				cout << "*";
			} else if (s >= first_pos + cons_start and j < contigExtension->low_rep().size() and s == first_pos + contigExtension->low_rep().at(j).first) {
				cout << "*";
				j++;
			} else {
				cout << "X";
			}
			s++;
		}
		cout << endl;
		cout << "number of low represented characters: " << contigExtension->low_rep().size() << endl;
	}
}

void
Contig::printContig(ofstream& ostrFile)
{
	ostrFile << ">contig" << _contigNum;
	switch (_ContigStatus) {
	case MatePairFound :		ostrFile << "_pairFound" << endl; 		break;
	case LengthExceed :			ostrFile << "_LengthExceed" << endl;	break;
	case RepetitiveSequence :	ostrFile << "_repSeq" << endl;			break;
	case NoMoreExtensions :		ostrFile << "_noMoreExt" << endl;		break;
	case ReadCycleFound :		ostrFile << "_ReadCycle" << endl;		break;
	default :					cout << "ERROR in contig " << _contigNum << endl;
								ostrFile << "_erroneusStatus" << endl;	break;
	}
	ostrFile << _contig << endl;
}

void
Contig::addContig (vector <string> & fasta_vector, unsigned int & size)
{
	string fasta (">contig");
	fasta.append ( static_cast<ostringstream*>( &(ostringstream() << _contigNum) )->str());
	switch (_ContigStatus) {
	case MatePairFound :		fasta.append ("_pairFound\n");		break;
	case LengthExceed :			fasta.append ("_LengthExceed\n");	break;
	case RepetitiveSequence : 	fasta.append ("_repSeq\n");			break;
	case NoMoreExtensions : 	fasta.append ("_noMoreExt\n");		break;
	case ReadCycleFound : 		fasta.append ("_ReadCycle\n");		break;
	default : 					cout << "ERROR in contig " << _contigNum << endl;
								fasta.append ("_erroneusStatus\n");	break;
	}
	fasta.append (_contig);
	fasta.append ("\n");
	size += fasta.size();
	fasta_vector.push_back (fasta);
}

// ONLY FOR DEBUG
void
Contig::printLayout (Hash * HashHead, ofstream & out)
{
	out << _contigNum << " ";
	for (vector <contigRead>::const_iterator c_iter = _reads.begin();
											 c_iter != _reads.end();
											 c_iter++) {
		readOriented r = HashHead->HASHvalues.at(c_iter->index);
		/*
		out << "(" << r.first << "," << c_iter->contPosition << ","
			<< HashHead->readsMulti.at(r.first).length()
			<< "," << r.second << ") ";
			*/
		out << "(" << r.id << "," << c_iter->contPosition << ","
			<< HashHead->readsMulti.at(r.id).length()
			<< "," << r.strand << ") ";
	}
	out << endl;
}


}  /* end of namespace contig */
