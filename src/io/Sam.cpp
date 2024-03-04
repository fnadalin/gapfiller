#include "Sam.h"

Sam::Sam () {
	_is_open = false;
	_b = new bam1_t;
	_b->data = (uint8_t*)malloc(1024);
}

Sam::Sam (const string & s) {

	unsigned short dot = s.find_last_of(".");
	string ext = s.substr(dot+1);
	_is_sam = (ext.compare("sam") == 0);

	if (_is_sam){
	    _filesam = samopen(s.c_str(), "r", &(_h));
	} else if (ext == "bam") {
	    _filesam = samopen(s.c_str(), "rb", &(_h));
	}

	_is_open = (_filesam != NULL);
	if (_is_open) {
		_b = new bam1_t;
		_b->data = (uint8_t*)malloc(1024);
	}
}


Sam::~Sam () {
	if (_is_open) samclose(_filesam);
	if (_b != NULL) {
		free(_b->data);
		delete _b;
	}
}


void Sam::open (const string & s) {

	unsigned short dot = s.find_last_of(".");
	if (dot >= s.size()) {
		return;
	}
	string ext = s.substr(dot+1);
	_is_sam = (ext.compare("sam") == 0);

	if (_is_sam){
		_filesam = samopen(s.c_str(), "r", &(_h));
	} else if (ext == "bam") {
		_filesam = samopen(s.c_str(), "rb", &(_h));
	}

	_is_open = (_filesam != NULL);
}


bool Sam::read () {
	return (samread (_filesam, _b) > 0);
}


string Sam::id () {
	string id (bam1_qname(_b));
	return id;
}


string Sam::sequence (bool strand) {

	if (_b == NULL) {
		return string();
	}
	uint8_t * seq =  bam1_seq(_b);
	int32_t l = _b->core.l_qseq;
	string ts;
	ts.reserve(l+1);
	for (int32_t i = 0; i < l; i++) {
		switch (bam1_seqi(seq, i)) {
		case 1:
			ts.push_back('A');
			break;
		case 2:
			ts.push_back('C');
			break;
		case 4:
			ts.push_back('G');
			break;
		case 8:
			ts.push_back('T');
			break;
		default:
			ts.push_back('N');
			break;
		}
	}
	if (strand) {
		return ts;
	} else {
		return reverse_complement_standalone_str (ts);
	}
}


bool Sam::is_first () {
	return (_b == NULL) ? false : _b->core.flag & BAM_FREAD1;
}


bool Sam::is_reverse () {
		return (_b == NULL) ? false : _b->core.flag & BAM_FREVERSE;
}

/*
bool Sam::is_strand_ok (bool pair_end) {
	bool innie = (is_first () and not is_reverse ()) or (not is_first() and is_reverse());
	return pair_end ? innie : not innie;
}
*/

