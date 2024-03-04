#include "options/Singleton.h"

Singleton * Singleton::_instance = NULL;

Singleton::Singleton() { }

Singleton *
Singleton::getInstance(void)
{
    if (_instance)
        return _instance;
    else
        return ( _instance = new Singleton() );
}

void
Singleton::printInstance(void)
{
	cout << endl;
	cout << "Bit size of layout fields:" << endl;
	cout << "contig_id = " <<_CONTIG_ID_MAX_SIZE << endl;
    cout << "number of sublists = " << _NUM_SUBLISTS_MAX_SIZE << endl;
    cout << "read ID = " << _READ_ID_MAX_SIZE << endl;
    cout << "distance between reads = " << _READ_DIST_MAX_SIZE << endl;
    cout << "read length = " << _READ_LEN_MAX_SIZE << endl;
}

unsigned int
Singleton::k(void)
{
	return _k;
}

unsigned long int
Singleton::q(void)
{
	return _q;
}

unsigned long int
Singleton::h(void)
{
	return _h;
}

unsigned short
Singleton::BL(void)
{
	return _BL;
}

unsigned short
Singleton::overlap(void)
{
	return _overlap;
}

unsigned short
Singleton::extensionLength(void)
{
	return _extensionLength;
}

void
Singleton::set_global_variables (unsigned int k, unsigned int b, unsigned int overlap, unsigned int max_read_length)
{
	_k = k;
	_q = 0;
	for (unsigned int i = 0; i < 2*k; i++) {
		_q |= 1 << i;
	}
	_BL = b;
	_h = 1;
	while (b > 1) {
		_h = (4*_h)%_q;
		b--;
	}
	_extensionLength = 2 * max_read_length;
}

void
Singleton::set_layout_variables (unsigned int numReads, unsigned int max_read_length)
{
	_READ_LEN_MAX_SIZE = 0; // 8
	while (pow(2,_READ_LEN_MAX_SIZE)-1 < max_read_length) {
		_READ_LEN_MAX_SIZE++;
	}
	_READ_ID_MAX_SIZE = 0;
	while (pow(2,_READ_ID_MAX_SIZE) < numReads) {
		_READ_ID_MAX_SIZE++;
	}
	_CONTIG_ID_MAX_SIZE = _READ_ID_MAX_SIZE;
	_READ_INDEX_MAX_SIZE = 1 + _READ_ID_MAX_SIZE;

	_NUM_SUBLISTS_MAX_SIZE = 8;

	_READ_DIST_MAX_SIZE = 0;
	while (pow(2,_READ_DIST_MAX_SIZE) < INFTY) {
		_READ_DIST_MAX_SIZE++;
	}
}

void
Singleton::set_layout_variables_for_vector (unsigned int read_id, unsigned int len)
{
	_READ_ID_MAX_SIZE = read_id;
	_READ_INDEX_MAX_SIZE = 1 + _READ_ID_MAX_SIZE;
	_READ_DIST_MAX_SIZE = 0;
	while (pow(2,_READ_DIST_MAX_SIZE) < INFTY) {
		_READ_DIST_MAX_SIZE++;
	}
	_READ_LEN_MAX_SIZE = len;
	_CONTIG_ID_MAX_SIZE = _READ_ID_MAX_SIZE;
	_NUM_SUBLISTS_MAX_SIZE = 8;
}

bool
Singleton::feasibility_check ()
{
	if (_CONTIG_ID_MAX_SIZE > 32 or _READ_ID_MAX_SIZE > 31 or
		_READ_DIST_MAX_SIZE > 16 or _READ_LEN_MAX_SIZE > 14) {
		// read id is stored in 31 bits in HashValues
		// read length is stored in 14 bits in Reads
		// read dist uses 16 bits in read_layout_list
		return false;
	}
	return true;
}
