#pragma once

#include <iostream>
#include <string>
#include <cmath>
using namespace std;

#define INFTY 10000

class Singleton
{
    unsigned int			_CONTIG_ID_MAX_SIZE;
    unsigned int			_NUM_SUBLISTS_MAX_SIZE;
    unsigned int			_READ_INDEX_MAX_SIZE;
    unsigned int			_READ_ID_MAX_SIZE;
    unsigned int 			_READ_DIST_MAX_SIZE;
    unsigned int 			_READ_LEN_MAX_SIZE;

    unsigned int 			_k;
    unsigned long int 		_q;
	unsigned long int 		_h;
	unsigned short	 		_BL;
	unsigned short			_overlap;
	unsigned short			_extensionLength;

	clock_t					_start;
	clock_t					_time1;
	clock_t					_time2;

    static Singleton*       _instance;

    Singleton();
    
public:

    static Singleton*			getInstance(void);

    void 						printInstance(void);

    inline unsigned int			CONTIG_ID_MAX_SIZE (void) { return _CONTIG_ID_MAX_SIZE; }
	inline unsigned int			NUM_SUBLISTS_MAX_SIZE(void) { return _NUM_SUBLISTS_MAX_SIZE; }
	inline unsigned int			READ_INDEX_MAX_SIZE (void) { return _READ_INDEX_MAX_SIZE; }
	inline unsigned int			READ_ID_MAX_SIZE(void) { return _READ_ID_MAX_SIZE; }
	inline unsigned int 		READ_DIST_MAX_SIZE(void) { return _READ_DIST_MAX_SIZE; }
	inline unsigned int 		READ_LEN_MAX_SIZE(void) { return _READ_LEN_MAX_SIZE; }

	inline unsigned int			k(void);
	inline unsigned long int 	q(void);
	inline unsigned long int 	h(void);
	inline unsigned short	 	BL(void);
	inline unsigned short		overlap(void);
	inline unsigned short		extensionLength(void);

	void						set_global_variables (unsigned int, unsigned int, unsigned int, unsigned int);
	void						set_layout_variables (unsigned int numReads, unsigned int);
	void						set_layout_variables_for_vector (unsigned int, unsigned int);
	bool 						feasibility_check ();

};
