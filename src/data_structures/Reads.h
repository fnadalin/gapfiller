#ifndef READS_H_
#define READS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <stdint.h>
using namespace std;

#include <common.h>

//#define returnBit(var,pos) (var & (1 << pos))
//#define setBit(var,pos)    (var |= (1 << pos))
//#define clearBit(var,pos)  (var &= ~(1 << pos))

class Reads {

public:

	virtual ~Reads() {}
	Reads() { };
	Reads(const string & s, bool cs = false);
	Reads(const Reads &r); // copy constructor

	Reads & 				operator=(const Reads & r);

	// inline uint8_t * 		get_read() { return read; } //* Return the read on +
	// inline uint8_t * 		get_readRV() { return readRC; } //* Return the read on -

	unsigned short 			length(); 		// Return the length of the sequence.
	inline bool				firstMate() { return (information & (1 << 14)); } 	// true is the read is the first mate
	// inline bool		 	used() { return (information & (1 << 17)); }		// has the read already been used
	inline bool 			seed() { return (information & (1 << 15)); }		// is the read a seed
	// inline void 			setUsed() { information |= (1 << 17); }				// set the read to used
	inline void 			setSeed() { information |= (1 << 14); }				// set the read to seed
	string  				toString(bool); // returns char representation of string
	// string  				toString_deb(bool); // returns char representation of string

	uint16_t 				information; // stores read length, paired position, used status, seed status

protected:

	uint8_t *				read; // read as given in input
	uint8_t *				readRC; // read as given in input

};

#endif /* READS_H_ */
