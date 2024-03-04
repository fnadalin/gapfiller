
#include "data_structures/Reads.h"

Reads::Reads(const string & s, bool cs)
{
	//save read length
	information = s.length(); // 14 bits

	// set bit for first mate
	if(cs) {
		information |= (1 << 14);
	} else {
		information &= ~(1 << 14);
	}
	//set read to unused
	// information &= ~(1 << 17);
	// set bit to not Seed
	information &= ~(1 << 15);

	// unsigned short int bufferREAD_length = ceil((s.length())/(float)16);
	unsigned short int bufferREAD_length = (s.length() + 3) / 4;
	read = new uint8_t [bufferREAD_length];
	memset(read, 0, bufferREAD_length);
	int actualByte = 0;
	int actualBit = 0;

	/*
	for (unsigned int i = 0; i < s.length() ; i++) {
		switch (s.at(i)) {
		case 'A' : case 'a' : read[actualByte] &= ~(1 << (31-actualBit)); read[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		case 'C' : case 'c' : read[actualByte] &= ~(1 << (31-actualBit)); read[actualByte] |= (1 << (31-(actualBit+1)));  break;
		case 'G' : case 'g' : read[actualByte] |= (1 << (31-actualBit)); read[actualByte] &= ~(1 << (31-(actualBit+1)));   break;
		case 'T' : case 't' : read[actualByte] |= (1 << (31-actualBit)); read[actualByte] |= (1 << (31-(actualBit+1)));  break;
		default : read[actualByte] &= ~(1 << (31-actualBit)); read[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		}
		actualBit+=2;
		if(actualBit == 32) {
			actualByte++;
			actualBit=0;
		}
	}
	*/

	//cout << s << endl;

	for (unsigned int i = 0; i < s.length() ; i++) {
		switch (s.at(i)) {
		case 'A' : case 'a' : read[actualByte] &= ~(1 << (7-actualBit)); read[actualByte] &= ~(1 << (7-(actualBit+1))); break;
		case 'C' : case 'c' : read[actualByte] &= ~(1 << (7-actualBit)); read[actualByte] |= (1 << (7-(actualBit+1)));  break;
		case 'G' : case 'g' : read[actualByte] |= (1 << (7-actualBit)); read[actualByte] &= ~(1 << (7-(actualBit+1)));   break;
		case 'T' : case 't' : read[actualByte] |= (1 << (7-actualBit)); read[actualByte] |= (1 << (7-(actualBit+1)));  break;
		default : read[actualByte] &= ~(1 << (7-actualBit)); read[actualByte] &= ~(1 << (7-(actualBit+1))); break;
		}
		actualBit += 2;
		if(actualBit == 8) {
			actualByte++;
			actualBit=0;
		}
	}

	/*
	for (unsigned int i = 0; i < bufferREAD_length; i++) {
		cout << (uint32_t) read[i] << endl;
	}
	*/

	string sequence = reverse_complement_standalone_str(s);
	readRC = new uint8_t [bufferREAD_length];
	memset(readRC, 0, bufferREAD_length);
	actualByte = 0;
	actualBit = 0;

	/*
	for (unsigned int i = 0; i < sequence.length(); i++) {
		switch (sequence.at(i)) {
		case 'A' : case 'a' : readRC[actualByte] &= ~(1 << (31-actualBit)); readRC[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		case 'C' : case 'c' : readRC[actualByte] &= ~(1 << (31-actualBit)); readRC[actualByte] |= (1 << (31-(actualBit+1)));  break;
		case 'G' : case 'g' : readRC[actualByte] |= (1 << (31-actualBit)); readRC[actualByte] &= ~(1 << (31-(actualBit+1)));   break;
		case 'T' : case 't' : readRC[actualByte] |= (1 << (31-actualBit)); readRC[actualByte] |= (1 << (31-(actualBit+1)));  break;
		default : readRC[actualByte] &= ~(1 << (31-actualBit)); readRC[actualByte] &= ~(1 << (31-(actualBit+1))); break;
		}
		actualBit+=2;
		if(actualBit == 32) {
			actualByte++;
			actualBit=0;
		}
	}
	*/

	for (unsigned int i = 0; i < sequence.length(); i++) {
		switch (sequence.at(i)) {
		case 'A' : case 'a' : readRC[actualByte] &= ~(1 << (7-actualBit)); readRC[actualByte] &= ~(1 << (7-(actualBit+1))); break;
		case 'C' : case 'c' : readRC[actualByte] &= ~(1 << (7-actualBit)); readRC[actualByte] |= (1 << (7-(actualBit+1)));  break;
		case 'G' : case 'g' : readRC[actualByte] |= (1 << (7-actualBit)); readRC[actualByte] &= ~(1 << (7-(actualBit+1)));   break;
		case 'T' : case 't' : readRC[actualByte] |= (1 << (7-actualBit)); readRC[actualByte] |= (1 << (7-(actualBit+1)));  break;
		default : readRC[actualByte] &= ~(1 << (7-actualBit)); readRC[actualByte] &= ~(1 << (7-(actualBit+1))); break;
		}
		actualBit += 2;
		if(actualBit == 8) {
			actualByte++;
			actualBit=0;
		}
	}
}

Reads::Reads(const Reads &r)
{
	read = r.read;
	readRC = r.readRC;
	information = r.information;
}

Reads &
Reads::operator=(const Reads & r)
{
	if (this != &r) { // protect against invalid self-assignment
		this->read = r.read;
		this->readRC = r.readRC;
		this->information = r.information;
	}
	return *this;
}

unsigned short
Reads::length()
{
	// return (unsigned short int)(  uint16_t(information >> 16) << 16 | uint16_t(information) );
	return (uint16_t) (((1 << 14) - 1) & information);
}

string
Reads::toString (bool P)
{
	uint8_t Mask = 0;
	unsigned int block = 0;
	Mask |= 1 << 7;
	Mask |= 1 << 6;

	string out;
	out.reserve(length());

	uint8_t R;
	P ? R = read[0] : R = readRC[0];

	for(unsigned int i = 0 ; i < length() ; i++) {
		if(i%4 == 0 and i > 0) {
			block++;
			P ? R = read[block] : R = readRC[block];
		}
		uint8_t c = (R & Mask) >> 6;
		R = R << 2;
 		switch(c) {
			case 0: out.append("A"); break;
			case 1: out.append("C"); break;
			case 2: out.append("G"); break;
			case 3: out.append("T"); break;
			default: cout << "strange" << endl;
		}
	}

	return out;
}

/*
string
Reads::toString_deb (bool P)
{
	unsigned int Mask = 0;
	unsigned int block= 0;
	Mask |= 1 << 31;
	Mask |= 1 << 30;

	//string out = "";
	string out;
	out.resize(length());

	unsigned int R;
	P ? R = read[0] : R = readRC[0];

	for(unsigned int i = 0 ; i < length() ; i++) {
		if(i%16 == 0 and i > 0) {
			block++;
			P ? R = read[block] : R = readRC[block];
		}
		unsigned int c = (R & Mask) >> 30;
		cout << "c = " << c << endl;
		R = R << 2;
 		switch(c) {
			case 0: out.at(i) = 'A'; break;
			case 1: out.at(i) = 'C'; break;
			case 2: out.at(i) = 'G'; break;
			case 3: out.at(i) = 'T'; break;
			default: cout << "strange " << cout << "\n";
		}
 		cout << out << endl;
	}
	cout << out << endl;
	return out;
}
*/

