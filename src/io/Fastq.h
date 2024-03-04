#ifndef FASTQ_H_
#define FASTQ_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
using namespace std;

#include "errors/Data_Exception.h"
#include "errors/Incorrect_Format.h"
using namespace errors;

#include "Auto_Unzip.h"
using namespace useful;

#include "common.h"
//#include "Masked_Sequence.h"

#define CCT_LETTERS 4

namespace fastq
{

// If no quality is specified, then we use value 't' (equivalent to PHREAD value 20
#define DEFAULT_QUALITY_CHAR 't'

#define DEFAULT_FASTA_COLUMNS 60

/**
 * \brief A simple class for handling fasta formatted files.
 *
 * This class can read a single sequence from a (multi)fasta file,
 * store information and using the sequence.
 * inline
 * At this point only nucleotide and protein sequences are allowed,
 * but quality file compatibility are planned.
 */
class Fastq
{
public:
	typedef vector<t_quality> t_quality_vector;
	enum FASTQ_type { standard, illumina };

	virtual ~Fastq();
	Fastq(bool cs = false);
	Fastq(const char *i, const char *s, bool cs = false);
	Fastq(const string &i, const string &s, bool cs = false);
	Fastq(const char *i, const char *s, const char *c, bool cs = false);
	Fastq(const string &i, const string &s, const string &c, bool cs = false);
	size_t length() const {	return sequence.length(); } //* Return the length of the sequence.
	void set_sequence(const char * str) { sequence = string(str); }
	void set_id(const char * str) { id = string(str); }
	void set_qual(bool v) {have_quality = v;}
	void set_quality(const char * str) { have_quality = true; quality = string(str); }
	const string & get_sequence() const { return sequence; } //* Return the sequence.
	const string & get_id() const { return id; } //* Return the id.
	const string & get_quality() const { return quality; } //* Return the quality.
	size_t size() const { return sequence.size(); }
	t_quality get_quality(size_t position) const throw (Data_Exception);
	t_quality_vector get_quality_vector() const;
	string description() const { return description(20); } //* Return a short rappresentation (20 characters) of the fasta.
	string description(size_t n) const;
	string reverse_complement() const;
	string reverse_complement(size_t start, size_t stop) const throw (Data_Exception);
	string reverse() const;
	string reverse(size_t start, size_t stop) const throw (Data_Exception);
//	void masking(Masked_Sequence& ms, char masking_character = masked_base);
	bool get_colorspace() const { return colorspace ; }
	void convertToColorspace();
	friend istream& operator>>(istream& buffer, Fastq& fasta) throw (Incorrect_Format);
	friend ostream& operator<<(ostream& buffer, const Fastq& fasta);

	bool get_mask_lowercase_flag() { return mask_lowercase; }
	void set_mask_lowercase_flag(bool flag = true) { mask_lowercase = flag; }


	size_t getColumns() const { return columns; }
	void setColumns(size_t columns) { this->columns = columns; }

	void loweringNs();

	static FASTQ_type check_FASTQ_type_file(const string &file);
	void set_input_type(FASTQ_type t) { input_type = t; }
	void set_input_type(const string &file) { input_type = check_FASTQ_type_file(file); }
	void set_output_type(FASTQ_type t) { output_type = t; }
	void set_output_type(const string &file) { output_type = check_FASTQ_type_file(file); }

	void setCoverage(unsigned int coverage) {this->cov = coverage;}
	unsigned int coverage() {return this->cov;}

protected:
	string id; ///< The id
	string sequence; ///< The sequence
	string quality; ///< Quality (compressed)
	string comment; ///< The comment
	bool   colorspace; ///< Data are in color space (SOLID technology)
	bool	have_quality;
	size_t	columns;
	bool	mask_lowercase;
	unsigned int cov;
	FASTQ_type input_type;
	FASTQ_type output_type;

	istream & read_from_fasta(istream & buffer) throw (Incorrect_Format);
	istream & read_from_fastq(istream & buffer) throw (Incorrect_Format);

};

}

#endif /*FASTQ_H_*/
