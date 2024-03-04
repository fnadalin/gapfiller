#include "options/Menu.h"

Menu::Menu ()
{
	stringstream ss;
	ss << fast_description() << endl << endl << "Allowed options";
	_head = ss.str();
}

void
Menu::add_options ()
{
	_menu.add_options()
	("help", "produce help message")
	("k", po::value<int>(), "length of the word used to hash (default: 12)")
	("block-length", po::value<unsigned int>(), "length of perfect match (default: 15)")
	("output-prefix", po::value<string>(), "output files prefix (default: \"GapFiller_output\")")
	("gz", "compress output with gzip")
	("bz2", "compress output with bzip2")
	("seed1", po::value<vector<string> >(), "seed1 fasta file (can be compressed with gzip or bzip2, or a pipe)")
	("seed2", po::value<vector<string> >(), "seed2 fasta file (can be compressed with gzip or bzip2, or a pipe)")
	("seed-sam", po::value<vector<string> > (), "seed sam file sorted by ID, with header (sam or bam format; can be repeated multiple times)")
	("query", po::value<vector<string> >(), "query fasta file: use different reads for extension instead of seeds (can be compressed with gzip or bzip2, or a pipe)")
	("query-sam", po::value<vector<string> >(), "query sam file: use different reads for extension instead of sam seeds (sam or bam format; can be repeated multiple times)")
	("seed-ins", po::value<unsigned int>(), "seed reads insert size")
	("seed-var", po::value<unsigned int>(), "seed reads insert variation")
	("store-layout", "store contigs layout (default: false)\n")
	("overlap", po::value<unsigned int>(), "minimum suffix-prefix overlap (default: 30)")
	("mismatch-rate", po::value<int>(), "maximum number of mismatches every 100 bp (default: 5)")
	("extThreshold", po::value<unsigned int>(), "number of reads needed to extend a contig (default: 2)")
	("limit", po::value<unsigned int>(), "limits the number of extended reads (useful for tests)")
	("no-read-cycle", "allow reads to be used multiple times within the same contig (default: false)")
	("mate-pairs", "default: paired-ends")
	("verbose", "print a lot of information! Use with --limit option");
}

void
Menu::parse_options (int argc, char* argv [])
{
	try {
		po::store(po::parse_command_line(argc, argv, _menu), _vm);
		po::notify(_vm);
	} catch (boost::program_options::error & error) {
		cerr << error.what() << endl;
		cerr << "Try \"--help\" for help" << endl;
		exit(2);
	}
}
