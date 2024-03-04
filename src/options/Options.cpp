#include "options/Options.h"

Options::Options ()
{
	help = HELP;
	k = K;
	blockLength = BLOCK_LENGTH;
	layout = LAYOUT;
	outputPrefix = OUTPUT_PREFIX;
	gzip = GZIP;
	bzip2 = BZIP2;
	seed_fasta = SEED_FASTA;
	seed_sam = SEED_SAM;
	query_fasta = QUERY_FASTA;
	query_sam = QUERY_SAM;
	seed_as_query = SEED_AS_QUERY;
	overlap = OVERLAP;
	globalMismatch = GLOBAL_MISMATCH;
	extThreshold = EXT_THRESHOLD;
	limitedExtension = LIMITED_EXTENSION;
	no_cycle = NO_CYCLE;
	reverse_seed = REVERSE_SEED;
	reverse_mate = REVERSE_MATE;
	verbose = VERBOSE;
}

void
Options::set_values (Menu & values)
{
	if (values._vm.count("help")) {
		cout << values._head << endl;
		cout << values._menu << endl;
		exit(0);
	}

	if (values._vm.count("k")) {
		k = values._vm["k"].as<int> ();
		if (k > 15 or k < 1) {
			cerr << "--k parameter must be a positive number less or equal to 15" << endl;
			cerr << "Try \"--help\" for help" << endl;
			exit(1);
		}
	}

	if (values._vm.count("block-length")) {
		blockLength = values._vm["block-length"].as<unsigned int> ();
	}

	if (values._vm.count("store-layout")) {
		if (values._vm.count("query")){
			cerr << "Cannot use both --query and --store-layout" << endl;
			cerr << "Try \"--help\" for help" << endl;
			exit(1);
		}
		layout = true;
	}

	if (values._vm.count("output-prefix")) {
		outputPrefix = values._vm["output-prefix"].as<string> ();
	}

	if (values._vm.count("gz") and values._vm.count("bz2")) {
		cerr << "Cannot use both --gz and --bz2" << endl;
		exit(1);
	}
	if (values._vm.count("gz")) {
		gzip = true;
	}
	if (values._vm.count("bz2")) {
		bzip2 = true;
	}

	if (not ((values._vm.count("seed1") and values._vm.count("seed2")) or values._vm.count("seed-sam"))) {
		cerr << "--seed1 and --seed2 options are required; otherwise, use --seed-sam option" << endl;
		cerr << "Try \"--help\" for help" << endl;
		exit(1);
	} else if (values._vm.count("seed1") and values._vm.count("seed2") and values._vm.count("seed-sam")) {
		cerr << "Cannot use both --seed1/--seed2 and --seed-sam" << endl;
		cerr << "Try \"--help\" for help" << endl;
		exit(1);
	} else {
		// set contig size
		if (not values._vm.count("seed-ins") or not values._vm.count("seed-var")) {
			cout << "--seed-ins and --seed-var options are required" << endl;
			cerr << "Try \"--help\" for help" << endl;
			exit(1);
		} else {
			seed_ins = values._vm["seed-ins"].as<unsigned int> ();
			seed_var = values._vm["seed-var"].as<unsigned int> ();
		}
		// store seed filenames
		if (values._vm.count("seed1") and values._vm.count("seed2")) {
			if (seed_fw_files.size() != seed_rv_files.size()) {
				cerr << "Different number of first and second seed reads files" << endl;
				exit(1);
			}
			seed_fasta = true;
			seed_fw_files = values._vm["seed1"].as<vector<string> > ();
			seed_rv_files = values._vm["seed2"].as<vector<string> > ();
		} else { // vm.count("seed-sam")
			seed_sam = true;
			seed_sam_files = values._vm["seed-sam"].as<vector<string> > ();
		}
	}

	if (not values._vm.count("query") and not values._vm.count("query-sam")) {
		cout << "Use seed reads for extension" << endl;
	} else {
		seed_as_query = false;
		if (values._vm.count("query")) {
			query_fasta = true;
			query_files = values._vm["query"].as<vector<string> > ();
		}
		if (values._vm.count("query-sam")) {
			query_sam = true;
			query_sam_files = values._vm["query-sam"].as<vector<string> > ();
		}
	}

	if (values._vm.count("overlap")) {
		overlap = values._vm["overlap"].as<unsigned int> ();
		if (overlap < blockLength) {
			overlap = blockLength;
			cerr << "--overlap set to " << blockLength << endl;
		}
	}

	if (values._vm.count("mismatch-rate")) {
		globalMismatch = values._vm["global-mismatch"].as<int> ();
	}

	if (values._vm.count("extThreshold")) {
		extThreshold = values._vm["extThreshold"].as<unsigned int> ();
	}

	if (values._vm.count("limit")) {
		limitedExtension = true;
		maxReadsToExtend = values._vm["limit"].as<unsigned int> ();
		cout << "Extend only " << maxReadsToExtend << " reads" << endl;
	}

	if (values._vm.count("no-read-cycle")) {
		no_cycle = true;
	}

	if (values._vm.count("mate-pairs")) {
		reverse_seed = true;
		reverse_mate = false;
	}

	if (values._vm.count("verbose")) {
		verbose = true;
	}
}
