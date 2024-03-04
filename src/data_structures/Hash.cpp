#include "data_structures/Hash.h"

void process_mem_usage(double& vm_usage, double& resident_set)
{
	using std::ios_base;
	using std::ifstream;
	using std::string;

	vm_usage = 0.0;
	resident_set = 0.0;

	ifstream stat_stream("/proc/self/stat", ios_base::in);

	// dummy vars for leading entries in stat that we don't care about
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;

	// the two fields we want
	unsigned long vsize;
	long rss;

	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
				>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
				>> utime >> stime >> cutime >> cstime >> priority >> nice >> O
				>> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}

Hash::Hash (unsigned int k_hash, unsigned int b, unsigned int o)
{
	k = k_hash;
	q = 0;
	for (unsigned int i = 0; i < 2*k; i++) {
		q |= 1 << i;
	}
	blockLength = b;
	h = 1;
	while (b > 1) {
		h = (4*h)%q;
		b--;
	}

	limit = 512*1024*1024;

	HASHcounter = new unsigned int [q];
	for(unsigned long int i = 0; i < q; i++) {
		HASHcounter[i] = 0;
	}

	overlap = o;
	numReads = 0;
	maxReadLength = 0;

	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "Hash created VM: " << vm_usage << "; RSS: " << rss_usage << endl;
}

void
Hash::delete_HASHcounter ()
{
	if (HASHcounter != NULL) {
		delete [] HASHcounter;
	}
}

Hash::~Hash()
{
	delete_HASHcounter ();
	cout << "Hash deleted" << endl;
}

void
Hash::store_reads (const Options & Opt)
{
	clock_t time1 = clock();
	if (Opt.seed_as_query) {
		if (Opt.seed_fasta) read_fastq (Opt.seed_fw_files, Opt.seed_rv_files);
		else if (Opt.seed_sam) read_sam (Opt.seed_sam_files);
	} else {
		if (Opt.query_fasta) read_fastq (Opt.query_files);
		if (Opt.query_sam) read_sam (Opt.query_sam_files);
	}
	clock_t time2 = clock();
	cout << endl;
	cout << "Time needed to read sequences: " << (double) (time2 - time1) / CLOCKS_PER_SEC << endl;

	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "Reads stored VM: " << vm_usage << "; RSS: " << rss_usage << endl;
	cout << "Total number of reads for extension: " << numReads << endl;
}

// fill HashValue computing fingerprints for each read as well as for its reverse-complement
// change block position within the read, depending on its strand
void
Hash::fill ()
{
	clock_t time1 = clock();
	HASHvalues.resize(numReads*2);
	unsigned long int fp;
	// compute # reads with fingerprint = fp
	for(unsigned int i = 0; i < numReads; i++) {
		fp = ComputeFingerprint(i, true, 0);  // compute fingerprint for original read
		HASHcounter[fp]++;
		fp = ComputeFingerprint(i, false, overlap-blockLength); // reverse complemented
		HASHcounter[fp]++;
	}
	// compute # reads with fingerprint < fp
	unsigned int t1 = HASHcounter[0];
	HASHcounter[0] = 0;
	for(unsigned int i = 1; i < q; i++) {
		int t2 = HASHcounter[i];
		HASHcounter[i] = HASHcounter[i-1] + t1;
		t1 = t2;
	}

	// assign read pointer and update # reads with fingerprint <= fp
	for(unsigned int i = 0; i < numReads; i++) {
		fp = ComputeFingerprint(i, true, 0);
		readOriented r1;
		r1.id = i; r1.strand = true;
		HASHvalues.at(HASHcounter[fp]) = r1;
		HASHcounter[fp]++;
		fp = ComputeFingerprint(i, false, overlap-blockLength);
		readOriented r2;
		r2.id = i; r2.strand = false;
		HASHvalues.at(HASHcounter[fp]) = r2;
		HASHcounter[fp]++;
	}

	// re-store # of reads with fingerprint < fp (so that it corresponds to the first position in HASHvalues)
	for(unsigned int i=q; i > 0 ; i--) {
		HASHcounter[i]= HASHcounter[i-1];
	}
	HASHcounter[0] = 0;

	double vm_usage, rss_usage;
	process_mem_usage(vm_usage, rss_usage);
	cout << "HASH->readsMulti populated VM: " << vm_usage << "; RSS: " << rss_usage << endl;
	clock_t time2 = clock();
	cout << "Time needed to populate HASH: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << endl;
}

unsigned long int
Hash::ComputeFingerprint (unsigned int pos, bool orientation, unsigned int start)
{
	string read = readsMulti.at(pos).toString(orientation);
	unsigned long int fingerprint = 0;
	for (unsigned int j = start; j < start + blockLength and j < read.length(); j++) {
		char c = 0;
		switch(read.at(j)) {
		case 'a' : case 'A' : c = 0; break;
		case 'c' : case 'C' : c = 1; break;
		case 'g' : case 'G' : c = 2; break;
		case 't' : case 'T' : c = 3; break;
		}
		fingerprint = ((fingerprint << 2) + c) % q;
	}
	return fingerprint;
}

void
Hash::read_fastq (const vector <string> & filename_for, const vector <string> & filename_rev)
{
	vector<string>::const_iterator iter_fw = filename_for.begin();
	vector<string>::const_iterator iter_rv = filename_rev.begin();
	while (iter_fw != filename_for.end() and iter_rv != filename_rev.end()) {

		Auto_Unzip input_file_fw(*iter_fw);
		Auto_Unzip input_file_rv(*iter_rv);
		istream & inputFileFW = input_file_fw.filtered();
		istream & inputFileRV = input_file_rv.filtered();

		Fastq read_fw;
		Fastq read_rv;
		Fastq::FASTQ_type type_fw = Fastq::check_FASTQ_type_file(iter_fw->c_str());
		Fastq::FASTQ_type type_rv = Fastq::check_FASTQ_type_file(iter_rv->c_str());
		read_fw.set_input_type(type_fw);
		read_rv.set_input_type(type_rv);

		while (not inputFileFW.eof() and not inputFileRV.eof()) {
			inputFileFW >> read_fw;
			inputFileRV >> read_rv;
			if (read_fw.length() > maxReadLength) maxReadLength = read_fw.length();
			if (read_rv.length() > maxReadLength) maxReadLength = read_rv.length();

			Reads r = Reads(read_fw.get_sequence(), true);
			readsMulti.push_back(r);
			r = Reads(read_rv.get_sequence(), false);
			readsMulti.push_back(r);
			if (readsMulti.size() % 100000 == 0) {
				cout << ".";
				cout.flush();
				if(readsMulti.size() % 10000000 == 0) {
					cout << endl;
				}
			}
		}
		iter_fw++;
		iter_rv++;
	}
	numReads += readsMulti.size();
}

void
Hash::read_fastq (const vector <string> & filename)
{
	vector<string>::const_iterator iter_query = filename.begin();
	while (iter_query != filename.end()) {
		Auto_Unzip input_file_query(*iter_query);
		istream & inputFilequery = input_file_query.filtered();

		Fastq read_query;
		Fastq::FASTQ_type type_query = Fastq::check_FASTQ_type_file(iter_query->c_str());
		read_query.set_input_type(type_query);

		while (not inputFilequery.eof()) {
			inputFilequery >> read_query;
			if (read_query.length() > maxReadLength) maxReadLength = read_query.length();

			Reads r = Reads(read_query.get_sequence(), true);
			readsMulti.push_back(r);
			if (readsMulti.size() % 100000 == 0) {
				cout << ".";
				cout.flush();
				if (readsMulti.size() % 10000000 == 0) {
					cout << endl;
				}
			}
		}
		iter_query++;
	}
	numReads += readsMulti.size();
}

void
Hash::read_sam (const vector <string> & filename) throw (Data_Exception)
{
	vector<string>::const_iterator iter_sam = filename.begin();
	while (iter_sam != filename.end()) {
		Sam samFile (*iter_sam);
		if (not samFile.is_open()) {
			throw Data_Exception (*iter_sam);
		}
		while (samFile.read ()) {
			if (samFile.sequence(true).length() > maxReadLength) {
				maxReadLength = samFile.sequence(true).length();
			}

			Reads r = Reads(samFile.sequence(not samFile.is_reverse()), samFile.is_first());
			readsMulti.push_back(r);
			if (readsMulti.size() % 100000 == 0) {
				cout << ".";
				cout.flush();
				if (readsMulti.size() % 10000000 == 0) {
					cout << endl;
				}
			}
		}
		iter_sam++;
	}
	numReads += readsMulti.size();
}
