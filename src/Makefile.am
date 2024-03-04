
ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=GapFiller 
GapFiller_SOURCES=GapFiller.cpp \
	data_structures/Hash.cpp data_structures/Hash.h \
	data_structures/Reads.cpp data_structures/Reads.h \
	data_structures/Contig.cpp data_structures/Contig.h \
	data_structures/Contig_layout.cpp data_structures/Contig_layout.h \
	data_structures/Layout_list.cpp data_structures/Layout_list.h \
	data_structures/Layout_vector.cpp data_structures/Layout_vector.h \
	data_structures/include.h \
	io/Auto_Unzip.cpp io/Auto_Unzip.h \
	io/Auto_Zip.cpp io/Auto_Zip.h \
	io/Fastq.cpp io/Fastq.h \
	io/Sam.cpp io/Sam.h \
	io/Binary.cpp io/Binary.h \
	errors/Data_Exception.cpp errors/Data_Exception.h \
	errors/Data_Not_Found.cpp errors/Data_Not_Found.h \
	errors/File_Not_Found.cpp errors/File_Not_Found.h \
	errors/Generic_Exception.cpp errors/Generic_Exception.h \
	errors/Incorrect_Format.cpp errors/Incorrect_Format.h \
	options/Menu.cpp options/Menu.h \
	options/Options.cpp options/Options.h \
	options/Singleton.cpp options/Singleton.h \
	modules/Assembly.cpp modules/Assembly.h \
	modules/Assembly_only_sequence.cpp modules/Assembly_only_sequence.h \
	modules/Assembly_store_layout.cpp modules/Assembly_store_layout.h \
	\
	samtools/bam_aux.c samtools/bam_import.c samtools/bam_pileup.c samtools/bam_index.c \
	samtools/bam.c samtools/bam.h \
	samtools/bgzf.c samtools/bgzf.h \
	samtools/faidx.c samtools/faidx.h \
	samtools/knetfile.c samtools/knetfile.h \
	samtools/kstring.c samtools/kstring.h \
	samtools/razf.c samtools/razf.h \
	samtools/sam.c samtools/sam.h \
	samtools/sam_header.c samtools/sam_header.h \
   	samtools/glf.h samtools/kaln.h \
   	samtools/khash.h samtools/klist.h \
   	samtools/ksort.h samtools/kstring.h \
   	samtools/kseq.h samtools/bam_endian.h \
	common.h

GapFiller_CXXFLAGS= -O3 -fpermissive $(BOOST_CPPFLAGS)
GapFiller_CFLAGS= -O3

GapFiller_LDFLAGS= $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_IOSTREAMS_LIBS) -lz
GapFiller_LIBS= $(GapFiller_LDFLAGS)

