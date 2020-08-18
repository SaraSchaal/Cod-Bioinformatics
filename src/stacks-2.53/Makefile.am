SUBDIRS = htslib
HTSLIB = htslib/libhts.a

AM_CPPFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/htslib
AM_LDFLAGS = $(OPENMP_CFLAGS)
LDADD = $(HTSLIB) libcore.a

noinst_LIBRARIES = libcore.a libclean.a libpop.a
bin_PROGRAMS = ustacks cstacks sstacks process_radtags process_shortreads \
	kmer_filter clone_filter populations phasedstacks \
	tsv2bam gstacks

libcore_a_SOURCES = \
	src/aln_utils.h src/aln_utils.cc \
	src/BamI.h src/BamI.cc \
	src/BamUnalignedI.h \
	src/BowtieI.h \
	src/BustardI.h \
	src/catalog_utils.h src/catalog_utils.cc \
	src/cmb.h src/cmb.cc \
	src/constants.h src/constants.cc \
	src/DNANSeq.h src/DNANSeq.cc \
	src/DNASeq4.h src/DNASeq4.cc \
	src/DNASeq.h src/DNASeq.cc \
	src/FastaI.h \
	src/FastqI.h \
	src/GappedAln.h \
	src/gzFasta.h src/gzFasta.cc \
	src/gzFastq.h \
	src/input.h src/input.cc \
	src/kmers.h src/kmers.cc \
	src/locus.h src/locus.cc \
	src/locus_readers.h \
	src/log_utils.h src/log_utils.cc \
	src/MetaPopInfo.h src/MetaPopInfo.cc \
	src/models.h src/models.cc \
	src/mstack.h src/mstack.cc \
	src/mst.h src/mst.cc \
	src/nucleotides.h src/nucleotides.cc \
	src/renz.h src/renz.cc \
	src/SamI.h \
	src/Seq.h src/Seq.cc \
	src/sql_utilities.h src/sql_utilities.cc \
	src/stacks.h src/stacks.cc \
	src/Tsv.h \
	src/utils.h src/utils.cc \
	src/Vcf.h src/Vcf.cc \
	$(HTSLIB)

libclean_a_SOURCES = \
	src/clean.h src/clean.cc \
	src/file_io.h src/file_io.cc \
	src/write.h src/write.cc

libpop_a_SOURCES = \
	src/bootstrap.h \
	src/export_formats.h src/export_formats.cc \
	src/mapping_utils.h src/mapping_utils.cc \
	src/genotype_dictionaries.h src/genotype_dictionaries.cc \
	src/ordered.h \
	src/PopMap.h \
	src/PopSum.h src/PopSum.cc \
	src/smoothing.h \
	src/smoothing_utils.h \
	src/Hwp.h src/Hwp.cc

phasedstacks_SOURCES = libcore.a src/phasedstacks.h src/phasedstacks.cc

process_radtags_SOURCES    = libcore.a libclean.a src/process_radtags.h src/process_radtags.cc
process_shortreads_SOURCES = libcore.a libclean.a src/process_shortreads.h src/process_shortreads.cc
kmer_filter_SOURCES        = libcore.a libclean.a src/kmer_filter.h src/kmer_filter.cc
clone_filter_SOURCES       = libcore.a libclean.a src/clone_filter.h src/clone_filter.cc

ustacks_SOURCES = libcore.a src/ustacks.h src/ustacks.cc
cstacks_SOURCES = libcore.a src/cstacks.h src/cstacks.cc
sstacks_SOURCES = libcore.a src/sstacks.h src/sstacks.cc
tsv2bam_SOURCES = libcore.a src/tsv2bam.cc
gstacks_SOURCES = libcore.a src/gstacks.h src/gstacks.cc src/debruijn.h src/debruijn.cc src/Alignment.h src/SuffixTree.h src/SuffixTree.cc
populations_SOURCES = libcore.a libpop.a src/populations.h src/populations.cc

process_radtags_LDADD    = $(LDADD) libclean.a
process_shortreads_LDADD = $(LDADD) libclean.a
clone_filter_LDADD       = $(LDADD) libclean.a
kmer_filter_LDADD        = $(LDADD) libclean.a
populations_LDADD        = $(LDADD) libpop.a

dist_bin_SCRIPTS = scripts/denovo_map.pl scripts/ref_map.pl \
	scripts/integrate_alignments.py scripts/count_fixed_catalog_snps.py \
	scripts/stacks-integrate-alignments scripts/stacks-dist-extract scripts/stacks-gdb \
	scripts/stacks-samtools-tview scripts/stacks-count-reads-per-sample-per-locus \
	scripts/stacks-hist2d-loci-samples-coverage # scripts/denovo_map.py

dist_noinst_SCRIPTS = autogen.sh scripts/extract_interpop_chars.pl scripts/convert_stacks.pl

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) $(abs_top_srcdir)/tests/tap-driver.sh

TESTS = tests/process_radtags.t tests/kmer_filter.t tests/ustacks.t tests/pstacks.t

EXTRA_DIST = LICENSE INSTALL README ChangeLog $(TESTS)

pkglocalstatedir = $(localstatedir)/$(PACKAGE)

debug:
	$(MAKE) all "CXXFLAGS=-g -Wall -DDEBUG -O0"

install-data-hook:
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_BINDIR_,$(bindir)/,g' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/denovo_map.pl > $(DESTDIR)$(bindir)/denovo_map.pl.subst
	mv $(DESTDIR)$(bindir)/denovo_map.pl.subst $(DESTDIR)$(bindir)/denovo_map.pl
	chmod +x $(DESTDIR)$(bindir)/denovo_map.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_BINDIR_,$(bindir)/,g' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/ref_map.pl > $(DESTDIR)$(bindir)/ref_map.pl.subst
	mv $(DESTDIR)$(bindir)/ref_map.pl.subst $(DESTDIR)$(bindir)/ref_map.pl
	chmod +x $(DESTDIR)$(bindir)/ref_map.pl
	sed -i.bkp -e 's,_VERSION_,$(VERSION),' $(DESTDIR)$(bindir)/stacks-integrate-alignments
	sed -i.bkp -e 's,_VERSION_,$(VERSION),' $(DESTDIR)$(bindir)/stacks-samtools-tview
	rm -f $(DESTDIR)$(bindir)/stacks-integrate-alignments.bkp
	rm -f $(DESTDIR)$(bindir)/stacks-samtools-tview.bkp
	chmod +x $(DESTDIR)$(bindir)/stacks-integrate-alignments
	chmod +x $(DESTDIR)$(bindir)/stacks-dist-extract
	chmod +x $(DESTDIR)$(bindir)/stacks-gdb
	##sed -e 's,_VERSION_,$(VERSION),' -e 's,_BINDIR_,$(DESTDIR)$(bindir)/,g' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/denovo_map.py > $(DESTDIR)$(bindir)/denovo_map.py.subst
	##mv $(DESTDIR)$(bindir)/denovo_map.py.subst $(DESTDIR)$(bindir)/denovo_map.py
	##chmod +x $(DESTDIR)$(bindir)/denovo_map.py

uninstall-hook:
	rm -rf $(DESTDIR)$(pkglocalstatedir)
