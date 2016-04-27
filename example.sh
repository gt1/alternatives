#! /bin/bash
BUILDDIR=$PWD
LIBMAUS2VERSION=2.0.175-release-20160426230024
BIOBAMBAM2VERSION=2.0.40-release-20160425112514
BWTB3MVERSION=0.0.37-release-20160426231355
ALTERNATIVESVERSION=0.0.2-release-20160427104438
MAFFTURL=http://mafft.cbrc.jp/alignment/software/mafft-7.213-without-extensions-src.tgz
BWAVERSION=0.7.13
MEMGB=8
SAMTOOLSVERSION=1.3.1
export GENE=BRCA2

# number of processors on machine
PROC=`cat /proc/cpuinfo | egrep "^processor[[:space:]]*:" | wc -l`

# get annotation file
if [ ! -f refFlat.txt.gz ] ; then
	wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
fi

# download human reference and rename chromosomes
if [ ! -f GRCH38.fa ] ; then
	( for i in `seq 1 22` X Y ; do
		wget -O - ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/Primary_Assembly/assembled_chromosomes/FASTA/chr${i}.fa.gz
	done ) | gzip -d -c | awk '/^>/ {V=$0 ; sub(/.*chromosome /,">chr",V) ; sub(/,.*/,"",V) ; print V } !/^>/ {print}' > GRCH38.fa
fi

if [ ! -e tools/bin/pkg-config ] ; then
	curl --location https://pkg-config.freedesktop.org/releases/pkg-config-0.29.1.tar.gz | tar xzvf -
	cd pkg-config-0.29.1
	./configure --with-internal-glib --prefix=${BUILDDIR}/tools
	make -j${PROC} install
	cd ..
	rm -fR pkg-config-0.29.1
fi

export PATH=${BUILDDIR}/tools/bin:$PATH

if [ ! -e tools/lib/libz.a ] ; then
	curl --location http://downloads.sourceforge.net/project/libpng/zlib/1.2.8/zlib-1.2.8.tar.gz | tar xzvf -
	cd zlib-1.2.8
	./configure --prefix=${BUILDDIR}/tools
	make install
	cd ..
fi

# build libmaus2
if [ ! -f libmaus2/lib/libmaus2.a ] ; then
	wget -O - https://github.com/gt1/libmaus2/archive/${LIBMAUS2VERSION}.tar.gz | tar xzf -
	mv libmaus2-${LIBMAUS2VERSION} libmaus2-${LIBMAUS2VERSION}-src
	mkdir libmaus2-${LIBMAUS2VERSION}-build
	cd libmaus2-${LIBMAUS2VERSION}-build
	../libmaus2-${LIBMAUS2VERSION}-src/configure --prefix=${BUILDDIR}/libmaus2 --disable-compile-testprograms
	make -j${PROC} install
	cd ..
	rm -fR libmaus2-${LIBMAUS2VERSION}-src libmaus2-${LIBMAUS2VERSION}-build
fi

# build samtools
if [ ! -f samtools-${SAMTOOLSVERSION}/samtools ] ; then
	wget -O - http://downloads.sourceforge.net/project/samtools/samtools/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 | tar xjf -
	cd samtools-${SAMTOOLSVERSION}
	make -j${PROC}
	cd ..
fi

# compute fa index
if [ ! -f GRCH38.fa.fai ] ; then
	samtools-${SAMTOOLSVERSION}/samtools faidx GRCH38.fa
fi

# build biobambam2
if [ ! -f biobambam2/bin/bamtofastq ] ; then
	wget -O - https://github.com/gt1/biobambam2/archive/${BIOBAMBAM2VERSION}.tar.gz | tar xzf -
	mv biobambam2-${BIOBAMBAM2VERSION} biobambam2-${BIOBAMBAM2VERSION}-src
	mkdir -p biobambam2-${BIOBAMBAM2VERSION}-build
	cd biobambam2-${BIOBAMBAM2VERSION}-build
	../biobambam2-${BIOBAMBAM2VERSION}-src/configure --with-libmaus2=${BUILDDIR}/libmaus2 --prefix=${BUILDDIR}/biobambam2 --enable-install-experimental
	make -j${PROC} install
	cd ..
	rm -fR biobambam2-${BIOBAMBAM2VERSION}-build biobambam2-${BIOBAMBAM2VERSION}-src
fi

# build bwtb3m
if [ ! -f bwtb3m/bin/bwtb3m ] ; then
	wget -O - https://github.com/gt1/bwtb3m/archive/${BWTB3MVERSION}.tar.gz | tar xzf -
	mv bwtb3m-${BWTB3MVERSION} bwtb3m-${BWTB3MVERSION}-src
	mkdir -p bwtb3m-${BWTB3MVERSION}-build
	cd bwtb3m-${BWTB3MVERSION}-build
	../bwtb3m-${BWTB3MVERSION}-src/configure --with-libmaus2=${BUILDDIR}/libmaus2 --prefix=${BUILDDIR}/bwtb3m
	make -j${PROC} install
	cd ..
	rm -fR bwtb3m-${BWTB3MVERSION}-build bwtb3m-${BWTB3MVERSION}-src
fi
BWTB3M=${BUILDDIR}/bwtb3m/bin/bwtb3m
BWTB3MTOBWA=${BUILDDIR}/bwtb3m/bin/bwtb3mtobwa

# build bwa
if [ ! -f bwa-${BWAVERSION}/bwa ] ; then
	wget -O - https://github.com/lh3/bwa/releases/download/v${BWAVERSION}/bwa-${BWAVERSION}.tar.bz2 | tar xjf -
	cd bwa-${BWAVERSION}
	make -j${PROC}
	cd ..
fi
BWA=${BUILDDIR}/bwa-${BWAVERSION}/bwa

# compute index for bwa
if [ ! -f GRCH38.fa.bwt ] ; then
	"${BWA}" fa2pac "GRCH38.fa"
	"${BWTB3M}" inputtype=pacterm mem="${MEMGB}g" outputfilename=GRCH38.fa.pac.bwt GRCH38.fa.pac
	"${BWTB3MTOBWA}" GRCH38.fa.pac.bwt GRCH38.fa.bwt GRCH38.fa.sa
	"${BWA}" bwtupdate "GRCH38.fa.bwt"
	rm -f GRCH38.fa.pac.*
	"${BWA}" fa2pac -f "GRCH38.fa"
fi


# source of program for extracting reference regions from FastA
function extractProgram
{
cat <<EOF
/*
    faextract
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/aio/InputStreamInstance.hpp>
#include <libmaus2/bambam/BamAlignmentDecoderFactory.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/util/ArgInfo.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		std::string const fafn = arginfo.restargs.at(0);
		std::string const faidxfn = fafn + ".fai";
		std::string const srange = arginfo.restargs.at(1);
		libmaus2::bambam::CramRange range(srange);
		libmaus2::aio::InputStreamInstance idxCIS(faidxfn);
		libmaus2::fastx::FastAIndex faidx(idxCIS);
		// idxCIS.close();

		std::map<std::string,uint64_t> seqtoid;
		for ( uint64_t i = 0; i < faidx.sequences.size(); ++i )
		{
			libmaus2::fastx::FastAIndexEntry const & ie = faidx.sequences[i];
			std::string const & seqname = ie.name;
			seqtoid[seqname] = i;
		}
		
		std::string const & rangeseq = range.rangeref;
		
		if ( seqtoid.find(rangeseq) == seqtoid.end() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Cannot find sequence " << rangeseq << std::endl;
			lme.finish();
			throw lme;
		}

		libmaus2::aio::InputStreamInstance faistr(fafn);
		libmaus2::autoarray::AutoArray<char> const seq = faidx.readSequence(faistr, seqtoid.find(rangeseq)->second);
		
		int64_t const zstart = static_cast<int64_t>(range.rangestart)-1;
		int64_t const zend = static_cast<int64_t>(range.rangeend)-1;
		int64_t const zlen = zend-zstart+1;

		if ( zlen < 0 || zstart < 0 || zend >= seq.size() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Invalid input range" << std::endl;
			lme.finish();
			throw lme;
		}
		
		std::cout << ">" << srange << "\n";
		std::cout.write(seq.begin()+zstart,zlen);
		std::cout.put('\n');
		std::cout.flush();
		
		// std::cerr << "rangeseq=" << rangeseq << " rangestart=" << range.rangestart << " rangeend=" << range.rangeend << std::endl;

		// rangestart,rangestart
		// libmaus2::autoarray::AutoArray<char> readSequence(std::istream & in, int64_t const seqid)
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
EOF
}

# build program for extracting reference regions
if [ ! -f faextract ] ; then
	extractProgram > faextract.cpp
	FALIBS=`PKG_CONFIG_PATH=${BUILDDIR}/libmaus2/lib/pkgconfig ${BUILDDIR}/tools/bin/pkg-config --libs libmaus2`
	c++ -std=c++0x -Ilibmaus2/include -Llibmaus2/lib faextract.cpp -ofaextract -lmaus2 -Wl,-rpath=${BUILDDIR}/libmaus2/lib ${FALIBS}
	rm -f faextract.cpp
fi

# generate random read pairs from FastA input 
# (nucleotide sequence needs to be in a single line with no white spaces)
function randomreads
{
	grep -v ">" | tr '\n' ' ' | perl -p -e "s/\s*//g" | \
	awk -v depth=$1 -v insertsize=150 -v insertsizemul=25 -v readlen=150 -v erate=0.00025 \
	'
	function reversecomplement(s)
	{
		sl=length(s);
		t="";
		for ( j=1; j <= sl; ++j )
		{
			t = substr(s,j,1) "" t;
		}
		gsub(/A/,"X",t);
		gsub(/C/,"Y",t);
		gsub(/G/,"C",t);
		gsub(/T/,"A",t);
		gsub(/Y/,"G",t);
		gsub(/X/,"T",t);
		return t;
	}
	function hn(n)
	{
		t="";
		for ( j=1; j<=n; ++j )
		{
			t= t "H";
		}
		return t;
	}
	function replace(q)
	{
		tr = rand();
		if ( q == "A" ) { if ( tr < 0.333 ) return "C"; if ( tr < 0.666 ) return "G"; return "T"; }
		if ( q == "C" ) { if ( tr < 0.333 ) return "A"; if ( tr < 0.666 ) return "G"; return "T"; }
		if ( q == "G" ) { if ( tr < 0.333 ) return "A"; if ( tr < 0.666 ) return "C"; return "T"; }
		if ( q == "T" ) { if ( tr < 0.333 ) return "A"; if ( tr < 0.666 ) return "C"; return "G"; }
		return "N";
	}
	function err(s)
	{
		sl=length(s);
		t="";
		for ( j=1; j <= sl; ++j )
		{
			if ( rand() <= erate )	
				t = t replace(substr(s,j,1));
			else
				t = t "" substr(s,j,1);
		}
		
		return t;
	}
	!/^>/ {
	V=$1 ; 
	l=length(V) ; 
	for ( i=0; i < depth*(l/readlen); ++i ) 
	{
		insertdif=int((rand()-0.5)*insertsizemul);
		linsertsize=insertsize+insertdif;
		tlen=2*readlen+linsertsize;
		left=int((l-tlen+1)*rand());
		s1=substr(V,left,readlen);
		s2=substr(V,left+readlen+linsertsize,readlen);
		r1=reversecomplement(s1);
		r2=reversecomplement(s2);
		if ( rand() >= 0.5 )
		{
			printf "@read" i "_" left "s/1\n" err(s1) "\n+\n" hn(readlen) "\n";
			printf "@read" i "_" left "s/2\n" err(s2) "\n+\n" hn(readlen) "\n";
		}
		else
		{
			printf "@read" i "_" left "r/1\n" err(r2) "\n+\n" hn(readlen) "\n";
			printf "@read" i "_" left "r/2\n" err(r1) "\n+\n" hn(readlen) "\n";
		}
	}
}'
}

# generate read set
function genemod
{
# exon 10
#./faextract GRCH38.fa chr13:32332272-32333387

# modified BRCA2 gene with tandem repeat in exon 10
( ./faextract GRCH38.fa chr13:32315480-32332271 ; 
  ./faextract GRCH38.fa chr13:32332272-32332277 ; 
  ./faextract GRCH38.fa chr13:32332278-32332325 ; 
  ./faextract GRCH38.fa chr13:32332278-32332325 ; 
  ./faextract GRCH38.fa chr13:32332326-32399668 ) \
  | randomreads 12 # depth 12

# original BRCA2 gene +- 1000
CRANGE=`zcat refFlat.txt.gz | awk -v GENE=${GENE} '{ if ( match($0,GENE "\t") ) print $3 ":" $5-1000 "-" $6+1000 }'`
./faextract GRCH38.fa ${CRANGE} | randomreads 20

# original BRCA1 gene +- 1000 to simulate off target data
CRANGE=`zcat refFlat.txt.gz | awk -v GENE=BRCA1 '{ if ( match($0,GENE "\t") ) print $3 ":" $5-1000 "-" $6+1000 }'`
./faextract GRCH38.fa ${CRANGE} | randomreads 20
}

# we generate reads representing BRCA2 as in the reference and
# BRCA2 with a 48 base tandem repeat in exon10
# in addition we add some off target data to be filtered out
# before computing alternatives
# the reads are mapped to the human reference using BWA mem
# the resulting SAM files is sorted by coordinate and converted to BAM
if [ ! -f ${GENE}mod.bam ] ; then
	genemod | bwa-${BWAVERSION}/bwa mem -p GRCH38.fa - | biobambam2/bin/bamsort inputformat=sam > ${GENE}mod.bam
fi

INPUT=${GENE}mod.bam
REMAPPED="${INPUT%.bam}_remapped.bam"
ASSEMBLY="${INPUT%.bam}_assembly.bam"
NAMES="${INPUT%.bam}.names"
FILTERED="${INPUT%.bam}_filtered.bam"
SUBFILTERED="${INPUT%.bam}_subfiltered.bam"
FQ=${FILTERED%.bam}.fastq
FQ1=${FILTERED%.bam}_1.fastq
FQ2=${FILTERED%.bam}_2.fastq
THREADS=$2
RANGE=`zcat refFlat.txt.gz | awk -v GENE=${GENE} '{ if ( match($0,GENE "\t") ) print $3 ":" $5-1000 "-" $6+1000 }'`

if [ ! -f ${FILTERED} ] ; then
	if [ ! -f ${REMAPPED} ] ; then
		# convert BAM to FastQ, modify names to hide pair information fro BWA,
		# map reads as single end
		# sort by coordinate
		biobambam2/bin/bamtofastq O=/dev/null O2=/dev/null S=/dev/null < ${INPUT} | \
			awk 'NR%4==1 {sub(/\/1$/,"_one",$0) ; sub(/\/2$/,"_two",$0) ; print} NR%4!=1 {print}' |\
			${BWA} mem -t ${PROC} GRCH38.fa - | \
			biobambam2/bin/bamsort index=1 indexfilename=${REMAPPED}.bai inputformat=sam calmdnm=1 fixmates=1 calmdnmreference=GRCH37.fa sortthreads=${PROC} >${REMAPPED}
	fi

	# extract range around target and save names where at least one end map to the target region
	biobambam2/bin/bamcollate2 ranges="${RANGE}" I=${REMAPPED} | samtools-${SAMTOOLSVERSION}/samtools view -h - | awk -F '\t' '!/@/{sub(/_one|_two/,"",$1) ; print $1}' | sort -u > ${NAMES}
	# filter reads with matching names from BAM file
	biobambam2/bin/bamfilternames <${INPUT} names=${NAMES} outputthreads=${PROC} >${FILTERED}
	rm -f ${NAMES} ${REMAPPED} ${REMAPPED}.bai
fi

if [ ! -f ${FILTERED}.fa.gz ] ; then
	# extract reads as FastA
	biobambam2/bin/bamtofastq < "${FILTERED}" fasta=1 O=/dev/null O2=/dev/null S=/dev/null | gzip -9 > ${FILTERED%.bam}.fa.gz
fi

if [ ! -f ${FILTERED%.bam}.compact ] ; then
	# convert FastA to compact packed representation
	bwtb3m/bin/fagzToCompact verbose=0 outputfilename=${FILTERED%.bam}.compact ${FILTERED%.bam}.fa.gz
fi

if [ ! -f ${FILTERED%.bam}.bwt ] ; then
	# compute BWT of read set
	${BWTB3M} inputtype=compactstream sasamplingrate=32 isasamplingrate=32 outputfilename=${FILTERED%.bam}.bwt numthreads=1 ${FILTERED%.bam}.compact
fi

# compile alternatives program
if [ ! -f alternatives/bin/alternatives ] ; then
	wget -O - https://github.com/gt1/alternatives/archive/${ALTERNATIVESVERSION}.tar.gz | tar xzf -
	mv alternatives-${ALTERNATIVESVERSION} alternatives-${ALTERNATIVESVERSION}-src
	mkdir -p alternatives-${ALTERNATIVESVERSION}-build
	cd alternatives-${ALTERNATIVESVERSION}-build
	../alternatives-${ALTERNATIVESVERSION}-src/configure --with-libmaus2=${BUILDDIR}/libmaus2 --prefix=${BUILDDIR}/alternatives
	make -j${PROC} install
	cd ..
	rm -fR alternatives-${ALTERNATIVESVERSION}-src alternatives-${ALTERNATIVESVERSION}-build
fi

if [ ! -f ${FILTERED%.bam}.dot.gz ] ; then
	# compute assembly (overlap) graph and extract alternatives
	alternatives/bin/alternatives k=20 minreqscore=125 ${FILTERED%.bam}.hwt | gzip -9 > ${FILTERED%.bam}.dot.gz
fi

if [ -e `which neato` ] ; then
	if [ ! -f ${FILTERED%.bam}.svg.gz ] ; then
		# produce layout using graphviz
		zcat ${FILTERED%.bam}.dot.gz | sed 's|graph {|graph { layout="sfdp";|' | neato -Tsvg | gzip -9 > ${FILTERED%.bam}.svg.gz
	fi
else
	echo "Graphviz program neato not found, not generating visualisation"
fi

# compile mafft
if [ ! -f mafft/bin/mafft ] ; then
	if [ ! -f mafft-7.213-without-extensions-src.tgz ] ; then
		wget -c ${MAFFTURL}
	fi
	
	MAFFTINSTDIR=${PWD}/mafft
	tar xzf mafft-7.213-without-extensions-src.tgz
	cd mafft-7.213-without-extensions/core
	make PREFIX=${MAFFTINSTDIR} install
	cd ../..
	rm -fR mafft-7.213-without-extensions
fi

# check for alternatives with indels
for i in ${FILTERED%.bam}_bubbles_multi_indel*.fasta ; do
	# compute multiple alignment using mafft
	if [ ! -f ${i%.fasta}_clustal.fasta ] ; then
		mafft/bin/mafft --clustalout $i > ${i%.fasta}.clustal
	fi
	# align contigs to reference using bwa mem
	if [ ! -f ${i%.fasta}.bam ] ; then
		bwa-${BWAVERSION}/bwa mem GRCH38.fa $i | biobambam2/bin/bamsort inputformat=sam > ${i%.fasta}.bam
	fi
done
