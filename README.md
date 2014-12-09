alternatives
============

Detection of alternatives by local assembly of short reads

Source
------

The alternatives source code is hosted on github:

	git@github.com:gt1/alternatives.git

Compilation of alternatives
---------------------------

alternatives needs libmaus [https://github.com/gt1/libmaus] . When libmaus
is installed in ${LIBMAUSPREFIX} then alternatives can be compiled and
installed in ${HOME}/alternatives using

	- autoreconf -i -f
	- ./configure --with-libmaus=${LIBMAUSPREFIX} \
		--prefix=${HOME}/alternatives
	- make install

Running the program
-------------------

The program expects a .hwt file produced by bwtb3m [https://github.com/gt1/bwtb3m] 
as input. This .hwt file can be produced using

	fagzToCompact verbose=0 outputfilename=reads.compact reads.fa.gz
	bwtb3m inputtype=compactstream sasamplingrate=32 isasamplingrate=32 \
		outputfilename=reads.bwt numthreads=1 reads.compact
	alternatives reads.hwt

The program will then construct an overlap graph and extract path alternatives
in files prefixed by reads_bubbles_multi_indel_. These files are stored in the FastA
format. A complete example can be found in the script example.sh in the source tree .


