/*
    alternatives
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
#include <config.h>

#include <libmaus2/alignment/SimpleLocalAligner.hpp>
#include <libmaus2/util/OutputFileNameTools.hpp>
#include <libmaus2/lcs/HammingOverlapDetection.hpp>
#include <libmaus2/parallel/PosixSpinLock.hpp>
#include <libmaus2/util/NumberSerialisation.hpp>
#include <libmaus2/graph/StronglyConnectedComponents.hpp>
#include <libmaus2/graph/TopologicalSorting.hpp>
#include <libmaus2/graph/IdentityTargetProjector.hpp>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <set>

// find sequence id for position in text
static uint64_t findSequence(std::vector<uint64_t> const & seqstarts, uint64_t const pos)
{
	// binary search
	std::vector<uint64_t>::const_iterator ita = std::lower_bound(seqstarts.begin(),seqstarts.end(),pos);

	// behind last start
	if ( ita == seqstarts.end() )
		return (ita-seqstarts.begin())-1;
	else
	{
		// on a start?
		if ( pos == *ita )
			return ita-seqstarts.begin();
		// not on a sequence start
		else
			return (ita-seqstarts.begin())-1;
	}
}

// get length of sequence for sequence id
static uint64_t getSeqLen(std::vector<uint64_t> const & seqstarts, uint64_t const seqid)
{
	uint64_t const seqid2 = seqid>>1;
	return seqstarts.at(2*seqid2+1)	- seqstarts.at(2*seqid2) - 1;
}

// read container
struct ReadContainer
{
	typedef libmaus2::fm::BidirectionalDnaIndexImpCompactHuffmanWaveletTree index_type;
	typedef index_type::lf_type lf_type;
	index_type const & index;
	libmaus2::fm::SampledISA<lf_type> const & ISA;
	std::vector<uint64_t> const seqstart;

	ReadContainer(
		index_type const & rindex,
		libmaus2::fm::SampledISA<lf_type> const & rISA
	) : index(rindex), ISA(rISA), seqstart(index.getSeqStartPositions())
	{

	}

	// get read for id
	std::string getRead(uint64_t const readid) const
	{
		uint64_t const start  = seqstart[2*readid];
		uint64_t const len    = getSeqLen(seqstart,2*readid);
		uint64_t const r0     = ISA[start];
		std::string const uref = index.getTextUnmapped(r0,len);
		return uref;
	}
};

// return true if edge orientations represent an a b c type overlap
static bool isABCOverlap(
        ::libmaus2::lcs::OverlapOrientation::overlap_orientation const abo,
        ::libmaus2::lcs::OverlapOrientation::overlap_orientation const aco,
        ::libmaus2::lcs::OverlapOrientation::overlap_orientation const bco
        )
{
	return (
                // abc
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back  && // ab
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back  && // ac
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back     // bc
                )
                // Abc
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back && // Ab
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back && // Ac
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back     // bc
                )
                // aBc
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front && // aB
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back  && // ac
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back    // Bc
                )
                // ABc
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front  && // AB
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back && // Ac
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back    // Bc
                )
                // abC
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back  && // ab
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front && // aC
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front    // bC
                )
                // AbC
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back && // Ab
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front  && // AC
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front    // bC
                )
                // aBC
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front && // aB
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front && // aC
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front     // BC
                )
                // ABC
                ||
                (
                        abo==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front  && // AB
                        aco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front  && // AC
                        bco==libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front     // BC
                )
        );
}

// graph edge without edge source id
struct OverlapEntry
{
	// overlap type
	libmaus2::lcs::OverlapOrientation::overlap_orientation orientation;
	// b overhang behind overlap of a and b
	uint64_t overhang;
	// score of match between a and b
	int64_t score;
	// overlap between a and b
	uint64_t overlap;
	// edge target
	uint64_t target;

	OverlapEntry()
	{
	}

	OverlapEntry(
		libmaus2::lcs::OverlapOrientation::overlap_orientation rorientation,
		uint64_t roverhang,
		int64_t rscore,
		uint64_t roverlap,
		uint64_t rtarget
	) : orientation(rorientation), overhang(roverhang), score(rscore), overlap(roverlap), target(rtarget) {}

	libmaus2::lcs::OverlapOrientation::overlap_orientation decodeOverlapOrientation(uint64_t const v)
	{
		switch ( v )
		{
			case libmaus2::lcs::OverlapOrientation::overlap_cover_complete: return libmaus2::lcs::OverlapOrientation::overlap_cover_complete;
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front: return libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front;
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back: return libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back;
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front: return libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front;
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back: return libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back;
			case libmaus2::lcs::OverlapOrientation::overlap_a_covers_b:    return libmaus2::lcs::OverlapOrientation::overlap_a_covers_b;
			case libmaus2::lcs::OverlapOrientation::overlap_b_covers_a:    return libmaus2::lcs::OverlapOrientation::overlap_b_covers_a;
			case libmaus2::lcs::OverlapOrientation::overlap_ar_covers_b:   return libmaus2::lcs::OverlapOrientation::overlap_ar_covers_b;
			case libmaus2::lcs::OverlapOrientation::overlap_b_covers_ar:   return libmaus2::lcs::OverlapOrientation::overlap_b_covers_ar;
			case libmaus2::lcs::OverlapOrientation::overlap_a_complete_b:  return libmaus2::lcs::OverlapOrientation::overlap_a_complete_b;
			case libmaus2::lcs::OverlapOrientation::overlap_ar_complete_b: return libmaus2::lcs::OverlapOrientation::overlap_ar_complete_b;
			default:
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "cannot decode " << v << " as libmaus2::lcs::OverlapOrientation::overlap_orientation" << std::endl;
				lme.finish();
				throw lme;
			}
		}
	}

	void serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,static_cast<uint64_t>(orientation));
		libmaus2::util::NumberSerialisation::serialiseNumber(out,static_cast<uint64_t>(overhang));
		libmaus2::util::NumberSerialisation::serialiseSignedNumber(out,static_cast<uint64_t>(score));
		libmaus2::util::NumberSerialisation::serialiseNumber(out,static_cast<uint64_t>(overlap));
		libmaus2::util::NumberSerialisation::serialiseNumber(out,static_cast<uint64_t>(target));
	}
	OverlapEntry(std::istream & in)
	:
		orientation(decodeOverlapOrientation(libmaus2::util::NumberSerialisation::deserialiseNumber(in))),
		overhang(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		score(libmaus2::util::NumberSerialisation::deserialiseSignedNumber(in)),
		overlap(libmaus2::util::NumberSerialisation::deserialiseNumber(in)),
		target(libmaus2::util::NumberSerialisation::deserialiseNumber(in))
	{}
};

// projector to edge target vertex id
struct OverlapEntryTargetProjector
{
	uint64_t operator()(OverlapEntry const & OE) const
	{
		return OE.target;
	}
};

// comparator for sorting in non ascending overlap order
struct OverlapEntryOverlapComparator
{
	bool operator()(OverlapEntry const & A, OverlapEntry const & B) const
	{
		return A.overlap > B.overlap;
	}
};

// output operator for OverlapEntry
std::ostream & operator<<(std::ostream & out, OverlapEntry const & oe)
{
	return out << "OverlapEntry("
		<< oe.orientation << ","
		<< "overhang=" << oe.overhang << ","
		<< "score=" << oe.score << ","
		<< "overlap=" << oe.overlap << ","
		<< "target=" << oe.target << ")";
}

// load set of edges and duplicate information from file primedgesname
void loadEdges(std::map< uint64_t,std::vector<OverlapEntry> > & preedges, std::set<uint64_t> & dup, std::string const & primedgesname)
{
	std::cerr << "loading serialised edges...";
	libmaus2::aio::InputStreamInstance CIS(primedgesname);

	// number of source vertices
	uint64_t const numentries = libmaus2::util::NumberSerialisation::deserialiseNumber(CIS);

	for ( uint64_t i = 0; i < numentries; ++i )
	{
		// source vertex index
		uint64_t const s = libmaus2::util::NumberSerialisation::deserialiseNumber(CIS);
		// number of edges for s
		uint64_t const n = libmaus2::util::NumberSerialisation::deserialiseNumber(CIS);

		for ( uint64_t j = 0; j < n; ++j )
		{
			OverlapEntry OE(CIS);
			preedges[s].push_back(OE);
		}
	}

	// number of dup vertices
	uint64_t const numdup = libmaus2::util::NumberSerialisation::deserialiseNumber(CIS);
	for ( uint64_t i = 0; i < numdup; ++i )
		dup.insert(libmaus2::util::NumberSerialisation::deserialiseNumber(CIS));

	std::cerr << "done." << std::endl;
}

// save set of edges and duplicate information from primedgesname
void saveEdges(std::map< uint64_t,std::vector<OverlapEntry> > const & preedges, std::set<uint64_t> const & dup, std::string const & primedgesname)
{
	std::cerr << "serialising edges...";
	libmaus2::aio::OutputStreamInstance COS(primedgesname);
	libmaus2::util::NumberSerialisation::serialiseNumber(COS,preedges.size());
	for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(COS,ita->first);
		libmaus2::util::NumberSerialisation::serialiseNumber(COS,ita->second.size());
		for ( uint64_t i = 0; i < ita->second.size(); ++i )
			ita->second[i].serialise(COS);
	}
	libmaus2::util::NumberSerialisation::serialiseNumber(COS,dup.size());
	for ( std::set<uint64_t>::const_iterator ita = dup.begin(); ita != dup.end(); ++ita )
		libmaus2::util::NumberSerialisation::serialiseNumber(COS,*ita);
	std::cerr << "done." << std::endl;
}

// simple pile up class
struct Pileup
{
	// piles
	std::map< uint64_t,std::vector<char> > M;

	// add symbol c in column i
	void add(uint64_t const i, char const c)
	{
		M[i].push_back(c);
	}

	// add string s starting from column i
	void add(uint64_t const i, std::string const & s)
	{
		for ( uint64_t j = 0; j < s.size(); ++j )
			add(i+j,s[j]);
	}

	// get consensus for column i and erase column i
	std::pair<char,uint64_t> consensus(uint64_t const i)
	{
		char cons = 'N';
		uint64_t support = 0;

		if ( M.find(i) != M.end() )
		{
			std::map<char,uint64_t> C;
			std::vector<char> const & V = M.find(i)->second;
			for ( uint64_t j = 0; j < V.size(); ++j )
				C[V[j]]++;

			std::vector < std::pair<uint64_t,char> > VV;
			for ( std::map<char,uint64_t>::const_iterator ita = C.begin(); ita != C.end(); ++ita )
			{
				VV.push_back(std::pair<uint64_t,char>(ita->second,ita->first));

				if (
					ita->first == 'a' ||
					ita->first == 'A' ||
					ita->first == 'c' ||
					ita->first == 'C' ||
					ita->first == 'g' ||
					ita->first == 'G' ||
					ita->first == 't' ||
					ita->first == 'T'
				)
					support += ita->second;
			}
			std::sort(VV.begin(),VV.end());

			if ( VV.size() == 1 )
				cons = VV[0].second;
			else
			{
				if ( VV[VV.size()-1].first == VV[VV.size()-2].first )
					cons = 'N';
				else
					cons = VV[VV.size()-1].second;
			}

			M.erase(M.find(i));
		}

		return std::pair<char,uint64_t>(cons,support);
	}
};

struct CurEndBase
{
	// current end on vertex
	enum cur_end_type { cur_end_left, cur_end_right };
};

std::ostream & operator<<(std::ostream & out, CurEndBase::cur_end_type const cur_end)
{
	switch ( cur_end )
	{
		case CurEndBase::cur_end_left:
			return out << "cur_end_left";
		case CurEndBase::cur_end_right:
		default:
			return out << "cur_end_right";
	}
}

// set of edges
struct EdgeSet : public CurEndBase
{
	// edge store
	std::map< uint64_t,std::vector<OverlapEntry> > & preedges;

	EdgeSet(std::map< uint64_t,std::vector<OverlapEntry> > & rpreedges) : preedges(rpreedges) {}

	// print set of edges
	std::ostream & printEdges(std::ostream & out) const
	{
		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		{
			std::vector<OverlapEntry> const & V = ita->second;

			for ( uint64_t i = 0; i < V.size(); ++i )
				out << std::setw(6) << ita->first << std::setw(0) << " -> " << V[i] << "\n";
		}

		return out;
	}


	static cur_end_type decodeCurEndType(uint64_t const i)
	{
		switch ( i )
		{
			case cur_end_left:
				return cur_end_left;
			case cur_end_right:
				return cur_end_right;
			default:
			{
				libmaus2::exception::LibMausException lme;
				lme.getStream() << "Unknown value in decodeCurEndType\n";
				lme.finish();
				throw lme;
			}
		}
	}

	// null vertex container
	struct NullContainer
	{
		void operator()(uint64_t const) {}
	};

	// front adding vertex container
	struct DequeFrontContainer
	{
		std::deque<uint64_t> & D;

		DequeFrontContainer(std::deque<uint64_t> & rD) : D(rD) {}

		void operator()(uint64_t const v)
		{
			D.push_front(v);
		}
	};

	// back adding vertex container
	struct DequeBackContainer
	{
		std::deque<uint64_t> & D;

		DequeBackContainer(std::deque<uint64_t> & rD) : D(rD) {}

		void operator()(uint64_t const v)
		{
			D.push_back(v);
		}
	};

	// result of follow unique path operation
	struct FollowUniqueResult
	{
		// end vertex
		uint64_t v;
		// number of edges used
		uint64_t c;
		// current end were path could not be continued
		cur_end_type cur_end;

		FollowUniqueResult() : v(0), c(0), cur_end(cur_end_left) {}
		FollowUniqueResult(
			uint64_t const rv,
			uint64_t const rc,
			cur_end_type const rcur_end
		)
		: v(rv), c(rc), cur_end(rcur_end)
		{

		}
	};

	// returns true if OE is a left edge (outgoing edge is from front of read)
	static bool isLeftEdge(OverlapEntry const & OE)
	{
		return
			OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front
			||
			OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back;
	}

	// get set of edges pointing into v by using edges leaving v and looking back at v from the target vertices
	std::vector< std::pair<uint64_t,OverlapEntry> > getReverseIncoming(uint64_t const v)
	{
		std::vector< std::pair<uint64_t,OverlapEntry> > V;
		std::set<uint64_t> S;

		if ( preedges.find(v) != preedges.end() )
		{
			std::vector<OverlapEntry> const & V = preedges.find(v)->second;
			for ( uint64_t i = 0; i < V.size(); ++i )
				S.insert(V[i].target);
		}

		for ( std::set<uint64_t>::const_iterator ita = S.begin(); ita != S.end(); ++ita )
		{
			uint64_t const rv = *ita;

			if ( preedges.find(rv) != preedges.end() )
			{
				std::vector<OverlapEntry> const & RV = preedges.find(rv)->second;

				for ( uint64_t i = 0; i < RV.size(); ++i )
					if ( RV[i].target == v )
						V.push_back(std::pair<uint64_t,OverlapEntry>(rv, RV[i]));
			}
		}

		return V;
	}

	// returns true if OE is a right edge (outgoing edge is from the back of the read)
	static bool isRightEdge(OverlapEntry const & OE)
	{
		return
			OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front
			||
			OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back;
	}

	// returns true if i is a border vertex (has edges either only and the left or only on the right)
	bool isBorderVertex(uint64_t const i) const
	{
		if  ( preedges.find(i) == preedges.end() )
			return false;

		return
			(countLeftEdges(i) == 0 && countRightEdges(i) != 0)
			||
			(countLeftEdges(i) != 0 && countRightEdges(i) == 0);
	}

	// count number of left edges in vector V
	static uint64_t countLeftEdges(std::vector<OverlapEntry> const & V)
	{
		uint64_t c = 0;
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( isLeftEdge(V[i]) )
				++c;

		return c;
	}

	// count number of right edges in vector V
	static uint64_t countRightEdges(std::vector<OverlapEntry> const & V)
	{
		uint64_t c = 0;
		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( isRightEdge(V[i]) )
				++c;

		return c;
	}

	// check whether graph is symmetric, returns true iff it is
	bool checkSymmetry() const
	{
		bool ok = true;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		{
			uint64_t const s = ita->first;
			std::vector<OverlapEntry> const & T = ita->second;

			for ( uint64_t i = 0; i < T.size(); ++i )
			{
				OverlapEntry const t = T[i];
				bool const okf = hasEdge(s,t.target);
				bool const okr = hasReverseEdge(s,t.target);

				assert ( okf );

				if ( ! okr )
				{
					std::cerr << "reverse edge for " << s << " -> " << getEdge(s,t.target) << " is missing" << std::endl;
				}

				ok = ok && okf && okr;
			}
		}

		return ok;
	}

	/**
	 * return true if graph contains the reverse edge for (s,t) or if (s,t) does not exist
	 **/
	bool hasReverseEdge(uint64_t const s, uint64_t const t) const
	{
		bool const have_s_t = hasEdge(s,t);
		bool const have_t_s = hasEdge(t,s);

		if ( have_s_t != have_t_s )
			return false;
		else if ( !have_s_t )
			return true;
		else
		{
			OverlapEntry const forward = getEdge(s,t);
			OverlapEntry const reverse = getEdge(t,s);

			switch ( forward.orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_cover_complete:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_cover_complete;
				case libmaus2::lcs::OverlapOrientation::overlap_a_covers_b:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_b_covers_a;
				case libmaus2::lcs::OverlapOrientation::overlap_b_covers_a:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_covers_b;
				case libmaus2::lcs::OverlapOrientation::overlap_ar_covers_b:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_b_covers_ar;
				case libmaus2::lcs::OverlapOrientation::overlap_b_covers_ar:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_ar_covers_b;
				case libmaus2::lcs::OverlapOrientation::overlap_a_complete_b:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_complete_b;
				case libmaus2::lcs::OverlapOrientation::overlap_ar_complete_b:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_ar_complete_b;
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back;
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back;
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front;
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
					return reverse.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front;
				default:
					return false;
			}
		}
	}

	// extract tip starting at border edge preseed
	// stores encountered vertices in C
	// number of edges used is stored in FUR
	// last vertex visited is also stored in FUR
	// returns true if preseed is at the end of a tip structure
	bool extractTip(uint64_t const preseed, std::deque<uint64_t> & C, FollowUniqueResult & FUR, uint64_t const maxfollow = std::numeric_limits<uint64_t>::max()) const
	{
		if ( isLeftTipSeed(preseed) )
		{
			C.push_back(preseed);
			FUR = followUniqueLeftRecord(preseed,C,maxfollow);
		}
		else if ( isRightTipSeed(preseed) )
		{
			C.push_back(preseed);
			FUR = followUniqueRightRecord(preseed,C,maxfollow);
		}

		// if tip seed and possible tip is a most maxfollow vertices long and next
		// vertex on chain has a link back
		if (
			C.size()
			&&
			C.size() <= maxfollow
			&&
			(
				(
					FUR.cur_end == cur_end_left
					&&
					countLeftEdges(FUR.v) == 1
					&&
					hasReverseEdge(FUR.v,getFirstLeftEdge(FUR.v).target)

					&&

					(
						(
						isLeftEdge(getEdge(getFirstLeftEdge(FUR.v).target,FUR.v))
						&&
						countLeftEdges(getFirstLeftEdge(FUR.v).target) > 1
						)
						||
						(
						isRightEdge(getEdge(getFirstLeftEdge(FUR.v).target,FUR.v))
						&&
						countRightEdges(getFirstLeftEdge(FUR.v).target) > 1
						)
					)
				)
				||
				(
					FUR.cur_end == cur_end_right
					&&
					countRightEdges(FUR.v) == 1
					&&
					hasReverseEdge(FUR.v,getFirstRightEdge(FUR.v).target)

					&&

					(
						(
						isLeftEdge(getEdge(getFirstRightEdge(FUR.v).target,FUR.v))
						&&
						countLeftEdges(getFirstRightEdge(FUR.v).target) > 1
						)
						||
						(
						isRightEdge(getEdge(getFirstRightEdge(FUR.v).target,FUR.v))
						&&
						countRightEdges(getFirstRightEdge(FUR.v).target) > 1
						)
					)

				)
			)
		)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	// comparator projecting onto second component of a pair
	template<typename type_a, typename type_b, template<typename> class base_comparator>
	struct PairSecondComparator
	{
		base_comparator<type_b> comparator;

		PairSecondComparator() : comparator()
		{

		}

		bool operator()(
			std::pair<type_a,type_b> const & A,
			std::pair<type_a,type_b> const & B
		) const
		{
			return comparator(A.second,B.second);
		}
	};

	// remove (single) edge
	void removeEdge(uint64_t const from, uint64_t const to)
	{
		if ( preedges.find(from) != preedges.end() )
		{
			std::vector<OverlapEntry> & edges = preedges.find(from)->second;

			uint64_t o = 0;
			for ( uint64_t i = 0; i < edges.size(); ++i )
				if ( edges[i].target != to )
					edges[o++] = edges[i];

			edges.resize(o);

			if ( ! o )
				preedges.erase(preedges.find(from));
		}
	}

	static uint64_t computeId(uint64_t const vertex, cur_end_type cur_end)
	{
		return (vertex << 1) | ((cur_end==cur_end_left) ? 0 : 1);
	}

	static uint64_t computeId(std::pair<uint64_t,cur_end_type> const & P)
	{
		return computeId(P.first,P.second);
	}

	static std::pair<uint64_t,cur_end_type> decodeId(uint64_t const vertex)
	{
		return std::pair<uint64_t,cur_end_type>(vertex>>1,decodeCurEndType(vertex&1));
	}

	std::pair< std::vector< uint64_t >, std::vector< uint64_t > > getStronglyConnectedComponents(
		uint64_t const root, cur_end_type rcur_end
	)
	{
		// extract edges
		std::map< uint64_t,std::vector<uint64_t> > const subedgesshift = extractConnectedSubEdgesShift(root,rcur_end);

		std::pair< std::vector< uint64_t >, std::vector< uint64_t > > SCCshift =
			libmaus2::graph::StronglyConnectedComponents::strongConnectContract<uint64_t,libmaus2::graph::IdentityTargetProjector>(
			subedgesshift,computeId(root,rcur_end)
		);

		std::set < std::set<uint64_t> > S;
		for ( uint64_t i = 1; i < SCCshift.second.size(); ++i )
		{
			std::set<uint64_t> R;

			for ( uint64_t j = SCCshift.second[i-1]; j < SCCshift.second[i]; ++j )
				R.insert(SCCshift.first[j] >> 1);

			S.insert(R);
		}

		std::vector< uint64_t > V0, V1;
		for ( std::set < std::set<uint64_t> >::const_iterator sita = S.begin(); sita != S.end(); ++sita )
		{
			std::set<uint64_t> const & subS = *sita;
			for ( std::set<uint64_t>::const_iterator ita = subS.begin(); ita != subS.end(); ++ita )
				V0.push_back(*ita);
			V1.push_back(subS.size());
		}

		uint64_t acc = 0;
		for ( uint64_t i = 0; i < V1.size(); ++i )
		{
			uint64_t const t = V1[i];
			V1[i] = acc;
			acc += t;
		}
		V1.push_back(acc);

		return std::pair< std::vector< uint64_t >, std::vector< uint64_t > >(V0,V1);
	}

	// extract component of end annotated vertices
	std::map< uint64_t,std::vector<uint64_t> > extractConnectedSubEdgesShift(uint64_t const root, cur_end_type const rcur_end)
	{
		// extracted edges map
		std::map< uint64_t,std::vector<uint64_t> > subedges;

		// set of annotated vertices seen so far
		std::set< std::pair<uint64_t,cur_end_type> > seen;
		// todo stack
		std::stack< std::pair<uint64_t,cur_end_type> > S;
		// push root
		S.push(std::pair<uint64_t,cur_end_type>(root,rcur_end));

		// extract set of 1d directed edges
		while ( !S.empty() )
		{
			// get node
			std::pair<uint64_t,cur_end_type> const P = S.top(); S.pop();
			// makr it as seen
			seen.insert(P);

			// compute id
			uint64_t const sid = computeId(P);

			// see if there are any outgoing edges
			if ( preedges.find(P.first) != preedges.end() )
			{
				// get edges
				std::vector<OverlapEntry> const & ledges = preedges.find(P.first)->second;

				// iterate over edges
				for ( uint64_t i = 0; i < ledges.size(); ++i )
				{
					OverlapEntry const & OE = ledges[i];

					// consider only edges on the correct side
					if (
						( (P.second == cur_end_left) && isLeftEdge(OE) )
						||
						( (P.second == cur_end_right) && isRightEdge(OE) )
					)
					{
						// compute new end of path continuation
						cur_end_type target_end = cur_end_left;

						switch ( OE.orientation )
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
								target_end = cur_end_right;
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								target_end = cur_end_left;
								break;
							default:
								break;
						}

						// target pair
						std::pair<uint64_t,cur_end_type> const T(OE.target,target_end);
						// compute target pair id
						uint64_t const tid = computeId(T);

						// insert edge
						subedges[sid].push_back(tid);

						// mark as todo if not already finished
						if ( seen.find(T) == seen.end() )
							S.push(T);
					}
				}
			}
		}

		return subedges;
	}

	// extract edges starting from root and store the employed edges (without their reverse edges if not used)
	std::map< uint64_t,std::vector<OverlapEntry> > extractConnectedSubEdges(uint64_t const root, cur_end_type rcur_end = cur_end_left)
	{
		std::map< uint64_t,std::vector<OverlapEntry> > subedges;

		std::set<uint64_t> seen;
		std::stack< std::pair<uint64_t,cur_end_type> > S;
		S.push(std::pair<uint64_t,cur_end_type>(root,rcur_end));

		// extract set of 1d directed edges
		while ( !S.empty() )
		{
			std::pair<uint64_t,cur_end_type> const P = S.top(); S.pop();
			seen.insert(P.first);

			if ( preedges.find(P.first) != preedges.end() )
			{
				std::vector<OverlapEntry> const & ledges = preedges.find(P.first)->second;

				for ( uint64_t i = 0; i < ledges.size(); ++i )
					if ( (P.second == cur_end_left) && isLeftEdge(ledges[i]) )
					{
						subedges[P.first].push_back(ledges[i]);

						switch (ledges[i].orientation)
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
								if ( seen.find(ledges[i].target) == seen.end() )
									S.push(std::pair<uint64_t,cur_end_type>(ledges[i].target,cur_end_right));
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
								if ( seen.find(ledges[i].target) == seen.end() )
									S.push(std::pair<uint64_t,cur_end_type>(ledges[i].target,cur_end_left));
								break;
							default:
								assert(0);
								break;
						}
					}
					else if ( (P.second == cur_end_right) && isRightEdge(ledges[i]) )
					{
						subedges[P.first].push_back(ledges[i]);

						switch (ledges[i].orientation)
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
								if ( seen.find(ledges[i].target) == seen.end() )
									S.push(std::pair<uint64_t,cur_end_type>(ledges[i].target,cur_end_right));
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								if ( seen.find(ledges[i].target) == seen.end() )
									S.push(std::pair<uint64_t,cur_end_type>(ledges[i].target,cur_end_left));
								break;
							default:
								assert(0);
								break;
						}
					}
			}
		}

		return subedges;
	}

	// returns true iff component containing root has a border vertex with only left edges
	static bool hasBorderVertexWithLeftEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root
	)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				(
					countLeftEdges(edges.find(*ita)->second)
					&&
					(!countRightEdges(edges.find(*ita)->second))
				)
			)
				return true;

		return false;
	}

	// returns true iff component containing root has a border vertex with only left edges
	static bool hasBorderVertexWithRightEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root
	)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				(
					countRightEdges(edges.find(*ita)->second)
					&&
					(!countLeftEdges(edges.find(*ita)->second))
				)
			)
				return true;

		return false;
	}

	// return (numerically) smallest left border vertex of component
	static uint64_t getBorderVertexWithLeftEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				countLeftEdges(edges.find(*ita)->second)
				&&
				!countRightEdges(edges.find(*ita)->second)
			)
				return *ita;

		libmaus2::exception::LibMausException ex;
		ex.getStream() << "getBorderVertexWithLeftEdge called for component without border vertices\n";
		ex.finish();
		throw ex;
	}

	// return all left border vertices of component
	static std::set<uint64_t> getBorderVerticesWithLeftEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		std::set<uint64_t> R;
		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				countLeftEdges(edges.find(*ita)->second)
				&&
				!countRightEdges(edges.find(*ita)->second)
			)
				R.insert(*ita);

		return R;
	}

	// return (numerically) smallest right border vertex of component
	static uint64_t getBorderVertexWithRightEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				countRightEdges(edges.find(*ita)->second)
				&&
				!countLeftEdges(edges.find(*ita)->second)
			)
				return *ita;

		libmaus2::exception::LibMausException ex;
		ex.getStream() << "getBorderVertexWithRightEdge called for component without border vertices\n";
		ex.finish();
		throw ex;
	}

	// return all right border vertices of component
	static std::set<uint64_t> getBorderVerticesWithRightEdge(
		std::map< uint64_t,std::vector<OverlapEntry> > const & edges,
		uint64_t const root)
	{
		std::set<uint64_t> component;
		std::stack<uint64_t> S;
		S.push(root);
		component.insert(root);

		while ( ! S.empty() )
		{
			uint64_t const v = S.top();
			S.pop();

			if ( edges.find(v) != edges.end() )
			{
				std::vector<OverlapEntry> const & V = edges.find(v)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( component.find(V[i].target) == component.end() )
					{
						component.insert(V[i].target);
						S.push(V[i].target);
					}
			}
		}

		std::set<uint64_t> R;
		for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			if (
				edges.find(*ita) != edges.end()
				&&
				countRightEdges(edges.find(*ita)->second)
				&&
				!countLeftEdges(edges.find(*ita)->second)
			)
				R.insert(*ita);

		return R;
	}

	// get set of vertices in connected component containing root
	static std::set<uint64_t> getConnectedComponent(std::map< uint64_t,std::vector<OverlapEntry> > const & edges, uint64_t const root)
	{
		// set of vertices
		std::set<uint64_t> R;
		// stack
		std::stack<uint64_t> S;
		S.push(root);
		R.insert(root);

		while ( !S.empty() )
		{
			uint64_t const s = S.top(); S.pop();

			if ( edges.find(s) != edges.end() )
			{
				std::vector<OverlapEntry> const & T = edges.find(s)->second;
				for ( uint64_t i = 0; i < T.size(); ++i )
				{
					OverlapEntry const & OE = T[i];
					if ( R.find(OE.target) == R.end() )
					{
						R.insert(OE.target);
						S.push(OE.target);
					}
				}
			}
		}

		return R;
	}

	void depthFirstLoops(
		std::pair<uint64_t,cur_end_type> const root,
		std::vector< std::pair<uint64_t,cur_end_type> > & S,
		std::set< std::pair<uint64_t,cur_end_type> > & onstack,
		std::vector< std::vector< std::pair<uint64_t,cur_end_type> > > & loops
	)
	{
		std::pair<uint64_t,cur_end_type> const s = S.back();

		if ( preedges.find(s.first) != preedges.end() )
		{
			std::vector<OverlapEntry> const & E = preedges.find(s.first)->second;

			for ( uint64_t i = 0; i < E.size(); ++i )
				if
				(
					(s.second == cur_end_left && isLeftEdge(E[i]))
					||
					(s.second == cur_end_right && isRightEdge(E[i]))
				)
				{
					cur_end_type target_end = cur_end_left;
					switch ( E[i].orientation )
					{
						case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
						case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
							target_end = cur_end_right;
							break;
						case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
						case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
							target_end = cur_end_left;
							break;
						default:
							break;
					}

					std::pair<uint64_t,cur_end_type> const T(E[i].target,target_end);

					if ( T == root )
					{
						S.push_back(root);
						loops.push_back(S);
						S.pop_back();
					}
					else if ( onstack.find(T) == onstack.end() )
					{
						onstack.insert(T);
						S.push_back(T);

						depthFirstLoops(root,S,onstack,loops);

						S.pop_back();
						onstack.erase(onstack.find(T));
					}
				}
		}
	}

	// break all strongly connected components and store deleted edges in removedEdges
	void breakStronglyConnectedComponents(std::map<uint64_t,std::vector<OverlapEntry> > & removedEdges, ReadContainer const * RC = 0)
	{
		uint64_t d_id = 0;
		bool changed = true;

		while ( changed )
		{
			changed = false;

			// get set of all vertices with outgoing edges
			std::set<uint64_t> componentunused;
			for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
				componentunused.insert(ita->first);

			while ( componentunused.size() )
			{
				// choose (numerically) smallest remaining vertex
				uint64_t const root0 = *(componentunused.begin());
				// get component
				std::set<uint64_t> const lcomponent = getConnectedComponent(preedges,root0);

				// mark vertices as used
				for ( std::set<uint64_t>::const_iterator ita = lcomponent.begin(); ita != lcomponent.end(); ++ita )
				{
					assert ( componentunused.find(*ita) != componentunused.end() );
					componentunused.erase(componentunused.find(*ita));
				}

				// try to find a border vertex or use root0 if is is on a circle
				std::vector< std::pair<uint64_t,cur_end_type> > root1V;
				// collect left border edges
				if ( hasBorderVertexWithLeftEdge(preedges,root0) )
				{
					std::set<uint64_t> R = getBorderVerticesWithLeftEdge(preedges,root0);
					for ( std::set<uint64_t>::const_iterator ita = R.begin(); ita != R.end(); ++ita )
						root1V.push_back(
							std::pair<uint64_t,cur_end_type>(*ita,cur_end_left)
						);
				}
				// collect right border edges
				if ( hasBorderVertexWithRightEdge(preedges,root0) )
				{
					std::set<uint64_t> R = getBorderVerticesWithRightEdge(preedges,root0);
					for ( std::set<uint64_t>::const_iterator ita = R.begin(); ita != R.end(); ++ita )
						root1V.push_back(
							std::pair<uint64_t,cur_end_type>(*ita,cur_end_right)
						);
				}
				// if there are no border edges then use any node
				if ( ! root1V.size() )
				{
					root1V.push_back(
						std::pair<uint64_t,cur_end_type>(root0,cur_end_right)
					);
				}

				// iterate over all border edges
				for ( uint64_t zroot1 = 0; zroot1 < root1V.size(); ++zroot1 )
				{
					// extract 1d edges
					uint64_t root1 = root1V[zroot1].first;
					cur_end_type root_cur_end = root1V[zroot1].second;

					std::cerr << "root0=" << root0 << " root1=" << root1 << std::endl;

					#if 0
					std::map< uint64_t,std::vector<OverlapEntry> > subedges = extractConnectedSubEdges(root1,root_cur_end);

					std::map< uint64_t,std::vector<uint64_t> > subedgesshift = extractConnectedSubEdgesShift(root1,root_cur_end);

					if ( ! subedges.size() )
						continue;


					// compute strongly connected components based on 1d edges
					std::pair< std::vector< uint64_t >, std::vector< uint64_t > > SCC =
						libmaus2::graph::StronglyConnectedComponents::strongConnectContract<OverlapEntry,OverlapEntryTargetProjector>(
							subedges,root1
						);
					#endif

					// compute strongly connected components based on 1d edges
					std::pair< std::vector< uint64_t >, std::vector< uint64_t > > SCC = getStronglyConnectedComponents(root1,root_cur_end);

					#if 0
					for ( uint64_t i = 1; i < SCCshift.second.size(); ++i )
						if ( SCCshift.second[i]-SCCshift.second[i-1] > 1 )
						{
							std::cerr << "graph* contains a non-trivial strongly connected component of size " << SCCshift.second[i]-SCCshift.second[i-1] << " ";
						}
					#endif

					// look for non trivial strongly connected components
					for ( uint64_t i = 1; i < SCC.second.size(); ++i )
						if ( SCC.second[i]-SCC.second[i-1] > 1 )
						{
							// std::cerr << "graph contains a non-trivial strongly connected component of size " << SCC.second[i]-SCC.second[i-1] << " ";

							for ( uint64_t j = SCC.second[i-1]; j < SCC.second[i]; ++j )
								std::cerr << SCC.first[j] << ";";

							if ( RC )
							{
								std::set<uint64_t> inside(
									SCC.first.begin() + SCC.second[i-1],
									SCC.first.begin() + SCC.second[i]
								);

								{
									std::ostringstream fnostr;
									fnostr << "loopy_" << d_id++;
									std::string const fnprefix = fnostr.str();
									std::string const dotfn = fnprefix + ".dot";
									libmaus2::aio::OutputStreamInstance COS(dotfn.c_str());
									printSubGraph(COS,inside);
									COS.flush();

									std::ostringstream comstr;
									std::string const svgfn = fnprefix + ".svg.gz";
									comstr << "neato -Tsvg " << dotfn << " | gzip -9 > " << svgfn;
									int const r = system(comstr.str().c_str());
									if ( r != EXIT_SUCCESS )
									{

									}

									remove(dotfn.c_str());
								}

								std::set<uint64_t> outsideConnected;
								for ( uint64_t j = SCC.second[i-1]; j < SCC.second[i]; ++j )
								{
									uint64_t const v = SCC.first[j];
									assert ( preedges.find(v) != preedges.end() );
									std::vector<OverlapEntry> const & E = preedges.find(v)->second;
									for ( uint64_t z = 0; z < E.size(); ++z )
										if ( inside.find(E[z].target) == inside.end() )
											outsideConnected.insert(v);
								}

								if ( outsideConnected.size() )
								{
									std::cerr << " outsideConnected.size()=" << outsideConnected.size() << std::endl;

									for ( std::set<uint64_t>::const_iterator ita = outsideConnected.begin();
										ita != outsideConnected.end(); ++ita )
										std::cerr << *ita << ";";
									std::cerr << " ";

									uint64_t const vroot = *(outsideConnected.begin());

									std::vector< std::vector< std::pair<uint64_t,cur_end_type> > > loops;

									if ( ! loops.size() )
									{
										std::vector< std::pair<uint64_t,cur_end_type> > S(1,std::pair<uint64_t,cur_end_type>(vroot,cur_end_left));
										std::set< std::pair<uint64_t,cur_end_type> > onstack;
										depthFirstLoops(S.front(),S,onstack,loops);
									}
									if ( ! loops.size() )
									{
										std::vector< std::pair<uint64_t,cur_end_type> > S(1,std::pair<uint64_t,cur_end_type>(vroot,cur_end_right));
										std::set< std::pair<uint64_t,cur_end_type> > onstack;
										depthFirstLoops(S.front(),S,onstack,loops);
									}

									std::sort(loops.begin(),loops.end());

									std::cerr << "number of loops " << loops.size() << "\n";

									for ( uint64_t i = 0; i < loops.size(); ++i )
									{
										std::vector<uint64_t> C;
										for ( uint64_t j = 0; j < loops[i].size(); ++j )
										{
											std::cerr << loops[i][j].first << ";";
											C.push_back(loops[i][j].first);
										}
										std::cerr << std::endl;

										std::pair<std::string,double> const contig = createLinearContig(C,*RC);
										std::cerr << contig.first << std::endl;
									}
								}

								#if 0
								std::pair<std::string,double> const contig = createLinearContig(
									std::vector<uint64_t>(SCC.first.begin() + SCC.second[i-1],SCC.first.begin()+SCC.second[i]),
									*RC);
								std::cerr << " " << contig.first;
								#endif
							}

							// extract the strongly connected sub component
							std::set<uint64_t> subcomp(SCC.first.begin() + SCC.second[i-1],SCC.first.begin() + SCC.second[i]);

							// minimum number of left and right edges anywhere in the component
							uint64_t minleftcnt = std::numeric_limits<uint64_t>::max();
							uint64_t minrightcnt = std::numeric_limits<uint64_t>::max();

							for ( std::set<uint64_t>::const_iterator ita = subcomp.begin(); ita != subcomp.end(); ++ita )
							{
								assert ( preedges.find(*ita) != preedges.end() );
								std::vector<OverlapEntry> const & V = preedges.find(*ita)->second;

								uint64_t lc = 0, rc = 0;
								for ( uint64_t j = 0; j < V.size(); ++j )
									if ( subcomp.find(V[j].target) != subcomp.end() )
									{
										if ( isLeftEdge(V[j]) && (subcomp.find(V[j].target) != subcomp.end()) )
											++lc;
										else if ( isRightEdge(V[j]) && (subcomp.find(V[j].target) != subcomp.end()) )
											++rc;
									}

								minleftcnt = std::min(lc,minleftcnt);
								minrightcnt = std::min(rc,minrightcnt);
							}

							// get minimum number of left or right edges with maximum outside links
							uint64_t maxoutlinks = 0;
							uint64_t maxoutlinknode = 0;
							for ( std::set<uint64_t>::const_iterator ita = subcomp.begin(); ita != subcomp.end(); ++ita )
							{
								assert ( preedges.find(*ita) != preedges.end() );
								std::vector<OverlapEntry> const & V = preedges.find(*ita)->second;

								uint64_t lc = 0, rc = 0;
								for ( uint64_t j = 0; j < V.size(); ++j )
									if ( subcomp.find(V[j].target) != subcomp.end() )
									{
										if ( isLeftEdge(V[j])  && (subcomp.find(V[j].target) != subcomp.end()) )
											++lc;
										else if ( isRightEdge(V[j])  && (subcomp.find(V[j].target) != subcomp.end()) )
											++rc;
									}

								uint64_t outcnt = 0;
								for ( uint64_t j = 0; j < V.size(); ++j )
									if ( subcomp.find(V[j].target) == subcomp.end() )
										++outcnt;

								if (
									((minleftcnt <= minrightcnt && lc == minleftcnt)
									||
									(minrightcnt < minleftcnt && rc == minrightcnt))
									&&
									( outcnt >= maxoutlinks )
								)
								{
									maxoutlinks = outcnt;
									maxoutlinknode = *ita;
								}
							}

							// count number of left and right outgoing edges to outside of component
							uint64_t outleft = 0;
							uint64_t outright = 0;
							assert ( preedges.find(maxoutlinknode) != preedges.end() );
							std::vector<OverlapEntry> const & V = preedges.find(maxoutlinknode)->second;
							for ( uint64_t j = 0; j < V.size(); ++j )
								if ( (subcomp.find(V[j].target) == subcomp.end()) && isLeftEdge(V[j]) )
									outleft++;
								else if ( (subcomp.find(V[j].target) == subcomp.end()) && isRightEdge(V[j]) )
									outright++;

							// count number of left and right incoming edges from inside sub component
							std::vector< std::pair<uint64_t,OverlapEntry> > const RIV = getReverseIncoming(maxoutlinknode);
							uint64_t inleft = 0;
							uint64_t inright = 0;
							for ( uint64_t j = 0; j < RIV.size(); ++j )
								if ( subcomp.find(RIV[j].first) != subcomp.end() )
									switch ( RIV[j].second.orientation )
									{
										case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
										case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
											inleft++;
											break;
										case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
										case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
											inright++;
											break;
										default:
											break;
								}

							// std::cerr << "inleft=" << inleft << " inright=" << inright << " outleft=" << outleft << " outright=" << outright << std::endl;

							#if 0
							for ( uint64_t j = 0; j < V.size(); ++j )
								std::cerr << V[j] << " " << ((subcomp.find(V[j].target) != subcomp.end())?'+':'-') << std::endl;
							#endif

							// remove all edges into component (and their reverse edges) for chosen vertex in either left or right direction
							std::vector<uint64_t> deltargets;
							if ( (inleft * outright) >= (inright * outleft) )
							{
								for ( uint64_t j = 0; j < V.size(); ++j )
									if (
										(subcomp.find(V[j].target) != subcomp.end())
										&&
										isRightEdge(V[j])
									)
										deltargets.push_back(V[j].target);
							}
							// either we did not use the branch above or it yielded no edges
							if ( deltargets.size() == 0 )
							{
								for ( uint64_t j = 0; j < V.size(); ++j )
									if (
										(subcomp.find(V[j].target) != subcomp.end())
										&&
										isLeftEdge(V[j])
									)
										deltargets.push_back(V[j].target);
							}
							// either we did not use the branches above or they yielded no edges
							if ( deltargets.size() == 0 )
							{
								for ( uint64_t j = 0; j < V.size(); ++j )
									if (
										(subcomp.find(V[j].target) != subcomp.end())
										&&
										isRightEdge(V[j])
									)
										deltargets.push_back(V[j].target);
							}

							// std::cerr << "number of edge targets removed is " << deltargets.size() << std::endl;

							for ( uint64_t j = 0; j < deltargets.size(); ++j )
							{
								uint64_t const src = maxoutlinknode;
								uint64_t const tgt = deltargets[j];

								if ( hasEdge(src,tgt) )
								{
									removedEdges[src].push_back(getEdge(src,tgt));
									removeEdge(src,tgt);
								}
								if ( hasEdge(tgt,src) )
								{
									removedEdges[tgt].push_back(getEdge(tgt,src));
									removeEdge(tgt,src);
								}

								std::cerr << " rm edge " << src << " <-> " << tgt;

								changed = true;
							}

							std::cerr << std::endl;
						}
				}
			}
		}
	}

	void disconnectTips(std::map<uint64_t,std::vector<OverlapEntry> > & removedEdges)
	{
		// remove long tips
		bool changed = true;

		while ( changed )
		{
			changed = false;

			// remove tips (of any length) from contig until there are no more tips
			bool foundTip = false;

			std::set<uint64_t> component;
			for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
				component.insert(ita->first);

			for ( std::set<uint64_t>::const_iterator ita = component.begin(); (!foundTip) && (ita != component.end()); ++ita )
			{
				// vertex with edges in one direction only
				if ( isTipSeed(*ita) )
				{
					// extract the tip (if any)
					std::deque<uint64_t> C;
					FollowUniqueResult FUR;
					foundTip = extractTip(*ita, C, FUR);

					if ( foundTip )
					{
						// std::cerr << "Found tip at " << *ita << std::endl;

						// get the vertex of the split (one beyond the end of the tip)
						uint64_t tipsplitvertex = 0;
						OverlapEntry revedge;

						if ( FUR.cur_end == cur_end_left )
						{
							assert ( countLeftEdges(FUR.v) == 1 );
							tipsplitvertex = getFirstLeftEdge(FUR.v).target;
						}
						else
						{
							assert ( countRightEdges(FUR.v) == 1 );
							tipsplitvertex = getFirstRightEdge(FUR.v).target;
						}

						std::cerr << "found tip at " << *ita << " tipsplitvertex " << tipsplitvertex << std::endl;

						assert ( hasEdge(FUR.v,tipsplitvertex) );
						assert ( hasReverseEdge(FUR.v,tipsplitvertex) );
						// get the reverse edge from the split vertex to the end of the tip
						revedge = getEdge(tipsplitvertex,FUR.v);

						// std::cerr << "revedge " << revedge << std::endl;

						// get direction of that edge
						cur_end_type revcurend = cur_end_left;
						switch ( revedge.orientation )
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								revcurend = cur_end_right;
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
								revcurend = cur_end_left;
								break;
							default:
								break;
						}

						// get set of all edges from splitting vertex
						std::vector<OverlapEntry> const & revedges = preedges.find(tipsplitvertex)->second;

						// std::cerr << "number of revedges is " << revedges.size() << std::endl;

						std::vector< std::pair<uint64_t,uint64_t> > tipSeedVector;

						for ( uint64_t i = 0; i < revedges.size(); ++i )
							if (
								(revcurend == cur_end_left && isLeftEdge(revedges[i]))
								||
								(revcurend == cur_end_right && isRightEdge(revedges[i]))
							)
							{
								OverlapEntry const & edge = revedges[i];

								// std::cerr << "revedges[]=" << edge << std::endl;

								FollowUniqueResult subFUR;
								std::deque<uint64_t> subC;

								switch ( edge.orientation )
								{
									case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
									case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
										subC.push_back(edge.target);
										subFUR = followUniqueRightRecord(edge.target,subC);
										break;
									case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
									case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
										subC.push_back(edge.target);
										subFUR = followUniqueLeftRecord(edge.target,subC);
										break;
									default:
										break;
								}

								// std::cerr << "subC.size()=" << subC.size() << std::endl;

								assert ( preedges.find(subFUR.v) != preedges.end() );

								bool istip = false;

								// if end is at back of subC
								if (
									subFUR.v == subC.back()
									&&
									subC.size() > 1
								)
								{
									assert ( hasEdge(subC[subC.size()-2],subC[subC.size()-1]) );
									assert ( hasEdge(subC[subC.size()-1],subC[subC.size()-2]) );
									assert ( hasReverseEdge(subC[subC.size()-2],subC[subC.size()-1]) );
									assert ( hasReverseEdge(subC[subC.size()-1],subC[subC.size()-2]) );

									if (
										isLeftEdge(getEdge(subC[subC.size()-1],subC[subC.size()-2]))
										&&
										countLeftEdges(subC[subC.size()-1]) == 1
										&&
										countRightEdges(subC[subC.size()-1]) == 0
									)
									{
										istip = true;
									}
									else if (
										isRightEdge(getEdge(subC[subC.size()-1],subC[subC.size()-2]))
										&&
										countRightEdges(subC[subC.size()-1]) == 1
										&&
										countLeftEdges(subC[subC.size()-1]) == 0

									)
									{
										istip = true;
									}
								}
								// otherwise end is at front
								else if (
									subFUR.v == subC.front()
									&&
									subC.size() > 1
								)
								{
									assert ( subC[0] == subFUR.v );
									assert ( hasEdge(subC[1],subC[0]) );
									assert ( hasEdge(subC[0],subC[1]) );
									assert ( hasReverseEdge(subC[1],subC[0]) );
									assert ( hasReverseEdge(subC[0],subC[1]) );

									if (
										isLeftEdge(getEdge(subC[0],subC[1])) &&
										countLeftEdges(subC[0]) == 1 &&
										countRightEdges(subC[0]) == 0
									)
									{
										istip = true;
									}
									else if (
										isRightEdge(getEdge(subC[0],subC[1])) &&
										countRightEdges(subC[0]) == 1 &&
										countLeftEdges(subC[0]) == 0
									)
									{
										istip = true;
									}
								}
								else
								{
									assert ( subC.size() == 1 );
									assert ( hasEdge(subC[0],tipsplitvertex) );
									assert ( hasEdge(tipsplitvertex,subC[0]) );

									if (
										isLeftEdge(getEdge(subC[0],tipsplitvertex)) &&
										countLeftEdges(subC[0]) == 1 &&
										countRightEdges(subC[0]) == 0
									)
									{
										istip = true;
									}
									else if (
										isRightEdge(getEdge(subC[0],tipsplitvertex)) &&
										countRightEdges(subC[0]) == 1 &&
										countLeftEdges(subC[0]) == 0
									)
									{
										istip = true;
									}
								}

								if ( istip )
								{
									tipSeedVector.push_back(std::pair<uint64_t,uint64_t>(subFUR.v,subC.size()));
								}
							}

						std::sort(tipSeedVector.begin(),tipSeedVector.end(),PairSecondComparator<uint64_t,uint64_t,std::greater>());

						assert ( tipSeedVector.size() );

						// check whether all tip seeds actually are tips
						for ( uint64_t i = 0; foundTip && i < tipSeedVector.size(); ++i )
						{
							std::deque<uint64_t> exC;
							FollowUniqueResult exFUR;
							bool const tipok = extractTip(tipSeedVector[i].first, exC, exFUR);

							if ( ! tipok )
								foundTip = false;
						}

						if ( foundTip )
						{
							// remove all tips but the longest one if there are multiple tips
							// if there is only one then remove the tip (means the tip is "off" a repeat)
							for ( uint64_t i = (tipSeedVector.size() > 1) ? 1 : 0; i < tipSeedVector.size(); ++i )
							{
								std::deque<uint64_t> exC;
								FollowUniqueResult exFUR;
								bool const tipok = extractTip(tipSeedVector[i].first, exC, exFUR);

								if ( ! tipok )
								{
									std::cerr << "extractTip failed for " << tipSeedVector[i].first << std::endl;
									return;
									// assert ( tipok );
								}

								uint64_t const linkvert = exFUR.v;

								std::cerr << "disconnecting tip of length " << tipSeedVector[i].second << " starting at " << linkvert << " from split node " << tipsplitvertex << std::endl;

								if ( hasEdge(tipsplitvertex,linkvert) )
								{
									removedEdges[tipsplitvertex].push_back(getEdge(tipsplitvertex,linkvert));
									removeEdge(tipsplitvertex,linkvert);
								}
								if ( hasEdge(linkvert,tipsplitvertex) )
								{
									removedEdges[linkvert].push_back(getEdge(linkvert,tipsplitvertex));
									removeEdge(linkvert,tipsplitvertex);
								}

								#if 0
								std::cerr << "Removing tip at " << tipSeedVector[i].first << " length " << tipSeedVector[i].second << " branching from " << tipsplitvertex << std::endl;
								removeEdges(std::vector<uint64_t>(exC.begin(),exC.end()));
								for ( uint64_t i = 0; i < exC.size(); ++i )
									component.erase(component.find(exC[i]));
								#endif

								changed = true;
							}
						}
					}
				}
			}
		}

	}

	// compute a linear contig based on a simple consensus of the reads
	std::pair<std::string,double> createLinearContig(std::vector<uint64_t> const & vertices, ReadContainer const & RC) const
	{
		if ( ! vertices.size() )
			return std::pair<std::string,double>(std::string(),0);
		else if ( vertices.size() == 1 )
			return std::pair<std::string,double>(RC.getRead(vertices[0]),1);
		else
		{
			assert ( preedges.find(vertices[0]) != preedges.end() );
			assert ( hasEdge(vertices[0],vertices[1]) );
			OverlapEntry const firstedge = getEdge(vertices[0],vertices[1]);
			cur_end_type cur_end = isLeftEdge(firstedge) ? cur_end_left : cur_end_right;
			#if 0
			std::vector<OverlapEntry> const & V0 = preedges.find(vertices[0])->second;
			assert ( V0.size() != 0 );
			assert ( V0.size() == 1 );
			#endif

			std::ostringstream ostr;
			Pileup pile;

			// length of previous read
			uint64_t prevreadlen = RC.getRead(vertices[0]).size();
			// uint64_t prevoverhang = 0;
			uint64_t offset = 0;
			uint64_t support = 0;

			// add first read to pile
			switch ( cur_end )
			{
				case cur_end_left:
					// std::cerr << "adding " << libmaus2::fastx::reverseComplementUnmapped(RC.getRead(vertices[0])) << " to pile at offset " << 0 << " for vertex rc " << vertices[0] << std::endl;
					pile.add ( 0, libmaus2::fastx::reverseComplementUnmapped(RC.getRead(vertices[0])) );
					break;
				case cur_end_right:
					// std::cerr << "adding " << (RC.getRead(vertices[0])) << " to pile at offset " << 0 << " for vertex " << vertices[0] << std::endl;
					pile.add ( 0, RC.getRead(vertices[0]) );
					break;
			}

			// add rest of the reads on the path
			for ( uint64_t i = 1; i < vertices.size(); ++i )
			{
				// source vertex
				uint64_t const sourceid = vertices[i-1];
				// connecting edge
				OverlapEntry const & edge = getEdge(sourceid,vertices[i]);
				// sanity check
				assert ( edge.target == vertices[i] );
				assert (
					(cur_end == cur_end_left && isLeftEdge(edge))
					||
					(cur_end == cur_end_right && isRightEdge(edge))
				);
				// overhang over previous read
				// uint64_t const overhang = edge.overhang;
				// overlap between the two reads
				uint64_t const overlap  = edge.overlap;

				// compute consensus for bases which are no longer covered by the current read
				for ( uint64_t j = 0; j < prevreadlen-overlap; ++j )
				{
					std::pair<char,uint64_t> const cons = pile.consensus(offset++);
					ostr.put ( cons.first );
					support += cons.second;
				}

				// add new read to pile and set current extension end according to edge
				switch ( edge.orientation )
				{
					case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
					case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
						// ostr << RC.getRead(vertices[i]).substr(overlap,overhang);
						// std::cerr << "adding " << (RC.getRead(vertices[i])) << " to pile at offset " << offset << " for vertex " << vertices[i] << std::endl;
						pile.add(offset,RC.getRead(vertices[i]));
						cur_end = cur_end_right;
						break;
					case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
					case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
						// ostr << libmaus2::fastx::reverseComplementUnmapped(RC.getRead(vertices[i])).substr(overlap,overhang);
						// std::cerr << "adding " << libmaus2::fastx::reverseComplementUnmapped(RC.getRead(vertices[i])) << " to pile at offset " << offset << " for vertex rc " << vertices[i] << std::endl;
						pile.add(offset,libmaus2::fastx::reverseComplementUnmapped(RC.getRead(vertices[i])));
						cur_end = cur_end_left;
						break;
					default:
						break;
				}

				// length of previous read
				prevreadlen = RC.getRead(vertices[i]).size();
				// overhang of previous read
				// prevoverhang = overhang;
			}

			// process remaining data
			// for ( uint64_t j = 0; j < prevoverhang; ++j )
			for ( uint64_t j = 0; j < prevreadlen; ++j )
			{
				std::pair<char,uint64_t> const cons = pile.consensus(offset++);
				ostr.put( cons.first );
				support += cons.second;
			}

			return std::pair<std::string,double>(ostr.str(),static_cast<double>(support)/ostr.str().size() );
		}
	}

	bool isLeftTipSeed(uint64_t const i) const
	{
		// does vertex have any edges?
		if ( preedges.find(i) == preedges.end() )
			return false;

		// get edge list
		std::vector<OverlapEntry> const & V = preedges.find(i)->second;

		// not a tip seed if it has more than one edge
		if ( V.size() != 1 )
			return false;

		// is the edge a left edge?
		switch ( V[0].orientation )
		{
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
				return true;
			default:
				return false;
		}
	}

	bool isRightTipSeed(uint64_t const i) const
	{
		if ( preedges.find(i) == preedges.end() )
			return false;

		std::vector<OverlapEntry> const & V = preedges.find(i)->second;

		if ( V.size() != 1 )
			return false;

		switch ( V[0].orientation )
		{
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
				return true;
			default:
				return false;
		}
	}

	bool isTipSeed(uint64_t const i) const
	{
		return isLeftTipSeed(i) || isRightTipSeed(i);
	}

	bool hasEdge(uint64_t const from, uint64_t const to) const
	{
		if ( preedges.find(from) == preedges.end() )
			return false;

		std::vector<OverlapEntry> const & V = preedges.find(from)->second;

		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i].target == to )
				return true;

		return false;
	}

	OverlapEntry getEdge(uint64_t const from, uint64_t const to) const
	{
		if ( ! hasEdge(from,to) )
		{
			std::ostringstream ostr;
			ostr << "getEdge called for non-existant edge " << from << " -> " << to << std::endl;

			if ( preedges.find(from) != preedges.end() )
			{
				std::vector<OverlapEntry> const & V = preedges.find(from)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					ostr << from << " -> " << V[i] << std::endl;;
			}
			if ( preedges.find(to) != preedges.end() )
			{
				std::vector<OverlapEntry> const & V = preedges.find(to)->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					ostr << to << " -> " << V[i] << std::endl;;
			}

			fail(ostr.str());
		}

		std::vector<OverlapEntry> const & V = preedges.find(from)->second;

		for ( uint64_t i = 0; i < V.size(); ++i )
			if ( V[i].target == to )
				return V[i];

		fail("getEdge called for non-existent edge");
		std::terminate();
	}

	std::vector < std::vector<uint64_t> > findShortContigs(uint64_t const maxfollow = std::numeric_limits<uint64_t>::max()) const
	{
		std::vector < std::vector<uint64_t> > R;
		std::set<uint64_t> removed;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
		{
			if ( removed.find(ita->first) == removed.end() )
			{
				// set of vertices
				std::deque<uint64_t> C;
				// information about last vertex on unique chain
				FollowUniqueResult FUR;

				if ( isLeftTipSeed(ita->first) )
				{
					C.push_back(ita->first);
					FUR = followUniqueLeftRecord(ita->first,C,maxfollow);
				}
				else if ( isRightTipSeed(ita->first) )
				{
					C.push_back(ita->first);
					FUR = followUniqueRightRecord(ita->first,C,maxfollow);
				}

				if (
					C.size()
					&&
					C.size() <= maxfollow
					&&
					(
						(
							FUR.cur_end == cur_end_left
							&&
							countLeftEdges(FUR.v) == 0
						)
						||
						(
							FUR.cur_end == cur_end_right
							&&
							countRightEdges(FUR.v) == 0
						)
					)
				)
				{
					for ( uint64_t i = 0; i < C.size(); ++i )
						removed.insert(C[i]);

					R.push_back(std::vector<uint64_t>(C.begin(),C.end()));
				}
			}

		}

		return R;
	}

	std::vector < std::vector<uint64_t> > findStraightContigs() const
	{
		std::vector < std::vector<uint64_t> > R;
		std::set<uint64_t> seen;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		{
			uint64_t const v = ita->first;

			if ( seen.find(v) == seen.end() )
			{
				// set of vertices
				std::deque<uint64_t> C;

				// edges only on one side
				if (
					(countLeftEdges(v) == 1 && countRightEdges(v) == 0)
					||
					(countLeftEdges(v) == 0 && countRightEdges(v) == 1)
				)
				{
					cur_end_type cur_end = countLeftEdges(v) ? cur_end_left : cur_end_right;

					uint64_t vv = v;
					C.push_back(vv);

					while (
						(cur_end == cur_end_left  && countLeftEdges(vv) == 1  && countRightEdges(vv) <= 1)
						||
						(cur_end == cur_end_right && countRightEdges(vv) == 1 && countLeftEdges(vv) <=1)
					)
					{
						OverlapEntry const OE = (cur_end == cur_end_left) ? getFirstLeftEdge(vv) : getFirstRightEdge(vv);

						switch ( OE.orientation )
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								cur_end = cur_end_left;
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
								cur_end = cur_end_right;
								break;
							default:
								break;
						}

						vv = OE.target;
						C.push_back(vv);
					}

					for ( uint64_t i = 0; i < C.size(); ++i )
						seen.insert(C[i]);

					if (
						(cur_end == cur_end_left && countLeftEdges(vv) == 0)
						||
						(cur_end == cur_end_right && countRightEdges(vv) == 0)
					)
					{
						R.push_back(std::vector<uint64_t>(C.begin(),C.end()));
					}
				}
			}
		}

		return R;
	}

	static std::set<uint64_t> flipLowBit(std::set<uint64_t> const & S)
	{
		std::set<uint64_t> R;
		for ( std::set<uint64_t>::const_iterator ita = S.begin(); ita != S.end(); ++ita )
			R.insert ( *ita ^ 1 );
		return R;
	}

	static std::map< uint64_t,std::vector<uint64_t> > getReverseEdgeSet(std::map< uint64_t,std::vector<uint64_t> > const & M)
	{
		std::map< uint64_t,std::vector<uint64_t> > R;

		for ( std::map< uint64_t,std::vector<uint64_t> >::const_iterator ita = M.begin(); ita != M.end(); ++ita )
		{
			uint64_t const v = ita->first;
			std::vector<uint64_t> const & V = ita->second;

			for ( uint64_t i = 0; i < V.size(); ++i )
				R[V[i]].push_back(v);
		}

		return R;
	}

	void findBubbles(std::map<uint64_t,std::vector<OverlapEntry> > & /* removedEdges */, ReadContainer const * RC = 0, std::string const * prefix = 0)
	{
		// set of unused vertices
		std::set<uint64_t> unused;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
			unused.insert(ita->first);


		// component id
		for ( uint64_t componentid = 0 ; unused.size(); ++componentid )
		{
			// get next (minimum id) unused node
			uint64_t const root = *(unused.begin());
			std::stack<uint64_t> S;
			std::set<uint64_t> component;
			S.push(root);
			component.insert(root);

			// get connected component in undirected graph
			while ( !S.empty() )
			{
				uint64_t const v = S.top();
				S.pop();

				if ( preedges.find(v) != preedges.end() )
				{
					std::vector<OverlapEntry> const & V = preedges.find(v)->second;
					for ( uint64_t i = 0; i < V.size(); ++i )
						if ( component.find(V[i].target) == component.end() )
						{
							S.push(V[i].target);
							component.insert(V[i].target);
						}
				}
			}

			// erase all vertices in the component
			for ( std::set<uint64_t>::const_iterator ita = component.begin(); ita != component.end(); ++ita )
			{
				assert ( unused.find(*ita) != unused.end() );
				unused.erase(unused.find(*ita));
			}

			// extract 1d edges
			uint64_t const root0 = *(component.begin());
			#if 0
			uint64_t root1;
			cur_end_type root_cur_end;
			#endif
			std::map< uint64_t,std::vector<OverlapEntry> > subedges;

			#if 0
			if ( hasBorderVertexWithLeftEdge(preedges,root0) )
			{
				root1 = getBorderVertexWithLeftEdge(preedges,root0);
				root_cur_end = cur_end_left;
			}
			else if ( hasBorderVertexWithRightEdge(preedges,root0) )
			{
				root1 = getBorderVertexWithRightEdge(preedges,root0);
				root_cur_end = cur_end_right;
			}
			else
			{
				root1 = root0;
				root_cur_end = cur_end_left;
			}
			#endif

			// try to find a border vertex or use root0 if is is on a circle
			std::vector< std::pair<uint64_t,cur_end_type> > root1V;
			// collect left border edges
			if ( hasBorderVertexWithLeftEdge(preedges,root0) )
			{
				std::set<uint64_t> R = getBorderVerticesWithLeftEdge(preedges,root0);
				for ( std::set<uint64_t>::const_iterator ita = R.begin(); ita != R.end(); ++ita )
					root1V.push_back(
						std::pair<uint64_t,cur_end_type>(*ita,cur_end_left)
					);
			}
			// collect right border edges
			if ( hasBorderVertexWithRightEdge(preedges,root0) )
			{
				std::set<uint64_t> R = getBorderVerticesWithRightEdge(preedges,root0);
				for ( std::set<uint64_t>::const_iterator ita = R.begin(); ita != R.end(); ++ita )
					root1V.push_back(
						std::pair<uint64_t,cur_end_type>(*ita,cur_end_right)
					);
			}
			// if there are no border edges then use any node
			if ( ! root1V.size() )
			{
				root1V.push_back(
					std::pair<uint64_t,cur_end_type>(root0,cur_end_right)
				);
			}

			std::cerr << std::string(80,'-') << std::endl;
			std::set < std::set<uint64_t> > vertsetsseen;

			for ( uint64_t root1Vi = 0; root1Vi < root1V.size(); ++root1Vi )
			{
				// id of bubble
				uint64_t bubbleid = 0;

				std::pair<uint64_t,cur_end_type> const root = root1V[root1Vi];

				std::cerr << "root " << root.first << "," << root.second << std::endl;

				// extract component of end annotated vertices
				std::map< uint64_t,std::vector<uint64_t> > const subedges = extractConnectedSubEdgesShift(root.first,root.second);
				std::map< uint64_t,std::vector<uint64_t> > const rsubedges = getReverseEdgeSet(subedges);

				std::set<uint64_t> vertset;
				for ( std::map< uint64_t,std::vector<uint64_t> >::const_iterator ita = subedges.begin(); ita != subedges.end(); ++ita )
				{
					vertset.insert(ita->first);
					std::vector<uint64_t> const & ledges = ita->second;
					for ( uint64_t i = 0; i < ledges.size(); ++i )
						vertset.insert(ledges[i]);
				}

				if ( vertsetsseen.find(flipLowBit(vertset)) == vertsetsseen.end() )
				{
					vertsetsseen.insert(vertset);

					// compute topological sorting starting from cnode
					std::pair<bool,std::map<uint64_t,uint64_t> > const TS = libmaus2::graph::TopologicalSorting::topologicalSorting(
						subedges,
						computeId(root),
						libmaus2::graph::IdentityTargetProjector()
					);

					if ( TS.first )
					{
						typedef std::pair<uint64_t,uint64_t> u8p;
						std::stack<u8p> S;
						S.push(u8p(computeId(root),0));
						std::vector< std::vector<uint64_t> > allpaths(1);

						while ( !S.empty() )
						{
							u8p const cur = S.top();
							S.pop();

							allpaths[cur.second].push_back(cur.first);

							if ( subedges.find(cur.first) != subedges.end() && subedges.find(cur.first)->second.size() )
							{
								std::vector<uint64_t> const & E = subedges.find(cur.first)->second;

								// continue on path
								S.push(u8p(E[0],cur.second));

								for ( uint64_t i = 1; i < E.size(); ++i )
								{
									uint64_t const t = E[i];
									// insert copy of path
									std::vector<uint64_t> const cpath = allpaths[cur.second];
									uint64_t nid = allpaths.size();
									allpaths.push_back(cpath);
									S.push(u8p(t,nid));
								}
							}
						}

						struct PathEndComparator
						{
							bool operator()(std::vector<uint64_t> const & A, std::vector<uint64_t> const & B) const
							{
								if ( ! A.size() && ! B.size() )
									return false;
								else if ( ! A.size() )
									return true;
								else if ( ! B.size() )
									return false;
								else
									return A.back() < B.back();
							}
						};


						std::sort(allpaths.begin(),allpaths.end(),PathEndComparator());

						std::cerr << "number of paths is " << allpaths.size() << std::endl;

						if ( RC )
						{
							uint64_t nextbubbleid = 0;

							uint64_t plow = 0;
							while ( plow != allpaths.size() )
							{
								uint64_t phigh = plow;
								while ( phigh != allpaths.size() && allpaths[phigh].back() == allpaths[plow].back() )
									++phigh;

								uint64_t minlen = std::numeric_limits<uint64_t>::max();
								uint64_t maxlen = std::numeric_limits<uint64_t>::min();
								for ( uint64_t i = plow; i < phigh; ++i )
								{
									std::vector<uint64_t> nv = allpaths[i];
									for ( uint64_t j = 0; j < nv.size(); ++j )
										nv[j] >>= 1;
									std::pair<std::string,double> const contig = createLinearContig(nv,*RC);
									minlen = std::min(minlen,contig.first.size());
									maxlen = std::max(maxlen,contig.first.size());
								}

								std::ostringstream indelstr;
								if ( maxlen != minlen )
									indelstr << "indel_" << maxlen-minlen << "_";

								std::ostringstream fnostr;
								if ( prefix )
									fnostr << *prefix << '_';
								fnostr << "bubbles_"
									<< ((phigh-plow>1)?"multi_":"single_")
									<< indelstr.str()
									<< componentid << "_" << root1Vi << "_" << bubbleid << ".fasta";

								libmaus2::aio::OutputStreamInstance COS(fnostr.str());

								std::cerr << std::string(80,'*') << std::endl;
								for ( uint64_t i = plow; i < phigh; ++i )
								{
									uint64_t const bubbleid = nextbubbleid++;

									#if 0
									std::cerr << "path[" << i << "]=";
									for ( uint64_t j = 0; j < allpaths[i].size(); ++j )
										std::cerr << allpaths[i][j] << ";";
									std::cerr << std::endl;
									#endif

									std::vector<uint64_t> nv = allpaths[i];
									for ( uint64_t j = 0; j < nv.size(); ++j )
										nv[j] >>= 1;

									#if 0
									for ( uint64_t j = 0; j < nv.size(); ++j )
										std::cerr << RC->getRead(j).size() << ";";
									std::cerr << std::endl;

									for ( uint64_t j = 0; j < nv.size(); ++j )
										std::cerr << RC->getRead(nv[j]) << "\n";
									#endif

									std::pair<std::string,double> const contig = createLinearContig(nv,*RC);
									std::ostringstream ostr;
									ostr << '>' << componentid << '_' << root1Vi << '_' << bubbleid << ' ' << contig.first.size() << ' ' << contig.second << '\n';
									ostr << contig.first << std::endl;

									COS << ostr.str();
								}

								COS.flush();

								plow = phigh;
							}
						}
					}
					else
					{
						std::cerr << "findBubbles: graph starting from " << root.first << "," << root.second << " contains loops" << std::endl;
					}

					#if 0
					std::pair<uint64_t,cur_end_type> cnode = root;

					while (
						subedges.find(computeId(cnode)) != subedges.end()
						&&
						subedges.find(computeId(cnode))->second.size()
					)
					{
						std::vector<uint64_t> const & pedges = subedges.find(computeId(cnode))->second;
						assert ( pedges.size() > 0 );

						if ( pedges.size() == 1 )
						{
							// follow linear chain
							std::vector< std::pair<uint64_t,cur_end_type> > lin;

							lin.push_back(cnode);

							while (
								subedges.find(computeId(cnode)) != subedges.end() &&
								subedges.find(computeId(cnode))->second.size()==1
							)
							{
								cnode = decodeId(subedges.find(computeId(cnode))->second.front());
								lin.push_back(cnode);
							}

							std::cerr << "linear chain of length " << lin.size() << std::endl;
						}
						else
						{
							// compute topological sorting starting from cnode
							std::pair<bool,std::map<uint64_t,uint64_t> > const TS = libmaus2::graph::TopologicalSorting::topologicalSorting(
								subedges,
								computeId(cnode),
								libmaus2::graph::IdentityTargetProjector()
							);

							if ( ! TS.first )
							{
								std::cerr << "Graph contains strongly connected components in findBubbles, root=" << root.first << "," << root.second << std::endl;
								break;
							}

							std::map<uint64_t,uint64_t> visitcnt;
							std::set<uint64_t> endnodes;
							std::stack<uint64_t> S;
							S.push(computeId(cnode));

							while ( S.size() )
							{
								uint64_t const node = S.top();
								S.pop();

								visitcnt[node]++;

								if (
									subedges.find(node) != subedges.end()
									&&
									subedges.find(node)->second.size() != 0
								)
								{
									std::vector<uint64_t> const & ledges = subedges.find(node)->second;
									for ( uint64_t i = 0; i < ledges.size(); ++i )
										S.push(ledges[i]);
								}
								else
								{
									endnodes.insert(node);
								}
							}

							std::cerr << "number of end nodes " << endnodes.size() << std::endl;

							#if 0
							struct QueueElement
							{
								uint64_t v;
								uint64_t r;
								uint64_t pathid;

								QueueElement() : v(0), r(0), pathid(0) {}
								QueueElement(uint64_t const rv, uint64_t const rr, uint64_t const rpathid) : v(rv), r(rr), pathid(rpathid) {}

								bool operator<(QueueElement const & O) const
								{
									return r > O.r;
								}
							};

							std::priority_queue < QueueElement > Q;

							std::map<uint64_t,uint64_t> seencnt;
							uint64_t activecnt;

							for ( uint64_t i = 0; i < pedges.size(); ++i )
							{
								Q.push(
									QueueElement(
										pedges[i],
										TS.second.find(pedges[i])->second,
										i
									)
								);
								activecnt++;
							}

							while ( Q.size() )
							{
								QueueElement const QE = Q.top();
								Q.pop();

								std::vector<uint64_t> const & ledges = subedges.find(QE.v)->second;
							}
							#endif

							break;
						}
					}
					#endif
				}
				else
				{
					std::cerr << "vertex set already seen" << std::endl;
				}
			}

			#if 0
			assert ( root1V.size() );
			uint64_t root1 = root1V.front().first;
			cur_end_type root_cur_end = root1V.front().second;

			subedges = extractConnectedSubEdges(root1,root_cur_end);

			// check that edges are all from edges in the component to edges in the component
			for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = subedges.begin(); ita != subedges.end(); ++ita )
			{
				assert ( component.find(ita->first) != component.end() );
				std::vector<OverlapEntry> const & V = ita->second;
				for ( uint64_t i = 0; i < V.size(); ++i )
					assert ( component.find(V[i].target) != component.end() );
			}

			if ( subedges.size() == 1 )
				continue;

			cur_end_type cur_end = root_cur_end;
			uint64_t curvert = root1;

			std::vector <
				std::vector < std::vector < std::pair<uint64_t,cur_end_type> > >
			> variantseqs;

			// while there are any outgoing edges in the relevant direction
			while (
				subedges.find(curvert) != subedges.end()
				&&
				(
					(cur_end == cur_end_left  && countLeftEdges (subedges.find(curvert)->second))
					||
					(cur_end == cur_end_right && countRightEdges(subedges.find(curvert)->second))
				)
			)
			{
				std::vector < std::pair<uint64_t,cur_end_type> > linearStrip;
				linearStrip.push_back(std::pair<uint64_t,cur_end_type>(curvert,cur_end));

				// follow linear chain until continuation is non existant or non unique
				while (
					subedges.find(curvert) != subedges.end()
					&&
					(
						(cur_end == cur_end_left  && countLeftEdges (subedges.find(curvert)->second) == 1)
						||
						(cur_end == cur_end_right && countRightEdges(subedges.find(curvert)->second) == 1)
					)
				)
				{
					std::vector<OverlapEntry> const & V = subedges.find(curvert)->second;
					assert ( V.size() );
					OverlapEntry OE = V[0];

					for ( uint64_t i = 0; i < V.size(); ++i )
						if (
							(isLeftEdge(V[i]) && (cur_end == cur_end_left))
							||
							(isRightEdge(V[i]) && (cur_end == cur_end_right))
						)
							OE = V[i];

					switch ( OE.orientation )
					{
						case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
						case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
							cur_end = cur_end_left;
							break;
						case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
						case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
							cur_end = cur_end_right;
							break;
						default:
							break;
					}
					curvert = OE.target;

					linearStrip.push_back(std::pair<uint64_t,cur_end_type>(curvert,cur_end));
				}

				#if 1
				std::cerr << "linear strip: ";
				for ( uint64_t i = 0; i < linearStrip.size(); ++i )
					std::cerr << linearStrip[i].first << ((i+1<linearStrip.size())?";":"");
				std::cerr << std::endl;
				#endif

				variantseqs.push_back(
					std::vector < std::vector < std::pair<uint64_t,cur_end_type> > >(
						1, linearStrip
					)
				);

				// if there are edges
				if (
					subedges.find(curvert) != subedges.end()
					&&
					(
						(cur_end == cur_end_left  && countLeftEdges (subedges.find(curvert)->second))
						||
						(cur_end == cur_end_right && countRightEdges(subedges.find(curvert)->second))
					)
				)
				{
					std::map< uint64_t,std::vector<uint64_t> > const subedgesshift = extractConnectedSubEdgesShift(curvert,cur_end);

					std::pair<bool,std::map<uint64_t,uint64_t> > TS =
						libmaus2::graph::TopologicalSorting::topologicalSorting<uint64_t,libmaus2::graph::IdentityTargetProjector>(
							subedgesshift,computeId(curvert,cur_end),libmaus2::graph::IdentityTargetProjector()
						);

					if ( TS.first )
					{
						std::map<uint64_t,uint64_t> TSr;

						for ( std::map<uint64_t,uint64_t>::const_iterator ita = TS.second.begin(); ita != TS.second.end(); ++ita )
							TSr[ita->first / 2] = ita->second;

						TS.second = TSr;
					}

					#if 0
					std::pair<bool,std::map<uint64_t,uint64_t> > TS =
						libmaus2::graph::TopologicalSorting::topologicalSorting<OverlapEntry,OverlapEntryTargetProjector>(
							subedges,curvert,OverlapEntryTargetProjector()
						);
					#endif

					if ( ! TS.first )
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "Input graph in findBubbles contains cycles" << std::endl;
						lme.finish();
						throw lme;
					}

					struct QueueElement
					{
						uint64_t v;
						uint64_t r;
						cur_end_type cur_end;
						uint64_t pathid;

						QueueElement() : v(0), r(0), cur_end(cur_end_left), pathid(0) {}
						QueueElement(uint64_t const rv, uint64_t const rr, cur_end_type rcur_end, uint64_t const rpathid) : v(rv), r(rr), cur_end(rcur_end), pathid(rpathid) {}

						bool operator<(QueueElement const & O) const
						{
							return r > O.r;
						}
					};

					std::priority_queue < QueueElement > Q;
					std::map < std::pair<uint64_t,cur_end_type>, uint64_t > activecnt;

					std::vector<OverlapEntry> const & RV = subedges.find(curvert)->second;

					std::vector< std::vector< std::pair<uint64_t,cur_end_type> > > bubblePaths;

					for ( uint64_t i = 0; i < RV.size(); ++i )
						if (
							(isLeftEdge(RV[i]) && cur_end == cur_end_left)
							||
							(isRightEdge(RV[i]) && cur_end == cur_end_right)
						)
						{
							OverlapEntry const & OE = RV[i];

							switch ( OE.orientation )
							{
								case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
								case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								{
									std::vector< std::pair<uint64_t,cur_end_type> > npath;
									npath.push_back(std::pair<uint64_t,cur_end_type>(curvert,cur_end));
									npath.push_back(std::pair<uint64_t,cur_end_type>(OE.target,cur_end_left));
									uint64_t const pathid = bubblePaths.size();
									bubblePaths.push_back(npath);

									assert ( TS.second.find(OE.target) != TS.second.end() );

									Q.push(
										QueueElement(
											OE.target,
											TS.second.find(OE.target)->second,
											cur_end_left,
											pathid
										)
									);


									activecnt[std::pair<uint64_t,cur_end_type>(OE.target,cur_end_left)]++;
									break;
								}
								case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
								case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
								{
									std::vector< std::pair<uint64_t,cur_end_type> > npath;
									npath.push_back(std::pair<uint64_t,cur_end_type>(curvert,cur_end));
									npath.push_back(std::pair<uint64_t,cur_end_type>(OE.target,cur_end_right));
									uint64_t const pathid = bubblePaths.size();
									bubblePaths.push_back(npath);

									assert ( TS.second.find(OE.target) != TS.second.end() );

									Q.push(
										QueueElement(
											OE.target,
											TS.second.find(OE.target)->second,
											cur_end_right,
											pathid
										)
									);
									activecnt[std::pair<uint64_t,cur_end_type>(OE.target,cur_end_right)]++;
									break;
								}
								default:
									break;
							}
						}

					std::set<uint64_t> pathkill;
					while ( Q.size() && activecnt.size() > 1 )
					{
						QueueElement QE = Q.top(); Q.pop();

						if ( ! (--activecnt[std::pair<uint64_t,cur_end_type>(QE.v,QE.cur_end)]) )
							activecnt.erase(activecnt.find(std::pair<uint64_t,cur_end_type>(QE.v,QE.cur_end)));

						if ( subedges.find(QE.v) != subedges.end() )
						{
							std::vector<OverlapEntry> const & RV = subedges.find(QE.v)->second;
							uint64_t outedges = 0;
							uint64_t pathid = QE.pathid;

							for ( uint64_t i = 0; i < RV.size(); ++i )
								if (
									(isLeftEdge(RV[i]) && QE.cur_end == cur_end_left)
									||
									(isRightEdge(RV[i]) && QE.cur_end == cur_end_right)
								)
								{
									switch ( RV[i].orientation )
									{
										case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
										case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
										{
											if ( outedges )
											{
												pathid = bubblePaths.size();
												std::vector< std::pair<uint64_t,cur_end_type> > opath = bubblePaths[QE.pathid];
												assert ( opath.size() );
												opath.pop_back();
												bubblePaths.push_back(opath);
											}

											bubblePaths[pathid].push_back(std::pair<uint64_t,cur_end_type>(RV[i].target,cur_end_left));
											Q.push(QueueElement(RV[i].target,TS.second.find(RV[i].target)->second,cur_end_left,pathid));
											activecnt[std::pair<uint64_t,cur_end_type>(RV[i].target,cur_end_left)]++;
											outedges++;
											break;
										}
										case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
										case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
										{
											if ( outedges )
											{
												pathid = bubblePaths.size();
												std::vector< std::pair<uint64_t,cur_end_type> > opath = bubblePaths[QE.pathid];
												assert ( opath.size() );
												opath.pop_back();
												bubblePaths.push_back(opath);
											}

											bubblePaths[pathid].push_back(std::pair<uint64_t,cur_end_type>(RV[i].target,cur_end_left));
											Q.push(QueueElement(RV[i].target,TS.second.find(RV[i].target)->second,cur_end_right,pathid));
											activecnt[std::pair<uint64_t,cur_end_type>(RV[i].target,cur_end_right)]++;
											outedges++;
											break;
										}
										default:
											break;
									}
								}

						}
						else
						{
							pathkill.insert(QE.pathid);

							std::cerr << "warning: input graph in findBubbles contains tips" << std::endl;

							#if 0
							libmaus2::exception::LibMausException lme;
							lme.getStream() << "Input graph in findBubbles contains tips" << std::endl;
							lme.finish();
							throw lme;
							#endif
						}
					}

					uint64_t o = 0;
					for ( uint64_t i = 0; i < bubblePaths.size(); ++i )
						if ( pathkill.find(i) == pathkill.end() )
							bubblePaths[o++] = bubblePaths[i];
					bubblePaths.resize(o);

					if ( activecnt.size() == 1 )
					{
						std::cerr << "bubble opening at vertex " << curvert << " closing at vertex " << activecnt.begin()->first.first
							<< " paths " << bubblePaths.size()
							<< std::endl;

						#if 0
						for ( uint64_t j = 0; j < bubblePaths.size(); ++j )
						{
							std::vector< std::pair<uint64_t,cur_end_type> > const & lpath = bubblePaths[j];

							std::cerr << "bubblePath[" << j << "]=";

							for ( uint64_t i = 0; i < lpath.size(); ++i )
								std::cerr << lpath[i].first << ((i+1<lpath.size())? ";":"");
							std::cerr << std::endl;
						}
						#endif

						variantseqs.push_back(bubblePaths);

						curvert = activecnt.begin()->first.first;
						cur_end = activecnt.begin()->first.second;
					}
					else
					{
						libmaus2::exception::LibMausException lme;
						lme.getStream() << "Unable to find closing node for bubble starting at node " << curvert << std::endl;
						lme.finish();
						throw lme;
					}
				}
			}

			std::cerr << "variant seqs for component id " << componentid << ": ";
			for ( uint64_t i = 0; i < variantseqs.size(); ++i )
				std::cerr << variantseqs[i].size() << ";";
			std::cerr << std::endl;

			uint64_t prod = 1;
			for ( uint64_t i = 0; i < variantseqs.size(); ++i )
				prod *= variantseqs[i].size();

			if ( RC && (prod > 1) )
			{
				std::ostringstream fnostr;

				if ( prefix )
					fnostr << *prefix << "_";

				fnostr << "bubbles_" << (bubbleid++) << ".fasta";
				std::string const fn = fnostr.str();
				libmaus2::aio::OutputStreamInstance COS(fn);


				for ( uint64_t z = 0; z < prod; ++z )
				{
					std::vector<uint64_t> idx;
					uint64_t zz = z;
					for ( uint64_t i = 0; i < variantseqs.size(); ++i )
					{
						idx.push_back(zz % variantseqs[i].size());
						zz /= variantseqs[i].size();
					}

					std::vector< std::pair<uint64_t,cur_end_type> > path;
					for ( uint64_t i = 0; i < idx.size(); ++i )
					{
						std::vector< std::pair<uint64_t,cur_end_type> > const & subpath =
							variantseqs[i][idx[i]];

						std::vector< uint64_t > subpathids;
						for ( uint64_t j = 0; j < subpath.size(); ++j )
							subpathids.push_back(subpath[j].first);

						std::pair<std::string,double> partcontigdepth = createLinearContig(subpathids,*RC);

						std::cerr << idx[i] << "{" << partcontigdepth.first.size() << "," << partcontigdepth.second << "}";

						if ( variantseqs[i].size() > 1 )
						{
							assert ( subpath.size() >= 3 );

							std::vector< uint64_t > subpathids;
							for ( uint64_t j = 1; j+1 < subpath.size(); ++j )
								subpathids.push_back(subpath[j].first);

							std::pair<std::string,double> partcontigdepth = createLinearContig(subpathids,*RC);

							std::cerr << "(" << partcontigdepth.second << ")";
						}

						std::cerr << ";";

						std::copy(
							subpath.begin() + (i==0?0:1),
							subpath.end(),
							std::back_insert_iterator< std::vector< std::pair<uint64_t,cur_end_type> > >(path)
						);
					}

					std::vector<uint64_t> vertices;
					for ( uint64_t i = 0; i < path.size(); ++i )
						vertices.push_back(path[i].first);

					std::pair<std::string,double> contigdepth = createLinearContig(vertices,*RC);
					std::string const contig = contigdepth.first;
					COS << ">contig_" << componentid << "_" << z << "_" << contig.size() << "_" << contigdepth.second << std::endl;
					COS << contig << std::endl;

					std::cerr << " reads " << path.size() << " contig length " << contig.size();
					std::cerr << std::endl;
				}

				COS.flush();
				COS.close();
			}
			#endif
		}

		#if 0
		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		{
			uint64_t const v = ita->first;

			if ( seen.find(v) == seen.end() )
			{
				// set of vertices
				std::deque<uint64_t> C;

				// edges only on one side
				if (
					(countLeftEdges(v) == 1 && countRightEdges(v) == 0)
					||
					(countLeftEdges(v) == 0 && countRightEdges(v) == 1)
				)
				{
					cur_end_type cur_end = countLeftEdges(v) ? cur_end_left : cur_end_right;

					uint64_t vv = v;
					C.push_back(vv);

					while (
						(cur_end == cur_end_left  && countLeftEdges(vv) == 1  && countRightEdges(vv) <= 1)
						||
						(cur_end == cur_end_right && countRightEdges(vv) == 1 && countLeftEdges(vv) <=1)
					)
					{
						OverlapEntry const OE = (cur_end == cur_end_left) ? getFirstLeftEdge(vv) : getFirstRightEdge(vv);

						switch ( OE.orientation )
						{
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
								cur_end = cur_end_left;
								break;
							case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
							case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
								cur_end = cur_end_right;
								break;
							default:
								break;
						}

						vv = OE.target;
						C.push_back(vv);
					}

					for ( uint64_t i = 0; i < C.size(); ++i )
						seen.insert(C[i]);

					if (
						(cur_end == cur_end_left && countLeftEdges(vv) == 0)
						||
						(cur_end == cur_end_right && countRightEdges(vv) == 0)
					)
					{
						R.push_back(std::vector<uint64_t>(C.begin(),C.end()));
					}
				}
			}
		}

		return R;
		#endif
	}


	std::vector < std::vector<uint64_t> > findTips(uint64_t const maxfollow = std::numeric_limits<uint64_t>::max()) const
	{
		std::vector < std::vector<uint64_t> > R;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
		{
			uint64_t const preseed = ita->first;
			// set of vertices
			std::deque<uint64_t> C;
			FollowUniqueResult FUR;

			if ( extractTip(preseed,C,FUR,maxfollow) )
			{
				// std::cerr << "tip starting at " << preseed << " length " << C.size() << std::endl;
				R.push_back(std::vector<uint64_t>(C.begin(),C.end()));
			}
		}

		return R;
	}

	void insertEdges(std::map<uint64_t,std::vector<OverlapEntry> > const & removedEdges)
	{
		for ( std::map<uint64_t,std::vector<OverlapEntry> >::const_iterator ita = removedEdges.begin(); ita != removedEdges.end(); ++ita )
		{
			std::vector<OverlapEntry> const & V = ita->second;

			for ( uint64_t i = 0; i < V.size(); ++i )
				if ( ! hasEdge(ita->first,V[i].target) )
					preedges[ita->first].push_back(V[i]);
		}
	}

	template<typename iterator>
	void removeCrossEdges(iterator ita, iterator ite, std::map<uint64_t,std::vector<OverlapEntry> > & removedEdges)
	{
		std::set<uint64_t> S(ita,ite);
		std::set<uint64_t> Ekill;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
			// node is not in S
			if ( S.find(ita->first) == S.end() )
			{
				// keep edges not pointing into S
				std::vector<OverlapEntry> & V = ita->second;
				uint64_t o = 0;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( S.find(V[i].target) == S.end() )
						V[o++] = V[i];
					else
						removedEdges[ita->first].push_back(V[i]);

				V.resize(o);

				if ( ! o )
					Ekill.insert(ita->first);
			}
			// node is in S
			else
			{
				// keep edges pointing into S
				std::vector<OverlapEntry> & V = ita->second;
				uint64_t o = 0;
				for ( uint64_t i = 0; i < V.size(); ++i )
					if ( S.find(V[i].target) != S.end() )
						V[o++] = V[i];
					else
						removedEdges[ita->first].push_back(V[i]);

				V.resize(o);

				if ( ! o )
					Ekill.insert(ita->first);
			}

		for ( std::set<uint64_t>::const_iterator ita = Ekill.begin(); ita != Ekill.end(); ++ita )
			preedges.erase(preedges.find(*ita));
	}

	void removeEdges(uint64_t const s, ReadContainer const * /* RC */ = 0)
	{
		#if 0
		if ( RC )
			std::cerr << "removing edges for read " << RC->getRead(s) << std::endl;
		#endif

		if ( preedges.find(s) != preedges.end() )
			preedges.erase(preedges.find(s));

		std::set<uint64_t> vecremove;

		for ( std::map< uint64_t,std::vector<OverlapEntry> >::iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
		{
			std::vector<OverlapEntry> & V = ita->second;
			uint64_t o = 0;

			for ( uint64_t i = 0; i < V.size(); ++i )
				if ( V[i].target != s )
					V[o++] = V[i];

			V.resize(o);

			if ( ! o )
				vecremove.insert(ita->first);
		}

		for ( std::set<uint64_t>::const_iterator ita = vecremove.begin(); ita != vecremove.end(); ++ita )
			preedges.erase(preedges.find(*ita));
	}

	void removeEdges(std::vector<uint64_t> const & V, ReadContainer const * RC = 0)
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			removeEdges(V[i],RC);
	}

	void removeEdges(std::vector< std::vector<uint64_t> > const & V, ReadContainer const * RC = 0)
	{
		for ( uint64_t i = 0; i < V.size(); ++i )
			removeEdges(V[i],RC);
	}

	FollowUniqueResult followUniqueLeft(
		uint64_t const i,
		uint64_t maxfollow = std::numeric_limits<uint64_t>::max()
	) const
	{
		NullContainer NC;
		return followUnique(i,cur_end_left,NC,maxfollow);
	}

	FollowUniqueResult followUniqueLeftRecord(
		uint64_t const i,
		std::deque<uint64_t> & C,
		uint64_t maxfollow = std::numeric_limits<uint64_t>::max()
	) const
	{
		DequeFrontContainer D(C);
		return followUnique(i,cur_end_left,D,maxfollow);
	}

	FollowUniqueResult followUniqueRight(
		uint64_t const i,
		uint64_t maxfollow = std::numeric_limits<uint64_t>::max()
	) const
	{
		NullContainer NC;
		return followUnique(i,cur_end_right,NC,maxfollow);
	}

	FollowUniqueResult followUniqueRightRecord(
		uint64_t const i,
		std::deque<uint64_t> & C,
		uint64_t maxfollow = std::numeric_limits<uint64_t>::max()
	) const
	{
		DequeBackContainer D(C);
		return followUnique(i,cur_end_right,D,maxfollow);
	}

	/*
	 * follow unique chain of edges until we reach a split node (not included) or the end of the chain
	 *
	 * returns a tuple consisting of the last node on the chain, the chain length in edges and
	 * the end on the last node which cannot be (uniquely) extended
	 */
	template<typename container_type>
	FollowUniqueResult followUnique(
		uint64_t const i,
		cur_end_type cur_end,
		container_type & container,
		uint64_t maxfollow = std::numeric_limits<uint64_t>::max()
	) const
	{
		uint64_t v = i;
		uint64_t c = 0;

		bool ok = true;

		while ( ok && c < maxfollow )
		{
			ok = false;

			switch ( cur_end )
			{
				case cur_end_left:
					if ( countLeftEdges(v) == 1 )
					{
						OverlapEntry OE = getFirstLeftEdge(v);

						// edge ends on the right side of the other read
						if ( OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back )
						{
							// check that following read has only one outgoing edge in direction and it leads us back to v
							if (
								countRightEdges(OE.target) == 1
								&&
								getFirstRightEdge(OE.target).target == v
							)
							{
								v = OE.target;
								cur_end = cur_end_left;
								ok = true;
								c += 1;

								container(v);
							}
						}
						// edge ends on the left side of the other read
						else
						{
							// check that following read has only one outgoing edge in direction and it leads us back to v
							if (
								countLeftEdges(OE.target) == 1
								&&
								getFirstLeftEdge(OE.target).target == v
							)
							{
								v = OE.target;
								cur_end = cur_end_right;
								ok = true;
								c += 1;

								container(v);
							}
						}
					}
					break;
				case cur_end_right:
					if ( countRightEdges(v) == 1 )
					{
						OverlapEntry OE = getFirstRightEdge(v);

						// edge ends on the right side of the other read
						if ( OE.orientation == libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back )
						{
							// check that following read has only one outgoing edge in direction and it leads us back to v
							if (
								countRightEdges(OE.target) == 1
								&&
								getFirstRightEdge(OE.target).target == v
							)
							{
								v = OE.target;
								cur_end = cur_end_left;
								ok = true;
								c += 1;

								container(v);
							}
						}
						// edge ends on the left side of the other read
						else
						{
							// check that following read has only one outgoing edge in direction and it leads us back to v
							if (
								countLeftEdges(OE.target) == 1
								&&
								getFirstLeftEdge(OE.target).target == v
							)
							{
								v = OE.target;
								cur_end = cur_end_right;
								ok = true;
								c += 1;

								container(v);
							}
						}
					}
					break;
			}
		}

		return FollowUniqueResult(v,c,cur_end);
	}

	uint64_t countLeftEdges(uint64_t const i) const
	{
		if ( preedges.find(i) == preedges.end() )
			return 0;

		std::vector<OverlapEntry> const & E = preedges.find(i)->second;

		uint64_t c = 0;
		for ( uint64_t i = 0; i < E.size(); ++i )
			switch ( E[i].orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
					c++;
					break;
				default:
					break;
			}

		return c;
	}

	void fail(std::string const & message) const
	{
		libmaus2::exception::LibMausException lme;
		lme.getStream() << message << "\n";
		lme.finish();
		throw lme;
	}

	OverlapEntry const & getFirstLeftEdge(uint64_t const i) const
	{
		if ( preedges.find(i) == preedges.end() )
			fail("getFirstLeftEdge() called on vertex with no left edges");

		std::vector<OverlapEntry> const & E = preedges.find(i)->second;

		for ( uint64_t i = 0; i < E.size(); ++i )
			switch ( E[i].orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
					return E[i];
					break;
				default:
					break;
			}

		fail("getFirstLeftEdge() called on vertex with no left edges");
		std::terminate();
	}

	OverlapEntry const & getFirstRightEdge(uint64_t const i) const
	{
		if ( preedges.find(i) == preedges.end() )
			fail("getFirstRightEdge() called on vertex with no right edges");

		std::vector<OverlapEntry> const & E = preedges.find(i)->second;

		for ( uint64_t i = 0; i < E.size(); ++i )
			switch ( E[i].orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
					return E[i];
					break;
				default:
					break;
			}

		fail("getFirstRightEdge() called on vertex with no right edges");
		std::terminate();
	}

	uint64_t countRightEdges(uint64_t const i) const
	{
		if ( preedges.find(i) == preedges.end() )
			return 0;

		std::vector<OverlapEntry> const & E = preedges.find(i)->second;

		uint64_t c = 0;
		for ( uint64_t i = 0; i < E.size(); ++i )
			switch ( E[i].orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
					c++;
					break;
				default:
					break;
			}

		return c;
	}

	void printSubGraph(std::ostream & out, std::set<uint64_t> const & vertices)
	{
		// get set of undirected edges
		std::set< std::pair<uint64_t,uint64_t> > P;
		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		{
			uint64_t s = ita->first;

			if ( vertices.find(s) != vertices.end() )
			{
				for ( uint64_t i = 0; i < ita->second.size(); ++i )
				{
					uint64_t t = ita->second[i].target;

					if ( vertices.find(t) != vertices.end() )
					{
						if ( s < t )
							P.insert(std::pair<uint64_t,uint64_t>(s,t));
						else
							P.insert(std::pair<uint64_t,uint64_t>(t,s));
					}
				}
			}
		}

		out << "graph {\n";
		for ( std::set< std::pair<uint64_t,uint64_t> >::const_iterator ita = P.begin(); ita != P.end(); ++ita )
		{
			OverlapEntry const edge = getEdge(ita->first,ita->second);
			// head: at edge target vertex
			std::string head = "none";
			// tail: at edge source vertex
			std::string tail = "none";

			switch ( edge.orientation )
			{
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
					head = "normal"; // right
					tail = "inv";    // left
					break;
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
					head = "inv";
					tail = "normal";
					break;
				case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
					head = "normal";
					tail = "normal";
					break;
				case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
					head = "inv";
					tail = "inv";
					break;
				default:
					break;
			}

			std::string edgestyle = std::string("[dir=both arrowhead=")+head+" arrowtail="+tail+" ]";

			out << "\t{ \"s"
				<< ita->first
				<< "\" -- \"s"
				<< ita->second << "\" "<<edgestyle<<" ;}\n";
		}
		out << "}\n";
	}
};

#include <libmaus2/lz/PlainOrGzipStream.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>

int alternatives(libmaus2::util::ArgInfo const arginfo)
{
	std::string const hwtname = arginfo.getUnparsedRestArg(0);
	std::string const isaname = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt") + ".isa";
	std::string const primedgesname = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt") + ".primedges";
	std::string const cleanededgesname = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt") + ".cleanededges";
	std::string const straightcontigsfa = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt") + "_straight_contigs.fa";
	std::string const straightcontigsfq = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt") + "_straight_contigs.fq";
	std::string const prefix = libmaus2::util::OutputFileNameTools::clipOff(hwtname,".hwt");
	uint64_t const k = arginfo.getValueUnsignedNumeric("k",12);
	uint64_t const minreqscore = arginfo.getValueUnsignedNumeric("minreqscore",40);
	uint64_t const numthreads = arginfo.getValueUnsignedNumeric("threads",1);

	libmaus2::fm::BidirectionalDnaIndexImpCompactHuffmanWaveletTree index(hwtname,numthreads);
	typedef libmaus2::fm::BidirectionalDnaIndexImpCompactHuffmanWaveletTree::lf_type lf_type;
	libmaus2::fm::SampledISA<lf_type>::unique_ptr_type const pISA(libmaus2::fm::SampledISA<lf_type>::load(index.LF.get(),isaname));
	libmaus2::fm::SampledISA<lf_type> const & ISA = *pISA;
	std::vector<uint64_t> const seqstart = index.getSeqStartPositions();
	assert ( seqstart.size() % 2 == 0 );

	std::vector<std::string> readnames;

	if ( 1 < arginfo.restargs.size() )
	{
		std::string const readsname = arginfo.getUnparsedRestArg(1);
		libmaus2::aio::PosixFdInputStream readsCIS(readsname);
		libmaus2::lz::PlainOrGzipStream readsPOGS(readsCIS);
		libmaus2::fastx::StreamFastAReaderWrapper SFARW(readsPOGS);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

		while ( SFARW.getNextPatternUnlocked(pattern) )
		{
			readnames.push_back(pattern.sid);
		}
	}

	ReadContainer const RC(index,ISA);

	libmaus2::lcs::HammingOverlapDetection const HOD;

	std::map< uint64_t,std::vector<OverlapEntry> > preedges;
	libmaus2::parallel::PosixSpinLock preedgeslock;
	std::set<uint64_t> dup;
	libmaus2::parallel::PosixSpinLock duplock;
	libmaus2::parallel::PosixSpinLock cerrlock;

	if ( libmaus2::util::GetFileSize::fileExists(cleanededgesname) )
	{
		loadEdges(preedges, dup, cleanededgesname);
	}
	else
	{
		if ( libmaus2::util::GetFileSize::fileExists(primedgesname) )
		{
			loadEdges(preedges, dup, primedgesname);
		}
		else
		{
			uint64_t const numreads = seqstart.size()/2;

			#if defined(_OPENMP)
			#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
			#endif
			for ( uint64_t readid = 0; readid < numreads; readid++ )
			{
				uint64_t const start  = seqstart[2*readid];
				uint64_t const len    = getSeqLen(seqstart,2*readid);
				uint64_t const r0     = ISA[start];
				std::string const ref = index.getText(r0,len);
				std::string const uref = index.getTextUnmapped(r0,len);

				std::set<uint64_t> overlapcand;
				for ( uint64_t i = 0; i < ref.size()-k+1; ++i )
				{
					std::string const sub = ref.substr(i,k);

					libmaus2::fm::BidirectionalIndexInterval const bint = index.biSearchBackward(sub.begin(), k);

					for ( uint64_t i = 0; i < bint.siz; ++i )
					{
						uint64_t const rr = bint.spf + i;
						uint64_t const p = (*index.SA)[rr];
						uint64_t const seq = findSequence(seqstart,p) >> 1;

						if ( seq != readid )
							overlapcand.insert(seq);
					}
				}

				std::vector<OverlapEntry> ledges;
				for ( std::set<uint64_t>::const_iterator ita = overlapcand.begin(); ita != overlapcand.end(); ++ita )
				{
					uint64_t const mseq   = *ita;
					uint64_t const mstart = seqstart[2*mseq];
					uint64_t const mlen   = getSeqLen(seqstart,2*mseq);
					uint64_t const mr     = ISA[mstart];
					std::string const mtext  = index.getText(mr,mlen);
					std::string const umtext = index.getTextUnmapped(mr,mlen);

					libmaus2::lcs::OverlapOrientation::overlap_orientation ori;
					uint64_t overhang = 0;
					int64_t maxscore = 0;

					bool const ok = HOD.detect(
						uref,umtext,5,ori,overhang,maxscore,false /* verbose */
					);

					if ( ok && (maxscore >= static_cast<int64_t>(minreqscore)) )
					{
						if ( overhang == 0 )
						{
							#if 0
							std::cerr << ori << std::endl;
							HOD.printOverlap(std::cerr,uref,umtext,ori,overhang);
							#endif

							if ( readid > mseq )
							{
								libmaus2::parallel::ScopePosixSpinLock llock(duplock);
								dup.insert(readid);
							}
						}
						else
						{
							#if 0
							std::cerr << ori << std::endl;
							HOD.printOverlap(std::cerr,uref,umtext,ori,overhang);
							#endif

							uint64_t const overlap = mlen - overhang;

							OverlapEntry OE(ori,overhang,maxscore,overlap,mseq);
							ledges.push_back(OE);
							// std::cerr << OE << std::endl;
						}
					}

					#if 0
					if ( maxscore >= 50 )
						HOD.printOverlap(std::cerr,uref,umtext,ori,overhang);
					#endif
				}

				std::sort(ledges.begin(),ledges.end(),OverlapEntryOverlapComparator());

				#if 0
				if ( ! ledges.size() )
				{
					libmaus2::parallel::ScopePosixSpinLock lcerrlock(cerrlock);
					std::cerr << "no edges for " << uref << std::endl;
				}
				#endif

				{
				libmaus2::parallel::ScopePosixSpinLock lpreedgeslock(preedgeslock);
				preedges[readid] = ledges;
				}

				if ( (readid) % 4096 == 0 )
				{
				libmaus2::parallel::ScopePosixSpinLock lcerrlock(cerrlock);
				std::cerr << readid << "/" << (seqstart.size()/2) << "\t" << len << "\t" << overlapcand.size() << "\t" << dup.size() << std::endl;
				}
			}

			saveEdges(preedges, dup, primedgesname);
		}

		// erase outgoing edges from duplicates
		for ( std::set<uint64_t>::const_iterator ita = dup.begin(); ita != dup.end(); ++ita )
			if ( preedges.find(*ita) != preedges.end() )
				preedges.erase(preedges.find(*ita));

		// erase edges to duplicates
		for ( std::map< uint64_t,std::vector<OverlapEntry> >::iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
		{
			std::vector<OverlapEntry> & V = ita->second;

			uint64_t o = 0;
			for ( uint64_t i = 0; i < V.size(); ++i )
				if ( dup.find(V[i].target) == dup.end() )
					V[o++] = V[i];

			V.resize(o);
		}

		std::vector<uint64_t> keys;
		for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin();
			ita != preedges.end(); ++ita )
			keys.push_back(ita->first);

		std::cerr << "Computing edge kill list...";
		std::map < uint64_t, std::set<uint64_t> > edgekill;
		libmaus2::parallel::PosixSpinLock edgekilllock;
		#if defined(_OPENMP)
		#pragma omp parallel for num_threads(numthreads)
		#endif
		for ( uint64_t ia = 0; ia < keys.size(); ++ia )
		{
			#if 0
			if ( (ia+1) % 1024 == 0 )
				std::cerr << (ia+1) << "/" << keys.size() << std::endl;
			#endif

			uint64_t const a = keys[ia];
			assert ( preedges.find(a) != preedges.end() );
			std::vector<OverlapEntry> const & aedges = preedges.find(a)->second;
			std::set<uint64_t> ledgekill; // = edgekill[a];

			for ( uint64_t i = 0; i < aedges.size(); ++i )
			{
				OverlapEntry const & abedge = aedges[i];
				uint64_t const b = abedge.target;
				libmaus2::lcs::OverlapOrientation::overlap_orientation aborientation = abedge.orientation;

				for ( uint64_t j = i+1; j < aedges.size(); ++j )
					if ( ledgekill.find(j) == ledgekill.end() )
					{
						OverlapEntry const & acedge = aedges[j];
						uint64_t const c = acedge.target;
						uint64_t const acoverlap = acedge.overlap;
						libmaus2::lcs::OverlapOrientation::overlap_orientation acorientation = aedges[j].orientation;

						if ( preedges.find(b) != preedges.end() )
						{
							std::vector<OverlapEntry> const & bedges = preedges.find(b)->second;

							for ( uint64_t z = 0; z < bedges.size(); ++z )
							{
								OverlapEntry const & bcedge = bedges[z];

								if (
									bcedge.target == c
									&&
									acoverlap < bcedge.overlap
									&&
									isABCOverlap(aborientation,acorientation,bedges[z].orientation)
								)
								{
									#if 0
									std::cerr << std::string(80,'-') << std::endl;
									// A -> B
									std::cerr << a << " -> " << aedges[i] << std::endl;
									// A -> C
									std::cerr << a << " -> " << aedges[j] << std::endl;
									// B -> C
									std::cerr << b << " -> " << bedges[z] << std::endl;
									#endif

									// kill edge A -> C
									ledgekill.insert(j);
								}
							}
						}
					}
			}

			{
			libmaus2::parallel::ScopePosixSpinLock ledgekilllock(edgekilllock);
			edgekill[a] = ledgekill;
			}

			#if 0
			for ( uint64_t i = 0; i < aedges.size(); ++i )
			{
				std::cerr << "\t[" << i << "]=" << aedges[i] << " " << ((ledgekill.find(i)!=ledgekill.end())?"K":"") << std::endl;

				#if 0
				uint64_t const c = aedges[i].target;

				for ( uint64_t j = 0; j < i; ++j )
				{
					uint64_t const b = aedges[j].target;
					if ( preedges.find(b) != preedges.end() )
					{
						std::vector<OverlapEntry> const & zedges = preedges.find(b)->second;
						for ( uint64_t z = 0; z < zedges.size(); ++z )
							if ( zedges[z].target == c )
							{
								std::cerr << "\t\t" << b << " -> " << zedges[z] << std::endl;

								std::string const read_a = RC.getRead(a);
								std::string const read_b = RC.getRead(b);
								std::string const read_c = RC.getRead(c);

								{
								libmaus2::lcs::OverlapOrientation::overlap_orientation ori;
								uint64_t overhang = 0;
								int64_t maxscore = 0;
								std::cerr << "A -> B" << std::endl;
								HOD.detect(read_a,read_b,5,ori,overhang,maxscore,true);
								std::cerr << "A -> C" << std::endl;
								HOD.detect(read_a,read_c,5,ori,overhang,maxscore,true);
								std::cerr << "B -> C" << std::endl;
								HOD.detect(read_b,read_c,5,ori,overhang,maxscore,true);
								}

							}
					}
				}
				#endif
			}
			#endif
		}
		std::cerr << "done.\n";

		std::cerr << "Performing transitive reduction...";
		libmaus2::util::Histogram hist;
		for ( uint64_t ia = 0; ia < keys.size(); ++ia )
		{
			uint64_t const a = keys[ia];

			if ( edgekill.find(a) != edgekill.end() )
			{
				std::set<uint64_t> const & ks = edgekill.find(a)->second;
				uint64_t o = 0;
				std::vector<OverlapEntry> & aedges = preedges.find(a)->second;

				for ( uint64_t i = 0; i < aedges.size(); ++i )
					if ( ks.find(i) == ks.end() )
						aedges[o++] = aedges[i];

				aedges.resize(o);

				if ( aedges.size() == 0 )
					preedges.erase(preedges.find(a));
				else
				{
					hist(aedges.size());
				}
			}
		}
		std::cerr << "done." << std::endl;
		hist.print(std::cerr);
		std::cerr << "number of duplicate reads " << dup.size() << std::endl;

		saveEdges(preedges, dup, cleanededgesname);
	}

	EdgeSet ES(preedges);

	std::cerr << "removing short tips." << std::endl;
	// look for graph tips (short branches)
	std::vector < std::vector<uint64_t> > tips = ES.findTips(5/* max follow */);
	// remove tips
	ES.removeEdges(tips,&RC);

	std::cerr << "removing short contigs." << std::endl;
	// look for short contigs
	std::vector < std::vector<uint64_t> > shorts = ES.findShortContigs(10);
	// remove short contigs
	ES.removeEdges(shorts,&RC);

	std::cerr << "finding straight contigs." << std::endl;
	std::vector < std::vector<uint64_t> > straight = ES.findStraightContigs();

	std::cerr << "creating linear contigs as fasta...";
	{
	libmaus2::aio::OutputStreamInstance COS(straightcontigsfa);
	for ( uint64_t i = 0; i < straight.size(); ++i )
	{
		std::string const s = ES.createLinearContig(straight[i],RC).first;
		COS << ">contig_" << i << "_" << s.size() << "\n";
		COS << s << "\n";
	}
	COS.flush();
	}
	std::cerr << "done." << std::endl;

	{
	libmaus2::aio::OutputStreamInstance COS(straightcontigsfq);
	for ( uint64_t i = 0; i < straight.size(); ++i )
	{
		std::string const s = ES.createLinearContig(straight[i],RC).first;
		COS << "@contig_" << i << "_" << s.size() << "\n";
		COS << s << "\n";
		COS << "+\n";
		COS << std::string(s.size(),'H') << "\n";
	}
	COS.flush();
	}

	std::cerr << "Removing straight contig edges...";
	ES.removeEdges(straight,&RC);
	std::cerr << "done." << std::endl;

	#if 1
	std::cerr << "Breaking strongly connected subgraphs...";
	std::map<uint64_t,std::vector<OverlapEntry> > removedStrongConnectedComponentEdges;
	ES.checkSymmetry();
	ES.breakStronglyConnectedComponents(removedStrongConnectedComponentEdges,&RC);
	std::cerr << "done." << std::endl;
	#endif

	#if 1
	std::cerr << "Disconnecting long tips...";
	std::map<uint64_t,std::vector<OverlapEntry> > removedLongTipEdges;
	ES.disconnectTips(removedLongTipEdges);
	std::cerr << "done." << std::endl;
	#endif

	#if 1
	std::cerr << "Checking for bubbles...";
	std::map<uint64_t,std::vector<OverlapEntry> > removedBubbleEdges;
	ES.findBubbles(removedBubbleEdges,&RC,&prefix);
	std::cerr << "done." << std::endl;
	#endif

	std::cerr << "Collecting unused vertex ids...";
	std::set<uint64_t> unusedvertices;
	for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
		unusedvertices.insert(ita->first);
	std::cerr << "done." << std::endl;

	std::map<uint64_t,uint64_t> cid;

	uint64_t curcid = 0;
	std::map<uint64_t,std::vector<OverlapEntry> > removedEdges;
	while ( unusedvertices.size() )
	{
		uint64_t const prim = *(unusedvertices.begin());
		EdgeSet::cur_end_type root1end = EdgeSet::cur_end_left;
		uint64_t root1;

		if ( ES.hasBorderVertexWithLeftEdge(ES.preedges,*(unusedvertices.begin())) )
		{
			root1 = ES.getBorderVertexWithLeftEdge(ES.preedges,*(unusedvertices.begin()));
			root1end = EdgeSet::cur_end_left;
		}
		else if ( ES.hasBorderVertexWithRightEdge(ES.preedges,*(unusedvertices.begin())) )
		{
			root1 = ES.getBorderVertexWithRightEdge(ES.preedges,*(unusedvertices.begin()));
			root1end = EdgeSet::cur_end_right;
		}
		else
		{
			root1 = (*(unusedvertices.begin()));
			root1end = EdgeSet::cur_end_left;
		}

		#if 0
		uint64_t const root1 =
			ES.hasBorderVertexWithLeftEdge(ES.preedges,*(unusedvertices.begin())) ?
				ES.getBorderVertexWithLeftEdge(ES.preedges,*(unusedvertices.begin()))
				: (*(unusedvertices.begin()));
		#endif

		#if 0
		std::pair< std::vector< uint64_t >, std::vector< uint64_t > > SCC =
			libmaus2::graph::StronglyConnectedComponents::strongConnectContract<OverlapEntry,OverlapEntryTargetProjector>(
				ES.extractConnectedSubEdges(prim),root1
				// ES.preedges,prim
			);
		#endif

		// compute strongly connected components based on 1d edges
		std::pair< std::vector< uint64_t >, std::vector< uint64_t > > const SCC = ES.getStronglyConnectedComponents(root1,root1end);


		bool hasloop = false;
		std::deque<uint64_t> contig;
		for ( uint64_t i = 1 ; i < SCC.second.size(); ++i )
			if ( SCC.second[i]-SCC.second[i-1] > 1 )
			{
				hasloop = true;

				for ( uint64_t j = SCC.second[i-1]; j < SCC.second[i] ; ++j )
					contig.push_back(SCC.first[j]);

				ES.removeCrossEdges(
					SCC.first.begin() + SCC.second[i-1],
					SCC.first.begin() + SCC.second[i],
					removedEdges
				);
			}


		if ( hasloop )
		{
			std::cerr << "extracting strongly connected components of size " << contig.size() << " from prim=" << prim << std::endl;
		}
		else
		{
			contig.push_back(prim);
			ES.followUniqueLeftRecord(prim,contig);
			ES.followUniqueRightRecord(prim,contig);
		}

		for ( uint64_t i = 0; i < contig.size(); ++i )
		{
			if ( unusedvertices.find(contig[i]) != unusedvertices.end() )
				unusedvertices.erase(unusedvertices.find(contig[i]));

			cid [ contig[i] ] = curcid;
		}

		#if 1
		std::cerr << "got contig id " << curcid <<" of length size " << contig.size() << " from prim " << prim;

		std::cerr << " left " << contig.front();
		if ( preedges.find(contig.front()) != preedges.end() )
			std::cerr << " edges " << preedges.find(contig.front())->second.size();
		std::cerr << " right " << contig.back();
		if ( preedges.find(contig.back()) != preedges.end() )
			std::cerr << " edges " << preedges.find(contig.back())->second.size();
		std::cerr << std::endl;
		#endif

		curcid++;
	}

	// reinsert removed edges
	ES.insertEdges(removedEdges);

	// get set of undirected edges
	std::set< std::pair<uint64_t,uint64_t> > P;
	for ( std::map< uint64_t,std::vector<OverlapEntry> >::const_iterator ita = preedges.begin(); ita != preedges.end(); ++ita )
	{
		uint64_t s = ita->first;

		for ( uint64_t i = 0; i < ita->second.size(); ++i )
		{
			uint64_t t = ita->second[i].target;
			if ( s < t )
				P.insert(std::pair<uint64_t,uint64_t>(s,t));
			else
				P.insert(std::pair<uint64_t,uint64_t>(t,s));
		}
	}

	std::cout << "graph {\n";
	std::cout << "\tsplines=\"line\";\n";
	std::cout << "\toverlap=scale;\n";
	for ( std::set< std::pair<uint64_t,uint64_t> >::const_iterator ita = P.begin(); ita != P.end(); ++ita )
	{
		OverlapEntry const edge = ES.getEdge(ita->first,ita->second);
		// head: at edge target vertex
		std::string head = "none";
		// tail: at edge source vertex
		std::string tail = "none";

		switch ( edge.orientation )
		{
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_front:
				head = "normal"; // right
				tail = "inv";    // left
				break;
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_back:
				head = "inv";
				tail = "normal";
				break;
			case libmaus2::lcs::OverlapOrientation::overlap_a_front_dovetail_b_front:
				head = "normal";
				tail = "normal";
				break;
			case libmaus2::lcs::OverlapOrientation::overlap_a_back_dovetail_b_back:
				head = "inv";
				tail = "inv";
				break;
			default:
				break;
		}

		std::string edgestyle = std::string("[dir=both arrowhead=")+head+" arrowtail="+tail+" ]";

		if ( ita->first < readnames.size() && ita->second < readnames.size() )
		{
			std::cout << "\t{ \"s"
				<< ita->first << "_" << cid[ita->first] << "_" << readnames.at(ita->first)
				<< "\" -- \"s"
				<< ita->second << "_" << cid[ita->second] << "_" << readnames.at(ita->second) << "\" "<<edgestyle<<" ;}\n";
		}
		else
		{
			std::cout << "\t{ \"s"
				<< ita->first << "_" << cid[ita->first]
				<< "\" -- \"s"
				<< ita->second << "_" << cid[ita->second] << "\" "<<edgestyle<<" ;}\n";
		}
	}
	std::cout << "}\n";

	// ES.printEdges(std::cerr);

	return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		
		std::cerr << "[V] this is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
		std::cerr << "[V] " << PACKAGE_NAME << " is distributed under version 3 of the GNU General Public License" << std::endl;
		
		if ( ! arginfo.restargs.size() )
		{
			std::cerr << "usage: " << argv[0] << " reads.hwt" << std::endl;
			std::cerr << "options (provide using key=value BEFORE reads.hwt):" << std::endl;
			std::cerr << "\tk: seed size for read overlap detection (default 12)" << std::endl;
			std::cerr << "\tminreqscore: minimum suffix/prefix alignment score for two reads to be considered overlapping (default 40)" << std::endl;
			std::cerr << "\tthreads: number of threads to use for overlap computation (default 1)" << std::endl;
			
			return EXIT_FAILURE;
		}
		
		return alternatives(arginfo);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
