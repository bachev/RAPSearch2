#ifndef __HASHSEARCH_H__
#define __HASHSEARCH_H_

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <boost/thread/thread.hpp>

#include "BlastStat.hpp"
#include "paras.h"
#include "cindex.hpp"
#include "hitUnit.hpp"
#include "Seg.hpp"


typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

typedef std::vector<char> POOL;
typedef POOL::iterator ITER;

typedef std::vector<uint> VUINT;
typedef std::vector<uchar> VUCHAR;
typedef std::vector<ushort> VUSHORT;
typedef std::vector<VUINT> MINDEX;
typedef std::vector<std::string> VNAMES;
typedef std::vector<VUSHORT> VCOMP;


// query package including sequences, length, and names
class CQrPckg
{
public:
	CQrPckg(VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames)
		: m_vSeqs(vSeqs), m_vLens(vLens), m_vNames(vNames) {}

public:
	VUCHAR& m_vSeqs;
	VUINT& m_vLens;
	VNAMES& m_vNames;
};


// query package including sequences, length, names, hash table, and statistical information
class CDbPckg
{
public:
	CDbPckg(MINDEX& vHash, VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames, VCOMP& vComp, std::vector<double>& vFreq, VUINT& vWordCnts, uint& unMedian)
		: m_vHash(vHash), m_vSeqs(vSeqs), m_vLens(vLens), m_vNames(vNames), m_vComp(vComp), m_vFreq(vFreq), m_vWordCnts(vWordCnts), m_unMedian(unMedian) {}

public:
	MINDEX& m_vHash;
	VUCHAR& m_vSeqs;
	VUINT& m_vLens;
	VNAMES& m_vNames;
	VCOMP& m_vComp;
	std::vector<double>& m_vFreq;
	VUINT& m_vWordCnts;
	uint m_unMedian;
};


// alignment package including sequence, length and the starting position of the seed
class CAlnPckg
{
public:
	CAlnPckg(uchar* p, uint unLen, uint unSeedBeg)
		: m_pSeq(p), m_unLen(unLen), m_unSeedBeg(unSeedBeg) {}

public:
	uchar* m_pSeq;
	uint m_unLen;
	uint m_unSeedBeg;
};


// structure including the information of a hit
typedef struct STALNMNT
{
	int nScore;
	int nMatch;
	int nQFwd;
	int nDFwd;
	int nQBwd;
	int nDBwd;
	std::vector<char> vMode;
	std::vector<int> vLen;
}STAlnmnt;


/* map for store hit results */
typedef std::pair<int, int> PIDX;
struct pairComp
{
	bool operator() (const PIDX& lhs, const PIDX& rhs) const
	{
		if (lhs.first != rhs.first)
		{
			return lhs.first < rhs.first;
		}
		else
		{
			return lhs.second < rhs.second;
		}
	}
};
typedef std::multimap<PIDX, CHitUnit, pairComp> MRESULT;
typedef MRESULT::iterator MIT;

/* the class for indexing, searching */
class CHashSearch
{
public:
	CHashSearch(int nThreadNum);
	~CHashSearch(void)
	{
		if (NULL != m_pBlastSig)
		{
			delete m_pBlastSig;
		}

		if (NULL != m_pComptor)
		{
			delete m_pComptor;
		}
	}

	// do the protein database search
	void Process(char* szDBFile, char* szQFile, char* szOFile, int nStdout, bool bEvalue, bool bLogE, double dThr, int nMaxOut, int nMaxM8, int nQueryTypeq, bool bPrintEmpty, bool bGapExt, bool bAcc, bool bHssp, int nMinLen, bool bXml, uint unDSize = 300000000, uint unQSize = 500000000, uint unMer = 6);
	// indexing the database
	void Process(char* szDBFile, char* szDbHash, bool bFullId, int nSplitNum = 0, uint unMer = 6);

private:
	// read file and build the hash table
	int BuildDHash(const char* szDbFile, std::string& sOutFile, int nSplitNum, bool bFullId); 
	// read file and build the hash table
	int BuildQHash(std::istream& input, int nQueryType, std::map<std::string,char>& mTransTable, std::map<char,char>& mComple, Seg* seg, Seg* segsht, std::vector<uchar>& vQSeqs, std::vector<uint>& vQLens, VNAMES& vQNames); 
	// probe that what is the type of query
	int GuessQueryType(POOL& vPool);
	// char -> compressed code
	void Encode(VUCHAR& v); 
	// char -> compressed code & count DB
	long int Encode(VUCHAR& v, std::vector<double>& vFreq); 
	// compressed code -> char
	void Decode(const std::vector<uchar>& v, std::string& sOut); 
	// char -> 10-base index
	int Tran2Ten(const std::vector<uchar>& v, uint nBeg); 
	// char -> 10-base index
	int Tran2Ten(CAlnPckg& QrAln, std::vector<char>& vValid); 

	// search all sequences in database
	void Search(std::string& sDbPre, int nSeqNum, std::vector<uchar>& vQSeqs, std::vector<uint>& vQLens, VNAMES& vQNames);
	// search for each search in database
	void Searching(int nQrIdx, CQrPckg& Query, CDbPckg& Db); 
	void ExtendSetPair(int nLen, CQrPckg& Query, CDbPckg& Db); 
	// search for each seed in a entry of database
	int  ExtendSeq2Set(int nSeed, uint unSeedLen, std::vector<uchar>& vExtra,
			int nQSeqIdx, CAlnPckg& QrAln, int nQOriLen, std::vector<char>& vValid,
			VUINT& vDSet, CDbPckg& Db,
			VNAMES& vQNames, VNAMES& vDNmaes,
			MRESULT& mRes, int nTreadID);

	// align two sequences
	bool AlignSeqs(int nSeed, CAlnPckg& QrAln, CAlnPckg& DbAln,  uint& unSeedLen, STAlnmnt& stAlnmnt, int nTreadID);
	// forward ungapped alignment
	int AlignFwd(uchar *queryseq, uchar *dataseq, uint len_queryseq, uint len_dataseq, int *extl, int *match, int score0);
	// backward ungapped alignment
	int AlignBwd(uchar *queryseq, uchar *dataseq, int pos1, int pos2, int *extl, int *match, int score0);
	// gapped alignment
	int AlignGapped(uchar *seq1, uchar *seq2, int M, int N, int *ext1, int *ext2, int *match_len, int *gap, std::vector<char>& vMode, std::vector<short>& vLen, int nTreadID);
	// retrieve a sequence according to a seed
	uchar* GetSeq(VUCHAR& vSeqs, VUINT& vLen, VNAMES& vNames, uint& unPos, uint& unLen, uint& unSeedBeg);
	
	// reset the struct
	void ResetResult(STAlnmnt& stAlnmnt);

	// calculate e-value and generate the subsequence 
	void CalRes(int nQIdx, uchar* pQ, int nQLen, uint unQSeedBeg, int nDIdx, uchar* pD, uint unDSeedBeg, CDbPckg& Db, uint unLocalSeedLen, STAlnmnt& stAlnmnt, MRESULT& mRes, int nTreadID);
	// calculate e-values of multi hits
	void SumEvalue(std::vector<CHitUnit>& v, int nSt, int nEd, int nLen, int nTreadID);
	// output the result
	void PrintRes(MRESULT& mRes, int nTreadID, CQrPckg& Query, CDbPckg& Db);

	// revise the size of database according to swift
	void GuessTotSeq(const char* szFile, long int& lnSeqNum, long int& lnAaNum);

	// merge the result files
	void MergeRes(int nDbBlockNum, VNAMES& vQNames, std::string& sDbPre);

	// init alignment parameters
	void InitAlignPara(bool bType, long int lnSLen, int nSNum, int nThreadNum);

	void PrintAln(std::vector<CHitUnit>& v, std::ostream& of);
	void PrintM8(std::vector<CHitUnit>& v, std::ostream& of);

	template<class T>
	void PrintXmlLine(const char* sTag, T s);
	void PrintXmlTag(const char* sTag);
	void PrintXmlTagR(const char* sTag);
	void PrintXmlBegin(std::string& sDb);
	void PrintXml(std::vector<CHitUnit>& v, int nIdx);
	void PrintXmlEnd();


private:
	uint m_unMer;
	uint m_unDSize;
	uint m_unQSize;
	bool m_bSeqType;
	// if dna, 6; if protein, 1
	int m_nIdxScl;
	uint m_unTotalIdx;
	int m_nQueryType;

	uchar m_uMask;
	uchar m_uSeg;
	uchar m_aChar2Code[256]; // char -> compressed code
	uchar m_aCode2Char[256]; // compressed code -> char
	uchar m_aCode2Ten[256]; // char -> 10-base

	int m_aSubMatrix[256][256];

	MRESULT m_mRes;

	/* for alignment and calculation */
	bool m_bEvalue;
	bool m_bLogE;
	double m_dThr;
	int GapIni;
	int GapExt;
	int MaxGap;
	HitComptor* m_pComptor;

	BlastStat* m_pBlastSig;
	double	GapExtSCutBits;
	double	GapExtSCut;
	double	UngapExtDropBits;
	double	UngapExtDrop;
	double	GapExtDropBits;
	double	GapExtDrop;
	double	UngapExtSCut;
	int MinMatch4Exp;

	std::vector<std::vector<std::vector<char> > > m_vTrace;
	std::vector<std::vector<std::vector<char> > > m_vETrace;
	std::vector<std::vector<std::vector<char> > > m_vDTrace;

	bool m_bFast;
	uint m_unMutSeedLen;
	VUINT m_vMutation;

	uint m_unTotalSeeds;
	uint m_unTotalQuery;
	uint m_unTotalSubj;

	/* for output */
	long long m_nMaxOut;
	long long m_nMaxM8;
	bool m_bPrintEmpty;
	bool m_bGapExt;
	std::string m_sStartTime;
	std::string m_sQFile;
	std::string m_sDFile;
	std::ofstream m_ofTemp;
	int m_nStdout;

	bool m_bXml;
	std::ofstream m_ofXml;
	uint m_unXmlSp;
	uint m_unXmlCnt;

	std::string m_sOutBase;
	// store the ouput
	std::string m_sOutput;
	std::string m_sM8;
	std::vector<CIndex> m_vOutIdx;
	std::vector<CIndex> m_vM8Idx;
	long long m_llOutCum;
	long long m_llM8Cum;
	int m_nSeqBase;
	std::string m_sLeft;

	// for multithread
	int m_nThreadNum;
	std::vector<BlastStat*> m_vpBlastSig;
	std::vector<int> m_vBlastPt; // -1 means available, 1 means used

	// for test on gap extension
	uint m_unGapExt;
	bool m_bAcc;
	bool m_bHssp;
	int m_nMinLen;
	long int m_lnSeqNum;
	long int m_lnTotalAa;

	// hssp criteria
	std::vector<int> m_vCriteria;
};


inline void CHashSearch::InitAlignPara(bool bType, long int lnSLen, int nSNum, int nThreadNUm)
{
	for (int i = 0; i < nThreadNUm; ++i)
	{
		BlastStat* pBlastSig = NULL;
		pBlastSig = new BlastStat(1, lnSLen, nSNum);
		m_vpBlastSig.push_back(pBlastSig);
	}

	GapIni = GAPINI;
	GapExt = GAPEXT;
	MaxGap = MAXGAP;
	// set up cutoff values
	// fix it!!!
	GapExtSCutBits = 25; //the dropoff (in bits) to invoke gapped alignment
	GapExtSCut = m_vpBlastSig[0]->Bits2RawScoreUngapped(GapExtSCutBits);

	UngapExtDropBits = 7; //in bits
	GapExtDropBits = 15; //in bits
	if (m_bAcc == false)
	{
		if(bType)
		{
			UngapExtSCut = 12; //blastx default for -f option
		}
		else
		{
			UngapExtSCut = 11; //blastp default
		}
	}
	else
	{
		if(bType)
		{
			UngapExtSCut = 37; //blastx default for -f option
		}
		else
		{
			UngapExtSCut = 11; //blastp default
		}
	}

	UngapExtDrop = m_vpBlastSig[0]->Bits2RawScoreUngapped(UngapExtDropBits);
	GapExtDrop = m_vpBlastSig[0]->Bits2RawScoreGapped(GapExtDropBits);

	MinMatch4Exp = 4;
}


inline int CHashSearch::Tran2Ten(const std::vector<uchar>& v, uint nBeg)
{
	if (nBeg >= v.size())
	{
		return -1;
	}
	uint un = 0;
	for (uint i = 0; i < m_unMer; ++i)
	{
		if (m_uMask == m_aCode2Ten[v[nBeg+i]])
		{
			return -1;
		}
		un = un*10 + m_aCode2Ten[v[nBeg+i]];
	}
	return un;
}


inline int CHashSearch::Tran2Ten(CAlnPckg& QrAln, std::vector<char>& vValid)
{
	if (QrAln.m_unSeedBeg >= QrAln.m_unLen-m_unMer+1)
	{
		return -1;
	}
	uint un = 0;
	for (uint i = 0; i < m_unMer; ++i)
	{
		if (m_uMask == vValid[QrAln.m_unSeedBeg+i])
		{
			return -1;
		}
		un = un*10 + m_aCode2Ten[QrAln.m_pSeq[QrAln.m_unSeedBeg+i]];
	}
	return un;
}


inline void CHashSearch::Decode(const std::vector<uchar>& v, std::string& sOut)
{
        // FIXME: BaCh 24.08.2016
        // reserve() ... really??? Shouldn't it be resize() or is there some other "magic"
        //  behind the scenes?
        // That doesn't feel right, check asap.
	sOut.reserve(v.size());
	for (uint i = 0; i < v.size(); ++i)
	{
		sOut += m_aCode2Char[v[i]];
	}
}


inline void CHashSearch::Encode(std::vector<uchar>& v)
{
	//cout << v.size() << endl;
	for (uint i = 0; i < v.size(); ++i)
	{
		v[i] = m_aChar2Code[v[i]];
	}
}


inline long int CHashSearch::Encode(std::vector<uchar>& v, std::vector<double>& vFreq)
{
	long int lnTotalAa = 0;
	//cout << v.size() << endl;
	for (uint i = 0; i < v.size(); ++i)
	{
		v[i] = m_aChar2Code[v[i]];
		if (m_uMask != v[i])
		{
			++vFreq[m_aCode2Ten[v[i]]];
			++lnTotalAa;
		}
	}

	return lnTotalAa;
}


inline uchar* CHashSearch::GetSeq(VUCHAR& vSeqs, VUINT& vLens, VNAMES& vNames, uint& unPos, uint& unLen, uint& unSeedBeg)
{
	uint unIdx = unPos >> 11;
	unSeedBeg = unPos & 0x000007FF;
	unLen = vLens[unIdx+1] - vLens[unIdx];
	return &vSeqs[0] + vLens[unIdx];
}


inline void CHashSearch::ResetResult(STAlnmnt& stAlnmnt)
{
	stAlnmnt.nScore = 0;
	stAlnmnt.nMatch = 0;
	stAlnmnt.nQFwd = 0;
	stAlnmnt.nDFwd = 0;
	stAlnmnt.nQBwd = 0;
	stAlnmnt.nDBwd = 0;
	stAlnmnt.vMode.clear();
	stAlnmnt.vLen.clear();
}




#endif // __HASHSEARCH_H_
