// this class is responsible for merge the temp results files
// by yongzhao
//
//
//


#include <algorithm>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/serialization/vector.hpp>
#include "mergeUnit.h"

using std::cout;
using std::cerr;
using std::endl;



CMergeUnit::CMergeUnit(const char* szFile)
{
	m_ifFile.open(szFile);
	if (!m_ifFile.good())
	{
		cout << "can not open temp file..." << endl;
	}
	
	m_sFile.assign(szFile);

	std::ifstream ifIdx((m_sFile+".idx").c_str());
	boost::archive::binary_iarchive ia(ifIdx);
	ia >> m_vIdx;
	ifIdx.close();
}


void CMergeUnit::Update(int nID, std::vector<CHitUnit>& v)
{
	if (m_vIdx[nID].m_nSize == 0)
	{
		return;
	}

	m_ifFile.seekg(m_vIdx[nID].m_llBeg, std::ios::beg);
	boost::archive::binary_iarchive ia(m_ifFile);
	std::vector<CHitUnit> vIn;
	ia >> vIn;

	v.insert(v.end(), vIn.begin(), vIn.end());
}

int CMergeUnit::GetLast()
{
	return m_vIdx.size();
}
