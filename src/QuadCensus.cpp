/*
 * QuadCensus.cpp
 *
 *  Created on: 14 Sep 2017
 *      Author: ortmann
 */

#include "QuadCensus.h"

#include <algorithm>

#include "Graph.h"

namespace oaqc {

QuadCensus::QuadCensus(const unsigned int n, const unsigned int m,
		const int * const edges) :
		_nodeOrbitCount(20), //
		_edgeOrbitCount(14), //
		_eTriCount(0), //
		_nTriCount(0), //
		_nNonIndC4Count(0), //
		_eNonIndC4Count(0), //
		_eOrbits(0), //
		_nOrbits(0), //
		_neighDeg(0), //
		_k3Count(0), //
		_2pathCount(0),//
		 _graph(n,m,edges){
	init();
	initCounts();
	calcK3K4C4();
	calcK3RelNonIndCounts();
	calcNonInducedFrequencies();
}

QuadCensus::~QuadCensus() {
	clear();
}

const unsigned long* QuadCensus::eOrbits() {
	return _eOrbits;
}

const unsigned long* QuadCensus::nOrbits() {
	return _nOrbits;
}

void QuadCensus::init() {
	const unsigned int n = _graph.n();
	const unsigned int m = _graph.m();
	_eTriCount = new unsigned long[m]();
	_nTriCount = new unsigned long[n]();
	_eNonIndC4Count = new unsigned long[m]();
	_nNonIndC4Count = new unsigned long[n]();
	_eOrbits = new unsigned long[_edgeOrbitCount * m]();
	_nOrbits = new unsigned long[_nodeOrbitCount * n]();
	_neighDeg = new unsigned long[n]();
}

void QuadCensus::clear() {
	delete[] _eTriCount;
	delete[] _nTriCount;
	delete[] _nNonIndC4Count;
	delete[] _eNonIndC4Count;
	delete[] _eOrbits;
	delete[] _nOrbits;
	delete[] _neighDeg;
}

void QuadCensus::initCounts() {
	for (unsigned int i = 0; i < _graph.n(); ++i) {
		const unsigned long deg = _graph.lastOutEdge(i) - _graph.firstInEdge(i);
		_2pathCount += choose2(deg);
		N_ORBIT(i, 11)= choose3 (deg);
		const unsigned int endInd = _graph.lastInEdge(i);
		for (unsigned int neigh = _graph.firstInEdge(i); neigh < endInd;
				++neigh) {
			const unsigned int neighInd = _graph.opInd(neigh);
			_neighDeg[i] += _graph.lastOutEdge(neighInd)
					- _graph.firstInEdge(neighInd);
			_neighDeg[neighInd] += deg;
		}
	}
}

void QuadCensus::calcK3K4C4() {
	int* const innerMark = new int[_graph.n()];
	std::fill_n(innerMark, _graph.n(), -1);
	int* const outerMark = new int[_graph.n()];
	std::fill_n(outerMark, _graph.n(), -1);

	unsigned int* const visCount = new unsigned int[_graph.n()]();
	unsigned int* const workCount = new unsigned int[_graph.n()]();

	for (unsigned int t4 = 1; t4 < _graph.n(); ++t4) {

		const unsigned int t4StartIncInd = _graph.firstInEdge(t4);
		const unsigned int t4EndIncInd = _graph.lastInEdge(t4);

		// mark the edges
		for (unsigned int i = t4StartIncInd; i < t4EndIncInd; i++) {
			outerMark[_graph.opInd(i)] = _graph.edgeInd(i);
		}
		// now find K3, K4, C4
		for (unsigned int t1EdgeInd = t4StartIncInd; t1EdgeInd < t4EndIncInd;
				++t1EdgeInd) {
			unsigned int t1 = _graph.opInd(t1EdgeInd);
			unsigned int l = outerMark[t1];
			outerMark[t1] = -1;
			unsigned int t1StartIncInd = _graph.firstInEdge(t1);
			unsigned int t1StartOutInd = _graph.firstOutEdge(t1);
			// for C4 counting
			for (; t1StartIncInd < t1StartOutInd; ++t1StartIncInd) {
				const unsigned int pos = _graph.opInd(t1StartIncInd);
				++visCount[pos];
				++workCount[pos];
			}
			unsigned int t1EndOutInd = t1StartIncInd;
			// for K4 counting
			for (unsigned int i = _graph.opInd(t1EndOutInd); i != t4;
					i = _graph.opInd(++t1EndOutInd)) {
				++visCount[i];
				++workCount[i];
				innerMark[i] = _graph.edgeInd(t1EndOutInd);
			}

			for (unsigned int t2EdgeInd = t1StartOutInd;
					t2EdgeInd < t1EndOutInd; ++t2EdgeInd) {
				const unsigned int t2 = _graph.opInd(t2EdgeInd);
				const unsigned int s1 = _graph.edgeInd(t2EdgeInd);
				innerMark[t2] = -1;
				if (outerMark[t2] == -1) {
					continue;
				}
				const unsigned int s2 = outerMark[t2];

				// this is an triangle
				++_eTriCount[s1];
				++_eTriCount[s2];
				++_eTriCount[l];

				++_nTriCount[t4];
				++_nTriCount[t1];
				++_nTriCount[t2];

				int t2StartOudInd = _graph.firstOutEdge(t2);
				const int t2EndOutInd = _graph.lastOutEdge(t2);
				while (t2StartOudInd < t2EndOutInd) {
					const int t3 = _graph.opInd(t2StartOudInd);
					const int s3 = _graph.edgeInd(t2StartOudInd);
					// this is a k4
					if (outerMark[t3] >= 0 && innerMark[t3] >= 0) {
						++E_ORBIT(s1, 13);
						++E_ORBIT(s2, 13);
						++E_ORBIT(s3, 13);
						++E_ORBIT(l, 13);
						++E_ORBIT(innerMark[t3], 13);
						++E_ORBIT(outerMark[t3], 13);

						++N_ORBIT(t1, 19);
						++N_ORBIT(t2, 19);
						++N_ORBIT(t3, 19);
						++N_ORBIT(t4, 19);
					}
					++t2StartOudInd;
				}
			}
		}
		for (unsigned int t1EdgeInd = t4StartIncInd; t1EdgeInd < t4EndIncInd;
				++t1EdgeInd) {
			const unsigned int t1 = _graph.opInd(t1EdgeInd);
			const unsigned int l = _graph.edgeInd(t1EdgeInd);
			unsigned int t1StartIncInd = _graph.firstInEdge(t1);
			unsigned int t2;
			while ((t2 = _graph.opInd(t1StartIncInd)) != t4) {
				const int c4count = visCount[t2] - 1;
				if (--workCount[t2] == 0) {
					const long val = choose2(c4count + 1);
					_nNonIndC4Count[t4] += val;
					_nNonIndC4Count[t2] += val;
					visCount[t2] = 0;
				}
				_nNonIndC4Count[t1] += c4count;
				_eNonIndC4Count[l] += c4count;
				_eNonIndC4Count[_graph.edgeInd(t1StartIncInd)] += c4count;
				++t1StartIncInd;
			}
		}
	}
	delete[] innerMark;
	delete[] outerMark;
	delete[] visCount;
	delete[] workCount;
}

void QuadCensus::calcK3RelNonIndCounts() {

	int* const mark = new int[_graph.n()];
	std::fill_n(mark, _graph.n(), -1);

	for (unsigned int t3 = 2; t3 < _graph.n(); ++t3) {
		unsigned int t3FirstInInd = _graph.firstInEdge(t3);
		unsigned int t3EndInInd = _graph.lastInEdge(t3);
		for (unsigned int i = t3FirstInInd; i < t3EndInInd; ++i) {
			mark[_graph.opInd(i)] = _graph.edgeInd(i);
		}
		for (unsigned int i = t3FirstInInd; i < t3EndInInd; ++i) {
			const unsigned int t1 = _graph.opInd(i);
			const int l = mark[t1];
			mark[t1] = -1;
			unsigned int t1StartOutInd = _graph.firstOutEdge(t1);
			unsigned int t2;
			while ((t2 = _graph.opInd(t1StartOutInd)) != t3) {
				if (mark[t2] >= 0) {
					++_k3Count;
					const unsigned
					int s1 = _graph.edgeInd(t1StartOutInd);
					const unsigned
					int s2 = mark[t2];

					const unsigned
					long s1T = _eTriCount[s1];
					const unsigned
					long s2T = _eTriCount[s2];
					const unsigned
					long lT = _eTriCount[l];

					const unsigned
					int t1D = _graph.lastOutEdge(t1) - _graph.firstInEdge(t1);
					const unsigned
					int t2D = _graph.lastOutEdge(t2) - _graph.firstInEdge(t2);
					const unsigned
					int t3D = _graph.lastOutEdge(t3) - _graph.firstInEdge(t3);

					// update diamond counts
					E_ORBIT(l, 11)+= s1T + s2T;
					E_ORBIT(s1, 11)+= s2T + lT;
					E_ORBIT(s2, 11)+= s1T + lT;

					N_ORBIT(t1, 17)+= s2T;
					N_ORBIT(t2, 17)+= lT;
					N_ORBIT(t3, 17)+= s1T;

					// update paw counts
					E_ORBIT(s1, 9)+= t3D;
					E_ORBIT(s2, 9)+= t1D;
					E_ORBIT(l, 9)+= t2D;

					N_ORBIT(t1, 14)+= t2D + t3D;
					N_ORBIT(t2, 14)+= t1D + t3D;
					N_ORBIT(t3, 14)+= t1D + t2D;
				}
				++t1StartOutInd;
			}
		}
	}
	delete[] mark;
}

void QuadCensus::calcNonInducedFrequencies() {
	// for each edge
	for (unsigned int i = 0; i < _graph.n(); i++) {
		const unsigned int endIndex = _graph.lastInEdge(i);
		const unsigned int srcDeg = _graph.lastOutEdge(i) - _graph.firstInEdge(i);
		for (unsigned int neighPos = _graph.firstInEdge(i); neighPos < endIndex;
				neighPos++) {
			const unsigned int tgtIndex = _graph.opInd(neighPos);
			const unsigned int eIndex = _graph.edgeInd(neighPos);
			const unsigned int tgtDeg = _graph.lastOutEdge(tgtIndex)
					- _graph.firstInEdge(tgtIndex);

			// set non-induced edge orbit counts
			E_ORBIT(eIndex, 12)= choose2(_eTriCount[eIndex] );
			E_ORBIT(eIndex, 11)= E_ORBIT(eIndex, 11) - 2 * _eTriCount[eIndex];
			E_ORBIT(eIndex, 10)=_eNonIndC4Count[eIndex];
			E_ORBIT(eIndex, 9)= E_ORBIT(eIndex, 9) - 2 * _eTriCount[eIndex];
			E_ORBIT(eIndex, 8)= _eTriCount[eIndex] * (srcDeg + tgtDeg - 4);
			E_ORBIT(eIndex, 7)=
			_nTriCount[i] + _nTriCount[tgtIndex] - 2 * _eTriCount[eIndex];
			E_ORBIT(eIndex, 6)= choose2(srcDeg - 1) + choose2(tgtDeg - 1);
			E_ORBIT(eIndex, 5)=(srcDeg - 1) * (tgtDeg - 1) - _eTriCount[eIndex];
			E_ORBIT(eIndex, 4)= _neighDeg[i] + _neighDeg[tgtIndex] - 2 * (srcDeg + tgtDeg) + 2 - 2 * _eTriCount[eIndex];
			E_ORBIT(eIndex, 3)= _eTriCount[eIndex] * (_graph.n() - 3);
			E_ORBIT(eIndex, 2)=((srcDeg - 1) + (tgtDeg - 1)) * (_graph.n() - 3);
			E_ORBIT(eIndex, 1)= _graph.m() - srcDeg - tgtDeg + 1;
			E_ORBIT(eIndex, 0)= choose2(_graph.n() - 2);

			// set non-induced node orbit counts
			N_ORBIT(i, 18)+= choose2(_eTriCount[eIndex]);
			N_ORBIT(i, 15)+= _nTriCount[tgtIndex] - _eTriCount[eIndex];
			N_ORBIT(i, 12)+= choose2(tgtDeg - 1);
			N_ORBIT(i, 10)+= _neighDeg[tgtIndex] - tgtDeg;
			N_ORBIT(i, 9)+= (srcDeg - 1) * (tgtDeg - 1) - _eTriCount[eIndex];
			N_ORBIT(i, 6)+= tgtDeg - 1;
			N_ORBIT(i, 3)+= _graph.m() - (srcDeg + tgtDeg - 1);

			N_ORBIT(tgtIndex, 18)+= choose2(_eTriCount[eIndex]);
			N_ORBIT(tgtIndex, 15)+= _nTriCount[i] - _eTriCount[eIndex];
			N_ORBIT(tgtIndex, 12)+= choose2(srcDeg - 1);
			N_ORBIT(tgtIndex, 10)+= _neighDeg[i] - srcDeg;
			N_ORBIT(tgtIndex, 9)+=
			(srcDeg - 1) * (tgtDeg - 1) - _eTriCount[eIndex];
			N_ORBIT(tgtIndex, 6)+= srcDeg - 1;
			N_ORBIT(tgtIndex, 3)+= _graph.m() - (srcDeg + tgtDeg - 1);
		}
	}
	// solve remaining non-induced node frequencies
	const unsigned
	long allPairs = choose3(_graph.n() - 1);
	for (unsigned int i = 0; i < _graph.n(); ++i) {
		const unsigned int deg = _graph.lastOutEdge(i) - _graph.firstInEdge(i);

		N_ORBIT(i, 17)-= _nTriCount[i];
		N_ORBIT(i, 16)= _nNonIndC4Count[i];
		N_ORBIT(i, 14)-= _nTriCount[i] * 4;
		N_ORBIT(i, 13)= _nTriCount[i] * (deg - 2);
		N_ORBIT(i, 11)= choose3(deg);
		N_ORBIT(i, 10)-= deg * (deg - 1) + 2 * _nTriCount[i];
		N_ORBIT(i, 8)= _k3Count - _nTriCount[i];
		N_ORBIT(i, 7)= _nTriCount[i] * (_graph.n() - 3);
		N_ORBIT(i, 5)= N_ORBIT(i, 6);
		N_ORBIT(i, 6)= _2pathCount - N_ORBIT(i, 6) - choose2(deg);
		N_ORBIT(i, 5)= N_ORBIT(i, 5) * (_graph.n() - 3);
		N_ORBIT(i, 4)= choose2(deg) * (_graph.n() - 3);
		N_ORBIT(i, 2)= (_graph.m() - deg) * (_graph.n() - 3);
		N_ORBIT(i, 1)= choose2(_graph.n() - 2) * deg;
		N_ORBIT(i, 0)= allPairs;
	}

}

void QuadCensus::calcInducedFrequencies() {
	for (unsigned int i = 0; i < _graph.n(); ++i) {
		const unsigned int endIndex = _graph.lastInEdge(i);
		for (unsigned int neighPos = _graph.firstInEdge(i); neighPos < endIndex;
				++neighPos) {
			const unsigned int eIndex = _graph.edgeInd(neighPos);

			E_ORBIT(eIndex, 12)-= E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 11)-= 4 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 10)-=
			E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 9)-=
			E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 8)-=
			E_ORBIT(eIndex, 11) + 4 * E_ORBIT(eIndex, 12)
			+ 4 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 7)-=
			E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 6)-=
			E_ORBIT(eIndex, 7) + E_ORBIT(eIndex, 8)
			+ E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 12)
			+ 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 5)-=
			E_ORBIT(eIndex, 8) + E_ORBIT(eIndex, 10)
			+ E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 12)
			+ 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 4)-=
			2 * E_ORBIT(eIndex, 7) + 2 * E_ORBIT(eIndex, 9)
			+ 2 * E_ORBIT(eIndex, 10)
			+ 3 * E_ORBIT(eIndex, 11)
			+ 4 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 3)-=
			E_ORBIT(eIndex, 8) + E_ORBIT(eIndex, 9)
			+ E_ORBIT(eIndex, 11) + 2 * E_ORBIT(eIndex, 12)
			+ 2 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 2)-=
			2 * E_ORBIT(eIndex, 3) + E_ORBIT(eIndex, 4)
			+ 2 * E_ORBIT(eIndex, 5)
			+ 2 * E_ORBIT(eIndex, 6)
			+ 2 * E_ORBIT(eIndex, 7)
			+ 3 * E_ORBIT(eIndex, 8)
			+ 2 * E_ORBIT(eIndex, 9)
			+ 2 * E_ORBIT(eIndex, 10)
			+ 3 * E_ORBIT(eIndex, 11)
			+ 4 * E_ORBIT(eIndex, 12)
			+ 4 * E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 1)-=
			E_ORBIT(eIndex, 4) + E_ORBIT(eIndex, 7)
			+ E_ORBIT(eIndex, 9) + E_ORBIT(eIndex, 10)
			+ E_ORBIT(eIndex, 11) + E_ORBIT(eIndex, 13);

			E_ORBIT(eIndex, 0)-=
			E_ORBIT(eIndex, 1) + E_ORBIT(eIndex, 2)
			+ E_ORBIT(eIndex, 3) + E_ORBIT(eIndex, 4)
			+ E_ORBIT(eIndex, 5) + E_ORBIT(eIndex, 6)
			+ E_ORBIT(eIndex, 7) + E_ORBIT(eIndex, 8)
			+ E_ORBIT(eIndex, 9) + E_ORBIT(eIndex, 10)
			+ E_ORBIT(eIndex, 11) + E_ORBIT(eIndex, 12)
			+ E_ORBIT(eIndex, 13);

		}
	}

	for (unsigned int i = 0; i < _graph.n(); ++i) {
		N_ORBIT(i, 18)-= 3 * N_ORBIT(i, 19);

		N_ORBIT(i, 17)-= 3 * N_ORBIT(i, 19);

		N_ORBIT(i, 16)-=
		3 * N_ORBIT(i, 19) + N_ORBIT(i, 18) + N_ORBIT(i, 17);

		N_ORBIT(i, 15)-= 3 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 17);

		N_ORBIT(i, 14)-=
		6 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18)
		+ 2 * N_ORBIT(i, 17);

		N_ORBIT(i, 13)-= 3 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18);

		N_ORBIT(i, 12)-=
		3 * N_ORBIT(i, 19) + N_ORBIT(i, 18) + 2 * N_ORBIT(i, 17)
		+ N_ORBIT(i, 15) + N_ORBIT(i, 14);

		N_ORBIT(i, 11)-=
		N_ORBIT(i, 19) + N_ORBIT(i, 18) + N_ORBIT(i, 13);

		N_ORBIT(i, 10)-=
		6 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18)
		+ 4 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ 2 * N_ORBIT(i, 15) + N_ORBIT(i, 14);

		N_ORBIT(i, 9)-=
		6 * N_ORBIT(i, 19) + 4 * N_ORBIT(i, 18)
		+ 2 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ N_ORBIT(i, 14) + 2 * N_ORBIT(i, 13);

		N_ORBIT(i, 8)-=
		N_ORBIT(i, 19) + N_ORBIT(i, 17) + N_ORBIT(i, 15);

		N_ORBIT(i, 7)-=
		3 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18) + N_ORBIT(i, 17)
		+ N_ORBIT(i, 14) + N_ORBIT(i, 13);

		N_ORBIT(i, 6)-=
		3 * N_ORBIT(i, 19) + N_ORBIT(i, 18) + 3 * N_ORBIT(i, 17)
		+ N_ORBIT(i, 16) + 3 * N_ORBIT(i, 15)
		+ N_ORBIT(i, 14) + N_ORBIT(i, 12) + N_ORBIT(i, 10)
		+ 3 * N_ORBIT(i, 8);

		N_ORBIT(i, 5)-=
		6 * N_ORBIT(i, 19) + 4 * N_ORBIT(i, 18)
		+ 4 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ 2 * N_ORBIT(i, 15) + 3 * N_ORBIT(i, 14)
		+ 2 * N_ORBIT(i, 13) + 2 * N_ORBIT(i, 12)
		+ N_ORBIT(i, 10) + N_ORBIT(i, 9)
		+ 2 * N_ORBIT(i, 7);

		N_ORBIT(i, 4)-=
		3 * N_ORBIT(i, 19) + 3 * N_ORBIT(i, 18) + N_ORBIT(i, 17)
		+ N_ORBIT(i, 16) + N_ORBIT(i, 14)
		+ 3 * N_ORBIT(i, 13) + 3 * N_ORBIT(i, 11)
		+ N_ORBIT(i, 9) + N_ORBIT(i, 7);

		N_ORBIT(i, 3)-=
		3 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18)
		+ 2 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ N_ORBIT(i, 15) + N_ORBIT(i, 14) + N_ORBIT(i, 13)
		+ N_ORBIT(i, 10) + N_ORBIT(i, 9);

		N_ORBIT(i, 2)-=
		3 * N_ORBIT(i, 19) + 2 * N_ORBIT(i, 18)
		+ 3 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ 3 * N_ORBIT(i, 15) + 2 * N_ORBIT(i, 14)
		+ N_ORBIT(i, 13) + 2 * N_ORBIT(i, 12)
		+ 2 * N_ORBIT(i, 10) + N_ORBIT(i, 9)
		+ 3 * N_ORBIT(i, 8) + N_ORBIT(i, 7)
		+ 2 * N_ORBIT(i, 6) + N_ORBIT(i, 5)
		+ N_ORBIT(i, 3);

		N_ORBIT(i, 1)-=
		3 * N_ORBIT(i, 19) + 3 * N_ORBIT(i, 18)
		+ 2 * N_ORBIT(i, 17) + 2 * N_ORBIT(i, 16)
		+ N_ORBIT(i, 15) + 2 * N_ORBIT(i, 14)
		+ 3 * N_ORBIT(i, 13) + N_ORBIT(i, 12)
		+ 3 * N_ORBIT(i, 11) + N_ORBIT(i, 10)
		+ 2 * N_ORBIT(i, 9) + 2 * N_ORBIT(i, 7)
		+ N_ORBIT(i, 5) + 2 * N_ORBIT(i, 4)
		+ N_ORBIT(i, 3);

		N_ORBIT(i, 0)-=
		N_ORBIT(i, 19) + N_ORBIT(i, 18) + N_ORBIT(i, 17)
		+ N_ORBIT(i, 16) + N_ORBIT(i, 15) + N_ORBIT(i, 14)
		+ N_ORBIT(i, 13) + N_ORBIT(i, 12) + N_ORBIT(i, 11)
		+ N_ORBIT(i, 10) + N_ORBIT(i, 9) + N_ORBIT(i, 8)
		+ N_ORBIT(i, 7) + N_ORBIT(i, 6) + N_ORBIT(i, 5)
		+ N_ORBIT(i, 4) + N_ORBIT(i, 3) + N_ORBIT(i, 2)
		+ N_ORBIT(i, 1);
	}
}

}
/* namespace oaqc */
