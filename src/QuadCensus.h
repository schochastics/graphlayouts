/*
 * QuadCensus.h
 *
 *  Created on: 14 Sep 2017
 *      Author: ortmann
 */

#ifndef QUADCENSUS_H_
#define QUADCENSUS_H_

#include "Graph.h"

namespace oaqc {

#define N_ORBIT(nIndex,orbit) _nOrbits[nIndex * _nodeOrbitCount + orbit]
#define E_ORBIT(eIndex,orbit) _eOrbits[eIndex * _edgeOrbitCount + orbit]

class QuadCensus {
public:
	QuadCensus(const unsigned int n, const unsigned int m,
			const int * const edges);
	virtual ~QuadCensus();

	const unsigned long getNOrbitCount() const{
		return _nodeOrbitCount;
	}
	const unsigned long getEOrbitCount() const {
		return _edgeOrbitCount;
	}

	const unsigned int* getMapping() const{
		return _graph.getMapping();
	}

	const unsigned long* eOrbits();
	const unsigned long* nOrbits();

	void calcInducedFrequencies();


private:


	const unsigned long _nodeOrbitCount;
	const unsigned long _edgeOrbitCount;


	unsigned long* _eTriCount;
	unsigned long* _nTriCount;
	unsigned long* _nNonIndC4Count;
	unsigned long* _eNonIndC4Count;
	unsigned long* _eOrbits;
	unsigned long* _nOrbits;
	unsigned long* _neighDeg;

	unsigned long long _k3Count;
	unsigned long long _2pathCount;

	const Graph _graph;

	void init();
	void clear();

	void initCounts();
	void calcK3K4C4();
	void calcK3RelNonIndCounts();
	void calcNonInducedFrequencies();

	inline unsigned long choose2(const unsigned long val) {
		if (val < 1) {
			return 0;
		}
		return (val * (val - 1)) / 2;
	}

	inline unsigned long choose3(const unsigned long val) {
		if (val < 3) {
			return 0;
		}
		return (val * (val - 1) * (val - 2)) / 6;
	}

};

} /* namespace oaqc */

#endif /* QUADCENSUS_H_ */
