/*
 * Graph.cpp
 *
 *  Created on: 14 Sep 2017
 *      Author: ortmann
 */

#include "Graph.h"
#include <algorithm>


namespace oaqc {

Graph::Graph(const unsigned int n, const unsigned int m, const int* const edges) :
		_edges(new Edge[2 * m]), //
		_n(n), //
		_m(m), //
		_inOffset(new unsigned int[n + 1]), //
		_outOffset(new unsigned int[n]), //
		_mapping(new unsigned int[n]) {
	// used to find the highest indexed outgoing edge
	_inOffset[_n] = 2 * _m;
	createGraph(edges);
}

Graph::~Graph() {
	delete[] _mapping;
	delete[] _outOffset;
	delete[] _inOffset;
	delete[] _edges;
}

void Graph::createGraph(const int* const edges) {
	bucketSort(edges);

	// store the graph
	for (unsigned int i = 0; i < _m; ++i) {
		const unsigned int n1 = _mapping[edges[i]];
		const unsigned int n2 = _mapping[edges[i + _m]];
		_edges[_outOffset[n1]].first = n2;
		_edges[_outOffset[n1]].second = i;
		++_outOffset[n1];
		_edges[_outOffset[n2]].first = n1;
		_edges[_outOffset[n2]].second = i;
		++_outOffset[n2];
	}

	// sort the edges according the opposite indices
	for (unsigned int i = 0; i < _n; ++i) {
		// this works since all ids are unique => second values are never compared
		// set the outOffset
		std::sort(_edges + _inOffset[i], _edges + _inOffset[i + 1]);
		for (unsigned int pos = _inOffset[i]; pos < _inOffset[i + 1]; ++pos) {
			if (_edges[pos].first > i) {
				_outOffset[i] = pos;
				break;
			}
		}
	}
}

void Graph::bucketSort(const int* const edges) {
	// calculate node degrees
	unsigned int* const deg = new unsigned int[_n]();
	for (unsigned int i = 0; i < _m; ++i) {
		++deg[edges[i]];
		++deg[edges[i + _m]];
	}
	// init the buckets array
	unsigned int maxDegree = 0;
	for (unsigned int i = 0; i < _n; ++i) {
		maxDegree = std::max(maxDegree, deg[i]);

	}
	unsigned int* const bucket = new unsigned int[maxDegree + 1]();

	// compute size of each bucket
	for (unsigned int i = 0; i < _n; ++i) {
		++bucket[deg[i]];
	}
	unsigned int first = 0;
	unsigned int size;
	// calculate the lower end point of each bucket
	for (unsigned int i = 0; i <= maxDegree; ++i) {
		size = bucket[i];
		bucket[i] = first;
		first += size;
	}
	unsigned int* const reverseMapping = new unsigned int[_n];
	for (unsigned int i = 0; i < _n; ++i) {
		const unsigned int degree = deg[i];
		unsigned int pos = bucket[degree];
		_mapping[i] = pos;
		reverseMapping[pos] = i;
		++bucket[degree];
	}

	unsigned int m = 0;
	for (unsigned int i = 0; i < _n; ++i) {
		_inOffset[i] = m;
		_outOffset[i] = m;
		m += deg[reverseMapping[i]];
	}
	delete[] bucket;
	delete[] deg;
	delete[] reverseMapping;
}

} /* namespace oaqc */
