/*
 * Graph.h
 *
 *  Created on: 14 Sep 2017
 *      Author: ortmann
 */
#ifndef GRAPH_H_
#define GRAPH_H_

#include <utility>

typedef std::pair< unsigned int, unsigned int> Edge;

namespace oaqc {

class Graph {
	public:
		Graph(const unsigned int n, const unsigned int m,
				const int * const edges);
		virtual ~Graph();

		Edge* const _edges; // TODO private

		inline const unsigned int* const getMapping() const{
			return _mapping;
		}

		inline unsigned int n() const {
			return _n;
		}

		inline unsigned int m() const {
			return _m;
		}

		inline unsigned int firstInEdge(const unsigned int v) const {
			return _inOffset[v];
		}

		inline unsigned int lastInEdge(const unsigned int v) const {
			return _outOffset[v];
		}

		inline unsigned int firstOutEdge(const unsigned int v) const {
			return _outOffset[v];
		}

		inline unsigned int lastOutEdge(const unsigned int v) const {
			return _inOffset[v + 1];
		}

		inline unsigned int opInd(const unsigned int pos) const {
			return _edges[pos].first;
		}

		inline unsigned int edgeInd(const unsigned int pos) const {
			return _edges[pos].second;
		}

	private:

		unsigned const int _n;
		unsigned const int _m;
		unsigned int* const _inOffset;
		unsigned int* const _outOffset;
		unsigned int* const _mapping;

		void createGraph(const int * const edges);
		void bucketSort(const int * const edges);

};

} /* namespace oaqc */

#endif /* GRAPH_H_ */
