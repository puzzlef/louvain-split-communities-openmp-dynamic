#pragma once
#include <utility>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "Graph.hxx"
#include "properties.hxx"
#include "louvain.hxx"
#include "split.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::vector;
using std::swap;
using std::max;




#pragma region ENVIRONMENT SETUP
#ifdef OPENMP
/**
 * Setup and perform the Louvain algorithm.
 * @param x original graph
 * @param o louvain options
 * @param fi initializing community membership and total vertex/community weights (vcom, vtot, ctot)
 * @param fm marking affected vertices (vaff, vcs, vcout, vcom, vtot, ctot)
 * @param fa is vertex allowed to be updated? (u)
 * @returns louvain result
 */
template <int SPLIT_CHUNK=2048, bool DYNAMIC=false, class G, class FI, class FM, class FA>
inline auto louvainSplitIterationInvokeOmp(const G& x, const LouvainOptions& o, FI fi, FM fm, FA fa) {
  using  K = typename G::key_type;
  using  W = LOUVAIN_WEIGHT_TYPE;
  using  B = char;
  // Options.
  double R = o.resolution;
  int    L = o.maxIterations, l = 0;
  int    P = o.maxPasses, p = 0;
  // Get graph properties.
  size_t X = x.size();
  size_t S = x.span();
  double M = edgeWeightOmp(x)/2;
  // Allocate buffers.
  int    T = omp_get_max_threads();
  vector<B> vaff(S);        // Affected vertex flag (any pass)
  vector<K> ucom, vcom(S);  // Community membership (first pass, current pass)
  vector<W> utot, vtot(S);  // Total vertex weights (first pass, current pass)
  vector<W> ctot;           // Total community weights (any pass)
  vector<K> bufk(T);        // Buffer for exclusive scan
  vector<size_t> bufs(T);   // Buffer for exclusive scan
  vector<vector<K>*> vcs(T);    // Hashtable keys
  vector<vector<W>*> vcout(T);  // Hashtable values
  if (!DYNAMIC) ucom.resize(S);
  if (!DYNAMIC) utot.resize(S);
  if (!DYNAMIC) ctot.resize(S);
  louvainAllocateHashtablesW(vcs, vcout, S);
  size_t Z = max(size_t(o.aggregationTolerance * X), X);
  size_t Y = max(size_t(o.aggregationTolerance * Z), Z);
  DiGraphCsr<K, None, None, K> cv(S, S);  // CSR for community vertices
  DiGraphCsr<K, None, W> y(S, Y);         // CSR for aggregated graph (input);  y(S, X)
  DiGraphCsr<K, None, W> z(S, Z);         // CSR for aggregated graph (output); z(S, X)
  // Data structures for splitting disconnected communities.
  vector<K> tcom(S);                      // Temporary community membership
  vector<vector<K>*> us(T), vs(T);        // Per-thread start, frontier vertices for BFS
  for (int t=0; t<T; ++t) {
    us[t] = new vector<K>();
    vs[t] = new vector<K>();
    (*us[t]).reserve(4*S/T);
    (*vs[t]).reserve(4*S/T);
  }
  // Perform Louvain algorithm.
  float tm = 0, ti = 0, tp = 0, tl = 0, ta = 0, ts = 0;  // Time spent in different phases
  float t  = measureDurationMarked([&](auto mark) {
    double E  = o.tolerance;
    auto   fc = [&](double el, int l) { return el<=E; };
    // Reset buffers, in case of multiple runs.
    fillValueOmpU(vaff, B());
    fillValueOmpU(ucom, K());
    fillValueOmpU(vcom, K());
    fillValueOmpU(utot, W());
    fillValueOmpU(vtot, W());
    fillValueOmpU(ctot, W());
    cv.respan(S);
    y .respan(S);
    z .respan(S);
    // Time the algorithm.
    mark([&]() {
      // Initialize community membership and total vertex/community weights.
      ti += measureDuration([&]() { fi(ucom, utot, ctot); });
      // Mark affected vertices.
      tm += measureDuration([&]() { fm(vaff, vcs, vcout, ucom, utot, ctot); });
      // Start timing first pass.
      auto t0 = timeNow(), t1 = t0;
      // Start local-moving, aggregation phases.
      // NOTE: In first pass, the input graph is a DiGraph.
      // NOTE: For subsequent passes, the input graph is a DiGraphCsr (optimization).
      for (l=0, p=0; M>0 && P>0;) {
        if (p==1) t1 = timeNow();
        bool isFirst = p==0;
        int m = 0;
        tl += measureDuration([&]() {
          if (isFirst) m = louvainMoveOmpW(ucom, ctot, vaff, vcs, vcout, x, utot, M, R, L, fc, fa);
          else         m = louvainMoveOmpW(vcom, ctot, vaff, vcs, vcout, y, vtot, M, R, L, fc);
        });
        ts += measureDuration([&]() {
          // TODO: Split only touched (which nodes leave) communities.
          if (isFirst) { splitDisconnectedCommunitiesBfsOmpW<SPLIT_CHUNK>(tcom, vaff, us, vs, x, ucom); swap(ucom, tcom); }
          else         { splitDisconnectedCommunitiesBfsOmpW<SPLIT_CHUNK>(tcom, vaff, us, vs, y, vcom); swap(vcom, tcom); }
        });
        l += max(m, 1); ++p;
        if (m<=1 || p>=P) break;
        size_t GN = isFirst? x.order() : y.order();
        size_t GS = isFirst? x.span()  : y.span();
        size_t CN = 0;
        if (isFirst) CN = louvainCommunityExistsOmpW(cv.degrees, x, ucom);
        else         CN = louvainCommunityExistsOmpW(cv.degrees, y, vcom);
        if (double(CN)/GN >= o.aggregationTolerance) break;
        if (isFirst) louvainRenumberCommunitiesOmpW(ucom, cv.degrees, bufk, x);
        else         louvainRenumberCommunitiesOmpW(vcom, cv.degrees, bufk, y);
        if (isFirst) {}
        else         louvainLookupCommunitiesOmpU(ucom, vcom);
        ta += measureDuration([&]() {
          cv.respan(CN); z.respan(CN);
          if (isFirst) louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, x, ucom);
          else         louvainCommunityVerticesOmpW(cv.offsets, cv.degrees, cv.edgeKeys, bufk, y, vcom);
          if (isFirst) louvainAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, x, ucom, cv.offsets, cv.edgeKeys);
          else         louvainAggregateOmpW(z.offsets, z.degrees, z.edgeKeys, z.edgeValues, bufs, vcs, vcout, y, vcom, cv.offsets, cv.edgeKeys);
        });
        swap(y, z);
        // fillValueOmpU(vcom.data(), CN, K());
        // fillValueOmpU(ctot.data(), CN, W());
        fillValueOmpU(vtot.data(), CN, W());
        fillValueOmpU(vaff.data(), CN, B(1));
        louvainVertexWeightsOmpW(vtot, y);
        louvainInitializeOmpW(vcom, ctot, y, vtot);
        E /= o.toleranceDrop;
      }
      if (p<=1) {}
      else      louvainLookupCommunitiesOmpU(ucom, vcom);
      if (p<=1) t1 = timeNow();
      tp += duration(t0, t1);
    });
  }, o.repeat);
  // Clean up data structures for splitting disconnected communities.
  for (int t=0; t<T; ++t) {
    delete us[t];
    delete vs[t];
  }
  louvainFreeHashtablesW(vcs, vcout);
  return LouvainResult<K, W>(ucom, utot, ctot, l, p, t, tm/o.repeat, ti/o.repeat, tp/o.repeat, tl/o.repeat, ta/o.repeat, ts/o.repeat, countValueOmp(vaff, B(1)));
}
#endif
#pragma endregion




#pragma region STATIC APPROACH
#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Static Louvain.
 * @param x original graph
 * @param o louvain options
 * @returns louvain result
 */
template <int SPLIT_CHUNK=2048, class G>
inline auto louvainSplitStaticOmp(const G& x, const LouvainOptions& o={}) {
  using B = char;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    louvainVertexWeightsOmpW(vtot, x);
    louvainInitializeOmpW(vcom, ctot, x, vtot);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueOmpU(vaff, B(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainSplitInvokeOmp<SPLIT_CHUNK, false>(x, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region NAIVE-DYNAMIC APPROACH
#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Naive-dynamic Louvain.
 * @param y updated graph
 * @param deletions edge deletions for this batch update (undirected)
 * @param insertions edge insertions for this batch update (undirected)
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <int SPLIT_CHUNK=2048, class G, class K, class V, class W>
inline auto louvainSplitNaiveDynamicOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  using B = char;
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot)  {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [ ](auto& vaff, const auto& vcom, const auto& vtot, const auto& ctot, auto& vcs,  auto& vcout) {
    fillValueOmpU(vaff, B(1));
  };
  auto fa = [ ](auto u) { return true; };
  return louvainSplitInvokeOmp<SPLIT_CHUNK, true>(y, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region DYNAMIC DELTA-SCREENING
#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Delta-screening Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <int SPLIT_CHUNK=2048, class G, class K, class V, class W>
inline auto louvainSplitDynamicDeltaScreeningOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  using  B = char;
  size_t S = y.span();
  double R = o.resolution;
  double M = edgeWeightOmp(y)/2;
  int    T = omp_get_max_threads();
  vector<B> vertices(S), neighbors(S), communities(S);
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesDeltaScreeningOmpW(vertices, neighbors, communities, vcs, vcout, y, deletions, insertions, vcom, vtot, ctot, M, R);
    copyValuesOmpW(vaff, vertices);
  };
  auto fa = [&](auto u) { return vertices[u] == B(1); };
  return louvainSplitInvokeOmp<SPLIT_CHUNK, true>(y, o, fi, fm, fa);
}
#endif
#pragma endregion




#pragma region DYNAMIC FRONTIER APPROACH
#ifdef OPENMP
/**
 * Obtain the community membership of each vertex with Dynamic Frontier Louvain.
 * @param y updated graph
 * @param deletions edge deletions in batch update
 * @param insertions edge insertions in batch update
 * @param q initial community each vertex belongs to
 * @param qvtot initial total edge weight of each vertex
 * @param qctot initial total edge weight of each community
 * @param o louvain options
 * @returns louvain result
 */
template <int SPLIT_CHUNK=2048, class G, class K, class V, class W>
inline auto louvainDynamicFrontierOmp(const G& y, const vector<tuple<K, K, V>>& deletions, const vector<tuple<K, K, V>>& insertions, const vector<K>& q, const vector<W>& qvtot, const vector<W>& qctot, const LouvainOptions& o={}) {
  vector2d<K> qs;
  vector2d<W> qvtots, qctots;
  louvainSetupInitialsW(qs, qvtots, qctots, q, qvtot, qctot, o.repeat);
  int  r  = 0;
  auto fi = [&](auto& vcom, auto& vtot, auto& ctot) {
    vcom = move(qs[r]);
    vtot = move(qvtots[r]);
    ctot = move(qctots[r]); ++r;
    louvainUpdateWeightsFromOmpU(vtot, ctot, y, deletions, insertions, vcom);
  };
  auto fm = [&](auto& vaff, auto& vcs, auto& vcout, const auto& vcom, const auto& vtot, const auto& ctot) {
    louvainAffectedVerticesFrontierOmpW(vaff, y, deletions, insertions, vcom);
  };
  auto fa = [ ](auto u) { return true; };
  return louvainSplitInvokeOmp<SPLIT_CHUNK, true>(y, o, fi, fm, fa);
}
#endif
#pragma endregion
#pragma endregion
