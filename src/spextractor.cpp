// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: OpenMS-the reference implementation project $
// --------------------------------------------------------------------------

#include "MzPeakStreamLoad.h"   // [mzpeak] streaming .mzpeak input (SPEXTRACTOR_WITH_MZPEAK)
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <sys/resource.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <unordered_map>
#include <cstring>
// spextractor::TdfMzCalibration comes from <OpenMS/FORMAT/TdfMzCalibration.h>, installed by the
// OpenMS patch and already included via the mzPeak loader; including the in-tree copy as well
// redefines the struct.
#include <OpenMS/FORMAT/TdfMzCalibration.h>
#include "TdfLoad.h"   // loadTdfCalibration: sqlite3, tool-only
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <type_traits>
#include <array>
#include <functional>
#include <cstdio>
#include <fstream>
#include <atomic>
#include <thread>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>     // getenv
#include <cstring>     // memcpy
#ifdef __GLIBC__
#include <malloc.h>    // mallinfo2: retained-vs-live at the milestones [mem]
#endif

// The three OpenMS patches in patches/ are build prerequisites. This is the compile-time check for the
// third: without defaulted moves an rvalue MassTrace binds to the copy constructor (not noexcept), and
// every "move" of a trace in the band/split gathering is a deep copy of its points.
static_assert(std::is_nothrow_move_constructible_v<OpenMS::MassTrace> && std::is_nothrow_move_assignable_v<OpenMS::MassTrace>,
              "OpenMS MassTrace has no noexcept move operations: apply patches/openms-masstrace-move.patch to the OpenMS tree");
#include <exception>   // exception_ptr
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SpeXtractor SpeXtractor

  @brief Extracts pseudo-MS/MS ("pseudo-DDA") spectra from diaPASEF (ion-mobility DIA) data.

  Reconstructs DDA-like MS2 spectra from ion-mobility DIA (timsTOF / diaPASEF) data by
  exploiting the fact that, in PASEF, fragments are recorded at the ion-mobility (1/K0)
  elution coordinate of their precursor. Precursor and fragment mass traces are detected
  IM-aware (@ref MassTraceDetection with an ion-mobility tolerance), then fragments are
  assigned to a precursor when they co-localize in ion mobility, retention time, and
  elution-profile (Pearson) correlation. The result is written as searchable mzML.

  This is the OpenMS-native, BSD-3 analogue of Nesvilab the reference implementation. It is a discovery
  front-end (enables open / semi-tryptic / PTM searches); it does not perform the database
  search, FDR, library building or quantification.

  See docs/dataset D-BASELINE.md for the measured decisions behind every default, and docs/reviews/
  for the adversarial reviews.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SpeXtractor.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SpeXtractor.html
*/

/// @cond TOPPCLASSES

namespace
{
  /// A detected mass trace reduced to the quantities we correlate/gate on.
  /// Current + peak resident set, read from /proc/self/status (Linux). Logged at phase boundaries so
  /// the memory peak is MEASURED rather than estimated: it tells us whether the peak sits in the
  /// raw-load phase or in the parallel per-window phase. Returns "" where unsupported. [mem]
  /// seconds since tool start, for per-phase bottleneck attribution [perf]
  static double phase_clock_()
  {
    static const auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  static std::string clk_() { char b[32]; snprintf(b, sizeof b, " [t=%.0fs]", phase_clock_()); return b; }

  /// [perf-instr] Per-phase wall/CPU/RSS accounting. Sparse [t=..] markers left ~53% of wall time
  /// unattributed, which is how you optimise the wrong thing. Every phase is timed; the table is
  /// printed unconditionally at the end so a benchmark log always says where the time went.
  struct PhaseStat { double wall = 0; double cpu = 0; long rss_end_mb = 0; int n = 0; };
  inline double& flush_wall_() { static double v = 0; return v; }   // [perf-load] pick inside LOAD
  inline std::map<std::string, PhaseStat>& phase_stats_()
  { static std::map<std::string, PhaseStat> m; return m; }
  inline std::vector<std::string>& phase_order_()
  { static std::vector<std::string> v; return v; }
  static double cpu_seconds_()
  {
    struct rusage ru;                       // RUSAGE_SELF sums ALL threads -> parallel efficiency
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
    return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6
         + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
  }
  static long rss_mb_()
  {
    std::ifstream f("/proc/self/status"); std::string k;   // deliberate: Linux-only, returns 0 elsewhere
    while (f >> k) { if (k == "VmRSS:") { long v; f >> v; return v / 1024; } f.ignore(1 << 20, '\n'); }
    return 0;
  }
  struct Phase
  {
    std::string name; double w0, c0;
    explicit Phase(std::string n) : name(std::move(n)), w0(phase_clock_()), c0(cpu_seconds_()) {}
    ~Phase()
    {
      auto& st = phase_stats_()[name];
      if (st.n == 0) phase_order_().push_back(name);
      st.wall += phase_clock_() - w0; st.cpu += cpu_seconds_() - c0;
      st.rss_end_mb = rss_mb_(); ++st.n;
    }
  };
  static void report_phases_(double total_wall)
  {
    if (phase_order_().empty()) return;
    OPENMS_LOG_INFO << "\n[perf] phase breakdown (wall s | % of total | CPU s | par x | RSS MB at end)\n";
    double acc = 0;
    for (const auto& n : phase_order_())
    {
      const auto& p = phase_stats_()[n];
      acc += p.wall;
      char b[220];
      snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  %9.1f  %5.1fx  %8ld%s\n", n.c_str(), p.wall,
               total_wall > 0 ? 100.0 * p.wall / total_wall : 0.0, p.cpu,
               p.wall > 0 ? p.cpu / p.wall : 0.0, p.rss_end_mb, p.n > 1 ? "  (summed)" : "");
      OPENMS_LOG_INFO << b;
    }
    char b[160];
    snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  (unattributed %.1f s)\n", "TOTAL(measured)", acc,
             total_wall > 0 ? 100.0 * acc / total_wall : 0.0, total_wall - acc);
    OPENMS_LOG_INFO << b;
  }

  bool log_overlap_ = false;   ///< [B0] gate:coelution == "logoverlap"
  /// [wavelet] Precursor-XIC smoothing scale in SECONDS, derived from the MEASURED MS1 FWHM
  /// (trace:wavelet_smooth x median FWHM). 0 = off. Set once, after MS1 tracing, before the
  /// window loop -- it cannot be read from a parameter directly because the FWHM is a measurement.
  double wl_scale_ = 0.0;
  bool rank_by_intensity_ = false;   ///< [rank] assembly:rank_by == "intensity" (the reference implementation-style)
  double im_weight_sigma_ = 0.0;     ///< [recall] Gaussian IM-proximity weight on emitted fragment intensity (0=off)
  double mono_guard_ = 0.0;          ///< [mono-guard] averagine slack on leftward isotope steps (0=off)
  bool mono_select_ = false;         ///< [mono-select] pick the mono by averagine fit over the isotope run
  std::atomic<long long> stat_slice_{0};    ///< fragments VISITED by the IM slice
  std::atomic<long long> stat_rtpass_{0};   ///< of those, how many pass the RT gate
  std::atomic<long long> stat_prec_{0};     ///< precursors scored
  bool var_support_ = false;         ///< [Q1] Pearson variance over union support, not full grid G (0-pad fix)
  double corr_power_ = 0.0;          ///< [Q2] emitted fragment intensity *= corr^corr_power (0=off; engine-agnostic)
  bool drop_prec_iso_ = std::getenv("SPEXTRACTOR_DROP_PREC_ISO") != nullptr;   ///< [C9 A/B] drop precursor M+1..M+3 from fragment lists

  /// glibc arena accounting at a milestone: what the allocator holds in arenas vs mmapped blocks
  /// and how much of it is free-but-retained. RSS alone cannot tell live data from retained pages.
  static std::string mem_()
  {
#ifdef __GLIBC__
    const struct mallinfo2 mi = mallinfo2();
    return " [malloc arena=" + std::to_string(mi.arena >> 20) + " MB mmap=" + std::to_string(mi.hblkhd >> 20)
         + " MB free=" + std::to_string(mi.fordblks >> 20) + " MB]";
#else
    return "";
#endif
  }
  static std::string rss_()
  {
    std::ifstream st("/proc/self/status");
    if (!st) return "";
    std::string line, cur = "?", peak = "?";
    while (std::getline(st, line))
    {
      auto val = [&line]() {
        std::string v = line.substr(line.find(':') + 1);
        size_t a = v.find_first_not_of(" \t");
        return a == std::string::npos ? std::string("?") : v.substr(a);
      };
      if (line.rfind("VmRSS:", 0) == 0) cur = val();
      else if (line.rfind("VmHWM:", 0) == 0) peak = val();
    }
    return " [RSS " + cur + ", peak " + peak + "]";
  }

  /// [dyn-mem] Bytes we may still allocate. MemAvailable is the kernel's own estimate of what
  /// can be handed out without swapping (it accounts for reclaimable page cache), which is the
  /// right quantity on a SHARED node -- MemFree would ignore cache and MemTotal would ignore
  /// the other users. Add our own RSS back: memory we already hold is part of our budget, not
  /// somebody else's. Returns 0 if unreadable, which the caller treats as "admit one at a time".
  static size_t availableBytes_()
  {
    std::ifstream mi("/proc/meminfo");
    if (!mi) return 0;
    std::string line;
    size_t avail_kb = 0;
    while (std::getline(mi, line))
      if (line.rfind("MemAvailable:", 0) == 0)
      { avail_kb = strtoull(line.c_str() + 13, nullptr, 10); break; }
    if (!avail_kb) return 0;
    std::ifstream st("/proc/self/status");
    size_t rss_kb = 0;
    while (st && std::getline(st, line))
      if (line.rfind("VmRSS:", 0) == 0)
      { rss_kb = strtoull(line.c_str() + 6, nullptr, 10); break; }
    return (avail_kb + rss_kb) * 1024ull;
  }

  /// THE RT AXIS IS GLOBAL AND DISCRETE. Every acquired frame has one retention time, and every
  /// point of every mass trace sits on one of those frames -- so an RT value never needs to be
  /// stored per point. It used to be: `vector<pair<double,double>>` cost 16 B per point and
  /// duplicated the same ~33,000 values across hundreds of millions of points, and aligning two
  /// profiles meant a binary search per point because the values had to be compared numerically.
  /// Here the axis is stored ONCE, a point is (frame index, intensity) at 8 B, and alignment is
  /// integer equality. Built once after loading, read-only thereafter.
  inline vector<double>& rtAxis() { static vector<double> v; return v; }

  /// Frame index of an RT. Used only while BUILDING traces: the value comes from a frame, so the
  /// match is exact, and the nearest-neighbour fallback only guards a last-ulp difference.
  inline uint32_t rtIndex(double rt)
  {
    const vector<double>& a = rtAxis();
    if (a.empty()) return 0;
    auto it = lower_bound(a.begin(), a.end(), rt);
    if (it == a.end()) return (uint32_t)(a.size() - 1);
    if (it != a.begin() && (*it - rt) > (rt - *(it - 1))) --it;
    return (uint32_t)(it - a.begin());
  }

  struct Trace
  {
    double mz = 0.0;
    double rt = 0.0;   ///< apex/centroid RT
    double im = 0.0;   ///< centroid ion mobility (1/K0)
    double intensity = 0.0;
    /// The flight-time bin this trace was measured on, and the calibration factor of its apex
    /// frame. The pair is the AUTHORITATIVE m/z: `mz` below is a cached convenience for the gates
    /// that still compare in m/z space, and the value actually written out is recomputed from these
    /// two at export. 0 means the trace came from the OpenMS detector, which has no bin.
    uint32_t tof = 0;
    double   b = 0.0;
    /// THE PROFILE IS A SPAN, NOT A POINT LIST. A trace is contiguous in its window's frame
    /// sequence (with gaps the detector tolerates), so it needs an entry frame and its
    /// intensities -- not a frame per point and not a bin per point, which is what the three
    /// parallel vectors this replaces stored (12 B/point + 3 mallocs/trace, measured 35.6 GB).
    /// The intensities live in the owning TraceStore's arena at [off, off+len); ZERO MEANS
    /// MISSING, which is valid because every observed intensity is > noise_threshold >= 0.
    /// Every consumer walks REAL points in frame order (skipping zeros), so the floating-point
    /// sums are visited in the same order as before: the change is pure representation.
    /// `frame0` is LOCAL to the store's frame table (a window's, or the MS1 map's) -- a
    /// window's consecutive frames are ~13 apart on the global RT axis, so a global entry
    /// frame would not be contiguous. See docs/TRACE-STRUCT-REDESIGN.md, revision 2.
    const struct TraceStore* st = nullptr;   ///< set once the store is immutable (after merges)
    uint32_t frame0 = 0;   ///< first frame of the span, store-local
    uint32_t off = 0;      ///< offset of the span in the store's arena
    uint16_t len = 0;      ///< frames spanned (a gradient is ~1,400; asserted < 65536)
    uint16_t npts = 0;     ///< REAL points in the span
    uint16_t apex = 0;     ///< position of the apex within the span
    size_t np() const { return npts; }
    size_t span() const { return len; }
    float  xv(size_t k) const;                       ///< intensity at span position k (0 = missing)
    bool   real(size_t k) const { return xv(k) > 0.0f; }
    uint32_t gframe(size_t k) const;                 ///< global RT-axis index of span position k
    double rtAtSpan(size_t k) const { return rtAxis()[gframe(k)]; }
    void freeProfile() { len = 0; npts = 0; }        ///< the arena is the store's; nothing to free
  };

  /// The frame table and the arenas that a set of traces point into: one per window, one for the
  /// MS1 map. `rt_index`/`b` are per FRAME (store-local index); `inten` holds every trace's span
  /// back to back (zero = missing); `bins` holds the per-frame flight-time bin for the integer
  /// path only and is freed once valley splitting has assigned each child its apex bin.
  struct TraceStore
  {
    vector<uint32_t> rt_index;   ///< frame -> global RT axis
    vector<double>   b;          ///< frame -> calibration factor (empty on the OpenMS path)
    vector<double>   frame_rt;   ///< frame -> RT value, for toTrace()'s frame lookup
    vector<float>    inten;      ///< arena of spans
    vector<uint32_t> bins;       ///< arena of per-frame bins (integer path, until EPD)
    size_t frames() const { return rt_index.size(); }
    void setFrames(const vector<uint32_t>& ri, const vector<double>* bb)
    {
      rt_index = ri; if (bb) b = *bb; else b.clear();
      frame_rt.resize(ri.size());
      for (size_t f = 0; f < ri.size(); ++f) frame_rt[f] = rtAxis()[ri[f]];
    }
    /// Store-local frame of an RT value (exact by construction; nearest guards a last-ulp slip).
    uint32_t frameOf(double rt) const
    {
      auto it = lower_bound(frame_rt.begin(), frame_rt.end(), rt);
      if (it == frame_rt.end()) return (uint32_t)(frame_rt.size() - 1);
      if (it != frame_rt.begin() && (*it - rt) > (rt - *(it - 1))) --it;
      return (uint32_t)(it - frame_rt.begin());
    }
    /// Append `other`'s arenas and return the offset they now start at (for rebasing).
    uint32_t absorb(TraceStore& other)
    {
      const size_t base = inten.size();
      if (base + other.inten.size() >= (size_t)std::numeric_limits<uint32_t>::max())
        throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      inten.insert(inten.end(), other.inten.begin(), other.inten.end());
      if (!other.bins.empty()) { bins.resize(base, 0u); bins.insert(bins.end(), other.bins.begin(), other.bins.end()); }
      // Once ANY source has contributed bins, keep bins index-parallel to inten: a later binless
      // source would otherwise leave bins short, and a bins[t.off + k] read would run off the end.
      // No current caller mixes the two, which is exactly why this is worth pinning down here.
      else if (!bins.empty()) bins.resize(inten.size(), 0u);
      vector<float>().swap(other.inten); vector<uint32_t>().swap(other.bins);
      return (uint32_t)base;
    }
  };
  inline float    Trace::xv(size_t k) const { return st->inten[off + k]; }
  inline uint32_t Trace::gframe(size_t k) const { return st->rt_index[frame0 + k]; }

  /// The real points of a trace, packed in frame order, for the consumers that were written over
  /// point lists (smoother, FWHM): identical sequences, so identical arithmetic.
  inline void packReal(const Trace& t, vector<double>& rt, vector<double>& v)
  {
    rt.clear(); v.clear();
    for (size_t k = 0; k < t.span(); ++k)
      if (t.real(k)) { rt.push_back(t.rtAtSpan(k)); v.push_back((double)t.xv(k)); }
  }

  /// TRIM a trace to at most `cap` seconds around its apex, AFTER detection and valley splitting
  /// have decided what the peak is. This is the safe form of a length limit:
  ///   * detection is untouched, so the greedy consume-and-mark ownership of peaks is unchanged --
  ///     unlike terminating extension, which changes which peaks belong to which trace;
  ///   * it works on MS1 as well as MS2, because it is a step on Trace objects rather than
  ///     something inside OpenMS's extension loop (MS1 is always traced by the OpenMS path, and
  ///     MS1 is where the long traces are: 15.2% of dataset D MS1 traces span >45 s, the longest 1,244 s
  ///     of a 1,860 s gradient, against a 375 s maximum for fragments);
  ///   * what it removes is the far tail of a blob that ElutionPeakDetection did not split, not
  ///     part of a chromatographic peak -- real peaks here are 5-30 s wide.
  /// On the span representation this is free: three integers move, no bytes are copied. The arena
  /// still holds the untrimmed bytes; reclaiming them needs a compaction pass (MS1 already has one).
  inline void trimToSpan(Trace& t, double cap)
  {
    if (cap <= 0.0 || t.span() == 0) return;
    const size_t a = t.apex;
    size_t lo = a, hi = a;
    // grow symmetrically outward from the apex while the span fits, preferring the earlier frame on
    // a tie so the result does not depend on which side is tested first
    for (;;)
    {
      const bool can_lo = lo > 0, can_hi = hi + 1 < t.span();
      if (!can_lo && !can_hi) break;
      const double d_lo = can_lo ? t.rtAtSpan(hi) - t.rtAtSpan(lo - 1) : 1e30;
      const double d_hi = can_hi ? t.rtAtSpan(hi + 1) - t.rtAtSpan(lo) : 1e30;
      if (d_lo <= d_hi) { if (!can_lo || d_lo > cap) { if (!can_hi || d_hi > cap) break; ++hi; } else --lo; }
      else              { if (!can_hi || d_hi > cap) { if (!can_lo || d_lo > cap) break; --lo; } else ++hi; }
    }
    // Zero-trim the ends. makeSpan starts and ends a span on a real point by construction, but
    // growing outward from the apex can stop on a gap, and a leading or trailing zero is pure
    // padding: it costs 4 B, it is skipped by every consumer, and it makes `span()` overstate the
    // trace. Interior zeros are the gaps and MUST stay -- they are the presence encoding.
    while (lo < hi && t.xv(lo) == 0.0f) ++lo;
    while (hi > lo && t.xv(hi) == 0.0f) --hi;
    if (lo == 0 && hi + 1 == t.span()) return;                  // already within the cap, ends real
    uint16_t np = 0;
    for (size_t k = lo; k <= hi; ++k) if (t.xv(k) > 0.0f) ++np;
    t.frame0 += (uint32_t)lo; t.off += (uint32_t)lo;
    t.len = (uint16_t)(hi - lo + 1); t.apex = (uint16_t)(a - lo); t.npts = np;
  }

  /// Drop the spans of the traces nobody references and compact the arena IN PLACE, walking the kept
  /// spans in ORIGINAL OFFSET order -- not container order. ms1_traces is sorted by m/z after its
  /// spans were appended in detection order, so container order is not offset order, and the
  /// previous walk (monotone write cursor over the container) overwrote spans it had not copied yet:
  /// ~1% of precursor XICs, deterministically. Spans are pairwise disjoint on this path (every
  /// child gets its own makeSpan append; trimToSpan only shrinks), and that is VALIDATED before a
  /// byte moves: in ascending-offset order the write cursor can never pass the next source.
  /// A freed trace keeps its scalars and loses its span. Returns the number of profiles released.
  inline Size compactUnreferenced(vector<Trace>& tr, TraceStore& store, const vector<bool>& needed)
  {
    vector<uint32_t> ids;
    Size freed = 0;
    for (size_t i = 0; i < tr.size(); ++i)
    {
      if (!needed[i]) { if (tr[i].np() > 0) ++freed; tr[i].freeProfile(); continue; }
      if (tr[i].len > 0) ids.push_back((uint32_t)i);     // a zero-length span owns no bytes
    }
    sort(ids.begin(), ids.end(), [&](uint32_t a, uint32_t b) { return tr[a].off < tr[b].off; });
    uint64_t prev_end = 0;
    for (uint32_t i : ids)
    {
      const uint64_t beg = tr[i].off, end = beg + (uint64_t)tr[i].len;
      if (end > store.inten.size() || beg < prev_end)
        throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "MS1 spans overlap or exceed the arena");
      prev_end = end;
    }
    size_t w = 0;
    for (uint32_t i : ids)
    {
      Trace& t = tr[i];
      if (w != t.off) std::memmove(store.inten.data() + w, store.inten.data() + t.off, (size_t)t.len * sizeof(float));
      t.off = (uint32_t)w; w += t.len;
    }
    store.inten.resize(w); store.inten.shrink_to_fit();
    return freed;
  }

  /// Build a span from (store-local frame, intensity) points into `store`, returning the trace
  /// with `st` UNSET: the caller sets it once the store is immutable. Points may arrive unsorted.
  inline Trace makeSpan(vector<pair<uint32_t, float>>& pts, TraceStore& store)
  {
    Trace t;
    if (pts.empty()) return t;
    sort(pts.begin(), pts.end());
    const uint32_t f0 = pts.front().first, f1 = pts.back().first;
    const size_t L = (size_t)f1 - f0 + 1;
    if (L >= 65536) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    t.frame0 = f0; t.len = (uint16_t)L; t.off = (uint32_t)store.inten.size();
    store.inten.resize(store.inten.size() + L, 0.0f);
    float best = -1.0f; uint16_t np = 0;
    for (const auto& q : pts)
    {
      const size_t k = q.first - f0;
      if (store.inten[t.off + k] == 0.0f) ++np;      // one point per frame; a repeat overwrites
      store.inten[t.off + k] = q.second;
      if (q.second > best) { best = q.second; t.apex = (uint16_t)k; }
    }
    t.npts = np;
    return t;
  }

  /// A precursor hypothesis after isotope/charge inference.
  struct Precursor_
  {
    double mono_mz = 0.0;
    uint32_t mono_tof = 0;   ///< the monoisotope's flight-time bin, and its frame's factor: the
    double   mono_b = 0.0;   ///  reported m/z is recomputed from these when the spectrum is written
    int charge = 0;          ///< 0 == unknown
    double rt = 0.0;
    double im = 0.0;
    size_t trace_idx = 0;    ///< index into ms1 traces
    bool guessed = false;    ///< [way-4] charge was DEFAULTED, not isotope-supported
    int n_isotopes = 0;      ///< [Q4] isotope peaks behind the mono call (0/1 = weak, mono ±1.00335 uncertain)
  };

  /// MassTraceDetection locates the per-peak IM array by the exact name
  /// Constants::UserParam::ION_MOBILITY ("Ion Mobility"). Converters name it via CV term and
  /// PeakPickerIM names it "Ion Mobility Centroid"; rename whichever exists so tracing is IM-aware.
  void ensureIMArrayName(MSSpectrum& s)
  {
    auto& fdas = s.getFloatDataArrays();
    if (s.containsIMData())
    {
      fdas[s.getIMData().first].setName(Constants::UserParam::ION_MOBILITY);
      return;
    }
    for (auto& fda : fdas)
    {
      if (fda.getName() == Constants::UserParam::ION_MOBILITY_CENTROID)
      {
        fda.setName(Constants::UserParam::ION_MOBILITY);
        return;
      }
    }
  }

  /// Reported-m/z estimator, set once from `trace:mz_estimator` before any trace is converted.
  /// A file-scope value because toTrace() is a free helper with no access to the tool's parameters.
  inline std::string& mzEstimator() { static std::string v = "apex"; return v; }

  Trace toTrace(const MassTrace& mt, TraceStore& store)
  {
    Trace t;
    // Reported m/z of a mass trace. "apex" (the default) takes the m/z of the most intense peak;
    // "mean" is the OpenMS intensity-weighted centroid; "median" the median m/z. Set by
    // `trace:mz_estimator` -- see docs/dataset D-BASELINE.md for the comparison that chose apex.
    const std::string& mz_est = mzEstimator();
    if (mz_est == "apex" && mt.getSize() > 0) t.mz = mt[mt.findMaxByIntPeak(false)].getMZ();
    else if (mz_est == "median" && mt.getSize() > 0) { MassTrace m2 = mt; m2.updateMedianMZ(); t.mz = m2.getCentroidMZ(); }
    else t.mz = mt.getCentroidMZ();
    t.rt = mt.getCentroidRT();
    t.im = mt.getCentroidIM();
    t.intensity = mt.getMaxIntensity(false);
    // One peak per frame (MassTraceDetection takes at most one candidate per spectrum), so the
    // profile is a span over the store's frame table with zeros where frames were missed.
    const Size np = mt.getSize();
    vector<pair<uint32_t, float>> pts;
    pts.reserve(np);
    for (Size i = 0; i < np; ++i) pts.emplace_back(store.frameOf(mt[i].getRT()), (float)mt[i].getIntensity());
    Trace sp = makeSpan(pts, store);
    sp.mz = t.mz; sp.rt = t.rt; sp.im = t.im; sp.intensity = t.intensity;
    return sp;
  }

  /// Per-window correlation grid: all fragment XICs are resampled onto ONE common RT axis
  /// (the sorted union of fragment-trace RTs) so a precursor can be correlated against every
  /// candidate fragment as a precomputed dot product instead of an O(P log F) per-pair search.
  /// This is the hot-path optimization: precursor XIC is prepared once per precursor (not per
  /// fragment), fragment means/norms are precomputed once per window. [perf]
  /// CSR layout: one `vector<vector<...>>` costs 24 B of header PER FRAGMENT before a single


  /// [cross-frame] Sliding-window aggregation of adjacent RT frames within ONE isolation window.
  ///
  /// the reference implementation sums "adjacent neighbor RT frames" of the same isolation window into a composite
  /// (m/z x 1/K0) matrix BEFORE peak picking, which raises per-point S/N by ~sqrt(N) and lets weak
  /// fragments form traces at all (Nat Commun 16:95 Methods). Because a diaPASEF isolation window is
  /// sampled exactly ONCE per cycle, "adjacent frames" == adjacent CYCLES.
  ///
  /// SLIDING (stride 1), not block: block summing would divide our RT point count by N, and peaks here
  /// span only 4-22 points, which would leave too few points to correlate at all. The cost of sliding
  /// is that adjacent output points share input frames and are therefore AUTOCORRELATED, which inflates
  /// Pearson r for signal and noise alike - the reason this must be benchmarked, not assumed.
  ///
  /// Output spectrum i = intensity-weighted (m/z, IM) clustering of the peaks of spectra
  /// [i-N/2, i+N/2], carrying the ORIGINAL RT of spectrum i so the RT axis is unchanged.
  void aggregateFrames_(PeakMap& wmap, int n_frames, double ppm, double im_tol)
  {
    if (n_frames <= 1 || wmap.empty()) return;
    const int half = n_frames / 2;
    const Size N = wmap.size();

    struct Pk { double mz; double im; double in; };
    PeakMap out;
    out.reserve(N);

    vector<Pk> buf;
    for (Size i = 0; i < N; ++i)
    {
      const Size lo = (i > (Size)half) ? i - half : 0;
      const Size hi = std::min(N - 1, i + (Size)half);

      buf.clear();
      for (Size k = lo; k <= hi; ++k)
      {
        const MSSpectrum& sp = wmap[k];
        if (!sp.containsIMData()) continue;
        const auto& imarr = sp.getFloatDataArrays()[sp.getIMData().first];
        for (Size p = 0; p < sp.size(); ++p) buf.push_back({sp[p].getMZ(), (double)imarr[p], sp[p].getIntensity()});
      }
      sort(buf.begin(), buf.end(), [](const Pk& a, const Pk& b) {
        if (a.mz != b.mz) return a.mz < b.mz;
        return a.im < b.im;
      });

      MSSpectrum agg;
      agg.setRT(wmap[i].getRT());                 // RT axis preserved (sliding, stride 1)
      agg.setMSLevel(wmap[i].getMSLevel());
      agg.setPrecursors(wmap[i].getPrecursors());
      OpenMS::DataArrays::FloatDataArray im_out;
      im_out.setName(Constants::UserParam::ION_MOBILITY);   // keep tracing IM-aware

      // greedy cluster in m/z (ppm) then IM: same ion across adjacent cycles collapses, intensities SUM
      for (size_t a = 0; a < buf.size(); )
      {
        const double mz0 = buf[a].mz;
        const double tol = mz0 * ppm * 1e-6;
        size_t b = a;
        double sum = 0.0, wmz = 0.0, wim = 0.0;
        while (b < buf.size() && buf[b].mz - mz0 <= tol)
        {
          if (fabs(buf[b].im - buf[a].im) <= im_tol)
          {
            sum += buf[b].in; wmz += buf[b].mz * buf[b].in; wim += buf[b].im * buf[b].in;
            buf[b].in = -1.0;                    // consumed
          }
          ++b;
        }
        if (sum > 0.0)
        {
          Peak1D pk; pk.setMZ(wmz / sum); pk.setIntensity((float)sum);
          agg.push_back(pk);
          im_out.push_back((float)(wim / sum));
        }
        // advance to next unconsumed peak
        while (a < buf.size() && buf[a].in < 0.0) ++a;
        if (a >= buf.size()) break;
      }
      agg.getFloatDataArrays().push_back(im_out);
      agg.sortByPosition();
      out.addSpectrum(std::move(agg));
    }
    wmap = std::move(out);
  }

  /// [compact] Quantised frame store for the per-window buffer.
  ///
  /// Holding all 24 isolation windows as full OpenMS PeakMaps costs ~89 GB resident before any window
  /// is consumed, which is what forces `max_concurrent_windows` down and caps our parallelism at ~7
  /// cores (the reference implementation uses ~43). A frame only needs (m/z, intensity, 1/K0) per peak plus one RT, so we
  /// store structure-of-arrays with quantised keys and materialise a PeakMap ONLY inside the worker
  /// that is about to trace it, freeing the compact copy immediately.
  ///
  /// Precision (deliberately far finer than every downstream tolerance, so tracing is unaffected):
  ///   m/z : uint32 at 1e-5 Da  -> 0.007 ppm at m/z 1400, vs a 20 ppm trace tolerance (2800x margin)
  ///   1/K0: uint16 over [0.4,1.8] -> 2.1e-5, vs a 0.01 IM tolerance (470x margin)
  ///   intensity: float32, unchanged (Peak1D already stores float)
  /// Cost: 10 B/peak (4+4+2) vs ~20 B/peak for Peak1D+FloatDataArray, plus no per-MSSpectrum overhead.
  constexpr double MZ_Q = 1e5;      ///< m/z quantum = 1e-5 Da

  /// THE M/Z AXIS IS THE INSTRUMENT'S FLIGHT-TIME INDEX. A timsTOF measures an integer TOF bin;
  /// m/z is a calibrated function of it, not a measured quantity. The compact store had been
  /// quantising the CALIBRATED m/z onto a uniform 1e-5 Da grid -- a second discretisation on top of
  /// the instrument's own, and a coarse approximation of the wrong axis. This tool only runs on
  /// ion-mobility TOF data, so the native axis is always available and there is no reason to
  /// approximate it.
  ///
  /// Why the native axis is better than a uniform m/z grid, not merely smaller:
  ///   * it is LOSSLESS -- the bin IS the measurement, so no quantum has to be justified;
  ///   * ~5e5 bins instead of ~1.7e8 grid points, so a bin index is a usable array subscript;
  ///   * the SAME ion gives the SAME bin in every frame, while its calibrated m/z wobbles with the
  ///     digitizer temperature -- so grouping in TOF space does not have to spend part of its ppm
  ///     tolerance absorbing calibration drift, which is exactly what the m/z-space tolerance does;
  ///   * a tolerance in m/z ppm becomes a handful of integer bins (m ~ tof^2, so a relative m/z
  ///     tolerance is about half of it in tof), which makes candidate lookup pure arithmetic.
  ///
  /// The per-frame factor `b` (from the digitizer temperature) is what turns a bin into an m/z, and
  /// it is applied at REPORT time, where it belongs -- the extraction never needs a calibrated m/z.
  using TofIdx = uint32_t;

  /// The flight-time axis for one run: the calibration model plus each frame's temperature factor.
  /// Fail-closed: without a calibration there is no TOF axis and the caller must not pretend there
  /// is one.
  struct TofAxis
  {
    spextractor::TdfMzCalibration cal;
    vector<double> b_by_frame;        ///< frame factor per FRAME ID (index 0 unused, as in the tdf)
    bool ok = false;

    double factor(size_t frame_id) const
    { return frame_id < b_by_frame.size() ? b_by_frame[frame_id] : cal.frameFactor(cal.T1_ref); }
    /// frame 0 is the unknown-frame sentinel; nothing may be calibrated against it
    bool mapped(uint32_t frame_id) const { return frame_id != 0 && factor(frame_id) > 0.0; }
    double mzOf(TofIdx tof, double b) const { return cal.tofToMz((double)tof, b); }
    TofIdx tofOf(double mz, double b) const
    {
      const double t = cal.mzToTof(mz, b);
      return t > 0.0 ? (TofIdx)llround(t) : 0u;
    }
    /// A relative m/z tolerance as a bin count, derived rather than assumed: convert both ends.
    TofIdx span(TofIdx tof, double b, double ppm) const
    {
      const double mz = mzOf(tof, b);
      if (!(mz > 0.0)) return 1u;
      const double t2 = cal.mzToTof(mz * (1.0 + ppm * 1e-6), b);
      const double d = fabs(t2 - (double)tof);
      return (TofIdx)std::max(1.0, d + 0.5);
    }
  };

  inline TofAxis& tofAxis() { static TofAxis a; return a; }

  /// Populate the flight-time axis from a Bruker analysis.tdf (the .d directory, or a sidecar).
  /// Returns false and leaves the axis unusable if anything is missing -- the detector that needs it
  /// refuses to run rather than substituting a reference temperature for every frame.
  inline bool loadTofAxis(const std::string& tdf_path, String& why)
  {
    TofAxis& ax = tofAxis();
    ax.ok = false;
    std::vector<double> t1;
    if (!spextractor::loadTdfCalibration(tdf_path, ax.cal, t1, why)) return false;
    ax.b_by_frame.assign(t1.size(), 0.0);
    for (size_t i = 1; i < t1.size(); ++i) ax.b_by_frame[i] = ax.cal.frameFactor(t1[i]);
    // Frames the tdf did not list keep the reference factor rather than a zero -- but index 0
    // is NOT such a frame. It is the "unknown frame" sentinel (the tdf's Frames.Id starts at 1),
    // and frameIdOf() returns 0 for any nativeID it cannot parse. Backfilling it to bref is what
    // let an mzML+sidecar input skip the fallback warning and calibrate every unmapped frame at
    // reference temperature -- the invented calibration this path is supposed to refuse.
    const double bref = ax.cal.frameFactor(ax.cal.T1_ref);
    for (size_t i = 1; i < ax.b_by_frame.size(); ++i) if (!(ax.b_by_frame[i] > 0.0)) ax.b_by_frame[i] = bref;
    ax.b_by_frame[0] = 0.0;
    ax.ok = true;
    return true;
  }
  constexpr double IM_LO = 0.4, IM_HI = 1.8;
  constexpr double IM_Q = 65535.0 / (IM_HI - IM_LO);

  struct CompactFrame
  {
    double rt = 0.0;
    uint32_t frame_id = 0;   ///< vendor frame Id: the key the per-frame calibration factor is on
    vector<uint32_t> mzq;
    vector<float> inten;
    vector<uint16_t> imq;
    size_t bytes() const { return mzq.size() * 4 + inten.size() * 4 + imq.size() * 2; }
  };

  /// Peaks that cannot be represented exactly are DROPPED and COUNTED, never silently altered.
  /// [adv-fix] The first version CLAMPED out-of-range 1/K0 to IM_LO/IM_HI, which would pin outliers
  /// from many frames onto two artificial, perfectly-aligned mobility values — manufacturing spurious
  /// traces out of nothing. NaN also passed straight through clamp into llround (UB on unsigned
  /// conversion). Dropping + counting is the only defensible behaviour.
  struct CompactStats
  {
    size_t no_im_array = 0, size_mismatch = 0, bad_mz = 0, bad_im = 0, kept = 0;
    size_t unmapped_frame = 0;   ///< nativeID carried no parseable "frame=": cannot be calibrated
    /// merge per-frame stats after a parallel batch [perf-load]
    CompactStats& operator+=(const CompactStats& o)
    {
      no_im_array += o.no_im_array; size_mismatch += o.size_mismatch;
      bad_mz += o.bad_mz; bad_im += o.bad_im; kept += o.kept;
      unmapped_frame += o.unmapped_frame;
      return *this;
    }
  };

  /// Convert one picked MSSpectrum into the compact form (drops everything tracing does not read).
  /// Vendor frame Id from a native identifier. Bruker writes "frame=<N> scan=<M>"; the mzPeak
  /// reader writes "mzpeak=<frame index> window=<w>". Anything else yields 0, which makes the
  /// calibration fall back to the reference temperature -- so the caller must treat 0 as unknown.
  inline uint32_t nidField(const String& nid, const char* key)
  {
    const size_t i = nid.find(key);
    if (i == std::string::npos) return 0;
    uint32_t v = 0; size_t k = i + strlen(key);
    for (; k < nid.size() && isdigit((unsigned char)nid[k]); ++k) v = v * 10 + (uint32_t)(nid[k] - '0');
    return v;
  }
  /// Only the vendor's own "frame=" key names a Frames.Id. The mzPeak reader's "mzpeak=" is an
  /// ARCHIVE index and must not be used as one; that input falls back to the reference factor.
  inline uint32_t frameIdOf(const String& nid) { return nidField(nid, "frame="); }
  /// "windowGroup=<WindowGroup>": which acquisition group this frame belongs to. A group holds each
  /// isolation m/z at most once, so (m/z, group) names one ion-mobility slice of the window; the
  /// scan range alone does not (two slices can share ScanNumBegin = 0 and differ in ScanNumEnd).
  /// 0 when absent (mzPeak input), which keeps that path keyed by m/z alone, as before.
  inline uint32_t windowGroupOf(const String& nid) { return nidField(nid, "windowGroup="); }

  CompactFrame compactify(const MSSpectrum& s, CompactStats& st)
  {
    CompactFrame f;
    f.rt = s.getRT();
    f.frame_id = frameIdOf(s.getNativeID());
    if (f.frame_id == 0) ++st.unmapped_frame;
    if (!s.containsIMData()) { st.no_im_array += s.size(); return f; }
    const auto& im = s.getFloatDataArrays()[s.getIMData().first];
    const Size n = s.size();
    if (im.size() != n) { st.size_mismatch += n; return f; }   // [adv-fix] no out-of-bounds read
    f.mzq.reserve(n); f.inten.reserve(n); f.imq.reserve(n);
    for (Size i = 0; i < n; ++i)
    {
      const double mz = s[i].getMZ();
      const double imv = (double)im[i];
      // [adv-fix] validate BEFORE the unsigned conversion (NaN/inf/negative -> UB otherwise)
      if (!std::isfinite(mz) || mz <= 0.0 || mz * MZ_Q > 4.2e9) { ++st.bad_mz; continue; }
      if (!std::isfinite(imv) || imv < IM_LO || imv > IM_HI)    { ++st.bad_im; continue; }
      f.mzq.push_back((uint32_t)llround(mz * MZ_Q));
      f.inten.push_back(s[i].getIntensity());
      f.imq.push_back((uint16_t)llround((imv - IM_LO) * IM_Q));
      ++st.kept;
    }
    return f;
  }


  /// [det] Order-insensitive 64-bit digest of a stage's numeric output, for localising
  /// run-to-run nondeterminism: sort the bit patterns, then FNV-1a. Logged under SPEXTRACTOR_DET=1.
  static uint64_t detDigest_(vector<uint64_t> bits)
  {
    std::sort(bits.begin(), bits.end());
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t b : bits) { h ^= b; h *= 1099511628211ULL; }
    return h;
  }
  static bool detOn_() { static const bool on = std::getenv("SPEXTRACTOR_DET") != nullptr; return on; }
  static uint64_t traceDigest_(const vector<Trace>& ts)
  {
    vector<uint64_t> b; b.reserve(ts.size() * 3);
    for (const auto& t : ts)
    {
      uint64_t u; double v;
      v = t.mz; std::memcpy(&u, &v, 8); b.push_back(u);
      v = t.rt; std::memcpy(&u, &v, 8); b.push_back(u);
      v = t.intensity; std::memcpy(&u, &v, 8); b.push_back(u);
    }
    return detDigest_(std::move(b));
  }

  /// One window's peaks in INDEX SPACE, frame-major. This is what the compact store already holds;
  /// the point is that it is never converted back to Peak1D doubles. Replaces the materialised
  /// PeakMap for the integer detector. See docs/MZ-AXIS-DESIGN.md.
  struct PeakSlab
  {
    vector<uint32_t> frame_off;    ///< size F+1: frame f owns [frame_off[f], frame_off[f+1])
    vector<uint32_t> rt_index;     ///< size F: the frame's index on the GLOBAL RT axis
    vector<TofIdx>   tof;          ///< flight-time bin, ascending within a frame
    vector<double>   b;            ///< size F: the frame's calibration factor (bin -> m/z)
    vector<float>    inten;
    vector<uint16_t> imq;          ///< the store's quantised 1/K0, kept for export
    size_t frames() const { return rt_index.size(); }
    size_t peaks()  const { return tof.size(); }
  };

  /// Ordered, incremental digest of a window's slab: frame table (rt index, factor bits, count) then
  /// every (tof, intensity bits, imq) tuple in frame order. Unlike detDigest_ it is NOT sorted, so a
  /// reordered equal-intensity peak or a tuple swapped between frames changes it -- which is what a
  /// representation change of the store must not do. [det]
  static uint64_t slabDigest_(const PeakSlab& sl)
  {
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](uint64_t v) { h ^= v; h *= 1099511628211ULL; };
    for (size_t f = 0; f < sl.frames(); ++f)
    {
      uint64_t bb; std::memcpy(&bb, &sl.b[f], 8);
      mix(sl.rt_index[f]); mix(bb); mix(sl.frame_off[f + 1] - sl.frame_off[f]);
      for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
      { uint32_t ib; std::memcpy(&ib, &sl.inten[k], 4); mix(sl.tof[k]); mix(ib); mix(sl.imq[k]); }
    }
    return h;
  }

  /// Move a window's compact frames into a slab, consuming them as it goes (the compact store has
  /// to shrink in step, exactly as materializeWindow does, or the two representations coexist).
  PeakSlab toSlab(vector<CompactFrame>& frames)
  {
    PeakSlab sl;
    const size_t F = frames.size();
    sl.frame_off.resize(F + 1, 0);
    sl.rt_index.resize(F);
    size_t total = 0;
    for (const auto& f : frames) total += f.mzq.size();
    sl.tof.reserve(total); sl.inten.reserve(total); sl.imq.reserve(total);
    sl.b.resize(F);
    const TofAxis& ax = tofAxis();
    for (size_t i = 0; i < F; ++i)
    {
      CompactFrame& f = frames[i];
      sl.rt_index[i] = rtIndex(f.rt);
      // The compact store keeps m/z at a 1e-5 Da quantum, ~400x finer than one TOF bin at m/z 600,
      // so inverting the calibration recovers the ORIGINAL bin exactly -- no approximation is
      // introduced here, and the loader does not have to change.
      const double b = ax.factor(f.frame_id);
      sl.b[i] = b;
      sl.frame_off[i] = (uint32_t)sl.tof.size();
      for (size_t k = 0; k < f.mzq.size(); ++k)
        sl.tof.push_back(ax.tofOf((double)f.mzq[k] / MZ_Q, b));
      sl.inten.insert(sl.inten.end(), f.inten.begin(), f.inten.end());
      sl.imq.insert(sl.imq.end(), f.imq.begin(), f.imq.end());
      vector<uint32_t>().swap(f.mzq);
      vector<float>().swap(f.inten);
      vector<uint16_t>().swap(f.imq);
    }
    sl.frame_off[F] = (uint32_t)sl.tof.size();
    // The extension loop steps frames by index and takes the RT span as last minus first, so the
    // frames must be in RT order -- the OpenMS path sorts its spectra for the same reason. The
    // loader delivers them ordered; this checks rather than assumes.
    for (size_t i = 1; i < F; ++i)
      if (sl.rt_index[i] < sl.rt_index[i - 1])
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "window frames are not in retention-time order at frame", std::to_string(i));
    return sl;
  }


  /// Mass-trace detection on the integer m/z axis.
  ///
  /// The OpenMS detector accepts a peak when it lies within `mass_error_ppm` of a trace's RUNNING
  /// centroid, which needs a tolerance comparison against every live trace and makes the cost
  /// superlinear in peaks per band -- measured on dataset D: 4 bands 18:25 wall / 163 GB against 12 bands
  /// 12:53 / 103 GB, i.e. the same work costs more when the bands are bigger. On a uniform integer
  /// axis the candidates for a peak are found by ARITHMETIC: a live trace can only accept indices
  /// within `ppmSpan` of its centroid, so bucketing by `mz >> shift` with a bucket wider than the
  /// largest span puts every candidate in the peak's own bucket or its two neighbours.
  ///
  /// This is a DIFFERENT ALGORITHM, not a reimplementation: it will not reproduce the OpenMS
  /// detector peak for peak. It is the default (`trace:detector=integer`), and must be judged on
  /// identified peptides with both search engines, never on a spectrum digest. On dataset D and dataset A the
  /// two detectors agree on ~85% of the union of identified peptides, and which one identifies more
  /// flips with both the search engine and the file; the case for it is ~40% less memory.
  ///
  /// Mobility is part of the identity, as in the IM-aware OpenMS path: a candidate must also be
  /// within `im_tol` of the trace's mobility. Valley splitting is NOT done here -- the existing
  /// ElutionPeakDetection step runs afterwards on the traces this produces.
  /// Mass-trace detection on integer arrays, following OpenMS `MassTraceDetection::run_` step for
  /// step -- written from that source, not from a description of it. The first version was written
  /// from a description and lost 92% of the peptides; the review that found why is in
  /// docs/reviews/integer-tracing-2026-09-03/ and the corrected semantics are:
  ///   * seeds in intensity order, ties later-index-first (stable ascending sort, reverse walk);
  ///   * the two directions extend INTERLEAVED, one frame down then one up, sharing one centroid;
  ///   * acceptance is the CLOSEST peak within +-3 sd of the running m/z centroid and within the
  ///     mobility tolerance of the running mobility centroid; if that closest peak is already
  ///     taken, the frame is a MISS -- the next-closest is not tried;
  ///   * sd starts at centroid*ppm and is RE-ESTIMATED from the accepted peaks (reestimate_mt_sd);
  ///   * a direction stops after more than 5 consecutive misses; an empty frame is not a miss;
  ///   * the trace is valid if its RT span is within [min,max] length and
  ///     size / (frames visited - trailing misses) >= min_sample_rate (0.5 by default);
  ///   * peaks are marked taken only when the trace is valid;
  ///   * the seed is counted twice in the centroid, as run_ does (init, then update again).
  /// Peaks stay integers (flight-time bin, quantised mobility); the running centroid and sd are
  /// doubles in m/z, exactly as in OpenMS, so the tolerance semantics are the same. The bin<->m/z
  /// conversion is one calibration call per frame step plus one per ACCEPTED peak -- never per
  /// visited peak. Valley splitting runs afterwards on the result, as it does for the OpenMS path.
  /// Per-window, read-only preparation shared by every band task: frame of each peak, the
  /// noise-membership mask, per-frame "has any member" (OpenMS deletes sub-noise peaks first and
  /// only then asks whether a frame is empty), and the intensity-ordered seed list. Built ONCE.
  /// The first corrected version rebuilt all of this inside each of the 12 band tasks -- twelve
  /// full-window sorts per window and twelve sets of P-sized arrays alive at once -- which is
  /// exactly where its extra time (9:28 vs 7:22) and extra memory (174 vs 166 GB) came from.
  struct TracePrep
  {
    vector<char>     frame_live;  ///< F: has at least one member
    vector<uint32_t> order;       ///< seeds, intensity desc, later-index-first among equals
  };
  TracePrep prepareTracing(const PeakSlab& sl, double noise, double snr)
  {
    TracePrep tp;
    const size_t F = sl.frames(), P = sl.tof.size();
    // No per-peak frame map and no membership byte: a peak's frame is the loop variable at every
    // hot site (and one binary search at the two that lack it), and membership IS
    // `inten > noise`. Measured: 100% of above-noise peaks are seeds at the shipped snr of 1.0, so
    // `order` was reserved at P/4 and grew by doubling to ~2P -- count first, reserve exactly.
    tp.frame_live.assign(F, 0);
    const double min_apex = snr * noise;
    size_t n_seed = 0;
    for (size_t f = 0; f < F; ++f)
      for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
        if ((double)sl.inten[k] > noise) { tp.frame_live[f] = 1; if ((double)sl.inten[k] > min_apex) ++n_seed; }
    tp.order.reserve(n_seed);
    for (uint32_t k = 0; k < P; ++k)
      if ((double)sl.inten[k] > noise && (double)sl.inten[k] > min_apex) tp.order.push_back(k);
    // OpenMS: stable_sort ascending by intensity, then iterate in reverse -> among equal
    // intensities the LATER peak seeds first.
    sort(tp.order.begin(), tp.order.end(), [&](uint32_t a, uint32_t b) {
      if (sl.inten[a] != sl.inten[b]) return sl.inten[a] > sl.inten[b];
      return a > b;
    });
    return tp;
  }

  vector<Trace> detectTracesInteger_(const PeakSlab& sl, const TracePrep& tp, TraceStore& store, double mass_ppm,
                                     double im_tol, double noise, double min_len_sec, double min_sample_rate,
                                     double max_len_sec, double max_span_sec, TofIdx band_lo = 0, TofIdx band_hi = 0)
  {
    vector<Trace> out;
    const size_t F = sl.frames(), P = sl.tof.size();
    if (F == 0 || P == 0) return out;
    const bool banded = band_hi > band_lo;
    const TofAxis& ax = tofAxis();
    const double msr = min_sample_rate >= 0.0 ? min_sample_rate : 0.5;   // -1 = OpenMS default
    const double imq_tol = im_tol * IM_Q;
    const Size max_missed = 5;                                            // trace_termination_outliers
    // a peak's frame, for the two sites that do not already have it in a loop variable
    auto frameOfPeak = [&](uint32_t k) {
      return (uint32_t)(std::upper_bound(sl.frame_off.begin(), sl.frame_off.end(), k) - sl.frame_off.begin() - 1);
    };
    const auto member = [&](uint32_t k) { return (double)sl.inten[k] > noise; };

    vector<uint64_t> visited((P + 63) / 64, 0);   // one BIT per peak: 8x less than the char array
    auto seen = [&](size_t k) { return (visited[k >> 6] >> (k & 63)) & 1u; };
    auto mark = [&](size_t k) { visited[k >> 6] |= (uint64_t)1 << (k & 63); };
    auto mzAt = [&](uint32_t k, size_t f) { return ax.mzOf(sl.tof[k], sl.b[f]); };
    auto iwm = [](double v, double w, double& c, double& cnt, double& den) {   // updateIterativeWeightedMean_
      const double ct = 1.0 + (w * v) / cnt, dt = 1.0 + w / den;
      c *= (ct / dt); cnt += w * v; den += w;
    };
    auto sdRobust = [](double mz, double w, double mean, double& sd, double& wsum) {   // updateWeightedSDEstimateRobust
      const double d1 = std::log(wsum) + 2.0 * std::log(sd);
      const double d2 = std::log(w) + 2.0 * std::log(std::abs(mz - mean));
      const double denom = std::sqrt(std::exp(d1) + std::exp(d2));
      const double ws = wsum + w;
      const double t = denom / std::sqrt(ws);
      if (t > std::numeric_limits<double>::epsilon()) sd = t;
      wsum = ws;
    };
    /// findBestPeak_ + isPeakAcceptable_: the closest peak in m/z within +-3 sd of the running
    /// centroid whose mobility is inside the gate; -1 for "no candidate", and `taken` set when the
    /// closest one is already used, which OpenMS counts as a miss rather than trying the next.
    /// The +-3 sd window is converted to a BIN range once per frame step, so the walk is bounded
    /// by two calibration calls and a binary search; only candidates inside that window are
    /// converted back to m/z for the comparison.
    auto best = [&](size_t f, double cmz, double sd, double cim, bool& taken) -> long {
      taken = false;
      const uint32_t a = sl.frame_off[f], b = sl.frame_off[f + 1];
      if (b <= a) return -1;
      const double bf = sl.b[f];
      const double lo_mz = cmz - 3.0 * sd, hi_mz = cmz + 3.0 * sd;
      const TofIdx lo_t = ax.tofOf(lo_mz > 0.0 ? lo_mz : 1e-9, bf);
      const TofIdx hi_t = ax.tofOf(hi_mz, bf);
      const uint32_t s0 = (uint32_t)(lower_bound(sl.tof.begin() + a, sl.tof.begin() + b, lo_t) - sl.tof.begin());
      long bk = -1; double bd = std::numeric_limits<double>::infinity();
      for (uint32_t k = s0; k < b && sl.tof[k] <= hi_t; ++k)
      {
        if (!member(k)) continue;
        if (std::abs((double)sl.imq[k] - cim) > imq_tol) continue;
        // "Closest" is measured in m/z, as OpenMS does, NOT in bins: I had argued the calibration
        // is monotone so the two agree, and that is true only for candidates on the same side of
        // the centroid -- across it, the map's curvature can flip which is nearer. Measured: the
        // bin metric cost 1.4% of peptides (12,650 -> 12,466) and bought no time (the speed came
        // from parallel valley splitting). The walk is still bounded in bins; only candidates
        // inside the window pay the conversion.
        const double d = std::abs(ax.mzOf(sl.tof[k], bf) - cmz);
        if (d < bd) { bd = d; bk = (long)k; }          // strict <: first of equals wins
      }
      if (bk >= 0 && seen((size_t)bk)) { taken = true; return -1; }
      return bk;
    };

    struct Dir { Size hits = 0, missed = 0, scans = 0; bool active = true; };
    vector<pair<uint32_t, uint32_t>> got;   // (frame, peak); marked only if the trace is valid
    for (uint32_t seed : tp.order)
    {
      if (banded && (sl.tof[seed] < band_lo || sl.tof[seed] >= band_hi)) continue;   // not this band's core
      if (seen(seed)) continue;
      const size_t f0 = frameOfPeak(seed);
      const double w0 = sl.inten[seed], m0 = mzAt(seed, f0);
      double cmz = m0, cnt = w0 * m0, den = w0;
      iwm(m0, w0, cmz, cnt, den);                             // seed counted twice, as run_ does
      double cim = sl.imq[seed], cnt_im = w0 * cim, den_im = w0;
      iwm((double)sl.imq[seed], w0, cim, cnt_im, den_im);
      double sd = (cmz / 1e6) * mass_ppm, wsum = w0;
      got.clear(); got.emplace_back((uint32_t)f0, seed);
      uint32_t sfmin = (uint32_t)f0, sfmax = (uint32_t)f0;    // running span, for the length cap
      Dir dn, up;
      long fd = (long)f0, fu = (long)f0;
      while ((fd > 0 && dn.active) || (fu + 1 < (long)F && up.active))
      {
        for (int side = 0; side < 2; ++side)
        {
          Dir& d = side == 0 ? dn : up;
          long& fi = side == 0 ? fd : fu;
          if (!d.active) continue;
          if (side == 0 ? !(fi > 0) : !(fi + 1 < (long)F)) continue;
          const size_t f = (size_t)(side == 0 ? fi - 1 : fi + 1);
          // "empty" AFTER sub-noise deletion, as OpenMS tests it. SPEXTRACTOR_INT_EMPTY_RAW=1 tests
          // raw emptiness instead: the A/B that attributes the 12,650 -> 12,466 move between the
          // two corrected versions, which is the only semantic change between them.
          static const bool empty_raw = std::getenv("SPEXTRACTOR_INT_EMPTY_RAW") != nullptr;
          if (empty_raw ? (sl.frame_off[f + 1] > sl.frame_off[f]) : (bool)tp.frame_live[f])
          {
            bool taken = false;
            const long k = best(f, cmz, sd, cim, taken);
            // The length cap TERMINATES this direction rather than discarding the trace: an
            // over-length trace that is thrown away also throws away every peak in it, which is
            // why trace:max_trace_length_sec is falsified. Here the trace simply stops growing and
            // ElutionPeakDetection still splits what was found.
            static const bool span_terminate = std::getenv("SPEXTRACTOR_SPAN_TERMINATE") != nullptr;
            if (span_terminate && k >= 0 && max_span_sec > 0.0)
            {
              const uint32_t nlo = std::min(sfmin, (uint32_t)f), nhi = std::max(sfmax, (uint32_t)f);
              if (rtAxis()[sl.rt_index[nhi]] - rtAxis()[sl.rt_index[nlo]] > max_span_sec)
              { d.active = false; fi += (side == 0 ? -1 : 1); ++d.scans; continue; }
            }
            if (k >= 0)
            {
              sfmin = std::min(sfmin, (uint32_t)f); sfmax = std::max(sfmax, (uint32_t)f);
              const double mz = mzAt((uint32_t)k, f), w = sl.inten[k];
              iwm(mz, w, cmz, cnt, den);
              iwm((double)sl.imq[k], w, cim, cnt_im, den_im);
              sdRobust(mz, w, cmz, sd, wsum);   // reestimate_mt_sd: the tool never turns it off
              got.emplace_back((uint32_t)f, (uint32_t)k);
              ++d.hits; d.missed = 0;
            }
            else ++d.missed;
          }
          fi += (side == 0 ? -1 : 1);
          ++d.scans;
          if (d.missed > max_missed) d.active = false;
        }
      }
      // isTraceValid_
      uint32_t fmin = got[0].first, fmax = fmin;
      for (const auto& g : got) { fmin = std::min(fmin, g.first); fmax = std::max(fmax, g.first); }
      const double span_s = rtAxis()[sl.rt_index[fmax]] - rtAxis()[sl.rt_index[fmin]];
      if (span_s < min_len_sec) continue;
      if (max_len_sec > 0.0 && span_s > max_len_sec) continue;   // the tool passes 0 for "off"
      const Size total = dn.scans + up.scans + 1;
      const Size adjusted = total > dn.missed + up.missed ? total - dn.missed - up.missed : 1;
      if ((double)got.size() / (double)adjusted < msr) continue;
      for (const auto& g : got) mark(g.second);               // valid: now they are taken

      // the span: (store-local frame, intensity) points, plus the per-frame bin arena the valley
      // splitter needs to give each child ITS measured apex bin
      vector<pair<uint32_t, float>> pts;
      pts.reserve(got.size());
      uint32_t kap = got[0].second, fap = got[0].first;
      for (const auto& g : got) { pts.emplace_back(g.first, sl.inten[g.second]);
                                  if (sl.inten[g.second] > sl.inten[kap]) { kap = g.second; fap = g.first; } }
      Trace t = makeSpan(pts, store);
      store.bins.resize(store.inten.size(), 0u);
      for (const auto& g : got) store.bins[t.off + (g.first - t.frame0)] = sl.tof[g.second];
      t.tof = sl.tof[kap]; t.b = sl.b[fap];
      t.mz = cmz;                                              // the weighted centroid, as OpenMS reports
      t.rt = rtAxis()[sl.rt_index[fap]];
      t.im = IM_LO + cim / IM_Q;
      t.intensity = sl.inten[kap];
      out.push_back(std::move(t));
    }
    return out;
  }

  /// Integer traces -> MassTrace objects for EPD -> Trace again, with the flight-time bin
  /// re-attached from the apex frame's calibration. Only the TRACES become OpenMS objects -- the
  /// peak arrays never do, which is the point of the integer path.
  ///
  /// CHUNKED AND PARALLEL, like the OpenMS path's own valley-splitting step (do_split over
  /// bands*4 chunks as tasks). The previous version converted every trace of a window into a
  /// MassTrace (a list<Peak2D>, ~48 B/point) on ONE thread and then ran ElutionPeakDetection over
  /// all of them at once, also on one thread: a serial section of a stage that is 95% of the window
  /// loop, and the whole window's list-based copy alive at the same time. That is consistent with
  /// what was measured -- occupancy 44x against 81x for the OpenMS path, and MORE peak memory,
  /// after the per-band sorts and the per-visit calibration calls had already been removed and
  /// changed nothing.
  vector<Trace> splitIntegerTraces(vector<Trace>& in, double split_valleys, TraceStore& wst, int nchunk)
  {
    if (split_valleys <= 0.0 || in.size() < 2) return std::move(in);   // OpenMS branch: EPD needs > 1
    const TofAxis& ax = tofAxis();
    nchunk = std::max(1, std::min<int>(nchunk, (int)in.size()));
    vector<vector<Trace>> parts(nchunk);
    vector<TraceStore> stores(nchunk);                         // chunk-private arenas, merged after
    const size_t n_in = in.size();
    #pragma omp taskloop grainsize(1) default(shared)
    for (int ci = 0; ci < nchunk; ++ci)
    {
      const size_t lo = (size_t)ci * n_in / nchunk, hi = (size_t)(ci + 1) * n_in / nchunk;
      TraceStore& cst = stores[ci];
      cst.setFrames(wst.rt_index, &wst.b);
      // STREAMED: one MassTrace in flight, not the chunk's worth. The two-phase version built every
      // input trace of the chunk as a MassTrace and kept them all alive across the split, so the
      // input payload and the whole split output coexisted -- the two largest ledger lines after the
      // slab. Per-trace is exactly equivalent: with width_filtering=off and masstrace_snr_filtering
      // =false the cross-trace routines are unreachable, and the lock-free patch's nested path is
      // itself `for (i..) detectElutionPeaks_(mt_vec[i], out)` appending in the same order.
      //
      // epd is hoisted to chunk scope deliberately. splitValleys constructed one and did
      // getParameters/setValue x4/setParameters on EVERY call; moving that inside a per-trace loop
      // would be ~4e4 Param round-trips per chunk and would trade the memory win for a slowdown.
      ElutionPeakDetection epd;
      epd.setLogType(ProgressLogger::NONE);
      {
        Param ep = epd.getParameters();
        ep.setValue("chrom_fwhm", split_valleys);
        ep.setValue("width_filtering", "off");
        ep.setValue("masstrace_snr_filtering", "false");
        epd.setParameters(ep);
      }
      vector<Trace>& out = parts[ci];
      out.reserve(hi - lo);   // ponytail: parents, not children -- a heavily split chunk regrows it.
                              // Measured no regression (dataset D got faster and lighter); revisit only if it shows up.
      vector<Peak2D> pk;                                      // reused: a vector, not a per-point list
      vector<MassTrace> split;                                // reused: detectPeaks clears it itself
      for (size_t q = lo; q < hi; ++q)
      {
        const Trace& t = in[q];
        pk.clear();
        for (size_t k = 0; k < t.span(); ++k)                 // REAL points only: EPD smooths what it is given
        {
          if (!t.real(k)) continue;
          Peak2D pt; pt.setRT(t.rtAtSpan(k)); pt.setIntensity(t.xv(k));
          const uint32_t bin = wst.bins.empty() ? 0u : wst.bins[t.off + k];
          pt.setMZ(ax.ok && bin ? ax.mzOf(bin, wst.b[t.frame0 + k]) : t.mz);
          pk.push_back(pt);
        }
        MassTrace mt(pk);
        mt.setCentroidIM(t.im);
        epd.detectPeaks(mt, split);
        // Each child carries its own measured m/z: the parent's per-frame bin was written into
        // every Peak2D above, so the child's APEX PEAK m/z is that bin's m/z and the bin comes back
        // by inverting the calibration -- which is what the pre-span code did.
        //
        // The first version of this walked the parent list to look the bin up by frame. That was
        // wrong: parents within a chunk are in IM/m-z order, not frame order, and their frame ranges
        // overlap freely, so the walk could attribute a child to the wrong parent and hand it another
        // trace's bin. It cost 3.3% of peptides on dataset D while the OpenMS arm stayed digest-identical,
        // which is what localised it here. [2026-09-04]
        for (const MassTrace& cm : split)
        {
          Trace c = toTrace(cm, cst);
          if (c.np() == 0) continue;
          const uint32_t cf = c.frame0 + c.apex;               // apex frame, window-local
          c.b = wst.b.empty() ? ax.factor(0) : wst.b[cf];
          const double amz = cm.getSize() ? cm[cm.findMaxByIntPeak(false)].getMZ() : c.mz;
          c.tof = ax.ok ? ax.tofOf(amz, c.b) : 0u;
          out.push_back(std::move(c));
        }
      }
    }
    // the parents' spans are dead: rebuild the window arena from the chunk arenas
    vector<Trace>().swap(in);
    vector<float>().swap(wst.inten); vector<uint32_t>().swap(wst.bins);
    size_t tot = 0, arena = 0; bool any_bins = false;
    for (int ci = 0; ci < nchunk; ++ci)
    { tot += parts[ci].size(); arena += stores[ci].inten.size(); any_bins |= !stores[ci].bins.empty(); }
    // absorb() appends without reserving and the window arena was just swapped empty, so without
    // this the 48 chunk appends regrow it geometrically. Note this trades regrowth for allocating
    // the whole destination while every source is still live: it removes the copies, but whether
    // it lowers the PEAK is a measurement, not a deduction.
    wst.inten.reserve(arena); if (any_bins) wst.bins.reserve(arena);
    vector<Trace> out; out.reserve(tot);
    for (int ci = 0; ci < nchunk; ++ci)
    {
      const uint32_t base = wst.absorb(stores[ci]);
      for (auto& t : parts[ci]) { t.off += base; out.push_back(std::move(t)); }
      vector<Trace>().swap(parts[ci]);
    }
    return out;
  }

  /// Rebuild a PeakMap for ONE window from its compact frames, for MassTraceDetection.
  /// Frames are consumed (cleared) as they are materialised so the compact store shrinks in step.
  PeakMap materializeWindow(vector<CompactFrame>& frames)
  {
    PeakMap m;
    m.reserve(frames.size());
    for (auto& f : frames)
    {
      MSSpectrum s;
      s.setRT(f.rt);
      s.setMSLevel(1); // MassTraceDetection traces MS1-level only
      s.reserve(f.mzq.size());
      OpenMS::DataArrays::FloatDataArray ima;
      ima.setName(Constants::UserParam::ION_MOBILITY);
      ima.reserve(f.mzq.size());
      for (size_t i = 0; i < f.mzq.size(); ++i)
      {
        Peak1D p;
        p.setMZ((double)f.mzq[i] / MZ_Q);
        p.setIntensity(f.inten[i]);
        s.push_back(p);
        ima.push_back((float)(IM_LO + (double)f.imq[i] / IM_Q));
      }
      s.getFloatDataArrays().push_back(std::move(ima));
      m.addSpectrum(std::move(s));
      // release this frame's compact storage immediately
      vector<uint32_t>().swap(f.mzq);
      vector<float>().swap(f.inten);
      vector<uint16_t>().swap(f.imq);
    }
    vector<CompactFrame>().swap(frames);
    return m;
  }

  /// [stream] Pick-and-compact consumer: removes the one-shot-read memory floor.
  ///
  /// MEASURED: `loadExperiment` reaches 90.3 GB RSS / 99.0 GB peak at t=0s, BEFORE any compaction
  /// exists — so the compact store (10 B/peak, 20.9 GB for 2.19e9 peaks) can never lower the peak.
  /// The floor is the whole run being resident at once. BrukerTimsFile::loadDIAStreaming hands us one
  /// frame at a time, so we peak-pick and compact each frame on arrival and never hold the raw run.
  ///
  /// Windows are keyed by the SAME (lo,hi) isolation-window key used downstream, taken from the
  /// spectrum's own precursor, so routing is identical to the old split - not by swath_nr, whose
  /// ordering is the reader's business and need not match ours.
  /// Isolation window identity: (lo*100, hi*100, WindowGroup). The group is part of it because
  /// one diaPASEF scheme (Meier 2020 "py3", PXD017703) acquires the SAME m/z window in two window
  /// groups with shifted, overlapping ion-mobility slices ~1.7 s apart in the cycle. Keyed by m/z
  /// alone the two halves land in one window as all of group A then all of group B (RT not
  /// monotone, and a precursor outside the overlap would see every other frame empty). Every
  /// other scheme has one slice per m/z, so for those this key orders and partitions exactly as
  /// the old (lo, hi) pair did.
  using WinKey = std::array<int, 3>;

  class PickCompactConsumer : public FullSwathFileConsumer
  {
  public:
    PickCompactConsumer(PeakPickerIM& picker, map<WinKey, vector<CompactFrame>>& windows,
                        PeakMap& ms1, CompactStats& stats,
                        std::function<WinKey(double, double, uint32_t)> keyfn)
      : picker_(picker), windows_(windows), ms1_(ms1), stats_(stats), keyfn_(std::move(keyfn)) {}

    size_t frames_seen = 0;

    /// [perf-load] Pick+compact a BATCH of buffered frames in parallel, then append in arrival
    /// order. OpenMS hands frames over a strictly serial `for (fid..) consumeSpectrum(spec)` loop,
    /// so before this the peak picking of 1.26e9 peaks ran on ONE thread and was 53% of total wall
    /// (measured: 1046 s wall, parallel factor 1.0). Buffering here parallelises it without
    /// touching the reader. Appending strictly by buffer index keeps the output order identical to
    /// the serial path; the emitted spectra matched it in m/z and differed in ~2% of intensities at
    /// the 1e-5..1e-4 relative level (Sage count unchanged) -- see docs/RUNTIME-PLAN [pick-det].
    void flush_()
    {
      const int n = (int)buf_.size();
      if (n == 0) return;
      const double _tf0 = phase_clock_();
      std::vector<CompactFrame> compacted((size_t)n);
      std::vector<CompactStats> st((size_t)n);          // per-frame stats, merged serially below
      PeakPickerIM picker = picker_;                    // Param is a deep-copy value type
      // SPEXTRACTOR_PICK_SERIAL=1 runs this loop on one thread: the A/B that separates picker
      // scheduling-sensitivity from downstream chaos (kimi, plan review). [pick-det]
      static const bool par_pick = std::getenv("SPEXTRACTOR_PICK_SERIAL") == nullptr;
      static const bool sched_set = []{ const bool st = std::getenv("SPEXTRACTOR_PICK_STATIC") != nullptr; omp_set_schedule(st ? omp_sched_static : omp_sched_dynamic, st ? 0 : 1); return true; }(); (void)sched_set;
      std::exception_ptr pick_err;                      // an exception escaping an OMP region aborts
      // chunk 1: with chunk 8 a 64-frame batch made only 8 chunks = 8 workers (codex, plan review)
      #pragma omp parallel for firstprivate(picker) schedule(runtime) if(par_pick)   // SPEXTRACTOR_PICK_STATIC=1 -> static (history-vs-race probe)
      for (int i = 0; i < n; ++i)
      {
        try
        {
          auto& s = buf_[(size_t)i];
          ensureIMArrayName(s);
          if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
          ensureIMArrayName(s);
          compacted[(size_t)i] = compactify(s, st[(size_t)i]);
          s.clear(true);                                // raw frame dies here, inside the loop
        }
        catch (...)
        {
          #pragma omp critical(pick_err)
          if (!pick_err) pick_err = std::current_exception();
        }
      }
      if (pick_err) std::rethrow_exception(pick_err);
      for (int i = 0; i < n; ++i)
      {
        stats_ += st[(size_t)i];
        windows_[key_[(size_t)i]].push_back(std::move(compacted[(size_t)i]));
      }
      flush_wall_() += phase_clock_() - _tf0;
      buf_.clear(); key_.clear();
    }

  protected:
    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      // MS1 stays a PeakMap (MassTraceDetection needs it whole), but its pick is NOT small: 1,343
      // frames x ~690k raw peaks = 37% of all peaks, and it used to run serially here (~400 s of
      // LOAD, measured 2026-09-02). Buffer and pick in parallel like MS2; order is kept by index.
      ms1_buf_.push_back(std::move(s));
      ++frames_seen;
      if (ms1_buf_.size() >= kBatch) flushMS1_();
    }

    void flushMS1_()
    {
      const int n = (int)ms1_buf_.size();
      if (n == 0) return;
      const double _tf0 = phase_clock_();
      PeakPickerIM picker = picker_;
      std::exception_ptr pick_err;
      #pragma omp parallel for firstprivate(picker) schedule(dynamic, 1)
      for (int i = 0; i < n; ++i)
      {
        try
        {
          auto& s = ms1_buf_[(size_t)i];
          ensureIMArrayName(s);
          if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
          ensureIMArrayName(s);
        }
        catch (...)
        {
          #pragma omp critical(pick_err_ms1)
          if (!pick_err) pick_err = std::current_exception();
        }
      }
      if (pick_err) std::rethrow_exception(pick_err);
      for (auto& s : ms1_buf_) ms1_.addSpectrum(std::move(s));
      ms1_buf_.clear();
      flush_wall_() += phase_clock_() - _tf0;
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t /*swath_nr*/) override
    {
      const auto& prec = s.getPrecursors();
      if (prec.empty()) return;                      // no isolation window -> not routable
      const double c = prec[0].getMZ();
      const double lo = c - prec[0].getIsolationWindowLowerOffset();
      const double hi = c + prec[0].getIsolationWindowUpperOffset();
      flushMS1_();                                   // MS1 arrives first; pick it before any MS2 work
      key_.push_back(keyfn_(lo, hi, windowGroupOf(s.getNativeID())));
      buf_.push_back(std::move(s));                  // pick+compact deferred to flush_()
      ++frames_seen;
      if (buf_.size() >= kBatch) flush_();
    }

    void ensureMapsAreFilled_() override { flushMS1_(); flush_(); }   // pick up the tails

  private:
    /// deliberate: 256 frames is ~1 GB of raw buffer at ~37k peaks/frame -- 2.5 work items per
    /// thread at 100 threads (64 gave 8 chunks of 8 = 8 workers), small next to the 12 GB compact
    /// store. SPEXTRACTOR_PICK_BATCH overrides for the sweep.
    static const size_t kBatch;
    std::vector<MapType::SpectrumType> buf_;
    std::vector<WinKey> key_;
    PeakPickerIM& picker_;
    std::vector<MapType::SpectrumType> ms1_buf_;   // MS1 frames awaiting the parallel pick
    map<WinKey, vector<CompactFrame>>& windows_;
    PeakMap& ms1_;
    CompactStats& stats_;
    std::function<WinKey(double, double, uint32_t)> keyfn_;
  };
  const size_t PickCompactConsumer::kBatch =
    std::getenv("SPEXTRACTOR_PICK_BATCH") ? (size_t)std::atoi(std::getenv("SPEXTRACTOR_PICK_BATCH")) : 256;

  /// [wavelet] Stationary ("a trous") B3-spline wavelet smoothing of an elution profile.
  ///
  /// WHY. the reference implementation correlates a SMOOTHED precursor XIC against RAW fragment XICs (Nat Commun
  /// 16:95, Methods: 2D Gaussian pre-filter, then Savitzky-Golay for elongated peaks). We copied
  /// its 0.3 correlation threshold WITHOUT copying the smoothing. Pearson r between two noisy
  /// profiles is attenuated toward 0 as S/N falls, so an unsmoothed 0.3 is a systematically weaker
  /// constraint than a smoothed 0.3 -- and the attenuation is worst exactly at faint precursors,
  /// which is where the coverage analysis says we lose peptides (MISS precursors are 3.2x dimmer).
  /// Raising the threshold to compensate was tried and produced only an FDR artefact; smoothing
  /// attacks the cause instead of the symptom.
  ///
  /// A-TROUS, not a decimated DWT: no downsampling, so the output stays sample-aligned with the
  /// input and can be scattered onto the correlation grid at the SAME RT positions. That matters
  /// for correctness here, not just convenience -- the overlap gate counts grid points where the
  /// precursor is non-zero, so a transform that moved or resampled points would silently change
  /// what `min_correlation_points` means.
  ///
  /// B3-spline scaling filter [1,4,6,4,1]/16, dilated by 2^(j-1) at level j; c_J is a low-pass of
  /// support ~2^J samples. J comes from the MEASURED MS1 FWHM (see atrousLevels), never a guess.
  void atrousSmooth(const vector<double>& in, int levels, vector<double>& out)
  {
    static const double H[5] = {1.0 / 16, 4.0 / 16, 6.0 / 16, 4.0 / 16, 1.0 / 16};
    const int n = (int)in.size();
    out = in;
    if (n < 3 || levels <= 0) return;
    vector<double> tmp(n);
    for (int j = 1; j <= levels; ++j)
    {
      const int step = 1 << (j - 1);
      for (int i = 0; i < n; ++i)
      {
        double acc = 0.0;
        for (int k = -2; k <= 2; ++k)
        {
          int idx = i + k * step;
          // MIRROR both ends. Zero-padding would pull an edge-truncated peak toward zero and
          // fabricate a downslope -- the same bias the FWHM measurement excludes truncated peaks for.
          if (idx < 0) idx = -idx;
          if (idx >= n) idx = 2 * (n - 1) - idx;
          idx = std::max(0, std::min(n - 1, idx));   // degenerate after double reflection (tiny n)
          acc += H[k + 2] * out[idx];
        }
        tmp[i] = acc;
      }
      out.swap(tmp);
    }
  }

  /// Number of a-trous levels whose smoothing support best matches @p scale_sec at this profile's
  /// OWN sampling. Sampling is read per-trace rather than assumed: diaPASEF cycles are regular, but
  /// a valley-split trace can carry gaps, and a fixed level count would then smooth different
  /// physical widths on different traces.
  int atrousLevels(const Trace& tr, double scale_sec)
  {
    if (tr.np() < 5 || scale_sec <= 0.0) return 0;
    thread_local vector<double> prt, pv;
    packReal(tr, prt, pv);                             // the smoother's gaps are between REAL points
    vector<double> d;
    d.reserve(prt.size());
    for (size_t i = 1; i < prt.size(); ++i)
    {
      const double g = prt[i] - prt[i - 1];
      if (g > 0) d.push_back(g);
    }
    if (d.empty()) return 0;
    nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
    const double dt = d[d.size() / 2];                 // median, not mean: gaps are heavy-tailed
    if (dt <= 0.0) return 0;
    const double want = scale_sec / dt;                // desired support, in SAMPLES
    if (want < 2.0) return 0;                          // finer than one step: nothing to smooth
    const int j = (int)std::lround(std::log2(want));
    return std::max(0, std::min(j, 4));                // level 5 spans ~64 samples, wider than any peak here
  }

  /// Per-window scoring statistics from IM-sorted fragment traces (see FragStats).
  /// What the scorer needs per window, now that a fragment's profile is a span in the window's
  /// own frame sequence: the Pearson denominator G (distinct frames carrying a REAL fragment
  /// point -- not the window's frame count), per-fragment mean and 1/norm over real points, and
  /// the precursor-to-fragment nearest-frame table. No copy of any profile: the arena is the grid.
  struct FragStats
  {
    vector<double> mean, invnorm;   ///< per fragment, over the full grid of G frames
    vector<int>    nearest_local;   ///< MS1-store frame -> window frame with real fragment data, or -1
    vector<char>   used;            ///< window frame carries >= 1 real fragment point
    size_t G = 0;
  };

  FragStats buildFragStats(const vector<Trace>& frags, const TraceStore& wst, const TraceStore& ms1st, double delta_rt)
  {
    FragStats g;
    const size_t F = wst.frames();
    g.used.assign(F, 0);
    for (const auto& f : frags)
      for (size_t k = 0; k < f.span(); ++k) if (f.real(k)) g.used[f.frame0 + k] = 1;
    for (char u : g.used) g.G += u ? 1 : 0;
    const double G = (double)g.G;
    // Precursor points live on MS1 frames; each maps to the NEAREST window frame that carries real
    // fragment data, within gate:delta_rt, ties to the LATER frame -- the rule the old grid used.
    vector<uint32_t> uf; uf.reserve(g.G);
    for (uint32_t f = 0; f < F; ++f) if (g.used[f]) uf.push_back(f);
    g.nearest_local.assign(ms1st.frames(), -1);
    size_t hi = 0;
    for (size_t m = 0; m < ms1st.frames(); ++m)
    {
      const double t = ms1st.frame_rt[m];
      while (hi < uf.size() && wst.frame_rt[uf[hi]] < t) ++hi;
      double bd = delta_rt; int gi = -1;
      if (hi > 0)          { const double d = fabs(wst.frame_rt[uf[hi - 1]] - t); if (d <= bd) { bd = d; gi = (int)uf[hi - 1]; } }
      if (hi < uf.size())  { const double d = fabs(wst.frame_rt[uf[hi]] - t);     if (d <= bd) { bd = d; gi = (int)uf[hi]; } }
      g.nearest_local[m] = gi;
    }
    g.mean.resize(frags.size()); g.invnorm.resize(frags.size());
    #pragma omp taskloop default(shared)
    for (long ii = 0; ii < (long)frags.size(); ++ii)
    {
      const Trace& f = frags[(size_t)ii];
      double sum = 0, sumsq = 0;
      for (size_t k = 0; k < f.span(); ++k) { const double v = f.xv(k); if (v > 0.0) { sum += v; sumsq += v * v; } }
      g.mean[ii] = sum / G;
      const double var = sumsq - G * g.mean[ii] * g.mean[ii];
      g.invnorm[ii] = var > 0 ? 1.0 / sqrt(var) : 0.0;  // 0 => constant/degenerate -> never correlates
    }
    return g;
  }
}

class TOPPSpeXtractor : public TOPPBase
{
public:
  TOPPSpeXtractor() :
    // official=false: SpeXtractor is a STANDALONE tool that links against OpenMS, not a tool IN
    // OpenMS, so it is deliberately absent from ToolHandler's official list. TOPPBase rejects an
    // unregistered name when official is left true -- which is how the rename surfaced that the
    // old name had been passing itself off as an official TOPP tool.
    TOPPBase("SpeXtractor", "Extracts pseudo-MS/MS spectra from diaPASEF (ion-mobility DIA) data.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input diaPASEF data (ion-mobility DIA; 1/K0 / VSSC): mzML, mzPeak, or a Bruker .d directory.");
#ifdef SPEXTRACTOR_WITH_MZPEAK
    setValidFormats_("in", {"mzML", "mzpeak", "d"});
#else
    setValidFormats_("in", {"mzML", "d"});      // stock OpenMS does not know the mzpeak format
#endif
    registerOutputFile_("out", "<file>", "", "Output pseudo-MS2 spectra. The FORMAT FOLLOWS THE "
                        "EXTENSION: .mzpeak (default) or .mzML. mzPeak is columnar and much smaller; "
                        "mzML is what DDA search engines read today, so pass an .mzML name (or "
                        "-out_type mzML) if the next step is a search.");
#ifdef SPEXTRACTOR_WITH_MZPEAK
    setValidFormats_("out", {"mzpeak", "mzML"});
#else
    setValidFormats_("out", {"mzML"});
#endif
    registerStringOption_("out_type", "<type>", "", "Force the output format instead of taking it "
                          "from the extension of -out. Empty = follow the extension; mzPeak if the "
                          "extension says nothing.", false);
    setValidStrings_("out_type", {"", "mzpeak", "mzML"});

    registerTOPPSubsection_("gate", "Precursor-to-fragment co-localization gate");
    // Conservative, documented defaults; these are placeholders pending a real-data sweep. [H-3]
    registerDoubleOption_("gate:delta_im", "<1/K0>", 0.01, "Max |precursor - fragment| ion mobility difference (1/K0). (the reference implementation deltaApexIM default 0.01)", false);
    registerDoubleOption_("gate:delta_rt", "<sec>", 3.0, "Max |precursor - fragment| apex RT difference (seconds).", false);
    registerDoubleOption_("gate:min_correlation", "<0..1>", 0.3, "Min Pearson correlation of elution profiles. (the reference implementation ms1MS2Corr default 0.3; raise for cleaner but sparser spectra)", false);
    registerIntOption_("gate:min_correlation_points", "<n>", 3, "Min overlapping XIC points required to accept a correlation.", false, true);

    registerTOPPSubsection_("assembly", "Pseudo-spectrum assembly");
    registerIntOption_("assembly:min_fragments", "<n>", 3, "Emit a pseudo-spectrum only if it has at least this many fragments.", false);
    registerIntOption_("assembly:max_fragments", "<n>", 500, "Keep at most this many (top-ranked) fragments per pseudo-spectrum.", false);
    registerFlag_("assembly:competitive", "Assign each fragment to only its single best-correlating precursor (de-chimerize). Shorthand for rp_max=1; the harder endpoint of the rp_max dial. [bench]");
    registerDoubleOption_("assembly:apportion", "<p>", 0.0, "[route-4] Split a shared fragment's INTENSITY across the precursors that claim it, weight w_i = corr_i^p / sum(corr^p), instead of copying it at full intensity into all of them (mean fan-out 6.45). Fragment COUNT is unchanged - only the weighting - so it sidesteps every falsified count knob (rank-pruning, competitive, min_corr, min_corr_pts). 0 = OFF. Try 1 (linear) or 2 (sharper); p->inf degenerates to competitive, which IS falsified. Share-all path only.", false);
    registerStringOption_("assembly:rank_by", "<mode>", "correlation", "[rank] Which score decides WHICH fragments survive the assembly:max_fragments cap. 'correlation' (current) keeps the top-N best-CORRELATING; 'intensity' keeps the top-N most INTENSE, as the reference implementation does (Nat Commun 16:95: 'only the top N highest intensity peaks (RF max, default 500)'). This matters because the search engine re-ranks BY INTENSITY: Sage keeps max_peaks=150 by intensity and MSFragger reads the full list. A high-intensity, mediocre-correlation fragment that our cap discards would have been scored by the engine, so a correlation-ranked cap can throw away exactly the peaks that get matched. 70-75% of our spectra sit AT the 500 cap, so the choice is live for most spectra, not a corner case.", false);
    setValidStrings_("assembly:rank_by", {"correlation", "intensity"});
    registerDoubleOption_("assembly:im_weight_sigma", "<1/K0>", 0.0, "[recall] Down-weight each emitted fragment's intensity by a Gaussian of its ion-mobility distance from the precursor, exp(-(dIM)^2 / 2 sigma^2). TIMS precedes fragmentation so a TRUE fragment shares the precursor's 1/K0 exactly; a co-isolating interferent at the edge of the gate:delta_im band is suppressed. WITHIN-spectrum (no fragment removed or redistributed across precursors) so it escapes the fan-out depletion that falsified competitive/apex/deconvolution, and directly un-buries faint true fragments under the engine's intensity top-150. 0 = off. Try ~0.005 (half the 0.01 IM gate). [bench]", false);
    registerDoubleOption_("assembly:corr_power", "<k>", 2.0, "[Q2] Multiply each emitted fragment's intensity by corr^k (k>=0), where corr is the fragment-precursor co-elution score. EMIT-ONLY (excluded from the max_fragments cap key), so it preserves peak membership and only re-ranks the intensities the search engine sees toward well-correlating (true co-eluting) fragments and away from bright interferents; WITHIN-spectrum so it escapes fan-out depletion. DEFAULT 2.0 since 2026-07-28: validated on dataset D (tuning) + dataset A + dataset B (held out) at frozen k=2 -- +8-10%% Sage / +5-8%% MSFragger over base on all three, both engines, at controlled empirical FDR; k=2 sits one step below the ~k=3 Sage peak (the safe side of the crossover where c^k starts crushing the faint tail). 0 = off (pre-2026-07-28 behaviour). Dimensionless + bounded (c in [0.3,1]) so it transfers across instruments, unlike assembly:im_weight_sigma. [bench]", false);
    registerIntOption_("assembly:rp_max", "<n>", 0, "FALSIFIED on dataset D (2026-07, all cross-precursor redistributions lost peptides; note the A/B ran with corr_power silently OFF until 2026-09-02 [E5]) -- keep 0. Soft RP rank-pruning (DIA-Umpire style): keep each fragment in only its top-N best-correlating precursors. 0 = unlimited (share-all, the default); 1 = competitive (winner-take-all). Small N (2-8) is the useful range given typical fan-out ~6. RFmax (per-precursor fragment cap) is 'max_fragments'. Runs emit a fan-out histogram so you can see whether a given rp_max actually bites. [bench]", false);
    registerFlag_("assembly:open_search_safe", "[way-4] ANNOTATION-SAFE output for BLIND/OPEN search. In open search a wrong charge is catastrophic in a way closed search hides: assigning z' instead of z shifts the neutral mass by dM=(z'-z)(m/z-1.0073), so ONE charge step at m/z 500 fabricates a ~499 Da 'modification'. Closed search barely notices (charge unset cost only 1.3%: 8,123->8,019) because the engine enumerates charge; an open engine trusts the annotation. This flag therefore: (a) emits charge UNSET for every precursor whose charge was GUESSED rather than isotope-supported, so the engine enumerates instead of trusting a default; (b) records the isotope-offset hypothesis as a meta value so delta-mass inference can correct BEFORE assigning a modification, rather than reporting phantom +-1.00335/+-2.0067 Da shifts. Use for open/variant search; leave off for closed-search peptide counting.");
    registerIntOption_("assembly:default_charge", "<z>", 2, "Charge for precursors with no isotope support: >0 assigns that charge (the reference implementation-style; 2 = most common); 0 = leave charge UNSET so the search engine searches its charge range and picks the right one (recovers 3+/4+ that a forced 2 would lose).", false);

    registerTOPPSubsection_("trace", "Mass-trace detection (see MassTraceDetection)");
    registerDoubleOption_("trace:mass_error_ppm", "<ppm>", 15.0, "m/z tolerance for trace detection.", false);
    registerDoubleOption_("trace:noise_threshold_int", "<int>", 100.0, "MS1 (precursor) noise intensity threshold.", false);
    registerDoubleOption_("trace:ms2_noise_threshold_int", "<int>", 10.0, "MS2 (fragment) noise threshold. Lower it (e.g. 50, 20) to recover weak fragments toward the reference implementation-like density; the grid-Pearson optimization keeps that affordable.", false);
    registerDoubleOption_("trace:min_length_sec", "<sec>", 3.0, "Minimum MS1 mass-trace length (seconds).", false, true);
    // the reference implementation runs MS2 with deliberately RELAXED constraints vs MS1 (Nat Commun 16:95). The OpenMS
    // apex gate is chrom_peak_snr*noise_threshold_int, so the stock snr=3.0 silently TRIPLES the
    // fragment threshold (a "noise 30" run really gated apices at 90). Decouple both for MS2. [recall]
    registerDoubleOption_("trace:ms1_chrom_peak_snr", "<x>", 3.0, "Apex signal-to-noise multiplier for MS1 traces (effective apex threshold = this x noise_threshold_int). OpenMS default 3.0.", false, true);
    registerDoubleOption_("trace:ms2_chrom_peak_snr", "<x>", 1.0, "Apex signal-to-noise multiplier for MS2 (fragment) traces. 1.0 = the ms2_noise_threshold_int you set IS the apex threshold; 3.0 would reproduce the old hidden 3x. Lower = more weak fragments recovered.", false);
    registerDoubleOption_("trace:ms2_min_length_sec", "<sec>", 0.0, "Minimum MS2 (fragment) mass-trace length (seconds). 0 = no length filter (the reference implementation-style relaxed MS2); MS1 keeps trace:min_length_sec.", false);
    registerStringOption_("trace:detector", "<mode>", "integer", "Mass-trace detector. 'integer' "
                          "(default) works on the global integer m/z axis: it finds a peak's candidate traces by "
                          "arithmetic instead of a tolerance search, and never converts the compact store back "
                          "to double m/z. It needs the vendor flight-time calibration and falls back to 'openms' "
                          "without it. 'openms' is OpenMS MassTraceDetection on a materialised PeakMap. These are "
                          "DIFFERENT algorithms, not reimplementations: they agree on about 85% of the union of "
                          "identified peptides, and which one identifies more depends on the search engine AND the "
                          "file. Integer costs ~40% less memory. See docs/MZ-AXIS-DESIGN.md.", false);
    setValidStrings_("trace:detector", {"openms", "integer"});
    registerStringOption_("trace:mz_estimator", "<mode>", "apex", "Reported m/z of a mass trace. "
                          "'apex' (default) uses the most intense peak's m/z; 'mean' the OpenMS intensity-weighted "
                          "centroid; 'median' the median. apex was chosen on measured identifications with both "
                          "search engines -- see docs/dataset D-BASELINE.md. Affects EVERY reported precursor and "
                          "fragment m/z.", false);
    setValidStrings_("trace:mz_estimator", {"apex", "mean", "median"});
    registerDoubleOption_("trace:ms2_split_valleys", "<chrom_fwhm>", 7.0, "[way-2] Split MS2 mass traces at chromatographic local minima via ElutionPeakDetection (which MassTraceDetection never does - it only terminates on outliers, so two peptides 25 s apart within the m/z tolerance MERGE). Value is chrom_fwhm in seconds; it sets both the Savitzky-Golay window and the local-extrema half-window. 0 = OFF. Try 6-8 (peaks here are 5-30 s). width_filtering is forced off so this only SPLITS, never deletes.", false);
    registerFlag_("diag:selftest_wavelet", "[wavelet] Run assertions on the a-trous smoother and exit. Tests the SHIPPED functions, not a copy: DC preservation (filter must sum to 1, or every correlation is rescaled), noise reduction, peak-height retention at the FWHM-derived scale, mirror edges (a truncated peak must not be dragged toward zero), and level selection from sampling. Exits non-zero on failure.", true);
    registerFlag_("diag:selftest_arena", "[arena] Run assertions on the MS1 arena compaction (compactUnreferenced) and exit: kept spans must survive a container order that differs from their offset order", true);
    registerDoubleOption_("trace:wavelet_smooth", "<x>", 0.0, "[wavelet] Smooth the PRECURSOR XIC with a stationary (a-trous) B3-spline wavelet before correlating fragments against it, at a scale of x TIMES THE MEASURED MS1 FWHM. the reference implementation smooths its precursor profile (2D Gaussian + Savitzky-Golay) and then applies a 0.3 correlation threshold; we copied the 0.3 but not the smoothing. Pearson r between two NOISY profiles is attenuated toward 0 as S/N falls, so an unsmoothed 0.3 is a systematically weaker constraint than a smoothed 0.3 -- and worst at faint precursors, which is where the coverage analysis says peptides are lost. Fragment XICs stay RAW (as in the reference implementation) and emitted peak intensities are unaffected. 0 = off. Try 0.25-0.5 (measured MS1 FWHM is 3.61 s on dataset D, so 0.5 ~ 1.8 s ~ 1.3 cycles).", false);
    registerDoubleOption_("trace:split_valleys_fwhm", "<x>", 0.0, "[fwhm] Set the valley-splitting window as a MULTIPLE OF THE MEASURED MS1 FWHM instead of absolute seconds. Overrides trace:ms2_split_valleys when > 0 (the MS2 window is the one it sets). Absolute seconds cannot be right across methods: 14.0 s is +0.0% peptides on a 31 min gradient and -1.06% on a 5.6 min one. Measured MS1 FWHM is 3.61 s on dataset D, so the old 7.0 s default is ~1.9x FWHM and 14.0 s is ~3.9x. 0 = off (use absolute seconds).", false);
    registerDoubleOption_("trace:ms1_split_valleys", "<chrom_fwhm>", 7.0, "Split MS1 mass traces at chromatographic valleys via ElutionPeakDetection. DEFAULT 7.0 s: MassTraceDetection never splits at a local minimum, so two peptides eluting ~25 s apart within the m/z tolerance merge into ONE multi-modal trace, producing a merged precursor with the wrong monoisotope and charge. Measured +36.7% peptides on dataset B and +37.4% on dataset A at peptide-level FDR. Costs ~1.7x wall time and ~25% peak RAM. 0 = off.", false);
    registerDoubleOption_("trace:max_span_sec", "<sec>", 120.0, "Trim a mass trace to at most this many "
                          "seconds around its apex, AFTER detection and valley splitting have decided what the peak "
                          "is. Applies to precursor and fragment traces alike. An unbounded trace has no locality: "
                          "15.2%% of dataset D precursor traces span >45 s and the longest spans 1,244 s of a 1,860 s "
                          "gradient, while real chromatographic peaks here are 5-30 s -- so what a 120 s cap removes "
                          "is the far tail of a blob ElutionPeakDetection did not split. Trimming AFTER the fact is "
                          "deliberate: it leaves detection's peak ownership untouched, unlike terminating extension "
                          "(SPEXTRACTOR_SPAN_TERMINATE=1, diagnostic), and unlike trace:max_trace_length_sec, which "
                          "DISCARDS an over-length trace and deletes all signal at that m/z (falsified). 0 = no cap. "
                          "CHANGES OUTPUT.", false);
    setMinFloat_("trace:max_span_sec", 0.0);
    registerDoubleOption_("trace:max_trace_length_sec", "<sec>", -1.0, "DANGER - NOT a guard. MassTraceDetection::isTraceValid_ RETURNS FALSE on over-length traces and peaks are marked visited only for VALID traces, so an over-length blob is rejected, re-seeded by the next apex, and rejected again - emitting NOTHING at that m/z. Setting this DELETES signal (the falsified direction). Kept only for diagnostics; leave at -1. Splitting requires ElutionPeakDetection, not this. [merged-trace]", false, true);
registerIntOption_("trace:frame_aggregation_ms1_n", "<n>", 1, "[cross-frame] Same as frame_aggregation_n but for the MS1 (PRECURSOR) map. This is the one that can add DEPTH: measurement shows we identify MORE spectra than the reference implementation (3.7%% vs 2.0%% of spectra) but cover FEWER distinct peptides (4.6 vs 1.46 PSMs/peptide), i.e. we miss low-abundance PRECURSORS. Aggregating MS2 alone cannot fix that. 1 = off.", false);
    registerIntOption_("trace:frame_aggregation_n", "<n>", 1, "[cross-frame] Number of ADJACENT RT frames (= adjacent diaPASEF CYCLES, since each isolation window is sampled once per cycle) summed into a composite (m/z x 1/K0) frame BEFORE MS2 peak/trace detection, the reference implementation-style. 1 = off (current behaviour). 5 = user request. Sliding window, stride 1, so the RT axis and point count are preserved; note adjacent output points then share input frames and are AUTOCORRELATED, which inflates Pearson r - benchmark, do not assume.", false);
    registerDoubleOption_("trace:ms2_min_sample_rate", "<r>", -1.0, "MS2 min fraction of visited scans that must contain a peak (MassTraceDetection min_sample_rate). OpenMS default 0.5 hard-deletes gappy fragment traces; MEASURED: 0.3 and 0.1 both LOSE peptides (8,865 -> 8,735 at 0.3), so -1 (= OpenMS 0.5) is the default. Lowering it produces longer, more-merged traces with wrong apexes.", false);
    registerDoubleOption_("trace:ms1_min_sample_rate", "<r>", -1.0, "MS1 min_sample_rate. -1 = OpenMS default 0.5 (strict, MS1 stays strict).", false, true);

    // NB: "threads" is reserved by TOPPBase (the standard -threads option), so this lives under "perf".
    // [route-1] feature-level consolidation: collapse spectra that are the SAME chromatographic
    // feature emitted more than once (RT slices / IM sub-ranges / charge hypotheses). Targets the
    // measured 4.6 vs 1.46 PSMs-per-peptide gap. OFF by default (delta_rt 0) - it is a deletion.
    registerTOPPSubsection_("consolidate", "Feature-level de-duplication of emitted spectra");
    registerTOPPSubsection_("merge", "Recombine split spectra (union, not selection)");
    registerDoubleOption_("merge:rt_window", "<sec>", 0.0, "FALSIFIED on dataset D (every merge/collapse variant lost peptides) -- keep 0 = OFF. [merge] Merge spectra sharing a precursor identity within this RT window into ONE spectrum by UNION of peaks; matching peaks take the MAX intensity by default (merge:sum_intensity sums instead). This is not consolidate: consolidate SELECTS the richest member and DISCARDS the rest (measured -12.7% peptides), while 41% of peaks are present in only one member. Measured MS1 FWHM is 3.61 s (2.61 cycles) and flat across the gradient, so ~3.6 s is one peak width -- the principled window. 0 = off.", false);
    registerDoubleOption_("merge:mass_ppm", "<ppm>", 20.0, "[merge] Precursor m/z tolerance for deciding two spectra share an identity.", false, true);
    registerDoubleOption_("merge:delta_im", "<1/K0>", 0.02, "[merge] Precursor ion-mobility tolerance.", false, true);
    registerDoubleOption_("merge:mz_tol", "<Da>", 0.01, "[merge] Peak m/z tolerance when summing the union. Peaks within this of the running intensity-weighted centroid are combined.", false, true);
    registerDoubleOption_("merge:min_cosine", "<x>", 0.0, "[merge] Require this normalised dot product between two spectra before merging them. Coordinate proximity is NOT the problem (the coordinate gate already admits only 1.8% different-peptide pairs); the problem is that only ~8% of emitted spectra identify at all, so a coordinate gate merges good spectra with unidentifiable ones whose peaks dilute the match. Measured reference: within-cycle pairs 0.845, same peptide across cycles 0.687, unrelated 0.007. 0 = off.", false);
    registerFlag_("merge:sum_intensity", "[merge] Sum intensities of matching peaks instead of taking the MAX. Summing distorts the merged spectrum: members share only ~37% of peaks (Jaccard 0.374), so shared peaks get Nx and unique peaks 1x, wrecking the fragment intensity ratios the engine scores on. MAX takes each ion's best single observation. Off (i.e. MAX) by default.");
    registerFlag_("merge:any_charge", "[merge] Also merge spectra whose precursor charges disagree. Off by default: 93.2% of split pairs already agree in charge, so merging across charges mostly pools genuinely different precursors.");
    registerDoubleOption_("consolidate:delta_rt", "<sec>", 0.0, "FALSIFIED on dataset D (-12.7% peptides at 5 s; every merge/collapse variant lost) -- keep 0 = OFF. Consolidate spectra whose precursors agree within this RT window (and mass_ppm/delta_im) into ONE feature, keeping the richest spectrum. Attacks emission multiplicity (ours 4.6 PSMs/peptide vs the reference implementation 1.46) without changing fragment content.", false);
    registerDoubleOption_("consolidate:mass_ppm", "<ppm>", 20.0, "Precursor m/z tolerance for feature consolidation.", false, true);
    registerDoubleOption_("consolidate:delta_im", "<1/K0>", 0.02, "Precursor IM tolerance for feature consolidation. Set larger than the 12 m/z x 2 IM sub-range split to merge a feature emitted from BOTH IM sub-ranges.", false, true);
    registerFlag_("consolidate:same_charge_only", "Only merge spectra that also agree in charge (safer; leaves genuine multi-charge observations intact).");

    registerStringOption_("perf:stream_load", "<true/false>", "true",
                          "[stream] Read the .d frame-by-frame via BrukerTimsFile::loadDIAStreaming, "
                          "peak-picking and compacting each frame on arrival, instead of loadExperiment() "
                          "which holds the entire run (MEASURED 90.3 GB RSS at t=0s -- that one-shot read IS "
                          "the memory floor and no downstream compaction can lower it, so with this off the "
                          "floor alone exceeds the whole tool's current peak). ON BY DEFAULT since 2026-09-04; "
                          "every benchmark figure on record was taken with it on. Bruker .d only -- ignored "
                          "for other inputs, which fall back to the one-shot reader. "
                          "CHANGES OUTPUT: this is not a pure performance switch. Measured on dataset D 2026-09-04, "
                          "true = 655,776 spectra in 5:11 at 78.4 GB, false = 656,371 in 10:30 at 121.8 GB, "
                          "spectrum lists NOT identical (0.09% more spectra with the one-shot reader). The two "
                          "readers pick peaks on different frame groupings; which is closer to the truth is "
                          "UNMEASURED. Everything this project has benchmarked used true.", false);
    setValidStrings_("perf:stream_load", {"true", "false"});
    registerIntOption_("trace:native_ms2_neighbors", "<n>", 0, "[stream] NATIVE cross-frame aggregation in the READER (BrukerTimsFile Config dia_ms2_n_neighbors): adjacent frames each side, 0=off, 1=3-frame sum, 2=5-frame sum. This sums RAW frames BEFORE centroiding, which is what the reference implementation does; the tool-level frame_aggregation_n sums AFTER centroiding and lost 26% of peptides. Requires perf:stream_load.", false);
    registerIntOption_("trace:native_ms1_neighbors", "<n>", 0, "[stream] Same, for MS1 frames (BrukerTimsFile Config ms1_n_neighbors). Requires perf:stream_load.", false);

    registerTOPPSubsection_("perf", "Concurrency / memory tradeoff");
    registerStringOption_("gate:coelution", "<mode>", "pearson", "[B0] Fragment-precursor co-elution statistic. 'pearson' takes variance over the FULL RT grid, so a fragment present in a handful of points is scored against ~1300 implicit zeros -- it rewards sparse fragments that spike where the precursor spikes, and its scale depends on grid size rather than evidence. 'logoverlap' is the AlphaDIA-style alternative: profile overlap in log space, summed only where the precursor has signal, normalised to [0,1]. The two are on DIFFERENT SCALES, so gate:min_correlation must be retuned when switching.", false);
    setValidStrings_("gate:coelution", {"pearson", "logoverlap"});
    registerFlag_("gate:variance_support", "[Q1] Compute the Pearson co-elution variance over the UNION SUPPORT of the precursor and fragment profiles (points where either is non-zero) rather than the full RT grid G (~1300). Fixes the [B0] zero-padding: a fragment present in a handful of grid points is no longer scored against ~1300 implicit zeros, and the score no longer depends on grid size. Only affects gate:coelution=pearson. Retune gate:min_correlation when enabling -- the score scale changes. [bench]");
    registerIntOption_("perf:ms1_trace_bands", "<n>", 12, "[parallel] Band-parallel MS1 mass-trace detection, same exact-partition scheme as perf:trace_bands. MS1 tracing was running SINGLE-THREADED (the call simply never passed a band count) while MS2 was banded, so the MS1 phase pinned one core. That is invisible at default thresholds because MS1 is only ~4% of frames, but it makes the MS1 SENSITIVITY sweep impossible: lowering trace:noise_threshold_int / ms1_chrom_peak_snr explodes MassTraceDetection's intensity-sorted apex list, and MassTraceDetection is a SERIAL greedy loop, so runtime blows up superlinearly on one core (measured: 40 min with no progress at 1.5 of 128 cores, vs 24 min for a whole default run). 1 = off (old behaviour).", false);
    registerIntOption_("perf:trace_bands", "<n>", 12, "[parallel] Split each window's fragment m/z range into N bands traced CONCURRENTLY. The diaPASEF method fixes the window count (24 here), which caps window-level parallelism far below the core count; banding lifts the ceiling to threads. Traces cannot span more than trace:mass_error_ppm in m/z, so banding with a halo is exact. 0 = auto (windows x bands >= threads), 1 = off.", false);
    registerDoubleOption_("perf:mem_fraction", "<f>", 0.75, "[dyn-mem] Fraction of currently-FREE RAM (/proc/meminfo MemAvailable + our own RSS) the window loop may commit. Concurrency is re-decided at every window admission, so the run adapts to other jobs on a shared node. One window is always admitted regardless, so a single oversized window cannot deadlock.", false);
    registerIntOption_("perf:max_concurrent_windows", "<n>", 0, "UPPER bound on isolation windows processed concurrently (memory admission may run fewer). Peak RAM scales with windows IN FLIGHT (each holds its frames + traces + grid), not with the total window count, so lowering this trades wall time for RAM. 0 = unlimited (use all threads).", false);

    registerTOPPSubsection_("diag", "Diagnostics ('debug' is reserved by TOPP)");
    registerStringOption_("diag:dump_ms1_tsv", "<prefix>", "", "[ms1-funnel] Write MS1 traces and inferred precursors to <prefix>.traces.tsv / <prefix>.precursors.tsv so precursor loss can be attributed to a stage (no trace / no hypothesis / no spectrum) rather than inferred from the final count. Empty = off.", false, true);

    registerIntOption_("max_charge", "<n>", 5, "Maximum precursor charge considered during isotope inference.", false);

    // NOTE: the subsection MUST be registered or printUsage_() throws ElementNotFound — which turns
    // any user parameter typo into a FATAL crash instead of a usage message. [bugfix]
    registerTOPPSubsection_("charge", "Isotope / charge-state inference");
    registerIntOption_("charge:min_charge", "<n>", 2, "Minimum precursor charge to emit; hypotheses below it "
                       "are dropped before extraction. Default 2, because tryptic peptides are essentially never singly "
                       "charged in ESI: z=1 hypotheses were ~30%% of emission and ~1.7%% of identified peptides, and removing "
                       "them raised peptide counts on both search engines on all three benchmark files while cutting runtime "
                       "17%%. Use 1 to emit them. 3 or more removes charge 2 and most of the result. Measurements: "
                       "docs/ISOTOPE-DUPLICATION-2026-09-03.md.", false);
    // [route-1] Only 25.7%% of the reference implementation's precursors are matched by us with the RIGHT m/z AND charge;
    // 17.2%% have the wrong charge and 22.9%% (decoy-corrected) the wrong monoisotope. A wrong charge or
    // mono gives a wrong NEUTRAL MASS, so the peptide can never be identified no matter how good the
    // fragments are. Partner-COUNTING cannot tell a real envelope from a coincidence; averagine shape
    // x isotope co-elution can.
    registerStringOption_("charge:scoring", "<mode>", "count", "Charge/monoisotope inference. 'count' is the partner-count heuristic (ties favour LOW charge, break-on-first-miss) and is the DEFAULT: measured against a DIA-NN reference on dataset B it gives 71.8% charge agreement vs 49.7% for 'envelope', and it eliminates the catastrophic 4->2 confusion (10,199 cases) that fabricates a ~-1198 Da phantom modification at m/z 600 in open search. 'envelope' scores each (charge, mono) hypothesis by averagine cosine x isotope co-elution and allows gapped envelopes; it was the default until 2026-07-21 and is kept for comparison. NOTE both are close to trivial: always answering z=2 scores 69.6%.", false);
    setValidStrings_("charge:scoring", {"envelope", "count"});
    registerDoubleOption_("charge:ambiguity_margin", "<x>", 0.0, "If a different-charge hypothesis scores within this margin of the best, emit it too (the reference implementation retains multiple charges when ambiguous). 0 = single charge per precursor. ONLY EFFECTIVE with charge:scoring=envelope -- the default 'count' mode returns before this is consulted, so the value is ignored there.", false);
    registerFlag_("assembly:dedup_precursors", "[route-4] Collapse precursor hypotheses that are the same species (same charge, m/z within ppm, co-eluting RT, co-located IM). Measured redundancy: 4.6 identified spectra per unique peptide vs the reference implementation's 1.46.");
    registerStringOption_("assembly:require_isotope_support", "<true/false>", "true",
                          "[emission] Drop precursor hypotheses with NO isotope partner -- the 'guessed' "
                          "singletons. ON BY DEFAULT since 2026-09-04, when the shipped default (off) was "
                          "measured for the first time and found to be 6.9x SLOWER and worse on both engines: "
                          "dataset D off = 1,255,577 spectra, 36:25 wall, 6.0x occupancy, Sage 12,212 / MSFragger "
                          "12,148; on = 655,776 spectra, 5:15, 65.7x, Sage 12,482 / MSFragger 12,337. An "
                          "earlier measurement (2026-07-26, at 10,036 peptides, before the charge gate, the "
                          "apex estimator and the integer detector) had put the cost at -0.78% peptides; the "
                          "sign has since flipped. Set false to keep the guessed singletons.", false);
    setValidStrings_("assembly:require_isotope_support", {"true", "false"});
    registerDoubleOption_("charge:mono_averagine_guard", "<slack>", 0.0, "[mono-guard] Reject a leftward isotope step (the walk that defines the MONOISOTOPE in charge:scoring=count) when the candidate peak is too weak to be a real lighter isotope under averagine. The partner search matches m/z, RT and IM but NEVER intensity, so it can latch onto a noise peak or a co-eluting species one isotope spacing below the true monoisotope -- making the reported precursor mass ~1.00335 Da TOO LIGHT. Closed search hides this (the engine enumerates isotope_errors and corrects it for free); OPEN/blind search turns it into a phantom -1.003 Da modification, and it is our single largest open-search delta-mass artefact on dataset D (8.20%% isotope-shifted PSMs vs the reference implementation's 3.46%%). The bound is averagine's LOOSEST: stepping from isotope m+1 to m the lighter peak is >= (m+1)/lambda times the heavier one, minimised at m=0 -> 1/lambda (lambda = 0.000594*neutral mass), so a value of 1.0 enforces exactly that minimum and never rejects a legitimate envelope; smaller values add slack for noise. 0 = off (old behaviour). [bench]", false);
    registerFlag_("charge:mono_averagine_select", "[mono-select] Choose the MONOISOTOPE by averagine envelope fit over the whole isotope run found by the charge:scoring=count walk, instead of taking the furthest-left peak reached. 'Furthest-left' is a STOPPING rule and can only trade over-reach (mono ~1.003 Da too LIGHT, a phantom -1 Da modification in open search) for under-reach (+1.003 Da) -- measured: an averagine STOPPING guard cut the -1 artefact 64%% but tripled +1, because a walk that halts early also leaves the M+1 peak unconsumed so it seeds its own too-heavy precursor. This is a SCORING decision instead: score every candidate mono position in the run by cosine against the averagine envelope and keep the best. Every found peak is still consumed, so de-isotoping and emission are unchanged -- only the reported precursor mass moves. [bench]");
    registerDoubleOption_("charge:iso_im_tolerance", "<1/K0>", 0.05, "IM tolerance for matching isotope partners. Default 0.05 (= 5x the fragment gate). A reviewer argued 0.01 (isotopes of one species share mobility, so 0.05 can admit cross-species partners), but 0.01 has NOT been measured on dataset D; change only behind the dataset D set gate.", false);

    // Bounds, not documentation: without these a negative count wraps through `Size` (min_fragments
    // -1 becomes ~1.8e19 and nothing is ever emitted) and a negative tolerance silently gates
    // everything out. Both fail as an empty-but-valid result, the worst failure mode here.
    // [adv-review codex/kimi 2026-09-03]
    setMinInt_("assembly:min_fragments", 1);
    setMinInt_("assembly:max_fragments", 1);
    setMinInt_("assembly:default_charge", 0);
    setMinInt_("max_charge", 1);
    setMinInt_("perf:trace_bands", 0);      // 0 = auto (documented)
    setMinFloat_("gate:delta_im", 0.0);
    setMinFloat_("gate:delta_rt", 0.0);
    setMinFloat_("gate:min_correlation", -1.0);
    setMaxFloat_("gate:min_correlation", 1.0);
    setMinFloat_("charge:iso_im_tolerance", 0.0);
    setMinFloat_("trace:ms1_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_noise_threshold_int", 0.0);
    setMinFloat_("trace:noise_threshold_int", 0.0);
    setMinInt_("perf:ms1_trace_bands", 0);      // 0 = auto (documented)
    setMinFloat_("perf:mem_fraction", 0.05);
    setMaxFloat_("perf:mem_fraction", 0.95);
    // deliberate: mass-defect filter default OFF (would silently drop PTM/nonspecific — the discovery cases). [H-7]
  }

  /// Configure an IM-aware MassTraceDetection and run it on @p map.
  /// MassTraceDetection wrapper. NOTE the apex gate inside MassTraceDetection is
  ///   peak_int > chrom_peak_snr * noise_threshold_int
  /// so the EFFECTIVE apex threshold is snr*noise, while trace MEMBERSHIP only needs > noise.
  /// Leaving snr at the OpenMS default 3.0 silently triples the fragment apex threshold.
  /// the reference implementation deliberately runs MS2 with "more relaxed constraints" than MS1 (Nat Commun 16:95
  /// Methods) — so MS1 and MS2 get INDEPENDENT snr / min-length here. [recall]
  /// `map` is taken by NON-const ref and CLEARED as soon as MassTraceDetection has consumed it:
  /// it is dead from that point on, but at ~20 B/peak (Peak1D 16 + IM float 4) it is ~1.1 GB per
  /// window that would otherwise stay resident through the valley-splitting peak. [compact]
  vector<Trace> detectTraces_(PeakMap& map, TraceStore& store, double delta_im, double noise,
                              double snr, double min_len, double min_sample_rate, double max_len,
                              double split_valleys, String* span_log = nullptr, int bands = 1)
  {
    MassTraceDetection mtd;
    mtd.setLogType(ProgressLogger::NONE); // quiet + thread-safe (called from the parallel window loop)
    Param p = mtd.getParameters();
    p.setValue("mass_error_ppm", getDoubleOption_("trace:mass_error_ppm"));
    p.setValue("noise_threshold_int", noise);
    p.setValue("chrom_peak_snr", snr);
    p.setValue("ion_mobility_tolerance", delta_im);        // IM-aware tracing [C-1]
    p.setValue("min_trace_length", min_len);
    // min_sample_rate (OpenMS default 0.5) requires a peak in >=50% of visited scans AND terminates
    // trace extension when occupancy drops below it — an intensity-INDEPENDENT deletion of gappy but
    // real fragment traces in sparse diaPASEF data. Relaxing it for MS2 is the sole surviving MS2
    // trace-validity lever. [recall]
    if (min_sample_rate >= 0.0) p.setValue("min_sample_rate", min_sample_rate);
    if (max_len > 0.0) p.setValue("max_trace_length", max_len);
    mtd.setParameters(p);

    // [epd-parallel] Reusable valley-split step. Called from a balanced parallel-for below, so 12
    // chunks split concurrently instead of once serially. Lock-free now that OpenMS EPD is patched
    // (patches/openms-epd-lockfree.patch); 96 concurrent callers (8 windows x 12 chunks) do not
    // contend. Output is order-independent downstream: MS2 traces are canonically re-sorted right
    // after detectTraces_, and ms1_traces are canonically sorted before inferPrecursors_.
    auto do_split = [&](vector<MassTrace>& in) -> vector<MassTrace> {
      if (split_valleys <= 0.0 || in.empty()) return std::move(in);
      ElutionPeakDetection epd;
      epd.setLogType(ProgressLogger::NONE);
      Param ep = epd.getParameters();
      ep.setValue("chrom_fwhm", split_valleys);
      ep.setValue("width_filtering", "off");
      ep.setValue("masstrace_snr_filtering", "false");
      epd.setParameters(ep);
      vector<MassTrace> out; out.reserve(in.size());
      epd.detectPeaks(in, out);
      return out;
    };

    vector<MassTrace> mts;
    if (bands > 1)
    {
      // [bands] The diaPASEF method fixes the window count at 24, so parallelising over
      // windows caps us at 24 threads however many cores exist. Fragment traces at
      // different m/z are INDEPENDENT, so a window's peaks can be partitioned by m/z and
      // traced concurrently. A trace can never span more than mass_error_ppm in m/z, so a
      // halo of 20x that tolerance makes the partition exact; each trace is kept only by the
      // band whose CORE contains its centroid, so halo duplicates are dropped rather than merged.
      const double ppm = getDoubleOption_("trace:mass_error_ppm");
      // Band edges from m/z QUANTILES, not uniform m/z: peak density varies by orders of
      // magnitude across a window, and uniform edges would leave one band holding most of the work.
      vector<double> sample;
      for (const auto& s : map)
        for (Size i = 0; i < s.size(); i += 97) sample.push_back(s[i].getMZ());
      sort(sample.begin(), sample.end());
      vector<double> edge(bands + 1);
      if (sample.empty()) { edge.assign(bands + 1, 0.0); }
      else
        for (int b = 0; b <= bands; ++b)
          edge[b] = sample[std::min(sample.size() - 1, (size_t)((double)b / bands * sample.size()))];
      edge.front() = 0.0; edge.back() = std::numeric_limits<double>::max();
      vector<double> ().swap(sample);

      // Partition ONCE (peaks are copied into exactly one core band plus any halo it falls in),
      // then release the source map so the halo duplicates are the only overhead.
      vector<PeakMap> sub(bands);
      for (auto& s : map)          // non-const: each spectrum is released once distributed
      {
        const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
        for (int b = 0; b < bands; ++b)
        {
          // Size each halo from ITS OWN edge: one h derived from edge[b] under-sizes the upper
          // side by the ppm growth across the band. [adv-review 2026-09-03]
          const double lo = edge[b] - edge[b] * ppm * 1e-6 * 20.0;
          const double hi = edge[b + 1] + edge[b + 1] * ppm * 1e-6 * 20.0;
          MSSpectrum t; t.setRT(s.getRT()); t.setMSLevel(1);
          OpenMS::DataArrays::FloatDataArray ia;
          ia.setName(Constants::UserParam::ION_MOBILITY);
          // [R3, 2026-09-02] spectra arrive m/z-sorted (loader sortByPosition), so the band's peak range is
          // two lower_bounds instead of a scan of every peak for every band (O(bands x peaks) -> O(bands log
          // peaks)). Same peaks in the same order -> byte-identical; the scan stays as the unsorted fallback.
          Size i0 = 0, i1 = s.size();
          if (s.isSorted())
          {
            i0 = (Size)(s.MZBegin(lo) - s.begin());
            i1 = (Size)(s.MZBegin(hi) - s.begin());
            t.reserve(i1 - i0); ia.reserve(i1 - i0);
            for (Size i = i0; i < i1; ++i) { t.push_back(s[i]); if (ima && i < ima->size()) ia.push_back((*ima)[i]); }
          }
          else
            for (Size i = 0; i < s.size(); ++i)
              if (s[i].getMZ() >= lo && s[i].getMZ() < hi)
              { t.push_back(s[i]); if (ima && i < ima->size()) ia.push_back((*ima)[i]); }
          if (!t.empty()) { t.getFloatDataArrays().push_back(std::move(ia)); sub[b].addSpectrum(std::move(t)); }
        }
        s.clear(true);   // distributed: release it now, so the source and the copies do not both
      }                  // stand at full size (the transient was ~2x the map)
      map.clear(true);
      vector<vector<MassTrace>> per(bands);
      #pragma omp taskloop grainsize(1) default(shared)
      for (int b = 0; b < bands; ++b)
      {
        MassTraceDetection m2; m2.setLogType(ProgressLogger::NONE); m2.setParameters(p);
        vector<MassTrace> t;
        sub[b].sortSpectra();
        if (sub[b].size() < 3) { sub[b].clear(true); continue; }   // MTD aborts below 3 spectra
        m2.run(sub[b], t);
        sub[b].clear(true);
        // keep only traces whose centroid lies in this band's CORE -> no duplicates from the halo
        for (auto& mt : t)
        {
          const double c = mt.getCentroidMZ();
          if (c >= edge[b] && c < edge[b + 1]) per[b].push_back(std::move(mt));
        }
      }
      size_t tot = 0; for (const auto& v : per) tot += v.size();
      mts.reserve(tot);
      for (auto& v : per) { for (auto& t : v) mts.push_back(std::move(t)); vector<MassTrace>().swap(v); }
    }
    else
    {
      if (map.size() >= 3) mtd.run(map, mts);   // MTD aborts below 3 spectra; a tiny input is not an error
      map.clear(true); // [compact] dead from here; ~1.1 GB/window freed before the EPD peak
    }

    // [way-2] SPLIT MASS TRACES AT CHROMATOGRAPHIC VALLEYS.
    // MassTraceDetection never splits at a local minimum: it terminates only on consecutive-miss
    // outliers or sample_rate, and max_trace_length only DISCARDS. So two peptides eluting 25 s apart
    // within the m/z tolerance merge into ONE multi-modal trace, whose XIC then correlates with
    // everything under it. ElutionPeakDetection is the OpenMS stage that actually smooths and splits
    // at local minima - we never ran it. width_filtering is forced to "off" because "fixed" is a
    // DELETION filter ([min_fwhm,max_fwhm]) and blanket deletion is on the falsified list.
    // [epd-parallel] Valley-splitting is the dominant serial phase (EPD-off dataset D 20:46 vs 28:40 for
    // MS2 alone). It is per-trace independent, so run it as a BALANCED parallel-for over the merged
    // traces -- more chunks than threads + dynamic schedule so no chunk's long traces stall the
    // rest. Lock-free (patched OpenMS), so the 8 concurrent windows above never contend. Trace SET
    // is identical to one serial split (per-trace, deterministic); only order differs, and both
    // consumers are order-independent (MS2 re-sorted after detectTraces_, MS1 sorted before
    // inferPrecursors_). Reverted once for a -0.8% "regression" that was actually the ms1_traces
    // determinism bug, now fixed.
    if (split_valleys > 0.0 && mts.size() > 1)
    {
      const int nchunk = std::max(1, bands * 4);      // 4x chunks vs threads -> dynamic load balance
      const size_t n_in = mts.size();
      vector<vector<MassTrace>> parts(nchunk);
      #pragma omp taskloop grainsize(1) default(shared)
      for (int ci = 0; ci < nchunk; ++ci)
      {
        const size_t lo = (size_t)ci * n_in / nchunk, hi = (size_t)(ci + 1) * n_in / nchunk;
        if (lo >= hi) continue;
        vector<MassTrace> chunk(std::make_move_iterator(mts.begin() + lo),
                                std::make_move_iterator(mts.begin() + hi));
        parts[ci] = do_split(chunk);
      }
      size_t tot = 0; for (auto& pth : parts) tot += pth.size();
      vector<MassTrace> out; out.reserve(tot);
      for (auto& pth : parts) { for (auto& t : pth) out.push_back(std::move(t)); vector<MassTrace>().swap(pth); }
      if (span_log) *span_log += " | valley-split " + String(n_in) + " -> " + String(out.size())
                               + " traces (" + String((int)(1000.0 * out.size() / std::max<size_t>(n_in,1)) / 10.0) + "%)";
      mts.swap(out);
    }
    // [merged-trace] instrumentation: is trace merging actually happening? Log the RT-span
    // distribution. Peaks are 5-30 s here, so a large p90/max means multi-modal merged traces.
    if (!mts.empty())
    {
      vector<double> spans;
      spans.reserve(mts.size());
      for (const auto& mt : mts)
      {
        double lo = 1e30, hi = -1e30;
        for (Size i = 0; i < mt.getSize(); ++i) { double r = mt[i].getRT(); lo = min(lo, r); hi = max(hi, r); }
        spans.push_back(hi - lo);
      }
      sort(spans.begin(), spans.end());
      size_t n = spans.size();
      size_t over45 = spans.end() - lower_bound(spans.begin(), spans.end(), 45.0);
      // the tail fractions the trace:max_span_sec cap is priced against
      size_t over120 = spans.end() - lower_bound(spans.begin(), spans.end(), 120.0);
      size_t over240 = spans.end() - lower_bound(spans.begin(), spans.end(), 240.0);
      if (span_log) *span_log = "traces=" + String(n) + " span_s median=" + String(spans[n/2])
                      + " p90=" + String(spans[(size_t)(0.9*n)]) + " p99=" + String(spans[(size_t)(0.99*n)])
                      + " max=" + String(spans.back())
                      + " frac>45s=" + String((int)(1000.0*over45/n)/10.0) + "%"
                      + " frac>120s=" + String((int)(10000.0*over120/n)/100.0) + "%"
                      + " frac>240s=" + String((int)(10000.0*over240/n)/100.0) + "%";
    }
    vector<Trace> out;
    out.reserve(mts.size());
    for (const auto& mt : mts) out.push_back(toTrace(mt, store));
    return out;
  }

  /// [route-1] Averagine isotope envelope, Poisson approximation: for neutral mass M the expected
  /// number of 13C is lambda ~ 0.000594*M, giving relative intensities [1, l, l^2/2, l^3/6, ...].
  /// This is what lets us pick the MONOISOTOPE: the mono is NOT always the most intense peak
  /// (above ~1800 Da M+1 exceeds M), so intensity-ranking alone mis-assigns it.
  static void averagineRatios(double neutral_mass, int n, vector<double>& out)
  {
    const double lam = std::max(0.10, 0.000594 * neutral_mass);
    out.assign(n, 0.0);
    double term = 1.0;
    for (int k = 0; k < n; ++k) { out[k] = term; term *= lam / (double)(k + 1); }
    double mx = 0.0; for (double v : out) mx = std::max(mx, v);
    if (mx > 0) for (double& v : out) v /= mx;
  }

  static double cosineSim(const vector<double>& a, const vector<double>& b)
  {
    double dot = 0, na = 0, nb = 0;
    for (size_t i = 0; i < a.size() && i < b.size(); ++i) { dot += a[i]*b[i]; na += a[i]*a[i]; nb += b[i]*b[i]; }
    return (na > 0 && nb > 0) ? dot / (sqrt(na) * sqrt(nb)) : 0.0;
  }

  /// [route-1] Pearson of two XICs over their OVERLAPPING RT points, properly mean-centred on the
  /// overlap (NOT zero-padded over a global grid). Isotopes of one peptide must co-elute; a
  /// coincidental peak at the right m/z spacing usually does not.
  /// Pearson correlation of two elution profiles over their shared frames. On a global axis the
  /// intersection is an integer merge -- no RT comparison, no tolerance, no binary search.
  static double xicCorr(const Trace& a, const Trace& b, double /*rt_tol*/)
  {
    // both are MS1 traces on the same frame table: the shared frames are the overlap of the two
    // spans, and a frame counts only where BOTH have a real point -- the old index intersection
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0; int m = 0;
    const uint32_t lo = std::max(a.frame0, b.frame0);
    const uint32_t hi = std::min(a.frame0 + a.len, b.frame0 + b.len);
    for (uint32_t f = lo; f < hi; ++f)
    {
      const double x = a.xv(f - a.frame0), y = b.xv(f - b.frame0);
      if (x > 0.0 && y > 0.0) { sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; ++m; }
    }
    if (m < 3) return -2.0;
    const double cov = sxy - sx * sy / m, vx = sxx - sx * sx / m, vy = syy - sy * sy / m;
    return (vx > 0 && vy > 0) ? cov / std::sqrt(vx * vy) : -2.0;
  }

  /// [route-1] Isotope/charge inference on MS1 traces. Scores each (charge, monoisotope) hypothesis
  /// by averagine-shape agreement x isotope-XIC co-elution, instead of counting partners.
  vector<Precursor_> inferPrecursors_(const vector<Trace>& ms1, int max_charge,
                                      double delta_rt, double iso_im_tol, double mass_ppm,
                                      bool envelope_scoring, double ambig_margin)
  {
    const double ISO = 1.0033548;
    const double PROTON_ = 1.00727646;
    const size_t N = ms1.size();

    // m/z-sorted index for O(log N) isotope-partner lookup (avoids O(N^2) over all MS1 traces). [H-8]
    vector<size_t> by_mz(N);
    for (size_t i = 0; i < N; ++i) by_mz[i] = i;
    sort(by_mz.begin(), by_mz.end(), [&](size_t a, size_t b) { return ms1[a].mz < ms1[b].mz; });
    vector<double> sorted_mz(N);
    for (size_t i = 0; i < N; ++i) sorted_mz[i] = ms1[by_mz[i]].mz;

    vector<size_t> order(N);
    for (size_t i = 0; i < N; ++i) order[i] = i;
    sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      if (ms1[a].intensity != ms1[b].intensity) return ms1[a].intensity > ms1[b].intensity;
      if (ms1[a].mz != ms1[b].mz) return ms1[a].mz < ms1[b].mz;
      return a < b; // total order so greedy isotope assignment is deterministic [code-review]
    });

    vector<bool> used(N, false);
    vector<Precursor_> out;
    for (size_t seed : order)
    {
      if (used[seed]) continue;
      const Trace& s = ms1[seed];
      const double tol = s.mz * mass_ppm * 1e-6;

      // Find an unused isotope partner near `target` m/z, co-localized in RT and IM (binary-search
      // the m/z index, scan only the ppm-window candidates); -1 if none. [H-8]
      auto findPartner = [&](double target) -> long {
        size_t klo = lower_bound(sorted_mz.begin(), sorted_mz.end(), target - tol) - sorted_mz.begin();
        size_t khi = upper_bound(sorted_mz.begin(), sorted_mz.end(), target + tol) - sorted_mz.begin();
        for (size_t k = klo; k < khi; ++k)
        {
          size_t j = by_mz[k];
          if (used[j] || j == seed) continue;
          if (fabs(ms1[j].rt - s.rt) <= delta_rt && fabs(ms1[j].im - s.im) <= iso_im_tol)
            return (long)j;
        }
        return -1;
      };

      // [baseline-fix v2] When envelope scoring is OFF, run the ORIGINAL algorithm end-to-end, not
      // just the original scoring. The previous attempt reverted only the score while keeping the new
      // candidate generation (gapped walk + lead 0..3), producing a THIRD algorithm: 2.52M MS1 traces
      // -> 1.21M precursors, 967,956 spectra, 7,308 peptides (baseline: 1.76M -> 956k, 741k, 9,103).
      // The original is: contiguous walk with BREAK-ON-FIRST-MISS, seed is always the starting point
      // (no lead offsets), monoisotope = furthest-LEFT partner reached, pick z by partner COUNT with
      // ties favouring the LOWEST z.
      if (!envelope_scoring)
      {
        int best_z = 0, best_n = 1;
        double best_mono = s.mz;
        vector<size_t> best_partners;
        for (int z = 1; z <= max_charge; ++z)
        {
          vector<size_t> partners;
          double mono = s.mz;
          double prev_int = s.intensity;
          vector<long> lefts;
          for (int k = 1; k <= 5; ++k)          // lighter partners define the monoisotope
          {
            long j = findPartner(s.mz - k * ISO / z);
            if (j < 0) break;                   // BREAK on first miss (contiguous only)
            // [mono-guard] A LEFTWARD step asserts the new peak is a LIGHTER ISOTOPE of the same
            // peptide, but findPartner only checks m/z, RT and IM -- never intensity. So the walk
            // happily latches onto a noise peak or a co-eluting species' peak that happens to sit
            // one isotope spacing below, and the monoisotope ends up one isotope TOO LIGHT.
            // Closed search hides this (the engine enumerates isotope_errors [-1,2] and corrects it
            // for free); OPEN search cannot, so it becomes a phantom -1.003 Da "modification" --
            // measured as our single largest open-search delta-mass artefact on dataset D.
            // Averagine bounds how weak a genuine lighter isotope may be: stepping from isotope m+1
            // to m, theo[m]/theo[m+1] = (m+1)/lambda, minimised at m=0 -> 1/lambda. Requiring the
            // candidate to reach that (times a slack factor) rejects the latch without touching a
            // real envelope.
            if (mono_guard_ > 0.0)
            {
              const double neutral = (s.mz - k * ISO / z) * z - z * PROTON_;
              const double lam = std::max(0.10, 0.000594 * neutral);
              if (ms1[j].intensity < mono_guard_ * prev_int / lam) break;
            }
            lefts.push_back(j);
            partners.push_back((size_t)j); mono = ms1[j].mz; prev_int = ms1[j].intensity;
          }
          vector<long> rights;
          for (int k = 1; k <= 5; ++k)          // heavier partners add confidence
          {
            long j = findPartner(s.mz + k * ISO / z);
            if (j < 0) break;
            rights.push_back(j);
            partners.push_back((size_t)j);
          }
          // [mono-select] Choosing the monoisotope by "furthest-left peak reached" is a STOPPING
          // rule, and a stopping rule can only trade over-reach (-1.003 Da) for under-reach (+1.003
          // Da) -- measured: the guard above cut the -1 artefact 64% but tripled +1, because a walk
          // that stops early also leaves the M+1 peak unconsumed, so it later seeds its OWN
          // precursor one isotope too heavy. Make it a SCORING decision instead: walk the full run,
          // then pick which peak IS the monoisotope by averagine fit of the envelope starting
          // there. Corrects both directions, and every found peak is still marked used below, so
          // de-isotoping and emission are untouched (only the reported mono m/z changes).
          if (mono_select_ && !lefts.empty())
          {
            const int nL = (int)lefts.size();
            vector<double> runI;                                  // ascending m/z: lefts, seed, rights
            for (int i = nL - 1; i >= 0; --i) runI.push_back(ms1[lefts[i]].intensity);
            runI.push_back(s.intensity);
            for (long j : rights) runI.push_back(ms1[j].intensity);
            double best_cos = -1.0; int best_c = nL;              // default: seed is the mono
            vector<double> theo, obs;
            for (int c = 0; c <= nL; ++c)                         // mono cannot be right of the seed
            {
              const int len = std::min(4, (int)runI.size() - c);
              if (len < 2) break;                                 // need >=2 peaks to score a shape
              const double neutral = (s.mz - (double)(nL - c) * ISO / z) * z - z * PROTON_;
              averagineRatios(neutral, len, theo);
              obs.assign(runI.begin() + c, runI.begin() + c + len);
              double mx = 0.0; for (double v : obs) mx = std::max(mx, v);
              if (mx <= 0.0) continue;
              for (double& v : obs) v /= mx;
              // [adv-fix] The SAME trap the envelope scorer documents below: the cosine of a short
              // vector is nearly always high, and `len` SHRINKS as the candidate mono moves right,
              // so raw cosine systematically prefers a too-HEAVY mono. Measured (first cut, raw
              // cosine): -1 bin 7,616 -> 2,314 as intended, but +1 exploded 4,697 -> 21,591 and the
              // unmodified fraction fell 31.2% -> 26.4%, at flat emission -- i.e. precursors that
              // were CORRECT were pushed one isotope right. Weight by evidence exactly as the
              // envelope scorer does: (n-1)/3 capped at 1 (4 isotopes 1.00, 3 -> 0.67, 2 -> 0.33).
              const double cs = cosineSim(obs, theo) * std::min(1.0, (double)(len - 1) / 3.0);
              if (cs > best_cos) { best_cos = cs; best_c = c; }
            }
            mono = (best_c < nL) ? ms1[lefts[nL - 1 - best_c]].mz : s.mz;
          }
          const int n = 1 + (int)partners.size();
          if (n > best_n) { best_n = n; best_z = z; best_mono = mono; best_partners = partners; }
        }
        Precursor_ pc0;
        pc0.trace_idx = seed;
        pc0.rt = s.rt; pc0.im = s.im;
        pc0.mono_mz = best_mono;
        // Carry the bin of whichever trace the monoisotope call landed on, so the exported mass is
        // calibrated from a measurement rather than from a number that passed through the gates.
        pc0.mono_tof = ms1[seed].tof; pc0.mono_b = ms1[seed].b;
        for (size_t j : best_partners)
          if (ms1[j].mz == best_mono) { pc0.mono_tof = ms1[j].tof; pc0.mono_b = ms1[j].b; break; }
        pc0.charge = best_z;                    // 0 => unknown, default_charge applied later
        pc0.n_isotopes = best_n;                // [Q4] envelope confidence for the mono call
        out.push_back(pc0);
        used[seed] = true;
        for (size_t j : best_partners) used[j] = true;   // de-isotope
        continue;
      }

      struct Hyp { int z; double mono; double score; vector<size_t> partners; };
      vector<Hyp> hyps;
      vector<double> theo, obs(4), obsn(4);

      for (int z = 1; z <= max_charge; ++z)
      {
        // [route-1] GAPPED envelopes: no break-on-first-miss (a sub-threshold M+1 with a present
        // M+2 previously killed the whole walk). Try the seed AND up to 3 isotopes below it as the
        // monoisotope candidate; averagine decides which one really is the mono.
        for (int lead = 0; lead <= 3; ++lead)
        {
          const double mono_mz = s.mz - lead * ISO / z;
          long mono_idx = (lead == 0) ? (long)seed : findPartner(mono_mz);
          if (lead > 0 && mono_idx < 0) continue;
          const double neutral = mono_mz * z - z * PROTON_;
          if (neutral < 400.0 || neutral > 6000.0) continue;

          vector<size_t> members;
          vector<long> iso_idx(4, -1);
          obs.assign(4, 0.0);
          obs[0] = ms1[mono_idx].intensity;
          iso_idx[0] = mono_idx;
          if ((size_t)mono_idx != seed) members.push_back((size_t)mono_idx);
          int found = 1;
          for (int k = 1; k <= 3; ++k)
          {
            long j = findPartner(mono_mz + k * ISO / z);
            if (j >= 0) { obs[k] = ms1[j].intensity; iso_idx[k] = j; ++found;
                          if ((size_t)j != seed) members.push_back((size_t)j); }
          }
          if (found < 2) continue;   // need >=2 isotopes before a shape can be scored

          averagineRatios(neutral, 4, theo);
          double mx = 0.0; for (double v : obs) mx = std::max(mx, v);
          if (mx <= 0) continue;
          for (int k = 0; k < 4; ++k) obsn[k] = obs[k] / mx;
          const double cos_sim = cosineSim(obsn, theo);

          double csum = 0.0; int cn = 0;
          for (int k = 1; k <= 3; ++k)
            if (iso_idx[k] >= 0) { csum += xicCorr(ms1[iso_idx[0]], ms1[iso_idx[k]], delta_rt); ++cn; }
          const double mean_corr = cn ? std::max(0.0, csum / cn) : 0.0;

          // Composite: shape AND co-elution must BOTH hold, WEIGHTED BY EVIDENCE.
          // [adv-fix] The naive product cos*corr + 0.01*n is statistically indefensible: a 2-isotope
          // hypothesis (cos 0.99 x corr 0.95 = 0.9605) BEATS 4 genuine isotopes (0.95 x 0.85 = 0.8475),
          // because the cosine of a 2-vector is nearly always high and the +0.01*n bonus (0.02 between
          // 2 and 4 isotopes) sits far below correlation noise. That biases us toward exactly the
          // sparse, weakly-supported envelopes we are trying to eliminate.
          // Fix: multiply by an evidence weight (n-1)/3 capped at 1 => 4 isotopes 1.00, 3 -> 0.67,
          // 2 -> 0.33. The counterexample now scores 0.808 vs 0.313 - correctly ordered.
          const double evidence = std::min(1.0, (found - 1) / 3.0);
          hyps.push_back({z, mono_mz, cos_sim * mean_corr * evidence, members});
        }
      }

      Precursor_ pc;
      pc.trace_idx = seed;
      pc.rt = s.rt; pc.im = s.im;

      if (hyps.empty())
      {
        pc.mono_mz = s.mz; pc.charge = 0;   // no envelope evidence at all -> unknown charge
      pc.mono_tof = s.tof; pc.mono_b = s.b;
        out.push_back(pc);
        used[seed] = true;
        continue;
      }

      // [baseline-fix] When envelope scoring is OFF we must still DE-ISOTOPE, otherwise every isotope
      // peak seeds its own precursor and the spectrum count explodes (measured: 741k -> 1,209,295
      // spectra, 9,103 -> 7,830 peptides). The off-branch therefore falls back to the ORIGINAL
      // partner-COUNT criterion with the original low-z tie-break, and still marks members used.
      if (!envelope_scoring)
      {
        for (auto& h : hyps) h.score = (double)h.partners.size();  // count, not shape
        sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) {
          if (a.score != b.score) return a.score > b.score;
          if (a.z != b.z) return a.z < b.z;                        // original: ties favour LOW z
          return a.mono < b.mono;
        });
      }
      else
      sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.z != b.z) return a.z > b.z;   // [route-1] ties now favour HIGHER z (old code favoured low z)
        return a.mono < b.mono;
      });
      const Hyp& best = hyps[0];
      pc.mono_mz = best.mono; pc.charge = best.z;
      pc.mono_tof = s.tof; pc.mono_b = s.b;
      for (size_t j : best.partners)
        if (ms1[j].mz == best.mono) { pc.mono_tof = ms1[j].tof; pc.mono_b = ms1[j].b; break; }
      // partners + the seed = envelope size. `members` excludes the seed, so this is the whole
      // envelope only because the seed is always one of its peaks. [adv-review kimi 2026-09-03]
      pc.n_isotopes = 1 + (int)best.partners.size();   // [Q4] envelope confidence for the mono call
      out.push_back(pc);
      used[seed] = true;
      for (size_t j : best.partners) used[j] = true;

      // [route-3] Ambiguous charge -> retain a second hypothesis (the reference implementation: "when the charge state
      // cannot be confidently determined, multiple values are retained"). 0 disables.
      if (ambig_margin > 0.0)
      {
        for (size_t h = 1; h < hyps.size(); ++h)
        {
          if (hyps[h].z == best.z) continue;
          if (best.score - hyps[h].score > ambig_margin) break;
          Precursor_ alt = pc;
          alt.mono_mz = hyps[h].mono; alt.charge = hyps[h].z;
          alt.n_isotopes = 1 + (int)hyps[h].partners.size();   // [Q4]
          out.push_back(alt);
          break;                            // at most one alternative
        }
      }
    }
    return out;
  }

  /// [route-4] Collapse precursor hypotheses that are the SAME species (same charge, m/z within ppm,
  /// co-eluting in RT, co-located in IM). Measured redundancy is 4.6 identified spectra per unique
  /// peptide vs the reference implementation's 1.46, so duplicates cost spectra, RAM and FDR budget without adding
  /// peptides. Keeps the most intense representative. Deterministic (total-order sort).
  static size_t dedupPrecursors_(vector<Precursor_>& pcs, const vector<Trace>& ms1,
                                 double mass_ppm, double rt_tol, double im_tol)
  {
    if (pcs.empty()) return 0;
    vector<size_t> idx(pcs.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
      if (pcs[a].charge != pcs[b].charge) return pcs[a].charge < pcs[b].charge;
      if (pcs[a].mono_mz != pcs[b].mono_mz) return pcs[a].mono_mz < pcs[b].mono_mz;
      if (pcs[a].rt != pcs[b].rt) return pcs[a].rt < pcs[b].rt;
      return a < b;
    });
    vector<bool> drop(pcs.size(), false);
    for (size_t a = 0; a < idx.size(); ++a)
    {
      size_t i = idx[a];
      if (drop[i]) continue;
      const double tol = pcs[i].mono_mz * mass_ppm * 1e-6;
      for (size_t b = a + 1; b < idx.size(); ++b)
      {
        size_t j = idx[b];
        if (pcs[j].charge != pcs[i].charge) break;
        if (pcs[j].mono_mz - pcs[i].mono_mz > tol) break;
        if (drop[j]) continue;
        if (fabs(pcs[j].rt - pcs[i].rt) <= rt_tol && fabs(pcs[j].im - pcs[i].im) <= im_tol)
        {
          if (ms1[pcs[j].trace_idx].intensity > ms1[pcs[i].trace_idx].intensity) { drop[i] = true; break; }
          drop[j] = true;
        }
      }
    }
    vector<Precursor_> keep;
    keep.reserve(pcs.size());
    size_t removed = 0;
    for (size_t i = 0; i < pcs.size(); ++i) { if (drop[i]) ++removed; else keep.push_back(pcs[i]); }
    pcs.swap(keep);
    return removed;
  }

  /// The fragment RT gate field. 99.4% of scoring visits die on |rt - pc.rt| > delta_rt, and
  /// reading an 8-byte double to reject each one is the last large bandwidth term in the scorer:
  /// 867e9 visits x 8 B = 6.9 TB per dataset D run, against a measured ~150-226 GB/s. `q` is the same RT
  /// quantised into buckets NEVER NARROWER than delta_rt, so two RTs within delta_rt of each other
  /// can never be more than one bucket apart; the reject then reads 2 bytes instead of 8 and the
  /// exact double test runs only on what survives it. Identical candidates in identical order --
  /// the exact test is the original one, unmoved, and the quantised test only precedes it. [soa2]
  ///
  /// The reject admits +/-2 buckets where +/-1 is provably enough in exact arithmetic. The extra
  /// bucket is insurance against a floor() landing on the far side of a boundary through rounding,
  /// which would silently DROP a candidate; it costs a pass rate of ~0.9% instead of ~0.6%.
  struct FragRt
  {
    vector<double> rt;
    vector<uint16_t> q;
    double rt0 = 0.0, inv = 0.0;
    double operator[](size_t i) const { return rt[i]; }
    /// Signed at the low end and clamped at both: a precursor may sit outside the window's
    /// fragment RT range. The test is written `!(x > 0.0)` rather than `x <= 0.0` so that a NaN
    /// coordinate lands in bucket 0 instead of reaching (int)NaN, which is undefined behaviour.
    int bucketOf(double v) const
    { const double x = (v - rt0) * inv; return !(x > 0.0) ? 0 : (x >= 65535.0 ? 65535 : (int)x); }
    /// inv == 0 is the DISABLED state, not an error: every bucket is then 0, every dq is 0, and
    /// the exact test alone decides. That is the fallback for every degenerate parameter below,
    /// and it costs no branch in the scan.
    void build(const vector<Trace>& fr, double delta_rt)
    {
      const size_t n = fr.size();
      rt.resize(n); q.assign(n, 0);
      rt0 = 0.0; inv = 0.0;
      double lo = 0.0, hi = 0.0;
      for (size_t i = 0; i < n; ++i) { rt[i] = fr[i].rt; if (!i || rt[i] < lo) lo = rt[i]; if (!i || rt[i] > hi) hi = rt[i]; }
      // A denormal delta_rt makes 1/delta_rt infinite, and an RT range spanning +/-DBL_MAX makes
      // the span infinite; both then produce inf*0 = NaN inside bucketOf. Neither survives real
      // acquisition, but both are reachable and both were UB, so every one degrades to disabled.
      if (!(delta_rt > 0.0) || !std::isfinite(delta_rt)) return;
      if (!std::isfinite(lo) || !std::isfinite(hi)) return;
      const double span = hi - lo;
      if (!std::isfinite(span)) return;
      double iv = 1.0 / delta_rt;
      // Buckets are delta_rt wide, but widened if the window's RT span would overflow uint16.
      // Wider is always safe: it can only put more fragments in reach of the exact test.
      if (span > 0.0 && span * iv > 65000.0) iv = 65000.0 / span;
      if (!std::isfinite(iv) || !(iv > 0.0)) return;
      rt0 = lo; inv = iv;
      for (size_t i = 0; i < n; ++i) q[i] = (uint16_t)bucketOf(rt[i]);
    }
  };

  /// Resample the precursor XIC onto the window grid (once), then for every IM+RT+overlap-passing
  /// candidate fragment compute the Pearson correlation and, if >= min_corr, call emit(fi, corr).
  /// Shared by the per-precursor path and the competitive (best-precursor) assignment. [bench]
  template <class Emit>
  void scoreCandidates_(const Precursor_& pc, const vector<Trace>& frag_traces,
                        const vector<double>& frag_im, const FragRt& frag_rt,
                        const vector<Trace>& ms1_traces,
                        const FragStats& fg, double delta_im, double delta_rt, double min_corr,
                        int min_corr_pts, vector<float>& pdense, Emit&& emit) const
  {
    const double G = (double)fg.G;
    const Trace& p_tr = ms1_traces[pc.trace_idx];
    // [wavelet] Correlate against a SMOOTHED precursor profile, as the reference implementation does. Only the
    // INTENSITIES change: RT positions, touched grid points and the overlap count are untouched,
    // so this cannot silently redefine the min_correlation_points gate. Fragment XICs stay RAW
    // (the reference implementation smooths only the precursor side), and emitted peak intensities are unaffected
    // because they come from the fragment traces, not from here.
    const vector<double>* p_int = nullptr;
    if (wl_scale_ > 0.0)
    {
      const int lv = atrousLevels(p_tr, wl_scale_);
      if (lv > 0)
      {
        thread_local vector<double> raw_i, sm_i, raw_rt;   // window loop is parallel; scratch is per-thread
        packReal(p_tr, raw_rt, raw_i);               // the smoother runs over the REAL points, packed
        atrousSmooth(raw_i, lv, sm_i);
        p_int = &sm_i;
      }
    }
    vector<int> touched;
    touched.reserve(p_tr.np());
    size_t pi = 0;                                    // index among REAL points (matches the packed smoother)
    for (size_t k = 0; k < p_tr.span(); ++k)
    {
      if (!p_tr.real(k)) continue;
      const double pval = p_int ? (*p_int)[pi] : (double)p_tr.xv(k);
      ++pi;
      const int gi = fg.nearest_local[p_tr.frame0 + k];   // MS1 frame -> window frame, precomputed
      if (gi >= 0) { if (pdense[gi] == 0.0f) touched.push_back(gi); pdense[gi] += (float)pval; }
    }
    double psum = 0, psumsq = 0, plogsum = 0;
    // plogsum normalises the precursor profile for gate:coelution=logoverlap. Constant across
    // every fragment of this precursor, so hoist it out of the inner loop.
    for (int gi : touched)
    {
      double v = pdense[gi];
      psum += v; psumsq += v * v;
      if (log_overlap_) plogsum += std::log1p(v);
    }
    double pmean = psum / G, pvar = psumsq - G * pmean * pmean;
    double pinv = pvar > 0 ? 1.0 / sqrt(pvar) : 0.0;
    if (pinv > 0.0)
    {
      size_t lo = lower_bound(frag_im.begin(), frag_im.end(), pc.im - delta_im) - frag_im.begin();
      size_t hi = upper_bound(frag_im.begin(), frag_im.end(), pc.im + delta_im) - frag_im.begin();
      long long n_rtpass = 0;
      stat_slice_ += (long long)(hi - lo); ++stat_prec_;
      const int pq = frag_rt.bucketOf(pc.rt);
      for (size_t fi = lo; fi < hi; ++fi)
      {
        const int dq = (int)frag_rt.q[fi] - pq;
        if (dq < -2 || dq > 2) continue;                                // 2 B reject, see FragRt
        // Read the gate field from a parallel array, and do NOT touch the 96-byte trace record
        // until the fragment has passed: 99.4% of visits end here, and at 8 B instead of 96 B the
        // scan moves ~12x fewer bytes. [soa]
        if (fabs(frag_rt[fi] - pc.rt) > delta_rt) continue;            // RT gate
        ++n_rtpass;
        if (fg.invnorm[fi] == 0.0) continue;                           // degenerate fragment
        const Trace& f = frag_traces[fi];
        if (fabs(f.mz - pc.mono_mz) < 0.01) continue;                  // exclude precursor peak [M-7]
        // [C9 A/B, 2026-09-02] 28% of emitted spectra carry the precursor's M+1 and 25% its M+2 as
        // "fragments" (unfragmented survivors). SPEXTRACTOR_DROP_PREC_ISO=1 drops M+1..M+3 (z known) too.
        if (drop_prec_iso_ && pc.charge > 0)
        {
          bool iso = false;
          for (int k = 1; k <= 3 && !iso; ++k) iso = fabs(f.mz - (pc.mono_mz + k * 1.00335 / pc.charge)) < 0.01;
          if (iso) continue;
        }
        double dot = 0; int overlap = 0;
        // [R2-late, 2026-09-02] overlap <= support, so a fragment with fewer support points than
        // min_corr_pts can never pass the overlap guard below: skip it before the dot product.
        // Output-identical by construction (the guard's continue is the only thing it would reach).
        if (f.np() < (size_t)std::max(min_corr_pts, 0)) continue;
        // two sequential streams: the fragment's span in the arena and pdense over window frames
        for (size_t k = 0; k < f.span(); ++k)
        { const float fv = f.xv(k); if (fv == 0.0f) continue;
          const float pv = pdense[f.frame0 + k]; if (pv != 0.0f) { dot += (double)pv * fv; ++overlap; } }
        if (overlap < min_corr_pts) continue;                          // overlap guard [H-4]
        // [B0] Pearson takes variance over the FULL grid G, so a fragment present in 5 of ~1300
        // grid points has its variance computed against ~1300 implicit zeros -- the zero-padding
        // defect the adversarial review flagged. It rewards sparse fragments that happen to spike
        // where the precursor spikes, and its scale depends on G rather than on the evidence.
        // AlphaDIA instead collapses fragments with sum(log(intensity+1)): a likelihood under
        // independence, where one missing point does not zero the score (unlike a product) and
        // one huge point cannot dominate (unlike a sum). Adapted to a (precursor, fragment) pair:
        // profile overlap in log space, summed ONLY where the precursor has signal, normalised to
        // [0,1]. No zero padding, no dependence on G.
        double c;
        if (var_support_ && !log_overlap_)                             // [Q1] pearson-only, per the flag doc
        {
          // [Q1] Pearson over the UNION SUPPORT (points where precursor OR fragment is non-zero),
          // not the full grid G. Removes the [B0] zero-padding: no ~1300 implicit zeros inflating
          // sparse fragments, no dependence on grid size. n = |touched U f.support|. fsum/fsumsq
          // over f.support are accumulated here (support is tiny, so a second pass is negligible).
          double fsum = 0.0, fsumsq = 0.0;
          for (size_t k = 0; k < f.span(); ++k)
          { const double fv = f.xv(k); if (fv > 0.0) { fsum += fv; fsumsq += fv * fv; } }
          double n = (double)(touched.size() + f.np() - (size_t)overlap);
          double pm = psum / n, fm = fsum / n;
          double pv2 = psumsq - n * pm * pm, fv2 = fsumsq - n * fm * fm;
          if (pv2 <= 0.0 || fv2 <= 0.0) continue;                      // degenerate on union support
          c = (dot - n * pm * fm) / std::sqrt(pv2 * fv2);
        }
        else if (log_overlap_)
        {
          // Normalise BOTH log-profiles to unit sum, then intersect. The first attempt divided
          // the intersection by the precursor's total log-intensity alone, which made the score
          // bounded by the fragment/precursor INTENSITY RATIO: a genuine fragment at 1% of
          // precursor intensity scored ~0.01 however perfectly it co-eluted, and was gated out.
          // Measured -18.9% peptides at gate 0.5 and -9.1% at 0.3 -- a magnitude statistic
          // standing in for a shape statistic. Unit-sum normalisation makes this a distributional
          // intersection: 1 iff the two normalised elution shapes agree, invariant to scale --
          // which is what Pearson gets right and the first version did not.
          double fsum = 0.0;
          for (size_t k = 0; k < f.span(); ++k)
          { const double fv = f.xv(k); if (fv > 0.0) fsum += std::log1p(fv); }
          if (fsum <= 0.0 || plogsum <= 0.0) continue;
          double inter = 0.0;
          for (size_t k = 0; k < f.span(); ++k)
          {
            const double fv = f.xv(k); if (fv <= 0.0) continue;
            inter += std::min(std::log1p(fv) / fsum,
                              std::log1p((double)pdense[f.frame0 + k]) / plogsum);
          }
          c = inter;
        }
        else c = (dot - G * pmean * fg.mean[fi]) * pinv * fg.invnorm[fi]; // Pearson
        if (c < min_corr) continue;                                    // correlation gate
        emit(fi, c);
      }
      stat_rtpass_ += n_rtpass;
    }
    for (int gi : touched) pdense[gi] = 0.0f;                          // reset scratch
  }

  /// [arena] Assertions on the SHIPPED compactUnreferenced(): the compaction must survive traces
  /// whose container order differs from their arena-offset order, which the previous walk did not.
  ExitCodes selftestArena_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (!ok) { writeLogError_(String("FAIL: ") + what); ++fail; }
      else writeLogInfo_(String("  ok  ") + what);
    };
    // four spans A(4) B(3) C(5) D(2) at offsets 0, 4, 7, 12, every value distinct; span k of trace T
    // holds T*100 + k + 1, with a gap (zero) inside D so interior zeros are exercised
    auto build = [](TraceStore& st, vector<Trace>& tr) {
      st = TraceStore(); tr.clear();
      const vector<pair<uint32_t, int>> spec = {{4, 1}, {3, 2}, {5, 3}, {3, 4}};
      for (const auto& [L, tag] : spec)
      {
        vector<pair<uint32_t, float>> pts;
        for (uint32_t k = 0; k < L; ++k) if (!(tag == 4 && k == 1)) pts.emplace_back(k, (float)(tag * 100 + k + 1));
        tr.push_back(makeSpan(pts, st));
      }
    };
    auto bytes = [](const TraceStore& st, const Trace& t) {
      return vector<float>(st.inten.begin() + t.off, st.inten.begin() + t.off + t.len);
    };
    TraceStore st0; vector<Trace> tr0; build(st0, tr0);
    vector<vector<float>> ref; for (const auto& t : tr0) ref.push_back(bytes(st0, t));
    chk(st0.inten.size() == 15 && tr0[3].len == 3 && tr0[3].npts == 2, "fixture: 15 floats, D spans 3 frames with a gap");

    // 1. container order D,B,A,C (offsets 12,4,0,7), keep D,B,A, drop C. The old container-order walk
    //    moved D to 0 and B to 2 before reading A from [0,4): A came back as D's and B's bytes.
    { vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]}; TraceStore s = st0;
      const Size freed = compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      chk(freed == 1, "1: one profile released");
      chk(s.inten.size() == 10, "1: arena holds exactly the kept spans (4+3+3)");
      chk(bytes(s, v[0]) == ref[3], "1: D survives (interior zero kept)");
      chk(bytes(s, v[1]) == ref[1], "1: B survives");
      chk(bytes(s, v[2]) == ref[0], "1: A survives -- the old walk overwrote it");
      chk(v[3].len == 0 && v[3].npts == 0, "1: C freed (span cleared)"); }
    // 2. offset order, everything needed: unchanged bytes, no move
    { vector<Trace> v = tr0; TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, true)) == 0, "2: nothing released");
      bool ok = s.inten.size() == 15; for (size_t i = 0; i < 4; ++i) ok = ok && bytes(s, v[i]) == ref[i];
      chk(ok, "2: offset-order input is a no-op"); }
    // 3. the freed trace first, then a trimmed span: the trim moved off/len/apex (as trimToSpan does)
    { vector<Trace> v = {tr0[2], tr0[0], tr0[1], tr0[3]}; TraceStore s = st0;
      v[1].off += 1; v[1].len -= 2; v[1].frame0 += 1; v[1].npts = 2;          // A trimmed to [1,3)
      const vector<float> a_trim(ref[0].begin() + 1, ref[0].begin() + 3);
      chk(compactUnreferenced(v, s, vector<bool>{false, true, true, true}) == 1, "3: C released first");
      chk(bytes(s, v[1]) == a_trim && bytes(s, v[2]) == ref[1] && bytes(s, v[3]) == ref[3], "3: trimmed A, B, D survive");
      chk(s.inten.size() == 2 + 3 + 3, "3: arena sized to the trimmed spans"); }
    // 4. nobody needed / a needed trace whose span was already released (len == 0) / duplicate references
    { vector<Trace> v = tr0; TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, false)) == 4 && s.inten.empty(), "4a: none needed -> empty arena");
      vector<Trace> w = tr0; TraceStore s2 = st0; w[0].freeProfile();
      chk(compactUnreferenced(w, s2, vector<bool>(4, true)) == 0 && s2.inten.size() == 11 && bytes(s2, w[1]) == ref[1], "4b: a needed len==0 trace owns no bytes"); }
    // 5. compaction is idempotent
    { vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]}; TraceStore s = st0;
      compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      const vector<float> once = s.inten;
      compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      chk(s.inten == once && bytes(s, v[2]) == ref[0], "5: second compaction is a no-op"); }
    // 6. overlapping spans are refused before any byte moves
    { vector<Trace> v = tr0; TraceStore s = st0; v[1].off = 2;                    // B now overlaps A
      bool threw = false; try { compactUnreferenced(v, s, vector<bool>(4, true)); } catch (const OpenMS::Exception::Precondition&) { threw = true; }
      chk(threw, "6: overlapping spans throw Precondition"); }

    writeLogInfo_(fail ? "[arena] SELFTEST FAILED" : "[arena] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /// [wavelet] Assertions on the SHIPPED a-trous smoother. Each check corresponds to a way the
  /// transform could be wrong while still "looking smoothed" in a plot.
  ExitCodes selftestWavelet_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (!ok) { writeLogError_(String("FAIL: ") + what); ++fail; }
      else writeLogInfo_(String("  ok  ") + what);
    };

    // 1. DC preservation. The B3 taps must sum to exactly 1: if they do not, smoothing RESCALES
    //    the precursor profile, and since Pearson is scale-invariant the error would be silent
    //    here but corrupt gate:coelution=logoverlap, which is NOT scale-invariant.
    {
      vector<double> flat(64, 7.5), out;
      atrousSmooth(flat, 4, out);
      double worst = 0.0;
      for (double v : out) worst = std::max(worst, fabs(v - 7.5));
      chk(worst < 1e-9, "DC preserved (constant signal unchanged)");
    }

    // 2. Noise reduction, and 3. peak retention -- the two must hold TOGETHER, else "smoothing"
    //    is either a no-op or a peak destroyer. Deterministic pseudo-noise (no RNG dependence).
    {
      const int n = 128;
      vector<double> sig(n), noisy(n), out;
      for (int i = 0; i < n; ++i)
      {
        const double x = (i - 64.0) / 8.0;               // Gaussian peak, sigma = 8 samples
        sig[i] = 1000.0 * std::exp(-0.5 * x * x);
        noisy[i] = sig[i] + 60.0 * std::sin(i * 2.7) * std::cos(i * 1.31);   // high-frequency
      }
      atrousSmooth(noisy, 2, out);
      double e_in = 0, e_out = 0;
      for (int i = 0; i < n; ++i) { e_in += fabs(noisy[i] - sig[i]); e_out += fabs(out[i] - sig[i]); }
      chk(e_out < 0.75 * e_in, "noise reduced (residual vs clean peak shrinks >25%)");
      chk(out[64] > 0.80 * sig[64], "peak height retained (>80% at apex)");
    }

    // 4. Mirror edges. Zero-padding would fabricate a downslope at a truncated peak -- the same
    //    bias the FWHM estimator avoids by excluding truncated traces.
    {
      vector<double> edge(32, 0.0), out;
      for (int i = 0; i < 32; ++i) edge[i] = 1000.0 - 25.0 * i;    // apex sits AT sample 0
      atrousSmooth(edge, 2, out);
      chk(out[0] > 0.85 * edge[0], "mirror boundary (edge apex not pulled toward zero)");
    }

    // 5. Level selection follows the SAMPLING, not a constant. Two traces of the same physical
    //    width but different cycle times must get different level counts.
    {
      // Profiles index the global axis, so the test builds an axis with the two spacings on it.
      vector<double>& ax = rtAxis();
      const vector<double> saved = ax;
      ax.clear();
      for (int i = 0; i < 40; ++i) ax.push_back(i * 0.5);      // dt = 0.5 s -> indices 0..39
      for (int i = 0; i < 40; ++i) ax.push_back(1000.0 + i * 2.0);   // dt = 2.0 s -> indices 40..79
      TraceStore ts;
      { vector<uint32_t> ri(80); for (int i = 0; i < 80; ++i) ri[i] = (uint32_t)i; ts.setFrames(ri, nullptr); }
      vector<pair<uint32_t, float>> pf, ps, pt;
      for (int i = 0; i < 40; ++i) pf.emplace_back((uint32_t)i, 100.0f);
      for (int i = 0; i < 40; ++i) ps.emplace_back((uint32_t)(40 + i), 100.0f);
      for (int i = 0; i < 3; ++i)  pt.emplace_back((uint32_t)i, 1.0f);
      Trace fast = makeSpan(pf, ts), slow = makeSpan(ps, ts), tiny = makeSpan(pt, ts);
      fast.st = &ts; slow.st = &ts; tiny.st = &ts;
      const int lf = atrousLevels(fast, 4.0), ls = atrousLevels(slow, 4.0);
      chk(lf > ls, "levels scale with sampling (finer dt -> more levels)");
      chk(atrousLevels(fast, 0.0) == 0, "scale 0 disables smoothing");
      chk(atrousLevels(tiny, 4.0) == 0, "too-short XIC -> 0 levels");
      ax = saved;
    }

    writeLogInfo_(fail ? "[wavelet] SELFTEST FAILED" : "[wavelet] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /// Build a pseudo-MS2 from a fragment list (m/z,intensity) + ranking scores: keeps top
  /// max_frags by score, drops if < min_frags, annotates the synthetic precursor. Deterministic.
  /// EXPORT-TIME CALIBRATION. On a flight-time axis the measured quantity is a bin; m/z is a
  /// calibrated function of it and of the frame's digitizer temperature. Applying that function
  /// ONCE, here, when the value is written, is what keeps the reported mass free of anything the
  /// intermediate gates did in cached m/z space -- and the apex frame's own factor is the most
  /// accurate one to use. A zero bin means the trace came from the OpenMS detector, which has no
  /// bin, and the value it computed stands.
  static double exportMz_(uint32_t tof, double b, double cached)
  {
    if (tof == 0 || !(b > 0.0) || !tofAxis().ok) return cached;
    const double v = tofAxis().mzOf(tof, b);
    return std::isfinite(v) && v > 0.0 ? v : cached;
  }

  void assembleFromList_(const Precursor_& pc, double win_lo, double win_hi,
                         vector<pair<double, double>>& frags, vector<double>& frag_scores,
                         Size min_frags, Size max_frags, MSSpectrum& out) const
  {
    if (frags.size() < min_frags) return;
    if (frags.size() > max_frags)
    {
      vector<size_t> idx(frags.size());
      for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
      // [rank] Rank by correlation (ours) or by intensity (the reference implementation's "top N highest intensity").
      // The engine re-ranks by INTENSITY downstream (Sage max_peaks=150), so a correlation-ranked
      // cap can discard high-intensity peaks the engine would have scored. Ties break on m/z then
      // index in BOTH modes so the output stays deterministic.
      const bool by_int = rank_by_intensity_;
      partial_sort(idx.begin(), idx.begin() + max_frags, idx.end(), [&](size_t a, size_t b) {
        const double ka = by_int ? frags[a].second : frag_scores[a];
        const double kb = by_int ? frags[b].second : frag_scores[b];
        if (ka != kb) return ka > kb;
        if (frags[a].first != frags[b].first) return frags[a].first < frags[b].first;
        return a < b;
      });
      vector<pair<double, double>> keep; keep.reserve(max_frags);
      for (size_t i = 0; i < max_frags; ++i) keep.push_back(frags[idx[i]]);
      frags.swap(keep);
    }
    out.setMSLevel(2);
    out.setType(SpectrumSettings::SpectrumType::CENTROID);
    out.setRT(pc.rt);
    for (const auto& fr : frags) out.emplace_back(fr.first, fr.second);
    out.sortByPosition();
    OpenMS::Precursor prec;
    prec.setMZ(exportMz_(pc.mono_tof, pc.mono_b, pc.mono_mz));   // calibrate on write
    if (pc.charge > 0) prec.setCharge(pc.charge);
    prec.setIsolationWindowLowerOffset(max(0.0, pc.mono_mz - win_lo));   // [C-5]
    prec.setIsolationWindowUpperOffset(max(0.0, win_hi - pc.mono_mz));
    prec.setDriftTime(pc.im);
    prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    out.setPrecursors({prec});
    // [provenance] Stamp the inference class onto the emitted spectrum so identified peptides can be
    // joined back to guessed (no isotope support) vs isotope-supported precursors -- the number that
    // sets the peptide cost of an isotope-support gate. Written as a <userParam> in the mzML.
    out.setMetaValue("spx_guessed", (int)(pc.guessed ? 1 : 0));
    // [Q4] Isotope-offset hypothesis for open/blind search: the number of isotope peaks behind the
    // mono call. 0/1 => the mono is weakly determined and could be off by +-k*1.00335 Da, so a
    // downstream delta-mass corrector should treat a near-integer-Da delta as an isotope miscall to
    // reassign rather than a real modification. This is the metadata assembly:open_search_safe's doc
    // promised but never wrote (only spx_guessed was stamped). Harmless in closed search (ignored).
    out.setMetaValue("spx_n_isotopes", pc.n_isotopes);
  }

  /// [E5, 2026-09-02] The ONE place the emitted-intensity weights live. Returns {cap-key intensity,
  /// emitted intensity}: the IM-proximity weight feeds both (IM-on-c, orthogonal to correlation); the
  /// [Q2] correlation-power weight feeds ONLY the emitted (engine-visible) intensity -- putting it in the
  /// cap key made frag_scores = base*c^(k+1) and double-hit borderline-c faint fragments at the 500-cap.
  /// Until today apportion and rp_max emitted RAW intensities, so every A/B of those flags silently
  /// compared corr_power=2 against corr_power=0.
  pair<double, double> weighted_(double base, double c, double dim) const
  {
    double inten = base;
    if (im_weight_sigma_ > 0.0) inten *= std::exp(-(dim * dim) / (2.0 * im_weight_sigma_ * im_weight_sigma_));
    double emit_inten = inten;
    if (corr_power_ > 0.0 && c > 0.0) emit_inten *= std::pow(c, corr_power_);
    return {inten, emit_inten};
  }

  /// Per-precursor assembly: this precursor claims every fragment passing its gate (a fragment may
  /// be shared across precursors -> chimeric). Leaves @p out empty on failure.
  void assembleOne_(const Precursor_& pc, double win_lo, double win_hi,
                    const vector<Trace>& frag_traces, const vector<double>& frag_im,
                    const FragRt& frag_rt,
                    const vector<Trace>& ms1_traces, const FragStats& fg,
                    double delta_im, double delta_rt, double min_corr, int min_corr_pts,
                    Size min_frags, Size max_frags, vector<float>& pdense, MSSpectrum& out) const
  {
    vector<pair<double, double>> frags;
    vector<double> frag_scores;
    scoreCandidates_(pc, frag_traces, frag_im, frag_rt, ms1_traces, fg, delta_im, delta_rt, min_corr,
                     min_corr_pts, pdense, [&](size_t fi, double c) {
      const auto [inten, emit_inten] = weighted_(frag_traces[fi].intensity, c, frag_im[fi] - pc.im);
      frags.emplace_back(exportMz_(frag_traces[fi].tof, frag_traces[fi].b, frag_traces[fi].mz), emit_inten);
      frag_scores.push_back(c * inten);
    });
    assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, out);
  }

  ExitCodes main_(int argc, const char** argv) override
  {
    phase_clock_();   // [perf-instr] fix the epoch HERE: it is a first-call static, and the
                      // streaming path first touched it only AFTER loading, so load time was
                      // invisible and every [t=..] was relative to post-load. [codex review]

    // -threads defaults to 1 in TOPPBase, a poor default for a tool whose window loop IS the
    // runtime. Default to every core. The command line is scanned rather than the value tested,
    // so that an EXPLICIT `-threads 1` still means one thread: testing `value == 1` alone would
    // silently override the one user who actually wants it serial.
    int n_threads_req = getIntOption_("threads");
    {
      bool given = false;
      for (int i = 1; i < argc && !given; ++i) given = (String(argv[i]) == "-threads");
      if (!given && n_threads_req == 1)
      {
        n_threads_req = std::max(1, (int)std::thread::hardware_concurrency());
        omp_set_num_threads(n_threads_req);
        writeLogInfo_("-threads not given: using all " + String(n_threads_req) + " cores. Pass "
                      "-threads <n> to limit it, which is what you want on a shared machine.");
      }
    }
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const double delta_im = getDoubleOption_("gate:delta_im");
    const double delta_rt = getDoubleOption_("gate:delta_rt");
    if (getFlag_("diag:selftest_wavelet")) return selftestWavelet_();
    if (getFlag_("diag:selftest_arena")) return selftestArena_();
    log_overlap_ = (getStringOption_("gate:coelution") == "logoverlap");
    rank_by_intensity_ = (getStringOption_("assembly:rank_by") == "intensity");
    im_weight_sigma_ = getDoubleOption_("assembly:im_weight_sigma");
    mono_guard_ = getDoubleOption_("charge:mono_averagine_guard");
    mono_select_ = getFlag_("charge:mono_averagine_select");
    var_support_ = getFlag_("gate:variance_support");
    corr_power_ = getDoubleOption_("assembly:corr_power");
    const double min_corr = getDoubleOption_("gate:min_correlation");
    const int min_corr_pts = getIntOption_("gate:min_correlation_points");
    const Size min_frags = (Size)getIntOption_("assembly:min_fragments");
    const Size max_frags = (Size)getIntOption_("assembly:max_fragments");
    const int max_charge = getIntOption_("max_charge");
    mzEstimator() = getStringOption_("trace:mz_estimator");   // before any toTrace() call
    const double mass_ppm = getDoubleOption_("trace:mass_error_ppm");

    PeakMap exp;
    // Accept mzML, mzPeak, or a Bruker .d directory; FileHandler auto-detects the format
    // (mzPeak -> MzPeakFile; .d -> BrukerTimsFile, requires WITH_OPENTIMS).
    PeakMap ms1_map;
    map<WinKey, vector<CompactFrame>> ms2_by_window;
    CompactStats cstat;
    auto winKey = [](double lo, double hi, uint32_t window_group) -> WinKey {
      return {(int)llround(lo * 100.0), (int)llround(hi * 100.0), (int)window_group};
    };
    const bool stream_load = (getStringOption_("perf:stream_load") == "true");
    bool streamed = false;

    if (stream_load)
    {
      // [stream] frame-by-frame: pick + compact on arrival, never hold the whole run.
      // Also the ONLY place the reader's NATIVE pre-centroid frame aggregation is reachable.
      BrukerTimsFile::Config cfg;
      cfg.dia_ms2_n_neighbors = getIntOption_("trace:native_ms2_neighbors");
      cfg.ms1_n_neighbors     = getIntOption_("trace:native_ms1_neighbors");
      PeakPickerIM spicker;
      Param sp = spicker.getParameters();
      sp.setValue("pickIMCluster:im_tolerance_cluster", getDoubleOption_("gate:delta_im"));
      sp.setValue("pickIMCluster:ppm_tolerance_cluster", getDoubleOption_("trace:mass_error_ppm"));
      if (const char* m = std::getenv("SPEXTRACTOR_PICK_MZ_MODE")) sp.setValue("pickIMCluster:mz_mode", std::string(m));   // [C2 pick-level A/B] weighted|seed|top3
      spicker.setParameters(sp);
      PickCompactConsumer consumer(spicker, ms2_by_window, ms1_map, cstat, winKey);
      try
      {
#ifdef SPEXTRACTOR_WITH_MZPEAK
        if (FileHandler::getTypeByFileName(in) == FileTypes::MZPEAK)
#else
        if (false)
#endif
        {
#ifdef SPEXTRACTOR_WITH_MZPEAK
          spx::loadMzPeakStreaming(in, consumer, n_threads_req);   // consumer counts frames itself
#else
          throw Exception::NotImplemented(__FILE__, __LINE__, "streaming .mzpeak input needs a build with -DMZPEAK_ROOT");
#endif
        }
        else BrukerTimsFile().loadDIAStreaming(in, consumer, cfg);
        streamed = true;
        phase_stats_()["LOAD(stream)"].wall = phase_clock_();      // [perf-instr] start-to-here
        phase_stats_()["LOAD(stream)"].cpu  = cpu_seconds_();
        phase_stats_()["LOAD(stream)"].rss_end_mb = rss_mb_();
        if (phase_stats_()["LOAD(stream)"].n++ == 0) phase_order_().push_back("LOAD(stream)");
        writeLogInfo_("[perf-load] pick(flush) wall=" + String(flush_wall_(), 1) + " s inside LOAD(stream)");
        if (std::getenv("SPEXTRACTOR_LOAD_ONLY")) { report_phases_(phase_clock_()); return EXECUTION_OK; }   // [perf-load] decode/hand-off/pick split only
        writeLogInfo_("[stream] frame-by-frame load done: " + String(consumer.frames_seen)
                      + " frames, MS1 " + String(ms1_map.size()) + ", windows "
                      + String(ms2_by_window.size()) + " (native ms2_neighbors="
                      + String(cfg.dia_ms2_n_neighbors) + ")" + clk_() + rss_() + mem_());
      }
      catch (const Exception::BaseException& e)
      {
        writeLogWarn_(String("[stream] streaming load failed (") + e.getName()
                      + "); falling back to loadExperiment. " + e.getMessage());
      }
    }

    if (!streamed)
    {
      Phase _ph("LOAD(full)");
      FileHandler().loadExperiment(in, exp,
#ifdef SPEXTRACTOR_WITH_MZPEAK
                                   {FileTypes::MZML, FileTypes::MZPEAK, FileTypes::BRUKER_TDF},
#else
                                   {FileTypes::MZML, FileTypes::BRUKER_TDF},
#endif
                                   log_type_);
    }

    //-------------------------------------------------------------
    // Input normalization / validation [C-4]
    //-------------------------------------------------------------
    // [stream] The streaming path fills ms1_map/ms2_by_window directly via an IM-aware consumer
    // (PickCompactConsumer calls ensureIMArrayName), leaving `exp` EMPTY -- so the exp-based IM check
    // only applies to the non-streaming load. Running it on the streamed (empty) `exp` wrongly aborted
    // with "no ion-mobility data at MS1", which is why native pre-centroiding aggregation could never
    // be benchmarked. Streamed data is guarded by the ms1_map.empty() check below instead.
    if (!streamed)
    {
      IMFormat imf = IMTypes::determineIMFormat(exp, 1);
      if (imf == IMFormat::NONE)
      {
        writeLogError_("Error: input has no ion-mobility data at MS1. SpeXtractor requires diaPASEF (IM DIA) input.");
        return ILLEGAL_PARAMETERS;
      }
    }
    // deliberate: assumes 1/K0 (VSSC). ms/CCS conversion and FAIMS rejection are Phase 2 [M-2].

    // Centroid IM frames that are still profile, keeping per-peak IM at a tolerance no coarser
    // than the trace IM tolerance [C-2]. Split MS1 vs MS2 into separate working maps.
    PeakPickerIM picker;
    Param pp = picker.getParameters();
    pp.setValue("pickIMCluster:im_tolerance_cluster", delta_im);
    pp.setValue("pickIMCluster:ppm_tolerance_cluster", mass_ppm);
    picker.setParameters(pp);

    // Peak-pick (centroid) all IM frames in PARALLEL. PeakPickerIM is firstprivate so each thread
    // has its own copy -> the mutable ccs_warning_shown_ flag is not shared (a shared const picker
    // would be a real data race). pickIMCluster mutates only its own spectrum; distinct indices are
    // independent. [par-Crit-1]
    writeLogInfo_("Loaded " + String(exp.size()) + " frames." + clk_() + rss_()); // raw residency peak? [mem]
    writeLogInfo_("Peak-picking " + String(exp.size()) + " frames (OpenMP)...");
    // Param is a deep-copy value type (mutable ParamNode root_, defaulted copy) so firstprivate
    // gives each thread a fully independent picker -> no refcount/COW race. pickIMCluster can throw
    // (e.g. unexpected IM format); an exception escaping an OMP region terminates the process, so
    // capture it and rethrow serially. [par-Crit-1, code-review]
    std::exception_ptr pick_err;
    #pragma omp parallel for firstprivate(picker) schedule(dynamic, 16)
    for (long i = 0; i < (long)exp.size(); ++i)
    {
      try
      {
        MSSpectrum& s = exp[i];
        if (!s.containsIMData()) continue;
        if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
        ensureIMArrayName(s); // so MassTraceDetection reads per-peak IM (IM-aware tracing) [C-1]
      }
      catch (...)
      {
        #pragma omp critical
        if (!pick_err) pick_err = std::current_exception();
      }
    }
    if (pick_err) std::rethrow_exception(pick_err);

    // Serial dispatch: MOVE (not copy) frames into the MS1 map / per-window MS2 maps to avoid
    // transient memory doubling; container inserts aren't thread-safe so this stays serial.
    if (!streamed)
    for (MSSpectrum& s : exp)
    {
      if (!s.containsIMData()) continue;
      if (s.getMSLevel() == 1)
      {
        ms1_map.addSpectrum(std::move(s));
      }
      else if (s.getMSLevel() == 2 && !s.getPrecursors().empty())
      {
        const OpenMS::Precursor& pr = s.getPrecursors()[0];
        double lo = pr.getMZ() - pr.getIsolationWindowLowerOffset();
        double hi = pr.getMZ() + pr.getIsolationWindowUpperOffset();
        ms2_by_window[winKey(lo, hi, windowGroupOf(s.getNativeID()))].push_back(compactify(s, cstat));
        s.clear(true);   // [mem] the picked frame is dead once compacted; waiting for exp.clear()
      }                  // below held the whole picked run (20 B/peak) beside the compact store
    }
    exp.clear(true); // frames moved out; release the container
    {
      size_t cb = 0, cp = 0;
      for (const auto& kv : ms2_by_window) for (const auto& f : kv.second) { cb += f.bytes(); cp += f.mzq.size(); }
      const size_t lost = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
      if (lost)
      {
        writeLogWarn_("[compact] DROPPED " + String(lost) + " peaks: no_im_array=" + String(cstat.no_im_array)
                      + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz)
                      + " bad_im=" + String(cstat.bad_im) + " (kept " + String(cstat.kept) + ") -- if this is"
                      + " nonzero the compact store is NOT lossless and results differ from the old path.");
      }
      size_t ms1_pk = 0; for (const auto& sp : ms1_map) ms1_pk += sp.size();
      writeLogInfo_("Split into MS1 + " + String(ms2_by_window.size()) + " windows (picked); compact store "
                    + String(cb / (1024ULL * 1024ULL)) + " MB for " + String(cp) + " peaks (~"
                    + String(cp ? (double)cb / cp : 0.0) + " B/peak); MS1 " + String(ms1_map.size()) + " frames / "
                    + String(ms1_pk) + " peaks (" + String(ms1_pk * 20 / (1024ULL * 1024ULL)) + " MB as PeakMap)." + clk_() + rss_() + mem_()); // [mem]
    }

    if (ms1_map.empty())
    {
      writeLogError_("Error: no MS1 ion-mobility frames found; cannot seed precursors.");
      return INCOMPATIBLE_INPUT_DATA;
    }

    // THE GLOBAL RT AXIS. Every frame contributes one retention time; every trace point sits on a
    // frame. Building it here, once, is what lets a trace point be 8 bytes (frame index +
    // intensity) instead of 16 (RT + intensity), and lets two profiles be aligned by integer
    // equality instead of a binary search per point. Frames are counted once even when an MS1 and
    // an MS2 frame share an RT, which is what makes the index a valid shared key.
    {
      vector<double>& ax = rtAxis();
      ax.clear();
      { size_t nf = ms1_map.size(); for (const auto& kv : ms2_by_window) nf += kv.second.size();
        ax.reserve(nf + 1024); }   // MS1 *and* every window's MS2 frames land here [mem]
      for (const auto& sp : ms1_map) ax.push_back(sp.getRT());
      for (const auto& kv : ms2_by_window) for (const auto& f : kv.second) ax.push_back(f.rt);
      sort(ax.begin(), ax.end());
      ax.erase(unique(ax.begin(), ax.end()), ax.end());
      writeLogInfo_("RT axis: " + String(ax.size()) + " distinct frame retention times"
                    + (ax.empty() ? "" : " spanning " + String(ax.front()) + "-" + String(ax.back()) + " s"));
      if (ax.empty())
      {
        writeLogError_("Error: no frame retention times; cannot build the RT axis.");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //-------------------------------------------------------------
    // MS1 traces + precursor/charge inference (shared, once) [H-8,H-10]
    //-------------------------------------------------------------
    const int agg_ms1 = getIntOption_("trace:frame_aggregation_ms1_n");
    ms1_map.sortSpectra();
    if (agg_ms1 > 1) { aggregateFrames_(ms1_map, agg_ms1, mass_ppm, delta_im); writeLogInfo_("MS1 cross-frame aggregation N=" + String(agg_ms1) + rss_()); }
    const double ms1_snr = getDoubleOption_("trace:ms1_chrom_peak_snr");
    const double ms2_snr = getDoubleOption_("trace:ms2_chrom_peak_snr");
    const double ms2_minlen = getDoubleOption_("trace:ms2_min_length_sec");
    const double ms2_msr = getDoubleOption_("trace:ms2_min_sample_rate");
    double ms2_split = getDoubleOption_("trace:ms2_split_valleys");
    const double max_trace_len = getDoubleOption_("trace:max_trace_length_sec");
    const double max_span = getDoubleOption_("trace:max_span_sec");
    const int agg_n = getIntOption_("trace:frame_aggregation_n");
    String ms1_span;
    // [parallel] MS1 tracing was UNBANDED: the call never passed a band count, so it defaulted to
    // 1 and MassTraceDetection -- a serial greedy loop over an intensity-sorted apex list -- ran on
    // ONE core while MS2 was banded across the machine. Invisible at default thresholds (MS1 is ~4%
    // of frames), fatal for the MS1 sensitivity sweep, which is the one lever the literature names
    // as the root cause of low-abundance loss in this tool class.
    int ms1_bands = getIntOption_("perf:ms1_trace_bands");
    if (ms1_bands < 1) ms1_bands = 1;
    double _t_ms1 = phase_clock_(), _c_ms1 = cpu_seconds_();   // [perf-instr]
    // detectTraces_ distributes its bands as TASKS, which requires an enclosing parallel region;
    // the window loop provides one, this phase must open its own. It also lets MS1 tracing use the
    // whole machine instead of the band count (it was measured at only 5.0x on 12 bands).
    vector<Trace> ms1_traces;
    // The MS1 traces' frame table is the MS1 map's own spectrum sequence (RT-sorted above); their
    // spans live in this store, which must outlive scoring.
    TraceStore ms1_store;
    { vector<uint32_t> ri; ri.reserve(ms1_map.size());
      for (const auto& sp : ms1_map) ri.push_back(rtIndex(sp.getRT()));
      ms1_store.setFrames(ri, nullptr); }
    #pragma omp parallel
    #pragma omp single
    ms1_traces = detectTraces_(ms1_map, ms1_store, delta_im, getDoubleOption_("trace:noise_threshold_int"),
                                             ms1_snr, getDoubleOption_("trace:min_length_sec"), getDoubleOption_("trace:ms1_min_sample_rate"), max_trace_len, getDoubleOption_("trace:ms1_split_valleys"), &ms1_span, ms1_bands);
    for (auto& t : ms1_traces) t.st = &ms1_store;          // the store is immutable from here
    { const double cap = getDoubleOption_("trace:max_span_sec");
      if (cap > 0.0) for (auto& t : ms1_traces) trimToSpan(t, cap); }
    if (detOn_()) writeLogInfo_("[det] MS1 traces n=" + String(ms1_traces.size()) + " digest=" + String(traceDigest_(ms1_traces)));
    { auto& st = phase_stats_()["MS1_TRACE"]; if (st.n++ == 0) phase_order_().push_back("MS1_TRACE");
      st.wall += phase_clock_() - _t_ms1; st.cpu += cpu_seconds_() - _c_ms1; st.rss_end_mb = rss_mb_(); }
    writeLogInfo_("MS1 " + ms1_span);   // [merged-trace] measure, do not assume

    // [fwhm] Measure the actual chromatographic peak width from the MS1 traces just built, and
    // use it to set the MS2 valley-splitting window. chrom_fwhm is an ABSOLUTE TIME, so a fixed
    // default cannot be right across methods: 14.0 s gave +0.0% peptides on a 31 min gradient
    // and -1.06% on a 5.6 min one, where the same 14 s is proportionally 5.5x wider.
    // Median, not mean: FWHM distributions are heavy-tailed and a mean is dominated by a handful
    // of merged or truncated traces. Traces that never fall to half maximum on both sides are
    // EXCLUDED rather than imputed -- a truncated peak is not a narrow one, and imputing it
    // would bias the estimate downward exactly where peaks are broadest.
    double ms1_fwhm_med = -1.0;
    {
      vector<double> fw;
      fw.reserve(ms1_traces.size() / 4 + 1);
      vector<double> prt, pv;
      for (const auto& t : ms1_traces)
      {
        if (t.np() < 5) continue;
        packReal(t, prt, pv);                         // walk REAL points: a gap is not a zero crossing
        size_t ap = 0;
        for (size_t i = 1; i < pv.size(); ++i) if (pv[i] > pv[ap]) ap = i;
        const double half = pv[ap] * 0.5;
        if (half <= 0.0) continue;
        double lo = -1.0, hi = -1.0;
        for (size_t i = ap; i > 0; --i)
          if (pv[i - 1] <= half)
          { const double d = pv[i] - pv[i - 1];
            const double f = d > 0 ? (pv[i] - half) / d : 0.0;
            lo = prt[i] - f * (prt[i] - prt[i - 1]); break; }
        for (size_t i = ap; i + 1 < pv.size(); ++i)
          if (pv[i + 1] <= half)
          { const double d = pv[i] - pv[i + 1];
            const double f = d > 0 ? (pv[i] - half) / d : 0.0;
            hi = prt[i] + f * (prt[i + 1] - prt[i]); break; }
        if (lo >= 0.0 && hi > lo) fw.push_back(hi - lo);
      }
      if (!fw.empty())
      {
        std::nth_element(fw.begin(), fw.begin() + fw.size() / 2, fw.end());
        ms1_fwhm_med = fw[fw.size() / 2];
        writeLogInfo_("Measured MS1 FWHM: median " + String(ms1_fwhm_med) + " s over "
                      + String(fw.size()) + " traces (of " + String(ms1_traces.size())
                      + "; truncated peaks excluded)");
      }
      else writeLogInfo_("Measured MS1 FWHM: NO usable traces -- falling back to absolute seconds");
    }
    // ms2_split is declared at :1324 but the FWHM is only known HERE, after MS1 tracing --
    // that ordering is why only the MS2 side can use the measured width for free. Doing the
    // same for MS1 needs a two-pass trace (measure with splitting off, then re-trace), which is
    // a real cost and is deliberately not done blind.
    // [wavelet] The smoothing scale is a MEASUREMENT, not a parameter: trace:wavelet_smooth is a
    // multiple of the median MS1 FWHM just measured above. A fixed absolute width cannot be right
    // across gradients (the same 14 s was +0.0% on a 31 min method and -1.06% on a 5.6 min one).
    // Refuse rather than guess if no FWHM could be measured -- silently falling back to some
    // absolute default is how the split_valleys work produced method-dependent results.
    const double wl_mult = getDoubleOption_("trace:wavelet_smooth");
    if (wl_mult > 0.0 && ms1_fwhm_med > 0.0)
    {
      wl_scale_ = wl_mult * ms1_fwhm_med;
      writeLogInfo_("[wavelet] precursor-XIC smoothing ON: " + String(wl_mult) + " x measured FWHM "
                    + String(ms1_fwhm_med) + " s = " + String(wl_scale_) + " s");
    }
    else if (wl_mult > 0.0)
    {
      wl_scale_ = 0.0;
      writeLogWarn_("[wavelet] trace:wavelet_smooth requested but NO MS1 FWHM could be measured; "
                    "smoothing stays OFF rather than guessing an absolute width");
    }
    const double sv_fwhm = getDoubleOption_("trace:split_valleys_fwhm");
    if (sv_fwhm > 0.0 && ms1_fwhm_med > 0.0)
    {
      const double old_split = ms2_split;
      ms2_split = sv_fwhm * ms1_fwhm_med;
      writeLogInfo_("MS2 valley-split window set from measured FWHM: " + String(sv_fwhm) + " x "
                    + String(ms1_fwhm_med) + " s = " + String(ms2_split) + " s (was "
                    + String(old_split) + " s absolute)");
    }
    else if (sv_fwhm > 0.0)
      writeLogWarn_("trace:split_valleys_fwhm requested but no FWHM could be measured; "
                    "keeping the absolute trace:ms2_split_valleys value");
    const double iso_im_tol = getDoubleOption_("charge:iso_im_tolerance");
    const double ms2_noise = getDoubleOption_("trace:ms2_noise_threshold_int");
    writeLogInfo_("MS2 trace gate: apex > " + String(ms2_snr * ms2_noise) + " (snr " + String(ms2_snr)
                  + " x noise " + String(ms2_noise) + "), membership > " + String(ms2_noise)
                  + ", min_length " + String(ms2_minlen) + " s");
    const bool env_scoring = (getStringOption_("charge:scoring") == "envelope");
    const double ambig_margin = getDoubleOption_("charge:ambiguity_margin");
    // [determinism] Canonically order the MS1 traces before precursor inference. inferPrecursors_
    // is a GREEDY seed-and-mark-used loop; its output depends on ms1_traces' CONTAINER ORDER, and
    // that order is non-deterministic -- MS1 detectTraces_ runs at top level so EPD's internal
    // parallel-for is active and appends split traces in thread-completion order. The greedy
    // seed sort already tie-breaks on index (a<b), but the index maps to a different physical
    // trace each run. Result: ~0.7% run-to-run peptide wobble, and the deterministic -0.8% when a
    // parallel change reshuffled the order. A total physical sort here mirrors the MS2-side sort
    // that already exists (fragment traces are canonicalised right after detectTraces_) and makes
    // the whole tool reproducible -- which also makes the emitted-count/peptide correctness gate
    // MEASURABLE. [adv-review kimi/codex 2026-07-26]
    // A non-finite coordinate makes the comparator below non-transitive, which is UB in sort()
    // and silently corrupts every lower_bound on the result. The fragment path already filtered;
    // this one did not. [adv-review codex/kimi 2026-09-03]
    {
      const size_t before = ms1_traces.size();
      ms1_traces.erase(remove_if(ms1_traces.begin(), ms1_traces.end(), [](const Trace& t) {
        return !std::isfinite(t.mz) || !std::isfinite(t.rt) || !std::isfinite(t.im)
               || !std::isfinite(t.intensity);
      }), ms1_traces.end());
      if (before != ms1_traces.size())
        writeLogInfo_("dropped " + String(before - ms1_traces.size()) + " MS1 traces with non-finite coordinates");
    }
    sort(ms1_traces.begin(), ms1_traces.end(), [](const Trace& a, const Trace& b) {
      if (a.mz != b.mz) return a.mz < b.mz;
      if (a.rt != b.rt) return a.rt < b.rt;
      if (a.im != b.im) return a.im < b.im;
      return a.intensity < b.intensity;
    });
    double _t_win = 0, _c_win = 0;                             // [perf-instr] window loop (function scope)
    double _t_inf = phase_clock_(), _c_inf = cpu_seconds_();   // [perf-instr]
    vector<Precursor_> precursors = inferPrecursors_(ms1_traces, max_charge, delta_rt, iso_im_tol,
                                                     mass_ppm, env_scoring, ambig_margin);
    if (detOn_())
    {
      vector<uint64_t> b;
      for (const auto& pc : precursors) { uint64_t u; double v; v = pc.rt; std::memcpy(&u, &v, 8); b.push_back(u); v = pc.im; std::memcpy(&u, &v, 8); b.push_back(u); }
      writeLogInfo_("[det] precursors n=" + String(precursors.size()) + " digest=" + String(detDigest_(std::move(b))));
    }
    if (getFlag_("assembly:dedup_precursors"))
    {
      const size_t before = precursors.size();
      const size_t removed = dedupPrecursors_(precursors, ms1_traces, mass_ppm, delta_rt, delta_im);
      writeLogInfo_("Precursor dedup: removed " + String(removed) + " of " + String(before)
                    + " hypotheses (" + String((int)(1000.0 * removed / std::max<size_t>(before, 1)) / 10.0) + "%)");
    }
    // Emulate the reference implementation, which assigns a charge to every precursor: give isotope-unsupported
    // precursors the default charge (0 = keep the strict isotope-only behavior and drop them).
    const int default_charge = getIntOption_("assembly:default_charge");
    // rp_max = per-fragment precursor cap. 0 = share-all (stream, memory-light). competitive flag = rp_max 1.
    int rp_max = getIntOption_("assembly:rp_max");
    const double apportion = getDoubleOption_("assembly:apportion");
    if (getFlag_("assembly:competitive") && rp_max <= 0) rp_max = 1;
    Size n_default = 0;
    const bool open_safe = getFlag_("assembly:open_search_safe");
    // [way-4] A GUESSED charge and an isotope-SUPPORTED charge must not be emitted identically.
    // Closed search hides the difference (engine enumerates); open search does not.
    for (auto& pc : precursors)
      if (pc.charge == 0) { pc.guessed = true; pc.charge = open_safe ? 0 : default_charge; ++n_default; }
    // [emission] Isotope-support gate: drop the guessed singletons rather than emit them. They are
    // coordinate-distinct (merge/dedup cannot cut them) and identify 8.8x less often; dropping all
    // guessed cost 0.78% peptides for -39.5% emission on dataset D. OFF by default (open-search caveat).
    if (getStringOption_("assembly:require_isotope_support") == "true")
    {
      const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [](const Precursor_& pc){ return pc.guessed; }), precursors.end());
      writeLogInfo_("require_isotope_support: dropped " + String(before - precursors.size())
                    + " guessed precursors, " + String(precursors.size()) + " isotope-supported remain");
    }
    // [min-charge, 2026-09-03] 29.3% of dataset D emission is SINGLY-charged precursor hypotheses (271,559
    // spectra) and they carry 259 unique peptides, 1.73%. They are identified on 0.42% of their own
    // spectra against 15.3% for z=2 -- 36x less often. Tryptic peptides are essentially never 1+ in
    // ESI, and charge inference breaks ties toward the LOW charge, so this is where the mis-assignments
    // land. Dropping them is a far better emission lever than any isotope collapse tested: -29.3%
    // emission for -1.73% peptides, against -15.3% for -4.6% at the best collapse setting.
    {
      const int zmin = getIntOption_("charge:min_charge"); const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [zmin](const Precursor_& pc){ return pc.charge > 0 && pc.charge < zmin; }), precursors.end());
      if (before != precursors.size())
        writeLogInfo_("charge:min_charge=" + String(zmin) + ": dropped " + String(before - precursors.size())
                      + " precursors below that charge, " + String(precursors.size()) + " remain");
    }
    // [step-02 emission-controlled arm, 2026-09-02] SPEXTRACTOR_MIN_ISOTOPES=k keeps only precursors whose
    // envelope has >= k isotope peaks (n_isotopes counts the mono): a precursor QUALITY gate that cuts
    // emission without touching fragment sharing. Pre-registered falsifier in dataset D-BASELINE.
    if (const char* mi = std::getenv("SPEXTRACTOR_MIN_ISOTOPES"))
    {
      const int k = std::atoi(mi); const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [k](const Precursor_& pc){ return pc.n_isotopes < k; }), precursors.end());
      writeLogInfo_("SPEXTRACTOR_MIN_ISOTOPES=" + String(k) + ": dropped " + String(before - precursors.size())
                    + " precursors with fewer isotope peaks, " + String(precursors.size()) + " remain");
    }
    // [step-02 sub-arm 2] SPEXTRACTOR_PRECURSOR_LIST=<tsv: rt_sec mz z, header line> keeps only precursors that
    // MATCH a listed one (|dRT| <= 10 s, |dm/z| <= 10 ppm, same z or listed z == 0). With the reference implementation's own
    // precursor list this is the precursor-matched sub-arm: identical precursor population, our spectra.
    if (const char* pl = std::getenv("SPEXTRACTOR_PRECURSOR_LIST"))
    {
      struct Ref { double mz, rt; int z; };
      std::vector<Ref> refs; std::ifstream f(pl); std::string line; std::getline(f, line);   // header
      while (std::getline(f, line)) { std::istringstream ls(line); Ref r{}; if (ls >> r.rt >> r.mz >> r.z) refs.push_back(r); }
      std::sort(refs.begin(), refs.end(), [](const Ref& a, const Ref& b){ return a.mz < b.mz; });
      const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(), [&](const Precursor_& pc) {
        const double tol = pc.mono_mz * 10e-6;
        auto lo = std::lower_bound(refs.begin(), refs.end(), pc.mono_mz - tol, [](const Ref& r, double v){ return r.mz < v; });
        for (auto it = lo; it != refs.end() && it->mz <= pc.mono_mz + tol; ++it)
          if (std::fabs(it->rt - pc.rt) <= 10.0 && (it->z == 0 || pc.charge == 0 || it->z == pc.charge)) return false;
        return true; }), precursors.end());
      writeLogInfo_("SPEXTRACTOR_PRECURSOR_LIST: " + String(refs.size()) + " reference precursors; kept "
                    + String(precursors.size()) + " of " + String(before) + " (matched within 10 s / 10 ppm / z)");
    }
    // [ms1-funnel] Measured recall against a DIA-NN truth set is 74.3% vs the reference implementation's 84.4%, but a
    // single number cannot say WHERE a precursor is lost: no MS1 trace at all, a trace that never
    // became a hypothesis, or a hypothesis that never produced a spectrum. Dump both upstream
    // stages so the loss can be attributed to one of them instead of guessed at.
    const String ms1_dump = getStringOption_("diag:dump_ms1_tsv");
    if (!ms1_dump.empty())
    {
      // [fwhm] FWHM straight from each trace's own XIC rather than from
      // ElutionPeakDetection, so the measurement does not depend on whether split_valleys ran.
      // Apex, then walk outwards to half maximum with linear interpolation between samples;
      // -1 when the trace never falls to half max on both sides (truncated peak), which must be
      // reported rather than silently imputed -- a truncated peak is not a narrow one.
      auto fwhm_of = [](const Trace& tr) -> double {
        if (tr.np() < 3) return -1.0;
        vector<double> prt, pv; packReal(tr, prt, pv);
        size_t ap = 0;
        for (size_t i = 1; i < pv.size(); ++i) if (pv[i] > pv[ap]) ap = i;
        const double half = pv[ap] * 0.5;
        double lo = prt[0], hi = prt[pv.size() - 1];
        for (size_t i = ap; i > 0; --i)
          if (pv[i - 1] <= half)
          {
            const double d = pv[i] - pv[i - 1];
            const double f = d > 0 ? (pv[i] - half) / d : 0.0;
            lo = prt[i] - f * (prt[i] - prt[i - 1]);
            break;
          }
        for (size_t i = ap; i + 1 < pv.size(); ++i)
          if (pv[i + 1] <= half)
          {
            const double d = pv[i] - pv[i + 1];
            const double f = d > 0 ? (pv[i] - half) / d : 0.0;
            hi = prt[i] + f * (prt[i + 1] - prt[i]);
            break;
          }
        return hi - lo;
      };
      ofstream ft((ms1_dump + ".traces.tsv").c_str());
      // [model] 12 significant digits, not the ostream default 6: at m/z 765 six digits quantise
      // to ~1 mDa = 1.3 ppm inside a 15 ppm isotope window, which both flips borderline partner
      // matches and destroys the ppm-residual feature a charge model needs. Every dump written
      // before 2026-07-23 carries that quantisation.
      ft.precision(12);
      ft << "mz\trt\tim\tintensity\tfwhm_s\tn_xic\n";
      for (const auto& t : ms1_traces)
        ft << t.mz << '\t' << t.rt << '\t' << t.im << '\t' << t.intensity << '\t'
           << fwhm_of(t) << '\t' << t.np() << '\n';
      ft.close();
      ofstream fp((ms1_dump + ".precursors.tsv").c_str());
      fp.precision(12);
      // [model] trace_idx makes the hypothesis -> seed-trace link explicit, so an offline replica
      // of inferPrecursors_ can be checked row-for-row rather than assumed to agree.
      fp << "mono_mz\tcharge\trt\tim\tguessed\ttrace_idx\n";
      for (const auto& pc : precursors)
        fp << pc.mono_mz << '\t' << pc.charge << '\t' << pc.rt << '\t' << pc.im << '\t'
           << (pc.guessed ? 1 : 0) << '\t' << pc.trace_idx << '\n';
      fp.close();
      writeLogInfo_("MS1 funnel dump: " + String(ms1_traces.size()) + " traces, "
                    + String(precursors.size()) + " precursors -> " + ms1_dump + ".{traces,precursors}.tsv");
    }
    // sort precursors by mono m/z so each window can binary-search its members instead of
    // rescanning all ~1M precursors per window (total order for determinism). [code-review]
    sort(precursors.begin(), precursors.end(), [](const Precursor_& a, const Precursor_& b) {
      if (a.mono_mz != b.mono_mz) return a.mono_mz < b.mono_mz;
      if (a.rt != b.rt) return a.rt < b.rt;
      if (a.im != b.im) return a.im < b.im;
      return a.trace_idx < b.trace_idx;
    });
    vector<double> prec_mz(precursors.size());
    for (size_t i = 0; i < precursors.size(); ++i) prec_mz[i] = precursors[i].mono_mz;
    // NOTE: ms1_map is already empty here. detectTraces_ takes it by non-const reference and
    // releases every spectrum as it distributes it into the m/z bands, then clears the container
    // itself; a clear() at this point frees nothing. The bytes it used to be credited with are
    // released one phase earlier -- into the picking threads' allocator arenas, which is why the
    // resident set does not fall there either. [mem]
    // MEMORY: scoring reads the profile of ms1_traces[pc.trace_idx], so ms1_traces must STAY -- but only the
    // traces an actual precursor points at are ever dereferenced. Free the xic of the rest
    // (2.5M traces vs 1.2M precursors => roughly half are dead weight held through scoring). [mem]
    {
      vector<bool> needed(ms1_traces.size(), false);
      for (const auto& pc : precursors)
        if (pc.trace_idx < ms1_traces.size()) needed[pc.trace_idx] = true;
      // The spans share one arena, so releasing a trace means COMPACTING, in place (a second arena
      // would double the MS1 arena at exactly the moment it is largest) and in OFFSET order -- see
      // compactUnreferenced() for why container order corrupted ~1% of the kept spans.
      const Size freed = compactUnreferenced(ms1_traces, ms1_store, needed);
      writeLogInfo_("Released XICs of " + String(freed) + " unreferenced MS1 traces." + rss_() + mem_());
    }
    writeLogInfo_(clk_() + " Detected " + String(ms1_traces.size()) + " MS1 traces -> " + String(precursors.size())
                  + " precursor hypotheses (" + String(n_default) + " assigned default charge)." + rss_());
    std::exception_ptr worker_err;

    //-------------------------------------------------------------
    // PARALLEL over windows (independent; with ~24 windows this fully uses the cores, and each
    // window's trace detection AND gating run on its own core). Each window frees its raw frames
    // right after trace extraction. Per-window buckets are merged in fixed window order then
    // canonically sorted, so output is thread-count-independent. Exceptions (e.g. MassTraceDetection
    // on a degenerate window) are captured and rethrown serially. [C-3,H-5,C-5,par-Crit-2/3/4]
    //-------------------------------------------------------------
    vector<pair<WinKey, vector<CompactFrame>*>> window_list;
    window_list.reserve(ms2_by_window.size());
    for (auto& kv : ms2_by_window) window_list.emplace_back(kv.first, &kv.second);
    // [heavy-first] Process the biggest windows FIRST. With n_conc slots and unequal windows,
    // a small window admitted late runs nearly alone in the tail; starting the largest first
    // gives the scheduler the whole run to overlap it. The widest tile holds ~23.7% of all
    // precursors, so tail order dominates wall time. Cost proxy = compact bytes (peaks to trace),
    // which correlates with both trace and score work. Output order is unaffected: a canonical
    // total-order sort over all emitted spectra runs after the loop, so reordering here is safe.
    {
      auto wbytes = [](const vector<CompactFrame>& fr) { size_t b = 0; for (const auto& f : fr) b += f.bytes(); return b; };
      { auto& st = phase_stats_()["PRECURSOR_INFER"]; if (st.n++ == 0) phase_order_().push_back("PRECURSOR_INFER");
        st.wall += phase_clock_() - _t_inf; st.cpu += cpu_seconds_() - _c_inf; st.rss_end_mb = rss_mb_(); }
      _t_win = phase_clock_(); _c_win = cpu_seconds_();          // [perf-instr] window loop
      std::stable_sort(window_list.begin(), window_list.end(),
                       [&](const auto& a, const auto& b) { return wbytes(*a.second) > wbytes(*b.second); });
    }
    vector<vector<MSSpectrum>> win_out(window_list.size());
    // Fan-out histogram (only filled when rp_max>0): buckets = #precursors each fragment gates,
    // BEFORE capping. Answers "does this rp_max actually bite?" Buckets: 1,2,3,4,5,6-10,11-25,26-50,51+.
    vector<uint64_t> fanout_hist(9, 0);
    Size win_done = 0;
    // MEMORY vs SPEED: every concurrently-processed window holds its own frames + traces + grid, so
    // peak RSS scales with the number of windows in flight, not with the window count. Capping
    // concurrency trades wall time for RAM. 0 = unlimited (previous behaviour). [mem]
    const int max_conc = getIntOption_("perf:max_concurrent_windows");
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    // [bands] The window count is FIXED BY THE ACQUISITION METHOD (24 tiles here), so
    // min(windows, threads) pinned us to 24 cores no matter the machine -- measured: 24
    // concurrent windows on a 224-core node ran at load ~20. The ceiling is now the thread
    // count, with each window split into m/z bands to supply the extra work units.
    // [bands] PER-WINDOW band count, sized from the window itself rather than one global
    // constant. Measured motivation: tracing is 69.3% of window-loop CPU (materialize 0.10,
    // trace 69.34, fraggrid 1.75, score 28.81, emit 0.01), and MassTraceDetection has no OpenMP
    // of its own -- banding is the only way to parallelise it.
    //
    // A global count is wrong because the HALO IS A FIXED WIDTH PER BAND (20x mass_error_ppm on
    // each side, which is what makes the partition exact). Its overhead therefore depends on how
    // wide the band is:
    //     tile 328-469 Da (141 wide):  24 bands ->  5% halo overhead
    //     tile 564-583 Da (19 wide):   24 bands -> 57% halo overhead
    // So a wide tile tolerates many bands and a narrow one does not. Sizing per window also
    // hands the most parallelism to the widest tile -- which holds 23.7% of all precursors and
    // is the Amdahl tail of this loop.
    //
    // min_band_da = 8x the halo keeps overhead near 12% at the narrow end. Bands are capped by
    // the thread budget so windows x bands cannot oversubscribe.
    int n_bands = getIntOption_("perf:trace_bands");
    if (n_bands <= 0)  // auto: enough bands that windows x bands covers every thread
      n_bands = (int)std::ceil((double)n_threads / std::max<size_t>(window_list.size(), 1));
    n_bands = std::max(1, n_bands);
    // TOTAL live threads are outer x inner. Capping only the inner team (as a first cut did)
    // gave 24 windows x 9 bands = 216 threads on a 224-core box that another job was already
    // loading to ~173 -- measured: 11 of 24 windows in 98 min against a 40 min single-band
    // reference. Divide the thread budget between the two levels instead of multiplying it.
    // [master/worker] n_conc no longer partitions the thread pool -- tasks do that. It is now
    // only the number of windows allowed IN FLIGHT, i.e. how many materialised windows may hold
    // memory at once, and the free-RAM admission gate below is the real bound. Deriving it from
    // threads/bands (the old 100/12 = 8) was what pinned the tool to a 96-thread grid.
    // MEASURED 2026-09-03: letting every window be in flight raised window-loop occupancy from
    // 64.8x to 86.2x -- and made the run SLOWER (15:50 against 12:53) because peak RSS went 103 ->
    // 152 GB and system time rose 75%: more resident windows means a bigger working set, more page
    // faults and worse locality, not more useful work. Occupancy and working set are separate
    // knobs and the task pool only fixes the first. So keep a default cap on windows IN FLIGHT --
    // the old threads/bands value, which was accidentally a reasonable working-set bound -- while
    // the pool goes on filling every thread from whatever is admitted.
    // I capped this earlier because the task pool alone had been SLOWER (15:50 vs 12:53) and I
    // read that as an oversized working set. The measurement says otherwise: once the scoring gate
    // stopped striding 96-byte records, uncapped is slightly faster (7:04 vs 7:15) at the SAME peak
    // memory (157 GB both ways). The pool's regression was the gate's bandwidth problem amplified
    // by concurrency, not a working-set problem of its own -- so the cap is not needed, and taking
    // it out costs nothing. `perf:max_concurrent_windows` still bounds it for a small machine.
    int n_conc = (int)window_list.size();
    if (max_conc > 0) n_conc = std::min(n_conc, max_conc);
    if (n_conc < 1) n_conc = 1;
#ifdef _OPENMP
    // outer loop over windows + inner loop over bands = two live levels
    if (n_bands > 1) omp_set_max_active_levels(2);
#endif
    writeLogInfo_("Parallelism: " + String(window_list.size()) + " windows x " + String(n_bands)
                  + " m/z bands = " + String((int)window_list.size() * n_bands) + " trace units; "
                  + String(n_threads) + " worker threads over one task pool; up to "
                  + String(n_conc) + " windows in flight (memory-bounded)");
    // [dyn-mem] Thread count is the UPPER bound; free RAM is the real one. A static cap either
    // wastes cores on a big machine or OOMs on a shared one, and this node is shared. Admit each
    // window only when its projected footprint fits in currently-free RAM, re-read from
    // /proc/meminfo at each admission so the run yields to other jobs instead of fighting them.
    // Always admit when nothing is in flight, so the loop cannot deadlock on a single huge window.
    bool integer_detector = (getStringOption_("trace:detector") == "integer");
    String why = "not attempted";
    if (integer_detector && !tofAxis().ok)
    {
      // The flight-time axis comes from the input's own analysis.tdf: a .d directory has one, and an
      // mzPeak archive can be accompanied by one (SPEXTRACTOR_MZPEAK_TDF, as the exact-calibration
      // path already uses).
      why = "no candidate tdf";
      std::vector<std::string> cand;
      cand.push_back(in + "/analysis.tdf");
      if (const char* sc = std::getenv("SPEXTRACTOR_MZPEAK_TDF")) cand.push_back(sc);
      for (const std::string& c : cand)
        if (loadTofAxis(c, why)) { writeLogInfo_("flight-time axis from " + c + ": "
              + String(tofAxis().b_by_frame.size()) + " frame factors"); break; }
    }
    // FAIL CLOSED on frames the calibration cannot address. A .d names every frame "frame=<Id>";
    // an mzML (or an archive index) may not, and frameIdOf() returns 0 for anything it cannot
    // parse. With index 0 backfilled to the reference factor those frames used to calibrate
    // SILENTLY at reference temperature, skipping the warning below entirely -- an mzML plus a
    // sidecar tdf loaded "fine" and produced masses that were never calibrated per frame. The
    // temperature term is small (0.669 ppm, tests/test_calibration.py) but it is exactly the
    // invented calibration this path exists to refuse, so an unmapped frame now costs the integer
    // detector rather than being papered over.
    if (integer_detector && tofAxis().ok && cstat.unmapped_frame > 0)
    {
      writeLogWarn_("trace:detector=integer needs a per-frame calibration factor for EVERY frame, and "
                    + String(cstat.unmapped_frame) + " frame(s) carry a nativeID with no parseable "
                    "\"frame=<Id>\" key, so they cannot be addressed in the tdf's Frames table. Falling "
                    "back to trace:detector=openms rather than calibrating them at reference "
                    "temperature. A Bruker .d always names its frames; an mzML converted from one may "
                    "not. Pass -trace:detector openms to silence this.");
      integer_detector = false;
    }
    if (integer_detector && !tofAxis().ok)
    {
      // FALL BACK, do not refuse. Failing closed is right when the alternative is inventing a
      // calibration -- that would put wrong masses in the output. Here there is a second detector
      // that needs no calibration at all, so refusing would only make the tool unusable on inputs
      // it handled before this became the default. Which one ran is recorded in the output as
      // spx:detector, so a file is always attributable.
      writeLogWarn_("trace:detector=integer needs the vendor flight-time calibration (MzCalibration "
                    "+ Frames.T1 from analysis.tdf) and it was not available for this input ("
                    + why + "). Falling back to trace:detector=openms; the two detectors agree on about "
                    "85% of the union of identified peptides, so this file is NOT comparable peptide-for-peptide "
                    "with one traced on the integer axis. Pass -trace:detector openms to silence this.");
      integer_detector = false;
    }
    const double mem_frac = getDoubleOption_("perf:mem_fraction");
    std::atomic<size_t> inflight{0};
    // Projected footprint per window, from its own compact size. RE-CALIBRATED 2026-09-03: the
    // 10x multiplier booked the average dataset A window at 8.9 GB against a measured 12-14 GB, i.e. it
    // under-booked by 35-55% in the UNSAFE direction -- it admitted more windows than fit. That was
    // survivable while the thread grid capped concurrency at 8 anyway; with the task pool, memory
    // is the only bound, so the estimate has to be honest. [adv-review kimi 2026-09-03]
    auto project = [&](const vector<CompactFrame>& fr) {
      size_t b = 0; for (const auto& f : fr) b += f.bytes();
      return (size_t)(b * 16.0) + (size_t)256 * 1024 * 1024;
    };
    writeLogInfo_("Processing " + String(window_list.size()) + " isolation windows (OpenMP over windows, <="
                  + String(n_conc) + " concurrent, admission bounded by " + String((int)(mem_frac * 100))
                  + "% of free RAM)..." + rss_() + mem_());

    // [timers] phase_clock_ is never called between the start of this loop and its end, so
    // EVERYTHING said about where the loop's 85.9% goes has been inference. Two performance
    // fixes were already attempted and reverted against a model that was never verified. These
    // are atomic nanosecond accumulators -- one add per stage per window, so the measurement
    // cost is unmeasurable against stages that run for minutes.
    std::atomic<long long> t_mat{0}, t_trace{0}, t_grid{0}, t_score{0}, t_emit{0};
    std::atomic<bool> span_logged{false};
    // [master/worker, 2026-09-03] The loop used to be `parallel for num_threads(n_conc)` with
    // nested `parallel for num_threads(n_bands)` inside, i.e. a RIGID n_conc x n_bands grid
    // (8 x 12 = 96 of 100 threads). Two costs followed from the rigidity: the grid never fills the
    // last 4 threads, and when fewer windows remain than slots the survivors keep only their own
    // 12-thread team while the rest of the machine idles -- measured 67.9x of 100 on dataset A, 71% of
    // the grid's own ceiling. Here ONE thread (the master) walks the windows heaviest-first and
    // throttles on memory, and every unit of work inside a window is an OpenMP TASK. Idle threads
    // steal across window boundaries, so the tail of one window is filled by another's bands.
    // Output is unaffected: every result is written to an index-addressed slot (win_out[wi],
    // pslot[pi-plo], per[b]), never appended in completion order.
    #pragma omp parallel num_threads(n_threads)
    #pragma omp single
    for (long wi = 0; wi < (long)window_list.size(); ++wi)
    {
      #pragma omp task default(shared) firstprivate(wi)
      {
      // A lambda, because `continue` is not a legal exit from a task's structured block and the
      // window body has one early exit (a window whose fragments all fail the finiteness filter).
      [&]() {
      try
      {
        const double win_lo = window_list[wi].first[0] / 100.0;
        const double win_hi = window_list[wi].first[1] / 100.0;
        // [dyn-mem] wait for room before materialising anything
        const size_t need = project(*window_list[wi].second);
        for (;;)
        {
          bool go = false;
          #pragma omp critical(admit)
          {
            const size_t budget = (size_t)(availableBytes_() * mem_frac);
            // "at least one": an empty pipeline always admits, however large the window
            if (inflight.load() == 0 || inflight.load() + need <= budget) { inflight += need; go = true; }
          }
          if (go) break;
          // The master is one of the workers: run queued tasks while waiting for memory rather
          // than sleeping on a thread the pool could be using.
          #pragma omp taskyield
          std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
        struct Release { std::atomic<size_t>& c; size_t n; ~Release() { c -= n; } } rel{inflight, need};
        // Per-window stage seconds and the completion line, printed on EVERY exit of the body
        // (empty window, any assembly arm, exception): what the profile could only infer. [mem]
        double w_prep = 0, w_trace = 0, w_split = 0, w_score = 0, w_emit = 0;
        auto secs = [](auto d) { return std::chrono::duration<double>(d).count(); };
        struct WinDone { std::function<void()> f; ~WinDone() { f(); } } win_done_guard{[&]() {
          #pragma omp critical(winlog)
          writeLogInfo_("  window " + String(++win_done) + "/" + String(window_list.size()) + " " + String(win_lo) + "-" + String(win_hi)
                        + ": prep " + String(w_prep) + " s | trace " + String(w_trace) + " s | split " + String(w_split)
                        + " s | score " + String(w_score) + " s | emit " + String(w_emit) + " s | spectra " + String(win_out[wi].size()) + clk_());
        }};
        // [compact] materialise ONLY this window, and free its compact frames as we go
        auto _t0 = std::chrono::steady_clock::now();
        // The integer detector reads the compact store directly; only the OpenMS path needs the
        // peaks converted back into Peak1D doubles, which is the largest per-window allocation.
        PeakMap wmap;
        PeakSlab wslab;
        if (integer_detector)
        {
          wslab = toSlab(*window_list[wi].second);
          if (detOn_())
          {
            #pragma omp critical
            writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " slab n=" + String(wslab.peaks()) + " digest=" + String(slabDigest_(wslab)));
          }
        }
        else
        {
          wmap = materializeWindow(*window_list[wi].second);
        }
        if (detOn_())   // [det] compact store digest BEFORE tracing: splits "pick" from "trace+assemble"
        {
          // materializeWindow consumed the frames; digest the materialised map instead (same content)
          vector<uint64_t> b;
          for (const auto& sp : wmap) for (const auto& pk : sp) { uint64_t u; double v = pk.getMZ(); std::memcpy(&u, &v, 8); b.push_back(u); float f = pk.getIntensity(); uint32_t w; std::memcpy(&w, &f, 4); b.push_back(w); }
          #pragma omp critical
          writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " compact n=" + String(b.size() / 2) + " digest=" + String(detDigest_(std::move(b))));
        }
        auto _t1 = std::chrono::steady_clock::now();
        t_mat += std::chrono::duration_cast<std::chrono::nanoseconds>(_t1 - _t0).count();
        for (MSSpectrum& s : wmap) s.setMSLevel(1); // MassTraceDetection traces MS1-level only
        wmap.sortSpectra();
        // [cross-frame] sum adjacent cycles for THIS window before detecting fragment traces
        if (agg_n > 1) aggregateFrames_(wmap, agg_n, mass_ppm, delta_im);
        String ms2_span;
        auto _t2 = std::chrono::steady_clock::now();
        // per-window band count: n_bands is the AUTO cap (threads/windows); a wide tile may exceed
        // it since bands are cheap relative to the trace work they unlock, a narrow one may not
        // reach it. Explicit -perf:trace_bands N pins every window to N.
        // MEASURED on dataset D: a flat 12 bands beats a per-window width-derived band count by 25% and beats
        // flat 24 as well (1,844 s vs 2,483 s vs 1,987 s). The heuristic under-provisioned
        // because it sized bands against HALO overhead, and the halo term turned out to be
        // dominated by a much larger effect: trace CPU falls 11,168 -> 5,378 -> 2,972 s across
        // 5/12/24 bands, i.e. banding REDUCES TOTAL WORK, not just wall time.
        // MassTraceDetection scans candidates per apex within the band's m/z range, so halving a
        // band halves the candidate list for every apex in it -- superlinear, and it swamps the
        // fixed per-band halo cost. 24 bands then loses to 12 on wall time despite half the CPU,
        // because 24 windows x 24 bands = 576 threads on 100 cores and scheduling overtakes the
        // algorithmic gain.
        // [auto-band-fix] Use the RESOLVED n_bands, not the raw option. In auto mode
        // (perf:trace_bands=0) n_bands was computed as ceil(threads/windows) at the top and
        // n_conc was sized from it, but the actual call here re-read the raw 0 and max(1,0)=1
        // traced every window SERIALLY while the log advertised multiple bands. Passing n_bands
        // keeps the advertised and actual band counts identical. Explicit values pass through
        // unchanged (n_bands == the option when the option is > 0).
        const int w_bands = n_bands;
        vector<Trace> frag_traces;
        TraceStore wst;                                   // this window's frame table + arenas; outlives scoring
        if (integer_detector)
        {
          wst.setFrames(wslab.rt_index, &wslab.b);
          // Banded like the OpenMS path. A band SEEDS only in its own flight-time core but may
          // extend anywhere, so no halo is needed and no trace is cut; the one approximation is
          // that two bands can take the same peak (each has its own `visited`), which is the same
          // class of error the OpenMS halo carries -- measured at 0.1% of peptides. Bands are
          // TASKS, so they fill the pool alongside every other window's work.
          if (agg_n > 1)
            throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "trace:frame_aggregation_n > 1 is not applied by trace:detector=integer");
          TofIdx tlo = 0xFFFFFFFFu, thi = 0;
          for (TofIdx t : wslab.tof) { tlo = std::min(tlo, t); thi = std::max(thi, t); }
          const int nb = std::max(1, w_bands);
          size_t n_seeds = 0, n_parents = 0;
          const size_t n_peaks = wslab.peaks();
          {
          auto _tp = std::chrono::steady_clock::now();
          const TracePrep tprep = prepareTracing(wslab, ms2_noise, ms2_snr);   // once per window
          w_prep = secs(std::chrono::steady_clock::now() - _tp);
          n_seeds = tprep.order.size();
          // Each band appends spans to a PRIVATE store; after the taskloop joins, the window
          // store absorbs them in band order and every offset is rebased (R4 sequencing).
          vector<TraceStore> bst(std::max(1, nb));
          if (thi <= tlo || nb == 1)
            frag_traces = detectTracesInteger_(wslab, tprep, bst[0], mass_ppm, delta_im, ms2_noise, ms2_minlen, ms2_msr, max_trace_len, max_span);
          else
          {
            vector<vector<Trace>> per(nb);
            const double wdt = (double)(thi - tlo + 1) / nb;
            #pragma omp taskloop grainsize(1) default(shared)
            for (int b = 0; b < nb; ++b)
            {
              const TofIdx blo = tlo + (TofIdx)(b * wdt);
              const TofIdx bhi = (b + 1 == nb) ? thi + 1 : tlo + (TofIdx)((b + 1) * wdt);
              per[b] = detectTracesInteger_(wslab, tprep, bst[b], mass_ppm, delta_im, ms2_noise, ms2_minlen, ms2_msr,
                                            max_trace_len, max_span, blo, bhi);
            }
            size_t tot = 0, arena = 0; bool any_bins = false;
            for (int b = 0; b < nb; ++b)
            { tot += per[b].size(); arena += bst[b].inten.size(); any_bins |= !bst[b].bins.empty(); }
            frag_traces.reserve(tot);
            wst.inten.reserve(arena); if (any_bins) wst.bins.reserve(arena);   // as in splitIntegerTraces
            writeLogInfo_("[mem] window " + String(win_lo) + "-" + String(win_hi) + " merge: " + String(nb) + " band arenas, "
                          + String(arena) + " floats (" + String(arena * (any_bins ? 8 : 4) / (1024ULL * 1024ULL)) + " MB) x2 while absorbed, " + String(tot) + " traces");
            for (int b = 0; b < nb; ++b)
            {
              const uint32_t base = wst.absorb(bst[b]);
              for (auto& t : per[b]) { t.off += base; frag_traces.push_back(std::move(t)); }
              vector<Trace>().swap(per[b]);
            }
          }
          if (thi <= tlo || nb == 1) { const uint32_t base = wst.absorb(bst[0]); for (auto& t : frag_traces) t.off += base; }
          for (auto& t : frag_traces) t.st = &wst;
          }                                                  // tprep dies here
          // The slab is dead too: valley splitting reads only the window store (which copied the
          // frame table) and the spans. Holding 12 GB of peaks plus the prep across the stage where
          // the MassTrace payloads peak was the largest lifetime defect the plan found.
          PeakSlab().frame_off.swap(wslab.frame_off);
          vector<TofIdx>().swap(wslab.tof); vector<float>().swap(wslab.inten);
          vector<double>().swap(wslab.b); vector<uint32_t>().swap(wslab.rt_index);
          vector<uint16_t>().swap(wslab.imq);
          // The same valley-splitting step the OpenMS path runs; skipping it was one of the
          // defects behind the first version's 92% peptide loss. It rebuilds the arena, so the
          // store pointers are set again afterwards; the per-frame bins are dead after it.
          n_parents = frag_traces.size();
          auto _ts = std::chrono::steady_clock::now();
          frag_traces = splitIntegerTraces(frag_traces, ms2_split, wst, std::max(1, w_bands) * 4);
          w_split = secs(std::chrono::steady_clock::now() - _ts);
          for (auto& t : frag_traces) t.st = &wst;
          if (max_span > 0.0) for (auto& t : frag_traces) trimToSpan(t, max_span);
          vector<uint32_t>().swap(wst.bins);
          writeLogInfo_("[mem] window " + String(win_lo) + "-" + String(win_hi) + ": seeds " + String(n_seeds) + " of " + String(n_peaks)
                        + " peaks (" + String((int)(1000.0 * n_seeds / std::max<size_t>(n_peaks, 1)) / 10.0) + "%); traces " + String(n_parents)
                        + " -> " + String(frag_traces.size()) + " (cap " + String(frag_traces.capacity()) + ", "
                        + String(frag_traces.capacity() * sizeof(Trace) / (1024ULL * 1024ULL)) + " MB); arena " + String(wst.inten.size()) + "/"
                        + String(wst.inten.capacity()) + " floats (" + String(wst.inten.capacity() * 4 / (1024ULL * 1024ULL)) + " MB)" + rss_());
        }
        else
        {
          { vector<uint32_t> ri; ri.reserve(wmap.size());          // wmap is RT-sorted above
            for (const auto& sp : wmap) ri.push_back(rtIndex(sp.getRT()));
            wst.setFrames(ri, nullptr); }
          frag_traces = detectTraces_(wmap, wst, delta_im, ms2_noise, ms2_snr, ms2_minlen, ms2_msr, max_trace_len, ms2_split, &ms2_span, w_bands);
          for (auto& t : frag_traces) t.st = &wst;
          if (max_span > 0.0) for (auto& t : frag_traces) trimToSpan(t, max_span);
        }
        // the window's spans and records, released with the window (the old add-only line was
        // cumulative over the run, not a peak)
        if (detOn_())
        {
          #pragma omp critical
          writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " frag traces n=" + String(frag_traces.size()) + " digest=" + String(traceDigest_(frag_traces)));
        }
        // (compact frames already released inside materializeWindow) [par-Crit-3]
        if (!span_logged.exchange(true)) { // [merged-trace] one window's span statistics, once
          #pragma omp critical
          writeLogInfo_("MS2 " + ms2_span);
        }

        frag_traces.erase(remove_if(frag_traces.begin(), frag_traces.end(), [](const Trace& t) {
          return !std::isfinite(t.im) || !std::isfinite(t.mz) || !std::isfinite(t.rt);
        }), frag_traces.end());
        if (frag_traces.empty()) return;   // nothing traceable in this window

        sort(frag_traces.begin(), frag_traces.end(), [](const Trace& a, const Trace& b) {
          if (a.im != b.im) return a.im < b.im;
          if (a.mz != b.mz) return a.mz < b.mz;
          if (a.rt != b.rt) return a.rt < b.rt;
          return a.intensity < b.intensity;
        });
        vector<double> frag_im(frag_traces.size());
        for (size_t i = 0; i < frag_traces.size(); ++i) frag_im[i] = frag_traces[i].im;
        FragRt frag_rt; frag_rt.build(frag_traces, delta_rt);

        // Precompute the common correlation grid once per window; pdense is reused across the
        // window's precursors (this window runs on one thread). [perf]
        auto _t3 = std::chrono::steady_clock::now();
        t_trace += std::chrono::duration_cast<std::chrono::nanoseconds>(_t3 - _t2).count();
        w_trace = secs(_t3 - _t2) - w_prep - w_split;
        FragStats fg = buildFragStats(frag_traces, wst, ms1_store, delta_rt);
        auto _t4 = std::chrono::steady_clock::now();
        t_grid += std::chrono::duration_cast<std::chrono::nanoseconds>(_t4 - _t3).count();
        vector<float> pdense(wst.frames(), 0.0f);   // indexed by window-local frame

        vector<MSSpectrum>& bucket = win_out[wi];
        // only this window's precursors (by mono m/z); NaN mono_mz sorts to the end, excluded.
        // (charge 0 is NOT dropped: emitted with charge UNSET so the engine searches its range.)
        size_t plo = lower_bound(prec_mz.begin(), prec_mz.end(), win_lo) - prec_mz.begin();
        size_t phi = upper_bound(prec_mz.begin(), prec_mz.end(), win_hi) - prec_mz.begin();
        if (apportion > 0.0 && rp_max <= 0)
        {
          // [route-4] INTENSITY APPORTIONMENT — the only assembly lever that is not a count knob.
          //
          // Route 3 measured that for ~80% (chance-corrected) of the peptides the reference implementation finds and we
          // miss, we DO have a precursor at the right (RT, m/z, IM): we emit a spectrum there, but the
          // right fragments never arrive with enough weight. Under share-all a fragment is copied at
          // FULL intensity into every gated precursor (mean fan-out 6.45), so a chimeric contribution
          // is indistinguishable from a private one and the true owner is diluted in rank.
          //
          // Every count-based remedy is falsified (rank-pruning, competitive, min_corr, min_corr_pts —
          // all monotonic losses). So: keep the fragment in ALL precursors (count unchanged) but split
          // its INTENSITY by relative correlation, w_i = corr_i^p / sum_j corr_j^p. The owner keeps
          // most of the signal, chimeric copies are down-weighted, and because max_fragments ranks on
          // score*intensity the ranking changes without any fragment being removed.
          // p = apportion (1 = linear, higher = sharper). p->inf degenerates to competitive, which is
          // falsified, so keep p modest.
          const size_t NF = frag_traces.size();
          vector<float> wsum(NF, 0.0f);
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            scoreCandidates_(pc, frag_traces, frag_im, frag_rt, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              wsum[fi] += (float)std::pow(std::max(c, 0.0), apportion);
            });
          }
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            vector<pair<double, double>> frags;
            vector<double> frag_scores;
            scoreCandidates_(pc, frag_traces, frag_im, frag_rt, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              const double w = wsum[fi] > 0.0f ? std::pow(std::max(c, 0.0), apportion) / wsum[fi] : 1.0;
              const auto [inten, emit_inten] = weighted_(frag_traces[fi].intensity * w, c, frag_im[fi] - pc.im);   // [E5]
              frags.emplace_back(exportMz_(frag_traces[fi].tof, frag_traces[fi].b, frag_traces[fi].mz), emit_inten);
              frag_scores.push_back(c * inten);
            });
            MSSpectrum ms2;
            assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, ms2);
            if (!ms2.empty()) bucket.push_back(std::move(ms2));
          }
        }
        else if (rp_max <= 0)
        {
          // Share-all (default, current best): each fragment goes to EVERY gated precursor.
          //
          // [perf-3] PRECURSOR-LEVEL PARALLELISM. The outer loop's unit of work is the isolation
          // window -- an ACQUISITION-METHOD CONSTANT (24 tiles) with one tile holding 23.7% of
          // all precursors, so Amdahl caps that loop at 4.2x regardless of thread count. /proc
          // sampling caught the tail running 1,993 s at 1.000 user cores with 21 of 22 threads
          // in futex_wait: one window left, scoring its precursors serially. Precursors within a
          // window are independent, so they are the correct unit.
          //
          // An earlier attempt at this was REVERTED for LOSING SPECTRA: it appended to the shared
          // `bucket`, and a race dropped emissions (caught by a 15 MB smaller mzML, not by the
          // peptide count). This writes into a PRE-SIZED slot per precursor -- no append, no
          // critical section, so a lost write is structurally impossible. `pdense` is the scratch
          // buffer assembleOne_ writes through, so it must be firstprivate: sharing it across
          // threads is precisely the corruption the reverted attempt suffered.
          //
          // Order-safety: a canonical total-order sort over all_out runs after the window loop
          // (RT, precursor m/z, charge, then fragment sequence), so emission order is already
          // thread-count-independent. Only slot CONTENTS matter here.
          // THREAD CAP -- this region is NESTED inside the window loop (num_threads(n_conc)).
          // Without a cap it requests the full team PER WINDOW, i.e. n_conc x omp_get_max_threads()
          // threads on a fixed core count. That is precisely the defect diagnosed in
          // ElutionPeakDetection.cpp:328 (a bare parallel-for with no num_threads, reached from
          // inside this same loop, measured at 43.5% kernel time) -- reproduced here by omission.
          // Divide the budget instead of multiplying it.
          auto _t5 = std::chrono::steady_clock::now();
          vector<MSSpectrum> pslot(phi > plo ? phi - plo : 0);
          // grainsize, not num_threads: the pool decides how many threads land here, so a window
          // that outlives its neighbours gets the whole machine instead of a fixed team of 12.
          #pragma omp taskloop grainsize(32) firstprivate(pdense) default(shared)
          for (long pii = (long)plo; pii < (long)phi; ++pii)
          {
            const size_t pi = (size_t)pii;
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            MSSpectrum ms2;
            assembleOne_(pc, win_lo, win_hi, frag_traces, frag_im, frag_rt, ms1_traces, fg,
                         delta_im, delta_rt, min_corr, min_corr_pts, min_frags, max_frags, pdense, ms2);
            if (!ms2.empty()) pslot[pi - plo] = std::move(ms2);
          }
          auto _t6 = std::chrono::steady_clock::now();
          t_score += std::chrono::duration_cast<std::chrono::nanoseconds>(_t6 - _t5).count();
          w_score = secs(_t6 - _t5);
          for (auto& s : pslot) if (!s.empty()) bucket.push_back(std::move(s));
          const auto _t7 = std::chrono::steady_clock::now();
          t_emit += std::chrono::duration_cast<std::chrono::nanoseconds>(_t7 - _t6).count();
          w_emit = secs(_t7 - _t6);
        }
        else
        {
          // Soft RP rank-pruning (rp_max>=1; 1 == competitive): keep each fragment in only its
          // top-rp_max best-correlating precursors. Pass A: per fragment, a bounded top-rp_max heap of
          // (corr,pi) + total gate count for the fan-out histogram. Pass B: assemble per winning precursor.
          const size_t NF = frag_traces.size();
          vector<vector<pair<double, long>>> topk(NF); // ascending by corr; back() = weakest kept
          vector<uint32_t> nfan(NF, 0);                 // total precursors that gated this fragment
          for (size_t pi = plo; pi < phi; ++pi)
          {
            const Precursor_& pc = precursors[pi];
            if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
            scoreCandidates_(pc, frag_traces, frag_im, frag_rt, ms1_traces, fg, delta_im, delta_rt, min_corr,
                             min_corr_pts, pdense, [&](size_t fi, double c) {
              ++nfan[fi];
              auto& tk = topk[fi];
              // insertion sort into a size-<=rp_max list kept ascending; deterministic tie-break on pi.
              if ((int)tk.size() < rp_max) {
                tk.emplace_back(c, (long)pi);
                for (size_t j = tk.size() - 1; j > 0 && (tk[j].first < tk[j-1].first ||
                     (tk[j].first == tk[j-1].first && tk[j].second < tk[j-1].second)); --j)
                  std::swap(tk[j], tk[j-1]);
              } else if (c > tk.front().first || (c == tk.front().first && (long)pi < tk.front().second)) {
                tk.front() = {c, (long)pi};
                for (size_t j = 1; j < tk.size() && (tk[j-1].first > tk[j].first ||
                     (tk[j-1].first == tk[j].first && tk[j-1].second > tk[j].second)); ++j)
                  std::swap(tk[j-1], tk[j]);
              }
            });
          }
          // fan-out histogram (before capping)
          vector<uint64_t> h(9, 0);
          for (size_t fi = 0; fi < NF; ++fi) {
            uint32_t n = nfan[fi];
            if (n == 0) continue;
            size_t b = n <= 5 ? n - 1 : n <= 10 ? 5 : n <= 25 ? 6 : n <= 50 ? 7 : 8;
            ++h[b];
          }
          // invert to per-precursor fragment lists (ordered by pi -> deterministic)
          map<long, vector<size_t>> won;
          for (size_t fi = 0; fi < NF; ++fi) for (auto& cp : topk[fi]) won[cp.second].push_back(fi);
          for (auto& kv2 : won)
          {
            const Precursor_& pc = precursors[kv2.first];
            vector<pair<double, double>> frags;
            vector<double> frag_scores;
            for (size_t fi : kv2.second)
            {
              // this fragment's corr for THIS precursor (from its kept top-k list)
              double c = -1.0;
              for (auto& cp : topk[fi]) if (cp.second == kv2.first) { c = cp.first; break; }
              const auto [inten, emit_inten] = weighted_(frag_traces[fi].intensity, c, frag_im[fi] - pc.im);   // [E5]
              frags.emplace_back(exportMz_(frag_traces[fi].tof, frag_traces[fi].b, frag_traces[fi].mz), emit_inten);
              frag_scores.push_back(c * inten);
            }
            MSSpectrum ms2;
            assembleFromList_(pc, win_lo, win_hi, frags, frag_scores, min_frags, max_frags, ms2);
            if (!ms2.empty()) bucket.push_back(std::move(ms2));
          }
          #pragma omp critical
          for (size_t b = 0; b < 9; ++b) fanout_hist[b] += h[b];
        }
      }
      catch (...)
      {
        #pragma omp critical
        if (!worker_err) worker_err = std::current_exception();
      }
      }();   // window body
      }      // omp task
    }
    if (worker_err) std::rethrow_exception(worker_err);
    // MEMORY: every MS1 structure is dead once the window loop has joined -- scoring was its only
    // reader (assembleOne_ -> scoreCandidates_ dereferences ms1_traces[pc.trace_idx] and its span).
    // They used to ride to the end of main_, i.e. through the canonical sort, the optional merge
    // and consolidate passes, the PeakMap assembly and the whole mzML write. [mem]
    vector<Trace>().swap(ms1_traces);
    ms1_store = TraceStore();
    vector<Precursor_>().swap(precursors);
    vector<double>().swap(prec_mz);
    writeLogInfo_("[mem] MS1 traces, spans, precursors released after the window loop." + rss_() + mem_());

    if (rp_max > 0)
    {
      uint64_t tot = 0; for (auto v : fanout_hist) tot += v;
      if (tot > 0)
      {
        const char* lbl[9] = {"1","2","3","4","5","6-10","11-25","26-50","51+"};
        String msg = "Fan-out histogram (precursors per fragment, before rp_max=" + String(rp_max) + " cap): ";
        for (size_t b = 0; b < 9; ++b)
          msg += String(lbl[b]) + "=" + String(fanout_hist[b]) + "(" + String((int)(1000.0 * fanout_hist[b] / tot) / 10.0) + "%) ";
        writeLogInfo_(msg);
      }
    }

    { auto& st = phase_stats_()["WINDOW_LOOP"]; if (st.n++ == 0) phase_order_().push_back("WINDOW_LOOP");
      st.wall += phase_clock_() - _t_win; st.cpu += cpu_seconds_() - _c_win; st.rss_end_mb = rss_mb_(); }
    writeLogInfo_("All windows done." + clk_() + rss_() + mem_()); // parallel-phase peak? [mem]
    {
      // These are CPU-seconds summed across threads, not wall time -- a stage that runs on 20
      // threads contributes 20x its wall. That is the right denominator for "where does the work
      // go", but it must not be read as wall time, so both are printed.
      const double S = 1e-9;
      const double tot = S * (t_mat + t_trace + t_grid + t_score + t_emit);
      auto pct = [&](long long v) { return tot > 0 ? 100.0 * S * v / tot : 0.0; };
      writeLogInfo_("Window-loop stage breakdown (each stage's own ELAPSED time, summed over "
                    "windows -- NOT CPU-seconds: stages of different windows run concurrently, so "
                    "this sum exceeds the loop's wall time; " + String(tot) + " s total):");
      writeLogInfo_("  materialize " + String(S * t_mat)   + " s (" + String(pct(t_mat))   + "%)");
      writeLogInfo_("  trace       " + String(S * t_trace) + " s (" + String(pct(t_trace)) + "%)");
      writeLogInfo_("  fraggrid    " + String(S * t_grid)  + " s (" + String(pct(t_grid))  + "%)");
      writeLogInfo_("  score       " + String(S * t_score) + " s (" + String(pct(t_score)) + "%)");
      writeLogInfo_("  emit        " + String(S * t_emit)  + " s (" + String(pct(t_emit))  + "%)");
      // Is the score stage doing arithmetic, or walking fragments the RT gate then discards?
      // One atomic add per precursor, so this costs nothing against a 10^11-iteration loop.
      const long long np = stat_prec_.load(), sl = stat_slice_.load(), rp = stat_rtpass_.load();
      if (np > 0)
        writeLogInfo_("  score gate: " + String(np) + " precursors visited " + String(sl)
                      + " fragments (" + String(sl / np) + "/precursor), "
                      + String(rp) + " passed the RT gate ("
                      + String(rp > 0 ? (double)sl / (double)rp : 0.0) + "x scanned per survivor)");
    }

    vector<MSSpectrum> all_out;
    {
      size_t n_tot = 0;
      for (const auto& bucket : win_out) n_tot += bucket.size();
      all_out.reserve(n_tot); // avoid repeated reallocation of ~1M spectrum objects
    }
    for (auto& bucket : win_out) { for (auto& s : bucket) all_out.push_back(std::move(s)); bucket.clear(); bucket.shrink_to_fit(); }

    // Canonical total-order sort (RT, precursor m/z, charge): thread-count-independent and free of
    // RT-tie ambiguity (RT-only sortSpectra is NOT canonical). [par-Crit-2]
    { Phase _ph("SORT(canonical)");
    sort(all_out.begin(), all_out.end(), [](const MSSpectrum& a, const MSSpectrum& b) {
      if (a.getRT() != b.getRT()) return a.getRT() < b.getRT();
      const double amz = a.getPrecursors().empty() ? 0.0 : a.getPrecursors()[0].getMZ();
      const double bmz = b.getPrecursors().empty() ? 0.0 : b.getPrecursors()[0].getMZ();
      if (amz != bmz) return amz < bmz;
      const int ac = a.getPrecursors().empty() ? 0 : a.getPrecursors()[0].getCharge();
      const int bc = b.getPrecursors().empty() ? 0 : b.getPrecursors()[0].getCharge();
      if (ac != bc) return ac < bc;
      // fragment-sequence tiebreaks so equal-precursor spectra have a canonical order [code-review]
      if (a.size() != b.size()) return a.size() < b.size();
      for (Size i = 0; i < a.size(); ++i)
      {
        if (a[i].getMZ() != b[i].getMZ()) return a[i].getMZ() < b[i].getMZ();
        if (a[i].getIntensity() != b[i].getIntensity()) return a[i].getIntensity() < b[i].getIntensity();
      }
      return false;
    }); }

    //-------------------------------------------------------------
    // [route-1] FEATURE-LEVEL CONSOLIDATION
    //
    // We emit 4.6 PSMs per unique peptide; the reference implementation emits 1.46. We are not producing WORSE spectra
    // (purity 8.3% vs their 4.0%) - we are producing the SAME peptide many times while missing others.
    // Emission is per precursor HYPOTHESIS (per trace), not per chromatographic FEATURE, so one real
    // feature can be emitted from several RT slices / IM sub-ranges / charge hypotheses.
    //
    // Consolidate: spectra whose precursors agree in m/z (ppm), RT and IM are ONE feature -> keep the
    // single best (most fragments, then most intense) and drop the rest. This changes emission
    // multiplicity WITHOUT touching fragment content, so it is orthogonal to every falsified
    // fragment-count knob. It should also relieve FDR competition among near-duplicate spectra.
    //
    // NOTE this is a DELETION, and deletion has repeatedly cost peptides here - hence default OFF and
    // a hard kill criterion (see TRACING.md): keep only if unique peptides rise >=3%.
    //-------------------------------------------------------------
    //-------------------------------------------------------------
    // [merge] Recombine SPLIT spectra into one joint spectrum -- union of peaks, intensities
    // SUMMED. This is not consolidate: consolidate SELECTS the richest member and discards the
    // rest, which measured -12.7% peptides. Measured basis for merging instead:
    //   * within-cycle pairs share a precursor (93.2% same charge, 7.5 ppm) and 41% of peaks
    //     are present in only ONE member; the union is 1.355x the larger member
    //   * ms1_split_valleys quadruples the affected groups (3,398 -> 14,379)
    // RT window: measured MS1 FWHM on dataset D is 3.61 s (2.61 cycles), flat across the gradient
    // (1.04x end to end), and our own spread is 2.66 cycles -- so ONE FWHM is the principled
    // window. It is 13x narrower than the 45.7 s at which an earlier analysis still saw
    // "benefit", which is the regime where two spectra cannot be one feature.
    //-------------------------------------------------------------
    const double merge_rt = getDoubleOption_("merge:rt_window");
    if (merge_rt > 0.0 && !all_out.empty())
    {
      const double mg_ppm = getDoubleOption_("merge:mass_ppm");
      const double mg_im  = getDoubleOption_("merge:delta_im");
      const double mg_mz  = getDoubleOption_("merge:mz_tol");
      const bool   mg_z   = !getFlag_("merge:any_charge");
      const bool   mg_sum = getFlag_("merge:sum_intensity");
      const double mg_cos = getDoubleOption_("merge:min_cosine");
      vector<MSSpectrum> kept;
      kept.reserve(all_out.size());
      vector<char> used(all_out.size(), 0);
      Size n_merged = 0, n_groups = 0;

      for (size_t i = 0; i < all_out.size(); ++i)
      {
        if (used[i]) continue;
        const auto& pi = all_out[i].getPrecursors();
        if (pi.empty()) { kept.push_back(std::move(all_out[i])); used[i] = 1; continue; }
        const double mz_i = pi[0].getMZ();
        const int    z_i  = pi[0].getCharge();
        const double im_i = pi[0].getDriftTime();

        vector<size_t> grp{i};
        // all_out is sorted by RT, so stop as soon as the window is exceeded
        for (size_t j = i + 1; j < all_out.size(); ++j)
        {
          if (all_out[j].getRT() - all_out[i].getRT() > merge_rt) break;
          if (used[j]) continue;
          const auto& pj = all_out[j].getPrecursors();
          if (pj.empty()) continue;
          if (fabs(pj[0].getMZ() - mz_i) > mz_i * mg_ppm * 1e-6) continue;
          if (mg_z && pj[0].getCharge() != z_i) continue;
          const double im_j = pj[0].getDriftTime();
          if (im_i > 0 && im_j > 0 && fabs(im_j - im_i) > mg_im) continue;
          // [merge-cos] CONTENT gate. Precursor coordinates are NOT the problem: measured on
          // identified pairs, SAME-peptide neighbours differ by a median 9.0 ppm and 0.0074
          // 1/K0 while DIFFERENT-peptide neighbours differ by 205,065 ppm and 0.0906, so the
          // coordinate gate already admits only 1.8% chimeras. The problem is that only 12,380
          // of 161,798 spectra identify at all, so a coordinate gate happily merges a good
          // spectrum with several spectra that carry no identifiable signal -- their peaks
          // dilute the matched fraction without contributing matches. Cosine sees that;
          // coordinates cannot. Measured reference distribution: within-cycle pairs 0.845,
          // same peptide across cycles 0.687, unrelated spectra 0.007.
          if (mg_cos > 0.0)
          {
            const MSSpectrum& A = all_out[i];
            const MSSpectrum& B = all_out[j];
            double dot = 0.0, na = 0.0, nb = 0.0;
            for (const auto& pa : A) na += (double)pa.getIntensity() * pa.getIntensity();
            for (const auto& pb : B) nb += (double)pb.getIntensity() * pb.getIntensity();
            if (na <= 0.0 || nb <= 0.0) continue;
            // both are m/z-sorted, so a merge-walk is O(|A|+|B|) rather than O(|A|*|B|)
            Size ia = 0, ib = 0;
            while (ia < A.size() && ib < B.size())
            {
              const double da = A[ia].getMZ() - B[ib].getMZ();
              if (fabs(da) <= mg_mz)
              { dot += (double)A[ia].getIntensity() * B[ib].getIntensity(); ++ia; ++ib; }
              else if (da < 0) ++ia;
              else ++ib;
            }
            if (dot / (sqrt(na) * sqrt(nb)) < mg_cos) continue;
          }
          grp.push_back(j);
        }
        if (grp.size() == 1) { kept.push_back(std::move(all_out[i])); used[i] = 1; continue; }

        // apex = highest TIC. Its RT/IM/precursor become the merged spectrum's coordinates, but
        // its CONTENT is the union of the whole group -- the point is not to pick the apex.
        size_t apex = grp[0];
        for (size_t g : grp) if (all_out[g].calculateTIC() > all_out[apex].calculateTIC()) apex = g;

        vector<pair<double, double>> pk;
        for (size_t g : grp)
          for (const auto& p : all_out[g]) pk.emplace_back(p.getMZ(), p.getIntensity());
        sort(pk.begin(), pk.end());
        MSSpectrum out = all_out[apex];
        out.clear(false);                       // keep metadata, drop peaks
        for (size_t a = 0; a < pk.size(); )
        {
          double mz_sum = pk[a].first * pk[a].second, in_sum = pk[a].second, in_max = pk[a].second;
          size_t b = a + 1;
          // combine every peak within mg_mz of the RUNNING centroid, so a dense cluster does not
          // chain arbitrarily far from where it started
          while (b < pk.size() && fabs(pk[b].first - mz_sum / in_sum) <= mg_mz)
          { mz_sum += pk[b].first * pk[b].second; in_sum += pk[b].second;
            in_max = std::max(in_max, pk[b].second); ++b; }
          Peak1D q; q.setMZ(mz_sum / in_sum);
          // [merge-fix] SUM systematically distorts the spectrum. Jaccard between members is
          // 0.374, so ~63% of peaks appear in only ONE member: summing gives shared peaks Nx and
          // unique peaks 1x, wrecking the fragment intensity ratios the search engine scores on
          // -- and penalising precisely the unique peaks that make merging worth doing. MAX
          // takes each ion's best single observation, which is what a more complete spectrum of
          // the same feature should look like.
          q.setIntensity(mg_sum ? in_sum : in_max);
          out.push_back(q);
          a = b;
        }
        // [merge-fix] RE-CAP. Merge runs AFTER assembleFromList_ has already capped each
        // spectrum at max_fragments, so merging N spectra yields up to N x max_fragments peaks
        // that are never re-ranked. ~47% of spectra sit AT the cap, so merged spectra routinely
        // land at 2-3x it -- far more random-match opportunities, diluting the score. Keep the
        // top max_frags by intensity, then restore m/z order (mzML requires ascending m/z).
        if (max_frags > 0 && out.size() > max_frags)
        {
          vector<Peak1D> v(out.begin(), out.end());
          std::nth_element(v.begin(), v.begin() + max_frags, v.end(),
                           [](const Peak1D& x, const Peak1D& y) { return x.getIntensity() > y.getIntensity(); });
          v.resize(max_frags);
          sort(v.begin(), v.end(), [](const Peak1D& x, const Peak1D& y) { return x.getMZ() < y.getMZ(); });
          out.clear(false);
          for (const auto& p : v) out.push_back(p);
        }
        for (size_t g : grp) used[g] = 1;
        n_merged += grp.size(); ++n_groups;
        kept.push_back(std::move(out));
      }
      writeLogInfo_("Spectrum merge: " + String(all_out.size()) + " -> " + String(kept.size())
                    + " spectra (" + String(n_groups) + " groups absorbed " + String(n_merged)
                    + " members, window " + String(merge_rt) + " s)" + rss_());
      all_out.swap(kept);
    }

    const double cons_rt = getDoubleOption_("consolidate:delta_rt");
    if (cons_rt > 0.0)
    {
      const double cons_ppm = getDoubleOption_("consolidate:mass_ppm");
      const double cons_im = getDoubleOption_("consolidate:delta_im");
      const bool cons_z = getFlag_("consolidate:same_charge_only");
      vector<MSSpectrum> kept;
      kept.reserve(all_out.size());
      vector<char> dropped(all_out.size(), 0);
      // all_out is already sorted by (RT, mz, charge); a feature is a run within cons_rt of each other
      for (size_t i = 0; i < all_out.size(); ++i)
      {
        if (dropped[i]) continue;
        size_t best = i;
        const auto& pi = all_out[i].getPrecursors();
        if (pi.empty()) { kept.push_back(std::move(all_out[i])); continue; }
        const double mz_i = pi[0].getMZ();
        const int z_i = pi[0].getCharge();
        // BUGFIX: IM is stamped on the PRECURSOR (:1137 prec.setDriftTime), not on the
        // MSSpectrum. getDriftTime() here returned -1 always, so "im_i > 0" was false and the
        // mobility gate NEVER FIRED -- every consolidate result to date merged on (RT, m/z,
        // charge) only.
        const double im_i = pi[0].getDriftTime();
        // scan forward while RT is still within the window (sorted by RT)
        for (size_t j = i + 1; j < all_out.size(); ++j)
        {
          if (all_out[j].getRT() - all_out[i].getRT() > cons_rt) break;
          if (dropped[j]) continue;
          const auto& pj = all_out[j].getPrecursors();
          if (pj.empty()) continue;
          if (fabs(pj[0].getMZ() - mz_i) > mz_i * cons_ppm * 1e-6) continue;
          if (cons_z && pj[0].getCharge() != z_i) continue;
          const double im_j = pj[0].getDriftTime();
          if (im_i > 0 && im_j > 0 && fabs(im_j - im_i) > cons_im) continue;
          // same feature: keep whichever carries more fragments (tie -> higher TIC)
          if (all_out[j].size() > all_out[best].size() ||
              (all_out[j].size() == all_out[best].size() &&
               all_out[j].calculateTIC() > all_out[best].calculateTIC()))
          {
            dropped[best] = 1; best = j;
          }
          else dropped[j] = 1;
        }
        kept.push_back(std::move(all_out[best]));
        dropped[best] = 1;
      }
      writeLogInfo_("Feature consolidation: " + String(all_out.size()) + " -> " + String(kept.size())
                    + " spectra (" + String((int)(1000.0 * kept.size() / std::max<size_t>(all_out.size(), 1)) / 10.0)
                    + "% kept)" + rss_());
      all_out.swap(kept);
    }

    PeakMap out_exp;
    { Phase _ph("ASSEMBLE(PeakMap)");
      for (auto& s : all_out) out_exp.addSpectrum(std::move(s)); }
    const Size n_out = out_exp.size();
    addDataProcessing_(out_exp, getProcessingInfo_(DataProcessing::DATA_PROCESSING));
    // [provenance] Record the emitted-intensity reweighting so an identical historical command line is
    // still reproducible after assembly:corr_power's default changed 0->2 (2026-07-28). Without this,
    // the same command produces different intensities across tool versions with no trace of why.
    // [provenance] WHICH m/z calibration produced this file. Without it an archived mzML is not
    // attributable to a calibration path and a silent fallback would be invisible afterwards.
    // The value travels ON THE LOADED EXPERIMENT (set by BrukerTimsFile), not through a header
    // static: an inline static can be duplicated across a shared-library boundary, in which case
    // the tool would read its own copy and always report "unset" [vibe/claude review 2026-09-01].
#ifdef SPEXTRACTOR_WITH_MZPEAK
    if (FileHandler::getTypeByFileName(in) == FileTypes::MZPEAK)
      out_exp.setMetaValue("spx:mz_calibration", spx::lastMzPeakCalibration());   // [codex #13] not BrukerTimsFile's state
    else
#endif
    out_exp.setMetaValue("spx:mz_calibration", BrukerTimsFile::lastMzCalibration());
    out_exp.setMetaValue("spx:detector", String(integer_detector ? "integer" : "openms"));
    out_exp.setMetaValue("spx:corr_power", corr_power_);
    out_exp.setMetaValue("spx:im_weight_sigma", im_weight_sigma_);
    out_exp.setMetaValue("spx:require_isotope_support",
                         (int)(getStringOption_("assembly:require_isotope_support") == "true"));
    // Output format. The extension is authoritative, because that is what a user typing
    // `-out pseudo.mzML` means; -out_type overrides it; mzPeak is the default when neither says.
    FileTypes::Type out_type = FileTypes::MZML;
    const String out_type_opt = getStringOption_("out_type");
#ifdef SPEXTRACTOR_WITH_MZPEAK
    if (!out_type_opt.empty())
      { String t = out_type_opt; t.toLower(); out_type = (t == "mzml") ? FileTypes::MZML : FileTypes::MZPEAK; }
    else
    {
      const FileTypes::Type by_name = FileHandler::getTypeByFileName(out);
      out_type = (by_name == FileTypes::MZML || by_name == FileTypes::MZPEAK) ? by_name
                                                                             : FileTypes::MZPEAK;
    }
    if (out_type == FileTypes::MZPEAK)
      writeLogInfo_("Writing mzPeak. NOTE: DDA search engines read mzML, not mzPeak -- if this file "
                    "is going straight into a search, write .mzML instead (or -out_type mzML).");
#else
    // Built without mzPeak support: mzML is the only container this binary can write. Refuse an
    // explicit request for mzPeak rather than silently writing something else.
    String out_type_lc = out_type_opt; out_type_lc.toLower();          // toLower() mutates
    if (!out_type_lc.empty() && out_type_lc != "mzml")
      throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "-out_type mzpeak needs a build with mzPeak support (-DMZPEAK_ROOT=...); this binary "
            "can only write mzML.");
#endif
    { Phase _ph(out_type == FileTypes::MZML ? "WRITE(mzML)" : "WRITE(mzPeak)");
      FileHandler().storeExperiment(out, out_exp, {out_type}, log_type_); }
    writeLogInfo_("Wrote " + String(n_out) + " pseudo-MS2 spectra to " + out
                  + " (" + String(FileTypes::typeToName(out_type)) + ")");
    report_phases_(phase_clock_());
    return EXECUTION_OK;
  }
};

/// @endcond

int main(int argc, const char** argv)
{
  TOPPSpeXtractor tool;
  return tool.main(argc, argv);
}
