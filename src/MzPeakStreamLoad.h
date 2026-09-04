// SpeXtract: streaming diaPASEF load from an .mzpeak archive through the mzPeak C++ library
// (github.com/OpenMS/mzpeak, fork okohlbacher/mzpeak-openms). Mirrors BrukerTimsFile::
// loadDIAStreaming's contract so PickCompactConsumer sees the same stream: every MS1 spectrum
// first (frame order), then one MSSpectrum per (MS2 frame, isolation window) holding only the
// peaks inside that window's 1/K0 band, handed over serially in frame order.
//
// mzPeak stores a diaPASEF frame as ONE spectrum carrying N precursors (one per window, each
// with an ion-mobility band) and a per-peak ion-mobility array; ims-compact archives carry raw
// TOF and the library reconstructs m/z with the archive's own two-point transform. Peaks arrive
// mobility-major, so each per-window spectrum is m/z-sorted here.
//
// ponytail: decode is parallel over contiguous frame ranges with one Index per thread (decodes
// serialise on a reader's mutex); the hand-off is serial. No caching beyond the library's own.
#pragma once
#ifdef SPEXTRACT_WITH_MZPEAK

#include <mzpeak/open.h>
#include <mzpeak/index.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/spectrum_metadata.h>

// [D1 exact m/z] the archive stores raw TOF behind a two-point transform; the exact ModelType-1 calibration
// lives only in the embedded vendor/analysis.tdf.gz. Recover tof = round((sqrt(mz) - c0) / c1) and re-apply
// TdfMzCalibration per frame. ON BY DEFAULT and fail-closed; SPEXTRACT_MZPEAK_EXACT=0 disables it and
// falls back to the archive's two-point transform, which costs ~12% of identified peptides.
#include <OpenMS/FORMAT/TdfMzCalibration.h>   // installed by the OpenMS patch (same header as src/TdfMzCalibration.h)
#include <zip.h>
#include <zlib.h>
#if __has_include(<SQLiteCpp/SQLiteCpp.h>)
#include <SQLiteCpp/SQLiteCpp.h>   // OpenMS' in-tree extern (src/openms/extern/SQLiteCpp); not installed by OpenMS
#define SPX_HAVE_SQLITECPP 1
#endif
#include <arrow/api.h>
#include <arrow/io/memory.h>
#include <parquet/file_reader.h>
#include <parquet/arrow/schema.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <cstring>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace spx
{
  struct MzPeakWin { double lo, hi, im_lo, im_hi; };

  /// Calibration provenance of the last mzPeak load (the archive's two-point tof transform, applied by the
  /// library; NOT the TDF table model). Stamped into the output instead of BrukerTimsFile's state.
  inline std::string& lastMzPeakCalibration() { static std::string s = "unset"; return s; }


  /// Exact TOF->m/z for an ims-compact archive: (c0, c1) from the peaks table's column metadata, the
  /// MzCalibration row and Frames.T1 from the embedded vendor tdf. Throws on anything missing (fail closed).
  struct MzPeakExactMz
  {
    double c0 = 0, c1 = 0;
    spextract::TdfMzCalibration cal;
    std::vector<double> t1_by_frame;          // index = Frames.Id
    bool enabled = false;

    static std::vector<char> readMember_(zip_t* z, const char* name, zip_int64_t from = 0, zip_int64_t len = -1)
    {
      zip_int64_t idx = zip_name_locate(z, name, 0);
      if (idx < 0) throw std::runtime_error(std::string("mzpeak: member missing: ") + name);
      zip_stat_t st; zip_stat_index(z, (zip_uint64_t)idx, 0, &st);
      if (len < 0) len = (zip_int64_t)st.size - from;
      zip_file_t* f = zip_fopen_index(z, (zip_uint64_t)idx, 0);
      if (!f) throw std::runtime_error(std::string("mzpeak: cannot open member ") + name);
      if (from > 0 && zip_fseek(f, from, SEEK_SET) != 0) { zip_fclose(f); throw std::runtime_error(std::string("mzpeak: cannot seek in ") + name + " (member not stored?)"); }
      std::vector<char> buf((size_t)len); zip_int64_t got = 0;
      while (got < len) { zip_int64_t n = zip_fread(f, buf.data() + got, (zip_uint64_t)(len - got)); if (n <= 0) break; got += n; }
      zip_fclose(f);
      if (got != len) throw std::runtime_error(std::string("mzpeak: short read of ") + name);
      return buf;
    }

    explicit MzPeakExactMz(const std::string& path)
    {
      int err = 0; zip_t* z = zip_open(path.c_str(), ZIP_RDONLY, &err);
      if (!z) throw std::runtime_error("mzpeak: zip_open failed on " + path);
      try
      {
        // (1) two-point transform params from the parquet footer of spectra_peaks.parquet (stored member)
        zip_int64_t idx = zip_name_locate(z, "spectra_peaks.parquet", 0);
        if (idx < 0) throw std::runtime_error("mzpeak: spectra_peaks.parquet missing");
        zip_stat_t st; zip_stat_index(z, (zip_uint64_t)idx, 0, &st);
        std::vector<char> tail = readMember_(z, "spectra_peaks.parquet", (zip_int64_t)st.size - 8, 8);
        uint32_t flen; std::memcpy(&flen, tail.data(), 4);
        if (std::string(tail.data() + 4, 4) != "PAR1") throw std::runtime_error("mzpeak: peaks member is not parquet");
        std::vector<char> footer = readMember_(z, "spectra_peaks.parquet", (zip_int64_t)st.size - 8 - flen, flen);
        uint32_t flen2 = flen;
        std::shared_ptr<parquet::FileMetaData> md = parquet::FileMetaData::Make(footer.data(), &flen2);
        std::shared_ptr<arrow::Schema> schema;
        parquet::ArrowReaderProperties props;
        PARQUET_THROW_NOT_OK(parquet::arrow::FromParquetSchema(md->schema(), props, md->key_value_metadata(), &schema));
        auto point = schema->GetFieldByName("point");
        if (!point) throw std::runtime_error("mzpeak: no 'point' column");
        auto tof = std::dynamic_pointer_cast<arrow::StructType>(point->type())->GetFieldByName("tof");
        if (!tof || !tof->metadata()) throw std::runtime_error("mzpeak: no tof column / metadata (not an ims-compact archive?)");
        auto tp = tof->metadata()->Get("mzpeak:transform_params");
        if (!tp.ok()) throw std::runtime_error("mzpeak: tof column has no mzpeak:transform_params");
        if (std::sscanf(tp->c_str(), "%lf,%lf", &c0, &c1) != 2 || !(c1 > 0)) throw std::runtime_error("mzpeak: bad transform_params " + *tp);

        // (2) the vendor tdf: gunzip vendor/analysis.tdf.gz to a temp file, read MzCalibration + Frames.T1
        // --no-vendor archives have no embedded tdf: accept a sidecar via SPEXTRACT_MZPEAK_TDF=<analysis.tdf.gz>
        std::vector<char> gz;
        if (const char* side = std::getenv("SPEXTRACT_MZPEAK_TDF"))
        {
          FILE* sf = std::fopen(side, "rb"); if (!sf) throw std::runtime_error(std::string("mzpeak: cannot read SPEXTRACT_MZPEAK_TDF ") + side);
          char b[1 << 16]; size_t n; while ((n = std::fread(b, 1, sizeof b, sf)) > 0) gz.insert(gz.end(), b, b + n); std::fclose(sf);
        }
        else gz = readMember_(z, "vendor/analysis.tdf.gz");
        char tmpl[] = "/tmp/spx_mzpeak_tdf_XXXXXX"; int fd = mkstemp(tmpl); if (fd < 0) throw std::runtime_error("mzpeak: mkstemp");
        close(fd);
        { z_stream zs{}; inflateInit2(&zs, 16 + MAX_WBITS); FILE* out = std::fopen(tmpl, "wb");
          zs.next_in = (Bytef*)gz.data(); zs.avail_in = (uInt)gz.size(); std::vector<char> ob(1 << 20); int rc;
          do { zs.next_out = (Bytef*)ob.data(); zs.avail_out = (uInt)ob.size(); rc = inflate(&zs, Z_NO_FLUSH);
               if (rc != Z_OK && rc != Z_STREAM_END) { std::fclose(out); inflateEnd(&zs); throw std::runtime_error("mzpeak: gunzip of analysis.tdf.gz failed"); }
               std::fwrite(ob.data(), 1, ob.size() - zs.avail_out, out); } while (rc != Z_STREAM_END);
          inflateEnd(&zs); std::fclose(out); }
#ifndef SPX_HAVE_SQLITECPP
        throw std::runtime_error("mzpeak exact m/z needs SQLiteCpp headers at build time (OpenMS extern); rebuild inside the OpenMS tree");
#else
        {
          SQLite::Database db(std::string(tmpl), SQLite::OPEN_READONLY);
          SQLite::Statement q(db, "SELECT ModelType, DigitizerTimebase, DigitizerDelay, C0, C1, C2, T1, dC1, dC2, C3, C4 FROM MzCalibration");
          if (!q.executeStep()) throw std::runtime_error("mzpeak: no MzCalibration row in the embedded tdf");
          cal.model_type = q.getColumn(0).getInt(); cal.digitizer_timebase = q.getColumn(1).getDouble(); cal.digitizer_delay = q.getColumn(2).getDouble();
          cal.C0 = q.getColumn(3).getDouble(); cal.C1 = q.getColumn(4).getDouble(); cal.C2 = q.getColumn(5).getDouble(); cal.T1_ref = q.getColumn(6).getDouble();
          cal.dC1 = q.getColumn(7).getDouble(); cal.dC2 = q.getColumn(8).getDouble(); cal.C3 = q.getColumn(9).getDouble(); cal.C4 = q.getColumn(10).getDouble();
          if (q.executeStep()) throw std::runtime_error("mzpeak: more than one MzCalibration row");
          if (!cal.isSupported()) throw std::runtime_error("mzpeak: " + cal.unsupportedReason());
          SQLite::Statement f(db, "SELECT Id, T1 FROM Frames");
          while (f.executeStep())
          {
            const long id = f.getColumn(0).getInt64(); if (id < 0 || id > 10000000) throw std::runtime_error("mzpeak: implausible Frames.Id");
            if ((size_t)id >= t1_by_frame.size()) t1_by_frame.resize((size_t)id + 1, cal.T1_ref);
            t1_by_frame[(size_t)id] = f.getColumn(1).isNull() ? cal.T1_ref : f.getColumn(1).getDouble();
          }
        }
#endif
        std::remove(tmpl);
        enabled = true;
      }
      catch (...) { zip_close(z); throw; }
      zip_close(z);
    }

    /// Recover the integer tof from the archive's two-point m/z and re-apply the exact model.
    void exact(std::vector<double>& mz, long frame_id) const
    {
      // Fail rather than fall back to the reference T1: substituting it for an unknown frame stamps
      // masses calibrated at the wrong temperature as "exact", which is silent wrongness the
      // provenance record would then vouch for. [adv-review codex 2026-09-03]
      if (frame_id < 0 || (size_t)frame_id >= t1_by_frame.size())
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "mzPeak exact calibration: no tdf T1 for frame", std::to_string(frame_id));
      const double b = cal.frameFactor(t1_by_frame[(size_t)frame_id]);
      for (double& v : mz)
      {
        if (!std::isfinite(v) || v <= 0.0)
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "mzPeak exact calibration: non-positive archive m/z", std::to_string(v));
        const double tof = std::llround((std::sqrt(v) - c0) / c1);
        const double e = cal.tofToMz(tof, b);
        if (!std::isfinite(e))
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "mzPeak exact calibration produced a non-finite m/z at tof", std::to_string(tof));
        v = e;
      }
    }
  };

  /// One decoded frame split per window (index-aligned with `wins` of that frame).
  inline void frameToSpectra_(const MzPeak::Spectrum& sp, double rt, std::size_t frame_idx,
                              const std::vector<MzPeakWin>& wins, std::vector<OpenMS::MSSpectrum>& out,
                              const MzPeakExactMz* exact = nullptr, long frame_id = -1)
  {
    using namespace OpenMS;
    std::vector<double> mz_exact;
    if (exact && exact->enabled) { mz_exact = sp.mz(); exact->exact(mz_exact, frame_id); }
    const std::vector<double>& mz = (exact && exact->enabled) ? mz_exact : sp.mz();
    const std::vector<float>& in = sp.intensity();
    const std::vector<double>& im = sp.ion_mobility_array();
    // Array lengths must agree. Truncating to the shorter one turns a malformed archive into a
    // smaller frame that looks valid, and the peaks that vanish are never reported anywhere.
    if (mz.size() != in.size())
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "mzPeak frame has mismatched m/z and intensity array lengths",
            std::to_string(mz.size()) + " vs " + std::to_string(in.size()));
    const std::size_t n = mz.size();
    const bool have_im = im.empty() || im.size() >= n;
    if (!im.empty() && im.size() < n)
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "mzPeak frame has a short ion-mobility array", std::to_string(im.size()));
    out.assign(wins.size(), MSSpectrum());
    for (std::size_t w = 0; w < wins.size(); ++w)
    {
      const MzPeakWin& win = wins[w];
      MSSpectrum& spec = out[w];
      spec.setRT(rt);
      spec.setMSLevel(wins.size() == 1 && !(win.lo < win.hi) ? 1 : 2);
      spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      spec.setNativeID("mzpeak=" + std::to_string(frame_idx) + " window=" + std::to_string(w));
      if (spec.getMSLevel() == 2)
      {
        Precursor prec;
        prec.setMZ(0.5 * (win.lo + win.hi));
        prec.setIsolationWindowLowerOffset(0.5 * (win.hi - win.lo));
        prec.setIsolationWindowUpperOffset(0.5 * (win.hi - win.lo));
        spec.setPrecursors({prec});
        spec.setMetaValue("ion mobility lower limit", win.im_lo);
        spec.setMetaValue("ion mobility upper limit", win.im_hi);
      }
      DataArrays::FloatDataArray im_array;
      IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
      spec.reserve(n / wins.size() + 64); im_array.reserve(n / wins.size() + 64);
      for (std::size_t k = 0; k < n; ++k)
      {
        const double k0 = have_im ? im[k] : 0.0;
        if (have_im && (k0 < win.im_lo || k0 > win.im_hi)) continue;
        if (!(in[k] > 0)) continue;
        spec.emplace_back(mz[k], in[k]);
        im_array.push_back(static_cast<float>(k0));
      }
      if (spec.empty()) continue;
      // mobility-major on disk -> sort by m/z, carrying the IM array along
      std::vector<std::size_t> ord(spec.size());
      for (std::size_t k = 0; k < ord.size(); ++k) ord[k] = k;
      std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return spec[a].getMZ() < spec[b].getMZ(); });
      std::vector<Peak1D> sorted; sorted.reserve(spec.size());
      DataArrays::FloatDataArray im_sorted; IMDataConverter::setIMUnit(im_sorted, DriftTimeUnit::VSSC); im_sorted.reserve(spec.size());
      for (std::size_t k : ord) { sorted.push_back(spec[k]); im_sorted.push_back(im_array[k]); }
      spec.clear(false);                                          // peaks only; RT/precursor/meta stay
      for (const Peak1D& pk : sorted) spec.push_back(pk);
      spec.getFloatDataArrays().push_back(std::move(im_sorted));
      spec.setIMPeakType(IMPeakType::IM_PROFILE);
    }
  }

  /// Stream an .mzpeak run into `consumer` (MS1 first, then MS2 per (frame, window), frame order).
  /// Returns the number of frames seen.
  inline std::size_t loadMzPeakStreaming(const std::string& path, OpenMS::FullSwathFileConsumer& consumer, int threads)
  {
    // 1) metadata sweep (no peak decode): RT, level, windows with IM bands, per frame
    MzPeak::Index index = MzPeak::open(path);
    MzPeak::Spectra spectra = index.spectra();
    const std::size_t n = spectra.size();
    lastMzPeakCalibration() = "mzpeak_two_point_transform (library-applied MS:1003825 from the archive; not the TDF MzCalibration model) archive=" + path.substr(path.find_last_of('/') + 1);
    std::vector<double> rt(n);
    std::vector<long> frame_id(n, -1);
    std::vector<std::vector<MzPeakWin>> wins(n);
    std::vector<std::size_t> ms1, ms2;
    std::unique_ptr<MzPeakExactMz> exact;
    // DEFAULT ON (2026-09-02 18:58): with the exact model recovered from the tdf, mzPeak input gives 12,082 Sage
    // peptides vs 10,785 with the archive's two-point transform (and 12,217 from .d): the transform was 91% of the
    // loss. Fail closed like the .d path: no tdf (embedded or SPEXTRACT_MZPEAK_TDF sidecar) -> error, unless
    // SPEXTRACT_MZPEAK_EXACT=0 explicitly accepts the two-point m/z.
    const char* ex = std::getenv("SPEXTRACT_MZPEAK_EXACT");
    if (!(ex && std::string(ex) == "0"))
    {
      try { exact = std::make_unique<MzPeakExactMz>(path); }
      catch (const std::exception& e)
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          std::string("mzPeak input: cannot recover the exact TDF calibration (") + e.what() + "). The archive's two-point transform is ~7 ppm off and costs ~12% peptides; set SPEXTRACT_MZPEAK_TDF=<analysis.tdf.gz> or SPEXTRACT_MZPEAK_EXACT=0 to accept it.", path);
      }
      lastMzPeakCalibration() = "tdf_table_modeltype1 (recovered from the archive's two-point transform + embedded vendor/analysis.tdf.gz) archive=" + path.substr(path.find_last_of('/') + 1);
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      MzPeak::Spectrum sp = spectra[i];
      const auto& md = sp.metadata();
      rt[i] = md.retention_time.value_or(0.0);
      { const std::string& id = md.id; const auto pos = id.find("frame="); if (pos != std::string::npos) frame_id[i] = std::atol(id.c_str() + pos + 6); }
      if (md.ms_level.value_or(1) == 1) { wins[i] = {MzPeakWin{0, 0, -1e9, 1e9}}; ms1.push_back(i); continue; }
      // Selected ions carry the per-window 1/K0 band. mzpeak-convert 0.9.x writes precursor_index
      // as NULL and the band as NAME-ONLY CV params; the reader then attaches every ion of the
      // frame to the FIRST precursor. So: flatten all ions of the spectrum and match each window
      // to the ion whose selected_ion_mz is its isolation target.
      struct Ion { double mz, lo, hi; };
      std::vector<Ion> ions;
      for (const auto& p : md.precursors)
        for (const auto& si : p.selected_ions)
        {
          Ion ion{si.selected_ion_mz.value_or(-1.0), -1e9, 1e9};
          if (si.ion_mobility_lower_limit && si.ion_mobility_upper_limit) { ion.lo = *si.ion_mobility_lower_limit; ion.hi = *si.ion_mobility_upper_limit; }
          for (const auto& cv : si.parameters)
          {
            if (!cv.name || !cv.value) continue;
            if (cv.name->find("ion mobility lower limit") != std::string::npos) ion.lo = std::stod(*cv.value);
            else if (cv.name->find("ion mobility upper limit") != std::string::npos) ion.hi = std::stod(*cv.value);
          }
          ions.push_back(ion);
        }
      for (const auto& p : md.precursors)
      {
        if (!p.isolation_window.target_mz) continue;
        const double c = *p.isolation_window.target_mz;
        const double lo = c - p.isolation_window.lower_offset.value_or(0.f), hi = c + p.isolation_window.upper_offset.value_or(0.f);
        const Ion* best = nullptr;
        for (const Ion& ion : ions) if (!best || std::fabs(ion.mz - c) < std::fabs(best->mz - c)) best = &ion;
        if (!best || std::fabs(best->mz - c) > 0.05 || best->lo < -1e8 || best->hi > 1e8)
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "mzPeak MS2 window without a matching ion-mobility band (spectrum " + std::to_string(i) + ")", std::to_string(c));
        wins[i].push_back(MzPeakWin{lo, hi, best->lo, best->hi});
      }
      if (!wins[i].empty()) ms2.push_back(i);
    }

    // 2) decode in parallel over contiguous ranges (row-group locality), hand off serially in order
    int nthr = std::max(1, threads);
#ifdef _OPENMP
    nthr = std::min(nthr, std::max(1, omp_get_max_threads()));
#endif
    auto run = [&](const std::vector<std::size_t>& ids)
    {
      const std::size_t per = 12;                                  // frames per thread per batch (~1 row group)
      const std::size_t batch = std::max<std::size_t>(per, per * (std::size_t)nthr);
      for (std::size_t b0 = 0; b0 < ids.size(); b0 += batch)
      {
        const std::size_t b1 = std::min(ids.size(), b0 + batch);
        const int nchunk = (int)((b1 - b0 + per - 1) / per);
        std::vector<std::vector<std::vector<OpenMS::MSSpectrum>>> out((std::size_t)nchunk);
        std::exception_ptr err;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nthr)
        for (int c = 0; c < nchunk; ++c)
        {
          try
          {
            const std::size_t c0 = b0 + (std::size_t)c * per, c1 = std::min(b1, c0 + per);
            // one reader per THREAD (not per work item): MzPeak::open re-reads the zip directory and
            // every parquet footer; ~1,450 opens cost ~1 h of system time on dataset D (measured 13:55).
            static thread_local std::unique_ptr<MzPeak::Index> tidx;   // Index is neither movable nor copyable:
            static thread_local std::string tidx_path;                  // construct it from the prvalue (elided)
            if (!tidx || tidx_path != path) { tidx.reset(new MzPeak::Index(MzPeak::open(path))); tidx_path = path; }
            MzPeak::Spectra tsp = tidx->spectra();
            std::vector<std::size_t> want(ids.begin() + (std::ptrdiff_t)c0, ids.begin() + (std::ptrdiff_t)c1);
            std::vector<MzPeak::Spectrum> got = tsp.get_spectra_batch(want);
            // A short batch result used to be clamped away with min(), leaving the tail frames as
            // empty spectra: fewer peaks, no warning, and a run that still reports success.
            // [adv-review codex 2026-09-03]
            if (got.size() != want.size())
              throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "mzPeak reader returned fewer frames than requested",
                    std::to_string(got.size()) + " of " + std::to_string(want.size()));
            auto& slot = out[(std::size_t)c]; slot.resize(want.size());
            for (std::size_t k = 0; k < want.size(); ++k)
              frameToSpectra_(got[k], rt[want[k]], want[k], wins[want[k]], slot[k], exact.get(), frame_id[want[k]]);
          }
          catch (...)
          {
            #pragma omp critical(mzpeak_err)
            if (!err) err = std::current_exception();
          }
        }
        if (err) std::rethrow_exception(err);
        for (auto& chunk : out) for (auto& frame : chunk) for (auto& spec : frame)
          if (!spec.empty()) consumer.consumeSpectrum(spec);
      }
    };
    run(ms1);
    run(ms2);
    return ms1.size() + ms2.size();
  }
} // namespace spx
#endif // SPEXTRACT_WITH_MZPEAK
