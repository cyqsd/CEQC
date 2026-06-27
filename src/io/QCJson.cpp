#include "ceqc/io/QCJson.hpp"
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ceqc::view {
namespace {

class json {
public:
  enum class Type { Null, Bool, Number, String, Object, Array };

  json() = default;
  json(std::nullptr_t) : type_(Type::Null) {}
  json(bool v) : type_(Type::Bool), bool_(v) {}
  json(const char* v) : type_(Type::String), string_(v ? v : "") {}
  json(const std::string& v) : type_(Type::String), string_(v) {}
  json(std::string&& v) : type_(Type::String), string_(std::move(v)) {}

  template <typename T, typename = std::enable_if_t<std::is_integral_v<T> && !std::is_same_v<std::decay_t<T>, bool>>>
  json(T v) : type_(Type::Number), number_(std::to_string(static_cast<long long>(v))) {}

  json(double v) : type_(Type::Number), number_(formatNumber(v)) {}

  static json object() {
    json j;
    j.type_ = Type::Object;
    return j;
  }

  static json array() {
    json j;
    j.type_ = Type::Array;
    return j;
  }


  json& operator[](const std::string& key) {
    ensureObject();
    return object_[key];
  }

  json& operator[](const char* key) {
    ensureObject();
    return object_[key ? key : ""];
  }

  void push_back(json v) {
    ensureArray();
    array_.push_back(std::move(v));
  }

  std::string dump(int indent = -1) const {
    std::ostringstream os;
    dumpTo(os, indent, 0);
    return os.str();
  }

private:
  Type type_ = Type::Null;
  bool bool_ = false;
  std::string number_;
  std::string string_;
  std::map<std::string, json> object_;
  std::vector<json> array_;

  static std::string formatNumber(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream os;
    os << std::setprecision(15) << v;
    return os.str();
  }

  static std::string escape(const std::string& s) {
    std::ostringstream os;
    os << '"';
    for (unsigned char c : s) {
      switch (c) {
      case '"': os << "\\\""; break;
      case '\\': os << "\\\\"; break;
      case '\b': os << "\\b"; break;
      case '\f': os << "\\f"; break;
      case '\n': os << "\\n"; break;
      case '\r': os << "\\r"; break;
      case '\t': os << "\\t"; break;
      default:
        if (c < 0x20) {
          os << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c)
             << std::dec << std::setfill(' ');
        } else {
          os << static_cast<char>(c);
        }
      }
    }
    os << '"';
    return os.str();
  }

  static void writeIndent(std::ostream& os, int count) {
    for (int i = 0; i < count; ++i) os.put(' ');
  }

  void ensureObject() {
    if (type_ != Type::Object) {
      type_ = Type::Object;
      object_.clear();
      array_.clear();
      string_.clear();
      number_.clear();
    }
  }

  void ensureArray() {
    if (type_ != Type::Array) {
      type_ = Type::Array;
      object_.clear();
      array_.clear();
      string_.clear();
      number_.clear();
    }
  }

  void dumpTo(std::ostream& os, int indent, int depth) const {
    switch (type_) {
    case Type::Null:
      os << "null";
      return;
    case Type::Bool:
      os << (bool_ ? "true" : "false");
      return;
    case Type::Number:
      os << (number_ == "null" ? "null" : number_);
      return;
    case Type::String:
      os << escape(string_);
      return;
    case Type::Object: {
      os << '{';
      bool first = true;
      for (const auto& kv : object_) {
        if (!first) os << ',';
        if (indent >= 0) {
          os << '\n';
          writeIndent(os, depth + indent);
        }
        os << escape(kv.first) << ':';
        if (indent >= 0) os << ' ';
        kv.second.dumpTo(os, indent, depth + indent);
        first = false;
      }
      if (!object_.empty() && indent >= 0) {
        os << '\n';
        writeIndent(os, depth);
      }
      os << '}';
      return;
    }
    case Type::Array: {
      os << '[';
      bool first = true;
      for (const auto& v : array_) {
        if (!first) os << ',';
        if (indent >= 0) {
          os << '\n';
          writeIndent(os, depth + indent);
        }
        v.dumpTo(os, indent, depth + indent);
        first = false;
      }
      if (!array_.empty() && indent >= 0) {
        os << '\n';
        writeIndent(os, depth);
      }
      os << ']';
      return;
    }
    }
  }
};

json statsJson(const ceqc::model::QCMetricStats& st) {
  json j = json::object();
  j["count"] = st.count;
  j["mean"] = st.mean;
  j["rms"] = st.rms;
  j["min"] = st.min;
  j["max"] = st.max;
  j["low_count"] = st.lowCount;
  j["jumps"] = st.jumps;
  return j;
}

json intArray(const std::vector<int>& v) {
  json a = json::array();
  for (int x : v) a.push_back(x);
  return a;
}

json stringArray(const std::vector<std::string>& v) {
  json a = json::array();
  for (const auto& x : v) a.push_back(x);
  return a;
}

json xyzArray(const double xyz[3]) {
  json a = json::array();
  a.push_back(xyz[0]);
  a.push_back(xyz[1]);
  a.push_back(xyz[2]);
  return a;
}

json mapStringInt(const std::map<std::string,int>& m) {
  json j = json::object();
  for (const auto& kv : m) j[kv.first] = kv.second;
  return j;
}

json mapStringDouble(const std::map<std::string,double>& m) {
  json j = json::object();
  for (const auto& kv : m) j[kv.first] = kv.second;
  return j;
}

json statMapObject(const std::map<std::string, ceqc::model::QCMetricStats>& m) {
  json j = json::object();
  for (const auto& kv : m) j[kv.first] = statsJson(kv.second);
  return j;
}

json statMapArray(const std::map<std::string, ceqc::model::QCMetricStats>& m, const std::string& keyName = "name") {
  json a = json::array();
  for (const auto& kv : m) {
    json row = statsJson(kv.second);
    row[keyName] = kv.first;
    a.push_back(row);
  }
  return a;
}

std::string systemName(const std::string& satOrSystem) {
  char s = satOrSystem.empty() ? '?' : satOrSystem[0];
  switch (s) {
  case 'G': return "GPS";
  case 'R': return "GLONASS";
  case 'E': return "Galileo";
  case 'C': return "BeiDou";
  case 'J': return "QZSS";
  case 'S': return "SBAS";
  case 'I': return "NavIC";
  default: return satOrSystem;
  }
}

std::string obsCodeSystem(const std::string& code) {
  auto p = code.find(':');
  return p == std::string::npos ? std::string{} : code.substr(0, p);
}

std::string obsCodeName(const std::string& code) {
  auto p = code.find(':');
  return p == std::string::npos ? code : code.substr(p + 1);
}

json timeplotRows(const ceqc::model::QCDerivedSummary& d) {
  json a = json::array();
  for (const auto& kv : d.satelliteTimeplot) {
    json r = json::object();
    r["satellite"] = kv.first;
    r["system"] = systemName(kv.first);
    r["plot"] = kv.second;
    a.push_back(r);
  }
  return a;
}

json skyplotArray(const ceqc::model::QCDerivedSummary& d) {
  json a = json::array();
  for (const auto& kv : d.satelliteMaxElevationDeg) {
    auto az = d.satelliteMaxElevationAzimuthDeg.find(kv.first);
    json r = json::object();
    r["satellite"] = kv.first;
    r["system"] = systemName(kv.first);
    r["max_elevation_deg"] = kv.second;
    if (az != d.satelliteMaxElevationAzimuthDeg.end()) r["azimuth_at_max_elevation_deg"] = az->second;
    else r["azimuth_at_max_elevation_deg"] = nullptr;
    a.push_back(r);
  }
  return a;
}

json riseSetArray(const std::vector<ceqc::model::QCRiseSetEvent>& evs) {
  json a = json::array();
  for (const auto& ev : evs) {
    json r = json::object();
    r["satellite"] = ev.satellite;
    r["system"] = systemName(ev.satellite);
    r["first"] = ev.first;
    r["last"] = ev.last;
    r["duration_hours"] = ev.durationHours;
    r["observation_count"] = ev.obsCount;
    r["has_ephemeris"] = ev.hasEphemeris;
    if (std::isfinite(ev.maxElevationDeg)) r["max_elevation_deg"] = ev.maxElevationDeg;
    else r["max_elevation_deg"] = nullptr;
    a.push_back(r);
  }
  return a;
}

json epochPositions(const std::vector<ceqc::model::QCEpochPosition>& eps) {
  json a = json::array();
  for (const auto& ep : eps) {
    json r = json::object();
    r["time"] = ep.time;
    r["status"] = ep.status;
    r["used_svs"] = ep.usedSVs;
    r["x_m"] = ep.x;
    r["y_m"] = ep.y;
    r["z_m"] = ep.z;
    r["latitude_deg"] = ep.latDeg;
    r["longitude_deg"] = ep.lonDeg;
    r["height_m"] = ep.heightM;
    r["clock_bias_m"] = ep.clockBiasM;
    r["postfit_rms_m"] = ep.postfitRmsM;
    r["max_residual_m"] = ep.maxResidualM;
    a.push_back(r);
  }
  return a;
}

json snrChartArray(const std::map<std::string, ceqc::model::QCMetricStats>& snr) {
  json a = json::array();
  for (const auto& kv : snr) {
    json r = statsJson(kv.second);
    r["series"] = kv.first;
    r["system"] = obsCodeSystem(kv.first);
    r["system_name"] = systemName(obsCodeSystem(kv.first));
    r["observation_code"] = obsCodeName(kv.first);
    a.push_back(r);
  }
  return a;
}

json multipathChartArray(const ceqc::model::QCDerivedSummary& d) {
  json a = json::array();
  auto add = [&](const std::string& name, double rms, int samples) {
    json r = json::object();
    r["combination"] = name;
    r["rms_m"] = rms;
    r["rms_cm"] = rms * 100.0;
    r["samples"] = samples;
    a.push_back(r);
  };
  int mp1n = 0, mp2n = 0;
  auto it1 = d.multipathMovingCount.find("MP1"); if (it1 != d.multipathMovingCount.end()) mp1n = it1->second;
  auto it2 = d.multipathMovingCount.find("MP2"); if (it2 != d.multipathMovingCount.end()) mp2n = it2->second;
  if (mp1n > 0) add("MP1", d.mp1Meters, mp1n);
  if (mp2n > 0) add("MP2", d.mp2Meters, mp2n);
  for (const auto& kv : d.multipathMovingRMS) {
    if (kv.first == "MP1" || kv.first == "MP2") continue;
    int n = 0; auto ni = d.multipathMovingCount.find(kv.first); if (ni != d.multipathMovingCount.end()) n = ni->second;
    add(kv.first, kv.second, n);
  }
  return a;
}

json residualSystemArray(const ceqc::model::ResidualStats& r) {
  json a = json::array();
  for (const auto& kv : r.biasRemovedBySystem) {
    json row = statsJson(kv.second);
    row["system"] = kv.first;
    row["system_name"] = systemName(kv.first);
    a.push_back(row);
  }
  return a;
}

json elevationBinsArray(const std::map<std::string, std::vector<int>>& bins) {
  json a = json::array();
  for (const auto& kv : bins) {
    json r = json::object();
    r["name"] = kv.first;
    r["bins"] = intArray(kv.second);
    a.push_back(r);
  }
  return a;
}

json warningsArray(const ceqc::model::QCSummary& s) {
  json a = json::array();
  std::set<std::string> seen;
  auto add = [&](const std::string& w) {
    if (seen.insert(w).second) a.push_back(w);
  };
  for (const auto& w : s.warnings) add(w);
  if (s.derived) {
    for (const auto& w : s.derived->thresholdWarnings) add(std::string("threshold: ") + w);
    for (const auto& w : s.derived->position.warnings) add(std::string("position: ") + w);
  }
  if (s.residuals) for (const auto& w : s.residuals->warnings) add(std::string("residual: ") + w);
  if (s.rtcm3) for (const auto& w : s.rtcm3->warnings) add(std::string("rtcm3: ") + w);
  if (s.ubx) for (const auto& w : s.ubx->warnings) add(std::string("ubx: ") + w);
  return a;
}

} // namespace

void writeQCJson(std::ostream& os, const ceqc::model::QCSummary& s) {
  json root = json::object();
  root["schema"] = "ceqc-qc-json-v1";
  root["generator"] = "ceqc";
  root["format_note"] = "Native CEQC QC JSON. It is not emitted by the +teqc compatibility renderer.";

  json src = json::object();
  src["path"] = s.sourcePath;
  src["rinex_version"] = s.version;
  src["kind"] = ceqc::model::toString(s.kind);
  src["nav_input_files"] = stringArray(s.navInputFiles);
  src["marker_name"] = s.markerName;
  src["marker_number"] = s.markerNumber;
  src["receiver_number"] = s.receiverNumber;
  src["receiver_type"] = s.receiverType;
  src["receiver_version"] = s.receiverVersion;
  src["antenna_number"] = s.antennaNumber;
  src["antenna_type"] = s.antennaType;
  root["source"] = src;

  json tw = json::object();
  tw["first_epoch"] = s.firstEpoch ? ceqc::model::formatUTC(*s.firstEpoch) : std::string{};
  tw["last_epoch"] = s.lastEpoch ? ceqc::model::formatUTC(*s.lastEpoch) : std::string{};
  tw["interval_seconds"] = s.estimatedIntervalS;
  root["time_window"] = tw;

  json counts = json::object();
  counts["epochs"] = s.epochCount;
  counts["observation_records"] = s.observationRecords;
  counts["observation_values_decoded"] = s.observationValues;
  counts["observation_empty_slots"] = s.observationBlankSlots + s.missingObservations;
  counts["observation_decode_failed"] = 0;
  counts["observation_slots_total"] = s.observationValues + s.missingObservations + s.observationBlankSlots;
  counts["navigation_records"] = s.navigationRecords;
  counts["navigation_values"] = s.navigationValues;
  counts["navigation_fields"] = s.navigationFields;
  counts["broadcast_ephemerides"] = s.broadcastEphemerides;
  counts["meteorological_records"] = s.meteorologicalRecords;
  counts["systems"] = mapStringInt(s.systemAppearance);
  counts["system_observation_records"] = mapStringInt(s.systemAppearance);
  counts["system_observation_values"] = mapStringInt(s.systemObservationValues);
  counts["system_blank_slots"] = mapStringInt(s.systemBlankSlots);
  counts["satellites"] = mapStringInt(s.satelliteAppearance);
  counts["navigation_satellites"] = mapStringInt(s.navigationSatelliteAppearance);
  root["counts"] = counts;

  if (s.rtcm3) {
    json r = json::object();
    r["bytes_read"] = static_cast<long long>(s.rtcm3->bytesRead);
    r["frames_total"] = s.rtcm3->frameCount;
    r["frames_good"] = s.rtcm3->goodFrameCount;
    r["bad_crc"] = s.rtcm3->badCRCCount;
    r["bad_length"] = s.rtcm3->badLengthCount;
    r["resync_count"] = s.rtcm3->resyncCount;
    r["decode_failed_messages"] = s.rtcm3->badCRCCount + s.rtcm3->badLengthCount;
    json decodeFailedByReason = json::object();
    decodeFailedByReason["bad_crc"] = s.rtcm3->badCRCCount;
    decodeFailedByReason["bad_length"] = s.rtcm3->badLengthCount;
    decodeFailedByReason["payload_decode_error"] = 0;
    r["decode_failed_by_reason"] = decodeFailedByReason;
    r["epochs"] = s.rtcm3->epochCount;
    r["satellite_observations"] = s.rtcm3->observationCount;
    r["observation_values"] = s.rtcm3->observationValueCount;
    json messageCounts = json::object();
    for (const auto& kv : s.rtcm3->messageCounts) messageCounts[std::to_string(kv.first)] = kv.second;
    r["message_counts"] = messageCounts;
    root["rtcm3"] = r;
  }
  if (s.ubx) {
    json u = json::object();
    u["bytes_read"] = static_cast<long long>(s.ubx->bytesRead);
    u["frames_total"] = s.ubx->frameCount;
    u["frames_good"] = s.ubx->goodFrameCount;
    u["bad_checksum"] = s.ubx->badChecksumCount;
    u["decode_failed_messages"] = s.ubx->badChecksumCount;
    u["rawx_messages"] = s.ubx->rawxCount;
    u["sfrbx_messages"] = s.ubx->sfrbxCount;
    u["epochs"] = s.ubx->epochCount;
    u["satellite_observations"] = s.ubx->observationCount;
    u["observation_values"] = s.ubx->observationValueCount;
    root["ubx"] = u;
  }

  if (s.derived) {
    const auto& d = *s.derived;
    json qc = json::object();
    qc["options_active"] = stringArray(d.optionsActive);
    qc["epoch_sv_min"] = d.epochSVMin;
    qc["epoch_sv_mean"] = d.epochSVMean;
    qc["epoch_sv_max"] = d.epochSVMax;
    qc["lli_count"] = d.lliCount;
    qc["clock_slip_count"] = d.clockSlipCount;
    qc["deleted_observations"] = d.deletedObservations;
    qc["possible_obs_above_horizon"] = d.possibleObsAboveHorizon;
    qc["possible_obs_above_mask"] = d.possibleObsAboveMask;
    qc["complete_obs_above_mask"] = d.completeObsAboveMask;
    qc["deleted_obs_above_mask"] = d.deletedObsAboveMask;
    qc["masked_obs_below_mask"] = d.maskedObsBelowMask;
    qc["unknown_elevation_obs"] = d.unknownElevationObs;
    qc["legacy_gps_glo_satellites"] = d.legacyGPSGLOSatellites;
    root["qc"] = qc;

    json snr = json::object();
    snr["by_observation"] = statMapObject(d.snrStats);
    snr["chart"] = snrChartArray(d.snrStats);
    root["snr"] = snr;

    json ion = json::object();
    ion["by_series"] = statMapObject(d.ionStats);
    ion["histograms"] = elevationBinsArray(d.elevationBins);
    root["ionosphere"] = ion;

    json mp = json::object();
    mp["mp1_rms_m"] = d.mp1Meters;
    mp["mp2_rms_m"] = d.mp2Meters;
    mp["moving_rms_m"] = mapStringDouble(d.multipathMovingRMS);
    mp["moving_samples"] = mapStringInt(d.multipathMovingCount);
    mp["skip_reason"] = d.multipathSkipReason;
    mp["chart"] = multipathChartArray(d);
    root["multipath"] = mp;

    json tp = json::object();
    tp["overall"] = d.timeplot;
    tp["obs"] = d.obsTimeplot;
    tp["nav"] = d.navTimeplot;
    tp["position"] = d.positionTimeplot;
    tp["obs_bin_counts"] = intArray(d.obsBinCounts);
    tp["satellites"] = timeplotRows(d);
    root["timeplots"] = tp;

    json sky = json::object();
    sky["points"] = skyplotArray(d);
    sky["description"] = "Each point is the azimuth/elevation sample at the satellite's maximum elevation in the QC window.";
    root["skyplot"] = sky;

    root["rise_set"] = riseSetArray(d.riseSetEvents);

    json pos = json::object();
    const auto& p = d.position;
    pos["attempted"] = p.attempted;
    pos["skipped_no_navigation"] = p.skippedNoNavigation;
    pos["candidate_epochs"] = p.candidateEpochs;
    pos["numeric_solutions"] = p.epochNumericSolutions;
    pos["accepted_solutions"] = p.epochSolutions;
    pos["skipped_insufficient_svs"] = p.skippedInsufficientSVs;
    pos["rejected_bad_residual"] = p.rejectedBadResidual;
    pos["rejected_bad_height"] = p.rejectedBadHeight;
    pos["rejected_bad_jump"] = p.rejectedBadJump;
    pos["rejected_bad_geometry"] = p.rejectedBadGeometry;
    pos["average_xyz_m"] = xyzArray(p.averageXYZ);
    pos["has_header_approx"] = p.hasApprox;
    pos["header_approx_xyz_m"] = xyzArray(p.approxXYZ);
    pos["used_svs_by_system"] = mapStringInt(p.usedSVsBySystem);
    pos["epoch_positions"] = epochPositions(d.epochPositions);
    root["position"] = pos;
  }

  if (s.residuals) {
    const auto& r = *s.residuals;
    json jr = json::object();
    jr["candidate_observations"] = r.candidateObservations;
    jr["evaluated"] = r.evaluated;
    jr["skipped_no_station"] = r.skippedNoStation;
    jr["skipped_no_ephemeris"] = r.skippedNoEphemeris;
    jr["skipped_no_pseudorange"] = r.skippedNoPseudorange;
    jr["epoch_centered_mean_m"] = r.meanMeters;
    jr["epoch_centered_rms_m"] = r.rmsMeters;
    jr["epoch_centered_max_abs_m"] = r.maxAbsMeters;
    jr["code_minus_range_no_clock_mean_m"] = r.rawMeanMeters;
    jr["code_minus_range_no_clock_rms_m"] = r.rawRmsMeters;
    jr["code_minus_range_no_clock_max_abs_m"] = r.rawMaxAbsMeters;
    jr["by_system"] = residualSystemArray(r);
    jr["raw_by_system"] = statMapArray(r.rawBySystem, "system");
    jr["skipped_no_ephemeris_by_system"] = mapStringInt(r.skippedNoEphemerisBySystem);
    jr["skipped_no_pseudorange_by_system"] = mapStringInt(r.skippedNoPseudorangeBySystem);
    root["residuals"] = jr;
  }

  root["warnings"] = warningsArray(s);

  json charts = json::object();
  if (s.derived) {
    charts["snr_by_observation"] = snrChartArray(s.derived->snrStats);
    charts["multipath_rms"] = multipathChartArray(*s.derived);
    charts["satellite_max_elevation"] = skyplotArray(*s.derived);
    json tb = json::array();
    for (size_t i = 0; i < s.derived->obsBinCounts.size(); ++i) { json r = json::object(); r["bin"] = static_cast<int>(i); r["count"] = s.derived->obsBinCounts[i]; tb.push_back(r); }
    charts["time_bin_observations"] = tb;
  }
  if (s.residuals) charts["residual_rms_by_system"] = residualSystemArray(*s.residuals);
  root["charts"] = charts;

  os << root.dump(2) << "\n";
}

} // namespace ceqc::view
