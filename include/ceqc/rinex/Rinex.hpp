#pragma once
#include "ceqc/core/Time.hpp"
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace ceqc::model {

enum class RinexKind { Unknown, Obs, Nav, Met };

inline std::string toString(RinexKind k) {
  switch (k) {
  case RinexKind::Obs: return "obs";
  case RinexKind::Nav: return "nav";
  case RinexKind::Met: return "met";
  default: return "unknown";
  }
}

struct HeaderLine {
  std::string raw;
  std::string value;
  std::string label;
};

struct RinexHeader {
  double version = 0.0;
  RinexKind kind = RinexKind::Unknown;
  std::string satelliteSystem;
  std::vector<HeaderLine> lines;
  std::map<std::string, std::vector<size_t>> byLabel;
  std::map<std::string, std::vector<std::string>> obsTypes;
  std::vector<std::string> metTypes;
};

struct ObservationEpoch {
  TimePoint time{};
  int flag = 0;
  std::vector<std::string> satellites;
  size_t lineIndex = 0;
  std::string rawLine;
};

struct ObservationValue {
  std::string type;
  std::optional<double> value;
  std::string lli;
  std::string ssi;
  std::string raw;
};

struct ObservationRecord {
  size_t epochIndex = 0;
  TimePoint time{};
  int flag = 0;
  std::string satellite;
  std::string system;
  std::vector<ObservationValue> values;
  std::vector<std::string> rawLines;
};

struct NavigationField {
  std::string name;
  std::string unit;
  double value = 0.0;
  size_t index = 0;
};

struct NavigationRecord {
  std::string satellite;
  std::string system;
  std::optional<TimePoint> epoch;
  std::string recordType = "EPH";
  std::string messageType;
  std::string messageSubtype;
  std::vector<double> values;
  std::map<std::string, NavigationField> fields;
  std::vector<std::string> rawLines;
  size_t lineIndex = 0;
};

struct MeteorologicalRecord {
  TimePoint time{};
  std::map<std::string, double> values;
  std::string rawLine;
  size_t lineIndex = 0;
};

struct RinexData {
  std::vector<ObservationEpoch> observationEpochs;
  std::vector<ObservationRecord> observationRecords;
  std::vector<NavigationRecord> navigationRecords;
  std::vector<MeteorologicalRecord> meteorologicalRecords;
  std::vector<std::string> parseWarnings;
};

struct RTCM3Station { int id=0; int itrf=0; double x=0,y=0,z=0,antennaHeight=0; bool hasPosition=false; };
struct RTCM3Ephemeris { int messageType=0; std::string system; std::string satellite; int count=0; };
struct RTCM3Summary {
  std::string sourcePath; long long bytesRead=0; int frameCount=0; int goodFrameCount=0; int badCRCCount=0; int badLengthCount=0; int resyncCount=0;
  std::map<int,int> messageCounts; std::map<int,int> stationCounts; std::map<int,RTCM3Station> stations;
  int epochCount=0; int observationCount=0; int observationValueCount=0; std::map<std::string,RTCM3Ephemeris> ephemeris; std::vector<std::string> warnings;
};
struct UBXSummary {
  std::string sourcePath; long long bytesRead=0; int frameCount=0; int goodFrameCount=0; int badChecksumCount=0; int nmeaLines=0; int rawxCount=0; int sfrbxCount=0;
  int epochCount=0; int observationCount=0; int observationValueCount=0; std::map<std::string,int> messageCounts; std::vector<std::string> warnings;
};

struct RinexFile {
  std::string path;
  RinexHeader header;
  std::vector<std::string> body;
  RinexData data;
  std::optional<RTCM3Summary> rtcm3;
  std::optional<UBXSummary> ubx;
};

struct ValidationIssue {
  std::string severity;
  std::string message;
  std::string label;
};
}
