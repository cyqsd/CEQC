#pragma once
#include "ceqc/rinex/Rinex.hpp"
#include <map>
#include <optional>
#include <array>
#include <string>
#include <vector>

namespace ceqc::model {
struct QCMetricStats { int count=0; double mean=0,rms=0,min=0,max=0; int lowCount=0; int jumps=0; };
struct QCGapEvent { std::string from; std::string to; double seconds=0; int epochIndex=0; };
struct QCSlipEvent { std::string time; std::string satellite; std::string type; std::string detail; };
struct QCPositionSummary { bool attempted=false; int epochSolutions=0; double averageXYZ[3]{0,0,0}; double approxXYZ[3]{0,0,0}; bool hasApprox=false; std::vector<std::string> warnings; };
struct QCDerivedSummary {
  std::vector<std::string> optionsActive; int epochSVMin=0; int epochSVMax=0; double epochSVMean=0;
  std::vector<QCGapEvent> gapEvents; int lliCount=0; int clockSlipCount=0; std::vector<QCSlipEvent> slipEvents;
  std::map<std::string,QCMetricStats> snrStats, ionStats, iodStats, multipathStats, pseudorangePhase; int codeBandCount=0; int deletedObservations=0;
  int possibleObsAboveHorizon=0, possibleObsAboveMask=0, completeObsAboveMask=0, deletedObsAboveMask=0, maskedObsBelowMask=0, unknownElevationObs=0;
  double mp1Meters=0.0, mp2Meters=0.0; std::map<std::string,double> multipathMovingRMS; std::map<std::string,int> multipathMovingCount; int iodSlipsBelowMask=0, iodSlipsAboveMask=0, iodOrMPSlipsBelowMask=0, iodOrMPSlipsAboveMask=0;
  std::map<std::string,std::vector<int>> elevationBins; std::string timeplot; std::string obsTimeplot; std::map<std::string,std::string> satelliteTimeplot; std::map<std::string,double> satelliteMaxElevationDeg; std::vector<int> obsBinCounts;
  std::vector<int> gpsSVsWithoutObs, gpsSVsWithoutNav, gpsUnhealthySVs, glonassSVsWithoutObs, glonassSVsWithoutNav; int freqTimeCodeCount=0, msecMpEventBins=0;
  int legacyGPSGLOSatellites=0, legacyPossibleObsAboveHorizon=0, legacyPossibleObsAboveMask=0, legacyCompleteObsAboveMask=0, legacyDeletedObsAboveMask=0, legacyMaskedObsBelowMask=0;
  std::vector<std::string> symbolLegend; QCPositionSummary position;
};
struct ResidualStats {
  int candidateObservations=0,evaluated=0,skippedNoStation=0,skippedNoEphemeris=0,skippedNoPseudorange=0;
  int gps=0, glonass=0, galileo=0, beidou=0; double meanMeters=0,rmsMeters=0,maxAbsMeters=0; std::map<std::string,int> bySystem; std::vector<std::string> warnings;
};
struct QCOptions {
  std::map<std::string,bool> noOrbitSystems, noPositionSystems; std::string orbitSpec;
  bool averagePosition=true, everyEpochPosition=false, everyEpochXYZ=false, everyEpochGeodetic=false, everyEpochDecimal=false, clockSlips=true;
  double codeSigmas=2.0; bool dataIndicators=true; double epsilon=1.387779e-17; bool horizon=true, unhealthy=false, ion=true, iod=true;
  int bins=18, ionBins=18; double ionJumpCM=1e35, iodJumpCMPerMin=400, gapMinMinutes=10, gapMaxNoNavMinutes=90; bool longReport=true, shortReport=true, report=true, lli=true, movingAverage=true, mask=false;
  std::map<std::string,int> minSNR{{"1",4},{"2",4},{"5",4},{"6",4},{"7",4},{"8",4}}; int minSVs=6; bool multipath=true, mpRaw=false, pseudorangePhase=false, plot=false, positionOnly=false, riseSet=true, snr=true, ssv=false, svpr=false, symbolCodes=false, allSymbols=false, yCode=false;
  std::map<std::string,double> mpRMSCM{{"12",65},{"21",65},{"15",65},{"51",65},{"16",80},{"61",80},{"17",90},{"71",90},{"18",100},{"81",100}};
  double mpSigmas=4; int mpWindow=50; double mpCAASPercent=100, msecTol=1e-2; std::vector<std::string> navFiles; double positionConvM=0.1, positionHMinM=-500, positionHMaxM=9000, positionJumpM=100000;
  std::string root; double setHorizonDeg=0, setMaskDeg=10, setComparisonDeg=25; bool slipsEnabled=false, slipsAppend=false; std::string slipsTarget; int snBins=18, mpBins=18, width=72;
};
struct QCSummary {
  std::string sourcePath; std::vector<std::string> navInputFiles; std::map<std::string,int> navigationSatelliteAppearance; std::string markerName; std::string markerNumber; std::string receiverNumber; std::string receiverType; std::string receiverVersion; std::string antennaNumber; std::string antennaType; double version=0; RinexKind kind=RinexKind::Unknown; int epochCount=0, observationRecords=0, observationValues=0, missingObservations=0;
  int navigationRecords=0, navigationValues=0, navigationFields=0, broadcastEphemerides=0, meteorologicalRecords=0, meteorologicalValues=0; std::optional<TimePoint> firstEpoch,lastEpoch; double estimatedIntervalS=0;
  std::map<std::string,int> satelliteAppearance, systemAppearance; std::vector<ValidationIssue> headerIssues; std::vector<std::string> warnings; std::optional<ResidualStats> residuals; std::optional<QCDerivedSummary> derived; std::optional<RTCM3Summary> rtcm3; std::optional<UBXSummary> ubx;
};
}
