// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright
// holders. All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief Step3 of the Strangeness tutorial
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram
// STEP 1
// Apply selections on topological variables of V0s
// STEP 2
// Apply PID selections on V0 daughter tracks
// STEP 3
// Check the MC information of the V0s

struct strangeness_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry rKzeroShort{
      "kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambda{
      "kLambda", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rGenParticles{
      "genParticles", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f,
                                 "Accepted z-vertex range (cm)"};

  // Configurable parameters for V0 selection
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1,
                                         "DCA V0 Daughters"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.06,
                                           "DCA Pos To PV"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.06,
                                           "DCA Neg To PV"};
  Configurable<double> v0setting_cospa{
      "v0setting_cospa", 0.98,
      "V0 CosPA"};  // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> v0setting_radius{"v0setting_radius", 0.5, "v0radius"};
  Configurable<float> k0MassLow{"k0MassLow", 0.45f, "Low mass cut for K0"};
  Configurable<float> k0MassHigh{"k0MassHigh", 0.55f, "High mass cut for K0"};
  Configurable<float> lambdaMassLow{"lambdaMassLow", 1.06f,
                                    "Low mass cut for Lambda"};
  Configurable<float> lambdaMassHigh{"lambdaMassHigh", 1.164f,
                                     "High mass cut for Lambda"};

  // Configurable parameters for PID selection
  Configurable<float> NSigmaTPCPion{"NSigmaTPCPion", 4, "NSigmaTPCPion"};
  int mRunNumber = 0;

  void init(InitContext const&) {
    // Axes
    AxisSpec K0ShortMassAxis = {nBins, 0.45f, 0.55f,
                                "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec ptAxis = {nBins, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec LambdaMassAxis = {nBins, 1.06f, 1.164f,
                               "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec",
                        {HistType::kTH1F, {vertexZAxis}});

    // K0s reconstruction
    // Mass
    rKzeroShort.add("hMassK0Short", "hMassK0Short",
                    {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelected", "hMassK0ShortSelected",
                    {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortSelectedTruePions",
                    "hMassK0ShortSelectedTruePions",
                    {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassK0ShortTrueRec", "hMassK0ShortTrueRec",
                    {HistType::kTH1F, {K0ShortMassAxis}});
    // Pt
    rKzeroShort.add("hPtK0ShortSelected", "hPtK0ShortSelected",
                    {HistType::kTH1F, {ptAxis}});
    rKzeroShort.add("hPtK0ShortSelected2d", "hPtK0ShortSelected2d",
                    {HistType::kTH2F, {ptAxis, K0ShortMassAxis}});
    rKzeroShort.add("hPtK0ShortTrueRec", "hPtK0ShortTrueRec",
                    {HistType::kTH1F, {{ptAxis}}});
    rKzeroShort.add("hPtK0ShortTrueRec2d", "hPtK0ShortTrueRec2d",
                    {HistType::kTH2F, {ptAxis, K0ShortMassAxis}});

    // K0s topological/PID cuts
    rKzeroShort.add("hDCAV0Daughters", "hDCAV0Daughters",
                    {HistType::kTH1F, {{55, 0.0f, 2.2f}}});
    rKzeroShort.add("hV0CosPA", "hV0CosPA",
                    {HistType::kTH1F, {{100, 0.95f, 1.f}}});
    rKzeroShort.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s",
                    {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    rKzeroShort.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s",
                    {HistType::kTH2F, {{100, -5.f, 5.f}, {ptAxis}}});
    // Lambda
    rLambda.add("hMassLambda", "hMassLambda",
                {HistType::kTH1F, {LambdaMassAxis}});
    rLambda.add("hMassLambdaSelected", "hMassLambdaSelected",
                {HistType::kTH1F, {LambdaMassAxis}});
    // rLambda.add("hMassLambdaSelectedTruePions",
    // "hMassK0ShortSelectedTruePions", {HistType::kTH1F, {{200, 0.45f,
    // 0.55f}}});
    rLambda.add("hMassLambdaTrueRec", "hMassLambdaTrueRec",
                {HistType::kTH1F, {LambdaMassAxis}});
    // Pt
    rLambda.add("hPtLambdaSelected", "hPtLambdaSelected",
                {HistType::kTH1F, {{ptAxis}}});
    rLambda.add("hPtLambdaTrueRec", "hPtLambdaTrueRec",
                {HistType::kTH1F, {{ptAxis}}});

    // Generated level histograms
    rEventSelection.add("hVertexZGen", "hVertexZGen",
                        {HistType::kTH1F, {vertexZAxis}});
    rGenParticles.add("hPtK0ShortGen", "hPtK0ShortGen",
                      {HistType::kTH1F, {{ptAxis}}});
    rGenParticles.add("hPtK0ShortGen2", "hPtK0ShortGen2",
                      {HistType::kTH2F, {ptAxis, K0ShortMassAxis}});
    rGenParticles.add("hPtLambdaGen", "hPtLambdaGen",
                      {HistType::kTH1F, {{ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutzvertex);

  // Filters on V0s
  // Cannot filter on dynamic columns
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0setting_dcapostopv &&
                        nabs(aod::v0data::dcanegtopv) > v0setting_dcanegtopv &&
                        aod::v0data::dcaV0daughters < v0setting_dcav0dau);
  // Defining the type of the daughter tracks
  using DaughterTracksD =
      soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi>;

  void processData(
      soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&
          collision,
      soa::Filtered<aod::V0Datas> const& V0s, DaughterTracksD const&) {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // V0s
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracksD>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracksD>();

      rKzeroShort.fill(HIST("hMassK0Short"), v0.mK0Short());
      rLambda.fill(HIST("hMassLambda"), v0.mLambda());
      rKzeroShort.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), v0.v0cosPA());

      // Cut on dynamic columns
      if (v0.v0cosPA() < v0setting_cospa) continue;
      if (v0.v0radius() < v0setting_radius) continue;

      if (TMath::Abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (k0MassLow < v0.mK0Short() && v0.mK0Short() < k0MassHigh) {
        // rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"),
        // posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
        // rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"),
        // negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
        rKzeroShort.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
        rKzeroShort.fill(HIST("hPtK0ShortSelected"), v0.pt());
        rKzeroShort.fill(HIST("hPtK0ShortSelected2d"), v0.pt(), v0.mK0Short());
      }
      if (lambdaMassLow < v0.mLambda() && v0.mLambda() < lambdaMassHigh) {
        rLambda.fill(HIST("hMassLambdaSelected"), v0.mLambda());
        rLambda.fill(HIST("hPtLambdaSelected"), v0.pt());
      }
    }
  }
  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra,
                                   aod::pidTPCPi, aod::McTrackLabels>;

  void processRecMC(
      soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&
          collision,
      soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>> const& V0s,
      DaughterTracks const&,  // no need to define a variable for tracks, if we
                              // don't access them directly
      aod::McParticles const&) {
    // Fire up CCDB
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // initCCDB(bc);
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    // V0s
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();

      rKzeroShort.fill(HIST("hMassK0Short"), v0.mK0Short());
      rLambda.fill(HIST("hMassLambda"), v0.mLambda());
      rKzeroShort.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
      rKzeroShort.fill(HIST("hV0CosPA"), v0.v0cosPA());

      // Cut on dynamic columns
      if (v0.v0cosPA() < v0setting_cospa) continue;
      if (v0.v0radius() < v0setting_radius) continue;

      if (TMath::Abs(posDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > NSigmaTPCPion) {
        continue;
      }
      if (lambdaMassLow < v0.mLambda() && v0.mLambda() < lambdaMassHigh) {
        rLambda.fill(HIST("hMassLambdaSelected"), v0.mLambda());
        rLambda.fill(HIST("hPtLambdaSelected"), v0.pt());
      }

      // Filling the PID of the V0 daughters in the region of the K0 peak
      if (k0MassLow < v0.mK0Short() && v0.mK0Short() < k0MassHigh) {
        rKzeroShort.fill(HIST("hNSigmaPosPionFromK0s"),
                         posDaughterTrack.tpcNSigmaPi(),
                         posDaughterTrack.tpcInnerParam());
        rKzeroShort.fill(HIST("hNSigmaNegPionFromK0s"),
                         negDaughterTrack.tpcNSigmaPi(),
                         negDaughterTrack.tpcInnerParam());
        rKzeroShort.fill(HIST("hMassK0ShortSelected"), v0.mK0Short());
        rKzeroShort.fill(HIST("hPtK0ShortSelected"), v0.pt());
        rKzeroShort.fill(HIST("hPtK0ShortSelected2d"), v0.pt(), v0.mK0Short());
      }
      // for lambdas this is more complicated
      if (posDaughterTrack.has_mcParticle() &&
          negDaughterTrack
              .has_mcParticle()) {  // Checking that the daughter tracks come
                                    // from particles and are not fake
        auto posParticle = posDaughterTrack.mcParticle();
        auto negParticle = negDaughterTrack.mcParticle();
        if (posParticle.pdgCode() == 211 &&
            negParticle.pdgCode() ==
                -211) {  // Checking that the daughter tracks are true pions
          rKzeroShort.fill(HIST("hMassK0ShortSelectedTruePions"),
                           v0.mK0Short());
        }
      }

      // Checking that the V0 is a true K0s
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (v0mcParticle.pdgCode() == 310) {
          rKzeroShort.fill(HIST("hMassK0ShortTrueRec"), v0.mK0Short());
          rKzeroShort.fill(
              HIST("hPtK0ShortTrueRec"),
              v0.pt());  // To mimic distribution after the signal extraction
          rKzeroShort.fill(HIST("hPtK0ShortTrueRec2d"), v0.pt(),
                           v0.mK0Short());  // To mimic distribution after the
                                            // signal extraction
        }
        if (TMath::Abs(v0mcParticle.pdgCode()) ==
            3122) {  // Lambda and AntiLambda
          rLambda.fill(HIST("hMassLambdaTrueRec"), v0.mLambda());
          rLambda.fill(
              HIST("hPtLambdaTrueRec"),
              v0.pt());  // To mimic distribution after the signal extraction
          // rLambda.fill(HIST("hPtLambdaTrueRec2"), v0.pt()); // To mimic
          // distribution after the signal extraction
        }
      }
    }
  }

  void processGenMC(
      soa::Filtered<aod::McCollisions>::iterator const& mcCollision,
      const soa::SmallGroups<soa::Join<
          o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>&
          collisions,
      aod::McParticles const& mcParticles) {
    if (collisions.size() < 1)  // to process generated collisions that've been
                                // reconstructed at least once
      return;
    rEventSelection.fill(HIST("hVertexZGen"), mcCollision.posZ());
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() == 310) {
        rGenParticles.fill(HIST("hPtK0ShortGen"), mcParticle.pt());
      }
      if (mcParticle.pdgCode() == 3122) {
        rGenParticles.fill(HIST("hPtLambdaGen"), mcParticle.pt());
      }
    }
  }

  PROCESS_SWITCH(strangeness_tutorial, processRecMC,
                 "Process Run 3 mc, reconstructed", true);
  PROCESS_SWITCH(strangeness_tutorial, processGenMC,
                 "Process Run 3 mc, generated", true);
  PROCESS_SWITCH(strangeness_tutorial, processData, "Process Run 3 data",
                 false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
