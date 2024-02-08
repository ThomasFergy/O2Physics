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
/// \brief this is a starting point for the Strangeness tutorial
/// \author Nepeivoda Roman (roman.nepeivoda@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)

#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 0
// Starting point: loop over all V0s and fill invariant mass histogram

struct strangeness_tutorial {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection",
                                    {},
                                    OutputObjHandlingPolicy::AnalysisObject,
                                    true,
                                    true};
  HistogramRegistry rKzeroShort{
      "kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f,
                                 "Accepted z-vertex range (cm)"};
  Configurable<float> k0MassLow{"k0MassLow", 0.45f, "Low mass cut for K0"};
  Configurable<float> k0MassHigh{"k0MassHigh", 0.55f, "High mass cut for K0"};
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
    rKzeroShort.add("hMassK0Short", "hMassK0Short",
                    {HistType::kTH1F, {K0ShortMassAxis}});
    rKzeroShort.add("hMassLambda", "hMassLambda",
                    {HistType::kTH1F, {LambdaMassAxis}});
    rKzeroShort.add("hMassAntiLambda", "hMassAntiLambda",
                    {HistType::kTH1F, {LambdaMassAxis}});
    // rKzeroShort.add("hPtK0ShortSelected", "hPtK0ShortSelected",
    //                 {HistType::kTH1F, {{ptAxis}}});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  void process(
      soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&
          collision,
      aod::V0Datas const& V0s) {
    // Fill the event counter
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    for (const auto& v0 : V0s) {
      rKzeroShort.fill(HIST("hMassK0Short"), v0.mK0Short());
      rKzeroShort.fill(HIST("hMassLambda"), v0.mLambda());
      rKzeroShort.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());
      // rKzeroShort.fill(HIST("hPtK0ShortSelected"), v0.pt());
      // rKzeroShort.fill(HIST("hPtLambdaSelected"), v0.pt());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{adaptAnalysisTask<strangeness_tutorial>(cfgc)};
}
