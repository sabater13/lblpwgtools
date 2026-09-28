// numu CC sample
//
// uses numuCCcuts.h

#include "TVector3.h"
#include "TFile.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

#include "Cuts/numuCCcuts.h"

using namespace ana;


// Input CAFs
//const std::string kInputCAF = "/pnfs/dune/persistent/physicsgroups/dunendsim/abooth/nd-production/MicroProdN4p1/run-cafmaker/MicroProdN4p1_NDComplex_FHC.caf.full.spineonly/CAF.flat/0002000/*/*.root";
//const std::string kInputCAF = "/pnfs/dune/persistent/physicsgroups/dunendproto/dunendsim/MiniProdN5-Small-Files/run-cafmaker/MiniProdN5p3_NDComplex_FHC.caf.full.sanddrift.spineonly/CAF.flat/*/*.root";
const std::string kInputCAF = "/pnfs/dune/persistent/physicsgroups/dunendproto/dunendsim/MiniProdN5/run-cafmaker/MiniProdN5p1_NDComplex_FHC.caf.full.spineonly.sanddrift/CAF.flat/*/*.root";

// Output
const std::string kOutputDir  = "Plotting";
const std::string kOutputName = "numuCCsample.root";
const std::string kOutputPath = kOutputDir + "/" + kOutputName;

// Selection parameters
const double kFVMargin = 25.;  // [cm] default FV margin on every wall (scan varies upstream Z only)
const double kMuonMinLen = 25.;  // [cm] min primary-muon track length in the base selection

// FV scans
const std::vector<double> kZUpstreamMargins = {0, 5, 10, 15, 20, 25, 30, 40};  // [cm]
const std::vector<double> kMuonLenThresholds = {0, 5, 10, 15, 20, 30, 50};  // [cm]


// =============================================================================
// HELPERS  (macro-specific, not in cuts header)
// =============================================================================
// Best-overlap truth interaction for a reco interaction (nullptr if none).
static const caf::SRTrueInteractionProxy* BestTruthMatch(const caf::SRInteractionProxy *ixn){
    if (ixn->truth.empty()) return nullptr;
    const caf::SRProxy *sr = ixn->Ancestor<caf::SRProxy>();
    size_t tidx = 0; float maxOverlap = 0;
    for (size_t i = 0; i < ixn->truthOverlap.size(); ++i)
        if (ixn->truthOverlap[i] > maxOverlap) { tidx = i; maxOverlap = ixn->truthOverlap[i]; }
    return caf::FindInteraction(sr->mc, ixn->truth[tidx]);
}

// Turn any "f(true interaction) -> double" into a Var evaluated on the
// truth-matched interaction. Returns `fallback` if there is no truth match.
template <typename F> static Var MatchedVar(F f, double fallback = -1.){
    return Var([f, fallback](const caf::SRInteractionProxy *ixn) -> double {
        const caf::SRTrueInteractionProxy *t = BestTruthMatch(ixn);
        return t ? f(t) : fallback;
    });
}

static bool TrueVtxInFV(const caf::SRTrueInteractionProxy *t, double zUpMargin = kFVMargin){
    return kTrueVtxX(t) > NDLArXLo + kFVMargin && kTrueVtxX(t) < NDLArXHi - kFVMargin
        && kTrueVtxY(t) > NDLArYLo + kFVMargin && kTrueVtxY(t) < NDLArYHi - kFVMargin
        && kTrueVtxZ(t) > NDLArZLo + zUpMargin  && kTrueVtxZ(t) < NDLArZHi - kFVMargin;
}

// angle with respect to the beam direction
static double TrueMuonAngleDeg(const caf::SRTrueInteractionProxy *t){
    for (const auto &p : t->prim)
        if (std::abs(p.pdg) == 13) {
            TVector3 dir(p.p.px, p.p.py, p.p.pz);
            return dir.Angle(beam_dir_truth) * 180.0 / TMath::Pi();
        }
    return -1.;
}

// y = 1 - E_mu/E_nu for CC (exact in the lab frame). bjorkenY does not exist in
// so build it from the primary muon.
static double TrueInelasticity(const caf::SRTrueInteractionProxy *t){
    const double enu = t->E;
    if (enu <= 0) return -1.;
    for (const auto &p : t->prim)
        if (std::abs(p.pdg) == 13) {
            const double pmag = std::sqrt(p.p.px*p.p.px + p.p.py*p.p.py + p.p.pz*p.p.pz);
            const double emu  = std::sqrt(pmag*pmag + 0.10566*0.10566); // muon mass [GeV]
            return 1.0 - emu / enu;
        }
    return -1.;
}

static double ModeCategory(const caf::SRTrueInteractionProxy *t){
    switch ((int)t->mode) {
        case 0:  return 0.; // QE
        case 10: return 1.; // MEC / 2p2h
        case 1:  return 2.; // RES
        case 2:  return 3.; // DIS
        case 3:  return 4.; // COH
        default: return 5.; // other
    }
}


// =============================================================================
// VARIABLES
// =============================================================================
// Reco energies
// Calorimetric neutrino energy: sum of all primary particles
const Var kRecoNuE_Calo([](const caf::SRInteractionProxy *ixn) -> double {
    double E = 0;
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p)) E += kPartE(&p);
    return E;
});

// Hadronic energy: all primaries except the muon
// Might need to also add another variable for visible hadronic energy
// i.e. no neutrons
const Var kRecoHadronicE([](const caf::SRInteractionProxy *ixn) -> double {
    double E = 0;
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p) && !kIsMuon(&p)) E += kPartE(&p);
    return E;
});

// Reco inelasticity y = Ehad / Etotal
const Var kRecoInelasticity([](const caf::SRInteractionProxy *ixn) -> double {
    double totalE = 0, hadE = 0;
    for (const auto &p : ixn->part.dlp) {
        if (!kIsPrimary(&p)) continue;
        const double e = kPartE(&p);
        totalE += e;
        if (!kIsMuon(&p)) hadE += e;
    }
    return (totalE > 0) ? hadE / totalE : -1.0;
});

// Truth variables
const TruthVar kTrueInelasticityVar([](const caf::SRTrueInteractionProxy *t) -> double { return TrueInelasticity(t); });
const TruthVar kTrueMuonAngle([](const caf::SRTrueInteractionProxy *t) -> double { return TrueMuonAngleDeg(t); });

// Truth-matched
const Var kTruthMatchedInelasticity = MatchedVar([](const caf::SRTrueInteractionProxy *t) { return TrueInelasticity(t); });
const Var kTruthMatchedTrueMuonAngle = MatchedVar([](const caf::SRTrueInteractionProxy *t) { return TrueMuonAngleDeg(t); });
const Var kTruthMatchedModeCategory = MatchedVar([](const caf::SRTrueInteractionProxy *t) { return ModeCategory(t); });

// (reco - true) / true
const Var kInelResidual([](const caf::SRInteractionProxy *ixn) -> double {
    const double trueY = kTruthMatchedInelasticity(ixn);
    if (trueY <= 0.01) return -999.; // guard against y~0 and missing truth
    return (kRecoInelasticity(ixn) - trueY) / trueY;
});

const Var kNuEResidual([](const caf::SRInteractionProxy *ixn) -> double {
    const double trueE = kTruthMatchedNuE(ixn);
    if (trueE <= 0) return -999.;
    return (kRecoNuE_Calo(ixn) - trueE) / trueE;
});

// Cut-flow dummy variables
const Var kDummyVar([](const caf::SRInteractionProxy*) -> double { return 0.5; });
const TruthVar kTruthDummyVar([](const caf::SRTrueInteractionProxy*) -> double { return 0.5; });

// Truth-match category for cut-flow breakdown
// 0 = no truth match
// 1 = truth-matched signal (numu CC, in FV, contained)
// 2 = truth-matched, but background
const Var kTruthMatchCategory([](const caf::SRInteractionProxy *ixn) -> double {
    if (ixn->truth.empty())        return 0.;
    if (kNuMuTrueSignalMatch(ixn)) return 1.;
    return 2.;
});

// Detailed truth-match category
// 0 = no truth match
// 1 = truth-matched signal
// 2 = numu CC, but true vertex out of FV ("OOFV")
// 3 = not numu CC at all (NC, nue CC, nutau CC)
// 4 = numu CC, in FV, failed some other signal criterion (e.g. containment)
const Var kTruthMatchCategoryDetailed([](const caf::SRInteractionProxy *ixn) -> double {
    if (ixn->truth.empty())        return 0.;
    if (kNuMuTrueSignalMatch(ixn)) return 1.;

    const caf::SRTrueInteractionProxy *t = BestTruthMatch(ixn);
    if (!t || !kTrueNumuCC(t))     return 3.;
    if (!TrueVtxInFV(t))           return 2.;
    return 4.;
});


// One cut-flow stage, filled three ways (plain / 3-category / detailed category)
struct CutFlowStage {
    std::string name;
    std::unique_ptr<Spectrum> plain, cat, detailed;
};
struct FVScanPoint {
    double margin;
    std::unique_ptr<Spectrum> num, den, all, oofvDen, oofvLeak;
};
struct LenScanPoint {
    double threshold;
    std::unique_ptr<Spectrum> num, all, trueMuonCand;
};


// =============================================================================
// MAIN
// =============================================================================
void numuCCsample(){

    SpectrumLoader loader(kInputCAF);

    // -------------------------------------------------------------------------
    // Binning
    // -------------------------------------------------------------------------
    const Binning ENuBins           = Binning::Simple(25, 0, 10);     // also used for the (square) migration matrix
    const Binning EvisBins          = Binning::Simple(25, 0, 10);
    const Binning MuonPBins         = Binning::Simple(20, 0, 6);
    const Binning MuonAngleBins     = Binning::Simple(18, 0, 180);
    const Binning ResidualBins      = Binning::Simple(100, -1, 1);
    const Binning InelBins          = Binning::Simple(20, 0, 1);
    const Binning HadEBins          = Binning::Simple(20, 0, 5);
    const Binning ModeBins          = Binning::Simple(7, -1.5, 5.5);  // -1 = unknown, 0..5 see ModeCategory
    const Binning BkgCatBins        = Binning::Simple(7, -0.5, 6.5);  // 0..6, see kBkgCategory
    const Binning OneBin            = Binning::Simple(1, 0, 1);
    const Binning MuonCandLenBins   = Binning::Simple(30, 0, 100);    // fine binning at short lengths where PID breaks down
    const Binning IsTrueMuonBins    = Binning::Simple(3, -1.5, 1.5);  // -1, 0, 1
    const Binning TruePDGCatBins    = Binning::Simple(8, -0.5, 7.5);
    const Binning MuonVtxDistBins   = Binning::Simple(500, 0, 500);
    const Binning TruthMatchCatBins         = Binning::Simple(3, -0.5, 2.5);
    const Binning TruthMatchCatDetailedBins = Binning::Simple(5, -0.5, 4.5);
    // Vertex-position efficiency binning: inside the FV only
    const Binning VtxXEffBins = Binning::Simple(30, NDLArXLo + kFVMargin, NDLArXHi - kFVMargin);
    const Binning VtxYEffBins = Binning::Simple(30, NDLArYLo + kFVMargin, NDLArYHi - kFVMargin);
    const Binning VtxZEffBins = Binning::Simple(30, NDLArZLo + kFVMargin, NDLArZHi - kFVMargin);


    // -------------------------------------------------------------------------
    // Axes
    // -------------------------------------------------------------------------
    // Neutrino energy
    HistAxis axRecoENu("Reco E_{#nu} [GeV]", ENuBins, kRecoNuE_Calo);
    HistAxis axMatchedNuE("Truth-matched true E_{#nu} [GeV]", ENuBins, kTruthMatchedNuE);
    TruthHistAxis axTrueENu("True E_{#nu} [GeV]", ENuBins, kTrueNuE);

    HistAxis axNuEResVsTrueNuE("True E_{#nu} [GeV]", ENuBins, kTruthMatchedNuE,
                               "(Reco - True)/True E_{#nu}", ResidualBins, kNuEResidual);
    HistAxis axNuERes1D("(Reco - True)/True E_{#nu}", ResidualBins, kNuEResidual);

    HistAxis axMigration("True E_{#nu} [GeV]", ENuBins, kTruthMatchedNuE,
                         "Reco E_{#nu} [GeV]", ENuBins, kRecoNuE_Calo);
    HistAxis axEvisMigration("True E_{vis} [GeV]", EvisBins, kTruthMatchedEvis,
                             "Reco E_{vis} [GeV]", EvisBins, kRecoNuE_Calo);

    // Muon kinematics
    TruthHistAxis axTrueMuonP("True Muon Momentum [GeV/c]", MuonPBins, kTrueMuonP);
    HistAxis axMatchedMuonP("True Muon Momentum [GeV/c]", MuonPBins, kTruthMatchedTrueMuonP);
    TruthHistAxis axTrueMuAngle("True Muon Angle [deg]", MuonAngleBins, kTrueMuonAngle);
    HistAxis axMatchedMuAngle("True Muon Angle [deg]", MuonAngleBins, kTruthMatchedTrueMuonAngle);

    // Vertex position (truth-matched reco vs. all true)
    HistAxis axMatchedVtxX("Truth-matched true vtx X [cm]", VtxXEffBins, kMatchedTrueVtxX);
    HistAxis axMatchedVtxY("Truth-matched true vtx Y [cm]", VtxYEffBins, kMatchedTrueVtxY);
    HistAxis axMatchedVtxZ("Truth-matched true vtx Z [cm]", VtxZEffBins, kMatchedTrueVtxZ);
    TruthHistAxis axTrueVtxX("True vtx X [cm]", VtxXEffBins, kTrueVtxX);
    TruthHistAxis axTrueVtxY("True vtx Y [cm]", VtxYEffBins, kTrueVtxY);
    TruthHistAxis axTrueVtxZ("True vtx Z [cm]", VtxZEffBins, kTrueVtxZ);

    // Hadronic energy and inelasticity
    HistAxis axHadE("Reco Hadronic Energy [GeV]", HadEBins, kRecoHadronicE);
    HistAxis axRecoInel("Reco Inelasticity y", InelBins, kRecoInelasticity);
    TruthHistAxis axTrueInel("True Inelasticity y", InelBins, kTrueInelasticityVar);
    HistAxis axInelResVsTrueInel("True Inelasticity y", InelBins, kTruthMatchedInelasticity,
                                 "(Reco - True)/True y", ResidualBins, kInelResidual);

    // Muon-candidate diagnostics
    HistAxis axMuonCandLenVsIsTrueMuon("Muon candidate track length [cm]", MuonCandLenBins, kMuonCandidateLen,
                                       "Is true muon", IsTrueMuonBins,  kMuonCandidateIsTrueMuon);
    HistAxis axMuonCandLenVsTruePDG("Muon candidate track length [cm]", MuonCandLenBins, kMuonCandidateLen,
                                    "True PDG category", TruePDGCatBins,  kTruePDGofRecoPrimaryMuon);
    HistAxis axMuonToVtxDist("Primary muon start - vertex distance [cm]", MuonVtxDistBins, kMuonStartToVtxDist);
    HistAxis axMuonVtxDistVsIsTrueMuon("Muon start - vertex distance [cm]", MuonVtxDistBins, kMuonStartToVtxDist,
                                       "Is true muon", IsTrueMuonBins,  kMuonCandidateIsTrueMuon);
    HistAxis axMuonVtxDistVsSignal("Primary muon start - vertex distance [cm]", MuonVtxDistBins, kMuonStartToVtxDist,
                                   "Truth-match category", TruthMatchCatBins, kTruthMatchCategory);

    // Interaction mode / background category
    HistAxis axMode("Interaction mode", ModeBins, kTruthMatchedModeCategory);
    HistAxis axModeVsENu ("Truth-matched E_{#nu} [GeV]", ENuBins, kTruthMatchedNuE,
                          "Interaction mode", ModeBins,  kTruthMatchedModeCategory);
    HistAxis axModeVsMuonP("True Muon Momentum [GeV/c]", MuonPBins, kTruthMatchedTrueMuonP,
                           "Interaction mode", ModeBins, kTruthMatchedModeCategory);
    HistAxis axBkgCatVsENu("Reco E_{#nu} [GeV]", ENuBins, kRecoNuE_Calo,
                           "Background category", BkgCatBins, kBkgCategory);
    HistAxis axBkgCatVsMuonP("True Muon Momentum [GeV/c]", MuonPBins, kTruthMatchedTrueMuonP,
                             "Background category", BkgCatBins, kBkgCategory);

    // Cut-flow axes
    HistAxis axDummy("", OneBin, kDummyVar);
    TruthHistAxis axTruthDummy("", OneBin, kTruthDummyVar);
    HistAxis axTruthMatchCat("Truth-match category", TruthMatchCatBins, kTruthMatchCategory);
    HistAxis axTruthMatchCatDetailed("Truth-match category (detailed)", TruthMatchCatDetailedBins, kTruthMatchCategoryDetailed);


    // -------------------------------------------------------------------------
    // Selections
    // -------------------------------------------------------------------------
    const Cut kBaseSelection = kPartLenInInteractionCut_LongestTrack
                            && kVtxInFV
                            && kEventContained
                            && kRecoNumuCC
                            && kPrimMuonTrackLengthCut(kMuonMinLen)
                            && kIsBestMatchForTruth;
    const Cut kSignalSelection = kBaseSelection && kNuMuTrueSignalMatch;

    // Signal split by what the muon candidate is / where it ends up
    const Cut kSigMuonContained   = kSignalSelection && kMuonCandContained;
    const Cut kSigMuonEscapesTMS  = kSignalSelection && kMuonCandEscapesToTMS;
    const Cut kSigMuonNotTrueMuon = kSignalSelection && kMuonCandIsNotTrueMuon;

    // Shorthand for the reco-interaction loader and the true-signal denominator
    auto reco = [&]() -> decltype(auto) { return loader.Interactions(RecoType::kDLP); };
    const TruthCut& kTrueSignal = kTrueNumuCCContainedNDLArAndEscapeToTMS;

    // =========================================================================
    // CUT FLOW  (cumulative; each stage filled plain / 3-cat / detailed)
    // =========================================================================
    std::vector<CutFlowStage> cutFlow;

    auto addStage = [&](const std::string &name, const Cut *cut) {
        CutFlowStage s;
        s.name = name;
        if (cut) {
            s.plain    = std::make_unique<Spectrum>(reco()[*cut], axDummy);
            s.cat      = std::make_unique<Spectrum>(reco()[*cut], axTruthMatchCat);
            s.detailed = std::make_unique<Spectrum>(reco()[*cut], axTruthMatchCatDetailed);
        } else { // no cut: all reco interactions
            s.plain    = std::make_unique<Spectrum>(reco(), axDummy);
            s.cat      = std::make_unique<Spectrum>(reco(), axTruthMatchCat);
            s.detailed = std::make_unique<Spectrum>(reco(), axTruthMatchCatDetailed);
        }
        cutFlow.push_back(std::move(s));
    };

    addStage("AllReco", nullptr);
    {
        Cut running = kPartLenInInteractionCut_LongestTrack;             // longest track > 1 mm
        addStage("HasParticles", &running);
        running = running && kVtxInFV;                                   // vertex in FV
        addStage("InFV", &running);
        running = running && kEventContained;                            // non-muon particles contained
        addStage("Contained", &running);
        running = running && kRecoNumuCC;                                // has primary muon candidate
        addStage("HasMuon", &running);
        running = running && kPrimMuonTrackLengthCut(kMuonMinLen);       // muon candidate long enough
        addStage("HasMuonLongerThan25cm", &running);
        addStage("BestMatch", &kBaseSelection);                          // best truth match (= full base selection)
    }

    // Reference counts for headline efficiency/purity
    Spectrum sCF_Signal    (reco()[kSignalSelection],           axDummy);
    Spectrum sCF_TrueSignal(loader.NuTruths()[kTrueSignal],     axTruthDummy);

    // =========================================================================
    // NEUTRINO ENERGY ESTIMATOR
    // =========================================================================
    Spectrum sRecoNuE       (reco()[kSignalSelection], axRecoENu);
    Spectrum sNuERes2D      (reco()[kSignalSelection], axNuEResVsTrueNuE);
    Spectrum sNuERes1D      (reco()[kSignalSelection], axNuERes1D);
    Spectrum sMigrationMatrix(reco()[kSignalSelection], axMigration);
    Spectrum sEvisMigration (reco()[kSignalSelection], axEvisMigration);

    // Split: muon contained / escapes to TMS / candidate is not a true muon
    Spectrum sNuERes2D_Contained    (reco()[kSigMuonContained],   axNuEResVsTrueNuE);
    Spectrum sNuERes2D_EscapesTMS   (reco()[kSigMuonEscapesTMS],  axNuEResVsTrueNuE);
    Spectrum sNuERes2D_NotTrueMuon  (reco()[kSigMuonNotTrueMuon], axNuEResVsTrueNuE);
    Spectrum sMigrationMatrix_Contained  (reco()[kSigMuonContained],   axMigration);
    Spectrum sMigrationMatrix_EscapesTMS (reco()[kSigMuonEscapesTMS],  axMigration);
    Spectrum sMigrationMatrix_NotTrueMuon(reco()[kSigMuonNotTrueMuon], axMigration);

    // Efficiency / purity vs. Enu
    Spectrum sTrueNumuCC_ENu       (loader.NuTruths()[kTrueSignal], axTrueENu);
    Spectrum sRecoNumuCC_All_ENu   (reco()[kBaseSelection],         axMatchedNuE);
    Spectrum sRecoNumuCC_Matched_ENu(reco()[kSignalSelection],      axMatchedNuE);

    // =========================================================================
    // HADRONIC ENERGY AND INELASTICITY
    // =========================================================================
    Spectrum sRecoHadE       (reco()[kSignalSelection],       axHadE);
    Spectrum sRecoInelasticity(reco()[kSignalSelection],      axRecoInel);
    Spectrum sTrueInelasticity(loader.NuTruths()[kTrueSignal], axTrueInel);
    Spectrum sInelRes2D      (reco()[kSignalSelection],       axInelResVsTrueInel);

    // =========================================================================
    // MUON-CANDIDATE DIAGNOSTICS
    // =========================================================================
    Spectrum sMuonCandLenVsIsTrueMuon     (reco()[kBaseSelection],   axMuonCandLenVsIsTrueMuon);
    Spectrum sMuonCandLenVsTruePDG        (reco()[kBaseSelection],   axMuonCandLenVsTruePDG);
    Spectrum sMuonToVtxDist               (reco()[kSignalSelection], axMuonToVtxDist);
    Spectrum sMuonVtxDistVsIsTrueMuon     (reco()[kBaseSelection],   axMuonVtxDistVsIsTrueMuon);
    Spectrum sMuonVtxDistVsIsTrueMuon_signal(reco()[kSignalSelection], axMuonVtxDistVsIsTrueMuon);
    Spectrum sMuonVtxDistVsSignal         (reco()[kBaseSelection],   axMuonVtxDistVsSignal);

    // =========================================================================
    // INTERACTION MODE / BACKGROUND BREAKDOWN  (selected sample, no truth-match requirement)
    // =========================================================================
    Spectrum sModeVsENu_selected  (reco()[kBaseSelection], axModeVsENu);
    Spectrum sModeVsMuonP_selected(reco()[kBaseSelection], axModeVsMuonP);
    Spectrum sBkgCatVsENu         (reco()[kBaseSelection], axBkgCatVsENu);
    Spectrum sBkgCatVsMuonP       (reco()[kBaseSelection], axBkgCatVsMuonP);
    Spectrum sModeSignalOnly      (reco()[kSignalSelection],                        axMode);
    Spectrum sModeBkgOnly         (reco()[kBaseSelection && !kNuMuTrueSignalMatch], axMode);

    // =========================================================================
    // EFFICIENCY VS. MUON ANGLE
    // =========================================================================
    Spectrum sRecoNumuCC_MuonAngle    (reco()[kSignalSelection],       axMatchedMuAngle);
    Spectrum sTrueNumuCC_MuonAngle    (loader.NuTruths()[kTrueSignal], axTrueMuAngle);
    Spectrum sRecoNumuCC_All_MuonAngle(reco()[kBaseSelection],         axMatchedMuAngle);

    // =========================================================================
    // EFFICIENCY VS. TRUE VERTEX POSITION (spatial uniformity)
    // Flat efficiency = no spatial bias; drops near the FV edge tell you
    // whether the FV margin is appropriate.
    // =========================================================================
    Spectrum sRecoNumuCC_VtxX(reco()[kSignalSelection], axMatchedVtxX);
    Spectrum sTrueNumuCC_VtxX(loader.NuTruths()[kTrueSignal], axTrueVtxX);
    Spectrum sRecoNumuCC_VtxY(reco()[kSignalSelection], axMatchedVtxY);
    Spectrum sTrueNumuCC_VtxY(loader.NuTruths()[kTrueSignal], axTrueVtxY);
    Spectrum sRecoNumuCC_VtxZ(reco()[kSignalSelection], axMatchedVtxZ);
    Spectrum sTrueNumuCC_VtxZ(loader.NuTruths()[kTrueSignal], axTrueVtxZ);

    // =========================================================================
    // FV SCAN -- vary the upstream Z margin (other walls fixed at kFVMargin)
    //   num     : signal passing all cuts with this FV
    //   den     : true signal in this FV
    //   all     : all selected reco events (purity denominator)
    //   oofvDen : true numu CC outside this FV
    //   oofvLeak: selected reco events truth-matched to an OOFV numu CC
    // =========================================================================
    std::vector<FVScanPoint> fvScan;
    for (double margin : kZUpstreamMargins) {
        const Cut recoFV = kVtxX > NDLArXLo + kFVMargin && kVtxX < NDLArXHi - kFVMargin
                        && kVtxY > NDLArYLo + kFVMargin && kVtxY < NDLArYHi - kFVMargin
                        && kVtxZ > NDLArZLo + margin    && kVtxZ < NDLArZHi - kFVMargin;

        const TruthCut trueFV = kTrueVtxX > NDLArXLo + kFVMargin && kTrueVtxX < NDLArXHi - kFVMargin
                             && kTrueVtxY > NDLArYLo + kFVMargin && kTrueVtxY < NDLArYHi - kFVMargin
                             && kTrueVtxZ > NDLArZLo + margin    && kTrueVtxZ < NDLArZHi - kFVMargin;

        const TruthCut trueSignal = kTrueNumuCC && trueFV && kAllTrueContainedExceptMuonDownstream;
        const TruthCut trueOOFV   = kTrueNumuCC && !trueFV;

        const Cut truthMatch = MakeNuMuTrueSignalMatch(trueFV);

        // NOTE: unlike kBaseSelection this omits kPrimMuonTrackLengthCut -- see notes.
        const Cut fullSel = kPartLenInInteractionCut_LongestTrack
                         && recoFV
                         && kEventContained
                         && kRecoNumuCC
                         && kIsBestMatchForTruth;

        // Reco events truth-matched to a numu CC whose true vertex is outside this FV
        const Cut oofvMatch = Cut([margin](const caf::SRInteractionProxy *ixn) -> bool {
            const caf::SRTrueInteractionProxy *t = BestTruthMatch(ixn);
            return t && kTrueNumuCC(t) && !TrueVtxInFV(t, margin);
        });

        FVScanPoint p;
        p.margin  = margin;
        p.num     = std::make_unique<Spectrum>(reco()[fullSel && truthMatch], axDummy);
        p.den     = std::make_unique<Spectrum>(loader.NuTruths()[trueSignal], axTruthDummy);
        p.all     = std::make_unique<Spectrum>(reco()[fullSel],               axDummy);
        p.oofvDen = std::make_unique<Spectrum>(loader.NuTruths()[trueOOFV],   axTruthDummy);
        p.oofvLeak= std::make_unique<Spectrum>(reco()[fullSel && oofvMatch],  axDummy);
        fvScan.push_back(std::move(p));
    }

    // =========================================================================
    // MUON TRACK-LENGTH THRESHOLD SCAN
    // Pick the value where purity improves most without killing efficiency.
    // =========================================================================
    const Cut candIsTrueMuon = Cut([](const caf::SRInteractionProxy *ixn) -> bool {
        return kMuonCandidateIsTrueMuon(ixn) > 0.5;
    });

    std::vector<LenScanPoint> lenScan;
    for (double lenThr : kMuonLenThresholds) {
        const Cut hasMuonLen = Cut([lenThr](const caf::SRInteractionProxy *ixn) -> bool {
            for (const auto &p : ixn->part.dlp)
                if (kIsPrimary(&p) && kIsMuon(&p) && kPartLen(&p) > lenThr)
                    return true;
            return false;
        });

        const Cut selWithLen = kPartLenInInteractionCut_LongestTrack
                            && kVtxInFV
                            && kEventContained
                            && hasMuonLen
                            && kIsBestMatchForTruth;

        LenScanPoint p;
        p.threshold    = lenThr;
        p.num          = std::make_unique<Spectrum>(reco()[selWithLen && kNuMuTrueSignalMatch], axDummy);
        p.all          = std::make_unique<Spectrum>(reco()[selWithLen],                         axDummy);
        p.trueMuonCand = std::make_unique<Spectrum>(reco()[selWithLen && candIsTrueMuon],       axDummy);
        lenScan.push_back(std::move(p));
    }

    Spectrum sLenScanDen(loader.NuTruths()[kTrueSignal], axTruthDummy);


    // =========================================================================
    // RUN
    // =========================================================================
    loader.Go();

    // ── Headline efficiency and purity (stdout only) ─────────────────────────
    const Spectrum &sCF_BestMatch = *cutFlow.back().plain;
    const double nSelected   = sCF_BestMatch.Integral(sCF_BestMatch.POT());
    const double nSignal     = sCF_Signal.Integral(sCF_Signal.POT());
    const double nTrueSignal = sCF_TrueSignal.Integral(sCF_TrueSignal.POT());
    std::cout << "==============================================\n"
              << "Headline efficiency = " << nSignal / nTrueSignal
              << "  (" << nSignal << " / " << nTrueSignal << ")\n"
              << "Headline purity     = " << nSignal / nSelected
              << "  (" << nSignal << " / " << nSelected << ")\n"
              << "==============================================\n";

    // =========================================================================
    // SAVE
    // =========================================================================
    TFile fout(kOutputPath.c_str(), "RECREATE");
    auto save = [&fout](Spectrum &s, const std::string &name) { s.SaveTo(&fout, name.c_str()); };

    // Cut flow: sCF_<stage>, sCF_<stage>_cat, sCF_<stage>_catDetailed
    for (auto &st : cutFlow) {
        save(*st.plain,    "sCF_" + st.name);
        save(*st.cat,      "sCF_" + st.name + "_cat");
        save(*st.detailed, "sCF_" + st.name + "_catDetailed");
    }
    save(sCF_Signal,     "sCF_Signal");
    save(sCF_TrueSignal, "sCF_TrueSignal");

    // Neutrino energy
    save(sRecoNuE,                      "sRecoNuE");
    save(sNuERes2D,                     "sNuERes2D");
    save(sNuERes1D,                     "sNuERes1D");
    save(sMigrationMatrix,              "sMigrationMatrix");
    save(sEvisMigration,                "sEvisMigration");
    save(sNuERes2D_Contained,           "sNuERes2D_Contained");
    save(sNuERes2D_EscapesTMS,          "sNuERes2D_EscapesTMS");
    save(sNuERes2D_NotTrueMuon,         "sNuERes2D_NotTrueMuon");
    save(sMigrationMatrix_Contained,    "sMigrationMatrix_Contained");
    save(sMigrationMatrix_EscapesTMS,   "sMigrationMatrix_EscapesTMS");
    save(sMigrationMatrix_NotTrueMuon,  "sMigrationMatrix_NotTrueMuon");
    save(sTrueNumuCC_ENu,               "sTrueNumuCC_ENu");
    save(sRecoNumuCC_All_ENu,           "sRecoNumuCC_All_ENu");
    save(sRecoNumuCC_Matched_ENu,       "sRecoNumuCC_Matched_ENu");

    // Hadronic / inelasticity
    save(sRecoHadE,         "sRecoHadE");
    save(sRecoInelasticity, "sRecoInelasticity");
    save(sTrueInelasticity, "sTrueInelasticity");
    save(sInelRes2D,        "sInelRes2D");

    // Muon-candidate diagnostics
    save(sMuonCandLenVsIsTrueMuon,        "sMuonCandLenVsIsTrueMuon");
    save(sMuonCandLenVsTruePDG,           "sMuonCandLenVsTruePDG");
    save(sMuonToVtxDist,                  "sMuonToVtxDist");
    save(sMuonVtxDistVsIsTrueMuon,        "sMuonVtxDistVsIsTrueMuon");
    save(sMuonVtxDistVsIsTrueMuon_signal, "sMuonVtxIsTrueMuon_signal"); // name kept as-is for the plotting script
    save(sMuonVtxDistVsSignal,            "sMuonVtxDistVsSignal");

    // Interaction mode / background category
    save(sModeVsENu_selected,   "sModeVsENu_selected");
    save(sModeVsMuonP_selected, "sModeVsMuonP_selected");
    save(sModeSignalOnly,       "sModeSignalOnly");
    save(sModeBkgOnly,          "sModeBkgOnly");
    save(sBkgCatVsENu,          "sBkgCatVsENu");
    save(sBkgCatVsMuonP,        "sBkgCatVsMuonP");

    // Efficiency vs. muon angle
    save(sRecoNumuCC_MuonAngle,     "sRecoNumuCC_MuonAngle");
    save(sTrueNumuCC_MuonAngle,     "sTrueNumuCC_MuonAngle");
    save(sRecoNumuCC_All_MuonAngle, "sRecoNumuCC_All_MuonAngle");

    // Efficiency vs. vertex position
    save(sRecoNumuCC_VtxX, "sRecoNumuCC_VtxX");
    save(sTrueNumuCC_VtxX, "sTrueNumuCC_VtxX");
    save(sRecoNumuCC_VtxY, "sRecoNumuCC_VtxY");
    save(sTrueNumuCC_VtxY, "sTrueNumuCC_VtxY");
    save(sRecoNumuCC_VtxZ, "sRecoNumuCC_VtxZ");
    save(sTrueNumuCC_VtxZ, "sTrueNumuCC_VtxZ");

    // FV scan (margin as suffix)
    for (auto &p : fvScan) {
        const std::string tag = std::to_string((int)p.margin);
        save(*p.num,      "sFVScan_Num_ZUp"     + tag);
        save(*p.den,      "sFVScan_Den_ZUp"     + tag);
        save(*p.all,      "sFVScan_All_ZUp"     + tag);
        save(*p.oofvDen,  "sFVScan_OOFVDen_ZUp" + tag);
        save(*p.oofvLeak, "sFVScan_OOFVLeak_ZUp"+ tag);
    }

    // Muon length scan (threshold as suffix)
    save(sLenScanDen, "sLenScan_Den");
    for (auto &p : lenScan) {
        const std::string tag = std::to_string((int)p.threshold);
        save(*p.num,          "sLenScan_Num_"          + tag);
        save(*p.all,          "sLenScan_All_"          + tag);
        save(*p.trueMuonCand, "sLenScan_TrueMuonCand_" + tag);
    }

    fout.Close();
    std::cout << "Saved " << kOutputPath << std::endl;
    
}
