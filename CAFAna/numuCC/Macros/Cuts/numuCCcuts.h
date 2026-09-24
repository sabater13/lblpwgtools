#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/TruthMatching.h"


#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#define NDLArXLo -346.9
#define NDLArXHi  346.9
#define NDLArYLo -215.5
#define NDLArYHi   81.7
#define NDLArZLo  418.2
#define NDLArZHi  913.3


// Returns pointer to best-matched true interaction, or nullptr if no match.
const caf::SRTrueInteractionProxy* GetBestTruthMatchInt(const caf::SRInteractionProxy *ixn) {
    if (ixn->truth.empty()) return nullptr;
    const caf::SRProxy *sr = ixn->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float maxOverlap = -1.0;
    for (size_t i = 0; i < ixn->truthOverlap.size(); ++i)
        if (ixn->truthOverlap[i] > maxOverlap) { tidx = i; maxOverlap = ixn->truthOverlap[i]; }
    if (maxOverlap < 0.0) return nullptr;
    return caf::FindInteraction(sr->mc, ixn->truth[tidx]);
}
// Particle level truth-matching helper (like interaction level above)
const caf::SRTrueParticleProxy* GetBestTruthMatchPart(const caf::SRRecoParticleProxy *p) {
    if (p->truth.empty()) return nullptr;
    const caf::SRProxy *sr = p->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float maxOverlap = -1.0;
    for (size_t i = 0; i < p->truthOverlap.size(); ++i)
        if (p->truthOverlap[i] > maxOverlap) { tidx = i; maxOverlap = p->truthOverlap[i]; }
    if (maxOverlap < 0.0) return nullptr;
    return caf::FindParticle(sr->mc, p->truth[tidx]);
}



// VERTEX
// reco vertices
const Var kVtxX = SIMPLEVAR(vtx.x);
const Var kVtxY = SIMPLEVAR(vtx.y);
const Var kVtxZ = SIMPLEVAR(vtx.z);
// true vertices
const TruthVar kTrueVtxX = SIMPLETRUTHVAR(vtx.x);
const TruthVar kTrueVtxY = SIMPLETRUTHVAR(vtx.y);
const TruthVar kTrueVtxZ = SIMPLETRUTHVAR(vtx.z);
// truth matched vertices
const Var kMatchedTrueVtxX([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueVtxX(tixn);
});
const Var kMatchedTrueVtxY([](const caf::SRInteractionProxy *ixn) {
   auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueVtxY(tixn);
});
const Var kMatchedTrueVtxZ([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueVtxZ(tixn);
});

// START AND END POSITION
//Truth
const TruthPartVar kTruePartStartX = SIMPLETRUTHPARTVAR(start_pos.x);
const TruthPartVar kTruePartStartY = SIMPLETRUTHPARTVAR(start_pos.y);
const TruthPartVar kTruePartStartZ = SIMPLETRUTHPARTVAR(start_pos.z);
const TruthPartVar kTruePartEndX = SIMPLETRUTHPARTVAR(end_pos.x);
const TruthPartVar kTruePartEndY = SIMPLETRUTHPARTVAR(end_pos.y);
const TruthPartVar kTruePartEndZ = SIMPLETRUTHPARTVAR(end_pos.z);

// GENERAL
const RecoPartVar kPartEndX = SIMPLEPARTVAR(end.x);
const RecoPartVar kPartEndY = SIMPLEPARTVAR(end.y);
const RecoPartVar kPartEndZ = SIMPLEPARTVAR(end.z);
const RecoPartCut kIsPrimary = SIMPLEPARTVAR(primary) == 1;
const RecoPartVar kPartPDG = SIMPLEPARTVAR(pdg);
const RecoPartCut kIsMuon = kPartPDG == 13 || kPartPDG == -13;
const RecoPartCut kIsElectron = kPartPDG == 11 || kPartPDG == -11;
const RecoPartCut kIsPhoton = kPartPDG == 22;
const RecoPartCut kIsPrimaryMuon = kIsMuon && kIsPrimary;
const RecoPartCut kIsPrimaryPion = (kPartPDG == 211 || kPartPDG == -211) && kIsPrimary;
const RecoPartCut kIsPrimaryProton = (kPartPDG == 2212) && kIsPrimary;
const RecoPartCut kPartEscapesToTMS = kPartEndX < NDLArXHi - 25 && kPartEndX > NDLArXLo + 25
                                    && kPartEndY < NDLArYHi - 25 && kPartEndY > NDLArYLo + 25
                                    && kPartEndZ > NDLArZLo + 25;
const RecoPartCut kIsPrimaryMuonEscaping = kIsPrimaryMuon && kPartEscapesToTMS;

// CONTAINMENT
// for now use is_contained, variable from SPINE
// Looks like it is defined as vtx and all enegy contained
// But don't know which volume it is defined with, maybe ND-LAr active volume?
const RecoPartCut kPartNDLArContained = SIMPLEPARTVAR(contained) == 1;
const RecoPartCut kPartNDLArContainedExceptMuonDownstream = kPartNDLArContained || kIsPrimaryMuonEscaping;
// all non-primary particles from interaction contianed
const Cut kAllPartContained([](const caf::SRInteractionProxy *ixn){
    for(const auto &p: ixn->part.dlp){
        if((!kIsPrimary && !kPartNDLArContained)(&p)){
            return false;
        }
    }
    return true;
});
// vertex contained in ND-LAr active volume
const Cut kVtxContained = kVtxX > NDLArXLo && kVtxX < NDLArXHi
                        && kVtxY > NDLArYLo && kVtxY < NDLArYHi
                        && kVtxZ > NDLArZLo && kVtxZ < NDLArZHi;
// full event contained
const Cut kEventContained([](const caf::SRInteractionProxy * ixn){
    if (!kVtxContained(ixn)) return false;
    for(const auto &p: ixn->part.dlp){
        if ((!kPartNDLArContainedExceptMuonDownstream)(&p)){  // kPartNDLArContainedExceptMuonDownstream OR kPartNDLArContained
            return false;
        }
    }
    return true;
});
// true containement
const TruthPartCut kTruePartContained = kTruePartEndX < NDLArXHi - 5 && kTruePartEndX > NDLArXLo + 5
                                     && kTruePartEndY < NDLArYHi - 5 && kTruePartEndY > NDLArYLo + 5
                                     && kTruePartEndZ < NDLArZHi - 5 && kTruePartEndZ > NDLArZLo + 5;
const TruthPartCut kTruePartContainedButCanEscapeToTMS = kTruePartEndX < NDLArXHi - 25 && kTruePartEndX > NDLArXLo + 25
                                                        && kTruePartEndY < NDLArYHi - 25 && kTruePartEndY > NDLArYLo + 25
                                                        && kTruePartEndZ > NDLArZLo + 25;

const TruthPartCut kTrueMuon = SIMPLETRUTHPARTVAR(pdg) == 13;
const TruthPartCut kTrueMuonEscapesToTMS = kTruePartEndX < NDLArXHi - 25 && kTruePartEndX > NDLArXLo + 25
                                        && kTruePartEndY < NDLArYHi - 25 && kTruePartEndY > NDLArYLo + 25
                                        && kTruePartEndZ > NDLArZLo + 25;
const TruthPartCut kTrueMuonEscaping = kTrueMuon && kTrueMuonEscapesToTMS;
const TruthPartCut kTrueContainedExceptMuonDownstream = kTruePartContained || kTrueMuonEscaping;
const TruthCut kAllTrueContained([](const caf::SRTrueInteractionProxy *truth){
  for(const auto &p: truth->prim) {
    if(kTrueMuonEscaping(&p)) continue;
    if(!kTruePartContained(&p)) return false;
  }
  return true;
});
const TruthCut kAllTrueContained_exceptMuons([](const caf::SRTrueInteractionProxy *truth){
  for(const auto &p: truth->prim) {
    if(kTrueMuonEscaping(&p)) continue;
    if(!kTruePartContained(&p)) return false;
  }
  return true;
});
const TruthCut kAllTrueContainedExceptMuonDownstream([](const caf::SRTrueInteractionProxy *truth){
  for(const auto &p: truth->prim){
    if(!kTrueContainedExceptMuonDownstream(&p)) return false;
  }
  return true;
});

const TruthCut kAllTrueContainedNDLArAndEscapeToTMS([](const caf::SRTrueInteractionProxy *truth){
  for(const auto &p: truth->prim)
    if(!kTruePartContainedButCanEscapeToTMS(&p)) return false;
  return true;
});

// DIFFERENT VOLUMES
// start with original FV cut (i.e. 25 cm form wall edge) and then vary Z axis
const Cut kVtxInFV = kVtxX > NDLArXLo + 25 && kVtxX < NDLArXHi - 25
                    && kVtxY > NDLArYLo + 25 && kVtxY < NDLArYHi - 25
                    && kVtxZ > NDLArZLo + 25 && kVtxZ < NDLArZHi - 25;
const TruthCut kTrueVtxInFV = kTrueVtxX > NDLArXLo + 25 && kTrueVtxX < NDLArXHi - 25
                             && kTrueVtxY > NDLArYLo + 25 && kTrueVtxY < NDLArYHi - 25
                             && kTrueVtxZ > NDLArZLo + 25 && kTrueVtxZ < NDLArZHi - 25;

// qo AND q3 VARIABLES
// In truth for now 
const TruthVar kq0([](const caf::SRTrueInteractionProxy * sr){
    return sr->q0;
});
const TruthVar kq3([](const caf::SRTrueInteractionProxy * sr){
    return sr->modq;
});
const Var kTruthMatchedQ0([](const caf::SRInteractionProxy * ixn){
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kq0(tixn);
});
const Var kTruthMatchedQ3([](const caf::SRInteractionProxy * ixn){
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kq3(tixn);
});

// QUALITY CUT
const RecoPartVar kPartLen([](const caf::SRRecoParticleProxy * p) -> double {
    // if (isnan(p->start.x)) return -1.0;
    // if (isnan(p->start.y)) return -1.0;
    // if (isnan(p->start.z)) return -1.0;
    // if (isnan(p->end.x)) return -1.0;
    // if (isnan(p->end.y)) return -1.0;
    // if (isnan(p->end.z)) return -1.0;
    // if (isinf(p->start.x)) return -1.0;
    // if (isinf(p->start.y)) return -1.0;
    // if (isinf(p->start.z)) return -1.0;
    // if (isinf(p->end.x)) return -1.0;
    // if (isinf(p->end.y)) return -1.0;
    // if (isinf(p->end.z)) return -1.0;
    return std::hypot(p->end.x - p->start.x, p->end.y - p->start.y, p->end.z - p->start.z);
});
// ── Length of the tagged primary muon candidate specifically (not just longest track) ──
const Var kMuonCandidateLen([](const caf::SRInteractionProxy *ixn) -> double {
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p) && kIsMuon(&p))
            return kPartLen(&p);
    return -1.0;
});


// -- Distance [cm] between a reco particle's start point and the reco interaction vtx --
double PartStartToVtxDist(const caf::SRInteractionProxy * ixn, const caf::SRRecoParticleProxy * p){
    const double dx = (double)p->start.x - (double)ixn->vtx.x;
    const double dy = (double)p->start.y - (double)ixn->vtx.y;
    const double dz = (double)p->start.z - (double)ixn->vtx.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// -- Distance between the start of the primary muon candidate and the reco vertex --
const Var kMuonStartToVtxDist([](const caf::SRInteractionProxy * ixn){
    for (const auto &p : ixn->part.dlp){
        if (!kIsPrimary(&p) || !kIsMuon(&p)) continue;
        const double d = PartStartToVtxDist(ixn, &p);
        return std::isfinite(d) ? d : 1.0;
    }
    return -1.0; // no muon candidate
});


// ── Is the tagged muon candidate actually truth-matched to a muon? ──
//    1 = yes (correct PID), 0 = no (some other particle), -1 = no muon candidate/no truth
const Var kMuonCandidateIsTrueMuon([](const caf::SRInteractionProxy *ixn) -> double {
    for (const auto &p : ixn->part.dlp) {
        if (!kIsPrimary(&p) || !kIsMuon(&p)) continue;
        if (p.truth.empty()) return -1.;
        const caf::SRProxy *sr = ixn->Ancestor<caf::SRProxy>();
        size_t tidx = 0; float maxOverlap = 0;
        for (size_t i = 0; i < p.truthOverlap.size(); ++i)
            if (p.truthOverlap[i] > maxOverlap) { tidx = i; maxOverlap = p.truthOverlap[i]; }
        int pdg = std::abs(caf::FindParticle(sr->mc, p.truth[tidx])->pdg);
        return (pdg == 13) ? 1. : 0.;
    }
    return -1.; // no muon candidate at all
});

// ── Is the primary muon candidate contained in the ND-LAr active volume? ──
const Cut kMuonCandContained([](const caf::SRInteractionProxy *ixn) -> bool {
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p) && kIsMuon(&p))
            return kPartNDLArContained(&p);
    return false; // no muon candidate found
});

// ── Does the primary muon candidate escape into the TMS? ──
const Cut kMuonCandEscapesToTMS([](const caf::SRInteractionProxy *ixn) -> bool {
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p) && kIsMuon(&p))
            return kPartEscapesToTMS(&p);
    return false;
});

// ── Is the primary muon candidate actually NOT a true muon (mis-PID)? ──
//    kMuonCandidateIsTrueMuon: 1 = correct PID, 0 = mis-PID, -1 = no candidate/no truth.
const Cut kMuonCandIsNotTrueMuon([](const caf::SRInteractionProxy *ixn) -> bool {
    return kMuonCandidateIsTrueMuon(ixn) == 0.;
});


const Cut kPartLenInInteractionCut_LongestTrack([](const caf::SRInteractionProxy * ixn){
    if (ixn->part.dlp.empty()) return false;
    size_t longestIdx = 1e6;
    double longest = -1;
    for (size_t i = 0; i < ixn->part.dlp.size(); ++i){
        const auto &p = ixn->part.dlp[i];
        if (kIsPrimary(&p) && kPartLen(&p) > longest){
            longestIdx = i;
            longest = kPartLen(&p);
        }
    }
    if (longestIdx > ixn->part.dlp.size()) return false;
    return longest > 0.1;
});

// TRUE NUMU CC EVENT
const TruthCut kTrueNumuCC = SIMPLETRUTHVAR(iscc) == 1 && (SIMPLETRUTHVAR(pdg) == 14 || SIMPLETRUTHVAR(pdg) == -14);
// contained true numu CC
const TruthCut kTrueNumuCCContained = kTrueNumuCC && kTrueVtxInFV && kAllTrueContained;
const TruthCut kTrueNumuCCContainedNDLArAndEscapeToTMS = kTrueNumuCC && kTrueVtxInFV && kAllTrueContainedExceptMuonDownstream; //kAllTrueContained;

// TRUTH MATCHED NUMU CC EVENT
const Cut kNuMuTrueSignalMatch([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueNumuCCContainedNDLArAndEscapeToTMS(tixn);
});
const Cut kNuMuTrueSignalMatch_trueNumuCCInFVUncontained([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueNumuCC(tixn) && kTrueVtxInFV(tixn) && !kAllTrueContainedExceptMuonDownstream(tixn);
});
const Cut kNuMuTrueSignalMatch_trueNumuCCOutFV([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueNumuCC(tixn) && !kTrueVtxInFV(tixn);
});

Cut MakeNuMuTrueSignalMatch(const TruthCut& trueFV){
    return Cut([trueFV](const caf::SRInteractionProxy *ixn) -> bool {
        auto *tixn = GetBestTruthMatchInt(ixn);
        TruthCut signal = kTrueNumuCC && trueFV && kAllTrueContainedExceptMuonDownstream;
        return signal(tixn);
    });
}

Cut MakeNuMuTrueMatchNoContainment(const TruthCut& trueFV){
    return Cut([trueFV](const caf::SRInteractionProxy *ixn) -> bool {
        auto *tixn = GetBestTruthMatchInt(ixn);
        return (kTrueNumuCC && trueFV)(tixn);
    });
}

Cut OutsideFVMatch(const TruthCut& trueFV){
    return Cut([trueFV](const caf::SRInteractionProxy *ixn){
        auto *tixn = GetBestTruthMatchInt(ixn);
        return kTrueNumuCC(tixn) && !trueFV(tixn);
    });
}

TruthCut TrueOutsideFV_Z(const TruthCut& fullFV){
    return TruthCut([fullFV](const caf::SRTrueInteractionProxy* t){
        if(fullFV(t)) return false;

        // only fail Z
        bool failZ = !(kTrueVtxZ(t) > NDLArZLo + 25 && kTrueVtxZ(t) < NDLArZHi - 25);

        return failZ;
    });
}
Cut OutsideFVMatch_Z(const TruthCut& fullFV){
    TruthCut zCut = TrueOutsideFV_Z(fullFV);
    return Cut([zCut](const caf::SRInteractionProxy *ixn){
        auto *tixn = GetBestTruthMatchInt(ixn);
        return (kTrueNumuCC && zCut)(tixn);
    });
}

TruthCut TrueOutsideFV_Z_upstream(const TruthCut& fullFV){
    return TruthCut([fullFV](const caf::SRTrueInteractionProxy* t){
        if(fullFV(t)) return false;
        return kTrueVtxZ(t) < NDLArZLo + 25;
    });
}
Cut OutsideFVMatch_Z_upstream(const TruthCut& fullFV){
    TruthCut zCut = TrueOutsideFV_Z_upstream(fullFV);
    return Cut([zCut](const caf::SRInteractionProxy *ixn){
        auto *tixn = GetBestTruthMatchInt(ixn);
        return (kTrueNumuCC && zCut)(tixn);
    });
}

TruthCut TrueOutsideFV_Z_downstream(const TruthCut& fullFV){
    return TruthCut([fullFV](const caf::SRTrueInteractionProxy* t){
        if(fullFV(t)) return false;
        return kTrueVtxZ(t) > NDLArZHi - 25;
    });
}
Cut OutsideFVMatch_Z_downstream(const TruthCut& fullFV){
    TruthCut zCut = TrueOutsideFV_Z_downstream(fullFV);
    return Cut([zCut](const caf::SRInteractionProxy *ixn){
        auto *tixn = GetBestTruthMatchInt(ixn);
        return (kTrueNumuCC && zCut)(tixn);
    });
}

TruthCut TrueOutsideFV_XY(const TruthCut& fullFV){
    return TruthCut([fullFV](const caf::SRTrueInteractionProxy* t){
        if(fullFV(t)) return false;

        bool inZ = (kTrueVtxZ(t) > NDLArZLo + 25 && kTrueVtxZ(t) < NDLArZHi - 25);

        bool failX = !(kTrueVtxX(t) > NDLArXLo + 25 && kTrueVtxX(t) < NDLArXHi - 25);
        bool failY = !(kTrueVtxY(t) > NDLArYLo + 25 && kTrueVtxY(t) < NDLArYHi - 25);

        return inZ && (failX || failY);
    });
}
Cut OutsideFVMatch_XY(const TruthCut& fullFV){
    TruthCut xyCut = TrueOutsideFV_XY(fullFV);
    return Cut([xyCut](const caf::SRInteractionProxy *ixn){
        auto *tixn = GetBestTruthMatchInt(ixn);
        return (kTrueNumuCC && xyCut)(tixn);
    });
}

// RECO NUMU CC EVENT
const Cut kRecoNumuCC([](const caf::SRInteractionProxy * ixn){
    if (ixn->part.dlp.empty()) return false;
    for (size_t i = 0; i < ixn->part.dlp.size(); i++){
        const auto &p = ixn->part.dlp[i];
        if (kIsPrimary(&p) && kIsMuon(&p)){
            return true;
        }
    }
    return false;
});

const Cut kPrimMuonTrackLengthCut(double length){
    const Cut kRecoNumuCC_TLbt([length](const caf::SRInteractionProxy * ixn){
        if (ixn->part.dlp.empty()) return false;
        for (size_t i = 0; i < ixn->part.dlp.size(); i++){
            const auto &p = ixn->part.dlp[i];
            if (kIsPrimary(&p) && kIsMuon(&p)){
                return kPartLen(&p) >= length;
            }
        }
        return false;
    });
    return kRecoNumuCC_TLbt;
}


// TRUTH OVERLAP
const Cut OverlapCut(double overlapThreshold){
  const Cut kOverlapCut([overlapThreshold](const caf::SRInteractionProxy * sr)
  {
    if (sr->truthOverlap.empty()) return false;
    return *std::max_element(sr->truthOverlap.begin(), sr->truthOverlap.end()) > overlapThreshold;
  });
  return kOverlapCut;
}

const Cut kIsBestMatchForTruth([](const caf::SRInteractionProxy* ixn)
{
    if(ixn->truth.empty()) return false;
    const caf::SRProxy* sr = ixn->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float myOverlap = -1.0;
    for(size_t i = 0; i < ixn->truthOverlap.size(); ++i){
        if(ixn->truthOverlap[i] > myOverlap){
            tidx = i;
            myOverlap = ixn->truthOverlap[i];
        }
    }
    int myTruthID = ixn->truth[tidx];
    for(const auto& other : sr->common.ixn.dlp){
        if(&other == ixn) continue;
        if(other.truth.empty()) continue;
        size_t oidx = 0;
        float otherOverlap = -1.0;
        for(size_t j = 0; j < other.truthOverlap.size(); ++j){
            if(other.truthOverlap[j] > otherOverlap){
                otherOverlap = other.truthOverlap[j];
                oidx = j;
            }
        }
        if(other.truth[oidx] == myTruthID){
            if(otherOverlap > myOverlap)
                return false;
        }
    }
    return true;
});



// true PDG of reco primary particles that have a track length of less than 1 mm
const RecoPartVar kTruePDGofShortPrimTrackLength([](const caf::SRRecoParticleProxy * part){
    if (part->truth.empty()) return -1.0;
    const caf::SRProxy * sr = part->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float maxOverlap = 0;
    for (size_t i =0; i < part->truthOverlap.size(); i++){
        if (part->truthOverlap[i] > maxOverlap){
            tidx = i;
            maxOverlap = part->truthOverlap[i];
        }
    }
    if (!kIsPrimary(part)) return -1.0;
    if (kPartLen(part) < 0.1) return -1.0;
    int pdg = caf::FindParticle(sr->mc, part->truth[tidx])->pdg;
    switch(pdg){
        case 13: return 0.; // muon
        case -13: return 1.; // mu+
        case 11: return 2.; // electron
        case -11: return 3.; // positron
        case 22: return 4.; // photon
        case 211: return 5.; // pi+
        case -211: return 6.; // pi-
        case 111: return 7.; // pi0
        case 2212: return 8.; // proton
        case 321: return 9.; // K+
        case -321: return 10.; // K-
        case 130: return 11.; // K0
        case 310: return 11.; // K0
        case 311: return 11.; // K0
        case 2112: return 12.; // neutron
    }
    return 13.;
});
std::vector<std::string> PDGlabels = {"#mu^{-}", "#mu^{+}", "e^{-}", "e^{+}","#gamma","#pi^{+}","#pi^{-}","#pi^{0}","p","K^{+}","K^{-}","K^{0}", "n","other"};

// particle level containement variables (for distance to wall vs track length)
// distance to closest wall
double DistanceToWall(double x, double y, double z)
{
  double dx = std::min(x - NDLArXLo, NDLArXHi - x);
  double dy = std::min(y - NDLArYLo, NDLArYHi - y);
  double dz = std::min(z - NDLArZLo, NDLArZHi - z);

  return std::min({dx, dy, dz});
}
const RecoPartVar kDistToClosestWall_PartStart([](const caf::SRRecoParticleProxy * part){
    double distance = -1.0;
    double part_start_x = part->start.x;
    double part_start_y = part->start.y;
    double part_start_z = part->start.z;
    distance = DistanceToWall(part_start_x, part_start_y, part_start_z);
    return distance;
});
const RecoPartVar kDistToClosestWall_PartEnd([](const caf::SRRecoParticleProxy * part){
    double distance = -1.0;
    if (isnan(part->end.x)) return -1.0;
    double part_end_x = part->end.x;
    if (isnan(part->end.y)) return -1.0;
    double part_end_y = part->end.y;
    if (isnan(part->end.z)) return -1.0;
    double part_end_z = part->end.z;
    distance = DistanceToWall(part_end_x, part_end_y, part_end_z);
    return distance;
});


const RecoPartVar kMomentum([](const caf::SRRecoParticleProxy * part){
    double momentum = -1.0;
    return std::sqrt(part->p.x * part->p.x +
                     part->p.y * part->p.y +
                     part->p.z * part->p.z);
});

const auto beam_dir = TVector3(0.0,-0.05836,1.0);
TVector3 RecoPartDir(const caf::SRRecoParticleProxy * part){return TVector3(part->end.x, part->end.y, part->end.z) - TVector3(part->start.x, part->start.y, part->start.z);}
const RecoPartVar kPartAngle([](const caf::SRRecoParticleProxy * part){
    auto angle = RecoPartDir(part).Angle(beam_dir) * 180.0 / TMath::Pi();
    return angle;
});

double PDGMass(int pdg) {
    switch (std::abs(pdg)) {
        case 13:   return 0.10566; // muon
        case 211:  return 0.13957; // pion+/-
        case 2212: return 0.93827; // proton
        case 11:   return 0.000511; // electron
        case 22:   return 0.0;     // photon
        case 111:  return 0.13498; // pi0
        case 321:  return 0.49368; // kaon+/-
        case 2112: return 0.93957; // neutron
        default:   return 0.13957; // default to pion mass
    }
}
const RecoPartVar kPartKE([](const caf::SRRecoParticleProxy * part) -> double {
    int pdg = std::abs(part->pdg);
    double mass = PDGMass(pdg);
    double E = part->E;
    return E - mass;
});

const RecoPartVar kPartE = SIMPLEPARTVAR(E);
const Var kSumPrimaryPartE([](const caf::SRInteractionProxy *ixn) {
    double E = 0;
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p))
            E += kPartE(&p);
    return E;
});
const Var kTruthMatchedNuE([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return (double)tixn->E;
});
const TruthVar kTrueNuE([](const caf::SRTrueInteractionProxy *t) {
    return t->E;
});

// ── True visible energy: sum of true final-state primary energies,      ──
//    excluding neutrons (which typically leave no visible signal in LAr).
//    This is the truth-level analog of kRecoNuE_Calo, which sums the
//    energies of *reconstructed* (i.e. detector-visible) primary particles.
const TruthVar kTrueEvis([](const caf::SRTrueInteractionProxy *t) -> double {
    double E = 0;
    for (const auto &p : t->prim) {
        if (std::abs(p.pdg) == 2112) continue; // skip neutrons: invisible in LAr
        E += p.p.E - PDGMass(p.pdg); // Kinetic energy, consistent with kPartE = SIMPLEPARTVAR(E) used in kRecoNuE_Calo
    }
    return E;
});

// ── Truth-matched true visible energy, for use on the reco (Interactions) loop ──
const Var kTruthMatchedEvis([](const caf::SRInteractionProxy *ixn) -> double {
    auto *tixn = GetBestTruthMatchInt(ixn);
    return kTrueEvis(tixn);
});

// ── Background category for selected (kBaseSelection) events ──
//    0=NC, 1=nueCC, 2=nutauCC, 3=numuCC true-OOFV, 4=numuCC true-uncontained,
//    5=numuCC signal (sanity bin - should be ~0 if used only on background subset),
//    6=no truth match at all (spurious/unmatched reco, e.g. cosmic or noise)
const Var kBkgCategory([](const caf::SRInteractionProxy *ixn) -> double {
    if (ixn->truth.empty()) return 6.;
    auto *tixn = GetBestTruthMatchInt(ixn);

    if (!tixn->iscc)                 return 0.; // NC
    if (std::abs(tixn->pdg) == 12)   return 1.; // nue CC
    if (std::abs(tixn->pdg) == 16)   return 2.; // nutau CC
    // true numu CC from here on
    if (!kTrueVtxInFV(tixn))                          return 3.; // OOFV
    if (!kAllTrueContainedExceptMuonDownstream(tixn))  return 4.; // uncontained
    return 5.; // true signal (sanity check bin)
});
std::vector<std::string> BkgCategoryLabels = {
    "NC", "#nu_{e} CC", "#nu_{#tau} CC",
    "#nu_{#mu} CC OOFV", "#nu_{#mu} CC Uncontained",
    "#nu_{#mu} CC Signal", "No truth match"
};

const Var kTruthMatchedTrueMuonP([](const caf::SRInteractionProxy *ixn) {
    auto *tixn = GetBestTruthMatchInt(ixn);
    for (const auto &p : tixn->prim)
        if (std::abs(p.pdg) == 13)
            return (double)std::sqrt(p.p.px*p.p.px + p.p.py*p.p.py + p.p.pz*p.p.pz);
    return -1.0;
});
const TruthVar kTrueMuonP([](const caf::SRTrueInteractionProxy *t) {
    for (const auto &p : t->prim)
        if (std::abs(p.pdg) == 13)
            return (double)std::sqrt(p.p.px*p.p.px + p.p.py*p.p.py + p.p.pz*p.p.pz);
    return -1.0;
});

const TruthPartVar kTruePartMomentum([](const caf::SRTrueParticleProxy *p) -> double {
    return std::sqrt(p->p.px*p->p.px + p->p.py*p->p.py + p->p.pz*p->p.pz);
});

const auto beam_dir_truth = TVector3(0.0, -0.05836, 1.0);
const TruthPartVar kTruePartAngle([](const caf::SRTrueParticleProxy *p) -> double {
    TVector3 dir(p->p.px, p->p.py, p->p.pz);
    return dir.Angle(beam_dir_truth) * 180.0 / TMath::Pi();
});

const TruthPartVar kTruePartKE([](const caf::SRTrueParticleProxy *p) -> double {
    double px = p->p.px, py = p->p.py, pz = p->p.pz;
    double pmag = std::sqrt(px*px + py*py + pz*pz);
    double mass = PDGMass(std::abs(p->pdg));
    double E = std::sqrt(pmag*pmag + mass*mass);
    return E - mass;
});

const TruthPartCut kTruePartIsMuon   = SIMPLETRUTHPARTVAR(pdg) == 13 || SIMPLETRUTHPARTVAR(pdg) == -13;
const TruthPartCut kTruePartIsPion   = SIMPLETRUTHPARTVAR(pdg) == 211 || SIMPLETRUTHPARTVAR(pdg) == -211;
const TruthPartCut kTruePartIsProton = SIMPLETRUTHPARTVAR(pdg) == 2212;

const Var kTruePDGofRecoPrimaryMuon([](const caf::SRInteractionProxy *ixn) -> double {
    for (const auto &p : ixn->part.dlp) {
        if (!kIsPrimary(&p) || !kIsMuon(&p)) continue;
        if (p.truth.empty()) return 13.; // "other" bin
        const caf::SRProxy *sr = ixn->Ancestor<caf::SRProxy>();
        size_t tidx = 0; float maxOverlap = 0;
        for (size_t i = 0; i < p.truthOverlap.size(); ++i)
            if (p.truthOverlap[i] > maxOverlap) { tidx = i; maxOverlap = p.truthOverlap[i]; }
        int pdg = std::abs(caf::FindParticle(sr->mc, p.truth[tidx])->pdg);
        switch(pdg){
            case 13:   return 0.;
            case 211:  return 1.;  // pi+/-
            case 2212: return 2.;  // proton
            case 11:   return 3.;  // electron
            case 22:   return 4.;  // photon
            case 2112: return 5.;  // neutron
            case 321:  return 6.;  // kaon
        }
        return 7.; // other
    }
    return -1.; // no reco muon found
});
std::vector<std::string> RecoPDGlabels = {"#mu", "#pi^{#pm}", "p", "e", "#gamma", "n", "K", "other"};
const RecoPartVar kPartTruePDGofRecoPrimaryMuon([](const caf::SRRecoParticleProxy* p) -> double {
    if (p->truth.empty()) return -1.;
    const caf::SRProxy* sr = p->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float maxOverlap = 0;
    for (size_t i = 0; i < p->truthOverlap.size(); ++i) {
        if (p->truthOverlap[i] > maxOverlap) {
            tidx = i;
            maxOverlap = p->truthOverlap[i];
        }
    }
    int pdg = std::abs(caf::FindParticle(sr->mc, p->truth[tidx])->pdg);
    switch (pdg) {
        case 13:   return 0.; // muon
        case 211:  return 1.; // pion
        case 2212: return 2.; // proton
        case 11:   return 3.; // electron
        case 22:   return 4.; // photon
        case 2112: return 5.; // neutron
        case 321:  return 6.; // kaon
    }
    return 7.; // other
});

const Var kTrueIntMode([](const caf::SRInteractionProxy *ixn) -> double {
    auto *tixn = GetBestTruthMatchInt(ixn);
    // iscc==0 -> NC, iscc==1 -> CC; pdg tells neutrino flavour
    if (!tixn->iscc)                                 return 0.; // NC (any flavour)
    if (std::abs(tixn->pdg) == 12)                   return 1.; // nue CC
    if (std::abs(tixn->pdg) == 16)                   return 2.; // nutau CC
    // numu CC but outside signal definition (wrong FV or containment)
    return 3.; // numu CC out-of-FV / uncontained
});
std::vector<std::string> IntModeLabels = {"NC", "#nu_{e} CC", "#nu_{#tau} CC", "#nu_{#mu} CC OOF/uncontained"};

const TruthVar kTrueMuonEndZ([](const caf::SRTrueInteractionProxy *t) -> double {
    for (const auto &p : t->prim)
        if (std::abs(p.pdg) == 13)
            return (double)p.end_pos.z;
    return -9999.;
});
const TruthVar kTrueMuonEndX([](const caf::SRTrueInteractionProxy *t) -> double {
    for (const auto &p : t->prim)
        if (std::abs(p.pdg) == 13)
            return (double)p.end_pos.x;
    return -9999.;
});

// Cut on truth-matched neutrino energy [GeV]
Cut MakeNuEWindow(double eMin, double eMax) {
    return Cut([eMin, eMax](const caf::SRInteractionProxy *ixn) -> bool {
        double enu = kTruthMatchedNuE(ixn);
        return enu > eMin && enu < eMax;
    });
}
const Cut kNuE_2p5_3p0 = MakeNuEWindow(2.5, 3.0);

const Var kShortestPrimaryTrackLen([](const caf::SRInteractionProxy *ixn) -> double {
    double shortest = 1e6;
    bool found = false;
    for (const auto &p : ixn->part.dlp) {
        if (!kIsPrimary(&p)) continue;
        double len = kPartLen(&p);
        if (len < 0.1) continue;  // same threshold as your existing quality cut
        if (len < shortest) {
            shortest = len;
            found = true;
        }
    }
    return found ? shortest : -1.0;
});
const Var kLongestPrimaryTrackLen([](const caf::SRInteractionProxy *ixn) -> double {
    double longest = -1.0;
    bool found = false;
    for (const auto &p : ixn->part.dlp) {
        if (!kIsPrimary(&p)) continue;
        double len = kPartLen(&p);
        if (len < 0.1) continue;  // same threshold as your existing quality cut
        if (len > longest) {
            longest = len;
            found = true;
        }
    }
    return found ? longest : -1.0;
});

const Var kNPrimaryRecoParticles([](const caf::SRInteractionProxy *ixn) -> double {
    int n = 0;
    for (const auto &p : ixn->part.dlp)
        if (kIsPrimary(&p) && kPartLen(&p) > 0.1) ++n;
    return n;
});

const Var kNRecoParticles([](const caf::SRInteractionProxy *ixn) -> double {
    int n = 0;
    for (const auto &p : ixn->part.dlp)
        n = n + 1;
    return n;
});

const TruthVar kTrueNuMomentum([](const caf::SRTrueInteractionProxy * tixn) -> double {
    return std::sqrt(tixn->momentum.x * tixn->momentum.x + tixn->momentum.y * tixn->momentum.y + tixn->momentum.z * tixn->momentum.z);
});

const TruthVar kTrueNuAngle([](const caf::SRTrueInteractionProxy *tixn) -> double {
    TVector3 dir(tixn->momentum.x, tixn->momentum.y, tixn->momentum.z);
    return dir.Angle(beam_dir_truth) * 180.0 / TMath::Pi();
});

const TruthPartVar kTruePDG = SIMPLETRUTHPARTVAR(pdg);
const TruthPartVar kTruePartE = SIMPLETRUTHPARTVAR(p.E);
// true muon KE with matched reconstructed muon
const RecoPartVar kMatchedTruePartKE_TrueMuonCut([](const caf::SRRecoParticleProxy * part){
    if(part->truth.empty()) return -1.0;
    const caf::SRProxy *sr = part->Ancestor<caf::SRProxy>();
    size_t tidx = 0;
    float maxOverlap = 0;
    for(size_t i = 0; i < part->truthOverlap.size(); ++i) {
        if(part->truthOverlap[i] > maxOverlap) {
        tidx = i;
        maxOverlap = part->truthOverlap[i];
        }
    }
    const caf::SRTrueParticleProxy *tixn = caf::FindParticle(sr->mc, part->truth[tidx]);
    if (kTruePDG(tixn) != 13) return -1000.0;
    return kTruePartE(tixn) - 0.1056583755; // substarcting muon mass as reconstructed energy is KE and true energy is total energy
});

const RecoPartCut PartOverlapCut(double overlapThreshold){
  const RecoPartCut kPartOverlapCut([overlapThreshold](const caf::SRRecoParticleProxy * part)
  {
    if (part->truthOverlap.empty()) return false;
    return *std::max_element(part->truthOverlap.begin(), part->truthOverlap.end()) > overlapThreshold;
  });
  return kPartOverlapCut;
}

const RecoPartVar kPartStartX = SIMPLEPARTVAR(start.x);
const RecoPartVar kPartStartY = SIMPLEPARTVAR(start.y);
const RecoPartVar kPartStartZ = SIMPLEPARTVAR(start.z);
const RecoPartCut kPartNDLArContainedFV_StartAndEnd = kPartEndX < NDLArXHi - 25 && kPartEndX > NDLArXLo + 25
                                                   && kPartEndY < NDLArYHi - 25 && kPartEndY > NDLArYLo + 25
                                                   && kPartEndZ < NDLArZHi - 25 && kPartEndZ > NDLArZLo + 25
                                                   && kPartStartX < NDLArXHi - 25 && kPartStartX > NDLArXLo + 25
                                                   && kPartStartY < NDLArYHi - 25 && kPartStartY > NDLArYLo + 25
                                                   && kPartStartZ < NDLArZHi - 25 && kPartStartZ > NDLArZLo + 25;  

const RecoPartVar kMatchedTruePartAngle_TrueMuonCut([](const caf::SRRecoParticleProxy * part) -> float
{
    if (part->truth.empty()) return 0.0;
    const caf::SRProxy * sr = part->Ancestor<caf::SRProxy>();
    bool contained = kPartNDLArContained(part);
    size_t tidx = 0;
    float maxOverlap = 0;
    for (size_t i = 0; i < part->truthOverlap.size(); i++){
        if (part->truthOverlap[i] > maxOverlap){
            tidx = i;
            maxOverlap = part->truthOverlap[i];
        }
    }
    if (maxOverlap < 0.0 || kTruePDG(caf::FindParticle(sr->mc, part->truth[tidx])) != 13) return -1000.0;
    return kTruePartAngle(caf::FindParticle(sr->mc, part->truth[tidx]));
});

// Absolute residuals in cm
const Var kVtxResX = kVtxX - kMatchedTrueVtxX;
const Var kVtxResY = kVtxY - kMatchedTrueVtxY;
const Var kVtxResZ = kVtxZ - kMatchedTrueVtxZ;
// 3D displacement
const Var kVtxRes3D([](const caf::SRInteractionProxy *ixn) -> double {
    double dx = kVtxX(ixn) - kMatchedTrueVtxX(ixn);
    double dy = kVtxY(ixn) - kMatchedTrueVtxY(ixn);
    double dz = kVtxZ(ixn) - kMatchedTrueVtxZ(ixn);
    if (dx < -999 || dy < -999 || dz < -999) return -9999.;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
});
