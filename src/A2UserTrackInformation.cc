// A2UserTrackInformation
// Author: Dominik Werthmueller, 2017

#include "G4Track.hh"

#include "A2UserTrackInformation.hh"

//______________________________________________________________________________
A2UserTrackInformation::A2UserTrackInformation()
{
    // Constructor.
    //comment from Yannis in 2026, to clear up confussion: a better name would be something like fAncestorTrackID, 
    //because all descendants get the track id of their  related primary track. It is infact not a unique id,
    //to identify each track. The same applies for fPartID
    fTrackID = 0;
    fPartID = -1;
}

//______________________________________________________________________________
A2UserTrackInformation::A2UserTrackInformation(const G4Track* aTrack)
{
    // Constructor.

    fTrackID = aTrack->GetTrackID();
    fPartID = -1;
}

//______________________________________________________________________________
A2UserTrackInformation::A2UserTrackInformation(const A2UserTrackInformation* aTrackInfo)
{
    // Constructor.

    fTrackID = aTrackInfo->GetTrackID();
    fPartID = aTrackInfo->GetPartID();
}

//______________________________________________________________________________
void A2UserTrackInformation::Print() const
{
    // Print track information.

    G4cout << "A2UserTrackInformation: "
              "fTrackID: " << fTrackID << "  "
              "fPartID: "  << fPartID  << G4endl;
}

