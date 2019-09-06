#include "stdafx.h"
#include "ctmsnet.h"
#include <iostream>
#include <iomanip>

PCtmsNet TCtmsNet::CopyNet() const {
	PCtmsNet Net = TCtmsNet::New();
	for (TSignNet::TNodeI NI = BegNI(); NI < EndNI(); NI++) { Net->AddNode(NI); }
	for (TSignNet::TEdgeI EI = BegEI(); EI < EndEI(); EI++) { Net->AddEdge(EI); }
	if (*Net == TCtmsNet(*this)) { return Net; }
	return NULL;
}

// Returns a new sub-network with equal number of positive and negative edges chosen randomly.
PCtmsNet TCtmsNet::GetBalSgnProbSubNet() {
	TVec<TEdgeI> PosEVec, NegEVec;
	PCtmsNet NewNet = TCtmsNet::New();

	for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) {
		if (EI() == 1) { PosEVec.Add(EI); }
		else { NegEVec.Add(EI); }
	}
	PosEVec.Shuffle(TInt::Rnd);
	TInt::Rnd.Randomize();
	NegEVec.Shuffle(TInt::Rnd);
	for (int e = 0; e<NegEVec.Len() && e<PosEVec.Len(); e++) {
		if (!NewNet->IsNode(PosEVec[e].GetSrcNId())) { NewNet->AddNode(PosEVec[e].GetSrcNId()); }
		if (!NewNet->IsNode(PosEVec[e].GetDstNId())) { NewNet->AddNode(PosEVec[e].GetDstNId()); }
		if (!NewNet->IsNode(NegEVec[e].GetSrcNId())) { NewNet->AddNode(NegEVec[e].GetSrcNId()); }
		if (!NewNet->IsNode(NegEVec[e].GetDstNId())) { NewNet->AddNode(NegEVec[e].GetDstNId()); }
		IAssert(!NewNet->IsEdge(PosEVec[e].GetSrcNId(), PosEVec[e].GetDstNId()));
		IAssert(!NewNet->IsEdge(NegEVec[e].GetSrcNId(), NegEVec[e].GetDstNId()));
		NewNet->AddEdge(PosEVec[e]);
		NewNet->AddEdge(NegEVec[e]);
	}
	return NewNet;
}

// Returns a new sub-network with minimum Embeddedness of MinEmValue
PCtmsNet TCtmsNet::GetMinEmbeddedSubNet(const int MinEmValue) {
	PCtmsNet thisNet = PCtmsNet(this);
	PCtmsNet NewNet = TCtmsNet::New();
	int addE = 0, totE = 0;
	for (TEdgeI EI = thisNet->BegEI(); EI < thisNet->EndEI(); EI++) {
		TIntV NbV;
		totE++;
		TSnap::GetCmnNbrs(thisNet, EI.GetSrcNId(), EI.GetDstNId(), NbV);
		if (NbV.Len() < MinEmValue) { continue; }
		if (!NewNet->IsNode(EI.GetSrcNId())) { NewNet->AddNode(EI.GetSrcNId()); }
		if (!NewNet->IsNode(EI.GetDstNId())) { NewNet->AddNode(EI.GetDstNId()); }
		NewNet->AddEdge(EI);
		addE++;
	}
	int delE = 0, i = 0;
	bool chg;
	do {
		chg = false; i++;
		for (TSignNet::TEdgeI EI = NewNet->BegEI(); EI < NewNet->EndEI(); EI++) {
			TIntV NbV;
			TSnap::GetCmnNbrs(NewNet, EI.GetSrcNId(), EI.GetDstNId(), NbV);
			if (NbV.Len() >= MinEmValue) { continue; }
			NewNet->DelEdge(EI.GetSrcNId(), EI.GetDstNId());
			delE++;
			chg = true;
		}
	} while (chg == true);
	printf("\nTotal Edges:%d Deleted:%d Left:%d\n", totE, (totE - (addE - delE)), (addE - delE));
	return NewNet;
}

bool TCtmsNet::operator==(const TCtmsNet& SgnNet) const {
	if (Empty() && SgnNet.Empty()) { return true; } //No Node exists
	if (GetNodes() != SgnNet.GetNodes() || GetEdges() != SgnNet.GetEdges()) { return false; }
	for (TNodeI NI = SgnNet.BegNI(); NI < SgnNet.EndNI(); NI++) {
		if (!IsNode(NI.GetId())) { return false; }
		if (NI.GetDat() != GetNDat(NI.GetId())) { return false; }
		if (NI.GetDeg() != GetNode(NI.GetId()).GetDeg()) { return false; }
	}
	for (TEdgeI EI = SgnNet.BegEI(); EI < SgnNet.EndEI(); EI++) {
		if (!IsEdge(EI.GetSrcNId(), EI.GetDstNId())) { return false; }
		if (EI.GetDat() != GetEDat(EI.GetSrcNId(), EI.GetDstNId())) { return false; }
	}
	return true;
}

bool TCtmsNet::operator!=(const TCtmsNet& SgnNet) const {
	return !(*this == SgnNet);
}

// Modified Load function of Snap signnet.
PCtmsNet TCtmsNet::LoadEpinionsModified(const TStr& FNm) {
	TSsParser Ss(FNm, ssfWhiteSep, true, true, true);
	PCtmsNet Net = TCtmsNet::New();
	while (Ss.Next()) {
		const int src = Ss.GetInt(0);
		const int dst = Ss.GetInt(1);
		const int sgn = Ss.GetInt(2);
		if (src == dst) { continue; } // skip self edges
		if (!Net->IsNode(src)) {
			Net->AddNode(src);
		}
		if (!Net->IsNode(dst)) {
			Net->AddNode(dst);
		}
		Net->AddEdge(src, dst, sgn);
	}
	return Net;
}

PCtmsNet TCtmsNet::LoadSignedNet(const TStr& InFNm, const int& SrcColId, const int& DstColId, const int& SignColId) {
	TSsParser Ss(InFNm, ssfWhiteSep, true, true, true);
	PCtmsNet Graph = TCtmsNet::New();
	int SrcNId, DstNId, Sgn;
	while (Ss.Next()) {
		if (!Ss.GetInt(SrcColId, SrcNId) || !Ss.GetInt(DstColId, DstNId) || !Ss.GetInt(SignColId, Sgn)) { continue; }
		if (!Graph->IsNode(SrcNId)) { Graph->AddNode(SrcNId); }
		if (!Graph->IsNode(DstNId)) { Graph->AddNode(DstNId); }
		Graph->AddEdge(SrcNId, DstNId, Sgn);
	}
	Graph->Defrag();
	return Graph;
}

//Overloading of SNAP native GenRewire of PNGraph for TSignNet
PCtmsNet TCtmsNet::GenRewire(const int& NSwitch, TRnd Rnd) const {
	const int Nodes = GetNodes();
	const int Edges = GetEdges();
	PCtmsNet NetPt = TCtmsNet::New();
	TCtmsNet& Network = *NetPt;
	Network.Reserve(Nodes, Edges);
	TExeTm ExeTm;
	// generate a Network that satisfies the constraints
	printf("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("Randomizing edges (%d, %d)...\n", Nodes, Edges);
	//Abtin: Instead of TIntPairSet we use TIntTripleSet for Edge and its data. //TIntPrSet EdgeSet(Edges);
	typedef THashSet<TIntTr> TIntTrSet;
	TIntTrSet EdgeSet(Edges);
	THashSet<TInt> EdgeDatSet; //This can hold multiple signs (instead of only + and -) so we can rewire multiple signed subnets
	for (TNodeI NI = BegNI(); NI < EndNI(); NI++) {
		const int NId = NI.GetId();
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			EdgeSet.AddKey(TIntTr(NId, NI.GetOutNId(e), NI.GetOutEDat(e)));
			if (!EdgeDatSet.IsKey(NI.GetOutEDat(e))) { EdgeDatSet.AddKey(NI.GetOutEDat(e)); }
		}
		Network.AddNode(NI.GetId(), NI());
	}
	// edge switching
	uint skip = 0;
	//Abtin
	uint succsSwps = 0, skip2 = 0, skip3 = 0;
	for (uint swps = 0; swps < 2 * uint(Edges)*uint(NSwitch); swps++) {
		const int keyId1 = EdgeSet.GetRndKeyId(Rnd);
		//Abtin
		Rnd.Randomize();
		const int keyId2 = EdgeSet.GetRndKeyId(Rnd);
		if (keyId1 == keyId2) { skip++; continue; }
		const TIntTr& E1 = EdgeSet[keyId1];
		const TIntTr& E2 = EdgeSet[keyId2];
		//Abtin added next line. Don't swap + and - edges
		if (E1.Val3 != E2.Val3) { skip++; skip2++; continue; }
		//Abtin (a,b,sign) & (c,d,sign) --(swap)--> (a,d,sign) & (c,b,sign)
		TIntTr NewE1(E1.Val1, E2.Val2, E1.Val3), NewE2(E2.Val1, E1.Val2, E2.Val3);
		//Abtin: Looks if there is a similar Edge in the Network with any Edge Data (Sign). Generalizing only +/- form
		bool IsAlreadyEdge = false;
		for (int EDSetKeyId = 0; EDSetKeyId < EdgeDatSet.Len(); EDSetKeyId++) {
			if (EdgeSet.IsKey(TIntTr(E1.Val1, E2.Val2, EdgeDatSet[EDSetKeyId]))) { IsAlreadyEdge = true; break; } //NewE1 with all possible Edge Data
			if (EdgeSet.IsKey(TIntTr(E2.Val1, E1.Val2, EdgeDatSet[EDSetKeyId]))) { IsAlreadyEdge = true; break; } //NewE2 with all possible Edge Data
		}
		if (IsAlreadyEdge) { skip++; skip3++; continue; }
		/* Replaced following code with above to support multiple (more than 2) signs
		//Abtin added next 4 line. If new edges is in the hash list so we can't swap.
		if (EdgeSet.IsKey(TIntTr(E1.Val1, E2.Val2, E1.Val3))) { skip++; skip3++; continue; } //NewE1
		if (EdgeSet.IsKey(TIntTr(E1.Val1, E2.Val2, -E1.Val3))) { skip++; skip3++; continue; } //NewE1 converted sign
		if (EdgeSet.IsKey(TIntTr(E2.Val1, E1.Val2, E2.Val3))) { skip++; skip3++; continue; } //NewE2
		if (EdgeSet.IsKey(TIntTr(E2.Val1, E1.Val2, -E2.Val3))) { skip++; skip3++; continue; } //NewE2 converted sign
		*/
		//Abtin added the next line.
		if (NewE1.Val1 == NewE1.Val2 || NewE2.Val1 == NewE2.Val2) { skip++; }
		//Abtin changed the next line --> PREVIOUS: NewE1.Val1!=NewE2.Val1 && NewE1.Val2!=NewE2.Val1 && NewE1.Val2!=NewE2.Val1 && NewE1.Val2!=NewE2.Val2
		else if (NewE1.Val1 != NewE2.Val1 && NewE1.Val1 != NewE2.Val2 && NewE1.Val2 != NewE2.Val1 && NewE1.Val2 != NewE2.Val2 && !EdgeSet.IsKey(NewE1) && !EdgeSet.IsKey(NewE2)) {
			EdgeSet.DelKeyId(keyId1);  EdgeSet.DelKeyId(keyId2);
			EdgeSet.AddKey(TIntTr(NewE1));
			EdgeSet.AddKey(TIntTr(NewE2));
			succsSwps++;
		}
		else { skip++; }
		if (swps % Edges == 0) {
			printf("\r %uk/%uk: %uk skip %uk collision %uk InHash [%s]", swps / 1000u, 2 * uint(Edges)*uint(NSwitch) / 1000u, skip / 1000u, skip2 / 1000u, skip3 / 1000u, ExeTm.GetStr());
			if (ExeTm.GetSecs() > 2 * 3600) { printf(" *** Time limit!\n"); break; } // time limit 2 hours
		}
	}
	printf("\n  %uk / %uk was Successful.\n  %uk Skiped.\n  %uk Sign Collision.\n  %uk was already in Hash\n  Total Time [%s]\n\n",
		succsSwps / 1000u, 2 * uint(Edges)*uint(NSwitch) / 1000u, skip / 1000u, skip2 / 1000u, skip3 / 1000u, ExeTm.GetStr());
	for (int e = 0; e < EdgeSet.Len(); e++) {
		Network.AddEdge(EdgeSet[e].Val1, EdgeSet[e].Val2, EdgeSet[e].Val3);
	}
	printf("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	return NetPt;
}

//Rewriting the GenRewire. This function Does the GenRewire in different hashes for positive and negative edges (Positive Hash + Negative Hash)
typedef THashSet<TIntTr> TIntTrSet;
void SwitchEdges(TIntTrSet& EdgeSet, TIntTrSet& EdgeSetOther, const int& NSwitch, TRnd& Rnd) {
	// edge switching
	printf("\n********************\n");
	printf("Randomizing Subnet edges (%d edge)...\n", EdgeSet.Len());
	TExeTm ExeTm;
	uint skip = 0;
	//Abtin
	uint succsSwps = 0, skip3 = 0; //skip2 was for Sign collision which is not meaningful here
	for (uint swps = 0; swps < 2 * uint(EdgeSet.Len())*uint(NSwitch); swps++) {
		const int keyId1 = EdgeSet.GetRndKeyId(Rnd);
		const int keyId2 = EdgeSet.GetRndKeyId(Rnd);
		if (keyId1 == keyId2) { skip++; continue; }
		const TIntTr& E1 = EdgeSet[keyId1];
		const TIntTr& E2 = EdgeSet[keyId2];
		//Abtin (a,b,sign) & (c,d,sign) --(swap)--> (a,d,sign) & (c,b,sign)
		TIntTr NewE1(E1.Val1, E2.Val2, E1.Val3), NewE2(E2.Val1, E1.Val2, E2.Val3);
		//Abtin added next 4 line. If new edges is in the hash list so we can't swap.
		if (EdgeSet.IsKey(TIntTr(E1.Val1, E2.Val2, E1.Val3))) { skip++; skip3++; continue; } //NewE1
		if (EdgeSet.IsKey(TIntTr(E2.Val1, E1.Val2, E2.Val3))) { skip++; skip3++; continue; } //NewE2
		if (EdgeSetOther.IsKey(TIntTr(E1.Val1, E2.Val2, -E1.Val3))) { skip++; skip3++; continue; } //NewE1 converted sign
		if (EdgeSetOther.IsKey(TIntTr(E2.Val1, E1.Val2, -E2.Val3))) { skip++; skip3++; continue; } //NewE2 converted sign
																								   //Abtin added the next line.
		if (NewE1.Val1 == NewE1.Val2 || NewE2.Val1 == NewE2.Val2) { skip++; }
		//Abtin changed the next line --> PREVIOUS: NewE1.Val1!=NewE2.Val1 && NewE1.Val2!=NewE2.Val1 && NewE1.Val2!=NewE2.Val1 && NewE1.Val2!=NewE2.Val2
		else if (NewE1.Val1 != NewE2.Val1 && NewE1.Val1 != NewE2.Val2 && NewE1.Val2 != NewE2.Val1 && NewE1.Val2 != NewE2.Val2 && !EdgeSet.IsKey(NewE1) && !EdgeSet.IsKey(NewE2)) {
			EdgeSet.DelKeyId(keyId1);  EdgeSet.DelKeyId(keyId2);
			EdgeSet.AddKey(TIntTr(NewE1));
			EdgeSet.AddKey(TIntTr(NewE2));
			succsSwps++;
		}
		else { skip++; }
		if (swps % EdgeSet.Len() == 0) {
			printf("\r %uk/%uk: %uk skip %uk InHash [%s]", swps / 1000u, 2 * uint(EdgeSet.Len())*uint(NSwitch) / 1000u, skip / 1000u, skip3 / 1000u, ExeTm.GetStr());
			if (ExeTm.GetSecs() > 2 * 3600) { printf(" *** Time limit!\n"); break; } // time limit 2 hours
		}
	}
	printf("\n  %uk / %uk was Successful.\n  %uk Skiped.\n  %uk was already in Hash\n  Total Time [%s]\n",
		succsSwps / 1000u, 2 * uint(EdgeSet.Len())*uint(NSwitch) / 1000u, skip / 1000u, skip3 / 1000u, ExeTm.GetStr());
}
PCtmsNet TCtmsNet::GenRewire2Hash(const int& NSwitch, TRnd Rnd) const {
	const int Nodes = GetNodes();
	const int Edges = GetEdges();
	PCtmsNet NetPt = TCtmsNet::New();
	TCtmsNet& Network = *NetPt;
	Network.Reserve(Nodes, Edges);

	printf("\n***************************************************************\n");
	printf("Randomizing edges (%d, %d)...\n", Nodes, Edges);
	typedef THashSet<TIntTr> TIntTrSet;
	TIntTrSet PosEdgeSet(Edges), NegEdgeSet(Edges);
	for (TNodeI NI = BegNI(); NI < EndNI(); NI++) {
		const int NId = NI.GetId();
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			if (NI.GetOutEDat(e) == 1) { PosEdgeSet.AddKey(TIntTr(NId, NI.GetOutNId(e), NI.GetOutEDat(e))); }
			if (NI.GetOutEDat(e) == -1) { NegEdgeSet.AddKey(TIntTr(NId, NI.GetOutNId(e), NI.GetOutEDat(e))); }
		}
		Network.AddNode(NI.GetId(), NI());
	}

	SwitchEdges(PosEdgeSet, NegEdgeSet, NSwitch, Rnd);
	SwitchEdges(NegEdgeSet, PosEdgeSet, NSwitch, Rnd);

	// rebuilding the network using two hashes
	for (int e = 0; e < PosEdgeSet.Len(); e++) {
		Network.AddEdge(PosEdgeSet[e].Val1, PosEdgeSet[e].Val2, PosEdgeSet[e].Val3);
	}
	for (int e = 0; e < NegEdgeSet.Len(); e++) {
		Network.AddEdge(NegEdgeSet[e].Val1, NegEdgeSet[e].Val2, NegEdgeSet[e].Val3);
	}
	printf("\n***************************************************************\n");
	return NetPt;
}

PCtmsNet TCtmsNet::RewireNetwork2(const int &NSwitch, const bool separateHash) {
	return (separateHash ? GenRewire2Hash(NSwitch) : GenRewire(NSwitch));
}

//Overloaded signnet RewireNetwork() which returns a new Network instead of changing the original one.
PCtmsNet TCtmsNet::RewireNetwork(const int &NSwitch) {
	PCtmsNet Net = GetEdgeSubNet(1, -1);
	TIntH NIdDatH;
	for (TNodeI NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
		NIdDatH.AddDat(NI.GetId(), NI());
	}
	// rewire plus and minus networks
	PNGraph PlusG = TSnap::ConvertGraph<PNGraph>(Net->GetSignSubNet(+1));
	PNGraph MinusG = TSnap::ConvertGraph<PNGraph>(Net->GetSignSubNet(-1));
	PlusG = TSnap::GenRewire(PlusG, NSwitch);
	MinusG = TSnap::GenRewire(MinusG, NSwitch);
	// create network
	Net->Clr(false);
	for (TNGraph::TNodeI NI = PlusG->BegNI(); NI < PlusG->EndNI(); NI++) {
		Net->AddNode(NI.GetId()), NIdDatH.GetDat(NI.GetId());
	}
	for (TNGraph::TNodeI NI = MinusG->BegNI(); NI < MinusG->EndNI(); NI++) {
		if (!Net->IsNode(NI.GetId())) {
			Net->AddNode(NI.GetId(), NIdDatH.GetDat(NI.GetId()));
		}
	}
	for (TNGraph::TEdgeI EI = PlusG->BegEI(); EI < PlusG->EndEI(); EI++) {
		Net->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), +1);
	}
	for (TNGraph::TEdgeI EI = MinusG->BegEI(); EI < MinusG->EndEI(); EI++) {
		if (Net->IsEdge(EI.GetSrcNId(), EI.GetDstNId())) {
			Net->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), 0);
		}
		else Net->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), -1);
	}
	return Net;
}
//Modified SNAP signnet PermuteEdgeSigns() which returns a new Network instead of changing the original one.
PCtmsNet TCtmsNet::PermuteEdgeSigns2() {
	PCtmsNet Net = TCtmsNet::New();
	TIntV EDatV(GetEdges(), 0);
	for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) {
		EDatV.Add(EI());
		const int SrcNId = EI.GetSrcNId();
		const int DstNId = EI.GetDstNId();
		if (!Net->IsNode(SrcNId)) { Net->AddNode(SrcNId); }
		if (!Net->IsNode(DstNId)) { Net->AddNode(DstNId); }
		Net->AddEdge(EI);
	}
	TInt::Rnd.Randomize();
	EDatV.Shuffle(TInt::Rnd);
	int i = 0;
	for (TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) {
		EI() = EDatV[i++];
	}
	return Net;
}
//This function is dummy.. can be deleted
int TCtmsNet::GetStrictSgnPermSets(PCtmsNet& Net, THash<TIntPr, TIntVPr>& PermutableGroups) {
	TVec<TIntPrPr> PermutableSets;
	PCtmsNet StandardNet = Net->CopyNet();

	for (TSignNet::TEdgeI EI = StandardNet->BegEI(); EI < StandardNet->EndEI(); EI++) {
		if (EI.GetSrcNId() == EI.GetDstNId()) { StandardNet->DelEdge(EI.GetSrcNId(), EI.GetDstNId()); }
		if (EI() != 1 && EI() != -1) { printf("!"); StandardNet->DelEdge(EI.GetSrcNId(), EI.GetDstNId()); }
	}
	PCtmsNet PosSubNet = StandardNet->GetSignSubNet(1), NegSubNet = StandardNet->GetSignSubNet(-1);

	// All nodes must be present in both positive and negative networks
	PosSubNet->SortNIdById(); NegSubNet->SortNIdById();
	for (TSignNet::TNodeI NI = PosSubNet->BegNI(); NI < PosSubNet->EndNI(); NI++) {
		if (!NegSubNet->IsNode(NI.GetId())) { PosSubNet->DelNode(NI.GetId()); }
	}
	for (TSignNet::TNodeI NI = NegSubNet->BegNI(); NI < NegSubNet->EndNI(); NI++) {
		if (!PosSubNet->IsNode(NI.GetId())) { NegSubNet->DelNode(NI.GetId()); }
	}

	printf("Extracting sign-permutable sets\n");
	int c = 0, Decile = int(NegSubNet->GetNodes() / 100) + 1;
	for (TSignNet::TNodeI NI = NegSubNet->BegNI(); NI < NegSubNet->EndNI(); NI++) {
		const int BaseNdId = NI.GetId();
		THash<TInt, TVec<TInt>> Dist2NbrsNegPos;
		for (int eon = 0; eon < NI.GetOutDeg(); eon++) {
			const int MiddleNdId = NI.GetOutNId(eon);
			TSignNet::TNode MiddleNd = PosSubNet->GetNode(MiddleNdId); // Switch to Positive SubNet
			for (int eip = 0; eip < MiddleNd.GetInDeg(); eip++) {
				const int DestNdId = MiddleNd.GetInNId(eip);
				if (!Dist2NbrsNegPos.IsKey(DestNdId)) { Dist2NbrsNegPos.AddKey(DestNdId); }
				Dist2NbrsNegPos(DestNdId).Add(MiddleNdId);
			}
		}
		Dist2NbrsNegPos.SortByKey(false);
		const TSignNet::TNode BaseNodeP = PosSubNet->GetNode(BaseNdId);
		for (int eop = 0; eop < BaseNodeP.GetOutDeg(); eop++) {
			const int MiddleNdId = BaseNodeP.GetOutNId(eop);
			TSignNet::TNode MiddleNd = NegSubNet->GetNode(MiddleNdId); // Switch to Negative SubNet
			for (int i = 0; i < Dist2NbrsNegPos.Len(); i++) {
				const int DestNdId = Dist2NbrsNegPos.GetKey(i);
				if (BaseNdId > DestNdId) { break; } //This was inserted before
				if (!MiddleNd.IsInNId(DestNdId)) { continue; }
				IAssertR(BaseNdId != DestNdId, "BaseNode == DestNode");
				TIntPr SetKey(BaseNdId, DestNdId);
				if (!PermutableGroups.IsKey(SetKey)) {
					PermutableGroups.AddKey(SetKey);
					PermutableGroups(SetKey).Val1 = Dist2NbrsNegPos[i];
				}
				PermutableGroups(SetKey).Val2.Add(MiddleNdId);
			}
		}
		if (++c % Decile == 0) { printf("\r%d%%", c / Decile); }
	}
	printf("\n");
	for (int j = 0; j < PermutableGroups.Len(); j++) {
		const int baseNd = PermutableGroups.GetKey(j).Val1;
		const int destNd = PermutableGroups.GetKey(j).Val2;
		for (int s1 = 0; s1 < PermutableGroups[j].Val1.Len(); s1++) {
			for (int s2 = 0; s2 < PermutableGroups[j].Val2.Len(); s2++) {
				const int mdNd1 = PermutableGroups[j].Val1[s1];
				const int mdNd2 = PermutableGroups[j].Val2[s2];
				PermutableSets.Add(TIntPrPr(TIntPr(baseNd, destNd), TIntPr(mdNd1, mdNd2)));
				//printf("\n%s", PermutableSets[j].GetStr().CStr());
			}
		}
	}
	return PermutableSets.Len();
}
//This function is dummy.. can be deleted
PCtmsNet TCtmsNet::PermuteEdgeSignsStrict2() {
	PCtmsNet Net = CopyNet();
	int SuccCnt = 0, fail = 0;
	THash<TIntPr, TIntVPr> PermGrps;
	int sz = GetStrictSgnPermSets(Net, PermGrps);
	printf("Permutable-set Size: %d\n", sz);
	TInt::Rnd.Randomize();
	for (int i = 0; i < PermGrps.Len(); i++) {
		const int bNdId = PermGrps.GetKey(i).Val1;
		const int dNdId = PermGrps.GetKey(i).Val2;
		IAssert(bNdId != dNdId);
		for (int v1c = 0; v1c < PermGrps[i].Val1.Len(); v1c++) {
			for (int v2c = 0; v2c < PermGrps[i].Val2.Len(); v2c++) {
				if (TInt::Rnd.GetUniDev() < 0.5) { continue; }
				const int mNdId1 = PermGrps[i].Val1[v1c];
				const int mNdId2 = PermGrps[i].Val2[v2c];
				TSignNet::TEdgeI EOut1, EOut2, EIn1, EIn2;
				EOut1 = Net->GetEI(bNdId, mNdId1);
				EOut2 = Net->GetEI(bNdId, mNdId2);
				EIn1 = Net->GetEI(dNdId, mNdId1);
				EIn2 = Net->GetEI(dNdId, mNdId2);
				// Reversign the sign
				if (EOut1() != EOut2() && EIn1() != EIn2() && EOut1() != EIn1()) {
					EOut1() = EOut1() * (-1);
					EOut2() = EOut2() * (-1);
					EIn1() = EIn1() * (-1);
					EIn2() = EIn2() * (-1);
					SuccCnt++;
				}
				else { fail++; }
			}
		}
	}
	printf("Successful Permutations: %d\n", SuccCnt);
	printf("Failed Permutations: %d\n", fail);
	return Net;
}

int TCtmsNet::GetAllEdgeSgnPermSets(PCtmsNet& Net, TSignNet::TEdgeI& EI, TVec<TIntPrPr>& ESgnPermSets) {
	const int ESign = EI();
	TSignNet::TNode SrcBaseNd = Net->GetNode(EI.GetSrcNId()), DstBaseNd = Net->GetNode(EI.GetDstNId());
	//Looking for these sets: {A->C, A+>D, B+>C, B->D} which A->C or A+>D is EI
	for (int i = 0; i < SrcBaseNd.GetOutDeg(); i++) {
		const int DstMidNdId = SrcBaseNd.GetOutNId(i);
		if (Net->GetEI(SrcBaseNd.GetId(), DstMidNdId)() == ESign) { continue; }
		for (int j = 0; j < DstBaseNd.GetInDeg(); j++) {
			const int SrcMidNdId = DstBaseNd.GetInNId(j);
			if (Net->GetEI(SrcMidNdId, DstBaseNd.GetId())() == ESign) { continue; }
			if (!Net->IsEdge(SrcMidNdId, DstMidNdId)) { continue; }
			if (Net->GetEI(SrcMidNdId, DstMidNdId)() == (-1)*ESign) { continue; }
			//Adding New Sign Permutable Set
			ESgnPermSets.Add(TIntPrPr(TIntPr(SrcBaseNd.GetId(), SrcMidNdId), TIntPr(DstBaseNd.GetId(), DstMidNdId)));
		}
	}
	return ESgnPermSets.Len();
}

double TCtmsNet::GetEdgePosSgnProb(const int& SrcNdId, const int& DstNdId) {
	PCtmsNet thisNet = PCtmsNet((TCtmsNet *)this);
	TSignNet::TEdgeI Edge = GetEI(SrcNdId, DstNdId);
	const int sgn = Edge();
	TVec<TIntPrPr> EPermSets, temp;
	int RevCnt, OpRevCnt = 0;
	RevCnt = GetAllEdgeSgnPermSets(thisNet, Edge, EPermSets);
	if (RevCnt == 0) {
		if (sgn == 1) { return 1; }
		else { return 0; }
	}
	for (int rnd = 0; rnd < EPermSets.Len(); rnd++) {
		//const int rnd = TInt::Rnd.GetUniDevInt(0, EPermSets.Len()-1);
		TSignNet::TEdgeI PEdge = thisNet->GetEI(EPermSets[rnd].Val1.Val1, EPermSets[rnd].Val2.Val2);
		TSignNet::TEdgeI OEdge = thisNet->GetEI(EPermSets[rnd].Val1.Val2, EPermSets[rnd].Val2.Val1);
		TSignNet::TEdgeI QEdge = thisNet->GetEI(EPermSets[rnd].Val1.Val2, EPermSets[rnd].Val2.Val2);
		Edge() = (-1) * Edge();
		PEdge() = (-1) * PEdge();
		OEdge() = (-1) * OEdge();
		QEdge() = (-1) * QEdge();
		OpRevCnt += GetAllEdgeSgnPermSets(thisNet, Edge, temp);
		//Restoring the edge signs
		Edge() = (-1) * Edge();
		PEdge() = (-1) * PEdge();
		OEdge() = (-1) * OEdge();
		QEdge() = (-1) * QEdge();
	}
	OpRevCnt /= EPermSets.Len();
	if (sgn == -1) { return (double)RevCnt / (RevCnt + OpRevCnt); }
	else { return (double)OpRevCnt / (RevCnt + OpRevCnt); }
	/*
	TSignNet::TNode SrcBaseNd = thisNet->GetNode(SrcNdId), DstBaseNd = thisNet->GetNode(DstNdId);
	//Counting these sets: {(A,C)(A,D)(B,C)(B,D)}
	for (int i = 0; i < SrcBaseNd.GetOutDeg(); i++) {
	const int DstMidNdId = SrcBaseNd.GetOutNId(i);
	for (int j = 0; j < DstBaseNd.GetInDeg(); j++) {
	const int SrcMidNdId = DstBaseNd.GetInNId(j);
	if (!thisNet->IsEdge(SrcMidNdId, DstMidNdId)) {continue;}
	totCnt++;
	}
	}
	if (sgn == -1) {return double(sgnRevCnt)/totCnt;}//when edge is -1 and it gets reversed it is +1;
	else {return double(totCnt - sgnRevCnt)/totCnt;}//we want the probability of positive sign for the edge
	*/
}

bool TCtmsNet::GetRndESgnPermSet(PCtmsNet& Net, TSignNet::TEdgeI& EI, TIntPrPr& ERndPermSet) {
	const int ESign = EI();
	TSignNet::TNode SrcBaseNd = Net->GetNode(EI.GetSrcNId()), DstBaseNd = Net->GetNode(EI.GetDstNId());
	TVec<int> SrcNbrIdV, DstNbrIdV;
	//Creating vector of all valid neighbours
	for (int i = 0; i < SrcBaseNd.GetOutDeg(); i++) {
		const int DstMidNdId = SrcBaseNd.GetOutNId(i);
		if (Net->GetEI(SrcBaseNd.GetId(), DstMidNdId)() == ESign) { continue; }
		DstNbrIdV.Add(DstMidNdId);
	}
	for (int j = 0; j < DstBaseNd.GetInDeg(); j++) {
		const int SrcMidNdId = DstBaseNd.GetInNId(j);
		if (Net->GetEI(SrcMidNdId, DstBaseNd.GetId())() == ESign) { continue; }
		SrcNbrIdV.Add(SrcMidNdId);
	}
	// Shuffling vectors to get a random order of neighbours
	TRnd rnd1, rnd2;
	rnd1.Randomize(); rnd2.Randomize();
	SrcNbrIdV.Shuffle(rnd1);
	DstNbrIdV.Shuffle(rnd2);
	//Looking for these sets: {A->C, A+>D, B+>C, B->D} which A->C or A+>D is EI
	for (int s = 0; s < SrcNbrIdV.Len(); s++) {
		for (int d = 0; d < DstNbrIdV.Len(); d++) {
			if (!Net->IsEdge(SrcNbrIdV[s], DstNbrIdV[d])) { continue; }
			if (Net->GetEI(SrcNbrIdV[s], DstNbrIdV[d])() == (-1)*ESign) { continue; } //this edge must have the same sign as EI.
																					  //found one random set
			ERndPermSet = TIntPrPr(TIntPr(SrcBaseNd.GetId(), SrcNbrIdV[s]), TIntPr(DstBaseNd.GetId(), DstNbrIdV[d]));
			return true;
		}
	}
	return false;
}

PCtmsNet TCtmsNet::PermuteEdgeSignsStrict() {
	PCtmsNet Net = CopyNet();
	TVec<TSignNet::TEdgeI> EdgesV;
	// shuffling edge order to get a more accurate result
	for (TSignNet::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) { EdgesV.Add(EI); }
	TRnd rnd; rnd.Randomize();
	EdgesV.Shuffle(rnd);
	IAssert(EdgesV.Len() == Net->GetEdges());
	printf("Permuting network (network structure and node-sign-distribution is preserved)\n");
	int permCnt = 0;
	TInt::Rnd.Randomize();
	int c = 0, Decile = int(Net->GetEdges() / 100) + 1;
	for (int v = 0; v < EdgesV.Len(); v++) {
		if (c++ % Decile == 0) { printf("\r%d%%", c / Decile); }
		TSignNet::TEdgeI EI = EdgesV[v];
		IAssert(EI() == 1 || EI() == -1);
		if (TInt::Rnd.GetUniDev() < 0.5) { continue; } //Randomly decide to permute an edge or not		
		TIntPrPr EPSet;
		if (!GetRndESgnPermSet(Net, EI, EPSet)) { continue; }
		TSignNet::TEdgeI CurE, BsdE, OpE1, OpE2;
		CurE = Net->GetEI(EPSet.Val1.Val1, EPSet.Val2.Val1); IAssert(CurE == EI);
		BsdE = Net->GetEI(EPSet.Val1.Val1, EPSet.Val2.Val2);
		OpE1 = Net->GetEI(EPSet.Val1.Val2, EPSet.Val2.Val1);
		OpE2 = Net->GetEI(EPSet.Val1.Val2, EPSet.Val2.Val2);
		// Reversign Sign
		CurE() = (-1) * CurE();
		BsdE() = (-1) * BsdE();
		OpE1() = (-1) * OpE1();
		OpE2() = (-1) * OpE2();
		permCnt++;
	}
	printf("\n%d permutations done.\n", permCnt);
	return Net;
}

void TCtmsNet::CompareNetWith(const PCtmsNet& PNet2) const {
	using namespace std;
	const TCtmsNet& Net2 = *PNet2();
	cout << "\n------- Comparing Neworks -------\n" << endl;
	if (*this == Net2) { cout << "Networks are the SAME\n" << endl; }
	else { cout << "Networks are DIFFERENT\n" << endl; }

	cout << "This Network vs. Compared Network:  ";
	cout << "(" << GetNodes() << ", " << GetEdges() << ") vs. "
		<< "(" << Net2.GetNodes() << ", " << Net2.GetEdges() << ")" << endl;

	int NodeMissMatch = 0, NodeDataMissMatch = 0, InDegreeMissMatch = 0, OutDegreeMissMatch = 0, PosDMissMatch = 0, NegDMissMatch = 0;
	for (TNodeI NI = Net2.BegNI(); NI < Net2.EndNI(); NI++) {
		if (!IsNode(NI.GetId())) { NodeMissMatch++; continue; }
		if (NI.GetDat() != GetNDat(NI.GetId())) { NodeDataMissMatch++; }
		if (NI.GetInDeg() != GetNode(NI.GetId()).GetInDeg()) { InDegreeMissMatch++; }
		if (NI.GetOutDeg() != GetNode(NI.GetId()).GetOutDeg()) { OutDegreeMissMatch++; continue; }

		int thisNodePosDeg = 0, thisNodeNegDeg = 0, otherNodePosDeg = 0, otherNodeNegDeg = 0;
		for (int i = 0; i < NI.GetOutDeg(); i++) {
			(NI.GetOutEDat(i) == +1) ? otherNodePosDeg++ : otherNodeNegDeg++;
			(GetNode(NI.GetId()).GetOutEDat(i) == +1) ? thisNodePosDeg++ : thisNodeNegDeg++;
		}
		if (thisNodePosDeg != otherNodePosDeg) { PosDMissMatch++; }
		if (thisNodeNegDeg != otherNodeNegDeg) { NegDMissMatch++; }
	}

	int  EdgeMissMatch = 0, EdgeDataMissMatch = 0;
	for (TEdgeI EI = Net2.BegEI(); EI < Net2.EndEI(); EI++) {
		if (!IsEdge(EI.GetSrcNId(), EI.GetDstNId())) { EdgeMissMatch++; continue; }
		if (EI.GetDat() != GetEDat(EI.GetSrcNId(), EI.GetDstNId())) { EdgeDataMissMatch++; }
	}
	cout << "# Node Miss-Matches : " << GetNodes() - Net2.GetNodes() << endl;
	cout << "# Edge Miss-Matches : " << GetEdges() - Net2.GetEdges() << endl;
	cout << "Node ID Miss-Matches : " << NodeMissMatch << endl;
	cout << "Edge ID Miss-Matches : " << EdgeMissMatch << endl;
	cout << "Node Data Miss-Matches : " << NodeDataMissMatch << endl;
	cout << "Edge Data Miss-Matches : " << EdgeDataMissMatch << endl;
	cout << "Node In-Degree Miss-Matches: " << InDegreeMissMatch << endl;
	cout << "Node Out-Degree Miss-Matches: " << OutDegreeMissMatch << endl;
	cout << "Positive Degree Miss-Matches : " << PosDMissMatch << endl;
	cout << "Negative Degree Miss-Matches : " << NegDMissMatch << endl;
	cout << "\n------- ----------------- -------\n" << endl;
	return;
}

TVec<TIntTr> TCtmsNet::GetEmbDistrb() const {
	TVec<TIntTr> distriputionVec;
	THash<TInt, TIntPr> distribH;
	const PSignNet PThis = new TSignNet(*this);
	for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) {
		const TInt emb = TSnap::GetCmnNbrs(PThis, EI.GetSrcNId(), EI.GetDstNId());
		if (!distribH.IsKey(emb))
			distribH.AddDat(emb, TIntPr(0, 0));
		distribH.GetDat(emb).Val1++;
		if (EI() == 1)
			distribH.GetDat(emb).Val2++;
	}
	distribH.SortByKey(true);
	for (int i = 0; i < distribH.Len(); i++)
		distriputionVec.Add(TIntTr(distribH.GetKey(i), distribH[i].Val1, distribH[i].Val2));
	return distriputionVec;
}

bool TCtmsNet::IsSameTriad2(const PCtmsNet& Net1, const PCtmsNet& Net2) {
	static THash<TInt, TVec<TIntV> > PermH;
	if (PermH.Empty()) {
		PermH.AddDat(2).Add(TIntV::GetV(0, 1));
		PermH.AddDat(2).Add(TIntV::GetV(1, 0));
		PermH.AddDat(3).Add(TIntV::GetV(0, 1, 2));
		PermH.AddDat(3).Add(TIntV::GetV(0, 2, 1));
		PermH.AddDat(3).Add(TIntV::GetV(1, 0, 2));
		PermH.AddDat(3).Add(TIntV::GetV(1, 2, 0));
		PermH.AddDat(3).Add(TIntV::GetV(2, 0, 1));
		PermH.AddDat(3).Add(TIntV::GetV(2, 1, 0));
	}
	if (Net1->GetNodes() != Net2->GetNodes() || Net1->GetEdges() != Net2->GetEdges()) { return false; }
	IAssert(Net1->GetNodes() == 2 || Net1->GetNodes() == 3);

	const TVec<TIntV>& PermV = PermH.GetDat(Net1->GetNodes());
	for (int p = 0; p < PermV.Len(); p++) {
		const TIntV& Perm = PermV[p];
		PCtmsNet NewNet2 = TCtmsNet::New();
		for (int i = 0; i < Net2->GetNodes(); i++) {
			NewNet2->AddNode(i, Net2->GetNDat(Perm[i]));
		}

		for (TSignNet::TEdgeI EI2 = Net2->BegEI(); EI2 < Net2->EndEI(); EI2++) {
			NewNet2->AddEdge(Perm.SearchForw(EI2.GetSrcNId(), 0), Perm.SearchForw(EI2.GetDstNId(), 0), EI2.GetDat());
		}
		if (NewNet2 == Net1) { return true; }
		//NewNet2->Clr();
	}
	return false;
}

int TCtmsNet::EdgeSig(int e1, int e2) {
	if (e1 == 0 && e2 == 0) { return 0; }
	if (e1 == 1 && e2 == 0) { return 1; }
	if (e1 == -1 && e2 == 0) { return 2; }
	if (e1 == 0 && e2 == 1) { return 3; }
	if (e1 == 0 && e2 == -1) { return 4; }
	if (e1 == 1 && e2 == 1) { return 5; }
	if (e1 == 1 && e2 == -1) { return 6; }
	if (e1 == -1 && e2 == 1) { return 7; }
	if (e1 == -1 && e2 == -1) { return 8; }
	Fail; return -1;
}

TChA TCtmsNet::GetTriadStr(int a[2], int b[2], int c[2]) {
	const int e1 = EdgeSig(a[0], a[1]);
	const int e2 = EdgeSig(b[0], b[1]);
	const int e3 = EdgeSig(c[0], c[1]);
	const TIntTr MnTr1 = TMath::Mn(TIntTr(e1, e2, e3), TIntTr(e2, e3, e1), TIntTr(e3, e1, e2));
	const int e1a = EdgeSig(c[1], c[0]);
	const int e2a = EdgeSig(b[1], b[0]);
	const int e3a = EdgeSig(a[1], a[0]);
	const TIntTr MinTr = TMath::Mn(MnTr1, TMath::Mn(TIntTr(e1a, e2a, e3a), TIntTr(e2a, e3a, e1a), TIntTr(e3a, e1a, e2a)));

	const TInt EdgNum[] = { MinTr.Val1, MinTr.Val2, MinTr.Val3 };
	TChA EdgStr[3];
	for (int i = 0; i < 3; i++) {
		switch (EdgNum[i]) {
		case 1: EdgStr[i] = "+o"; break;
		case 2: EdgStr[i] = "-o"; break;
		case 3: EdgStr[i] = "o+"; break;
		case 4: EdgStr[i] = "o-"; break;
		case 5: EdgStr[i] = "++"; break;
		case 6: EdgStr[i] = "+-"; break;
		case 7: EdgStr[i] = "-+"; break;
		case 8: EdgStr[i] = "--"; break;
		default: EdgStr[i] = TChA();
		}
	}
	return EdgStr[0] += TChA("|") += EdgStr[1] += TChA("|") += EdgStr[2];
}

TChA TCtmsNet::GetTriadStr(const PCtmsNet& Nt, int srcId, int dstId, int nbrId, bool IsSigned) {
	int SDSgn[2] = { 0, 0 }, DNSgn[2] = { 0, 0 }, NSSgn[2] = { 0, 0 };
	if (IsSigned) {
		if (Nt->IsEdge(srcId, dstId)) {
			SDSgn[0] = Nt->GetEDat(srcId, dstId);
		}
		if (Nt->IsEdge(dstId, srcId)) {
			SDSgn[1] = Nt->GetEDat(dstId, srcId);
		}
		if (Nt->IsEdge(dstId, nbrId)) {
			DNSgn[0] = Nt->GetEDat(dstId, nbrId);
		}
		if (Nt->IsEdge(nbrId, dstId)) {
			DNSgn[1] = Nt->GetEDat(nbrId, dstId);
		}
		if (Nt->IsEdge(nbrId, srcId)) {
			NSSgn[0] = Nt->GetEDat(nbrId, srcId);
		}
		if (Nt->IsEdge(srcId, nbrId)) {
			NSSgn[1] = Nt->GetEDat(srcId, nbrId);
		}
	}
	else {
		if (Nt->IsEdge(srcId, dstId)) {
			SDSgn[0] = 1;
		}
		if (Nt->IsEdge(dstId, srcId)) {
			SDSgn[1] = 1;
		}
		if (Nt->IsEdge(dstId, nbrId)) {
			DNSgn[0] = 1;
		}
		if (Nt->IsEdge(nbrId, dstId)) {
			DNSgn[1] = 1;
		}
		if (Nt->IsEdge(nbrId, srcId)) {
			NSSgn[0] = 1;
		}
		if (Nt->IsEdge(srcId, nbrId)) {
			NSSgn[1] = 1;
		}
	}
	//return GetTriadStr(SDSgn, DNSgn, NSSgn);
	int *a = SDSgn, *b = DNSgn, *c = NSSgn;
	const int e1 = EdgeSig(a[0], a[1]);
	const int e2 = EdgeSig(b[0], b[1]);
	const int e3 = EdgeSig(c[0], c[1]);
	const TIntTr MnTr1 = TMath::Mn(TIntTr(e1, e2, e3), TIntTr(e2, e3, e1), TIntTr(e3, e1, e2));
	const int e1a = EdgeSig(c[1], c[0]);
	const int e2a = EdgeSig(b[1], b[0]);
	const int e3a = EdgeSig(a[1], a[0]);
	const TIntTr MinTr = TMath::Mn(MnTr1, TMath::Mn(TIntTr(e1a, e2a, e3a), TIntTr(e2a, e3a, e1a), TIntTr(e3a, e1a, e2a)));
	//return TCtmsNet::GetTriadStr(MinTr);
	const TInt EdgNum[] = { MinTr.Val1, MinTr.Val2, MinTr.Val3 };
	TChA EdgStr[3];
	for (int i = 0; i < 3; i++) {
		switch (EdgNum[i]) {
		case 1: EdgStr[i] = "+o"; break;
		case 2: EdgStr[i] = "-o"; break;
		case 3: EdgStr[i] = "o+"; break;
		case 4: EdgStr[i] = "o-"; break;
		case 5: EdgStr[i] = "++"; break;
		case 6: EdgStr[i] = "+-"; break;
		case 7: EdgStr[i] = "-+"; break;
		case 8: EdgStr[i] = "--"; break;
		default: EdgStr[i] = TChA();
		}
	}
	return EdgStr[0] += TChA("|") += EdgStr[1] += TChA("|") += EdgStr[2];
}

TChA TCtmsNet::GetTriadStr() const {
	TCtmsNet::TNodeI NI = this->BegNI();
	const int srcNdId = NI.GetId(); NI++;
	const int dstNdId = NI.GetId(); NI++;
	const int nbrNdId = NI.GetId();	
	return TCtmsNet::GetTriadStr(this->GetTriad(srcNdId, dstNdId, nbrNdId), srcNdId, dstNdId, nbrNdId);
}

void TCtmsNet::GenTriadEquivClasses(TTriadEqClasH& TriadGroups, const bool BiDirEdgeSide, const bool Signed, const bool ZeroEdgeSide) {
	PCtmsNet Triad;
	TVec<PCtmsNet> TriadIdV;
	TIntH TriadIdCntH;
	THash<TInt, int> NodeState;
	enum EdgeTypes { NoE = 0, FP, FN, BP, BN, FPBP, FNBN, FPBN, FNBP };
	NodeState.AddDat(0); NodeState.AddDat(1); NodeState.AddDat(2);

	TVec<int> EdgeState;
	if (Signed) {
		if (ZeroEdgeSide && BiDirEdgeSide) { EdgeState = TVec<int>::GetV(NoE, FP, FN, BP, BN, FPBP, FNBN, FPBN, FNBP); }
		else if (!ZeroEdgeSide && BiDirEdgeSide) { EdgeState = TVec<int>::GetV(FP, FN, BP, BN, FPBP, FNBN, FPBN, FNBP); }
		else if (ZeroEdgeSide && !BiDirEdgeSide) { EdgeState = TVec<int>::GetV(NoE, FP, FN, BP, BN); }
		else { EdgeState = TVec<int>::GetV(FP, FN, BP, BN); }
	}
	else {
		if (ZeroEdgeSide && BiDirEdgeSide) { EdgeState = TVec<int>::GetV(NoE, FP, BP, FPBP); }
		else if (!ZeroEdgeSide && BiDirEdgeSide) { EdgeState = TVec<int>::GetV(FP, BP, FPBP); }
		else if (ZeroEdgeSide && !BiDirEdgeSide) { EdgeState = TVec<int>::GetV(NoE, FP, BP); }
		else { EdgeState = TVec<int>::GetV(FP, BP); }
	}

	for (int i = 0; i < EdgeState.Len(); i++) {
		NodeState[0] = EdgeState[i];
		for (int j = 0; j < EdgeState.Len(); j++) {
			NodeState[1] = EdgeState[j];
			for (int k = 0; k < EdgeState.Len(); k++) {
				NodeState[2] = EdgeState[k];
				Triad = TCtmsNet::New();
				Triad->AddNode(NodeState.GetKey(0)); Triad->AddNode(NodeState.GetKey(1)); Triad->AddNode(NodeState.GetKey(2));
				for (TSignNet::TNodeI NI = Triad->BegNI(); NI < Triad->EndNI(); NI++) {
					TSignNet::TNodeI NxtNI(NI);	NxtNI++;
					if (NxtNI == Triad->EndNI()) { NxtNI = Triad->BegNI(); }
					switch (NodeState[NI.GetId()])
					{
					case FP://+o
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), 1);
						break;
					case FN://-o
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), -1);
						break;
					case BP://o+
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), 1);
						break;
					case BN://o-
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), -1);
						break;
					case FPBP://++
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), 1);
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), 1);
						break;
					case FNBN://--
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), -1);
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), -1);
						break;
					case FPBN://+-
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), 1);
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), -1);
						break;
					case FNBP://-+
						Triad->AddEdge(NI.GetId(), NxtNI.GetId(), -1);
						Triad->AddEdge(NxtNI.GetId(), NI.GetId(), 1);
						break;
					default:
						break;
					}
				}
				Triad->SortNIdById();
				// Count signed triad
				int TriadId = -1;
				for (int i = 0; i < TriadIdV.Len() && TriadId == -1; i++) {
					if (TCtmsNet::IsSameTriad(&*TriadIdV[i], &*Triad) != -1) { TriadId = i; break; }
				}
				if (TriadId == -1) {
					PCtmsNet NewTriad = Triad->CopyNet();
					TriadId = TriadIdV.Len();  TriadIdV.Add(NewTriad);
				}
				TriadIdCntH.AddDat(TriadId) += 1;
				// Deleting Triad
				Triad->Clr();
			}
		}
	}
	//Create (Class String Key, Triad Network, Number of Triads in the class)
	for (int i = 0; i < TriadIdCntH.Len(); i++) {
		TriadGroups.AddDat(TriadIdV[TriadIdCntH.GetKey(i)]->GetTriadStr(), TPair<PCtmsNet, int>(TriadIdV[TriadIdCntH.GetKey(i)], TriadIdCntH[i]));
	}
	return;
}

// Returns the probability of Triad in its Equivalency Class
double TCtmsNet::GetTriadProb2(const double& PlusProb) const {
	static TTriadEqClasH SignedTriadGroups, UnSignedTriadGroups;

	if (SignedTriadGroups.Empty()) { GenTriadEquivClasses(SignedTriadGroups); }
	if (UnSignedTriadGroups.Empty()) { GenTriadEquivClasses(UnSignedTriadGroups, true, false); }
	const int E = GetEdges();
	int P = 0;
	for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) {
		const int sign = EI();  IAssert(sign == 1 || sign == -1);
		if (sign == 1) { P++; }
	}
	if (GetNodes() == 3) {
		PCtmsNet Triad = CopyNet();
		Triad->SetAllEDat(1);
		const int UnSigEqCnt = UnSignedTriadGroups(Triad->GetTriadStr()).Val2;
		const int SigEqCnt = SignedTriadGroups(GetTriadStr()).Val2;
		return (SigEqCnt / UnSigEqCnt) * pow(PlusProb, P) * pow(1 - PlusProb, E - P);
	}
	return pow(PlusProb, P) * pow(1 - PlusProb, E - P);
}

void TCtmsNet::CountSignedTriads2(const TStr& OutFNm) const {
	printf("Counting signed triads");
	// permutations
	static TTriadEqClasH SigEqClasses, UnSigEqClasses;
	if (SigEqClasses.Empty()) { GenTriadEquivClasses(SigEqClasses); }
	if (UnSigEqClasses.Empty()) { GenTriadEquivClasses(UnSigEqClasses, true, false); }
	for (int i = 0; i < SigEqClasses.Len(); i++) { SigEqClasses[i].Val2 = 0; }
	for (int i = 0; i < UnSigEqClasses.Len(); i++) { UnSigEqClasses[i].Val2 = 0; }

	TIntV NbrV;
	PSignNet ThisPt = PSignNet((TSignNet*) this);
	double AllPlusE = 0, AllE = GetEdges();
	//Abtin added
	int AllTriadsCnt = 0;

	int c = 0, Decile = int(AllE / 100) + 1; //Abtin: +1 is to prevent docile get zero.
	for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) {
		TSnap::GetCmnNbrs(ThisPt, EI.GetSrcNId(), EI.GetDstNId(), NbrV);
		for (int n = 0; n < NbrV.Len(); n++) {
			PCtmsNet TriadNet = GetTriad(EI.GetSrcNId(), EI.GetDstNId(), NbrV[n]);
			// count signed triad
			SigEqClasses(TriadNet->GetTriadStr()).Val2 += 1;
			// count unsigned triads
			TriadNet->SetAllEDat(1);
			UnSigEqClasses(TriadNet->GetTriadStr()).Val2 += 1;
		}
		if (EI() == 1) { AllPlusE += 1; }
		if (++c % Decile == 0) { printf("."); }
	}
	const double PlusProb = AllPlusE / AllE;
	int SigTriadsClasPresent = 0, UnSigTriadsClasPresent = 0;
	//Abtin: Each Triad is counted according to the number of edges present in it
	for (int t = 0; t < SigEqClasses.Len(); t++) {
		const int E = SigEqClasses[t].Val1->GetEdges();
		SigEqClasses[t].Val2 /= E;
		AllTriadsCnt += SigEqClasses[t].Val2;
		if (SigEqClasses[t].Val2 != 0) { SigTriadsClasPresent++; }
	}
	//Abtin: Each Triad is counted according to the number of edges present in it
	for (int u = 0; u < UnSigEqClasses.Len(); u++) {
		const int E = UnSigEqClasses[u].Val1->GetEdges();
		UnSigEqClasses[u].Val2 /= E;
		if (UnSigEqClasses[u].Val2 != 0) { UnSigTriadsClasPresent++; }
		printf("\n%d. Unsigned triad: %s\t%d", u + 1, UnSigEqClasses[u].Val1->GetTriadStr().CStr(), UnSigEqClasses[u].Val2);
	}
	// draw
	FILE *T = fopen(TStr::Fmt("%s-SignTriad.tab", OutFNm.CStr()).CStr(), "wt");
	fprintf(T, "AB\tBC\tCA\tCount\tE[Count]\tSurprise\tTriadProb\n");
	printf("\n%d unsigned triads equivalent classes", UnSigTriadsClasPresent);
	printf("\n%d signed triads equivalent classes", SigTriadsClasPresent);
	printf("\nplus prob %d / %d = %f", int(AllPlusE), int(AllE), PlusProb);
	for (int t = 0; t < SigEqClasses.Len(); t++) {
		int PlusE = 0;
		PCtmsNet TriadNet = SigEqClasses[t].Val1->CopyNet();
		PCtmsNet TriadNetUnsig = TriadNet->CopyNet(); TriadNetUnsig->SetAllEDat(1);
		TStr FNm = TStr::Fmt("%s-SignTriad-%02d", OutFNm.CStr(), t + 1);
		FILE *F = fopen(TStr(FNm + ".dot").CStr(), "wt");
		fprintf(F, "digraph G {\n");
		fprintf(F, "  graph [splines=true, overlap=false]\n  node  [shape=ellipse, width=0.3, height=0.3 label=\"\"]\n");
		for (TEdgeI EI = TriadNet->BegEI(); EI < TriadNet->EndEI(); EI++) {
			fprintf(F, "  n%d -> n%d [label=\"%s\" len=2];\n", EI.GetSrcNId(), EI.GetDstNId(), EI() == 1 ? "+" : "-");
			if (EI() == 1) { PlusE++; }
		}
		const double TriadCnt = SigEqClasses[t].Val2;
		const double UnSignCnt = UnSigEqClasses(TriadNetUnsig->GetTriadStr()).Val2;
		const double TriadProbInEqCls = TriadNet->GetTriadProb2(PlusProb);
		const double TriadProb = TriadProbInEqCls * (UnSignCnt / AllTriadsCnt);
		// number of all isorphic nonsigned triads
		const double ExpCnt = TriadProbInEqCls * UnSignCnt; // is equal to: TriadProb * AllTriadsCnt
		const double Surp = (TriadCnt - ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0 - TriadProb));
		fprintf(F, "  label = \"T=%d, E[T]=%d, S=%.1f\";\n}\n", int(TriadCnt), int(ExpCnt), Surp);
		fclose(F);
		//TGraphViz::DoLayout(FNm+".dot", FNm+".gif",  gvlNeato);
		fprintf(T, "%s\t%d\t%.2f\t%.2f\t%f\n", TriadNet->GetTriadStr().CStr(), int(TriadCnt), ExpCnt, Surp, TriadProb);
	}
	fclose(T);
	printf("\nall triads Count: %d\n", AllTriadsCnt);
}

void TCtmsNet::CountSignedTriads3(const TStr& OutFNm) const {
	TTriadEqClasH NetSigTriDistrib;
	TTriadEqClasH PermSigTriDisrib, RewireSigTriDistrib;
	PCtmsNet ThisPt = PCtmsNet((TCtmsNet*) this), PermNet, RewiredNet;
	double NetAllTriads, PermAllTriads, RewAllTriads;
	PermNet = ThisPt->PermuteEdgeSigns2();
	RewiredNet = ThisPt->RewireNetwork2();

	if (NetSigTriDistrib.Empty()) { NetAllTriads = TCtmsNet::CalcNetTriadsDistrib(NetSigTriDistrib, ThisPt); }
	PermAllTriads = TCtmsNet::CalcNetTriadsDistrib(PermSigTriDisrib, PermNet);
	RewAllTriads = TCtmsNet::CalcNetTriadsDistrib(RewireSigTriDistrib, RewiredNet);

	FILE *T = fopen(TStr::Fmt("%s-TriadCountClasses.txt", OutFNm.CStr()).CStr(), "wt");
	fprintf(T, "AB\tBC\tCA\tCount\tProb\tPermCnt\tPermPr\tRewCnt\tRewPr\n");
	for (int i = 0; i < NetSigTriDistrib.Len(); i++)
	{
		TChA TriClassKey = NetSigTriDistrib[i].Val1->GetTriadStr();
		fprintf(T, "%s\t%d\t%f\t%d\t%f\t%d\t%f\n", TriClassKey.CStr(),
			NetSigTriDistrib(TriClassKey).Val2, NetSigTriDistrib(TriClassKey).Val2 / NetAllTriads,
			PermSigTriDisrib(TriClassKey).Val2, PermSigTriDisrib(TriClassKey).Val2 / PermAllTriads,
			RewireSigTriDistrib(TriClassKey).Val2, RewireSigTriDistrib(TriClassKey).Val2 / RewAllTriads);
	}
	fclose(T);
	return;
}

int TCtmsNet::CalcNetTriadsDistrib(TTriadEqClasH& NetSigTriDistrib, PCtmsNet& Net) {
	printf("Counting Networks' signed triads in each Equivalency Class\n");
	TIntV NbrV;
	const int AllE = Net->GetEdges();
	int c = 0, Decile = int(AllE / 100) + 1, AllTriads = 0;

	if (NetSigTriDistrib.Len() != 96) {
		NetSigTriDistrib.Clr();
		GenTriadEquivClasses(NetSigTriDistrib);
	}
	for (int i = 0; i < NetSigTriDistrib.Len(); i++) { NetSigTriDistrib[i].Val2 = 0; }

	//counting Net Triads in each class
	for (TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) {
		TSnap::GetCmnNbrs(Net, EI.GetSrcNId(), EI.GetDstNId(), NbrV);
		for (int n = 0; n < NbrV.Len(); n++) {
			PCtmsNet Tri = Net->GetTriad(EI.GetSrcNId(), EI.GetDstNId(), NbrV[n]);
			NetSigTriDistrib(Tri->GetTriadStr()).Val2 += 1;
		}
		if (++c % Decile == 0) { printf("."); }
	}
	for (int i = 0; i < NetSigTriDistrib.Len(); i++) {
		NetSigTriDistrib[i].Val2 /= NetSigTriDistrib[i].Val1->GetEdges();
		AllTriads += NetSigTriDistrib[i].Val2;
	}
	printf("\n> Counting Results:\t%d Edges\t%d Triads\n", AllE, AllTriads);
	return AllTriads;
}