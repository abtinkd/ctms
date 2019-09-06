#ifndef ctms_signed_network_h
#define ctms_signed_network_h

#include <Snap.h>
#include <mxdag.h>
#include <signnet.h>
#include <network.h>

class TCtmsNet;
typedef TPt<TCtmsNet> PCtmsNet;
typedef TPair<TIntPr, TIntPr> TIntPrPr;
typedef TPair<TIntV, TIntV> TIntVPr;
typedef THash<TChA, TPair<PCtmsNet, int>> TTriadEqClasH;

/////////////////////////////////////////////////
// Extended version of Signed network
class TCtmsNet : public TSignNet {
public:
	TCtmsNet() { }
	TCtmsNet(TSIn& SIn) : TSignNet(SIn) { }
	TCtmsNet(TSignNet& TSNet) : TSignNet(TSNet) {}
	TCtmsNet(PSignNet PSNet): TSignNet(*PSNet) {}
	static PCtmsNet New() { return new TCtmsNet;}
	static PCtmsNet Load(TSIn& SIn) { return new TCtmsNet(SIn); }
	PCtmsNet GetEdgeSubNet(const int& EDat1, const int& EDat2 = TInt::Mn, const int& EDat3 = TInt::Mn) const 
	{		
		return new TCtmsNet(TSignNet::GetEdgeSubNet(EDat1, EDat2, EDat3));
	}

	PCtmsNet GetSubNet(const int& MinEdgeWgt) const {
		return new TCtmsNet(TSignNet::GetSubNet(MinEdgeWgt));
	}
	PCtmsNet GetSignSubNet(const int& Sign) const {
		return new TCtmsNet(TSignNet::GetSignSubNet(Sign));
	}

	PCtmsNet GetTriad(const int& NId1, const int& NId2, const int& NId3) const {
		return new TCtmsNet(TSignNet::GetTriad(NId1, NId2, NId3));
	}
	
	static PCtmsNet LoadEpinions(const TStr& FNm) { return new TCtmsNet(TSignNet::LoadEpinions(FNm)); }
	static PCtmsNet LoadSlashdot(const TStr& FNm) { return new TCtmsNet(TSignNet::LoadSlashdot(FNm)); }
	static PCtmsNet GetSmallNet() { return new TCtmsNet(TSignNet::GetSmallNet()); }
	friend class TPt<TCtmsNet>;
	
	// implemented in ctmsnet.cpp	
	static PCtmsNet LoadSignedNet(const TStr& InFNm, const int& SrcColId = 0, const int& DstColId = 1, const int& SignColId = 2);
	PCtmsNet LoadEpinionsModified(const TStr& FNm);
	PCtmsNet CopyNet() const;
	PCtmsNet GetBalSgnProbSubNet();
	PCtmsNet GetMinEmbeddedSubNet(const int MinEmValue);
	bool operator==(const TCtmsNet& SgnNet) const;
	bool operator!=(const TCtmsNet& SgnNet) const;	
	PCtmsNet GenRewire(const int& NSwitch, TRnd Rnd = TRnd()) const;
	PCtmsNet GenRewire2Hash(const int& NSwitch, TRnd Rnd = TRnd()) const;
	PCtmsNet RewireNetwork2(const int &NSwitch = 0, const bool separateHash = true);
	PCtmsNet RewireNetwork(const int &NSwitch);
	PCtmsNet PermuteEdgeSigns2();
	PCtmsNet PermuteEdgeSignsStrict2();
	int GetStrictSgnPermSets(PCtmsNet& Net, THash<TIntPr, TIntVPr>& PermutableGroups);
	int GetAllEdgeSgnPermSets(PCtmsNet& Net, TSignNet::TEdgeI& EI, TVec<TIntPrPr>& ESgnPermSets);
	double GetEdgePosSgnProb(const int& SrcNdId, const int& DstNdId);
	bool GetRndESgnPermSet(PCtmsNet& Net, TSignNet::TEdgeI& EI, TIntPrPr& ERndPermSet);
	PCtmsNet PermuteEdgeSignsStrict();
	void CompareNetWith(const PCtmsNet& PNet2) const;
	TVec<TIntTr> GetEmbDistrb() const;
	bool IsSameTriad2(const PCtmsNet& Net1, const PCtmsNet& Net2);
	static int EdgeSig(int FwSgn, int BwSgn);
	static TChA GetTriadStr(int a[2], int b[2], int c[2]);
	static TChA GetTriadStr(const PCtmsNet& Nt, int srcId, int dstId, int nbrId, bool IsSigned = true);
	TChA GetTriadStr() const;
	void static GenTriadEquivClasses(TTriadEqClasH& TriadGroups, const bool BiDirEdgeSide = true, const bool Signed = true, const bool ZeroEdgeSide = false);
	double GetTriadProb2(const double& PlusProb) const;
	void CountSignedTriads2(const TStr& OutFNm) const;
	void CountSignedTriads3(const TStr& OutFNm) const;
	int static CalcNetTriadsDistrib(TTriadEqClasH& NetSigTriDistrib, PCtmsNet& Net);
};

#endif
