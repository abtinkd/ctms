#ifndef SIGN_PREDICTION_H
#define SIGN_PREDICTION_H

#include "ctmsnet.h"
#include <Snap.h>
#include "ml.h"
#define OutputDir "results/"

struct TPredictionResult {
	double accuracy, truePositive, trueNegative;
	long zeroFeaVCnt, incorrectZeroPredCnt, criticalPointCnt;
	TPredictionResult(double a=0.0, double tp=0.0, double tn=0.0, long z=0, long incz=0, long c=0) {
		accuracy = a; truePositive = tp; trueNegative = tn;
		zeroFeaVCnt = z; incorrectZeroPredCnt = incz; criticalPointCnt = c;
	}
};

class TSignPrediction {
protected:
	TStr NetName;
	PCtmsNet bNet;
	TIntTrV focusedNetEdges;
	THashSet<TChA> Features; //TChA is Feature string. Its corresponding Int ID is hash KeyId (int KeyId -> TChA Feature String)
	
public:
	static TPredictionResult getAccuracy(const TIntV& YV, const TFltV& ResultsV, const TFlt resValWhenZeroFeatV, const bool printDetailedResults = true);
};

class TSignPredictionNaive: public TSignPrediction {
protected:
	TFltV InDegBased, OutDegBased, weightedMeanBased, combinedDegBased, weightedCombDegBased;
	TIntV YV;
public:
	TSignPredictionNaive(const PCtmsNet& Net, const TStr NtNm, TIntTrV& selectedEdges) {
		NetName = NtNm;
		bNet = Net;
		for (int i = 0; i < selectedEdges.Len(); i++) {
			focusedNetEdges.Add(selectedEdges[i]);
			YV.Add(selectedEdges[i].Val3);
		}
		InDegBased.Gen(focusedNetEdges.Len());			InDegBased.PutAll(0);
		OutDegBased.Gen(focusedNetEdges.Len());			OutDegBased.PutAll(0);
		weightedMeanBased.Gen(focusedNetEdges.Len());	weightedMeanBased.PutAll(0);
		combinedDegBased.Gen(focusedNetEdges.Len());	combinedDegBased.PutAll(0);
		weightedCombDegBased.Gen(focusedNetEdges.Len());weightedCombDegBased.PutAll(0);
		CalcProbs();
	}

	void CalcProbs();
	double GetInDegBasedAcc(const bool printDetails = true) {return GetAccFor(InDegBased, "In-Deg-Based", printDetails);}
	double GetOutDegBasedAcc(const bool printDetails = true) {return GetAccFor(OutDegBased, "Out-Deg-Based", printDetails);}
	double GetweightedMeanBasedAcc(const bool printDetails = true) {return GetAccFor(weightedMeanBased, "Weighted-Mean-Based", printDetails);}
	double GetcombinedDegBasedAcc(const bool printDetails = true) {return GetAccFor(combinedDegBased, "Combined-Deg-Based", printDetails);}
	double GetweightedCombDegBasedAcc(const bool printDetails = true) {return GetAccFor(weightedCombDegBased, "Weighted-Comb-Deg-Based", printDetails);}

private:
	double GetAccFor(TFltV& probV, const TStr title, const bool printDetails) {
		TPredictionResult result = getAccuracy(YV, probV, 0.5, false);
		if (printDetails) {
			printf("\n_________________________________\n\n");
			printf("%s: %.4f\nACC2: %.4f\nTPR: %.4f TNR: %.4f Z: %d IZ: %d CP: %d\n", title.CStr(),
				result.accuracy, (result.truePositive + result.trueNegative)/2.0,
				result.truePositive, result.trueNegative,
				result.zeroFeaVCnt, result.incorrectZeroPredCnt, result.criticalPointCnt);
		}
		return result.accuracy;
	}
};

class TCrossValidation: public TSignPrediction {
private:
	TPredictionResult CalcCVFinalResult(TVec<TPredictionResult>& foldResV, char* title, const bool printDetails = true);
protected:
	virtual TPredictionResult TrainTestUseFeaV(TIntV& TrnIndxs, TIntV& TstIndxs) = 0;
	virtual TPredictionResult TrainTestUseNet(PCtmsNet& net, TIntV& TstIndxs) = 0;
	double RunValidation(const int NFold = 10, const bool alsoRunInSubNetMode = true);
	PCtmsNet GetNet(const TIntV& Indexes);
};

class TSignPredictionLearn: public TCrossValidation {
private:
	THash<TStr, TPair<TInt, TFltV>> DataSet;
	TVec<TFltV> XV;
	TFltV YV;
	TFltV Pred;
	PLogRegPredict LRModel;
	bool useNewtonMethod;
private:
	// My features
	void CreateFeatureV(const bool Signed = true, const bool BiDirEdgeSide = true);
	// Lescovek features
	void CreateFeatureVPp();
	// Abtin: Get all edges' corresponding triads with parameters: sign and birdirection
	void ExtractDataSet(const bool Signed = true, const bool BiDirEdgeSide = true);
	// Lescovek features
	void ExtractDataSetPp();
	void ScaleFeature();
	void SaveTheta(const TStr& Suffx = TStr());
protected:
	TPredictionResult TrainTestUseFeaV(TIntV& TrnIndxs, TIntV& TstIndxs);
	TPredictionResult TrainTestUseNet(PCtmsNet& net, TIntV& TstIndxs);
	void FitPredictionModel(const bool ScaleF = true, const bool Newton = true);
public:
	TSignPredictionLearn(const PCtmsNet& Net, const TStr NtNm) {
		NetName = NtNm;
		bNet = Net;
		for (TSignNet::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++)
			focusedNetEdges.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), EI()));
		ExtractDataSetPp();
		//ExtractDataSet();
	}
	TSignPredictionLearn(const PCtmsNet& Net, const TStr NtNm, TIntTrV& selectedEdges) {
		NetName = NtNm;
		bNet = Net;
		for (int i = 0; i < selectedEdges.Len(); i++) {focusedNetEdges.Add(selectedEdges[i]);}
		ExtractDataSetPp();
		//ExtractDataSet();
	}
	static TPredictionResult Test(const PLogRegPredict& LogRegM, const TVec<TFltV>& X, const TFltV& Y, TFltV& outCfys);
	static TPredictionResult Test(const TFltV& theta, const TVec<TFltV>& X, const TFltV& Y);
	const TPredictionResult Test(const PLogRegPredict& LogRegM);
	const TPredictionResult Test(const TFltV& theta) const;
	void FeatureSelection (const int FeCnt = 2);
	// Calculates and returns learned model theta
	void CalcTheta(TFltV& th, bool saveTheta = false);
	double CrossValidTest(const bool ScaleF = false, const bool Newton = true, const int NFold = 10);
	// Saves Extracted Dataset
	void SaveDataset(const TStr& FNameSuffx = TStr());
	void SaveYXV(const TStr& FNameSuffx = TStr());
	void SavePredictions(const TStr& FNameSuffx = TStr());
};

class TSignPredicNoLrn: public TCrossValidation {
private:
	THash<TIntTr, THash<TInt, TInt>> DataSet; // THash<TInt, TInt> : 1st TInt is Features KeyId (DataSet[x].Key <-> Features.GetKeyId(+o|-o|++)
	TVec<THash<TInt, TInt>> PXV, NXV; //instead of saving all feature values, here we save nonzero features
	TVec<TInt> YV;
	TVec<TFlt> NaivePred;
	TVec<TFltV> Pred;
	THash<TChA, TFlt> theta;
	double bNetPosProb;
private:
	static int EdgeSig(int FwSgn, int BwSgn);
	static TChA GetTriadStr(int a[2], int b[2], int c[2]);
	static TChA GetTriadStr(const PCtmsNet& Nt, int srcId, int dstId, int nbrId, bool IsSigned = true);
	static double CalcEffectiveEsPosProb(const PCtmsNet& Nt);
	TFlt GetEdgeNaivePrediction(const TIntPr& E) const;
	void CreateFeatureV();
	void ExtractDataSet();
	static void SaveTheta(const THash<TChA, TFlt>& thet, const TStr fileName);
	void SaveTheta();
	void SaveYXV();
	void SavePredictions();
protected:
	// Gets the probability of each given edge (specified by IndexV) being positive by OutV. The return value of function is the probabilty
	//  of an edge with no feature value to be positive
	double GetCfy(const TIntV& IndexV, TFltV& OutV, const THash<TChA, TFlt>& NewTheta);
	TPredictionResult Test(const THash<TChA, TFlt>& Th, const TIntV& TstEsIndexes);
	TPredictionResult TrainTestUseFeaV(TIntV& TrnIndxs, TIntV& TstIndxs);
	TPredictionResult TrainTestUseNet(PCtmsNet& net, TIntV& TstIndxs);
public:
	TSignPredicNoLrn(const PCtmsNet& Net, const TStr NtNm) {
		CreateFeatureV();
		NetName = NtNm;
		bNet = Net;
		for (TSignNet::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++)
			focusedNetEdges.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), EI()));
		ExtractDataSet();
	}
	TSignPredicNoLrn(const PCtmsNet& Net, const TStr NtNm, TIntTrV& selectedEdges) {
		CreateFeatureV();
		NetName = NtNm;
		bNet = Net;
		for (int i = 0; i < selectedEdges.Len(); i++) {focusedNetEdges.Add(selectedEdges[i]);}
		ExtractDataSet();
	}
	void SaveDataset();
	bool Empty() {return DataSet.Empty();}
	//Get whole network Theta
	void GetTheta(THash<TChA, TFlt>& thta);
	//Get Theta restricted to given edges
	void GetTheta(const TIntV& TrnIndices, THash<TChA, TFlt>& thta);
	//Similar to previous function but it considers edge embeddedness to calculate p0 and count CTMSs.
	//It performs better on balanced datasets and has higher precision.
	void GetThetaEmBased(const TIntV& TrnIndices, THash<TChA, TFlt>& thta);
	//Get nets Theta
	static void GetTheta(const PCtmsNet& Net, THash<TChA, TFlt>& thta, TStr saveFileName = "");
	static void normalize(THash<TChA, TFlt>& Th);
	static void GetFeatureV(THashSet<TChA>& Featu);
	//Test on whole network
	TPredictionResult Test(const THash<TChA, TFlt>& Th);
	double CrossValidTest(const bool CalSubNet = true, const int NFold = 10);
};
#endif // !SIGN_PREDICTION_H

