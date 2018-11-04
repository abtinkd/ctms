#include "sign_prediction.h"
#include <Snap.h>
#include <signnet.h>

TPredictionResult TSignPrediction::getAccuracy(const TIntV& YV, const TFltV& ResultsV, const TFlt resValWhenZeroFeatV, const bool printDetailedResults) {
	IAssert(YV.Len() == ResultsV.Len());
	for (int i = 0; i < ResultsV.Len(); i++) {
		IAssert(ResultsV[i] >= 0.0 && ResultsV[i] <= 1.0);
		IAssert(YV[i] == -1 || YV[i] == +1);
	}

	int correct = 0, negCor = 0, posCor = 0, Totneg = 0, Totpos = 0, zeroCnt = 0, incorrectZero = 0, criticalPoint = 0;

	for (int i = 0;  i < YV.Len(); i++) {
		if ((YV[i] == -1 && ResultsV[i] < 0.5) || (YV[i] == +1 && ResultsV[i] >= 0.5)) {correct++;}

		if (YV[i] == +1) {
			Totpos++;
			if (ResultsV[i] >= 0.5) {posCor++;}
		}
		if (YV[i] == -1) {
			Totneg++;
			if (ResultsV[i] < 0.5) {negCor++;}
		}

		if (ResultsV[i] == 0.5) {criticalPoint++;}
		if (ResultsV[i] == resValWhenZeroFeatV) {
			zeroCnt++;
			if ((YV[i] == +1 && ResultsV[i] < 0.5) || (YV[i] == -1 && ResultsV[i] >= 0.5)) {incorrectZero++;}
		}
	}

	TPredictionResult result;
	result.accuracy = double(correct)/(double)YV.Len();
	result.truePositive = (double)posCor/(double)Totpos;
	result.trueNegative = (double)negCor/(double)Totneg;
	result.zeroFeaVCnt = zeroCnt;
	result.incorrectZeroPredCnt = incorrectZero;
	result.criticalPointCnt = criticalPoint;

	if (printDetailedResults) {
		printf("\nTrue Positive ( %.2f ), True Negative ( %.2f ).\n", result.truePositive, result.trueNegative);
		printf("Empty FeaVec( %d ) | Incorrect zero estimate( %d ) | CRITICAL POINTS( %d ).\n",
			result.zeroFeaVCnt, result.incorrectZeroPredCnt, result.criticalPointCnt);
	}

	return result;
}

//////////////////////////////////////////////////////////////////////
//TCrossValidation
TPredictionResult TCrossValidation::CalcCVFinalResult(TVec<TPredictionResult>& foldResV, char* title, const bool printDetails) {
	const double fn = foldResV.Len();
	TPredictionResult result(0.0,0.0,0.0,0,0,0);
	for (int l = 0;  l < foldResV.Len(); l++) {
		result.accuracy += foldResV[l].accuracy;
		result.truePositive += foldResV[l].truePositive;
		result.trueNegative += foldResV[l].trueNegative;
		result.zeroFeaVCnt += foldResV[l].zeroFeaVCnt;
		result.incorrectZeroPredCnt += foldResV[l].incorrectZeroPredCnt;
		result.criticalPointCnt += foldResV[l].criticalPointCnt;
	}
	result.accuracy /= fn;
	result.truePositive /= fn;
	result.trueNegative /= fn;
	const double acc2 = (result.truePositive + result.trueNegative)/2.0;
	if (printDetails) {
		printf("\n_________________________________\n\n");
		printf("Mean Accuracy (%s): %.4f\nACC2: %.4f\nTPR: %.4f TNR: %.4f Z: %d IZ: %d CP: %d\n", title,
			result.accuracy, acc2,
			result.truePositive, result.trueNegative,
			result.zeroFeaVCnt, result.incorrectZeroPredCnt, result.criticalPointCnt);
	}
	return result;
}

double TCrossValidation::RunValidation(const int NFold, const bool alsoRunInSubNetMode) {
	TExeTm Tm;
	IAssert(!focusedNetEdges.Empty());

	TIntV Seq(focusedNetEdges.Len());
	for (int i = 0; i < Seq.Len(); i++) {Seq[i] = i;} //so Seq[0]==0, Seq[1]==1, Seq[2]==2, ...
	TInt::Rnd.Randomize();
	Seq.Shuffle(TInt::Rnd); //so Seq[0]!=0 (Rnd), Seq[1]!=1 (Rnd), Seq[2]!=2 (Rnd), ...

	TVec<TPredictionResult> CVRes, CVRes2;
	const int TestSz = focusedNetEdges.Len() / NFold;
	IAssertR(focusedNetEdges.Len() > NFold, "insufficient data for cross validation");
	int TstSrtIndx, TstEndIndx;

	for (int n = 0; n < NFold; n++) {
		printf("\nFold %d   --- --- --- --- --- --- --- --- --- --- --- ---  Time: %s\n", n+1, Tm.GetStr());
		// Setting NewXV NewYV
		TIntV TrainIndexes, TestIndexes;
		TstSrtIndx = n * TestSz;
		TstEndIndx = (n+1) * TestSz;
		printf("Recreating train & test data...");
		// Setting Test Data
		for (int t = TstSrtIndx; t < TstEndIndx; t++)
			TestIndexes.Add(Seq[t]);
		// Setting Train Data
		for (int t = 0; t < TstSrtIndx; t++)
			TrainIndexes.Add(Seq[t]);
		for (int t = TstEndIndx; t < Seq.Len(); t++)
			TrainIndexes.Add(Seq[t]);
		printf(" COMPLETED\n");
		if (alsoRunInSubNetMode)
		{
			printf("\nBuilding %d fold's Network ...\n",n+1);
			PSignNet foldNet = GetNet(TrainIndexes);
			printf("Using subnet ---> ");
			TPredictionResult res1 = TrainTestUseNet(foldNet, TestIndexes);
			CVRes.Add(res1);
			printf("Accuracy: %.4f\n", res1.accuracy);
		}
		{
			printf("\nOriginal data ---> ");
			TPredictionResult res2 = TrainTestUseFeaV(TrainIndexes, TestIndexes);
			CVRes2.Add(res2);
			printf("Accuracy: %.4f\n", res2.accuracy);
		}
	}
	if (alsoRunInSubNetMode)
		const TPredictionResult finalResult = CalcCVFinalResult(CVRes, "SubNet");
	const TPredictionResult finalResult2 = CalcCVFinalResult(CVRes2, "Origin");
	return finalResult2.accuracy;
}

PSignNet TCrossValidation::GetNet(const TIntV& Indexes) {
	PSignNet net = TSignNet::New();
	for (int i = 0; i < Indexes.Len(); i++) {
		const TIntTr CurE = focusedNetEdges[Indexes[i]];
		if (!net->IsNode(CurE.Val1)) {net->AddNode(CurE.Val1);}
		if (!net->IsNode(CurE.Val2)) {net->AddNode(CurE.Val2);}
		net->AddEdge(CurE.Val1, CurE.Val2, CurE.Val3);
	}
	return net;
}
//////////////////////////////////////////////////////////////////////
//TSignPredictionNaive
void TSignPredictionNaive::CalcProbs() {
	for (int v = 0; v < focusedNetEdges.Len(); v++) {
		const TIntTr Edge = focusedNetEdges[v];
		int SrcPosOutDeg = 0, DesPosInDeg = 0;
		const TSignNet::TNodeI SrcNI = bNet->GetNI(Edge.Val1);
		const TSignNet::TNodeI DesNI = bNet->GetNI(Edge.Val2);
		const int SrcOutDeg = SrcNI.GetOutDeg() - 1; // -1 for ignoring current edge
		const int DesInDeg  = DesNI.GetInDeg()  - 1; // -1 for ignoring current edge

		if (SrcOutDeg == 0 && DesInDeg == 0) {
			InDegBased[v] = OutDegBased[v] = weightedMeanBased[v] = combinedDegBased[v] = weightedCombDegBased[v] = 0.5;
			continue;
		}

		for (int i = 0;  i < SrcOutDeg+1; i++)
			if (bNet->GetEDat(SrcNI.GetId(), SrcNI.GetOutNId(i)) == +1)
				SrcPosOutDeg++;
		for (int i = 0;  i < DesInDeg+1; i++)
			if (bNet->GetEDat(DesNI.GetInNId(i), DesNI.GetId()) == +1)
				DesPosInDeg++;

		if(Edge.Val3 == +1) {
			SrcPosOutDeg--;
			DesPosInDeg--;
		}
		OutDegBased[v] = (float)(SrcPosOutDeg + 1) / (float)(SrcOutDeg + 2);
		InDegBased[v]  = (float)(DesPosInDeg + 1) / (float)(DesInDeg + 2);

		double w1 = 1 + sqrt(SrcOutDeg) * pow(OutDegBased[v]-0.5, 2);
		double w2 = 1 + sqrt(DesInDeg) * pow(InDegBased[v]-0.5, 2);
		weightedMeanBased[v] = (float)(w1*OutDegBased[v] + w2*InDegBased[v]) / (w1 + w2);

		combinedDegBased[v] = (float)(SrcPosOutDeg + DesPosInDeg + 1) / (float)(SrcOutDeg + DesInDeg + 2);

		const double alpha = sqrt((float)(SrcOutDeg + 1) / (float)(DesInDeg + 1)); // if we use sqrt then weightedCombDegBased ~ Mean
		weightedCombDegBased[v] = (float)((1/alpha)*SrcPosOutDeg + alpha*DesPosInDeg + 1) / 
			(float)((1/alpha)*SrcOutDeg + alpha*DesInDeg + 2);

		IAssertR(OutDegBased[v]>=0.0 && OutDegBased[v]<=1.0, OutDegBased[v].GetStr());
		IAssertR(InDegBased[v]>=0.0 && InDegBased[v]<=1.0, InDegBased[v].GetStr());
		IAssertR(weightedMeanBased[v]>=0.0 && weightedMeanBased[v]<=1.0, weightedMeanBased[v].GetStr());
		IAssertR(combinedDegBased[v]>=0.0 && combinedDegBased[v]<=1.0, combinedDegBased[v].GetStr());
		IAssertR(weightedCombDegBased[v]>=0.0 && weightedCombDegBased[v]<=1.0, weightedCombDegBased[v].GetStr());
	}
	return;
}

//////////////////////////////////////////////////////////////////////
//TSignPredictionLearn

/*
Describtion:
Ex. (A -N- B ) where (A,B) is the edge we want to estimate its sign and N is the neighbour
For making feature vector, we must first create all kinds of distinct triads with properties: Sign & Bidirection.
These traids can have zero edge side, because after removing the specified edge from triad, one side of triad may get lost; if it's
not bidirected!
Then we assure that that triad can happen in the features:
1. If the triad is complete (6 edges), it is not a feature, since we can't have a 6-edge triad feature after removing one
2. {A,N} and {B,N} must be present since we look for features in the (A,B)'s neighbourhood.
*/
void TSignPredictionLearn::CreateFeatureV(const bool Signed, const bool BiDirEdgeSide) {
	TTriadEqClasH fe;
	TSignNet::GenTriadEquivClasses(fe, BiDirEdgeSide, Signed, true);
	for (int c = 0;  c < fe.Len(); c++) {
		const PSignNet ClTriad = fe[c].Val1;
		if (ClTriad->GetEdges() == 6) {continue;} //refers to 1
		int SideECnt = 0;
		if (ClTriad->IsEdge(0, 1 , false)) {SideECnt++;} //refers to 2
		if (ClTriad->IsEdge(1, 2 , false)) {SideECnt++;} //refers to 2
		if (ClTriad->IsEdge(0, 2 , false)) {SideECnt++;} //refers to 2
		if (SideECnt != 0 && SideECnt != 1) {
			Features.AddKey(fe.GetKey(c));} //just corresponding strings is saved as features
	}
	return;
}

void TSignPredictionLearn::CreateFeatureVPp() {
	const TStr ES[] = {"F", "B"}, SS[] = {"p", "n"};
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					Features.AddKey(ES[i]+ES[j]+SS[k]+SS[l]);
				}
			}
		}
	}
	return;
}

void TSignPredictionLearn::FeatureSelection(const int FeCnt) {
	if (FeCnt < 1 || FeCnt >= Features.Len()) {return;}
	THash<TChA, TFlt> FeaturesCnt;
	for (int l = 0;  l < Features.Len(); l++) {FeaturesCnt.AddDat(Features.GetKey(l), 0);}
	for (int i = 0;  i < DataSet.Len(); i++) {
		IAssertR(DataSet[i].Val2.Len() == FeaturesCnt.Len(), "FeatureVec miss-match");
		for (int l = 0; l < DataSet[i].Val2.Len(); l++)	{
			FeaturesCnt(Features.GetKey(l)) += DataSet[i].Val2[l];}}
	FeaturesCnt.SortByDat(false);
	//Creating new feature vector
	THashSet<TChA> NwFeat;
	for (int cnt = 0;  cnt < FeCnt; cnt++) {
		NwFeat.AddKey(FeaturesCnt.GetKey(cnt));}

	THash<TStr, TPair<TInt, TFltV>> NewDataSet;
	XV.Clr();
	for (int i = 0;  i < DataSet.Len(); i++) {
		TFltV FeNwVals = TFltV(FeCnt);
		for (int l = 0; l < NwFeat.Len(); l++)	{
			const int FeMapId = Features.GetKeyId(NwFeat.GetKey(l));
			FeNwVals[l] = DataSet[i].Val2[FeMapId];
		}
		NewDataSet.AddDat(DataSet.GetKey(i), TPair<TInt, TFltV> (DataSet[i].Val1, FeNwVals));
		XV.Add(FeNwVals);
	}
	IAssert(XV.Len() == DataSet.Len()); IAssert(NewDataSet.Len() == DataSet.Len());
	Features.Clr();
	DataSet.Clr();
	Features = NwFeat;
	DataSet = NewDataSet;
	return;
}

void TSignPredictionLearn::CalcTheta(TFltV& th, bool saveTheta) {
	if (LRModel.Empty())
		FitPredictionModel(false, true);
	LRModel->GetTheta(th);
	if (saveTheta)
		SaveTheta();
	return;
}

void TSignPredictionLearn::ExtractDataSet(const bool Signed, const bool BiDirEdgeSide) {
	TIntV NbrsV;
	PSignNet SgnNet = bNet->CopyNet();

	CreateFeatureV(Signed, BiDirEdgeSide);
	if (!Signed) {SgnNet->SetAllEDat(1);}
	printf("Extracting DataSet. Signed:%d  BiDir:%d\n", Signed, BiDirEdgeSide);
	int c = 0, ECnt = (SgnNet->GetEdges()/100) +1;
	for (TSignNet::TEdgeI EI = SgnNet->BegEI(); EI < SgnNet->EndEI(); EI++) {
		const int ESign = bNet->GetEDat(EI.GetSrcNId(), EI.GetDstNId()); //(ESign instead of EI()) This line added for using unsigned features
		const TStr edge = TStr::Fmt("(%d, %d)", EI.GetSrcNId(), EI.GetDstNId());
		const int curId = DataSet.AddKey(edge);
		TFltV FeVals = TFltV(Features.Len()); FeVals.PutAll(0);
		TSnap::GetCmnNbrs(SgnNet, EI.GetSrcNId(), EI.GetDstNId(), NbrsV);
		for (int n = 0; n < NbrsV.Len(); n++) {
			PSignNet Tri = SgnNet->GetTriad(EI.GetSrcNId(), EI.GetDstNId(), NbrsV[n]);
			Tri->DelEdge(0, 1); //Deleting current edge
			if (BiDirEdgeSide) {
				const int FeId = Features.GetKeyId(Tri->GetTriadStr());
				IAssertR(FeId != -1, "Error! Found a triad not present in the features vector!");
				if (FeId != -1) {FeVals[FeId]++;}			
			}
			else //This condition is when bidirected edges is not allowed in features so we must consider all (A->N->B, A<-N->B, A->N<-B, A<-N<-B)
			{
				int E0_2[]={0, 2, 2, 0}, E1_2[]={1, 2, 2, 1};
				for (int i = 0; i < 3; i+=2) {
					for (int j = 0; j < 3; j+=2) {
						PSignNet tempTri = Tri->CopyNet();
						if (Tri->IsEdge(E0_2[i], E0_2[i+1]) && Tri->IsEdge(E0_2[i+1], E0_2[i])) {
							tempTri->DelEdge(E0_2[i], E0_2[i+1]);}
						if (Tri->IsEdge(E1_2[j], E1_2[j+1]) && Tri->IsEdge(E1_2[j+1], E1_2[j])) {
							tempTri->DelEdge(E1_2[j], E1_2[j+1]);}
						const int FeId = Features.GetKeyId(tempTri->GetTriadStr());
						IAssertR(FeId != -1, "Error! Found a triad not present in the features vector!");
						if (FeId != -1) {FeVals[FeId]++;}
					}
				}
			}
		}
		DataSet[curId] = TPair<TInt, TFltV> (ESign, FeVals);
		XV.Add(FeVals);
		YV.Add((ESign==-1)? 0 : 1);
		if (c++%ECnt == 0) {printf("\r%d%%",c/ECnt);}
	}
	printf("\n");
	SaveDataset();
	SaveYXV();
	Pred.Gen(focusedNetEdges.Len());
	return;
}

void TSignPredictionLearn::ExtractDataSetPp() {
	TExeTm Tm;
	TIntV NbrsV;

	CreateFeatureVPp();
	printf("Extracting DataSet (Only 16 Triads).\n");
	int c = 0, ECnt = (focusedNetEdges.Len()/100) +1;
	for (int i = 0; i < focusedNetEdges.Len(); i++) {
		const TIntTr CurE(focusedNetEdges[i]);
		const TStr edge = TStr::Fmt("(%d, %d)", CurE.Val1, CurE.Val2);
		const int curId = DataSet.AddKey(edge);
		TFltV FeVals = TFltV(Features.Len()); FeVals.PutAll(0);
		TSnap::GetCmnNbrs(bNet, CurE.Val1, CurE.Val2, NbrsV);
		for (int n = 0; n < NbrsV.Len(); n++) {
			// ->-> Fwd[0]Fwd[1] <--> Bwd[0]Fwd[1] ...
			bool Fwd[2] = {false, false}, Bwd[2] = {false, false};
			bool FwdSgnP[2] = {false, false}, BwdSgnP[2] = {false, false};
			if (bNet->IsEdge(CurE.Val1, NbrsV[n])) {Fwd[0] = true; 
			if(bNet->GetEDat(CurE.Val1, NbrsV[n]) == 1) {FwdSgnP[0] = true;}}
			if (bNet->IsEdge(NbrsV[n], CurE.Val1)) {Bwd[0] = true;
			if(bNet->GetEDat(NbrsV[n], CurE.Val1) == 1) {BwdSgnP[0] = true;}}
			if (bNet->IsEdge(NbrsV[n], CurE.Val2)) {Fwd[1] = true;
			if(bNet->GetEDat(NbrsV[n], CurE.Val2) == 1) {FwdSgnP[1] = true;}}
			if (bNet->IsEdge(CurE.Val2, NbrsV[n])) {Bwd[1] = true;
			if(bNet->GetEDat(CurE.Val2, NbrsV[n]) == 1) {BwdSgnP[1] = true;}}
			if (Fwd[0] && Fwd[1]) {
				TChA Feat = "FF";
				FwdSgnP[0] ? Feat += "p" : Feat += "n";
				FwdSgnP[1] ? Feat += "p" : Feat += "n";
				FeVals[Features.GetKeyId(Feat)]++;
			}
			if (Fwd[0] && Bwd[1]) {
				TChA Feat = "FB";
				FwdSgnP[0] ? Feat += "p" : Feat += "n";
				BwdSgnP[1] ? Feat += "p" : Feat += "n";
				FeVals[Features.GetKeyId(Feat)]++;
			}
			if (Bwd[0] && Fwd[1]) {
				TChA Feat = "BF";
				BwdSgnP[0] ? Feat += "p" : Feat += "n";
				FwdSgnP[1] ? Feat += "p" : Feat += "n";
				FeVals[Features.GetKeyId(Feat)]++;
			}
			if (Bwd[0] && Bwd[1]) {
				TChA Feat = "BB";
				BwdSgnP[0] ? Feat += "p" : Feat += "n";
				BwdSgnP[1] ? Feat += "p" : Feat += "n";
				FeVals[Features.GetKeyId(Feat)]++;
			}
		}
		DataSet[curId] = TPair<TInt, TFltV> (CurE.Val3, FeVals);
		XV.Add(FeVals);
		YV.Add((CurE.Val3==1)? 1.0 : -1.0);
		if (c++%ECnt == 0) {printf("\r%d%%",c/ECnt);}
	}
	printf("\r100%%  [%s]\n", Tm.GetStr());
	SaveDataset();
	SaveYXV();
	Pred.Gen(focusedNetEdges.Len());
	return;
}

void TSignPredictionLearn::ScaleFeature() {
	printf("Scaling feature Vector\n");
	TFltV MaxVal(Features.Len());
	MaxVal.PutAll(1);
	for (int l = 0; l < XV.Len(); l++) {
		for (int k = 0; k < XV[l].Len(); k++) {
			if (XV[l][k] > MaxVal[k]) {MaxVal[k] = XV[l][k];}
		}
	}

	for (int l = 0; l < XV.Len(); l++) {
		for (int k = 0; k < XV[l].Len(); k++) {
			XV[l][k] = XV[l][k]/MaxVal[k];
		}
	}
	return;
}

void TSignPredictionLearn::FitPredictionModel(const bool ScaleF, const bool Newton) {
	printf("Fitting Prediction Model...");
	if (ScaleF) {ScaleFeature();}
	TLogRegFit LRFit;
	if (Newton) {LRModel = LRFit.CalcLogRegNewton(XV, YV);}
	else {LRModel = LRFit.CalcLogRegGradient(XV, YV);}
	printf(" COMPLETED\n");
	return;
}

TPredictionResult TSignPredictionLearn::Test(const PLogRegPredict& LogRegM, const TVec<TFltV>& X, const TFltV& Y, TFltV& outCfys) {
	IAssertR(X.Len() == Y.Len(), "X and Y mismatch");
	TIntV YIntV(Y.Len());
	for (int i = 0; i < Y.Len(); i++)
		if (Y[i] > 0.0)
			YIntV[i] = 1;
		else
			YIntV[i] = -1;

	TFltV theta;
	LogRegM->GetTheta(theta);
	const double valueWhenZeroFeatV = (1 / (1 + exp(-theta.Last()))); // theta.Last() <-> theta0

	TFltV Results;
	LogRegM->GetCfy(X, Results);
	IAssert(X.Len() == Results.Len());
	outCfys = Results;
	return getAccuracy(YIntV, Results, valueWhenZeroFeatV);
}

TPredictionResult TSignPredictionLearn::Test(const TFltV& theta, const TVec<TFltV>& X, const TFltV& Y) {
	IAssertR(X.Len() == Y.Len(), "X and Y mismatch");
	const int feaSize = X[0].Len();
	IAssertR(theta.Len() == feaSize || theta.Len() == (feaSize+1), "Theta mismatch");
	for (int i = 1; i < X.Len(); i++)
		IAssertR(X[i].Len() == feaSize, "X vector incompatible");

	TIntV YIntV(Y.Len());
	for (int i = 0; i < Y.Len(); i++)
		if (Y[i] >= 0.0)
			YIntV[i] = 1;
		else
			YIntV[i] = -1;

	TFltV Results;
	TLogRegPredict::GetCfy(X, Results, theta);
	const double valueWhenZeroFeatV = (1 / (1 + exp(-theta.Last()))); // theta.Last() <-> theta0
	return getAccuracy(YIntV, Results, valueWhenZeroFeatV);
}

const TPredictionResult TSignPredictionLearn::Test(const PLogRegPredict& LogRegM) {
	IAssert(!LogRegM.Empty());
	return TSignPredictionLearn::Test(LogRegM, XV, YV, Pred);
}

const TPredictionResult TSignPredictionLearn::Test(const TFltV& theta) const {
	IAssert(!theta.Empty());
	return TSignPredictionLearn::Test(theta, XV, YV);
}

double TSignPredictionLearn::CrossValidTest(const bool ScaleF, const bool Newton, const int NFold) {
	TExeTm Tm;
	IAssert(!XV.Empty());
	IAssert(XV.Len() == YV.Len());
	if (ScaleF)
		ScaleFeature();
	useNewtonMethod = Newton;
	double acc = RunValidation(NFold, false);
	SavePredictions();
	return acc;
}

TPredictionResult TSignPredictionLearn::TrainTestUseFeaV(TIntV& TrnIndxs, TIntV& TstIndxs) {
	TVec<TFltV> TestXV, TrainXV;
	TFltV TestYV, TrainYV;
	for (int i = 0; i < TrnIndxs.Len(); i++) {
		TrainXV.Add(XV[TrnIndxs[i]]);
		TrainYV.Add(YV[TrnIndxs[i]]);
	}
	for (int i = 0; i < TstIndxs.Len(); i++) {
		TestXV.Add(XV[TstIndxs[i]]);
		TestYV.Add(YV[TstIndxs[i]]);
	}
	printf("Fitting prediction model ...\n");
	TLogRegFit LRFit;
	LRModel.Clr();
	if (useNewtonMethod)
		LRModel = LRFit.CalcLogRegNewton(TrainXV, TrainYV);
	else
		LRModel = LRFit.CalcLogRegGradient(TrainXV, TrainYV);

	TFltV cfyRes;
	TPredictionResult res = Test(LRModel, TestXV, TestYV, cfyRes);
	SaveTheta();

	IAssert(cfyRes.Len() == TstIndxs.Len());
	for (int i = 0; i < cfyRes.Len(); i++)
		Pred[TstIndxs[i]] = cfyRes[i];
	return res;
}

TPredictionResult TSignPredictionLearn::TrainTestUseNet(PSignNet& net, TIntV& TstIndxs) {
	return TPredictionResult();
}

void TSignPredictionLearn::SaveDataset(const TStr& Suffx) {
	FILE * S = fopen((OutputDir+NetName+"DSet"+Suffx+".data").CStr(), "wt");
	fprintf(S, "#Edge\tLabel\t");
	for (int f = 0; f < Features.Len(); f++) {fprintf(S, "{%s}\t", (Features[f]).CStr());}
	for (int d = 0; d < DataSet.Len(); d++) {
		fprintf(S,"\n");
		fprintf(S, "%s\t%d\t", DataSet.GetKey(d).CStr(), DataSet[d].Val1);
		for (int f = 0; f < DataSet[d].Val2.Len(); f++) {fprintf(S, "%.0f\t", DataSet[d].Val2[f]);}
	}
	fclose(S);
}

void TSignPredictionLearn::SaveYXV(const TStr& Suffx) {
	FILE * S = fopen((OutputDir+NetName+"Y_X"+Suffx+".data").CStr(), "wt");
	fprintf(S, "#Label\t");
	for (int f = 0; f < Features.Len(); f++) {fprintf(S, "{%s}\t", (Features[f]).CStr());}
	for (int d = 0; d < XV.Len(); d++) {
		fprintf(S,"\n%.0f\t", YV[d]);
		for (int f = 0; f < XV[d].Len(); f++) {fprintf(S, "%.1f\t", XV[d][f]);}
	}
	fclose(S);
}

void TSignPredictionLearn::SaveTheta(const TStr& Suffx) {
	FILE * S = fopen((OutputDir+NetName+"Theta"+Suffx+".data").CStr(), "a+t");
	TFltV thet; LRModel->GetTheta(thet);
	fprintf(S, "#Theta0\t");
	for (int f = 0; f < Features.Len(); f++) {fprintf(S, "{%s}\t", (Features[f]).CStr());}
	fprintf(S, "\n%f\t", thet[thet.Len()-1]); //theta0
	for (int t = 0; t < thet.Len()-1; t++) {fprintf(S, "%f\t", thet[t]);}
	fclose(S);
}

void TSignPredictionLearn::SavePredictions(const TStr& Suffx) {
	FILE * S = fopen((OutputDir+NetName+"Predictions"+Suffx+".data").CStr(), "wt");
	fprintf(S, "#Edge\tLabel\tPrediction");
	for (int p = 0; p < Pred.Len(); p++)
		fprintf(S,"\n(%d, %d)\t%d\t%f", focusedNetEdges[p].Val1, focusedNetEdges[p].Val2, focusedNetEdges[p].Val3, Pred[p]);
	fclose(S);
}

//////////////////////////////////////////////////////////////////////
//TSignPredicNoLrn

int TSignPredicNoLrn::EdgeSig(int e1, int e2) {
	if (e1==0 && e2==0) { return 0; }
	if (e1==1 && e2==0) { return 1; }
	if (e1==-1 && e2==0) { return 2; }
	if (e1==0 && e2==1) { return 3; }
	if (e1==0 && e2==-1) { return 4; }
	if (e1==1 && e2==1) { return 5; }
	if (e1==1 && e2==-1) { return 6; }
	if (e1==-1 && e2==1) { return 7; }
	if (e1==-1 && e2==-1) { return 8; }
	Fail; return -1;
}

TChA TSignPredicNoLrn::GetTriadStr(int a[2], int b[2], int c[2]) {
	const int e1 = EdgeSig(a[0], a[1]);
	const int e2 = EdgeSig(b[0], b[1]);
	const int e3 = EdgeSig(c[0], c[1]);
	const TIntTr MnTr1 = TMath::Mn(TIntTr(e1, e2, e3), TIntTr(e2, e3, e1), TIntTr(e3, e1, e2));
	const int e1a = EdgeSig(c[1], c[0]);
	const int e2a = EdgeSig(b[1], b[0]);
	const int e3a = EdgeSig(a[1], a[0]);
	const TIntTr MinTr = TMath::Mn(MnTr1, TMath::Mn(TIntTr(e1a, e2a, e3a), TIntTr(e2a, e3a, e1a), TIntTr(e3a, e1a, e2a)));

	const TInt EdgNum[] = {MinTr.Val1, MinTr.Val2, MinTr.Val3};
	TChA EdgStr[3];
	for (int i = 0; i < 3; i++) {
		switch(EdgNum[i]) {
		case 1 : EdgStr[i] = "+o"; break;
		case 2 : EdgStr[i] = "-o"; break;
		case 3 : EdgStr[i] = "o+"; break;
		case 4 : EdgStr[i] = "o-"; break;
		case 5 : EdgStr[i] = "++"; break;
		case 6 : EdgStr[i] = "+-"; break;
		case 7 : EdgStr[i] = "-+"; break;
		case 8 : EdgStr[i] = "--"; break;
		default : EdgStr[i] = TChA();
		}
	}
	return EdgStr[0]+"|"+EdgStr[1]+"|"+EdgStr[2];
}

TChA TSignPredicNoLrn::GetTriadStr(const PSignNet& Nt, int srcId, int dstId, int nbrId, bool IsSigned) {
	int SDSgn[2] = {0, 0}, DNSgn[2] = {0, 0}, NSSgn[2] = {0, 0};
	if (IsSigned) {
		if (Nt->IsEdge(srcId, dstId)) {
			SDSgn[0] = Nt->GetEDat(srcId, dstId);}
		if (Nt->IsEdge(dstId, srcId)) {
			SDSgn[1] = Nt->GetEDat(dstId, srcId);}
		if (Nt->IsEdge(dstId, nbrId)) {
			DNSgn[0] = Nt->GetEDat(dstId, nbrId);}
		if (Nt->IsEdge(nbrId, dstId)) {
			DNSgn[1] = Nt->GetEDat(nbrId, dstId);}
		if (Nt->IsEdge(nbrId, srcId)) {
			NSSgn[0] = Nt->GetEDat(nbrId, srcId);}
		if (Nt->IsEdge(srcId, nbrId)) {
			NSSgn[1] = Nt->GetEDat(srcId, nbrId);}
	}
	else {
		if (Nt->IsEdge(srcId, dstId)) {
			SDSgn[0] = 1;}
		if (Nt->IsEdge(dstId, srcId)) {
			SDSgn[1] = 1;}
		if (Nt->IsEdge(dstId, nbrId)) {
			DNSgn[0] = 1;}
		if (Nt->IsEdge(nbrId, dstId)) {
			DNSgn[1] = 1;}
		if (Nt->IsEdge(nbrId, srcId)) {
			NSSgn[0] = 1;}
		if (Nt->IsEdge(srcId, nbrId)) {
			NSSgn[1] = 1;}
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
	//return TSignNet::GetTriadStr(MinTr);
	const TInt EdgNum[] = {MinTr.Val1, MinTr.Val2, MinTr.Val3};
	TChA EdgStr[3];
	for (int i = 0; i < 3; i++) {
		switch(EdgNum[i]) {
		case 1 : EdgStr[i] = "+o"; break;
		case 2 : EdgStr[i] = "-o"; break;
		case 3 : EdgStr[i] = "o+"; break;
		case 4 : EdgStr[i] = "o-"; break;
		case 5 : EdgStr[i] = "++"; break;
		case 6 : EdgStr[i] = "+-"; break;
		case 7 : EdgStr[i] = "-+"; break;
		case 8 : EdgStr[i] = "--"; break;
		default : EdgStr[i] = TChA();
		}
	}
	return EdgStr[0]+"|"+EdgStr[1]+"|"+EdgStr[2];
}

double TSignPredicNoLrn::CalcEffectiveEsPosProb(const PSignNet& Nt) {
	int baseNetNONZeroEmbEs = 0;
	int posE = 0;
	for (TSignNet::TEdgeI EI = Nt->BegEI(); EI < Nt->EndEI(); EI++) { // instead of using TSnap::GetCmnNbrs() to increase performance
		TSignNet::TNodeI srcNI = Nt->GetNI(EI.GetSrcNId());
		TSignNet::TNodeI desNI = Nt->GetNI(EI.GetDstNId());
		for (int i = 0; i < srcNI.GetDeg(); i++) {
			if (desNI.IsNbrNId(srcNI.GetNbrNId(i))) { // so feature vector has at least 1 nonzero entry
				baseNetNONZeroEmbEs++;
				if (EI()==1)
					posE++;
				break;
			}
		}
	}
	return posE/(double)baseNetNONZeroEmbEs;
}

TFlt TSignPredicNoLrn::GetEdgeNaivePrediction(const TIntPr& Edge) const{
	int SrcPosOutDeg = 0, DesPosInDeg = 0;

	const TSignNet::TNodeI SrcNI = bNet->GetNI(Edge.Val1);
	const TSignNet::TNodeI DesNI = bNet->GetNI(Edge.Val2);
	const int SrcOutDeg = SrcNI.GetOutDeg() - 1; // -1 for ignoring current edge
	const int DesInDeg  = DesNI.GetInDeg()  - 1; // -1 for ignoring current edge

	for (int i = 0;  i < SrcOutDeg+1; i++)
		if (bNet->GetEDat(SrcNI.GetId(), SrcNI.GetOutNId(i)) == +1)
			SrcPosOutDeg++;
	for (int i = 0;  i < DesInDeg+1; i++)
		if (bNet->GetEDat(DesNI.GetInNId(i), DesNI.GetId()) == +1)
			DesPosInDeg++;

	if(bNet->GetEDat(Edge.Val1, Edge.Val2) == +1) {
		SrcPosOutDeg--;
		DesPosInDeg--;
	}

	const double OutDegBased = (float)(SrcPosOutDeg + 1) / (float)(SrcOutDeg + 2);
	const double InDegBased  = (float)(DesPosInDeg + 1) / (float)(DesInDeg + 2);

	double w1 = 1 + sqrt(SrcOutDeg) * pow(OutDegBased-0.5, 2);
	double w2 = 1 + sqrt(DesInDeg) * pow(InDegBased-0.5, 2);
	const double weightedMeanBased = (float)(w1*OutDegBased + w2*InDegBased) / (w1 + w2);

	IAssert(weightedMeanBased >=0.0 && weightedMeanBased <=1.0);
	return weightedMeanBased ;
}

void TSignPredicNoLrn::CreateFeatureV() {
	GetFeatureV(Features);
}

void TSignPredicNoLrn::GetFeatureV(THashSet<TChA>& Featu) {
	static TTriadEqClasH fe;
	if (fe.Empty())
		TSignNet::GenTriadEquivClasses(fe);
	Featu.AddKey("theta0");
	for (int c = 0;  c < fe.Len(); c++) //id <- c+1, triad string <- fe.GetKey(c)
		Featu.AddKey(fe.GetKey(c));
	return;
}

void TSignPredicNoLrn::GetTheta(THash<TChA, TFlt>& thta) {
	TIntV Idx;
	if (theta.Empty()) {
		Idx.Gen(focusedNetEdges.Len());
		for (int i = 0; i < focusedNetEdges.Len(); i++)
			Idx[i] = i;
		GetTheta(Idx, theta);
	}
	thta = theta;
	return;
}

TIntPr GetClassPosNegECount(const TChA TriClassStr, const THashSet<TChA>& Feas) {
	static THash<TChA, TIntPr> AllClassesPosNegEs;
	if (AllClassesPosNegEs.Len() != Feas.Len()) {
		for (int f = 0; f < Feas.Len(); f++) {
			const TChA clssStr = Feas.GetKey(f);
			TIntPr PNCount(0,0);
			for(int i = 0; i < clssStr.Len(); i++) {
				if (clssStr[i] == '+') PNCount.Val1++;
				if (clssStr[i] == '-') PNCount.Val2++;				
			}
			AllClassesPosNegEs.AddKey(clssStr);
			AllClassesPosNegEs(clssStr).Val1 = PNCount.Val1;
			AllClassesPosNegEs(clssStr).Val2 = PNCount.Val2;
		}
	}

	return AllClassesPosNegEs(TriClassStr);
}

void TSignPredicNoLrn::GetTheta(const TIntV& TrnIndices, THash<TChA, TFlt>& thta) {
	TExeTm Tm;
	printf("Calculating Theta using new function ...\n");
	THash<TChA, TInt> actualCTMSsCnt, CTMSsIfPos, CTMSsIfNeg;

	for (int i = 0; i < Features.Len(); i++) {
		const TChA feaStr = Features.GetKey(i);
		thta.AddDat(feaStr, 0.0);
		actualCTMSsCnt.AddDat(feaStr, 0);
		CTMSsIfPos.AddDat(feaStr, 0);
		CTMSsIfNeg.AddDat(feaStr, 0);
	}

	TIntV NbrV;
	int PosEListSize = 0, AllTriadsCnt = 0;
	int c=0, Decile=int(TrnIndices.Len()/100.0)+1; //Abtin: +1 is to prevent docile get zero.
	for (int i = 0; i < TrnIndices.Len(); i++) {
		const TIntTr CurE(focusedNetEdges[TrnIndices[i]]);
		IAssertR(DataSet.IsKey(CurE), "Edge does not exist in the DataSet!");
		const int EId = DataSet.GetKeyId(CurE);
		if (CurE.Val3 == 1) {PosEListSize++;}

		const THash<TInt, TInt>& curEFeaVals = DataSet[EId];
		for (int j = 0; j < curEFeaVals.Len(); j++) {
			const TChA feaStr = Features[curEFeaVals.GetKey(j)];
			actualCTMSsCnt(feaStr) += curEFeaVals[j];
			AllTriadsCnt += curEFeaVals[j];
		}

		const THash<TInt, TInt>& posEFeaVals = PXV[TrnIndices[i]];
		for (int j = 0; j < posEFeaVals.Len(); j++) {
			const TChA feaStr = Features[posEFeaVals.GetKey(j)];
			CTMSsIfPos(feaStr) += posEFeaVals[j];
		}

		const THash<TInt, TInt>& negEFeaVals = NXV[TrnIndices[i]];
		for (int j = 0; j < negEFeaVals.Len(); j++) {
			const TChA feaStr = Features[negEFeaVals.GetKey(j)];
			CTMSsIfNeg(feaStr) += negEFeaVals[j];
		}

		if (++c % Decile ==0) { printf("\r%.0f%%",(double)c/Decile); }
	}
	double p0 = PosEListSize/(double)TrnIndices.Len();

	for (int t = 0; t < thta.Len(); t++) {
		const TChA feaStr = thta.GetKey(t);
		const double ExpCnt = 1 + p0 * CTMSsIfPos(feaStr) + (1-p0) * CTMSsIfNeg(feaStr);
		const double TriadProb = ExpCnt / (double)AllTriadsCnt;
		const double Surp = (actualCTMSsCnt(feaStr)-ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0-TriadProb));
		thta[t] = Surp;
	}
	thta("theta0") = 0.0;
	SaveTheta(thta, NetName);
	normalize(thta);
	SaveTheta(thta, NetName + "Normal");
	printf("[%s]\n", Tm.GetStr());
	return;
}

//Similar to previous function but it considers edge embeddedness to calculate p0 and count CTMSs.
void TSignPredicNoLrn::GetThetaEmBased(const TIntV& TrnIndices, THash<TChA, TFlt>& thta) {
	TExeTm Tm;
	printf("Calculating Theta using new function(2) ...\n");
	THash<TChA, TInt> actualCTMSsCnt, CTMSs;
	THash<TInt, TIntPr> edgeEmDistrib;
	THash<TInt, THash<TChA, TInt>> EmCTMSsIfPos, EmCTMSsIfNeg;

	for (int i = 0; i < Features.Len(); i++) {
		const TChA feaStr = Features.GetKey(i);
		thta.AddDat(feaStr, 0.0);
		actualCTMSsCnt.AddDat(feaStr, 0);
		CTMSs.AddDat(feaStr, 0);
	}

	TIntV NbrV;
	int PosEListSize = 0, AllTriadsCnt = 0;
	int c=0, Decile=int(TrnIndices.Len()/100.0)+1; //Abtin: +1 is to prevent docile get zero.
	for (int i = 0; i < TrnIndices.Len(); i++) {
		const TIntTr CurE(focusedNetEdges[TrnIndices[i]]);
		IAssertR(DataSet.IsKey(CurE), "Edge does not exist in the DataSet!");
		const int EId = DataSet.GetKeyId(CurE);

		const THash<TInt, TInt>& curEFeaVals = DataSet[EId];
		int emb = 0;
		for (int j = 0; j < curEFeaVals.Len(); j++) {
			const TChA feaStr = Features[curEFeaVals.GetKey(j)];
			actualCTMSsCnt(feaStr) += curEFeaVals[j];
			AllTriadsCnt += curEFeaVals[j];
			emb += curEFeaVals[j];
		}

		if (!edgeEmDistrib.IsKey(emb)) {
			edgeEmDistrib.AddDat(emb, TIntPr(0, 0));
			EmCTMSsIfPos.AddDat(emb, CTMSs);
			EmCTMSsIfNeg.AddDat(emb, CTMSs);
		}

		if (CurE.Val3 == 1) {
			edgeEmDistrib(emb).Val1 += 1;}
		edgeEmDistrib(emb).Val2 += 1;

		const THash<TInt, TInt>& posEFeaVals = PXV[TrnIndices[i]];
		for (int j = 0; j < posEFeaVals.Len(); j++) {
			const TChA feaStr = Features[posEFeaVals.GetKey(j)];
			EmCTMSsIfPos(emb)(feaStr) += posEFeaVals[j];
		}

		const THash<TInt, TInt>& negEFeaVals = NXV[TrnIndices[i]];
		for (int j = 0; j < negEFeaVals.Len(); j++) {
			const TChA feaStr = Features[negEFeaVals.GetKey(j)];
			EmCTMSsIfNeg(emb)(feaStr) += negEFeaVals[j];
		}

		if (++c % Decile ==0) { printf("\r%.0f%%",(double)c/Decile); }
	}

	for (int t = 0; t < thta.Len(); t++) {
		const TChA feaStr = thta.GetKey(t);
		double ExpCnt = 1;
		for (int m = 0; m < edgeEmDistrib.Len(); m++) {
			const TInt em = edgeEmDistrib.GetKey(m);
			double p0 = edgeEmDistrib(em).Val1/(double)edgeEmDistrib(em).Val2;
			ExpCnt += p0 * EmCTMSsIfPos(em)(feaStr) + (1-p0) * EmCTMSsIfNeg(em)(feaStr);
		}
		const double TriadProb = ExpCnt / (double)AllTriadsCnt;
		const double Surp = (actualCTMSsCnt(feaStr)-ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0-TriadProb));
		thta[t] = Surp;
	}
	thta("theta0") = 0.0;
	SaveTheta(thta, NetName);
	normalize(thta);
	SaveTheta(thta, NetName + "Normal");
	printf("[%s]\n", Tm.GetStr());
	return;
}
//

void TSignPredicNoLrn::GetTheta(const PSignNet& Net, THash<TChA, TFlt>& thta, TStr saveFileName) {
	TExeTm Tm;
	THashSet<TChA> Features;
	GetFeatureV(Features);
	printf("Calculating Theta ...\n");
	static TTriadEqClasH SgnTriCnt, UnsgnTriCnt;
	if (SgnTriCnt.Empty()) {TSignNet::GenTriadEquivClasses(SgnTriCnt);}
	if (UnsgnTriCnt.Empty()) {TSignNet::GenTriadEquivClasses(UnsgnTriCnt, true, false);}
	for (int i = 0; i < Features.Len(); i++) {
		thta.AddDat(Features.GetKey(i), 0);}
	for (int i = 0; i < SgnTriCnt.Len(); i++) {
		SgnTriCnt[i].Val2 = 0;}
	for (int i = 0; i < UnsgnTriCnt.Len(); i++) {
		UnsgnTriCnt[i].Val2 = 0;}

	TIntV NbrV;
	double AllPlusE=0, AllE = Net->GetEdges();
	int  emptyFeVEs = 0, emptyFeVPosEs = 0; // for calculating theta0
	//Abtin added
	int AllTriadsCnt = 0;
	int c=0, Decile=int(AllE/100)+1; //Abtin: +1 is to prevent docile get zero.
	for (TSignNet::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) {
		TSnap::GetCmnNbrs(Net, EI.GetSrcNId(), EI.GetDstNId(), NbrV);
		if (NbrV.Len() == 0) { // all features were zero so this edge has no triad around and no common neighbour
			emptyFeVEs++;
			if (EI() == 1) {emptyFeVPosEs++;}
		}
		for (int n = 0; n < NbrV.Len(); n++) {
			const TChA TriStr = GetTriadStr(Net, EI.GetSrcNId(), EI.GetDstNId(), NbrV[n]);
			// count signed triad
			SgnTriCnt(TriStr).Val2 += 1;
			// count unsigned triads
			const TChA TriStrUnsgn = GetTriadStr(Net, EI.GetSrcNId(), EI.GetDstNId(), NbrV[n], false);
			UnsgnTriCnt(TriStrUnsgn).Val2 +=1;
		}
		if (EI() == 1) { AllPlusE += 1; }
		if (++c % Decile ==0) { printf("\r%.0f%%",(double)c/Decile); }
	}
	const double PlusProb = AllPlusE / AllE;
	int SigTriadsClasPresent = 0, UnSigTriadsClasPresent = 0;
	//Abtin: Each Triad is counted according to the number of edges present in it
	for (int t = 0; t < SgnTriCnt.Len(); t++) {
		const int E = SgnTriCnt[t].Val1->GetEdges();
		SgnTriCnt[t].Val2 /= E; 
		AllTriadsCnt += SgnTriCnt[t].Val2;
	}
	//Abtin: Each Triad is counted according to the number of edges present in it
	for (int u = 0; u < UnsgnTriCnt.Len(); u++) {
		const int E = UnsgnTriCnt[u].Val1->GetEdges();
		UnsgnTriCnt[u].Val2 /= E;
	}
	for (int t = 0; t < SgnTriCnt.Len(); t++) {
		//PSignNet TriadNetUnsig = SgnTriCnt[t].Val1->CopyNet(); TriadNetUnsig->SetAllEDat(1);
		TIntV TriNdV;
		SgnTriCnt[t].Val1->GetNIdV(TriNdV);
		const TChA Sgn2UnTriStr = GetTriadStr(SgnTriCnt[t].Val1, TriNdV[0], TriNdV[1], TriNdV[2], false);
		const double TriadCnt = SgnTriCnt[t].Val2;
		const double UnSignCnt = 1 + UnsgnTriCnt(Sgn2UnTriStr).Val2; //optional 1 added
		const double TriadProbInEqCls = SgnTriCnt[t].Val1->GetTriadProb2(PlusProb);
		const double TriadProb = TriadProbInEqCls * (UnSignCnt/(AllTriadsCnt + 1)); //optional 1 added
		const double ExpCnt = TriadProbInEqCls * UnSignCnt; // is equal to: TriadProb * AllTriadsCnt
		const double Surp = (TriadCnt-ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0-TriadProb));
		thta(SgnTriCnt.GetKey(t)) = Surp;
	}
	if (emptyFeVEs != 0)
		thta("theta0") = (emptyFeVPosEs/(double)emptyFeVEs) - 0.5;
	else
		thta("theta0") = PlusProb - 0.5;
	if(saveFileName != "") 
		SaveTheta(thta, saveFileName);
	normalize(thta);
	if(saveFileName != "") 
		SaveTheta(thta, saveFileName + "Normal");
	printf("  [%s]\n", Tm.GetStr());
	return;
}

void TSignPredicNoLrn::ExtractDataSet() {
	TExeTm Tm;
	TIntV NbrsV;
	printf("Extracting DataSet.. \n");
	PXV.Gen(focusedNetEdges.Len());
	NXV.Gen(focusedNetEdges.Len());
	YV.Gen(focusedNetEdges.Len());
	NaivePred.Gen(focusedNetEdges.Len());
	Pred.Gen(focusedNetEdges.Len());
	int c = 0, ECnt = (focusedNetEdges.Len()/100) +1;
	for (int i = 0; i < focusedNetEdges.Len(); i++) {
		const TIntTr CurE(focusedNetEdges[i]);
		THash<TInt, TInt> FeVals;
		THash<TInt, TInt> PosFeVals;
		THash<TInt, TInt> NegFeVals;
		TSnap::GetCmnNbrs(bNet, CurE.Val1, CurE.Val2, NbrsV);
		for (int n = 0; n < NbrsV.Len(); n++) {
			PSignNet ctms = bNet->GetTriad(CurE.Val1, CurE.Val2, NbrsV[n]); //E(CurE.Val1, CurE.Val2) --> E(0,1)
			const TChA TriStr = GetTriadStr(ctms, 0, 1, 2);
			TChA PTriStr, NTriStr;
			if (CurE.Val3 == 1) {
				PTriStr = TriStr;
				ctms->SetEDat(0, 1, -1); //Make edge sign negative
				NTriStr = GetTriadStr(ctms, 0, 1, 2);
			}
			else //CurE.Val3 == -1
			{
				NTriStr = TriStr;
				ctms->SetEDat(0, 1, 1); //Make edge sign positive
				PTriStr = GetTriadStr(ctms, 0, 1, 2);
			}

			FeVals.IsKey(Features.GetKeyId(TriStr)) ? FeVals(Features.GetKeyId(TriStr))++ : FeVals.AddDat(Features.GetKeyId(TriStr), 1);
			PosFeVals.IsKey(Features.GetKeyId(PTriStr)) ? PosFeVals(Features.GetKeyId(PTriStr))++ : PosFeVals.AddDat(Features.GetKeyId(PTriStr), 1);
			NegFeVals.IsKey(Features.GetKeyId(NTriStr)) ? NegFeVals(Features.GetKeyId(NTriStr))++ : NegFeVals.AddDat(Features.GetKeyId(NTriStr), 1);
		}
		DataSet.AddDat(CurE, FeVals);
		PXV[i] = PosFeVals;
		NXV[i] = NegFeVals;
		YV[i] = (CurE.Val3 == +1? +1: -1);
		NaivePred[i] = GetEdgeNaivePrediction(TIntPr(CurE.Val1, CurE.Val2));
		if (c++%ECnt == 0) {printf("\r%d%%",c/ECnt);}
	}
	bNetPosProb = CalcEffectiveEsPosProb(bNet);
	printf("\r100%%  [%s]\n", Tm.GetStr());
	SaveDataset();
	SaveYXV();
	return;
}

double TSignPredicNoLrn::GetCfy(const TIntV& IndexV, TFltV& OutV, const THash<TChA, TFlt>& NewTheta) {
	for (int i = 0;  i < IndexV.Len(); i++) {
		double pos = 0, neg = 0;
		long nonZeros = 0;
		const TInt curEIndex = IndexV[i];
		IAssert(curEIndex < PXV.Len());
		IAssert(curEIndex < NXV.Len());

		for (int j = 0; j < PXV[curEIndex].Len(); j++) {
			const TChA FeaStrToIntMap = Features[PXV[curEIndex].GetKey(j)];
			pos += NewTheta.GetDat(FeaStrToIntMap) * (double)PXV[curEIndex][j];
			nonZeros += PXV[curEIndex][j];
		}
		for (int j = 0; j < NXV[curEIndex].Len(); j++) {
			const TChA FeaStrToIntMap = Features[NXV[curEIndex].GetKey(j)];
			neg += NewTheta.GetDat(FeaStrToIntMap) * (double)NXV[curEIndex][j];
		}
		// -0.5< theta <0.5 => -1 < {(pos-neg)/nonZeros = pos/nonZeros - neg/nonZeros} < 1
		const double res = NewTheta.GetDat("theta0") + (pos - neg) * (5.0/(nonZeros!=0? nonZeros : 1));
		double mu = 1.0 / (1.0 + exp(-res));
		if (nonZeros == 0)
			mu = NaivePred[curEIndex];
		//	mu = (mu*nonZeros + NaivePred[curEIndex]) /(nonZeros + 1);
		Pred[curEIndex].Add(mu);
		OutV.Add(mu);
	}
	const double resWhenEmptyFeaV = 1.0 / (1.0 + exp(-NewTheta.GetDat("theta0")));
	return resWhenEmptyFeaV;
}

// Maps theta values to (-0.5, 0.5):	min -> -0.5		 max -> 0.5
void TSignPredicNoLrn::normalize(THash<TChA, TFlt>& Th) {
	TFlt min = 0, max = 0, median = 0, mean = 0;
	for (int i = 0; i < Th.Len(); i++) {
		if (Th[i] < min)
			min = Th[i];
		if (Th[i] > max)
			max = Th[i];
	}
	median = (max + min)/2;
	const TFlt theta0 = Th.GetDat("theta0");
	for (int i = 0; i < Th.Len(); i++) {
		Th[i] = (Th[i] - median)/(max - min);
		mean += Th[i];
	}
	Th("theta0") = mean / Th.Len();
	return;
}

TPredictionResult TSignPredicNoLrn::Test(const THash<TChA, TFlt>& Th, const TIntV& TstEsIndexes) {
	IAssertR(!Th.Empty(), "Theta is empty");
	IAssertR(!TstEsIndexes.Empty(), "Test vector is empty");

	TFltV Results;
	const TFlt valWhenEmptyFeatV = GetCfy(TstEsIndexes, Results, Th);
	TIntV YIntV(TstEsIndexes.Len());
	for (int i = 0; i < TstEsIndexes.Len(); i++)
		YIntV[i] = YV[TstEsIndexes[i]];

	return getAccuracy(YIntV, Results, valWhenEmptyFeatV);
}

TPredictionResult TSignPredicNoLrn::Test(const THash<TChA, TFlt>& Th) {
	TIntV seq(focusedNetEdges.Len());
	for (int i = 0; i < focusedNetEdges.Len(); i++)
		seq[i] = i;
	return Test(Th, seq);
}

double TSignPredicNoLrn::CrossValidTest(const bool CalSubNet, const int NFold) {
	// you can add normalization and feature selection here
	double acc = RunValidation(NFold, CalSubNet);
	SavePredictions();
	return acc;
}

TPredictionResult TSignPredicNoLrn::TrainTestUseFeaV(TIntV& TrnIndxs, TIntV& TstIndxs) {
	THash<TChA, TFlt> th;
	GetTheta(TrnIndxs, th);
	TPredictionResult res = Test(th, TstIndxs);
	return res;
}

TPredictionResult TSignPredicNoLrn::TrainTestUseNet(PSignNet& net, TIntV& TstIndxs) {
	THash<TChA, TFlt> th;
	GetTheta(net, th, NetName+"SUBNETUsed");
	TPredictionResult res = Test(th, TstIndxs);
	return res;
}

void TSignPredicNoLrn::SaveTheta(const THash<TChA, TFlt>& thet, const TStr fileName){
	FILE * S = fopen((OutputDir+fileName+"Theta.data").CStr(), "a+t");
	for (int f = 0; f < thet.Len(); f++) {fprintf(S, "{%s}\t", thet.GetKey(f).CStr());}
	fprintf(S, "\n");
	for (int t = 0; t < thet.Len(); t++) {fprintf(S, "%f\t", thet[t]);}
	fclose(S);
}

void TSignPredicNoLrn::SaveTheta() {
	SaveTheta(theta, NetName);
}

void TSignPredicNoLrn::SaveYXV() {
	FILE * S = fopen((OutputDir+NetName+"Y_XP.data").CStr(), "wt");
	FILE * L = fopen((OutputDir+NetName+"Y_XN.data").CStr(), "wt");
	fprintf(S, "#Edge\tLabel\t");
	fprintf(L, "#Edge\tLabel\t");
	for (int f = 0; f < Features.Len(); f++) {
		fprintf(S, "{%s}\t", (Features[f]).CStr());
		fprintf(L, "{%s}\t", (Features[f]).CStr());
	}
	for (int d = 0; d < PXV.Len(); d++) {
		fprintf(S,"\n(%d, %d)\t%d\t", focusedNetEdges[d].Val1, focusedNetEdges[d].Val2, focusedNetEdges[d].Val3);
		for (int f = 0; f < Features.Len(); f++) {
			const int feVal = (PXV[d].IsKey(f) ? PXV[d](f) : 0);
			fprintf(S, "%d\t", feVal);
		}
	}
	fclose(S);
	for (int d = 0; d < NXV.Len(); d++) {
		fprintf(L,"\n(%d, %d)\t%d\t", focusedNetEdges[d].Val1, focusedNetEdges[d].Val2, focusedNetEdges[d].Val3);
		for (int f = 0; f < Features.Len(); f++) {
			const int feVal = (NXV[d].IsKey(f) ? NXV[d](f) : 0);
			fprintf(L, "%d\t", feVal);
		}
	}
	fclose(L);
	return;
}

void TSignPredicNoLrn::SaveDataset() {
	FILE * S = fopen((OutputDir+NetName+"DSet.data").CStr(), "wt");
	fprintf(S, "#Edge\tLabel\t");
	for (int f = 0; f < Features.Len(); f++)
		fprintf(S, "{%s}\t", (Features[f]).CStr());
	for (int d = 0; d < DataSet.Len(); d++) {
		fprintf(S,"\n");
		fprintf(S, "(%d, %d)\t%d\t", DataSet.GetKey(d).Val1, DataSet.GetKey(d).Val2, DataSet.GetKey(d).Val3);
		for (int f = 0; f < Features.Len(); f++) {
			const int feVal = (DataSet[d].IsKey(f) ? DataSet[d](f) : 0);
			fprintf(S, "%d\t", feVal);
		}
	}
	fclose(S);
}

void TSignPredicNoLrn::SavePredictions() {
	FILE * S = fopen((OutputDir+NetName+"Predictions.data").CStr(), "wt");
	fprintf(S, "#Edge\tLabel\tPredictions...");
	for (int p = 0; p < Pred.Len(); p++) {
		fprintf(S,"\n(%d, %d)\t%d\t", focusedNetEdges[p].Val1, focusedNetEdges[p].Val2, focusedNetEdges[p].Val3);
		for (int v = 0; v < Pred[p].Len(); v++) {
			fprintf(S, "%f\t", Pred[p][v]);
		}
	}
	fclose(S);
}