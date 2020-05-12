#include "sign_prediction2.h"

// Maps theta values to (-0.5, 0.5):	min -> -0.5		 max -> 0.5
void normalize(THash<TChA, TFlt>& Th) {
	TFlt min = 0, max = 0, median = 0, mean = 0;
	for (int i = 0; i < Th.Len(); i++) {
		if (Th[i] < min)
			min = Th[i];
		if (Th[i] > max)
			max = Th[i];
	}
	median = (max + min) / 2;
	const TFlt theta0 = Th.GetDat("theta0");
	for (int i = 0; i < Th.Len(); i++) {
		Th[i] = (Th[i] - median) / (max - min);
		mean += Th[i];
	}
	Th("theta0") = mean / Th.Len();
	return;
}

void TNaivePredictor::build() {
	return;
}

void TBalanceBasedPredictor::build() {
	return;
}

void TStatusBasedPredictor::build() {
	return;
}

void TLogisticRegression::mapFeatureToIndex() {
	const TStr ES[] = { "F", "B" }, SS[] = { "p", "n" };
	TInt index = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					const TStr feature = ES[i] + ES[j] + SS[k] + SS[l];
					fe2ix.AddKey(feature);
					fe2ix.AddDat(feature, index);
					ix2fe.AddKey(index);
					ix2fe.AddDat(index, feature);
					index++;
				}
			}
		}
	}
	return;
}

void TLogisticRegression::build() {
	TVec<TFltV> XV;
	TFltV YV;
	printf("Extracting feature values...\n");
	int c = 0, ECnt = (network->GetEdges() / 100) + 1;
	for (TSignNet::TEdgeI EI = network->BegEI(); EI < network->EndEI(); EI++) {
		const TIntPr edge(EI.GetSrcNId(), EI.GetDstNId());
		THash<TChA, TInt> edgeFeaValues;
		extractFeatures(edge, edgeFeaValues);
		TFltV values = TFltV(fe2ix.Len()); values.PutAll(0.0);
		for (int i = 0; i < edgeFeaValues.Len(); i++) {			
			const TChA feature = edgeFeaValues.GetKey(i);
			values[fe2ix.GetDat(feature)] = edgeFeaValues.GetDat(feature);
		}
		XV.Add(values);
		const int sign = network->GetEDat(edge.Val1, edge.Val2);		
		YV.Add((sign == -1) ? 0 : 1);
		if (c++%ECnt == 0) {printf("\r%d%%", c / ECnt);}
	}
	printf("\r%d%%\n", 100);
	printf("Fitting predictive model...\n");	
	TLogRegFit LRFit;
	PLogRegPredict LRModel = LRFit.CalcLogRegNewton(XV, YV);
	printf(" Done.\n");
	
	LRModel->GetTheta(theta);	
	return;
}

void TCTMSProbabilisticInference::build() {	
	TTriadEqClasH fe;
	TCtmsNet::GenTriadEquivClasses(fe);
	THash<TChA, TInt> actualCTMSsCnt, CTMSsIfPos, CTMSsIfNeg;
	theta.AddDat("theta0", 0.0);
	for (int c = 0; c < fe.Len(); c++) { //id <- c+1, triad string <- fe.GetKey(c)
		const TChA feaStr = fe.GetKey(c);
		theta.AddDat(feaStr, 0.0);
		actualCTMSsCnt.AddDat(feaStr, 0);
		CTMSsIfPos.AddDat(feaStr, 0);
		CTMSsIfNeg.AddDat(feaStr, 0);
	}
	int c = 0, ECnt = (network->GetEdges() / 100) + 1;
	int PosEListSize = 0, AllTriadsCnt = 0;
	for (TSignNet::TEdgeI EI = network->BegEI(); EI < network->EndEI(); EI++) {
		const TIntPr edge(EI.GetSrcNId(), EI.GetDstNId());
		const TInt sign = EI.GetDat();
		if (sign == +1) PosEListSize++;
		for (int i = 0; i < 3; i++) {
			if (i == 1)
				network->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), +1);
			else if (i == 2)
				network->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), -1);
			THash<TChA, TInt> feaValues;
			extractFeatures(edge, feaValues);
			for (THashKeyDatI<TChA, TInt> fv = feaValues.BegI(); fv < feaValues.EndI(); fv++) {
				if (i == 0) {
					AllTriadsCnt += fv.GetDat();
					actualCTMSsCnt(fv.GetKey()) += fv.GetDat();
				}
				else if (i == 1)
					CTMSsIfPos(fv.GetKey()) += fv.GetDat();
				else if (i == 2)
					CTMSsIfNeg(fv.GetKey()) += fv.GetDat();
			}
			network->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), sign);
		}
		if (c++%ECnt == 0) { printf("\r%d%%", c / ECnt); }
	}
	double p0 = PosEListSize / (double)network->GetEdges();
	for (int t = 0; t < theta.Len(); t++) {
		const TChA feaStr = theta.GetKey(t);
		const double ExpCnt = 1 + p0 * CTMSsIfPos(feaStr) + (1 - p0) * CTMSsIfNeg(feaStr);
		const double TriadProb = ExpCnt / (double)AllTriadsCnt;
		const double Surp = (actualCTMSsCnt(feaStr) - ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0 - TriadProb));
		theta[t] = Surp;
	}
	normalize(theta);	
	return;
}

void TCTMSProbabilisticInferenceLocal::build() {
	return;
}

int TNaivePredictor::predict(const TInt srcNId, const TInt desNId) {
	int sign;
	const TSignNet::TNodeI SrcNI = network->GetNI(srcNId);
	const TSignNet::TNodeI DesNI = network->GetNI(desNId);
	const int SrcOutDeg = SrcNI.GetOutDeg();
	const int DesInDeg = DesNI.GetInDeg();

	int SrcPosOutDeg = 0, DesPosInDeg = 0;
	for (int i = 0; i < SrcOutDeg; i++)
		if (network->GetEDat(SrcNI.GetId(), SrcNI.GetOutNId(i)) == +1)
			SrcPosOutDeg++;
	for (int i = 0; i < DesInDeg; i++)
		if (network->GetEDat(DesNI.GetInNId(i), DesNI.GetId()) == +1)
			DesPosInDeg++;

	const double gnr = (SrcPosOutDeg + 1.0) / (SrcOutDeg + 2.0);
	const double rcp = (DesPosInDeg + 1.0) / (DesInDeg + 2.0);	
	IAssertR(gnr >= 0.0 && gnr <= 1.0, "invalid gnr");
	IAssertR(rcp >= 0.0 && rcp <= 1.0, "invalid rcp");
	if (pType == "gnr"){				
		sign = (gnr >= 0.5) ? 1 : -1;
	}
	else if (pType == "rcp") {				
		sign = (rcp >= 0.5) ? 1 : -1;
	}
	else if (pType == "cmp") {
		const double cmp = (SrcPosOutDeg + DesPosInDeg + 1.0) / (SrcOutDeg + DesInDeg + 2.0);
		IAssertR(cmp >= 0.0 && cmp <= 1.0, "invalid cmp");
		sign = (cmp >= 0.5) ? 1 : -1;
	}
	else if (pType == "wgr") {
		const double w1 = 1 + sqrt(SrcOutDeg) * pow(gnr - 0.5, 2);
		const double w2 = 1 + sqrt(DesInDeg) * pow(rcp - 0.5, 2);
		const double wgr = (float)(w1*gnr + w2*rcp) / (w1 + w2);
		IAssertR(wgr >= 0.0 && wgr <= 1.0, "invalid wgr");
		sign = (wgr >= 0.5) ? 1 : -1;
	}
	else {
		TInt::Rnd.Randomize();
		sign = (TInt::Rnd.GetUniDevInt(2) == 1) ? 1 : -1;
	}
	return sign;
}

int TBalanceBasedPredictor::predict(const TInt srcNId, const TInt desNId) {
	// ref: snap-exp\signnet.h > PredictBalance
	TIntV NbrV;
	TSnap::GetCmnNbrs(network, srcNId, desNId, NbrV);
	int Bal = 0;
	for (int n = 0; n < NbrV.Len(); n++) {		
		const int E1 = network->IsEdge(srcNId, NbrV[n]) ? network->GetEDat(srcNId, NbrV[n]) : network->GetEDat(NbrV[n], srcNId);
		const int E2 = network->IsEdge(desNId, NbrV[n]) ? network->GetEDat(desNId, NbrV[n]) : network->GetEDat(NbrV[n], desNId);
		if (E1*E2 == 1) { Bal++; }
	}
	return (Bal >= NbrV.Len() / 2) ? +1 : -1;	
}

int TStatusBasedPredictor::predict(const TInt srcNId, const TInt desNId) {
	// ref: snap-exp\signnet.h > PredictStatus
	TIntV NbrV;
	TSnap::GetCmnNbrs(network, srcNId, desNId, NbrV);
	double SrcS = 0, DstS = 0;
	for (int n = 0; n < NbrV.Len(); n++) {
		PSignNet Network = &*network;
		PSignNet TriadNet = Network->GetTriad(srcNId, desNId, NbrV[n]); // (src,dst) == (0,1)
		SrcS += TSignMicroEvol::GetStatus(TriadNet, 0) >= 0 ? 1 : -1;
		DstS += TSignMicroEvol::GetStatus(TriadNet, 1) >= 0 ? 1 : -1;
	}
	return (SrcS > DstS) ? -1 : +1;	
}

int TLogisticRegression::predict(const TInt srcNId, const TInt desNId) {
	if (!network->IsNode(srcNId) || !network->IsNode(desNId))
		return 0;
	const TIntPr edge(srcNId, desNId);
	THash<TChA, TInt> edgeFeaValues;
	extractFeatures(edge, edgeFeaValues);
	TFltV values = TFltV(fe2ix.Len()); values.PutAll(0.0);
	for (int i = 0; i < edgeFeaValues.Len(); i++) {
		const TChA feature = edgeFeaValues.GetKey(i);
		values[fe2ix.GetDat(feature)] = edgeFeaValues.GetDat(feature);
	}
	TVec<TFltV> X;
	X.Add(values);
	TFltV Result;
	TLogRegPredict::GetCfy(X, Result, theta);
	return (Result[0] < 0.50000000) ? -1 : +1;
}

int TCTMSProbabilisticInference::predict(const TInt srcNId, const TInt desNId) {
	if (!network->IsNode(srcNId) || !network->IsNode(desNId))
		return 0;
	const TIntPr edge(srcNId, desNId);
	network->AddEdge(srcNId, desNId, +1);
	THash<TChA, TInt> fvIfPositive;
	extractFeatures(edge, fvIfPositive);
	network->SetEDat(srcNId, desNId, -1);
	THash<TChA, TInt> fvIfNegative;
	extractFeatures(edge, fvIfNegative);
	network->DelEdge(srcNId,desNId);
	
	double posTrend = theta.GetDat("theta0"), negTrend = theta.GetDat("theta0");
	for (THashKeyDatI<TChA, TInt> fv = fvIfPositive.BegI(); fv < fvIfPositive.EndI(); fv++)
		posTrend += theta.GetDat(fv.GetKey()) * fv.GetDat();
	for (THashKeyDatI<TChA, TInt> fv = fvIfNegative.BegI(); fv < fvIfNegative.EndI(); fv++)
		negTrend += theta.GetDat(fv.GetKey()) * fv.GetDat();
	return posTrend >= negTrend ? 1 : -1;
}

int TCTMSProbabilisticInferenceLocal::predict(const TInt srcNId, const TInt desNId) {
	if (!network->IsNode(srcNId) || !network->IsNode(desNId))
		return 0;
	
	// calculating theta based on local data
	TIntPrV edges;
	const TSignNet::TNodeI SrcNI = network->GetNI(srcNId);
	const TSignNet::TNodeI DesNI = network->GetNI(desNId);
	for (int i = 0; i < SrcNI.GetOutDeg(); i++)
		edges.Add(TIntPr(SrcNI.GetId(), SrcNI.GetOutNId(i)));
	for (int i = 0; i < DesNI.GetOutDeg(); i++)
		edges.Add(TIntPr(DesNI.GetId(), DesNI.GetOutNId(i)));
	for (int i = 0; i < SrcNI.GetInDeg(); i++)
		edges.Add(TIntPr(SrcNI.GetInNId(i), SrcNI.GetId()));
	for (int i = 0; i < DesNI.GetInDeg(); i++)
		edges.Add(TIntPr(DesNI.GetInNId(i), DesNI.GetId()));

	static TTriadEqClasH fe;
	if (fe.Empty())
		TCtmsNet::GenTriadEquivClasses(fe);	
	THash<TChA, TInt> actualCTMSsCnt, CTMSsIfPos, CTMSsIfNeg;
	THash<TChA, TFlt> theta;
	theta.AddDat("theta0", 0.0);
	for (int c = 0; c < fe.Len(); c++) { //id <- c+1, triad string <- fe.GetKey(c)
		const TChA feaStr = fe.GetKey(c);
		theta.AddDat(feaStr, 0.0);
		actualCTMSsCnt.AddDat(feaStr, 0);
		CTMSsIfPos.AddDat(feaStr, 0);
		CTMSsIfNeg.AddDat(feaStr, 0);
	}	
	int PosEListSize = 0, AllTriadsCnt = 0;
	for (int i = 0; i < edges.Len(); i++) {
		const TIntPr edge(edges[i].Val1, edges[i].Val2);
		const TInt sign = network->GetEDat(edge.Val1, edge.Val2);
		if (sign == +1) PosEListSize++;
		for (int i = 0; i < 3; i++) {
			if (i == 1)
				network->SetEDat(edge.Val1, edge.Val2, +1);
			else if (i == 2)
				network->SetEDat(edge.Val1, edge.Val2, -1);
			THash<TChA, TInt> feaValues;
			extractFeatures(edge, feaValues);
			for (THashKeyDatI<TChA, TInt> fv = feaValues.BegI(); fv < feaValues.EndI(); fv++) {
				if (i == 0) {
					AllTriadsCnt += fv.GetDat();
					actualCTMSsCnt(fv.GetKey()) += fv.GetDat();
				}
				else if (i == 1)
					CTMSsIfPos(fv.GetKey()) += fv.GetDat();
				else if (i == 2)
					CTMSsIfNeg(fv.GetKey()) += fv.GetDat();
			}
			network->SetEDat(edge.Val1, edge.Val2, sign);
		}		
	}
	double p0 = PosEListSize / (double)edges.Len();
	for (int t = 0; t < theta.Len(); t++) {
		const TChA feaStr = theta.GetKey(t);
		const double ExpCnt = 1 + p0 * CTMSsIfPos(feaStr) + (1 - p0) * CTMSsIfNeg(feaStr);
		const double TriadProb = ExpCnt / (double)AllTriadsCnt;
		const double Surp = (actualCTMSsCnt(feaStr) - ExpCnt) / sqrt(AllTriadsCnt*TriadProb*(1.0 - TriadProb));
		theta[t] = Surp;
	}
	normalize(theta);
	// end of calculation

	const TIntPr edge(srcNId, desNId);
	network->AddEdge(srcNId, desNId, +1);
	THash<TChA, TInt> fvIfPositive;
	extractFeatures(edge, fvIfPositive);
	network->SetEDat(srcNId, desNId, -1);
	THash<TChA, TInt> fvIfNegative;
	extractFeatures(edge, fvIfNegative);
	network->DelEdge(srcNId, desNId);

	double posTrend = theta.GetDat("theta0"), negTrend = theta.GetDat("theta0");
	for (THashKeyDatI<TChA, TInt> fv = fvIfPositive.BegI(); fv < fvIfPositive.EndI(); fv++)
		posTrend += theta.GetDat(fv.GetKey()) * fv.GetDat();
	for (THashKeyDatI<TChA, TInt> fv = fvIfNegative.BegI(); fv < fvIfNegative.EndI(); fv++)
		negTrend += theta.GetDat(fv.GetKey()) * fv.GetDat();
	return posTrend >= negTrend ? 1 : -1;
}

void TLogisticRegression::extractFeatures(const TIntPr& edge,
	THash<TChA, TInt>& edgeFeaValues) {	
	TIntV NbrsV;
	TSnap::GetCmnNbrs(network, edge.Val1, edge.Val2, NbrsV);
	for (int n = 0; n < NbrsV.Len(); n++) {
		// ->-> Fwd[0]Fwd[1] <--> Bwd[0]Fwd[1] ...
		bool Fwd[2] = { false, false }, Bwd[2] = { false, false };
		bool FwdSgnP[2] = { false, false }, BwdSgnP[2] = { false, false };
		if (network->IsEdge(edge.Val1, NbrsV[n])) {
			Fwd[0] = true;
			if (network->GetEDat(edge.Val1, NbrsV[n]) == 1) { FwdSgnP[0] = true; }
		}
		if (network->IsEdge(NbrsV[n], edge.Val1)) {
			Bwd[0] = true;
			if (network->GetEDat(NbrsV[n], edge.Val1) == 1) { BwdSgnP[0] = true; }
		}
		if (network->IsEdge(NbrsV[n], edge.Val2)) {
			Fwd[1] = true;
			if (network->GetEDat(NbrsV[n], edge.Val2) == 1) { FwdSgnP[1] = true; }
		}
		if (network->IsEdge(edge.Val2, NbrsV[n])) {
			Bwd[1] = true;
			if (network->GetEDat(edge.Val2, NbrsV[n]) == 1) { BwdSgnP[1] = true; }
		}
		TChA Feat;
		if (Fwd[0] && Fwd[1]) {
			Feat = "FF";
			FwdSgnP[0] ? Feat += "p" : Feat += "n";
			FwdSgnP[1] ? Feat += "p" : Feat += "n";
		}
		if (Fwd[0] && Bwd[1]) {
			Feat = "FB";
			FwdSgnP[0] ? Feat += "p" : Feat += "n";
			BwdSgnP[1] ? Feat += "p" : Feat += "n";
		}
		if (Bwd[0] && Fwd[1]) {
			Feat = "BF";
			BwdSgnP[0] ? Feat += "p" : Feat += "n";
			FwdSgnP[1] ? Feat += "p" : Feat += "n";
		}
		if (Bwd[0] && Bwd[1]) {
			Feat = "BB";
			BwdSgnP[0] ? Feat += "p" : Feat += "n";
			BwdSgnP[1] ? Feat += "p" : Feat += "n";
		}
		edgeFeaValues.IsKey(Feat) ?
			edgeFeaValues(Feat)++ :
			edgeFeaValues.AddDat(Feat, 1);
	}
	return;
}

void TCTMSProbabilisticInference::extractFeatures(const TIntPr& edge, 
	THash<TChA, TInt>& edgeFeaValues) {	
	TIntV NbrsV;
	TSnap::GetCmnNbrs(network, edge.Val1, edge.Val2, NbrsV);
	for (int n = 0; n < NbrsV.Len(); n++) {
		PCtmsNet ctms = network->GetTriad(edge.Val1, edge.Val2, NbrsV[n]);
		const TChA TriStr = TCtmsNet::GetTriadStr(ctms, 0, 1, 2);
		edgeFeaValues.IsKey(TriStr) ? 
			edgeFeaValues(TriStr)++ : 
			edgeFeaValues.AddDat(TriStr, 1);
	}
}

void TCTMSProbabilisticInferenceLocal::extractFeatures(const TIntPr& edge,
	THash<TChA, TInt>& edgeFeaValues) {
	TIntV NbrsV;
	TSnap::GetCmnNbrs(network, edge.Val1, edge.Val2, NbrsV);
	for (int n = 0; n < NbrsV.Len(); n++) {
		PCtmsNet ctms = network->GetTriad(edge.Val1, edge.Val2, NbrsV[n]);
		const TChA TriStr = TCtmsNet::GetTriadStr(ctms, 0, 1, 2);
		edgeFeaValues.IsKey(TriStr) ?
			edgeFeaValues(TriStr)++ :
			edgeFeaValues.AddDat(TriStr, 1);
	}
}
