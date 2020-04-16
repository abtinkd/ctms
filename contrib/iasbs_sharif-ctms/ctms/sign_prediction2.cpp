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

}

void TStatusBasedPredictor::build() {

}

void TLogisticRegression::build() {

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
	int PosEListSize = 0, AllTriadsCnt = 0;
	for (TSignNet::TEdgeI EI = network->BegEI(); EI < network->EndEI(); EI++) {
		for (int i = 0; i < 3; i++) {
			const TIntPr edge(EI.GetSrcNId(), EI.GetDstNId());
			const TInt sign = EI.GetDat();			
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
					if (sign == +1) PosEListSize++;
				}
				else if (i == 1)
					CTMSsIfPos(fv.GetKey()) += fv.GetDat();
				else if (i == 2)
					CTMSsIfNeg(fv.GetKey()) += fv.GetDat();
			}
			network->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), sign);
		}		
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
	return 1;
}

int TStatusBasedPredictor::predict(const TInt srcNId, const TInt desNId) {
	return 1;
}

int TLogisticRegression::predict(const TInt srcNId, const TInt desNId) {
	return 1;
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

void::TLogisticRegression::extractFeatures() {

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
