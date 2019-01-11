#include <iostream>
#include <fstream>
#include "ctmsnet.h"
#include "ml.h"
#include "sign_prediction.h"
#include "convert_wikipedia.h"

// Define drand48 in MS Visual Studio as (0<rnd<1)
#define drand48() ((double)rand()/(RAND_MAX+1))
// Input File name and Location
#define LogDir "../Experimental Results/Temp/"
#define FileNEpinions "soc-sign-epinions.txt"
#define FileNSlashdot "soc-sign-Slashdot090221.txt"
#define FileNWikipedia "wikipedia.txt"
#define FileNWikipediaRawFormat "wikiElec.ElecBs3.txt"
#define FileNTest "test.txt"
#define DATE "94.05.19"
#define ConsoleToFile true

using namespace std;

// Print network statistics
void PrintNStats(const char NetName[], const PCtmsNet SignedNet) {
	cout << endl;
	cout << "Printing Networks Statistics (PrintNStats) ---------- " << NetName << endl;
	cout << "Nodes Count: " << SignedNet->GetNodes() << endl;
	cout << "Edges Count: " << SignedNet->GetEdges() << endl;
	cout << "Zero Nodes: " << SignedNet->GetNodes() - TSnap::CntNonZNodes(SignedNet) << endl;
	cout << "Self Edges: " << TSnap::CntSelfEdges(SignedNet) << endl;
	cout << "Unique Bi-direction Edges: " << TSnap::CntUniqBiDirEdges(SignedNet) << endl;
	//cout << "Average Clustering Coefficient: " << TSnap::GetClustCf(SignedNet) << endl;
	SignedNet->EdgeSignStat();
	cout << "------------------------------------------ ---------- -----------------" << endl;
}

// Test the default constructor
void DefaultConstructor() {
	int NNodes = 119217;
	int NEdges = 841200;
	int x,y;
	PCtmsNet SignedNet;

	SignedNet = TCtmsNet::New();
	
	for(int i = 0; i < NNodes; i++) {
		SignedNet->AddNode(i);
	}

	for(int i = 0; i < NEdges;) {
		x = (long) (drand48() * NNodes);
		y = (long) (drand48() * NNodes);
		if(x != y && !SignedNet->IsEdge(x,y)) {
			SignedNet->AddEdge(x,y);
			i++;
		}
	}

	PrintNStats("DefaultConstructor:Net", SignedNet);
}

void RemoveSelfEdges(PCtmsNet& sNet) {
	for(TSignNet::TEdgeI EI = sNet->BegEI(); EI < sNet->EndEI(); EI++) {
		int src = EI.GetSrcNId(), dst = EI.GetDstNId();
		if (src == dst) {
			sNet->DelEdge(src, dst);
		}
	}
	for(TSignNet::TNodeI NI = sNet->BegNI(); NI < sNet->EndNI(); NI++) {
		if (NI.GetDeg() == 0) {
			sNet->DelNode(NI.GetId());
		}
	}
}

void LoadSNet(PCtmsNet& signedNet, const char FLocation[]) {
	//TFIn FIn(FileLocation); signedNet = TCtmsNet::Load(FIn);
	//signedNet = TCtmsNet::LoadEpinions(FileLocation);
	//signedNet = TCtmsNet::LoadEpinionsModified(FileLocation);
	signedNet = TCtmsNet::LoadSignedNet(FLocation);
	//RemoveSelfEdges(signedNet);
	//PrintNStats("LoadEpinions", signedNet);
}

// Prints on the consol all kinds of distinct triads (signed or not, bidirected or not, include zero edge or not)
void PrintTriadEquiClasses() {
	bool Zero[] = {false, true}, BiDir[] = {false, true}, Signed[] = {false, true};
	for(bool b : BiDir)
		for(bool s : Signed)
			for(bool z : Zero)
			{
				cout << (s ? "Signed-" : "UnSigned-") << (b ? "BiDirection-" : "MonoDirection-") << (z ? "Zero Edge" : "No Zero Edge") << endl;
				TTriadEqClasH TriadGroups;
				int totalTriads = 0;
				TCtmsNet::GenTriadEquivClasses(TriadGroups, b, s, z);
				if (totalTriads == 0) {
					for (int j = 0; j < TriadGroups.Len(); j++) {totalTriads += TriadGroups[j].Val2;}}
				for (int j = 0; j < TriadGroups.Len(); j++) {
					std::cout << j+1 << ".\t" << TriadGroups[j].Val2 << "\t" << (double)TriadGroups[j].Val2/totalTriads << "\t" << TriadGroups[j].Val1->GetTriadStr().CStr() << std::endl;
				}
				std::cout << "Distinct Triads/Total Triads: " << TriadGroups.Len() << " / " << totalTriads << std::endl;
				getchar();
			}
}

void TraverseNodeSignedNetwork(PCtmsNet signedNet) {
	int countLine = 1;
	for(TSignNet::TNodeI NI = signedNet->BegNI(); NI < signedNet->EndNI(); NI++,countLine++) {
		cout << countLine << ".Node " << NI.GetId() << "   out-degree: " << NI.GetOutDeg() << "   in-degree: " << NI.GetInDeg() << endl;
	}
	cout << endl;
	cout << "NNNNNNNNNNNNNNNNNNNNNNN"  << endl;
}

void TraverseEdgeSignedNetwork(PCtmsNet signedNet) {
	for(TSignNet::TEdgeI EI = signedNet->BegEI(); EI < signedNet->EndEI(); EI++) {
		cout << " [" << EI.GetSrcNId() << ", " << EI.GetDstNId() << "]  " << EI.GetDat() << endl;
	}
	cout << endl;
	cout << "EEEEEEEEEEEEEEEEEEEEEEE"  << endl;
}

/*
Gets all Balanced or All positive and negative edges in the original network with
the minimum embeddedness of MinEmb
and maximum embeddedness of MaxEmb
but not necessarilly keep the MinEmb in the new network
priorities: 1.MinEmb/MaxEmb 2.Balance
*/
void GetPridictionEdges (const PCtmsNet& OrigNet, TIntTrV& outputEdges, const int MinEmb, const bool Balanced, const int MaxEmb = TInt::Mx) {
	TVec<TSignNet::TEdgeI> PosEVec, NegEVec;
	int pos = 0;

	for (TSignNet::TEdgeI EI = OrigNet->BegEI(); EI < OrigNet->EndEI(); EI++) {
		int nbrsCnt = TSnap::GetCmnNbrs(OrigNet, EI.GetSrcNId(), EI.GetDstNId());
		if (MinEmb <= nbrsCnt && nbrsCnt <= MaxEmb) {
			if (!Balanced) {
				outputEdges.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), EI()));
				if(EI() == 1) {pos++;}
			}
			else {
				if (EI() == 1) {PosEVec.Add(EI);}
				else {NegEVec.Add(EI);}
			}
		}
	}

	int minLen;
	if (Balanced) {
		PosEVec.Shuffle(TInt::Rnd);
		TInt::Rnd.Randomize();
		NegEVec.Shuffle(TInt::Rnd);
		if (NegEVec.Len() <= PosEVec.Len()) {minLen = NegEVec.Len();}
		else {minLen = PosEVec.Len();}
		for (int i = 0; i < minLen; i++) {
			outputEdges.Add(TIntTr(PosEVec[i].GetSrcNId(), PosEVec[i].GetDstNId(), PosEVec[i].GetDat()));
			outputEdges.Add(TIntTr(NegEVec[i].GetSrcNId(), NegEVec[i].GetDstNId(), NegEVec[i].GetDat()));
		}
	}

	printf("\n Total network edges: %d\n Prediction edges left: %d\n", OrigNet->GetEdges(), outputEdges.Len());
	printf(" Positive edge probability: %.4f\n", (Balanced? (double)minLen/outputEdges.Len() : (double)pos/outputEdges.Len()));
	return;
}

// runs algorithms for sign prediciton here. set algorithmsEnabled = 0 to run all algorithms
void runAllSignPredictionMethods(const PCtmsNet& Network, TIntTrV& Edges, const TStr outputFNm,
								 const bool SbNet, const int algorithmsEnabled = 0)
{
	if(algorithmsEnabled == 1 || algorithmsEnabled == 0) {//naive method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Naive method\n" << endl;
		TSignPredictionNaive naive(Network, outputFNm + "Naive.", Edges);
		naive.GetOutDegBasedAcc();
		naive.GetInDegBasedAcc();
		naive.GetweightedMeanBasedAcc();
		naive.GetcombinedDegBasedAcc();
		naive.GetweightedCombDegBasedAcc();
	}

	if(algorithmsEnabled == 2 || algorithmsEnabled == 0) {//run Abtin's method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Abt method (No learning; 96 triad classes)\n" << endl;
		TSignPredicNoLrn NoLearn(Network, outputFNm + "Abt.", Edges);
		NoLearn.CrossValidTest(SbNet);

	}

	if(algorithmsEnabled == 3 || algorithmsEnabled == 0) {//run Lescoves's method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Les method (16 triad logistic regression)\n" << endl;
		TSignPredictionLearn Learn(Network, outputFNm + "Les.", Edges);
		Learn.CrossValidTest();
	}

	if(algorithmsEnabled == 4 || algorithmsEnabled == 0 || algorithmsEnabled == 6) {//run Heuristic Balance method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Heuristic Balance\n" << endl;
		TEdgeSignPred heuristicPredition;
		for (TSignNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) { 
			heuristicPredition.AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat()); }
		for (int i = 0; i < Edges.Len(); i++)
			heuristicPredition.PredictBalance(Edges[i].Val1, Edges[i].Val2, Edges[i].Val3);

		heuristicPredition.PrintRes();
		const double TP = heuristicPredition.ResH.GetDat("Balance").Val1;
		const double FN = heuristicPredition.ResH.GetDat("Balance").Val2;
		const double FP = heuristicPredition.ResH.GetDat("Balance").Val3;
		const double TN = heuristicPredition.ResH.GetDat("Balance").Val4;

		const double acc = ((double)(TP+TN)/(TP+FP+TN+FN));
		const double tpr = (double)TP/(TP+FN);
		const double tnr = (double)TN/(TN+FP);
		const double acc2 = (tpr+tnr)/2.0;
		printf("\n_________________________________\n\n");
		printf("Accuracy: %.4f\nACC2: %.4f\nTPR: %.4f TNR: %.4f\n", acc, acc2, tpr, tnr);
	}

	if(algorithmsEnabled == 5 || algorithmsEnabled == 0 || algorithmsEnabled == 6) {//run Heuristic Status method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Heuristic Status\n" << endl;
		TEdgeSignPred heuristicPredition;
		for (TSignNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
			heuristicPredition.AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat());
		}
		for (int i = 0; i < Edges.Len(); i++)
			heuristicPredition.PredictStatus(Edges[i].Val1, Edges[i].Val2, Edges[i].Val3);

		heuristicPredition.PrintRes();
		const TStr met[] = {"SrcStat", "DstStat", "StatDif"};
		for (int i = 0; i < 3; i++) {
			const double TP = heuristicPredition.ResH.GetDat(met[i]).Val1;
			const double FN = heuristicPredition.ResH.GetDat(met[i]).Val2;
			const double FP = heuristicPredition.ResH.GetDat(met[i]).Val3;
			const double TN = heuristicPredition.ResH.GetDat(met[i]).Val4;

			const double acc = ((double)(TP+TN)/(TP+FP+TN+FN));
			const double tpr = (double)TP/(TP+FN);
			const double tnr = (double)TN/(TN+FP);
			const double acc2 = (tpr+tnr)/2.0;
			printf("\n%s _________________________________\n\n", met[i].CStr());
			printf("Accuracy: %.4f\nACC2: %.4f\nTPR: %.4f TNR: %.4f\n", acc, acc2, tpr, tnr);
		}		
	}
	return;
}

// preprocess network and analyze it
void analyzeNet(const char FileLocation[], const TStr outputFNmPrefix,
					   const int MinEm, const int Bal, const int Simp, const bool SbNet, const int runAlgorithms)
{
	PCtmsNet Net;
	TIntTrV predictionEdges;

	cout << "Loading network from file: " << FileLocation << " ...";
	LoadSNet(Net, FileLocation);
	cout << " DONE" << endl;

	if (Simp == 1) {
		cout << "Simplifying the network. (Removing bi-dir edges) ...";
		Net->SimplifyNet();
		cout << " DONE" << endl;
	}

	cout << "Getting sub-network" << " ...";
	GetPridictionEdges(Net, predictionEdges, MinEm, (Bal==0 ? false : true));
	cout << "DONE" << endl;

	runAllSignPredictionMethods(Net, predictionEdges, outputFNmPrefix, SbNet, runAlgorithms);
	return;
}

// preprocess network and analyze it. OLD STRICT version
void analyzeNetStrictOld(const char FileLocation[], const TStr outputFNmPrefix,
					   const int MinEm, const int Bal, const int Simp, const bool SbNet, const int runAlgorithms)
{
	PCtmsNet Net, MinEmNet;
	TIntTrV predictionEdges;

	cout << "Loading network from file: " << FileLocation << " ...";
	LoadSNet(Net, FileLocation);
	cout << " DONE" << endl;

	if (Simp == 1) {
		cout << "Simplifying the network. (Removing bi-dir edges) ...";
		Net->SimplifyNet();
		cout << " DONE" << endl;
	}

	cout << "Generating a subnetwork with minimum embeddedness:" << MinEm << " ...";
	MinEmNet = Net->GetMinEmbeddedSubNet(MinEm);
	cout << " DONE" << endl;
	if (Bal == 1) {
		cout << "Generating a subnetwork with balanced (50% positive) edges ...";
		Net = MinEmNet->GetBalSgnProbSubNet();
		cout << " DONE" << endl;
	}
	else {
		Net = MinEmNet;
	}

	int pos = 0;
	//Adding all new net edges to predictionEdges
	for(TSignNet::TEdgeI EI = Net->BegEI(); EI < Net->EndEI(); EI++) {
		predictionEdges.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), EI()));
		if(EI() == 1) {pos++;}
	}
	printf("\nPositive edge probability: %.4f\n", (double)pos/Net->GetEdges());

	runAllSignPredictionMethods(Net, predictionEdges, outputFNmPrefix, SbNet, runAlgorithms);
	return;
}

//Epinions entrance
void EpinionsAnalysis(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
					  const bool calcThetaSub = true, const int runAlgorithms = 0) {	
	cout << "*****  V." << DATE << " - EPINIONS SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr OutFNm = TStr::Fmt("%sEpinions.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));
	
	analyzeNet(FileNEpinions, OutFNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

void EpinionsAnalysisOld(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
						 const bool calcThetaSub = true, const int runAlgorithms = 0) {
	cout << "*****  V.Strict - EPINIONS SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr NtNm = TStr::Fmt("%sStrictEpinions.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));

	analyzeNetStrictOld(FileNEpinions, NtNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

//Slashdot entrance
void SlashdotAnalysis(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
					  const bool calcThetaSub = true, const int runAlgorithms = 0) {
	cout << "*****  V." << DATE << " - SLASHDOT SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr NtNm = TStr::Fmt("%sSlashdot.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));

	analyzeNet(FileNSlashdot, NtNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

void SlashdotAnalysisOld(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
						 const bool calcThetaSub = true, const int runAlgorithms = 0) {
	cout << "*****  V.Strict - SLASHDOT SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr NtNm = TStr::Fmt("%sStrictSlashdot.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));

	analyzeNetStrictOld(FileNSlashdot, NtNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

//Wikipedia entrance
void WikipediaAnalysis(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
					   const bool calcThetaSub = true, const int runAlgorithms = 0) {
	cout << "*****  V." << DATE << " - WIKIPEDIA SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr NtNm = TStr::Fmt("%sWikipedia.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));

	analyzeNet(FileNWikipedia, NtNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

void WikipediaAnalysisOld(const int MinEm = 25, const int Bal = 1, const int Simp = 0,
						  const bool calcThetaSub = true, const int runAlgorithms = 0) {
	cout << "*****  V.Strict - WIKIPEDIA SIGN PREDICTION (MinEm=" << MinEm << " Balanced:"  << 
		((Bal==1)? "Y":"N") << (Simp==1? " Simplified":"") << ")  *****\n" << endl;
	const TStr NtNm = TStr::Fmt("%sStrictWikipedia.%d.%s.",(Simp==1? "_smpl":""), MinEm, ((Bal==1)? "BL":""));

	analyzeNetStrictOld(FileNWikipedia, NtNm, MinEm, Bal, Simp, calcThetaSub, runAlgorithms);
	return;
}

void TestAnalysis() {
	PCtmsNet testNet, permutedNet, tempNet, rewiredNet, rewiredNet2, balancedNet;
	
	LoadSNet(testNet, FileNTest);
	PrintNStats("Test", testNet);
	balancedNet = testNet->GetBalSgnProbSubNet();
	TSignPredictionLearn Tpredictor(balancedNet, "BalTestSgn");
	Tpredictor.FeatureSelection(3);
	Tpredictor.SaveDataset("F");
	Tpredictor.SaveYXV("F");

	//Tpredictor.CrossValidTest(false, false);
	/*
	//PrintNStats("Test Network", testNet);
	permutedNet = testNet->PermuteEdgeSignsStrict();
	testNet->CompareNetWith(permutedNet);
	for (TSignNet::TEdgeI EI = testNet->BegEI(); EI < testNet->EndEI(); EI++) {
		cout << "Prob(" << EI.GetSrcNId() << ", " << EI.GetDstNId() << "):" <<
			testNet->GetEdgePosSgnProb(EI.GetSrcNId(), EI.GetDstNId()) << endl;
	}
	*/
	//testNet->CountSignedTriads2("../Test Triads/Test");
	//testNet->CountSignedTriads("Test");
	//testNet->CountSignedTriads3("TestNet");
	/*
	permutedNet = testNet->PermuteEdgeSigns2();
	PrintNStats("Test Network Permuted", permutedNet);
	permutedNet->CountSignedTriads2("..\\Test Signed Triads Modified2\\PermutedTestNet-M2");
	rewiredNet = testNet->RewireNetwork2();
	PrintNStats("Test Network Permuted", rewiredNet);
	rewiredNet->CountSignedTriads2("..\\Test Signed Triads Modified2\\RewiredTestNet-M2");
	*/
}

PCtmsNet GetNet(const TIntTrV& EL) {
	PCtmsNet net = TCtmsNet::New();
	for (int i = 0; i < EL.Len(); i++) {
		const TIntTr CurE = EL[i];
		if (!net->IsNode(CurE.Val1)) {net->AddNode(CurE.Val1);}
		if (!net->IsNode(CurE.Val2)) {net->AddNode(CurE.Val2);}
		net->AddEdge(CurE.Val1, CurE.Val2, CurE.Val3);
	}
	return net;
}

void CrossNetPredictionNOLrn(const char NN, const int MinEm = 1, const bool Bl = false, const bool useNewAlg = false) {
	const TStr info = TStr::Fmt("%d.%s.Abt", MinEm, (Bl? "Bl": ""));
	const TStr EpName = (TStr::Fmt("EpCross") + info), SlName = (TStr::Fmt("SlCross") + info),
		WiName = (TStr::Fmt("WiCross") + info), TstName = (TStr::Fmt("TstCross") + info);
	static PCtmsNet epinions, slashdot, wikipedia;
	static TIntTrV EpSelecEdgs, SlSelecEdgs, WiSelecEdgs, TsSelecEdgs;
	PCtmsNet testNet;
	cout << "***********  " << DATE << "  ABT" << info.CStr() << "CROSS NETWORK SIGN PREDICTION    ***********\n" << endl;
	cout << "Loading TestNet... ";	LoadSNet(testNet, FileNTest); cout << "DONE" << endl;
	if(epinions.Empty())
		cout << "Loading Epinions... ";	LoadSNet(epinions, FileNEpinions); cout << "DONE" << endl;
	if(slashdot.Empty())
		cout << "Loading Slashdot... ";	LoadSNet(slashdot, FileNSlashdot); cout << "DONE" << endl;
	if(wikipedia.Empty())
		cout << "Loading Wikipedia... "; LoadSNet(wikipedia, FileNWikipedia); cout << "DONE" << endl;
	
	if (EpSelecEdgs.Empty())
		GetPridictionEdges(epinions, EpSelecEdgs, MinEm, Bl);
	if (SlSelecEdgs.Empty())
		GetPridictionEdges(slashdot, SlSelecEdgs, MinEm, Bl);
	if (WiSelecEdgs.Empty())
		GetPridictionEdges(wikipedia, WiSelecEdgs, MinEm, Bl);
	GetPridictionEdges(testNet, TsSelecEdgs, MinEm, Bl);

	static TSignPredicNoLrn EpSPL(epinions, EpName.CStr(), EpSelecEdgs),
		SlSPL(slashdot, SlName, SlSelecEdgs),
		WiSPL(wikipedia, WiName, WiSelecEdgs);
	TSignPredicNoLrn TstSPL(testNet, TstName, TsSelecEdgs);
	IAssert(!epinions.Empty()); IAssert(!slashdot.Empty()); IAssert(!wikipedia.Empty());
	THash<TChA, TFlt> theTa;
	switch (NN)
	{
	case 'e':
		cout << "Setting Theta on Epinions...\n";
		if (useNewAlg)
			EpSPL.GetTheta(theTa);
		else
			TSignPredicNoLrn::GetTheta(GetNet(EpSelecEdgs), theTa, "EpinionsSubNet");
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 's':
		cout << "Setting Theta on Slashdot...\n";
		if (useNewAlg)
			SlSPL.GetTheta(theTa);
		else
			TSignPredicNoLrn::GetTheta(GetNet(SlSelecEdgs), theTa, "SlashdotSubNet");
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 'w':
		cout << "Setting Theta on Wikipedia...\n";
		if (useNewAlg)
			WiSPL.GetTheta(theTa);
		else
			TSignPredicNoLrn::GetTheta(GetNet(WiSelecEdgs), theTa, "WikipediaSubNet");
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 't':
		cout << "Setting Theta on TestNet...\n";
		if (useNewAlg)
			TstSPL.GetTheta(theTa);
		else
			TSignPredicNoLrn::GetTheta(GetNet(TsSelecEdgs), theTa, "TestnetSubNet");
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
	default:
		cout << "Invalid Network" << endl;
		break;
	}
	return;
}

void CrossNetPredictionLearn(const char NN, const int MinEm = 1, const bool Bl = false) {
	const TStr info = TStr::Fmt("%d.%s.Les", MinEm, (Bl? "Bl": ""));
	const TStr EpName = (TStr::Fmt("EpCross") + info), SlName = (TStr::Fmt("SlCross") + info),
		WiName = (TStr::Fmt("WiCross") + info), TstName = (TStr::Fmt("TstCross") + info);
	static PCtmsNet epinions, slashdot, wikipedia;
	static TIntTrV EpSelecEdgs, SlSelecEdgs, WiSelecEdgs, TsSelecEdgs;
	PCtmsNet testNet;
	cout << "***********  " << DATE << "  LES" << info.CStr() << "CROSS NETWORK SIGN PREDICTION    ***********\n" << endl;
	cout << "Loading TestNet... ";	LoadSNet(testNet, FileNTest); cout << "DONE" << endl;
	if(epinions.Empty())
		cout << "Loading Epinions... ";	LoadSNet(epinions, FileNEpinions); cout << "DONE" << endl;
	if(slashdot.Empty())
		cout << "Loading Slashdot... ";	LoadSNet(slashdot, FileNSlashdot); cout << "DONE" << endl;
	if(wikipedia.Empty())
		cout << "Loading Wikipedia... "; LoadSNet(wikipedia, FileNWikipedia); cout << "DONE" << endl;
	
	if (EpSelecEdgs.Empty())
		GetPridictionEdges(epinions, EpSelecEdgs, MinEm, Bl);
	if (SlSelecEdgs.Empty())
		GetPridictionEdges(slashdot, SlSelecEdgs, MinEm, Bl);
	if (WiSelecEdgs.Empty())
		GetPridictionEdges(wikipedia, WiSelecEdgs, MinEm, Bl);
	GetPridictionEdges(testNet, TsSelecEdgs, MinEm, Bl);

	static TSignPredictionLearn EpSPL(epinions, EpName.CStr(), EpSelecEdgs),
		SlSPL(slashdot, SlName, SlSelecEdgs),
		WiSPL(wikipedia, WiName, WiSelecEdgs);
	TSignPredictionLearn TstSPL(testNet, TstName, TsSelecEdgs);
	IAssert(!epinions.Empty()); IAssert(!slashdot.Empty()); IAssert(!wikipedia.Empty());
	TFltV theTa;
	switch (NN)
	{
	case 'e':
		cout << "Setting Theta on Epinions...\n";
		EpSPL.CalcTheta(theTa, true);
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 's':
		cout << "Setting Theta on Slashdot...\n";
		SlSPL.CalcTheta(theTa, true);
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 'w':
		cout << "Setting Theta on Wikipedia...\n";
		WiSPL.CalcTheta(theTa, true);
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
		break;
	case 't':
		cout << "Setting Theta on TestNet...\n";
		TstSPL.CalcTheta(theTa, true);
		cout << "DONE" << endl;
		cout << "Accuracy on Epinions: " <<	EpSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Slashdot: " << SlSPL.Test(theTa).accuracy << endl;
		cout << "Accuracy on Wikipedia: " << WiSPL.Test(theTa).accuracy << endl;
	default:
		cout << "Invalid Network" << endl;
		break;
	}
	return;
}

void crsbNt (int MinEm = 1, int nodes = 10) {
	PCtmsNet tstnt, abstrNt; LoadSNet(tstnt, FileNEpinions);
	tstnt->SortNIdById();
	abstrNt = TCtmsNet::New();
	TSignNet::TNodeI NI = tstnt->BegNI();
	for (int i = 0; i < nodes; i++) {
		if (!abstrNt->IsNode(NI.GetId())) {abstrNt->AddNode(NI);}
		for (int j = 0; j < NI.GetInDeg(); j++) {
			const int id = NI.GetInNId(j);
			if (!abstrNt->IsNode(id)) {abstrNt->AddNode(id);}
			if (!abstrNt->IsEdge(id, NI.GetId())) {abstrNt->AddEdge(id, NI.GetId(), tstnt->GetEDat(id, NI.GetId()));}
		}
		for (int j = 0; j < NI.GetOutDeg(); j++) {
			const int id = NI.GetOutNId(j);
			if (!abstrNt->IsNode(id)) {abstrNt->AddNode(id);}
			if (!abstrNt->IsEdge(NI.GetId(), id)) {abstrNt->AddEdge(NI.GetId(), id, tstnt->GetEDat(NI.GetId(), id));}
		}
		NI++;
	}
	abstrNt->SaveTxt("Net.txt");
	PCtmsNet subNet = abstrNt->GetMinEmbeddedSubNet(MinEm);
	subNet->SaveTxt("SubNet.txt");
	TIntTrV Edges; GetPridictionEdges(abstrNt, Edges, MinEm, false);
	PCtmsNet edgeNet = TCtmsNet::New();
	for (int i = 0;  i < Edges.Len(); i++) {
		if (!edgeNet->IsNode(Edges[i].Val1)) {edgeNet->AddNode(Edges[i].Val1);}
		if (!edgeNet->IsNode(Edges[i].Val2)) {edgeNet->AddNode(Edges[i].Val2);}
		if (!edgeNet->IsEdge(Edges[i].Val1, Edges[i].Val2)) {edgeNet->AddEdge(Edges[i].Val1, Edges[i].Val2, Edges[i].Val3);}
	}
	edgeNet->SaveTxt("EdgeNet.txt");
}

int main(int argc, char* argv[]) {
	FILE *outputStream;
	//DefaultConstructor();
	//PrintTriadEquiClasses();
	do{		
		char n;
		int cross;
		cout << "Same net cross validation(0) OR Cross network testing(1)? ";
		cin >> cross; cin.ignore(INT_MAX, '\n');
		cout << "(e)pinions or (s)lashdot or (w)ikipedia or (t)est ? ";
		cin >> n; cin.ignore(INT_MAX, '\n');
		if (cross == 1) {
			static int MinEm = -1, useLearn = 0;
			static bool Bl = false;
			if (MinEm == -1) {
				cout << "MinEm, BL(1), Learn(1): ";
				cin >> MinEm >> Bl >> useLearn;
				cin.ignore(INT_MAX, '\n');
			}
			
			if (ConsoleToFile) {
				const TStr logFile = TStr::Fmt("%s%s.%d.%s_CrossOn(%c).Log", LogDir,
					(useLearn?"Les":"Abt"), MinEm, (Bl?"Bl":""), n);
			if((outputStream = freopen(logFile.CStr(), "w", stdout)) == NULL)
				exit(-1);
			}

			if (useLearn == 1) {
				CrossNetPredictionLearn(n, MinEm, Bl);}
			else {
				CrossNetPredictionNOLrn(n, MinEm, Bl, true);}
		}
		else
		{
			int MnEm, Bl, Smpl;
			cout << "MinEm, BL(1), Simp(1): ";
			cin >> MnEm >> Bl >> Smpl;
			cin.ignore(INT_MAX, '\n');

			int strctVersion = 0;
			/*cout << "Paper version(0) OR Strict version(1)? ";
			cin >> strctVersion;
			cin.ignore(INT_MAX, '\n');*/

			bool calcThetaSub = false;
			cout << "Calc theta using new subnet (1)? ";
			cin >> calcThetaSub;
			cin.ignore(INT_MAX, '\n');

			int algoRun;
			cout << "ALGORITHMS - All(0) NAIVE(1) ABT(2) LES(3) BALANCE(4) STATUS(5) 4&5(6): ";
			cin >> algoRun; cin.ignore(INT_MAX, '\n');

			if (ConsoleToFile) {
			//redirecting console output to file
			TStr logFile = TStr::Fmt(LogDir);
			TStr algo;
			if(algoRun == 0) algo = "ALL";
			else if(algoRun == 1) algo = "NIV";
			else if(algoRun == 2) algo = "ABT";
			else if(algoRun == 3) algo = "LES";
			else if(algoRun == 4) algo = "BAL";
			else if(algoRun == 5) algo = "STS";
			switch (n)
			{
			case 'e':
				logFile += TStr::Fmt("%sEpinions.%d.%s_%s.Log",(Smpl==1? "_smpl":""), MnEm, ((Bl==1)? "BL":""), algo.CStr());
				break;
			case 's':
				logFile += TStr::Fmt("%sSlashdot.%d.%s_%s.Log",(Smpl==1? "_smpl":""), MnEm, ((Bl==1)? "BL":""), algo.CStr());
				break;
			case 'w':
				logFile += TStr::Fmt("%sWikipedia.%d.%s_%s.Log",(Smpl==1? "_smpl":""), MnEm, ((Bl==1)? "BL":""), algo.CStr());
				break;
			case 't':
				logFile += TStr::Fmt("%sTest.%d.%s_%s.Log",(Smpl==1? "_smpl":""), MnEm, ((Bl==1)? "BL":""), algo.CStr());
				break;
			default:
				logFile += TStr::Fmt("Log.Log");
				break;
			}
			if((outputStream = freopen(logFile.CStr(), "w", stdout)) == NULL)
				exit(-1);
			}

			//start analysis
			switch(n) {
			case 'e':
				(strctVersion == 1) ? EpinionsAnalysisOld(MnEm, Bl, Smpl, calcThetaSub, algoRun) :
					EpinionsAnalysis(MnEm, Bl, Smpl, calcThetaSub, algoRun);
				break;
			case 's':
				(strctVersion == 1) ? SlashdotAnalysisOld(MnEm, Bl, Smpl, calcThetaSub, algoRun) :
					SlashdotAnalysis(MnEm, Bl, Smpl, calcThetaSub, algoRun);
				break;
			case 'w':
				(strctVersion == 1) ? WikipediaAnalysisOld(MnEm, Bl, Smpl, calcThetaSub, algoRun) :
					WikipediaAnalysis(MnEm, Bl, Smpl, calcThetaSub, algoRun);
				break;
			case 't':
				TestAnalysis();
				break;
			default: cout << "Invalid network!" << endl;
			}
		}
		
		if (ConsoleToFile) 
		{
		//redirecting file output to console
		fclose(outputStream);
		outputStream = freopen("CON", "w", stdout);
		}

		cout << "____________________________________________________________________________" << endl;
		cout << "Continue? ";		
	} while(getchar() != 'n' && getchar() == '\n');
	
	return 0;
}