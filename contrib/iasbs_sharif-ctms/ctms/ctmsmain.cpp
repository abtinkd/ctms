#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "ctmsnet.h"
#include "sign_prediction.h"

using namespace std;


// Prints all kinds of distinct triads (signed/unsgined, mono/bidirected, with/without no edge sides)
void PrintTriadEquiClasses() {
	cout << "printing distribution of triad classes according to sign/direction/noEdge status..." << endl;
	bool Zero[] = { false, true }, BiDir[] = { false, true }, Signed[] = { false, true };
	for (bool b : BiDir)
		for (bool s : Signed)
			for (bool z : Zero)
			{
				cout << (s ? "Signed" : "UnSigned") << " - " <<
					(b ? "BiDirection" : "MonoDirection") << " - " <<
					(z ? "Zero Edge premitted" : "No Zero Edge permitted") << endl;
				cout << "#\tclass\t\tcount\tPr" << endl;
				TTriadEqClasH TriadGroups;
				int totalTriads = 0;
				TCtmsNet::GenTriadEquivClasses(TriadGroups, b, s, z);
				if (totalTriads == 0) {
					for (int j = 0; j < TriadGroups.Len(); j++) { totalTriads += TriadGroups[j].Val2; }
				}
				for (int j = 0; j < TriadGroups.Len(); j++) {
					cout << j + 1 << ".\t" << TriadGroups[j].Val1->GetTriadStr().CStr() <<
						"\t" << TriadGroups[j].Val2 << "\t" << (double)TriadGroups[j].Val2 / totalTriads << endl;
				}
				cout << "Distinct Triads/Total Triads: " << TriadGroups.Len() << " / " << totalTriads << endl;
				cout << "press Enter to continue..";
				getchar();
			}
}


/*
Gets all Balanced or All positive and negative edges in the original network with
the minimum embeddedness of MinEmb
and maximum embeddedness of MaxEmb
but not necessarilly keep the MinEmb in the new network
priorities: 1.MinEmb/MaxEmb 2.Balance
*/
void GetMaskEdges(const PCtmsNet& OrigNet, TIntTrV& outputEdges, const int MinEmb, const bool Balanced, const int MaxEmb = TInt::Mx) {
	TVec<TSignNet::TEdgeI> PosEVec, NegEVec;
	int pos = 0;

	for (TSignNet::TEdgeI EI = OrigNet->BegEI(); EI < OrigNet->EndEI(); EI++) {
		int nbrsCnt = TSnap::GetCmnNbrs(OrigNet, EI.GetSrcNId(), EI.GetDstNId());
		if (MinEmb <= nbrsCnt && nbrsCnt <= MaxEmb) {
			if (!Balanced) {
				outputEdges.Add(TIntTr(EI.GetSrcNId(), EI.GetDstNId(), EI()));
				if (EI() == 1) { pos++; }
			}
			else {
				if (EI() == 1) { PosEVec.Add(EI); }
				else { NegEVec.Add(EI); }
			}
		}
	}

	int minLen;
	if (Balanced) {
		PosEVec.Shuffle(TInt::Rnd);
		TInt::Rnd.Randomize();
		NegEVec.Shuffle(TInt::Rnd);
		if (NegEVec.Len() <= PosEVec.Len()) { minLen = NegEVec.Len(); }
		else { minLen = PosEVec.Len(); }
		for (int i = 0; i < minLen; i++) {
			outputEdges.Add(TIntTr(PosEVec[i].GetSrcNId(), PosEVec[i].GetDstNId(), PosEVec[i].GetDat()));
			outputEdges.Add(TIntTr(NegEVec[i].GetSrcNId(), NegEVec[i].GetDstNId(), NegEVec[i].GetDat()));
		}
	}	
	return;
}


void BuildNet(PCtmsNet& net, TIntTrV& maskE, const char FilePath[], const int MinEm, const bool isBalanced)
{	
	cout << "loading full network ..." << endl;
	net = TCtmsNet::LoadSignedNet(FilePath);		

	cout << "extracting sub network ..." << endl;
	GetMaskEdges(net, maskE, MinEm, isBalanced);
		
	cout << "  fullnet edge count: " << net->GetEdges() << endl;
	cout << "  sub-net edge count: " << maskE.Len() << endl;	

	int pos = 0;
	for (int i = 0; i < maskE.Len(); i++) {
		maskE[i].GetVal3() == 1 ? ++pos : 1;
	}	
	printf("  positive edge ratio:  %.2f\n", (double)pos / (double)maskE.Len());	
	
	return;
}


// runs algorithms for sign prediciton here. set algorithmsEnabled = 0 to run all algorithms
void runSignPredictionMethods(const PCtmsNet& Network, TIntTrV& Edges, const TStr outPrefix, const int algorithmsEnabled = 0)
{
	if (algorithmsEnabled == 1 || algorithmsEnabled == 0) {//naive method
		cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Naive method\n" << endl;
		TNaiveInference naive(Network, outPrefix + "naive.", Edges);
		naive.GetOutDegBasedAcc();
		naive.GetInDegBasedAcc();
		naive.GetweightedMeanBasedAcc();
		naive.GetcombinedDegBasedAcc();
		naive.GetweightedCombDegBasedAcc();
	}

	if (algorithmsEnabled == 2 || algorithmsEnabled == 0) {//run CTMS-based method
		cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: CTMS-based method (probabilistic inference, 96 CTMS classes)\n" << endl;
		TCTMSProbabilisticInference ctmsInference(Network, outPrefix, Edges, (outPrefix == "" ? false : true));
		ctmsInference.CrossValidTest(false);

	}

	if (algorithmsEnabled == 3 || algorithmsEnabled == 0) {//run LogReg-16triad method
		cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: LR method (logistic regression, 16 triad types)\n" << endl;
		TLogisticRegression lr16Inference(Network, outPrefix, Edges, (outPrefix == ""? false: true));
		lr16Inference.CrossValidTest();
	}

	if (algorithmsEnabled == 4 || algorithmsEnabled == 0) {//run Heuristic Balance method
		cout << "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Heuristic Balance\n" << endl;
		TEdgeSignPred heuristicPredition;
		heuristicPredition.Network = &*Network;
		/*for (TSignNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
			heuristicPredition.AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat()); }*/
		for (int i = 0; i < Edges.Len(); i++)
			heuristicPredition.PredictBalance(Edges[i].Val1, Edges[i].Val2, Edges[i].Val3);

		heuristicPredition.PrintRes();
		const double TP = heuristicPredition.ResH.GetDat("Balance").Val1;
		const double FN = heuristicPredition.ResH.GetDat("Balance").Val2;
		const double FP = heuristicPredition.ResH.GetDat("Balance").Val3;
		const double TN = heuristicPredition.ResH.GetDat("Balance").Val4;

		const double acc = ((double)(TP + TN) / (TP + FP + TN + FN));
		const double tpr = (double)TP / (TP + FN);
		const double tnr = (double)TN / (TN + FP);
		const double acc2 = (tpr + tnr) / 2.0;
		printf("\n_________________________________\n\n");
		printf("Acc: %.4f\nmean(TPR,TNR): %.4f\nTPR: %.4f TNR: %.4f\n", acc, acc2, tpr, tnr);
	}

	if (algorithmsEnabled == 5 || algorithmsEnabled == 0) {//run Heuristic Status method
		cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
		cout << "Sign prediction: Heuristic Status\n" << endl;
		TEdgeSignPred heuristicPredition;
		heuristicPredition.Network = &*Network;
		/*for (TSignNet::TEdgeI EI = Network->BegEI(); EI < Network->EndEI(); EI++) {
			heuristicPredition.AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat());
		}*/
		for (int i = 0; i < Edges.Len(); i++)
			heuristicPredition.PredictStatus(Edges[i].Val1, Edges[i].Val2, Edges[i].Val3);

		heuristicPredition.PrintRes();
		const TStr met[] = { "SrcStat", "DstStat", "StatDif" };
		for (int i = 0; i < 3; i++) {
			const double TP = heuristicPredition.ResH.GetDat(met[i]).Val1;
			const double FN = heuristicPredition.ResH.GetDat(met[i]).Val2;
			const double FP = heuristicPredition.ResH.GetDat(met[i]).Val3;
			const double TN = heuristicPredition.ResH.GetDat(met[i]).Val4;

			const double acc = ((double)(TP + TN) / (TP + FP + TN + FN));
			const double tpr = (double)TP / (TP + FN);
			const double tnr = (double)TN / (TN + FP);
			const double acc2 = (tpr + tnr) / 2.0;
			printf("\n%s _________________________________\n", met[i].CStr());
			printf("Acc: %.4f\nmean(TPR,TNR): %.4f\nTPR: %.4f \nTNR: %.4f\n", acc, acc2, tpr, tnr);
		}
	}
	return;
}


int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "Missing arguments for the main function. Instructions are in README file." << endl;
	}

	if (strcmp(argv[1], "EQC") == 0) {
		PrintTriadEquiClasses();
		return 0;
	}

	const char* filePath = argv[1];
	const int minEm = atoi(argv[2]);
	const bool isBalanced = (argc < 3 || atoi(argv[3]) == 0 ? false : true);

	cout << "input file: " << filePath << endl;
	cout << "minimum embeddedness: " << minEm << endl;
	cout << "limit to a random subnet of balanced edge signs: " << (isBalanced? "yes": "no") << endl;
	
	PCtmsNet net;
	TIntTrV maskE;
	BuildNet(net, maskE, filePath, minEm, isBalanced);
	
	do {
		
		int algoRun = -1;
		while (algoRun < 0 || algoRun > 5) {
			cout << "Methods: All(0) NAIVE(1) CTMS(2) LogReg(3) BALANCE(4) STATUS(5): ";
			cin >> algoRun; cin.ignore(INT_MAX, '\n');
		}

		string outPfx = "";
		cout << "Set a net name to enable logging in <./results>. (press Enter to disable it): ";
		getline(cin, outPfx);

		runSignPredictionMethods(net, maskE, outPfx.c_str(), algoRun);
		cout << "__________________________________________________________________________________\n";
		cout << "Use another method (y/n)? ";
	} while (getchar() != 'n' && getchar() == '\n');

	return 0;
}