#include <iostream>
#include <fstream>
#include <ctmsnet.h>

using namespace std;

/*
Gets all Balanced or All positive and negative edges in the original network with
the minimum embeddedness of MinEmb
and maximum embeddedness of MaxEmb
but not necessarilly keep the MinEmb in the new network
priorities: 1.MinEmb/MaxEmb 2.Balance
*/
void GetMaskEdges(const PCtmsNet& OrigNet, TIntTrV& outputEdges, const int MinEmb, const bool Balanced, const int MaxEmb) {
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

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "missing arguments:" << endl;
		cout << "arg1 -- source file path" << endl;
		cout << "arg2 -- output file path" << endl;
		cout << "arg3 -- minimum embeddedness" << endl;
		cout << "arg4 -- balanced or not?(1/0)" << endl;		
		return 1;
	}

	const char* sourceFilePath = argv[1];
	const char* outputFilePath = argv[2];
	const int minEm = atoi(argv[3]);
	const bool isBalanced = (argc < 5 || atoi(argv[4]) == 0 ? false : true);
	const TInt maxEm = (argc < 6 ? TInt::Mx : atoi(argv[5]));	

	cout << "source file: " << sourceFilePath << endl;
	cout << "output file: " << outputFilePath << endl;
	cout << "minimum embeddedness: " << minEm << endl;
	cout << "maximum embeddedness: " << maxEm << endl;
	cout << "is balanced? " << (isBalanced ? 'yes' : 'no') << endl;

	PCtmsNet sourceNet;	
	cout << "loading source-net ..." << endl;
	sourceNet = TCtmsNet::LoadSignedNet(sourceFilePath);

	cout << "extracting sub-net ..." << endl;
	TIntTrV maskE;
	GetMaskEdges(sourceNet, maskE, minEm, isBalanced, maxEm);	

	const long sourceEcount = sourceNet->GetEdges();
	const long subnetEcount = maskE.Len();
	long pos = 0;
	for (int i = 0; i < maskE.Len(); i++) {
		maskE[i].GetVal3() == 1 ? ++pos : 1;
	}	
	
	cout << "saving ..." << endl;
	ofstream outputFile;
	outputFile.open(outputFilePath);
	outputFile << "# source net : " << sourceFilePath << endl;
	outputFile << "# minimum embeddedness: " << minEm << endl;
	if (maxEm != TInt::Mx) {
		outputFile << "# maximum embeddedness: " << maxEm << endl;
	}
	outputFile << "# is balanced? " << (isBalanced? 'Y' : 'N') << endl;
	outputFile << "# edge count: " << subnetEcount << "/" << sourceEcount << endl;
	outputFile << "# positive edge count: " << pos << endl;
	outputFile << "# FromNodeId ToNodeId	Sign" << endl;
	for (int i = 0; i < maskE.Len(); i++) {
		outputFile << maskE[i].GetVal1() << "\t" << maskE[i].GetVal2() <<
			"\t" << maskE[i].GetVal3() << endl;
	}
	outputFile.close();
}