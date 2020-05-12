#include <iostream>
#include <string>
#include <fstream>
#include <ctmsnet.h>

using namespace std;

/*
Gets sampleFraction (balanced number ) of positive and negative edges
in the original network with
the minimum embeddedness of MinEmb
and maximum embeddedness of MaxEmb
but not necessarilly keep the MinEmb in the new network
priorities: 1.MinEmb/MaxEmb 2.Balance
*/
void GetMaskEdges(const PCtmsNet& OrigNet, TIntTrV& outputEdges, 
	const double sampleFraction, const int MinEmb, const bool Balanced, const int MaxEmb) {
	TVec<TSignNet::TEdgeI> PosEVec, NegEVec;

	for (TSignNet::TEdgeI EI = OrigNet->BegEI(); EI < OrigNet->EndEI(); EI++) {
		int nbrsCnt = TSnap::GetCmnNbrs(OrigNet, EI.GetSrcNId(), EI.GetDstNId());
		if (MinEmb <= nbrsCnt && nbrsCnt <= MaxEmb) {			
			if (EI() == 1) { PosEVec.Add(EI); }
			else { NegEVec.Add(EI); }
		}
	}
	
	int addCountPos, addCountNeg;
	if (Balanced) {
		const int minLen = NegEVec.Len() <= PosEVec.Len() ? NegEVec.Len() : PosEVec.Len();
		addCountPos = addCountNeg = minLen * sampleFraction;
	}
	else {
		addCountPos = PosEVec.Len() * sampleFraction;
		addCountNeg = NegEVec.Len() * sampleFraction;
	}
	TInt::Rnd.Randomize();
	PosEVec.Shuffle(TInt::Rnd);
	for (int i = 0; i < addCountPos; i++) {
		outputEdges.Add(TIntTr(PosEVec[i].GetSrcNId(), PosEVec[i].GetDstNId(), PosEVec[i].GetDat()));
	}
	TInt::Rnd.Randomize();
	NegEVec.Shuffle(TInt::Rnd);
	for (int i = 0; i < addCountNeg; i++) {
		outputEdges.Add(TIntTr(NegEVec[i].GetSrcNId(), NegEVec[i].GetDstNId(), NegEVec[i].GetDat()));
	}
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Following arguments are needed to extract the *test edges* from the network:" << endl;
		cout << "arg1 -- source file path" << endl;
		cout << "arg2 -- output file name" << endl;
		cout << "arg3 -- minimum embeddedness" << endl;
		cout << "arg4 -- balanced or not?(1/0) [optional]" << endl;
		cout << "arg5 -- 0 < sample size <= 100 [optional]" << endl;
		cout << "arg6 -- maximum embeddedness [optional]" << endl;
		return 1;
	}

	const string sourceFilePath = argv[1];
	const string outputFileName = argv[2];
	const int minEm = atoi(argv[3]);
	const bool isBalanced = (argc < 5 || atoi(argv[4]) == 0 ? false : true);
	const int sampleRatio = 
		(argc < 6 || atoi(argv[5]) <= 0 || atoi(argv[5]) > 100 ? 100 : atoi(argv[5]));
	const TInt maxEm = (argc < 7 ? TInt::Mx : atoi(argv[5]));
	
	const double sampleFraction = sampleRatio / 100.0;
	cout << "source file: " << sourceFilePath << endl;	
	cout << "sample ratio: " << sampleFraction << endl;
	cout << "minimum embeddedness: " << minEm << endl;
	cout << "maximum embeddedness: " << maxEm << endl;
	cout << "is balanced? " << (isBalanced ? "yes" : "no") << endl;

	PCtmsNet sourceNet;	
	cout << "loading source-net ..." << endl;
	sourceNet = TCtmsNet::LoadSignedNet(sourceFilePath.data());

	cout << "extracting sub-net ..." << endl;
	TIntTrV maskE;
	GetMaskEdges(sourceNet, maskE, sampleFraction, minEm, isBalanced, maxEm);

	const long sourceEcount = sourceNet->GetEdges();
	const long subnetEcount = maskE.Len();
	long pos = 0;
	for (int i = 0; i < maskE.Len(); i++) {
		maskE[i].GetVal3() == 1 ? ++pos : 1;
	}
	
	const string outputTrainFilePath = "out/train_" + outputFileName;
	const string outputLabelsFilePath = "out/lables_" + outputFileName;
	const string outputTestFilePath = "out/test_" + outputFileName;
	cout << "saving ..." << endl;
	cout << outputTrainFilePath << endl;
	cout << outputLabelsFilePath << endl;
	cout << outputTestFilePath << endl;
	
	ofstream outputTrainFile, outputLabelsFile, outputTestFile;
	outputTrainFile.open(outputTrainFilePath);
	outputLabelsFile.open(outputLabelsFilePath);
	outputTestFile.open(outputTestFilePath);
	outputTrainFile << "# source net : " << sourceFilePath << endl;
	outputLabelsFile << "# source net : " << sourceFilePath << endl;
	outputTestFile << "# source net : " << sourceFilePath << endl;
	outputTrainFile << "# complement net: " << outputLabelsFilePath << endl;
	outputLabelsFile << "# complement net: " << outputTrainFilePath << endl;
	outputTestFile << "# complement net: " << outputTrainFilePath << endl;
	
	outputLabelsFile << "# minimum embeddedness: " << minEm << endl;
	outputTestFile << "# minimum embeddedness: " << minEm << endl;
	if (maxEm != TInt::Mx) {		
		outputLabelsFile << "# maximum embeddedness: " << maxEm << endl;
		outputTestFile << "# maximum embeddedness: " << maxEm << endl;
	}
	outputLabelsFile << "# is balanced? " << (isBalanced? 'Y' : 'N') << endl;
	outputTestFile << "# is balanced? " << (isBalanced ? 'Y' : 'N') << endl;
	outputLabelsFile << "# sampling rate: " << sampleFraction << endl;
	outputTestFile << "# sampling rate: " << sampleFraction << endl;
	outputLabelsFile << "# edge count: " << subnetEcount << endl;
	outputTestFile << "# edge count: " << subnetEcount << endl;
	outputLabelsFile << "# positive edge ratio: " << 1.0 * pos / subnetEcount << endl;
	outputTestFile << "# positive edge ratio: " << 1.0 * pos / subnetEcount << endl;
	outputLabelsFile << "# FromNodeId ToNodeId	Sign" << endl;
	outputTestFile << "# FromNodeId ToNodeId" << endl;
	
	for (int i = 0; i < maskE.Len(); i++) {
		sourceNet->DelEdge(maskE[i].GetVal1(), maskE[i].GetVal2());
		outputLabelsFile << maskE[i].GetVal1() << "\t" << maskE[i].GetVal2() <<
			"\t" << maskE[i].GetVal3() << endl;
		outputTestFile << maskE[i].GetVal1() << "\t" << maskE[i].GetVal2() << endl;
	}
	outputLabelsFile.close();
	outputTestFile.close();

	for (TSignNet::TEdgeI EI = sourceNet->BegEI(); EI < sourceNet->EndEI(); EI++) {
		outputTrainFile << EI.GetSrcNId() << "\t" << EI.GetDstNId() <<
			"\t" << EI.GetDat() << endl;
	}
	outputTrainFile.close();	
}