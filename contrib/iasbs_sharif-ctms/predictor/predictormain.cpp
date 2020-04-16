#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <ctmsnet.h>
#include <sign_prediction2.h>

using namespace std;

void LoadEdges(TIntPrV& edges, const char* edgesFilePath) {
	ifstream inputFile(edgesFilePath);
	string line;
	if (inputFile.is_open()) {
		while (getline(inputFile, line)) {
			if (line.at(0) == '#') continue;
			stringstream ss(line);
			getline(ss, line, '\t');
			const TInt srcId = stoi(line);
			getline(ss, line, '\t');
			const TInt desId = stoi(line);
			edges.Add(TIntPr(srcId, desId));
		}
		inputFile.close();
	}
	else cout << edgesFilePath << " cannot be openned.\n";	
}

int main(int argc, char* argv[]) {
	if (argc < 5) {
		cout << "missing arguments:" << endl;
		cout << "arg1 -- algorithm:" << endl <<
			"   1 <- generative-based" << endl <<
			"   2 <- receptive-based" << endl <<
			"   3 <- compound-based" << endl <<
			"   4 <- weighted generative receptive combination" << endl <<
			"   5 <- heuristic balance" << endl <<
			"   6 <- heuristic status" << endl <<
			"   7 <- logistic regression" << endl <<
			"   8 <- CTMS" << endl <<
			"   9 <- Local CTMS (not included in the paper)" << endl;
		cout << "arg2 -- train file path" << endl;
		cout << "arg3 -- test file path" << endl;
		cout << "arg4 -- predictions file name" << endl;		
		return 1;
	}

	const int alg = atoi(argv[1]);
	const string networkFilePath = argv[2];
	const string edgesFilePath = argv[3];
	string outputFilePath = argv[4];

	cout << "train file: " << networkFilePath << endl;
	cout << "test file: " << edgesFilePath << endl;
	cout << "loading network ..." << endl;
	PCtmsNet network = TCtmsNet::LoadSignedNet(networkFilePath.data());
	cout << "loading edges ..." << endl;
	TIntPrV edges;
	LoadEdges(edges, edgesFilePath.data());
	
	string algorithm;
	TSignPredictor* predictor;
	switch (alg) {
	case 1: algorithm = "gnr";
		predictor = new TNaivePredictor(network, algorithm.data()); break;
	case 2: algorithm = "rcp";
		predictor = new TNaivePredictor(network, algorithm.data()); break;
	case 3: algorithm = "cmp";
		predictor = new TNaivePredictor(network, algorithm.data()); break;
	case 4: algorithm = "wgr";
		predictor = new TNaivePredictor(network, algorithm.data()); break;
	case 5: algorithm = "balance";
		predictor = new TBalanceBasedPredictor(network); break;
	case 6: algorithm = "status";
		predictor = new TStatusBasedPredictor(network); break;	
	case 7: algorithm = "logreg"; 
		predictor = new TLogisticRegression(network); break;
	case 8: algorithm = "ctms";
		predictor = new TCTMSProbabilisticInference(network); break;
	case 9: algorithm = "localctms";
		predictor = new TCTMSProbabilisticInferenceLocal(network); break;
	default:
		predictor = new TNaivePredictor(network, "rnd");
	}
	cout << "building model using \'" << algorithm << "\' algorithm..." << endl;
	predictor->build();	

	outputFilePath = "out/" + algorithm + "_" + outputFilePath;
	cout << "saving predictions in " << outputFilePath << endl;
	
	ofstream outputFile;
	outputFile.open(outputFilePath);
	outputFile << "# prediction results for: " << algorithm << endl;
	outputFile << "# network file: " << networkFilePath << endl;
	outputFile << "# edges file: " << edgesFilePath << endl;
	outputFile << "# FromNodeId ToNodeId	Sign" << endl;
	for (int i = 0; i < edges.Len(); i++) {
		if (!network->IsNode(edges[i].Val1)) continue;
		if (!network->IsNode(edges[i].Val2)) continue;
		const int sign = predictor->predict(edges[i].Val1, edges[i].Val2);
		outputFile << edges[i].Val1 << "\t" << edges[i].Val2 << "\t" << sign << endl;
	}
	outputFile.close();
}