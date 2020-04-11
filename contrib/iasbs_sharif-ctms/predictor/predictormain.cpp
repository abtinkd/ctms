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
		cout << "arg1 -- algorithm: NAIVE(1) BALANCE(2) STATUS(3) LogReg(4) CTMS(5)" << endl;
		cout << "arg2 -- input network file path" << endl;
		cout << "arg3 -- input unsigned edges file path" << endl;
		cout << "arg4 -- output predicted edges file path" << endl;		
		return 1;
	}

	const int alg = atoi(argv[1]);
	const char* networkFilePath = argv[2];
	const char* edgesFilePath = argv[3];
	const char* outputFilePath = argv[4];

	cout << "network file: " << networkFilePath << endl;
	cout << "unsigned edges file: " << edgesFilePath << endl;
	cout << "loading network ..." << endl;
	PCtmsNet network = TCtmsNet::LoadSignedNet(networkFilePath);
	cout << "loading edges ..." << endl;
	TIntPrV edges;
	LoadEdges(edges, edgesFilePath);
	
	TChA algorithm;
	TSignPredictor* predictor;
	switch (alg) {
	case 1: algorithm = "naive"; 
		predictor = new TNaivePredictor(network); break;
	case 2: algorithm = "balance";
		predictor = new TBalanceBasedPredictor(network); break;
	case 3: algorithm = "status";
		predictor = new TStatusBasedPredictor(network); break;	
	case 4: algorithm = "logreg"; 
		predictor = new TLogisticRegression(network); break;
	case 5: algorithm = "ctms";
		predictor = new TCTMSProbabilisticInference(network); break;
	default:
		predictor = new TNaivePredictor(network);
	}
	cout << "building model using \'" << algorithm.CStr() << "\' algorithm..." << endl;
	predictor->build();	

	ofstream outputFile;
	outputFile.open(outputFilePath);
	outputFile << "# network file: " << networkFilePath << endl;
	outputFile << "# edges file: " << edgesFilePath << endl;
	outputFile << "# FromNodeId ToNodeId	Sign" << endl;
	for (int i = 0; i < edges.Len(); i++) {
		const int sign = predictor->predict(edges[i].Val1, edges[i].Val2);
		outputFile << edges[i].Val1 << "\t" << edges[i].Val2 << "\t" << sign << endl;
	}
	outputFile.close();
}