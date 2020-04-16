#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>

using namespace std;

void read(const char* filePath, unordered_map<string, int>& listing) {
	ifstream file(filePath);
	if (file.is_open()) {
		string line;
		while (getline(file, line)) {
			if (line.at(0) == '#') continue;
			stringstream ss(line);
			string srcId, desId, sign;
			getline(ss, srcId, '\t');			
			getline(ss, desId, '\t');
			getline(ss, sign);			
			listing[srcId + " " + desId] = stoi(sign);
		}
	}
	else cout << "could not open file " << filePath << endl;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "missing arguments:" << endl;
		cout << "arg1 -- true lables file path" << endl;
		cout << "arg2 -- prediction file path" << endl;
		return 1;
	}

	const string filePath1 = argv[1];
	const string filePath2 = argv[2];

	unordered_map<string, int> actual;
	read(filePath1.data(), actual);
	unordered_map<string, int> inferred;
	read(filePath2.data(), inferred);
	
	int tp = 0, fp = 0, tn = 0, fn = 0, missing = 0, unknown = 0;
	for (auto entry : actual) {
		const string edge = entry.first;
		const int sign = entry.second;
		unordered_map<string, int>::const_iterator got = inferred.find(edge);
		if (got == inferred.end()) {
			missing++;
			continue;
		}
		if (sign == 1) {
			if (got->second == sign) tp++;
			else fn++;
		}
		else if (sign == -1) {
			if (got->second == sign) tn++;
			else fp++;
		}
		else unknown++;
	}
	cout << "TP: " << tp << endl;
	cout << "TN: " << tn << endl;
	cout << "FP: " << fp << endl;
	cout << "FN: " << fn << endl;
	cout << "missing: " << missing << endl;
	cout << "unknown: " << unknown << endl;
	cout << "----------------------" << endl;
	const double tpr = 1.0 * tp / (tp + fn);
	cout << "TPR: " << tpr << endl;
	const double tnr = 1.0 * tn / (tn + fp);
	cout << "TNR: " << tnr << endl;
	cout << "Accuracy: " << 1.0 * (tp + tn) / (tp + tn + fp + fn) << endl;
	cout << "Mean(TPR, TNR): " << (tpr + tnr) / 2.0 << endl;	
}