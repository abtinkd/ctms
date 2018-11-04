// convertWikipedia.cpp : Converts wikipedia standard format to our used format
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
int convertWiki(char fileLoc[], bool includeNeutralVotes)
{
	ifstream fin(fileLoc);
	ofstream fout;
	char* outFileLoc = (includeNeutralVotes ? "wikipedia_NeutralIncluded.txt" : "wikipedia.txt");
	fout.open(outFileLoc);
	if (!fin.is_open()) {cout << "ERROR Read"; return 1;}
	if (!fout.is_open()) {cout << "ERROR Write"; return 1;}
	fout << "# Directed graph: wikipedia neutral vote not included" << endl;
	fout << "# FromNodeId\tToNodeId\tSign" << endl;

	string line, FrmNId, ToNdId, Sgn;
	while (getline(fin, line)) {
		stringstream liness(line);
		string data;
		getline(liness, data, '\t');
		if (data == "U") {
			getline(liness, ToNdId, '\t');}
		else if (data == "V") {
			getline(liness, Sgn, '\t');
			getline(liness, FrmNId, '\t');
			if (Sgn != "0" || includeNeutralVotes) {
				fout << FrmNId << '\t' << ToNdId << '\t' << Sgn << '\n';}
		}
	}
	fin.close();
	fout.close();
	return 0;
}