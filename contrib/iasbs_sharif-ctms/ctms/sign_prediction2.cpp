#include "sign_prediction2.h"

void TNaivePredictor::build() {

}

void TBalanceBasedPredictor::build() {

}

void TStatusBasedPredictor::build() {

}

void::TLogisticRegression::extractFeatures() {

}

void TLogisticRegression::build() {

}

void TCTMSProbabilisticInference::extractFeatures() {

}

void TCTMSProbabilisticInference::build() {

}

int TNaivePredictor::predict(const TInt srcNId, const TInt desNId) {	
	if (pType == "gnr"){}
	else if (pType == "rcp") {}
	else if (pType == "cmp") {}
	else if (pType == "wgr") {}
	else {// to do random
	}
	return 1;
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
	return 1;
}