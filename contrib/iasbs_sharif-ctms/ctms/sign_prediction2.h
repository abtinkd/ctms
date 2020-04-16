#ifndef SIGN_PREDICTION2_H
#define SIGN_PREDICTION2_H

#include "ctmsnet.h"
#include <Snap.h>
#include "ml.h"

class TSignPredictor {
protected:	
	PCtmsNet network;	
public:
	TSignPredictor(const PCtmsNet& net) : network(net) {}
	virtual void build() = 0;
	virtual int predict(const TInt srcNId, const TInt desNId) = 0;
};

class TNaivePredictor : public TSignPredictor {
protected:
	const char* pType;
public:
	TNaivePredictor(const PCtmsNet& net, const char* type) : TSignPredictor(net), pType(type) {}
	void build();
	int predict(const TInt srcNId, const TInt desNId);
};

class TBalanceBasedPredictor : public TSignPredictor {
public:
	TBalanceBasedPredictor(const PCtmsNet& net) : TSignPredictor(net) {}
	void build();
	int predict(const TInt srcNId, const TInt desNId);
};

class TStatusBasedPredictor : public TSignPredictor {
public:
	TStatusBasedPredictor(const PCtmsNet& net) : TSignPredictor(net) {}
	void build();
	int predict(const TInt srcNId, const TInt desNId);
};

class TLogisticRegression : public TSignPredictor {
protected:
	THashSet<TChA> features; //TChA is feature string. Its corresponding Int ID is hash KeyId (int KeyId -> TChA Feature String)
	void extractFeatures();
public:
	TLogisticRegression(const PCtmsNet& net) : TSignPredictor(net) {}
	void build();
	int predict(const TInt srcNId, const TInt desNId);
};

class TCTMSProbabilisticInference : public TSignPredictor {
protected:
	THash<TChA, TFlt> theta;
	void extractFeatures(const TIntPr& edge, THash<TChA, TInt>& edgeFeaValues);
public:
	TCTMSProbabilisticInference(const PCtmsNet& net) : TSignPredictor(net) {}
	void build();
	int predict(const TInt srcNId, const TInt desNId);
};
#endif // !SIGN_PREDICTION2_H