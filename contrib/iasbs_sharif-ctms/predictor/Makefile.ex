## Main application file
MAIN = predictormain
CXXFLAGS += -std=c++11
LDFLAGS += -I../ctms
DEPH = ../$(EXSNAPEXP)/signnet.h ../ctms/ctmsnet.h ../ctms/sign_predictor.h ../ctms/ml.h 
DEPCPP = ../$(EXSNAPEXP)/signnet.cpp ../ctms/ctmsnet.cpp ../ctms/sign_predictor.cpp ../ctms/ml.cpp

