## Main application file
MAIN = predictormain
CXXFLAGS += -std=c++11
DEPH = ../$(EXSNAPEXP)/signnet.h ../ctms/ctmsnet.h ../ctms/sign_prediction.h ../ctms/ml.h 
DEPCPP = ../$(EXSNAPEXP)/signnet.cpp ../ctms/ctmsnet.cpp ../ctms/sign_prediction.cpp ../ctms/ml.cpp

