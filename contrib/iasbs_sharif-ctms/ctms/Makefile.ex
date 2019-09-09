#
#	configuration variables for the example

## Main application file
MAIN = ctmsmain
CXXFLAGS += -std=c++11
DEPH = ../$(EXSNAPEXP)/signnet.h ctmsnet.h sign_prediction.h ml.h convert_wikipedia.h
DEPCPP = ../$(EXSNAPEXP)/signnet.cpp ctmsnet.cpp sign_prediction.cpp ml.cpp convert_wikipedia.cpp

