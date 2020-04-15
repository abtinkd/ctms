## Main application file
MAIN = extractormain
CXXFLAGS += -std=c++11
LDFLAGS += -I../ctms
DEPH = ../$(EXSNAPEXP)/signnet.h ../ctms/ctmsnet.h
DEPCPP = ../$(EXSNAPEXP)/signnet.cpp ../ctms/ctmsnet.cpp

