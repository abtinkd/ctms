SNAPDIR=../../../snap
INCLUDES= -I $(SNAPDIR)/snap-core -I $(SNAPDIR)/glib-core -I $(SNAPDIR)/snap-adv -I $(SNAPDIR)/snap-exp
CC=g++
CFLAGS=-c -O3 

SOURCES= ctmsmain.cpp sign_prediction.cpp ml.cpp ctmsnet.cpp convert_wikipedia.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=ctmsmain
SNAPLIB=$(SNAPDIR)/snap-core/Snap.o

all:	$(SNAPLIB) $(SOURCES) $(EXECUTABLE)

$(SNAPLIB):	
	$(MAKE) -C $(SNAPDIR)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(SNAPLIB)

.cpp.o:
	$(CC) $(CFLAGS) $< $(INCLUDES) -o $@
clean:
	rm -f $(EXECUTABLE) *.o