DYNAMODIR=/home/mjki2mb2/dynamo
COILDIR=$(DYNAMODIR)/src/coil
MAGNETDIR=$(DYNAMODIR)/src/magnet
COILLIBDIR=$(DYNAMODIR)/lib/

TARGET = $(shell find . -name "*.cpp" | head -n 1 | sed 's/\.cpp$$//')
CXXFLAGS = -O3 -g -Wall -Werror -ansi -I$(COILDIR) -I$(MAGNETDIR) `pkg-config gtkmm-2.4 --cflags`
LINKFLAGS = -L $(COILLIBDIR) -lcoil `pkg-config gtkmm-2.4 --libs`

all: $(TARGET)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET) : $(TARGET).o 
	$(CXX) $(CXXFLAGS) $< $(LINKFLAGS)  -o $@

clean :	
	rm -f *.o $(TARGET)
