CC:=g++
CFLAGS:= -Wall -g -lmgl -Wno-return-type-c-linkage
LDFLAGS:=-std=gnu++11 -lmgl
INC:=-I /usr/local/include/
EXECUTABLE:=MrGBP

all: $(EXECUTABLE)

$(EXECUTABLE):
	$(CC) $(LDFLAGS) *.cpp -o $(EXECUTABLE) $(INC) $(LDFLAGS)

$(OBJECTS): %.o: %.cpp
	$(CC) $(CFLAGS) $< —o $@ $(INC)

clean:
	rm $(EXECUTABLE)
