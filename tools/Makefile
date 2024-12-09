CXX = g++
CXXFLAGS = -std=c++11 -O2 -g -lz

SRCS_BUILD = build_fai.cpp
OBJS_BUILD = $(SRCS_BUILD:.cpp=.o)
TARGET_BUILD = build_fai

SRCS_READ = read_fai.cpp
OBJS_READ = $(SRCS_READ:.cpp=.o)
TARGET_READ = read_fai

all: $(TARGET_BUILD) $(TARGET_READ)

$(TARGET_BUILD): $(OBJS_BUILD)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_BUILD)

$(TARGET_READ): $(OBJS_READ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_READ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS_BUILD) $(OBJS_READ) $(TARGET_BUILD) $(TARGET_READ)

.PHONY: all clean

