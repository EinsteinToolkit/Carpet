# Lock buffer_zone_width. One RK4 step with 3 ghosts is 9 zones.
# A 2:1 subcycled step is 21 zones.

CXX = g++
STD = -std=c++17

.PHONY: test clean

test: buffer_zones_test
	./buffer_zones_test

buffer_zones_test: buffer_zones_test.cc ../../src/buffer_zones.hh
	$(CXX) -Wall -Werror $(STD) -I../../src buffer_zones_test.cc -o $@

clean:
	rm -f buffer_zones_test
