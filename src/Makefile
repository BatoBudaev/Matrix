CC = g++
CFLAGS = -Wall -Wextra -Werror -std=c++17 -c -g
OS := $(shell uname)
TEST_FLAGS = $(shell pkg-config --cflags --libs gtest) -lm

all: clean s21_matrix_oop.a test

test: clean s21_test.o s21_matrix_oop.a
	$(CC) s21_test.o s21_matrix_oop.a -o test -L. $(TEST_FLAGS)
	./test

s21_test.o: s21_test.cpp
	$(CC) $(CFLAGS) s21_test.cpp -o s21_test.o

s21_matrix_oop.a: s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

s21_matrix_oop.o: s21_matrix_oop.cpp
	$(CC) $(CFLAGS) s21_matrix_oop.cpp -o s21_matrix_oop.o

clean:
	rm -rf *.o *.a *.out *.info *.html *.gc* test gcov_report s21_matrix_oop report

gcov_report:
	$(CC) -fprofile-arcs -ftest-coverage -std=c++17 s21_test.cpp s21_matrix_oop.cpp -o gcov_report $(TEST_FLAGS)
	./gcov_report
	lcov -t "test" -o test.info -c -d ./
	lcov -e test.info "s21_test.cpp" "s21_matrix_oop.cpp" -o filtered_test.info
	genhtml -o report filtered_test.info

open:
ifeq ($(OS),Darwin)
	open -a "Google Chrome" report/index.html
else
	xdg-open report/index.html
endif

clang:
	clang-format -n -style=google *.cpp *.h

clangI:
	clang-format -i -style=google *.cpp *.h

leak:
ifeq ($(OS),Darwin)
	CK_FORK=no leaks -atExit -- ./test
else
	valgrind --leak-check=full -s --track-origins=yes ./test
endif

install_gtest:
ifeq ($(OS),Darwin)
	brew update
	brew install googletest
else
	sudo apt-get update
	sudo apt install -y libgtest-dev
	sudo apt install -y cmake
endif
	
rebuild: all