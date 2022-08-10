W = -Wall -Werror -Wextra

all: clean s21_matrix.a test gcov_report

s21_matrix.a: s21_matrix.c
	@gcc -c $(W) s21_matrix.c -o s21_matrix.o
	@ar -rcs s21_matrix.a s21_matrix.o
	@rm -f *.o

test: test.c s21_matrix.a
	@gcc s21_matrix.a test.c -lcheck -o test

clean:
	@rm -rf *.o s21_matrix.a *.out report gcov_report *.gcda *.gcno *.info test Test

check:
	@cppcheck --enable=all --suppress=missingIncludeSystem --inconclusive --check-config *.c *.h
	@cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg  
	@python3 ../materials/linters/cpplint.py --extension=c *.c *.h
	@rm -rf CPPLINT.cfg

gcov_report: s21_matrix.a
	@gcc --coverage test.c -lcheck s21_matrix.c -o Test
	@./Test
	@lcov -t "Test" -o Test.info -c -d .
	@genhtml -o report Test.info
	@open ./report/index.html

install lcov:
	@curl -fsSL https://rawgit.com/kube/42homebrew/master/install.sh | zsh
	@brew install lcov