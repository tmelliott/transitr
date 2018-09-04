default: install

build:
	R CMD build .

FLAGS:=
CPPFLAGS:=
CXXFLAGS:=
install:
	R CMD INSTALL --configure-args='CPPFLAGS=$(FLAGS) $(CPPFLAGS) CXXFLAGS=$(FLAGS) $(CXXFLAGS)' .

check:
	R CMD build .
	R CMD check transitr_*

document:
	R -e "devtools::document()"

test:
	R -e "devtools::load_all(); devtools::test()"

clean:
	./cleanup

run:
	R --slave -f scripts/run_model.R

coverage:
	R -e "covr::report()"
