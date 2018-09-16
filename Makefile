default: install

build:
	R CMD build .

FLAGS:=
CPPFLAGS:=
CXXFLAGS:=
install: exports
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

startserver:
	cd simulations && yarn start &

SIM ?= sim000
simulation:
	R --slave -f scripts/run_simulation.R --args $(SIM)

view:
	R --slave -f scripts/track_simulations.R

coverage:
	R -e "covr::report()"


exports:
	R -e "Rcpp::compileAttributes()"
