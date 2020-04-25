default: install

build:
	R CMD build .

FLAGS:=
CPPFLAGS:=
CXXFLAGS:=
install: exports
	R CMD INSTALL --configure-args='CPPFLAGS="$(FLAGS) $(CPPFLAGS)" CXXFLAGS="$(FLAGS) $(CXXFLAGS)"' .

check:
	R CMD build .
	R CMD check transitr_*

document:
	R -e "devtools::document()"

FILTER ?= ".+"
test:
	R -e "devtools::load_all(); devtools::test(filter=\"${FILTER}\")"

clean:
	./cleanup

DEBUG ?= ""
ifeq ($(DEBUG), "")
	Rcmd = R
else
	Rcmd = R -d "$(DEBUG)"
endif

run:
	$(Rcmd) --slave -f scripts/run_model.R

startserver:
	cd simulations && yarn start &

SIM ?= sim000
simulation:
	$(Rcmd) --slave -f scripts/run_simulation.R --args $(SIM)

view:
	R --slave -f scripts/track_simulations.R

coverage:
	R -e "covr::report()"


exports:
	R -e "Rcpp::compileAttributes()"


DRY?=n
syncSims:
	rsync -avP$(DRY) --delete --exclude="sim1*" --exclude="sim0*" --exclude node_modules --exclude="*.zip" tell029@certellprd01.its.auckland.ac.nz:/data/transitr/simulations/ simulations


simprogress:
	@tail -n 1 simulations/$(SIM)/eta_state.csv | xargs -I{} R --slave -e "as.POSIXct(read.csv(textConnection(commandArgs(trailing=TRUE)[1]), header=F)[[3]], origin='1970-01-01')" --args {}
