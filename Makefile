default: install

build:
	R CMD build .

install:
	R CMD INSTALL .

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
	R -f scripts/run_model.R

coverage:
	R -e "covr::report()"
	