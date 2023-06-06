VERSION=$(shell grep Version DESCRIPTION | sed 's/Version: //')

default: install

install: build
	@echo "Installing package"
	@R CMD INSTALL transitr_$(VERSION).tar.gz

build: transitr_$(VERSION).tar.gz

transitr_$(VERSION).tar.gz: DESCRIPTION NAMESPACE R/*.R
	@echo "Building package"
	@R CMD build .

README.md: README.Rmd
	@echo "Building README.md"
	@R -e "rmarkdown::render('README.Rmd')"
	@rm README.html
