default: install

build:
	R CMD build .

README.md: README.Rmd
	@echo "Building README.md"
	@R -e "rmarkdown::render('README.Rmd')"
	@rm README.html
