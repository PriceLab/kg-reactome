default:
	@echo "targets:"
	@echo "   roxy: build roxygen documentation and NAMESPACE"
	@echo "   install"
	@echo "   demo"

roxy:
	R -e "devtools::document()"

install:
	R CMD INSTALL .

test:
	(cd inst/unitTests; R -f test_ReactionParser.R)


all: roxy install test
quick: install test
