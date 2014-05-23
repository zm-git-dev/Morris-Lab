
source: qarp_0.1-1.tar.gz
	mv qarp_0.1-1.tar.gz qarp/
	tar czf qarp.tar.gz qarp



# Build a binary R package for qarp.
# should everything go in there, including python code and shiny code?
qarp_0.1-1.tar.gz:
	(cd qarp/qarp && ../../roxygenize.R)
	R CMD build qarp/qarp


binary:
	(cd qarp/qarp && ../../roxygenize.R)
	R CMD INSTALL --build --preclean qarp/qarp

.PHONY: source binary
