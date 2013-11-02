
mds:
	(cd mds && ../roxygenize.R)
	R CMD build mds

binary:
	(cd mds && ../roxygenize.R)
	R CMD INSTALL --build --preclean mds

.PHONY: mds
