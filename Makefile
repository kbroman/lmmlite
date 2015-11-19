R_OPTS=--no-save --no-restore --no-init-file --no-site-file # --vanilla, but without --no-environ

all: doc extdata
.PHONY: doc extdata

# build package documentation
doc:
	R -e 'devtools::document()'

extdata: inst/extdata/recla.zip

inst/extdata/recla.zip: inst/extdata/create_example_data.R
	cd $(<D);R $(R_OPTS) -e "source('$(<F)')"
