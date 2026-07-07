RSCRIPT ?= Rscript

.PHONY: all run validate docs clean-final

all: run

run:
	$(RSCRIPT) run_all.R

validate:
	$(RSCRIPT) run_all.R
	$(RSCRIPT) -e 'stopifnot(file.exists("results/session/sessionInfo.txt"), file.exists("renv.lock"))'

docs:
	quarto render docs

clean-final:
	rm -rf data/processed/final results/final
	find results/tables results/figures results/logs results/session -mindepth 1 ! -name .gitkeep -exec rm -rf {} +
