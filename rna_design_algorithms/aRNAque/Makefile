################################################################################
#This script was taken from learna and edited on purpose
################################################################################


.DEFAULT_GOAL := show-help
SHELL := /bin/bash
PATH := $(PWD)/thirdparty/miniconda/miniconda/bin:$(PATH)


################################################################################
# Utility
################################################################################

## To clean project state
clean: clean-runtime clean-log-file clean-plots clean-thirdparty

## Remove runtime files
clean-runtime:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -rf {} +
	find . -name '.DS_Store' -exec rm -rf {} +

## Remove data files
clean-log-file:
	rm -rf data/Log/*


## Clean plots
clean-plots:
	rm -rf images/*.pdf
	rm -rf images/*.png
	rm -rf images/Eterna100/*
	rm -rf images/PseudoBase++/*
	rm -rf images/RFAM/*


## Remove thirdparty installs
clean-thirdparty:
	rm -rf thirdparty/miniconda/miniconda


################################################################################
# Setup LEARNA
################################################################################

## Install all dependencies
requirements:
	./thirdparty/miniconda/make_miniconda.sh
	conda env create -f environment.yml


## Plot reproduced results using pgfplots
plots:


################################################################################
# Help
################################################################################

# From https://drivendata.github.io/cookiecutter-data-science/
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
