include $(top_srcdir)/common.mk

install:
	cd $(top_srcdir)/python/casm && $(PYTHON) setup.py install --prefix=$(PYTHON_PREFIX) "$(GIT_HASH)"
	cd $(top_srcdir)/python/vasp && $(PYTHON) setup.py install --prefix=$(PYTHON_PREFIX) "$(GIT_HASH)"




##The following is enough to install the python module, bypassing setup.py. Seriously no uninstall or clean up from distutils? Weak.
#autopy_PYTHON=autopy/__init__.py autopy/noindent.py autopy/casm.py
#learn_PYTHON=autopy/learn/fit.py autopy/learn/linear_model.py autopy/learn/evolve.py autopy/learn/cross_validation.py autopy/learn/feature_selection.py autopy/learn/tools.py autopy/learn/__init__.py
#vaspwrapper_PYTHON=autopy/vaspwrapper/vaspwrapper.py autopy/vaspwrapper/converge.py autopy/vaspwrapper/relax.py autopy/vaspwrapper/__init__.py 
#project_PYTHON=autopy/project/project.py autopy/project/query.py autopy/project/selection.py autopy/project/__init__.py autopy/project/io.py 
#plotting_PYTHON=autopy/plotting/__init__.py autopy/plotting/plotting.py 
#
#autopydir=$(pythondir)/autopy
#learndir=$(autopydir)/learn
#vaspwrapperdir=$(autopydir)/vaspwrapper
#projectdir=$(autopydir)/project
#plottingdir=$(autopydir)/plotting
#
#
#dist_bin_SCRIPTS=scripts/auto.casm-learn scripts/auto.casm.view scripts/auto.vasp.kpoint_converge scripts/auto.vasp.relax scripts/auto.vasp.relax.report scripts/auto.vasp.setup 
#
##Delete empty directories when unistalling
#uninstall-hook:
#	if test -d ${learndir}; then rmdir ${learndir}; fi
#	if test -d ${vaspwrapperdir}; then rmdir ${vaspwrapperdir}; fi
#	if test -d ${projectdir}; then rmdir ${projectdir}; fi
#	if test -d ${plottingdir}; then rmdir ${plottingdir}; fi
#	if test -d ${autopydir}; then rmdir ${autopydir}; fi
