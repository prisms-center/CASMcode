import os, math, sys, json, re, warnings, shutil
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from configproperties import ConfigPropertiesBase

class Report(object):
    """The Report class contains functions for reporting the results of a selection.
    """
    def __init__(self, selection):
        """
        Construct a Report job object.

        Arguments
        ----------

            selection: casm.project.Selection object, default= yet to be implemented #todo
              Selection of all DiffTransConfigurations to submit a NEB calculation.
              default should be MASTER selection and yet to be implemented
        """
        print "Construct a casm.vaspwrapper.Report  instance:"
        self.selection = selection

class FinalizeBase(object):
    """Base class for Finalize
    """
    def __init__(self, config_obj, method = "relax"):
        """
        Construct a Finalize base object

        Arguments
        ---------
            config_obj: casm.vaspwrapper.Configproperties object

        """
        self.config_obj = config_obj
        self.method = method

    def report_status(self, status, failure_type=None):
        """Report calculation status to status.json file in configuration directory.

            Args:
                status: string describing calculation status. Currently used values are
                    not_submitted
                    submitted
                    complete
                    failed
                failure_type: optional string describing reason for failure. Currently used values are
                    unknown
                    electronic_convergence
                    run_limit"""

        output = dict()
        output["status"] = status
        if failure_type is not None:
            output["failure_type"] = failure_type

        outputfile = os.path.join(self.config_obj.calcdir, "status.json")
        with open(outputfile, 'w') as file:
            file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
        print "Wrote " + outputfile
        sys.stdout.flush()


    def finalize(self, properties_obj, super_poscarfile = None):
        if self.is_converged():
            # write properties.calc.json
            vaspdir = os.path.join(self.config_obj.calcdir, "run.final")
            speciesfile = self.config_obj.casm_directories.settings_path_crawl("SPECIES", self.config_obj.configname,
                                                                               self.config_obj.clex, self.config_obj.calc_subdir)
            output = properties_obj.properties(vaspdir, super_poscarfile, speciesfile)
            outputfile = os.path.join(self.config_obj.calcdir, "properties.calc.json")
            with open(outputfile, 'w') as file:
                file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
            print "Wrote " + outputfile
            sys.stdout.flush()
            self.report_status('complete')

    def is_converged(self):
        # Check for electronic convergence in completed calculations. Returns True or False.

        # Verify that the last relaxation reached electronic convergence
        available_methods = {'relax': vasp.Relax, 'neb': vasp.Neb}
        calculation = available_methods[self.method](self.config_obj.calcdir, self.config_obj.settings)
        for i in range(len(calculation.rundir)):
            try:
                #vrun = vasp.io.Vasprun(os.path.join(self.calcdir, calculation.rundir[-i-1], "vasprun.xml"))
                vrun_oszicar = vasp.io.Oszicar(os.path.join(self.config_obj.calcdir, calculation.rundir[-i-1],
                                                            self.config_obj.results_subdir, "OSZICAR"))
                vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(self.calcdir, calculation.rundir[-i-1]))
                if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
                if len(vrun_oszicar.num_elm[-1]) >= vrun_nelm:
                    print('The last relaxation run (' +
                          os.path.basename(relaxation.rundir[-i-1]) +
                          ') failed to achieve electronic convergence; properties.calc.json will not be written.\n')
                    self.report_status('failed', 'electronic_convergence')
                    return False
                break
            except:
                pass

        # Verify that the final static run reached electronic convergence
        #vrun = vasp.io.Vasprun(os.path.join(self.calcdir, "run.final", "vasprun.xml"))
        vrun_oszicar = vasp.io.Oszicar(os.path.join(self.config_obj.calcdir, "run.final",
                                                    self.config_obj.results_subdir, "OSZICAR"))
        vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(self.config_obj.calcdir, "run.final"))
        if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
        if vrun_oszicar.num_elm[-1] >= vrun_nelm:
            print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
            self.report_status('failed', 'electronic_convergence')
            return False

        return True

class PropertiesBase(object):

    def __init__(self, config_obj):
        self.config_obj = config_obj

    @staticmethod
    def properties(vaspdir, super_poscarfile=None, speciesfile=None):
        """ return a dict of output form a vasp directory"""

        output = dict()
        # load the OSZICAR and OUTCAR
        zcar = vasp.io.Oszicar(os.path.join(vaspdir, "OSZICAR"))
        ocar = vasp.io.Outcar(os.path.join(vaspdir, "OUTCAR"))

        # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

        if (super_poscarfile is not None) and (speciesfile is not None):
            species_settings = vasp.io.species_settings(speciesfile)
            super_poscar = vasp.io.Poscar(super_poscarfile, species_settings)
            unsort_dict = super_poscar.unsort_dict()
        else:
            # fake unsort_dict (unsort_dict[i] == i)
            super_poscar = vasp.io.Poscar(os.path.join(vaspdir, "POSCAR"))
            unsort_dict = dict(zip(range(0, len(super_poscar.basis)),
                                   range(0, len(super_poscar.basis))))
        super_contcar = vasp.io.Poscar(os.path.join(vaspdir, "CONTCAR"))

        # unsort_dict:
        #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
        #   unsorted_dict[sorted_index] == orig_index
        #   For example:
        #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR


        output["atom_type"] = super_poscar.type_atoms
        output["atoms_per_type"] = super_poscar.num_atoms
        output["coord_mode"] = super_poscar.coord_mode

        # as lists
        output["relaxed_forces"] = [None for i in range(len(ocar.forces))]
        for i, force in enumerate(ocar.forces):
            output["relaxed_forces"][unsort_dict[i]] = casm.NoIndent(force)

        output["relaxed_lattice"] = [casm.NoIndent(list(v)) for v in super_contcar.lattice()]
        output["relaxed_basis"] = [None for i in range(len(super_contcar.basis))]
        for i, ba in enumerate(super_contcar.basis):
            output["relaxed_basis"][unsort_dict[i]] = casm.NoIndent(list(ba.position))

        output["relaxed_energy"] = zcar.E[-1]

        if ocar.ispin == 2:
            output["relaxed_magmom"] = zcar.mag[-1]
            if ocar.lorbit in [1, 2, 11, 12]:
                output["relaxed_mag_basis"] = [None for i in range(len(super_contcar.basis))]
                for i, v in enumerate(super_contcar.basis):
                    output["relaxed_mag_basis"][unsort_dict[i]] = casm.NoIndent(ocar.mag[i])

        return output

class ContainerBase(object):
    """Contains the finalize and properties objects
    """
    def __init__(self, configname, ConfigProperties, Finalize, Properties):
        self.config_obj = ConfigProperties(configname)
        self.finalize_obj = Finalize(self.config_obj)
        self.properties_obj = Properties(self.config_obj)

    def get_objects(self):
        return self.config_obj, self.finalize_obj, self.properties_obj
    
