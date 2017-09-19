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

def FinalizeBase(object):
    """Base class fro Finalize
    """
    def __init__(self, config_obj):
        """
        Construct a Finalize base object

        Arguments
        ---------
            config_obj: casm.vaspwrapper.Configproperties object

        """
        self.config_obj = config_obj

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


    def finalize(self):
        if self.is_converged():
            # write properties.calc.json
            vaspdir = os.path.join(self.config_obj.calcdir, "run.final")
            speciesfile = self.config_obj.casm_directories.settings_path_crawl("SPECIES", self.config_obj.configname,
                                                                          self.config_obj.clex, self.config_obj.calc_subdir)
            output = self.config_obj.properties(vaspdir, True, speciesfile)
            outputfile = os.path.join(self.config_obj.calcdir, "properties.calc.json")
            with open(outputfile, 'w') as file:
                file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
            print "Wrote " + outputfile
            all_image_folders = [int(i.strip().split('_')[-1]) for i in os.listdir(self.config_obj.calcdir) if "N_images" in i]
            num_images = [int(self.config_obj.calcdir.strip().split('_')[-1])]
            if num_images == all_image_folders:
                shutil.copy(os.path.join(self.config_obj.calcdir, "properties.calc.json"),
                            os.path.join(os.path.split(self.config_obj.calcdir)[0], "properties.calc.json"))
                print "As the present run has highest number of images copied {0} to {1}".format(os.path.join(self.config_obj.calcdir,
                                                                                                              "properties.calc.json"),
                                                                                                 os.path.join(os.path.split(self.config_obj.calcdir)[0],
                                                                                                              "properties.calc.json"))
            sys.stdout.flush()
            report_status('complete', self.config_obj)

    def is_converged(self):
        # Check for electronic convergence in completed calculations. Returns True or False.

        # Verify that the last relaxation reached electronic convergence
        calculation = vasp.Neb(self.calcdir, self.run_settings())
        for i in range(len(calculation.rundir)):
            try:
                #vrun = vasp.io.Vasprun(os.path.join(self.calcdir, calculation.rundir[-i-1], "vasprun.xml"))
                vrun_oszicar = vasp.io.Oszicar(os.path.join(self.calcdir, calculation.rundir[-i-1], "01", "OSZICAR"))
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
        vrun_oszicar = vasp.io.Oszicar(os.path.join(self.calcdir, "run.final", "01", "OSZICAR"))
        vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(self.calcdir, "run.final"))
        if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
        if vrun_oszicar.num_elm[-1] >= vrun_nelm:
            print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
            self.report_status('failed', 'electronic_convergence')
            return False

        return True

    @staticmethod
    def properties(vaspdir, use_poscarfile=False, speciesfile=None):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""
        final_output = []
        num_images = vasp.io.get_incar_tag("IMAGES", vaspdir)
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            output = dict()
            #vrun = vasp.io.Vasprun( os.path.join(vaspdir, "vasprun.xml") )
            vrun_oszicar = vasp.io.Oszicar(os.path.join(vaspdir, img, "OSZICAR"))
            vrun_outcar = vasp.io.Outcar(os.path.join(vaspdir, img, "OUTCAR"))

            # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

            if (use_poscarfile is not None) and (speciesfile is not None):
                species_settings = vasp.io.species_settings(speciesfile)
                super_poscar = vasp.io.Poscar(os.path.join(vaspdir, img, "POSCAR"), species_settings)
                super_contcar = vasp.io.Poscar(os.path.join(vaspdir, img, "CONTCAR"), species_settings)
                unsort_dict = super_poscar.unsort_dict()
            else: # not implemented #TODO
                # fake unsort_dict (unsort_dict[i] == i)
                unsort_dict = dict(zip(range(0, len(vrun.basis)), range(0, len(vrun.basis))))
                super_poscar = vasp.io.Poscar(os.path.join(vaspdir, "POSCAR"))

            # unsort_dict:
            #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
            #   unsorted_dict[sorted_index] == orig_index
            #   For example:
            #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR

            output["Image_number"] = img
            output["atom_type"] = super_poscar.type_atoms
            output["atoms_per_type"] = super_poscar.num_atoms
            output["coord_mode"] = super_poscar.coord_mode

            # as lists
            output["relaxed_forces"] = [ None for i in range(len(vrun_outcar.forces))]
            for i, v in enumerate(vrun_outcar.forces):
                output["relaxed_forces"][unsort_dict[i] ] = casm.NoIndent(vrun_outcar.forces[i])

            output["relaxed_lattice"] = [casm.NoIndent(list(v)) for v in super_contcar.lattice()]
            output["relaxed_basis"] = [None for i in range(len(super_contcar.basis))]
            for i, v in enumerate(super_contcar.basis):
                output["relaxed_basis"][unsort_dict[i]] = casm.NoIndent(list(super_contcar.basis[i].position))

            output["relaxed_energy"] = vrun_oszicar.E[-1]
            final_output.append(output)

        return final_output
