"""Defines the neb module methods"""

import os
import sys
import json
import shutil
import vasp
import casm
import vaspwrapper
from casm.vaspwrapper.vasp_calculator_base import VaspCalculatorBase
from vasp import Neb as calculator

class Neb(VaspCalculatorBase):
    """
    The Neb class contains functions for setting up, executing, and parsing a VASP neb calculation.

    Attributes
    ----------
    selection : casm.project.Selection
        selection of configuration
    calctype : string
        calctype to setup and run the neb calculations
    auto : bool
    sort : bool

    Methods
    -------
    from_configuration_dir(configdir='string', calctype='string', bool, bool)
        returns a instance of the Neb class instantited with a single configuration
    config_properties(config_data=dict/Pandas.DataFrame)
        return a dict of the properties required to setup a configuration
    pre_setup
        creates folder and makes POS files for each image
    setup
        sets up the input vasp files for the selection
    config_setup
        sets up the input vasp files for a single configuration
    get_vasp_input_files(config_data=dict/Pandas.DataFrame, settings=dict)
        returns filenames of a vasp neb calculation
    submit
        submit a job for each configuration
    run
        runs the neb calcutation on the selection
    report
        reports results for the selection
    run_cmd(configdir='string', calctype='string')
        return a string of command to run a single configuration
    finalize(config_data=dict/pandas_data, super_poscarfile='string')
        checks convergnce and write a properties file for the selection
    properties(calcdir='string', super_poscarfile='string', speciesfile='string')
        return a dict containing all the relaxed properties for a configuration

    Notes
    -----
    The calculation creates the following directory structure for each configuration:
    config/
        calctype.name/
            N_images_`n_images`\
                run.0/
                ....

    'run.i' directories are only created when ready.

    This automatically looks for VASP settings files using:
    casm.project.DirectoryStructure.settings_path_crawl

    The class inhertes from VaspCalculatorbase and methods overload the functionality in the parent

    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """Construct a VASP neb job object."""
        print "Construct a casm.vaspwrapper.Neb instance:"
        VaspCalculatorBase.__init__(self, selection, calctype, auto, sort)
        self.results_subdir = '01'
        self.calculator = calculator

    def config_properties(self, config_data):
        """return configuration properties as a dict"""
        config_dict = super(Neb, self).config_properties(config_data)
        try:
            n_images = json.load(open(config_dict["setfile"]))["n_images"]
            endstate_calctype = json.load(open(config_dict["setfile"]))["endstate_calctype"]
        except:
            ## Error message if "n_images" not present in settings
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
        config_dict["n_images"] = n_images
        config_dict["endstate_calctype"] = endstate_calctype
        config_dict["calcdir"] = os.path.join(config_dict["calcdir"], "N_images_{}".format(n_images))
        return config_dict

    def pre_setup(self):
        """Makes settings file for interpolation and POS folder/files for each image"""
        proj = self.selection.proj
        try:
            os.mkdir(os.path.join(proj.path, ".casm/tmp"))
        except:
            pass
        sel_tmp = self.selection.saveas(os.path.join(proj.path, ".casm/tmp", "neb_interpolation_selection_tmp"), force=True)
        dict = {}
        dict["selection"] = os.path.join(proj.path, ".casm/tmp", "neb_interpolation_selection_tmp")
        for index, config_data in self.selection.data.iterrows():
            conf_dict = {"n_images" : config_data["n_images"],
                         "endstate_calctype" : config_data["endstate_calctype"]}
            dict[config_data["name"]] = conf_dict
            try:
                os.makedirs(config_data["calcdir"])
            except:
                pass
            try:
                for i in range(config_data["n_images"] + 2):
                    os.makedirs(os.path.join(config_data["calcdir"], 'poscars', str(i).zfill(2)))
            except:
                pass
        dict["n_images"] = conf_dict["n_images"]
        dict["endstate_calctype"] = conf_dict["endstate_calctype"]
        dict["calctype"]=self.calctype
        tmp_folder = os.path.join(proj.path, '.casm/tmp')
        filename = "neb_interpolation_settings.json"
        with open(os.path.join(tmp_folder, filename), 'w') as file:
            file.write(json.dumps(dict, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))

        ## write the selection Interpolation command
        args = "enum --method DiffTransConfigInterpolation -s {}".format(os.path.join(tmp_folder, filename))
        output = proj.command(args)
        print "Interpolation complete"
        #os.remove(os.path.join(tmp_folder, filename))

    def setup(self):
        """ Setup initial relaxation run for the selection"""
        super(Neb, self).setup()

    def config_setup(self, config_data):
        """ Setup initial relaxation run for a configuration

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config_calcdir/00 directory:
                POSCAR: sample structure of the configuration to be relaxed

        """
        super(Neb, self).config_setup(config_data)
        ## settings the images tag in incar file
        tmp_dict = {"images": config_data["n_images"]}
        vasp.io.set_incar_tag(tmp_dict, config_data["calcdir"])

    def get_vasp_input_files(self, config_data, settings):
        """returns filenames of a vasp neb calculation"""
        vaspfiles = super(Neb, self).get_vasp_input_files(config_data, settings)
        incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files = vaspfiles
        super_poscarfile = os.path.join(config_data["calcdir"], "poscars", "00", "POSCAR")
        return incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files

    def submit(self):
        """submit a job for each configuration"""
        super(Neb, self).submit()

    def run(self):
        """runs the neb calcutation on the selection"""
        super(Neb, self).run()

    def report(self):
        """reports results for the selection"""
        super(Neb, self).report()

    @staticmethod
    def run_cmd(configdir, calctype):
        """return a string of command to run a single configuration"""
        return "python -c \"import casm.vaspwrapper; obj = casm.vaspwrapper.Neb.from_configuration_dir('{0}', '{1}'); obj.run()\"\n".format(configdir,
                                                                                                                                            calctype)

    def finalize(self, config_data, super_poscarfile=None):
        """checks convergnce and write a properties file for the selection"""
        if super_poscarfile is None:
            super_poscarfile = os.path.join(config_data["calcdir"], "run.final/00/POSCAR")
        super(Neb, self).finalize(config_data, super_poscarfile)
        images_calculated = [int(i.strip().split('_')[-1]) for i in os.listdir(os.path.split(config_data["calcdir"])[0]) if "N_images" in i]
        num_images = config_data["n_images"]
        if num_images == max(images_calculated):
            shutil.copy(os.path.join(config_data["calcdir"], "properties.calc.json"),
                        os.path.join(os.path.split(config_data["calcdir"])[0],
                                     "properties.calc.json"))
            print "As the present run has highest number of images copied {0} to {1}".format(os.path.join(config_data["calcdir"],
                                                                                                          "properties.calc.json"),
                                                                                             os.path.join(os.path.split(config_data["calcdir"])[0],
                                                                                                          "properties.calc.json"))
        sys.stdout.flush()

    def properties(self, calcdir, super_poscarfile=None, speciesfile=None):
        """Make properties output as a list of dict of each image properties"""
        final_output = {}
        endpts = json.load(open(os.path.join(calcdir,"endpoint_props.json")))
        final_output["00"]=endpts["0"]
        num_images = vasp.io.get_incar_tag("IMAGES", calcdir)
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            vaspdir = calcdir + "/{}".format(img)
            output = super(Neb, self).properties(vaspdir, super_poscarfile, speciesfile)
            final_output[img] = output
        final_output[str(num_images+1).zfill(2)]=endpts[str(num_images+1)]
        return final_output
