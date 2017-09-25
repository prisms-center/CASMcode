import os, math, sys, json, re, warnings, shutil
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from configproperties import ConfigPropertiesBase
from report import FinalizeBase, PropertiesBase, ContainerBase

class ConfigProperties(ConfigPropertiesBase):
    """ Derived class with additional configuration properties and errors

    Additional Attributes
    --------------------

      n_images: int
        number of images in the calculation

    """

    def __init__(self, configname, calctype):
        """
             Add additional properties and error messages
        """
        ConfigPropertiesBase.__init__(self, configname, calctype)

        ## Error message if "n_images" not present in settings
        if not "n_images" in self.settings:
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()

        self.n_images = self.settings["n_images"]
        # append the n_images to calcdir
        self.calcdir = os.path.join(self.calcdir, "N_images_{}".format(self.n_images))
        self.results_subdir = '01'

class Neb(object):
    """The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              ....

        'run.i' directories are only created when ready.

        This automatically looks for VASP settings files using:
          casm.project.DirectoryStructure.settings_path_crawl

    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """
        Construct a VASP neb job object.

        Arguments
        ----------

            selection: casm.project.Selection object, default= yet to be implemented #todo
              Selection of all DiffTransConfigurations to submit a NEB calculation. 
              default should be MASTER selection and yet to be implemented

            sort: boolean, optional, default=True,
              Use True to sort atoms in POSCAR by type

        """
        print "Construct a casm.vaspwrapper.Relax instance:"

        self.selection = selection
        self.calctype = calctype
        self.auto = auto
        self.sort = sort

    def _pre_setup(self):
        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data["configname"], self.calctype)
            try:
                os.makedirs(config_obj.calcdir)
                for i in range(config_obj.n_images + 2):
                    os.makedirs(os.path.join(config_obj.calcdir, str(i).zfill(2)))
            except:
                pass

    def setup(self):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed

        """
        proj = self.selection.proj
        dict = {}
        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data["configname"], self.calctype)
            conf_dict = {"n_images" : config_obj.n_images,
                         "calctype" : config_obj.clex.calctype}
            dict[config_obj.configname] = conf_dict
            try:
                os.makedirs(config_obj.calcdir)
                for i in range(config_obj.n_images + 2):
                    os.makedirs(os.path.join(config_obj.calcdir, str(i).zfill(2)))
            except:
                pass

        filename = "_neb_tmp.json"
        with open(filename, 'w') as file:
            file.write(json.dumps(dict, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))

        ## write the selection Interpolation command
        args = "enum --method DiffTransConfigInterpolation -j {}".format(filename)
        output = proj.command(args)
        os.remove(filename)

        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data["configname"], self.calctype)
            # Find required input files in CASM project directory tree
            vaspfiles = casm.vaspwrapper.vasp_input_file_names(config_obj.casm_directories,
                                                               config_obj.configname,
                                                               config_obj.clex,
                                                               config_obj.calc_subdir)
            incarfile, prim_kpointsfile, prim_poscarfile, temp_poscarfile, speciesfile = vaspfiles
            # Find optional input files
            extra_input_files = []
            for s in config_obj.settings["extra_input_files"]:
                extra_input_files.append(config_obj.casm_directories.settings_path_crawl(s, config_obj.configname,
                                                                                         config_obj.clex, config_obj.calc_subdir))
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. Extra input file " + s + " not found in CASM project.")
            if config_obj.settings["initial"]:
                extra_input_files += [config_obj.casm_directories.settings_path_crawl(config_obj.settings["initial"],
                                                                                      config_obj.configname, config_obj.clex,
                                                                                      config_obj.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No initial INCAR file " + config_obj.settings["initial"] + " found in CASM project.")
            if config_obj.settings["final"]:
                extra_input_files += [config_obj.casm_directories.settings_path_crawl(config_obj.settings["final"],
                                                                                      config_obj.configname, config_obj.clex,
                                                                                      config_obj.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No final INCAR file " + config_obj.settings["final"] + " found in CASM project.")
            sys.stdout.flush()

            #make vasp input files
            sample_super_poscarfile = os.path.join(config_obj.calcdir, "00", "POSCAR")
            vasp.io.write_vasp_input(config_obj.calcdir, incarfile, prim_kpointsfile, prim_poscarfile,
                                     sample_super_poscarfile, speciesfile, self.sort, extra_input_files,
                                     config_obj.settings["strict_kpoints"])
            ## settings the images tag in incar file
            tmp_dict = {"images": config_obj.n_images}
            vasp.io.set_incar_tag(tmp_dict, config_obj.calcdir)

class Finalize(FinalizeBase):
    """Finalize object specific to this class"""

    def __init__(self, config_obj):
        """initialze the Finalize object"""
        FinalizeBase.__init__(self, config_obj)

    def finalize(self):
        if self.is_converged():
            properties_obj = Properties(self.config_obj)
            super_poscarfile = os.path.join(self.config_obj.calcdir, "00/POSCAR")
            super(Finalize, self).finalize(properties_obj, super_poscarfile)
            all_image_folders = [int(i.strip().split('_')[-1]) for i in os.listdir(self.config_obj.calcdir) if "N_images" in i]
            num_images = self.config_obj.n_images
            if num_images == max(all_image_folders):
                shutil.copy(os.path.join(self.config_obj.calcdir, "properties.calc.json"),
                            os.path.join(os.path.split(self.config_obj.calcdir)[0],
                                         "properties.calc.json"))
                print "As the present run has highest number of images copied {0} to {1}".format(os.path.join(self.config_obj.calcdir,
                                                                                                              "properties.calc.json"),
                                                                                                 os.path.join(os.path.split(self.config_obj.calcdir)[0],
                                                                                                              "properties.calc.json"))
            sys.stdout.flush()

class Properties(PropertiesBase):

    def __init__(self, config_obj):
        PropertiesBase.__init__(self, config_obj)

    def properties(self, super_poscarfile=None, speciesfile=None):
        """Make properties output as a list of dict of each image properties"""
        final_output = []
        num_images = vasp.io.get_incar_tag("IMAGES", self.config_obj.calcdir)
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            vaspdir = self.config_obj.calcdir + '0' + str(img)
            output = super(Properties, self).properties(vaspdir, super_poscarfile, speciesfile)
            output["Image_number"] = img
            final_output.append(output)
        return final_output

class Container(ContainerBase):

    def __init__(self, configname, calctype=None):
        ContainerBase.__init__(self, configname, calctype, ConfigProperties, Finalize, Properties)
