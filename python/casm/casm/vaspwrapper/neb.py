import os, math, sys, json, re, warnings, shutil
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from configproperties import ConfigPropertiesBase
from report import FinalizeBase, PropertiesBase, ContainerBase

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
        """
        print "Construct a casm.vaspwrapper.Neb instance:"
        VaspCalculatorBase.__init__(selection, calctype, auto, sort)
        self.append_selection_data()
        self.results_subdir = '01'

    def append_selection_data(self):
        """append configproperties to selection.data"""
        for config_data in self.selection.data:
            selection.data.append(self.config_properties(config_data))

    def config_properties(self, config_data):
        """configuration properties as a dict"""
        config_dict = super(Neb, self).config_properties(config_data)
        config_data["calcdir"] =  os.path.join(config_data["calcdir"], "N_images_{}".format(n_images))
        try:
            n_images = json.load(config_data["setfile"])["n_images"]
        except:
            ## Error message if "n_images" not present in settings
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()
        config_data["n_images"] = n_images
        return config_dict
        
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
            conf_dict = {"n_images" : config_data["n_images"],
                         "calctype" : config_data["calctype"]}
            dict[config_data["configname"] = conf_dict
            try:
                os.makedirs(config_data["calcdir"])
                for i in range(config_data["n_images"] + 2):
                    os.makedirs(os.path.join(config_data["calcdir"], str(i).zfill(2)))
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
            # Find required input files in CASM project directory tree
            vaspfiles = casm.vaspwrapper.vasp_input_file_names(self.casm_directories,
                                                               config_data["configname"],
                                                               self.clex,
                                                               self.calc_subdir)
            incarfile, prim_kpointsfile, prim_poscarfile, temp_poscarfile, speciesfile = vaspfiles
            settings = self.read_settings(config_data["setfile"])
            # Find optional input files
            extra_input_files = []
            for s in settings["extra_input_files"]:
                extra_input_files.append(self.casm_directories.settings_path_crawl(s, config_data["configname"],
                                                                                   self.clex, self.calc_subdir))
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. Extra input file " + s + " not found in CASM project.")
            if config_obj.settings["initial"]:
                extra_input_files += [self.casm_directories.settings_path_crawl(settings["initial"],
                                                                                config_data["configname"],
                                                                                self.clex, self.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No initial INCAR file " + settings["initial"] + " found in CASM project.")
            if config_obj.settings["final"]:
                extra_input_files += [config_obj.casm_directories.settings_path_crawl(settings["final"],
                                                                                      config_data["configname"],
                                                                                      self.clex, self.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No final INCAR file " + config_obj.settings["final"] + " found in CASM project.")
            sys.stdout.flush()

            #make vasp input files
            sample_super_poscarfile = os.path.join(config_obj.calcdir, "00", "POSCAR")
            vasp.io.write_vasp_input(config_data["calcdir"], incarfile, prim_kpointsfile, prim_poscarfile,
                                     sample_super_poscarfile, speciesfile, self.sort, extra_input_files,
                                     settings["strict_kpoints"])
            ## settings the images tag in incar file
            tmp_dict = {"images": config_obj.n_images}
            vasp.io.set_incar_tag(tmp_dict, config_data["calcdir"])

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
