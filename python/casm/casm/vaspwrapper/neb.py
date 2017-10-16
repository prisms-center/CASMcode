import os, math, sys, json, re, warnings, shutil
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from casm.vaspwrapper import VaspCalculatorBase
import vasp.Neb

class Neb(VaspCalculatorBase):
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
        self.results_subdir = '01'
        self.calculator = vasp.Neb

    @classmethod
    def neb(cls, configname, calctype, auto=True, sort=True):
        sel = Selection.selection_from_confignames([configname])
        obj = cls(sel, calctype, auto, sort)
        return obj

    def config_properties(self, config_data):
        """configuration properties as a dict"""
        config_dict = super(Neb, self).config_properties(config_data)
        try:
            n_images = json.load(config_data["setfile"])["n_images"]
        except:
            ## Error message if "n_images" not present in settings
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()
        config_data["n_images"] = n_images
        config_data["calcdir"] = os.path.join(config_data["calcdir"], "N_images_{}".format(n_images))
        return config_dict

    def pre_setup(self):
        proj = self.selection.proj
        dict = {}
        for config_data in self.selection.data:
            conf_dict = {"n_images" : config_data["n_images"],
                         "calctype" : config_data["calctype"]}
            dict[config_data["configname"]] = conf_dict
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

    def setup(self):
        """ Setup initial relaxation run for the selection
        """
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
        vaspfiles = super(Neb, self).get_vasp_input_files(config_data, settings)
        incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files = vaspfiles
        super_poscarfile = os.path.join(config_data["calcdir"], "00", "POSCAR")
        return incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files

    def submit(self):
        super(Neb, self).submit()

    def run(self):
        super(Neb, self).run()

    def report(self):
        super(Neb, self).report()

    @staticmethod
    def run_cmd(configdir, calctype):
        return "python -c \"import casm.vaspwrapper; obj = casm.vaspwrapper.Neb.neb('{0}', '{1}'); obj.run()\"\n".format(configdir, calctype)

    def finalize(self, config_data, super_poscarfile=None):
        if super_poscarfile is None:
            super_poscarfile = os.path.join(config_data["calcdir"], "00/POSCAR")
        super(Neb, self).finalize(config_data, super_poscarfile)
        all_image_folders = [int(i.strip().split('_')[-1]) for i in os.listdir(self.config_obj.calcdir) if "N_images" in i]
        num_images = config_data["n_images"]
        if num_images == max(all_image_folders):
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
        final_output = []
        num_images = vasp.io.get_incar_tag("IMAGES", calcdir)
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            vaspdir = calcdir + "{}".format(img)
            output = super(Neb, self).properties(vaspdir, super_poscarfile, speciesfile)
            output["Image_number"] = img
            final_output.append(output)
        return final_output
