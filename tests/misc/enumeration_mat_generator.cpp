#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/container/LinearAlgebra.hh"

#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/external/Eigen/Dense"
#include "casm/crystallography/Niggli.hh"

using namespace CASM;

void lat_generator(const std::string &pos_filename, int minvol, int maxvol)
{
    const boost::filesystem::path testdir("../unit/crystallography");

    const Structure test_struc(testdir/pos_filename);
    const Lattice test_lat=test_struc.lattice();
    const SymGroup effective_pg=test_struc.factor_group();

    Array<Lattice> enumerated_lats;
    test_lat.generate_supercells(enumerated_lats, effective_pg, minvol, maxvol, 3, Eigen::Matrix3i::Identity());

    jsonParser lat_dump;
    lat_dump["min_vol"]=minvol;
    lat_dump["max_vol"]=maxvol;
    lat_dump["source"]=pos_filename;
    lat_dump["lats"]=enumerated_lats;

    lat_dump.write(testdir.string()+"/"+pos_filename+"_"+std::to_string(minvol)+"_"+std::to_string(maxvol)+"_lats.json");

    return;
}

void mat_generator(const std::string &pos_filename, int minvol, int maxvol)
{
    const boost::filesystem::path testdir("../unit/crystallography");

    const Structure test_struc(testdir/pos_filename);
    const Lattice test_lat=test_struc.lattice();
    const SymGroup effective_pg=test_struc.factor_group();

    Array<Eigen::Matrix3i> enumerated_mats;

    SupercellEnumerator<Lattice> test_enumerator(test_lat, effective_pg, minvol, maxvol);

    for(auto it = test_enumerator.begin(); it != test_enumerator.end(); ++it)
    {
        enumerated_mats.push_back(it.matrix());
    }

    jsonParser mat_dump;
    mat_dump["min_vol"]=minvol;
    mat_dump["max_vol"]=maxvol;
    mat_dump["source"]=pos_filename;
    mat_dump["mats"]=enumerated_mats;

    mat_dump.write(testdir.string()+"/"+pos_filename+"_"+std::to_string(minvol)+"_"+std::to_string(maxvol)+"_mats.json");

    return;
}

int main()
{
    mat_generator("POS1",1,6);
    mat_generator("PRIM1",2,9);
    mat_generator("PRIM2",4,7);
    mat_generator("PRIM4",1,8);

    lat_generator("POS1", 2, 6);
    lat_generator("PRIM1", 2, 9);
    lat_generator("PRIM2", 3, 7);
    lat_generator("PRIM4", 1, 8);
    lat_generator("PRIM5", 1, 8);

    return 0;
}
