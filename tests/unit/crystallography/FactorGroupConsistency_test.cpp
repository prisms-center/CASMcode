#include "gtest/gtest.h"
#include <fstream>
#include <vector>
#include "autotools.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/eigen.hh"

using namespace CASM;

class CrystalGroupTest : public testing::Test
{
    protected:
        void add_primitive_structure(std::string prim_file_in_crystallography_testdir, int factor_group_size)
        {
            auto test_files_dir=autotools::abs_srcdir() + "/tests/unit/crystallography/";
            std::ifstream s;

            s.open(test_files_dir+prim_file_in_crystallography_testdir);
            primitive_structures.emplace_back(xtal::BasicStructure::from_poscar_stream(s));
            expected_primitive_factor_group_size.emplace_back(factor_group_size);
            s.close();

            return;
        }

        void SetUp() override
        {
            add_primitive_structure("POS1.txt", 16);
            add_primitive_structure("PRIM4.txt", 48);
            add_primitive_structure("hcp_mg.vasp", 24);

            Eigen::Matrix3l super_mat_0;
            super_mat_0 << 1,0,0,0,1,0,0,0,4;
            prim_to_super_trasnformation_matrices.emplace_back(std::move(super_mat_0));

            Eigen::Matrix3l super_mat_1;
            super_mat_1 << -1,1,1,1,-1,1,1,1,-1;
            prim_to_super_trasnformation_matrices.emplace_back(std::move(super_mat_1));

            for(const Eigen::Matrix3l& mat : prim_to_super_trasnformation_matrices)
            {
                std::vector<xtal::BasicStructure> superstructures_from_mat;
                for(const xtal::BasicStructure& prim : primitive_structures)
                {
                    superstructures_from_mat.emplace_back(xtal::make_superstructure(prim,mat));
                }
                superstructures.emplace_back(std::move(superstructures_from_mat));
            }

            return;
        }

        std::vector<xtal::BasicStructure> primitive_structures;
        std::vector<int> expected_primitive_factor_group_size;

        //Each inner vector is all the primitive structures with the same transformation applied
        //to each primitive structure
        std::vector<std::vector<xtal::BasicStructure>> superstructures;
        std::vector<Eigen::Matrix3l> prim_to_super_trasnformation_matrices;
};

TEST_F(CrystalGroupTest, FactorGroupSizes)
{
    for(int i=0; i<primitive_structures.size(); ++i)
    {
        auto factor_group=xtal::make_factor_group(primitive_structures[i]);
        int expected_size=expected_primitive_factor_group_size[i];
        EXPECT_EQ(factor_group.size(),expected_size);

        for(int j=0; j<prim_to_super_trasnformation_matrices.size(); ++j)
        {
            const xtal::BasicStructure& superstruc=superstructures[i][j];
            int prims_in_superstruc=prim_to_super_trasnformation_matrices[j].determinant();

            auto super_factor_group=xtal::make_factor_group(superstruc);
            EXPECT_EQ(super_factor_group.size(),expected_size*prims_in_superstruc);
        }
    }
}

/* { */
/*   xtal::SymOpVector fgroup; */
/*     xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc))); */
/*     auto sym_set = mapper.map_deformed_struc_impose_lattice(sstruc, xtal::Lattice(sstruc.lat_column_mat), 0, xtal::StrucMapping::big_inf(), 1e-3); */

/*     EXPECT_EQ(sym_set.size(), N); */
/*     fgroup = adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set); */
/* } */
