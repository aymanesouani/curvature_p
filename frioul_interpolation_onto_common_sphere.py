import commands

if __name__ == '__main__':
    db_paths = ['/hpc/scalp/data/REPRO_database/FS_database_OASIS_test_retest_FS5.3.0', '/hpc/scalp/data/REPRO_database/FS_database_KKI_test_retest_FS6.0']

    curv_types = ['mean_Rusinkiewicz', 'mean_Dong', 'mean_Maillot', 'mean_PetitJean', 'mean_Peyre', 'mean_Taubin',
                  'gauss_Dong', 'gauss_Maillot', 'gauss_PetitJean', 'gauss_Peyre', 'gauss_Rusinkiewicz', 'gauss_Taubin']

    for db_path in db_paths:
        for curv_type in curv_types:
            print(curv_type)
            cmd = "frioul_batch 'export PYTHONPATH=/hpc/meca/softs/dev/auzias/pyhon/slam/; /hpc/meca/users/auzias/miniconda3/envs/trimesh_dev/bin/python /hpc/meca/softs/dev/auzias/pyhon/curvatures_compare/interpolation_onto_common_sphere.py %s %s'" %(db_path, curv_type)
            print(cmd)
            a = commands.getoutput(cmd)
            print(a)