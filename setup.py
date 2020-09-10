import os

import versioneer
from setuptools import setup

# Create list of data files


def find_data_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join("..", path, filename))

    return paths


extra_files = find_data_files("nazgul/data")

stan_files = [f for f in find_data_files("nazgul/stan_models") if ".stan" in f]

extra_files.extend(stan_files)

setup(
    version=versioneer.get_version(),
    include_package_data=True,
    package_data={"": extra_files},
    cmdclass=versioneer.get_cmdclass(),

)
