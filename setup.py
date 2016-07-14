from setuptools import setup, find_packages

import os 
def locate_packages():
    packages = ['ld_cmc']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages

setup(
    name="CMC LD extraction",
    version="1.0",
    packages=locate_packages(),
    author="James Boocock and Eli Stahl",
    author_email="james.boocock@otago.ac.nz",
    description="CMC LD extraction software", 
    license="Mit",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'cmc_ld_extract = ld_cmc.main:main',
        ]
        },
    url="github.com/smilefreak/fine_mapping_pipelin",
    use_2to3=True,
)

