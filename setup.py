from setuptools import setup

setup(
    name='ws_py',
    version='1.0',
    packages=['IO', 'struc', 'generalutils', 'cmdline', 'Parsers', 'Plotter',
              'boltztrap', 'polyhedron', 'ElectronicStructure'],
    url='https://github.com/materials-theory/ws_py',
    license='GPL v3',
    author='Woosun Jang',
    author_email='jin890@yonsei.ac.kr',
    description='A set of python libraries to handle the VASP input / output files', install_requires=['pandas',
                                                                                                       'numpy']
)