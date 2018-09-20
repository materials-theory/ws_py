from setuptools import setup

setup(
    name='vw-master',
    version='1.0',
    packages=['tests', 'vw_py', 'vw_py.IO', 'vw_py.struc', 'vw_py.cmdline', 'vw_py.Parsers',
              'vw_py.Plotter', 'vw_py.boltztrap', 'vw_py.polyhedron', 'vw_py.generalutils',
              'vw_py.ElectronicStructure', 'examples'],
    url='https://github.com/materials-theory/ws_py',
    license='GPL 3.0',
    author='Woosun Jang',
    author_email='jin890@yonsei.ac.kr',
    description='Python script to post-process VASP generated files', install_requires=['numpy']
)
