from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
        name='dtailor',
        version='1.0',
        author='Joao Guimaraes',
        description='D-Tailor',
        packages=find_packages(),
        package_dir={'dtailor': 'dtailor'},
        install_requires=required
)
