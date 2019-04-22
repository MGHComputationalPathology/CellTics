"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import


from setuptools import setup, find_packages
exec(open('celltics/version.py').read())


def load_requirements():
    with open('requirements.txt', 'r') as f:
        return [l.strip(' ') for l in f if '://' not in l]

setup(
    name='celltics',
    version='__version__',
    py_modules=['celltics'],
    packages=find_packages(),
    url='https://github.com/MGHComputationalPathology/CellTics',
    license='BSD',
    author='Allison MacLeay',
    author_email='amacleay@mgh.harvard.edu',
    description='Center for Integrated Diagnostics at Mass General Hospital NGS tools',
    include_package_data=True,
    package_data={
        'celltics.tests.data.files': ['*.*'],
    },
    install_requires=load_requirements(),
    entry_points='''
        [console_scripts]
        celltics=celltics.main:cli
    ''',
)
