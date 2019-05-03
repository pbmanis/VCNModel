from setuptools import setup, find_packages
import os

# Use Semantic Versioning, http://semver.org/
version_info = (0, 1, 0, '')
__version__ = '%d.%d.%d%s' % version_info


setup(name='VCN',
      version=__version__,
      description='VCN SBEM Cell modeling',
      url='http://github.com/pbmanis/VCN_Model',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['src*']),
      install_requires=['matplotlib>=3.0', 'numpy>=1.1',
          ],
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'model_run=src.model_run:main',
               'allgbcivs=src.all_gbc_ivs:main',
               ],
          # 'gui_scripts': [
          #       'event_monger=src.event_monger:main',
          # ]
      },
      classifiers = [
             "Programming Language :: Python :: 3.6+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Manis Lab",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Computational Modeling :: Neuroscience",
             ],
    )
      