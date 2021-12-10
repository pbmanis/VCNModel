from setuptools import setup, find_packages

"""
This module is part of *vcnmodel*.

Support::

    NIH grants:
    DC R01DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Paul B. Manis, 2014-2022
"""

version = '0.9.1'

setup(name='vcnmodel',
      version=version,
      description='VCN SBEM Cell modeling',
      url='http://github.com/pbmanis/VCN_Model',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['vcnmodel*']),
      python_requires='>=3.7',
      # install_requires=['matplotlib>=3.0', 'numpy>=1.1',
#           ],
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'model_run=vcnmodel.model_run2:main',
               'allgbcivs=vcnmodel.all_gbc_iv:main',
               'show_swc=scnmodel.util.show_swc:main',
               'render=vcnmodel.util.render:main',
               'plot_sims=vcnmodel.plotters.plot_sims:main',
               'datatable=vcnmodel.DataTablesVCN:main',
               'hocswcmap = vcnmodel.util.hoc_swc_sectionmap:main',
               ],
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
      