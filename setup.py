from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
import numpy

"""
This setup.py is for *vcnmodel*.

Support::

    NIH grants:
    DC R01DC015901 (Spirou, Manis, Ellisman),
    DC R01 DC004551 (Manis, 2013-2019, Early development)
    DC R01 DC019053 (Manis, 2020-2025, Later development)

Paul B. Manis, 2014-2022
"""

version = '0.9.9'  # 30 August 2022
extensions = [
    Extension("sttc_cython",  ["vcnmodel/analyzers/sttc_cython.pyx"],
               include_dirs=[numpy.get_include()])
]
setup(name='vcnmodel',
      version=version,
      description='VCN SBEM Cell modeling',
      url='http://github.com/pbmanis/VCNModel',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['vcnmodel*']),
      ext_modules=cythonize(extensions),

      python_requires='>=3.9',
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'model_run=vcnmodel.model_run2:main',
               'allgbcivs=vcnmodel.all_gbc_iv:main',
               'show_swc=vcnmodel.util.show_swc:main',
               'render=vcnmodel.util.render:main',
               'plot_sims=vcnmodel.plotters.plot_sims:main',
               'datatable=vcnmodel.DataTablesVCN:main',
               'hocswcmap = vcnmodel.util.hoc_swc_sectionmap:main',
               ],
      },
      classifiers = [
             "Programming Language :: Python :: 3.9+",
             "Development Status ::  Beta",
             "Environment :: Console",
             "Intended Audience :: Manis Lab",
             "License :: MIT",
             "Operating System :: OS Independent",
             "Topic :: Software Development :: Tools :: Python Modules",
             "Topic :: Computational Modeling :: Neuroscience",
             ],
    )
