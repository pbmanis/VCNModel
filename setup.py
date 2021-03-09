from setuptools import setup, find_packages

version = '0.9.1'

setup(name='vcnmodel',
      version=version,
      description='VCN SBEM Cell modeling',
      url='http://github.com/pbmanis/VCN_Model',
      author='Paul B. Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['vcnmodel*']),
      python_requires='>=3.6',
      # install_requires=['matplotlib>=3.0', 'numpy>=1.1',
#           ],
      zip_safe=False,
      entry_points={
          'console_scripts': [
               'model_run=src.vcnmodel.model_run2:main',
               'allgbcivs=src.vcnmodel.all_gbc_iv:main',
               'show_swc=src.vcnmodel.util.show_swc:main',
               'render=src.vcnmodel.util.render:main',
               'plot_sims=src.vcnmodel.plotters.plot_sims:main',
               'datatable=src.vcnmodel.DataTablesVCN:main',
               'hocswcmap = src.vcnmodel.util.hoc_swc_sectionmap:main',
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
      