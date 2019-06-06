from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='mmodes',
      version='0.2.1',
      description='Metabolic Models based Ordinary Differential Equations Simulation',
      author='Jorge Carrasco Muriel',
      author_email='jorge.cmuriel@alumnos.upm.es',
      keywords='python metabolic-models bioinformatics systems-biology',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'cobra>=0.14.1', 'numpy', 'pandas', 'matplotlib', 'scipy', 'dill'
      ],
      zip_safe=False)
