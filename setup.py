from setuptools import setup
from setuptools import find_packages

with open("README.md","r") as fh:
    long_description = fh.read()
    
setup(name='NaivePolyCrys',
      version='0.1',
      description='Generate polycrystalline atomic configurations',
      long_description = long_description,
      long_description_content_type='text/markdown',
      url='',
      author='Stefan Bringuier',
      author_email='stefanbringuier@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          'numpy',
          'ase'],
      zip_safe=False,
      python_requires='>3.7.5')
