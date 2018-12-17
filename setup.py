from setuptools import setup, find_packages

setup(
    name='cta_transient_search',
    author='Lena Linhoff, Jana Moschner, Kai Brügge',
    author_email='lena.linhoff@tu-dortmund.de',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        'astropy',
        'numpy',
        'matplotlib',
        'click',
    ],
)
