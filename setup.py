from setuptools import find_packages, setup
setup(
    name='CodeA',
    packages=find_packages(),
    version='1.0',
    description='Dusty Particles Dynamics',
    author='Morris',
    license='MIT',
    install_requires=['numpy', 'scipy', 'matplotlib', 'tqdm']
)