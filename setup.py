import setuptools

import os.path

#readme = ""
#here = os.path.abspath(os.path.dirname(__file__))
#readme_path = os.path.join(here, "README.md")
#if os.path.exists(readme_path):
#    with open(readme_path, "rb") as stream:
#        readme = stream.read().decode("utf8")


setuptools.setup(
    name="SpIn",
    version="1.0.0",
    author="Camila Angulo, Oscar Chávez, Erica Sawczynec",
    author_email="spin@spin.com",
    description="A Python package for spectra information",
    url="https://github.com/ericasaw/SpIn",
    packages=setuptools.find_packages(),
    python_requires=">=3.8",
)