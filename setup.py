import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EBBEF2p",
    version="0.0.1",
    author="Bogdan Teodoru",
    author_email="bteodoru@gmail.com",
    description="Tool for Finite Element Analysis of beams on elastic foundations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bteodoru/ebbef2p",
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'matplotlib', 'itertools', 'math'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    python_requires='>=3.6',
)
