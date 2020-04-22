import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ebbef2p-python",
    version="0.0.2.dev1",
    author="Bogdan Teodoru",
    author_email="bteodoru@gmail.com",
    description="Tool for Finite Element Analysis of beams on elastic foundations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bteodoru/ebbef2p-python",
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'matplotlib'],
    classifiers=[
        # How mature is this project? Common values are
	#Development Status :: 1 - Planning
	#Development Status :: 2 - Pre-Alpha
	#Development Status :: 3 - Alpha
	#Development Status :: 4 - Beta
	#Development Status :: 5 - Production/Stable
	#Development Status :: 6 - Mature
	#Development Status :: 7 - Inactive
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    keywords='finite element method, math, numerics, elastic foundation',
    python_requires='>=3.6',
)
