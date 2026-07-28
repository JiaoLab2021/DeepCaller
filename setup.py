from setuptools import setup, find_packages
from deepcaller import __version__, __author__, __email__, __license__, __description__

setup(
    name="DeepCaller",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description=__description__,
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    license=__license__,
    url="https://github.com/JiaoLab2021/DeepCaller",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "pysam",
        "numpy",
        "pandas",
        "h5py",
        "tensorflow",
        "tensorflow-addons",
        "pyarrow",
        "setproctitle",
    ],
    entry_points={
        "console_scripts": [
            "deepcaller=deepcaller.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
    ],
)