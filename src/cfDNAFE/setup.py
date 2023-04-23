import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="CFFE",
    version="0.0.1",
    author="Wanxin Cui, Xiaoqing Peng",
    author_email="wxCui2020@csu.edu.cn",
    description="Fragment Feature Extraction For cfDNA Sequencing Data Analysis",
    license="Please see LICENSE.txt.",
    keywords=["cell free DNA", "WGS", "Fragmentation", 'mutation', "CNV", 'DNA methylation'],
    url="https:// ",
    packages=find_packages(),
    package_data={"cfDNAFFE": ["data/*", "temp/*"]},
    long_description=read("README.rst"),
    platforms="Linux/Unix",
)