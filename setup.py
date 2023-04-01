import os
from setuptools import setup, find_packages


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="cfDNAFE",
    version="0.1.0",
    author="Wanxin Cui, Xiaoqing Peng",
    author_email="xqpeng@csu.edu.cn",
    description="Fragment Feature Extraction For cfDNA Sequencing Data Analysis",
    license="Please see LICENSE.txt.",
    keywords=["cell free DNA", "WGS","WGBS", "Fragmentation", 'mutation', "CNV","methylation"],
    url="https:// ",
    packages=find_packages(),
    package_data={"cfDNAFE": ["data/*", "temp/*"]},
    long_description=read("README.rst"),
    platforms="Linux/Unix",
)