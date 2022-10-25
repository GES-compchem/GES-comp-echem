import setuptools


setuptools.setup(
    name="GES-comp-echem",
    version="0.2.1",
    description="",
    long_description="",
    packages=["compechem"],
    package_data={
        "compechem": ["wrappers/*", "functions/*", "tools/*"],
    },
    install_requires=[],
)
