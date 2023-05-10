import setuptools


setuptools.setup(
    name="GES-comp-echem",
    version="0.3.2",
    description="",
    long_description="",
    packages=["compechem"],
    package_data={
        "compechem": ["*", "wrappers/*", "engines/*", "functions/*", "tools/*", "core/*"],
    },
    install_requires=[],
)
