import setuptools


setuptools.setup(
    name="GES-comp-echem",
    version="0.2.2",
    description="",
    long_description="",
    packages=["compechem"],
    package_data={
        "compechem": ["*", "wrappers/*", "functions/*", "tools/*", "core/*"],
    },
    install_requires=[],
)
