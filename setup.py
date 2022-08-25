import setuptools


setuptools.setup(
    name="GES-comp-echem",
    version="0.2.0",
    description="",
    long_description="",
    packages=["compechem"],
    package_data={
        "compechem": ["algorithms/*", "wrappers/*", "functions/*", "tools/*"],
    },
    install_requires=[],
)
