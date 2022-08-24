import setuptools


setuptools.setup(
    name="GES-comp-echem",
    version="0.1.38a",
    description="",
    long_description="",
    packages=["compechem"],
    package_data={
        "compechem": ["algorithms/*", "wrappers/*", "functions/*", "tools/*"],
    },
    install_requires=[],
)
