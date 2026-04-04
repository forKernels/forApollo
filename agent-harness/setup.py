"""CLI-Anything agent harness for forApollo -- guidance, navigation, control."""

from setuptools import setup, find_namespace_packages

setup(
    name="cli-anything-forapollo",
    version="0.1.0",
    description="CLI-Anything agent harness for forApollo (guidance, navigation, control)",
    author="The Fantastic Planet",
    packages=find_namespace_packages(include=["cli_anything.*"]),
    install_requires=["click>=8.0", "numpy"],
    entry_points={
        "console_scripts": [
            "cli-anything-forapollo=cli_anything.forapollo.__main__:main",
        ],
    },
    python_requires=">=3.10",
)
