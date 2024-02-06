import setuptools

with open("README.md", "r", encoding="utf-8") as fhand:
    long_description = fhand.read()

setuptools.setup(
    name="mandoline",
    version="0.0.1",
    author="Olivier Chabot",
    description=("Fast approximate slices of AMReX plotfiles"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/olivecha/mandoline",
    project_urls={
        "Bug Tracker": "https://github.com/olivecha/mandoline/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["requests"],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "mandoline = mandoline.cli:main",
            "cuisine = mandoline.cook_cli:main",
        ]
    }
)
